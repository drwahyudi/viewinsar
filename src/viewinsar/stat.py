"""Block statistics along a line - the core of `viewinsar stat`.

Draw a line along a fault (as many vertices as needed). Every SPACING km along
it a point is placed, and at each point two circles of RADIUS km: one on each
side of the line, its edge GAP km (default 1) clear of the line. Block 1 is
the LEFT side of the line, walking from its first vertex to its last; block 2
the right. The values inside each circle describe that block at that point.

    the circle centre = point + / - (RADIUS + GAP) * the line's local normal

The normal is the line's direction over +/- RADIUS around the point, so a
hand-digitised trace with small wiggles still gives steady circles. On the
inside of a bend, the centre can come closer than RADIUS + GAP to ANOTHER part
of the line; the circle is then pushed further out along the same normal until
it clears every part of the line by GAP - or dropped if it cannot (or if the
way out crosses the line). So no circle ever reaches across, or within GAP of,
the line: every pixel used is at least GAP from it, on its own side.

Everything is computed in a metric frame of the data grid (X along the
columns, Y along -rows; MintPy's ground pixel size on a radar grid, a local
equirectangular frame on a lat/lon one), so the radius is a true distance. A
grid that is a mirror image of the map - every right-looking radar grid, which
is why view.py flips it - gets X along -columns instead, so "left" and "right"
are always those of the ground (and of the map as view.py shows it).
Numpy only: no MintPy, no matplotlib.
"""

import numpy as np

from .swath import get_pixel_size, is_projected

__all__ = ['grid_mirrored', 'polyline_frame', 'from_yx', 'to_yx', 'sample_points', 'distance_to_polyline', 'side_of_polyline',
           'ray_hits_polyline', 'place_circles', 'circle_pixels', 'circle_outline',
           'stat_line', 'summarize', 'kde_mode']

CHUNK = 2_000_000          # points x segments handled at once (memory bound)


##############################################################################
#                       the line, in metres                                  #
##############################################################################
def _lookup_det(lat, lon, max_points=250_000):
    """Sign source: det d(lon, lat) / d(col, row) over a radar grid, or None."""
    import warnings

    lat, lon = np.asarray(lat, dtype=np.float64), np.asarray(lon, dtype=np.float64)
    if lat.ndim != 2 or min(lat.shape) < 2:
        return None
    step = max(1, int(np.ceil(np.sqrt(lat.size / max_points))))
    la, lo = lat[::step, ::step], lon[::step, ::step]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)       # all-NaN rows / columns
        lon_c, lat_c = np.nanmedian(np.diff(lo, axis=1)), np.nanmedian(np.diff(la, axis=1))
        lon_r, lat_r = np.nanmedian(np.diff(lo, axis=0)), np.nanmedian(np.diff(la, axis=0))
    det = lon_c * lat_r - lon_r * lat_c
    return float(det) if np.isfinite(det) and det != 0 else None


def grid_mirrored(atr, lookup=None):
    """Is the (column, -row) frame of this grid a mirror image of the map?

    Geocoded: when X_STEP and Y_STEP have the same sign (north-up grids have
    X_STEP > 0 > Y_STEP). Radar-coded: from the latitude / longitude of
    geometryRadar.h5 (`lookup`, anything with .lat and .lon arrays) when
    given, else from the look side: the grid of a right-looking SAR is always
    a mirror image of the ground (ANTENNA_SIDE = -1, MintPy's default).
    """
    if 'Y_FIRST' in atr:
        return float(atr['X_STEP']) * float(atr['Y_STEP']) > 0
    if lookup is not None:
        det = _lookup_det(lookup.lat, lookup.lon)
        if det is not None:
            # (col, row) -> (X, Y) = (col, -row) flips; so the frame keeps the
            # orientation of the ground exactly when (col, row) -> (lon, lat) flips too
            return det > 0
    return str(atr.get('ANTENNA_SIDE', '-1')).strip() != '1'


def polyline_frame(atr, vertices_yx, mirrored=None):
    """Metric frame of a polyline drawn on a data grid.

    Parameters: atr         - dict, metadata of the data file
                vertices_yx - list of [row, col], full file, at least 2 distinct
                mirrored    - bool, the grid is a mirror image of the map
                              (default: grid_mirrored(atr))
    Returns:    dict with P (vertices in the frame), the pixel size, the
                segment vectors and lengths, the cumulative length, the
                frame -> km scale, and sx (-1 on a mirrored grid, else +1)
    """
    yx = np.asarray(vertices_yx, dtype=np.float64).reshape(-1, 2)
    keep = np.r_[True, np.any(np.diff(yx, axis=0) != 0, axis=1)]   # no repeated clicks
    yx = yx[keep]
    if len(yx) < 2:
        raise ValueError('the line needs at least two different points')

    lat_c = None
    if 'Y_FIRST' in atr and not is_projected(atr):
        lat_c = float(atr['Y_FIRST']) + float(atr['Y_STEP']) * (yx[:, 0].mean() + 0.5)
    x_step, y_step, unit = get_pixel_size(atr, lat_c=lat_c)
    if mirrored is None:
        mirrored = grid_mirrored(atr)
    sx = -1.0 if mirrored else 1.0
    pts = np.c_[sx * yx[:, 1] * x_step, -yx[:, 0] * y_step]
    seg = np.diff(pts, axis=0)
    seg_len = np.hypot(seg[:, 0], seg[:, 1])
    return dict(
        yx=yx, P=pts, seg=seg, seg_len=seg_len, cum=np.r_[0.0, np.cumsum(seg_len)],
        x_step=x_step, y_step=y_step, unit=unit, sx=sx, mirrored=bool(mirrored),
        scale=1e-3 if unit == 'm' else 1.0,                 # frame unit -> km
        dist_unit='km' if unit == 'm' else 'pixel',
    )


def from_yx(frame, rows, cols):
    """(row, col) of the full file -> frame coordinates, as an (..., 2) array."""
    rows, cols = np.asarray(rows, dtype=np.float64), np.asarray(cols, dtype=np.float64)
    return np.stack([frame['sx'] * cols * frame['x_step'], -rows * frame['y_step']], axis=-1)


def to_yx(frame, xy):
    """Frame coordinates -> fractional (row, col) of the full file."""
    xy = np.asarray(xy, dtype=np.float64)
    return -xy[..., 1] / frame['y_step'], xy[..., 0] / (frame['sx'] * frame['x_step'])


def position_along(frame, dist):
    """Frame coordinates of the points `dist` (frame units) along the line."""
    cum, seg, seg_len = frame['cum'], frame['seg'], frame['seg_len']
    dist = np.clip(np.atleast_1d(np.asarray(dist, dtype=np.float64)), 0.0, cum[-1])
    idx = np.clip(np.searchsorted(cum, dist, side='right') - 1, 0, len(seg) - 1)
    frac = (dist - cum[idx]) / np.maximum(seg_len[idx], 1e-300)
    return frame['P'][idx] + frac[:, None] * seg[idx]


def sample_points(frame, spacing_km, smooth_km=0.0):
    """Points every `spacing_km` along the line, centred on it.

    floor(L / spacing) + 1 points, spaced exactly `spacing_km` apart and
    centred, so the first and the last sit equally far from the line's ends.
    The normal at a point is that of the chord between the points `smooth_km`
    before and after it along the line (the local segment if 0; at a vertex,
    the bisector of both segments).

    Returns: dist_km - 1D array, distance of each point from the first vertex
             points  - (n, 2) array, in the frame
             normals - (n, 2) array, unit normals to the LEFT of the line
    """
    cum, seg, seg_len = frame['cum'], frame['seg'], frame['seg_len']
    total, step = cum[-1], spacing_km / frame['scale']
    if step <= 0:
        raise ValueError('the spacing must be positive')
    num = int(np.floor(total / step + 1e-9)) + 1
    dist = (total - (num - 1) * step) / 2.0 + np.arange(num) * step
    points = position_along(frame, dist)

    # the local direction: the segment, or the bisector on a vertex
    idx = np.clip(np.searchsorted(cum, dist, side='right') - 1, 0, len(seg) - 1)
    unit = seg / np.maximum(seg_len, 1e-300)[:, None]
    local = unit[idx].copy()
    for k, j in enumerate(idx):
        if j > 0 and abs(dist[k] - cum[j]) < 1e-9 * max(total, 1.0):
            both = unit[j - 1] + unit[j]
            if np.hypot(*both) > 1e-12:
                local[k] = both / np.hypot(*both)

    tangent = local
    half = max(float(smooth_km), 0.0) / frame['scale']
    if half > 0:
        chord = position_along(frame, dist + half) - position_along(frame, dist - half)
        size = np.hypot(chord[:, 0], chord[:, 1])
        # keep the local direction where the chord folds back (a hairpin)
        ok = (size > 1e-9 * half) & (np.einsum('ij,ij->i', chord, local) > 0)
        tangent = np.where(ok[:, None], chord / np.maximum(size, 1e-300)[:, None], local)
    normals = np.c_[-tangent[:, 1], tangent[:, 0]]
    return dist * frame['scale'], points, normals


def _nearest_on_segments(frame, xy):
    """t (unclipped), distance to every segment - for a chunk of points."""
    a, ab = frame['P'][:-1], frame['seg']
    ab2 = np.maximum(np.einsum('ij,ij->i', ab, ab), 1e-300)
    rel = xy[:, None, :] - a[None, :, :]
    t_raw = np.einsum('kij,ij->ki', rel, ab) / ab2
    t = np.clip(t_raw, 0.0, 1.0)
    diff = rel - t[..., None] * ab[None]
    return t_raw, np.hypot(diff[..., 0], diff[..., 1]), rel


def _chunks(frame, num):
    size = max(1, CHUNK // max(len(frame['seg']), 1))
    return [slice(i, min(i + size, num)) for i in range(0, num, size)]


def distance_to_polyline(frame, xy):
    """Shortest distance (frame units) from points to any part of the line."""
    xy = np.atleast_2d(np.asarray(xy, dtype=np.float64)).reshape(-1, 2)
    out = np.empty(len(xy))
    for sl in _chunks(frame, len(xy)):
        out[sl] = _nearest_on_segments(frame, xy[sl])[1].min(axis=1)
    return out


def side_of_polyline(frame, xy):
    """+1 for points left of the line (walking first -> last vertex), -1 right.

    The side of the nearest segment; where the nearest point is a vertex
    between two segments, the side of that vertex's pseudo-normal (the sum of
    both segment normals) - otherwise a point just outside a sharp bend could
    be put on the inside.
    """
    xy = np.atleast_2d(np.asarray(xy, dtype=np.float64)).reshape(-1, 2)
    ab, P = frame['seg'], frame['P']
    unit = ab / np.maximum(frame['seg_len'], 1e-300)[:, None]
    normal = np.c_[-unit[:, 1], unit[:, 0]]
    n_seg = len(ab)
    out = np.empty(len(xy), dtype=np.int64)
    for sl in _chunks(frame, len(xy)):
        t_raw, dist, rel = _nearest_on_segments(frame, xy[sl])
        k = np.arange(len(dist))
        j = np.argmin(dist, axis=1)
        side = ab[j, 0] * rel[k, j, 1] - ab[j, 1] * rel[k, j, 0]
        tj = t_raw[k, j]
        for on_vertex, vert, s0, s1 in (((tj >= 1) & (j < n_seg - 1), j + 1, j, j + 1),
                                        ((tj <= 0) & (j > 0), j, j - 1, j)):
            if on_vertex.any():
                pseudo = normal[s0[on_vertex]] + normal[s1[on_vertex]]
                w = xy[sl][on_vertex] - P[vert[on_vertex]]
                side[on_vertex] = np.einsum('ij,ij->i', w, pseudo)
        out[sl] = np.where(side >= 0, 1, -1)
    return out


def ray_hits_polyline(frame, start, direction):
    """Distance along a ray from `start` to the first point of the line it meets.

    `direction` is a unit vector; np.inf if the ray never meets the line.
    """
    a, ab = frame['P'][:-1], frame['seg']
    d = np.asarray(direction, dtype=np.float64)
    rel = a - np.asarray(start, dtype=np.float64)
    denom = d[0] * ab[:, 1] - d[1] * ab[:, 0]
    ok = np.abs(denom) > 1e-12 * np.maximum(np.hypot(ab[:, 0], ab[:, 1]), 1e-300)
    with np.errstate(divide='ignore', invalid='ignore'):
        s = (rel[:, 0] * ab[:, 1] - rel[:, 1] * ab[:, 0]) / denom
        u = (rel[:, 0] * d[1] - rel[:, 1] * d[0]) / denom
    hit = ok & (s >= 0) & (u >= -1e-12) & (u <= 1 + 1e-12)
    return float(s[hit].min()) if hit.any() else np.inf


def place_circles(frame, points, normals, radius_km, gap_km=1.0, max_push=3.0):
    """Circle centres on both sides of the line, each GAP clear of all of it.

    A centre starts RADIUS + GAP out from its point along the normal - to the
    left for block 1, to the right for block 2 - and moves out in steps of 2 %
    until (1) no part of the line is closer than RADIUS + GAP, (2) it is on its
    own side, and (3) the way out from the point does not cross the line. No
    such spot within `max_push` x (RADIUS + GAP): NaN, the circle is left out.

    Returns: centres   - (2, n, 2) array: [0] block 1 (left), [1] block 2
                         (right); NaN where no clear spot was found
             offset_km - (2, n) array, distance of each centre from its point
    """
    radius, gap = radius_km / frame['scale'], gap_km / frame['scale']
    need = radius + gap
    offs = need * np.linspace(1.0, max_push, int(round((max_push - 1.0) / 0.02)) + 1)
    centres = np.full((2, len(points), 2), np.nan)
    offsets = np.full((2, len(points)), np.nan)
    for b, sign in enumerate((1, -1)):
        for k in range(len(points)):
            out = sign * normals[k]
            start = points[k] + out * need * 1e-6
            reach = ray_hits_polyline(frame, start, out)
            cand = points[k] + offs[:, None] * out
            good = offs < reach
            if not good.any():
                continue
            good[good] &= distance_to_polyline(frame, cand[good]) >= need * (1 - 1e-9)
            if good.any():
                good[good] &= side_of_polyline(frame, cand[good]) == sign
            if good.any():
                i = int(np.argmax(good))
                centres[b, k], offsets[b, k] = cand[i], offs[i] * frame['scale']
    return centres, offsets


##############################################################################
#                       the pixels in a circle                               #
##############################################################################
def circle_pixels(data, box, frame, centre, radius_km, mask=None, keep_zero=False):
    """Valid pixels whose centre lies within `radius_km` of `centre`.

    Parameters: data   - 2D array, read within `box`
                box    - (x0, y0, x1, y1) of that window, in full-file pixels
                frame  - from polyline_frame()
                centre - (2,) frame coordinates; NaN gives nothing
    Returns:    values, rows, cols (full file) - 1D arrays; clipped - bool, the
                circle reaches past the window that was read
    """
    empty = np.array([]), np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    if not np.all(np.isfinite(centre)):
        return (*empty, False)
    radius = radius_km / frame['scale']
    row_c, col_c = to_yx(frame, centre)
    d_row, d_col = radius / frame['y_step'], radius / frame['x_step']
    length, width = data.shape[-2:]
    r0, r1 = int(np.floor(row_c - d_row)), int(np.ceil(row_c + d_row))
    c0, c1 = int(np.floor(col_c - d_col)), int(np.ceil(col_c + d_col))
    rr0, rr1 = max(r0, box[1]), min(r1, box[1] + length - 1)
    cc0, cc1 = max(c0, box[0]), min(c1, box[0] + width - 1)
    clipped = bool(rr0 > r0 or rr1 < r1 or cc0 > c0 or cc1 < c1)
    if rr1 < rr0 or cc1 < cc0:
        return (*empty, clipped)

    rows, cols = np.arange(rr0, rr1 + 1), np.arange(cc0, cc1 + 1)
    dx = frame['sx'] * cols * frame['x_step'] - centre[0]
    dy = -rows * frame['y_step'] - centre[1]
    inside = dy[:, None] ** 2 + dx[None, :] ** 2 <= radius ** 2
    idx = np.ix_(rows - box[1], cols - box[0])
    sub = data[idx]
    val = np.ma.getdata(sub).astype(np.float64)
    good = inside & np.isfinite(val)
    if np.ma.isMaskedArray(sub):
        good &= ~np.ma.getmaskarray(sub)
    if mask is not None:
        good &= np.asarray(mask)[idx].astype(bool)
    if not keep_zero:
        good &= val != 0
    ri, ci = np.nonzero(good)
    return val[ri, ci], rows[ri], cols[ci], clipped


def circle_outline(frame, centre, radius_km, num=73):
    """(rows, cols) of a circle's outline, for drawing it on the map."""
    ang = np.linspace(0, 2 * np.pi, num)
    radius = radius_km / frame['scale']
    xy = np.c_[centre[0] + radius * np.cos(ang), centre[1] + radius * np.sin(ang)]
    return to_yx(frame, xy)


def kde_mode(values, bins=2048):
    """The mode of a sample: the peak of its kernel density - the widest point
    of its violin.

    A Gaussian kernel with Scott's bandwidth (std * n ** -0.2), as in
    matplotlib's violinplot and scipy's gaussian_kde, evaluated on `bins`
    steps between the smallest and the largest value. Fewer than 3 values, or
    all equal: the median.
    """
    values = np.asarray(values, dtype=np.float64).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    lo, hi = values.min(), values.max()
    width = values.std(ddof=1) * values.size ** -0.2 if values.size > 1 else 0.0
    if values.size < 3 or hi <= lo or not width > 0:
        return float(np.median(values))
    hist, edges = np.histogram(values, bins=bins, range=(lo, hi))
    step = edges[1] - edges[0]
    half = int(np.ceil(4.0 * width / step))
    kernel = np.exp(-0.5 * (np.arange(-half, half + 1) * step / width) ** 2)
    density = np.convolve(hist, kernel, mode='full')[half:half + bins]
    return float((edges[:-1] + edges[1:])[np.argmax(density)] / 2.0)


def summarize(values):
    """count, median, mean, mode, std, q25, q75, min, max - NaN when empty."""
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return dict(count=0, median=np.nan, mean=np.nan, mode=np.nan, std=np.nan, q25=np.nan,
                    q75=np.nan, min=np.nan, max=np.nan)
    q25, med, q75 = np.percentile(values, [25, 50, 75])
    return dict(count=int(values.size), median=float(med), mean=float(values.mean()),
                mode=kde_mode(values), std=float(values.std()), q25=float(q25),
                q75=float(q75), min=float(values.min()), max=float(values.max()))


##############################################################################
def stat_line(data, atr, box, vertices_yx, spacing_km, radius_km, gap_km=1.0, mask=None,
              keep_zero=False, smooth_km=None, mirrored=None):
    """Everything `viewinsar stat` draws and writes, for one line.

    Returns: dict - frame, dist_km (n,), points (n, 2), normals (n, 2),
             centres (2, n, 2), offset_km (2, n), and per block b and point k:
             values[b][k], rows[b][k], cols[b][k] (1D arrays), stats[b][k]
             (dict, see summarize), clipped[b][k] (bool); plus the inputs.
             smooth_km: the half length of line the normals follow (default:
             the radius; 0 = the local segment). mirrored: see polyline_frame.
    """
    if radius_km <= 0:
        raise ValueError('the radius must be positive')
    if gap_km < 0:
        raise ValueError('the gap must not be negative')
    frame = polyline_frame(atr, vertices_yx, mirrored=mirrored)
    smooth_km = radius_km if smooth_km is None else smooth_km
    dist, points, normals = sample_points(frame, spacing_km, smooth_km=smooth_km)
    centres, offsets = place_circles(frame, points, normals, radius_km, gap_km)
    out = dict(frame=frame, dist_km=dist, points=points, normals=normals, centres=centres,
               offset_km=offsets, length_km=frame['cum'][-1] * frame['scale'],
               spacing_km=spacing_km, radius_km=radius_km, gap_km=gap_km, smooth_km=smooth_km,
               values=[[], []], rows=[[], []], cols=[[], []], stats=[[], []],
               clipped=[[], []])
    for b in range(2):
        for k in range(len(dist)):
            val, rows, cols, clipped = circle_pixels(data, box, frame, centres[b, k], radius_km,
                                                     mask=mask, keep_zero=keep_zero)
            out['values'][b].append(val)
            out['rows'][b].append(rows)
            out['cols'][b].append(cols)
            out['stats'][b].append(summarize(val))
            out['clipped'][b].append(clipped)
    return out

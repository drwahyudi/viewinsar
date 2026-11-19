"""Fault creep along a line - the core of `viewinsar stat`.

Draw a line along a fault (as many vertices as needed). Every SPACING km along
it a short profile crosses it, as Hussain et al. (2016) measured the creep of
the North Anatolian Fault: LENGTH km (default 5) either side of the line,
the pixels within WIDTH km (default 2.5) either side of the profile projected
onto it. Block 1 is the LEFT side of the line, walking from its first vertex
to its last; block 2 the right. One straight line through the pixels of each
side, and the creep rate is their offset at the line (their eq. 2):

    v = slope_1 x + b                (block 1, x < 0)
    v = slope_2 x + b + creep        (block 2, x > 0)

x the distance from the line along the profile, + towards block 2 - so the
creep is block 2 minus block 1 at the fault. Solved by least squares (the
BLUE) and by a percentile bootstrap that resamples each side on its own, n1
of n1 and n2 of n2 (model.fit_two_lines); Hussain et al. use the bootstrap.

A profile is perpendicular to the line's direction over +/- WIDTH around its
point, so a hand-digitised trace with small wiggles still gives steady
profiles. A pixel on the other side of the line than its x says - where the
line bends inside the window - is left out.

Everything is computed in a metric frame of the data grid (X along the
columns, Y along -rows; MintPy's ground pixel size on a radar grid, a local
equirectangular frame on a lat/lon one), so the distances are true ones. A
grid that is a mirror image of the map - every right-looking radar grid, which
is why view.py flips it - gets X along -columns instead, so "left" and "right"
are always those of the ground (and of the map as view.py shows it).
Numpy only: no MintPy, no matplotlib.
"""

import numpy as np

from .model import fit_two_lines
from .swath import get_pixel_size, is_projected

__all__ = ['grid_mirrored', 'polyline_frame', 'from_yx', 'to_yx', 'position_along',
           'sample_points', 'distance_to_polyline', 'side_of_polyline', 'profile_axes',
           'profile_outline', 'trace_offset', 'profile_pixels', 'creep_profile',
           'creep_profiles', 'creep_along_line']

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


##############################################################################
#                       the short profiles across the line                   #
##############################################################################
def profile_axes(normal):
    """(across, along): unit vectors of a profile at a point whose LEFT
    normal is `normal` - across the line towards block 2 (the right), and
    along it (its direction, first vertex -> last)."""
    n = np.asarray(normal, dtype=np.float64)
    return -n, np.array([n[1], -n[0]])


def profile_outline(frame, point, normal, length_km, width_km, x_from=None, x_to=None):
    """(rows, cols) of the closed outline of a profile window: from x_from to
    x_to km across the line (default -length .. +length, + towards block 2)
    and +/- width_km along it - for drawing on the map."""
    across, along = profile_axes(normal)
    x0 = -length_km if x_from is None else x_from
    x1 = length_km if x_to is None else x_to
    w = width_km
    corners = np.array([[x0, -w], [x1, -w], [x1, w], [x0, w], [x0, -w]]) / frame['scale']
    xy = np.asarray(point, dtype=np.float64) + corners[:, :1] * across + corners[:, 1:] * along
    return to_yx(frame, xy)


def _dense_line(frame):
    """The line sampled every quarter pixel (cached in the frame)."""
    if 'dense' not in frame:
        step = 0.25 * min(frame['x_step'], frame['y_step'])
        num = int(np.ceil(frame['cum'][-1] / step)) + 1
        frame['dense'] = position_along(frame, np.linspace(0.0, frame['cum'][-1], num))
    return frame['dense']


def trace_offset(frame, point, normal, length_km, width_km):
    """How far the line strays from the centre of a profile window: the
    largest |x| (km) of the line within +/- width_km of the profile - 0 on a
    straight line, more where it bends, or where another stretch of it comes
    back into the window."""
    across, along = profile_axes(normal)
    rel = _dense_line(frame) - np.asarray(point, dtype=np.float64)
    x, y = (rel @ across) * frame['scale'], (rel @ along) * frame['scale']
    inside = (np.abs(y) <= width_km) & (np.abs(x) <= length_km)
    return float(np.abs(x[inside]).max()) if inside.any() else 0.0


def profile_pixels(data, box, frame, point, normal, length_km, width_km, mask=None,
                   keep_zero=False, same_side=True):
    """Valid pixels of a profile window across the line at `point`: within
    length_km of the line along the profile, either side, and within width_km
    of the profile.

    Parameters: data   - 2D array, read within `box`
                box    - (x0, y0, x1, y1) of that window, in full-file pixels
                frame  - from polyline_frame()
                point  - (2,) frame coordinates of the profile on the line
                normal - (2,) unit normal to the LEFT of the line there
                same_side - leave out a pixel on the other side of the line
                         than its x says (where the line bends)
    Returns:    dict - values, rows, cols (full file), x_km (along the
                profile from the line, + towards block 2), y_km (along the
                line), clipped (the window reaches past the data read),
                dropped (pixels left out by same_side)
    """
    empty = dict(values=np.array([]), rows=np.array([], dtype=np.int64),
                 cols=np.array([], dtype=np.int64), x_km=np.array([]), y_km=np.array([]),
                 clipped=False, dropped=0)
    point = np.asarray(point, dtype=np.float64)
    if not (np.all(np.isfinite(point)) and np.all(np.isfinite(normal))):
        return empty
    rows_c, cols_c = profile_outline(frame, point, normal, length_km, width_km)
    length, width = data.shape[-2:]
    r0, r1 = int(np.floor(rows_c.min())), int(np.ceil(rows_c.max()))
    c0, c1 = int(np.floor(cols_c.min())), int(np.ceil(cols_c.max()))
    rr0, rr1 = max(r0, box[1]), min(r1, box[1] + length - 1)
    cc0, cc1 = max(c0, box[0]), min(c1, box[0] + width - 1)
    empty['clipped'] = clipped = bool(rr0 > r0 or rr1 < r1 or cc0 > c0 or cc1 < c1)
    if rr1 < rr0 or cc1 < cc0:
        return empty

    rows, cols = np.arange(rr0, rr1 + 1), np.arange(cc0, cc1 + 1)
    across, along = profile_axes(normal)
    dx = frame['sx'] * cols * frame['x_step'] - point[0]
    dy = -rows * frame['y_step'] - point[1]
    x = (dx[None, :] * across[0] + dy[:, None] * across[1]) * frame['scale']
    y = (dx[None, :] * along[0] + dy[:, None] * along[1]) * frame['scale']
    eps = 1e-9 * max(length_km, width_km)
    inside = (np.abs(x) <= length_km + eps) & (np.abs(y) <= width_km + eps) & (x != 0)
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
    out = dict(values=val[ri, ci], rows=rows[ri], cols=cols[ci], x_km=x[ri, ci],
               y_km=y[ri, ci], clipped=clipped, dropped=0)
    if same_side and ri.size:
        # only the pixels nearer the line than it strays can be on its other side
        reach = trace_offset(frame, point, normal, length_km, width_km)
        pixel = max(frame['x_step'], frame['y_step']) * frame['scale']
        near = np.abs(out['x_km']) <= reach + 2.0 * pixel
        if near.any():
            xy = from_yx(frame, out['rows'][near], out['cols'][near])
            wrong = np.zeros(ri.size, dtype=bool)
            # block 1 (x < 0) is the left of the line (+1), block 2 the right (-1)
            wrong[near] = side_of_polyline(frame, xy) != np.where(out['x_km'][near] < 0, 1, -1)
            if wrong.any():
                for key in ('values', 'rows', 'cols', 'x_km', 'y_km'):
                    out[key] = out[key][~wrong]
                out['dropped'] = int(wrong.sum())
    return out


def creep_profile(data, box, frame, point, normal, length_km=5.0, width_km=2.5, mask=None,
                  keep_zero=False, factor=1.0, n_boot=2000, seed=0, min_points=10):
    """The creep rate at one point of the line, as Hussain et al. (2016):
    profile_pixels(), then a straight line through the pixels of each block
    and their offset at the line (model.fit_two_lines) - block 2 minus block 1.

    Returns: dict of profile_pixels() (the values multiplied by `factor`),
             plus n_1, n_2 (pixels of block 1, block 2), trace_km (see
             trace_offset), and fit - fit_two_lines() of block 1 (side a) and
             block 2 (side b), or None with fewer than min_points a side.
    """
    out = profile_pixels(data, box, frame, point, normal, length_km, width_km, mask=mask,
                         keep_zero=keep_zero)
    out['values'] = out['values'] * factor
    a, b = out['x_km'] < 0, out['x_km'] > 0
    out.update(n_1=int(a.sum()), n_2=int(b.sum()), fit=None,
               trace_km=(trace_offset(frame, point, normal, length_km, width_km)
                         if np.all(np.isfinite(point)) else np.nan))
    if min(out['n_1'], out['n_2']) >= max(int(min_points), 3):
        out['fit'] = fit_two_lines(out['x_km'][a], out['values'][a], out['x_km'][b],
                                   out['values'][b], n_boot=n_boot, seed=seed,
                                   min_points=max(int(min_points), 3))
    return out


def creep_profiles(data, box, frame, points, normals, length_km=5.0, width_km=2.5, mask=None,
                   keep_zero=False, factor=1.0, n_boot=2000, seed=0, min_points=10):
    """creep_profile() at every point; the bootstrap of point k is seeded
    with (seed, k) - reproducible, and independent of the others. `seed` may
    be a tuple (seed, file index) for several files."""
    seed = tuple(np.atleast_1d(seed).astype(int).tolist())
    return [creep_profile(data, box, frame, points[k], normals[k], length_km, width_km,
                          mask=mask, keep_zero=keep_zero, factor=factor, n_boot=n_boot,
                          seed=[*seed, k], min_points=min_points)
            for k in range(len(points))]


def creep_along_line(data, atr, box, vertices_yx, spacing_km, length_km=5.0, width_km=2.5,
                     mask=None, keep_zero=False, factor=1.0, n_boot=2000, seed=0,
                     min_points=10, smooth_km=None, mirrored=None):
    """Everything `viewinsar stat` draws and writes, for one file.

    Points every spacing_km along the line (sample_points: centred on it), a
    profile across the line at each (creep_profile). smooth_km: the half
    length of line the profiles are perpendicular to (default: width_km; 0 =
    the local segment); mirrored: see polyline_frame.

    Returns: dict - frame, dist_km (n,), points (n, 2), normals (n, 2),
             profiles (a creep_profile() each), length_km (of the line), and
             the parameters.
    """
    if spacing_km <= 0:
        raise ValueError('the spacing must be positive')
    if length_km <= 0 or width_km <= 0:
        raise ValueError('the profiles need a positive length and width')
    frame = polyline_frame(atr, vertices_yx, mirrored=mirrored)
    smooth_km = width_km if smooth_km is None else smooth_km
    dist, points, normals = sample_points(frame, spacing_km, smooth_km=smooth_km)
    profiles = creep_profiles(data, box, frame, points, normals, length_km, width_km,
                              mask=mask, keep_zero=keep_zero, factor=factor, n_boot=n_boot,
                              seed=seed, min_points=min_points)
    return dict(frame=frame, dist_km=dist, points=points, normals=normals, profiles=profiles,
                length_km=frame['cum'][-1] * frame['scale'], spacing_km=spacing_km,
                profile_km=length_km, width_km=width_km, smooth_km=smooth_km, n_boot=n_boot,
                min_points=min_points)

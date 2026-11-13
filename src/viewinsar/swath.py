"""Buffer a line, linearly reference the pixels inside it, and bin them.

This module is the numerical core of viewinsar and deliberately depends on
numpy only, so it can be imported, tested and reused without MintPy or
matplotlib. The GIS equivalent of what it does is:

    1. buffer(line, half_width)                    - flat / square caps
    2. clip the raster with that polygon
    3. "Locate Features Along Routes" (ArcGIS) / LRS (QGIS) of the clipped
       pixels, i.e. project each one onto the line to get an along-profile
       "measure" and an off-profile distance
    4. bin the resulting (distance, value) cloud into a 2-D density grid
"""

import numpy as np

EARTH_RADIUS = 6.3781e6      # in meter

__all__ = [
    'EARTH_RADIUS', 'is_projected', 'get_pixel_size', 'center_incidence', 'line_geometry',
    'swath_extract', 'bilinear_sample', 'weighted_quantile', 'nice_step',
    'swath_density',
]


##############################################################################
#                              geometry                                      #
##############################################################################
def is_projected(atr):
    """Is this grid in a projected CRS (metres) rather than degrees?

    Deciding on UTM_ZONE alone is not enough: a GeoTIFF in a projected CRS
    usually carries EPSG / Y_UNIT but no UTM_ZONE, and treating its metre-sized
    pixels as degrees inflates every distance by a factor of ~111,320.
    """
    for key in ('Y_UNIT', 'X_UNIT'):
        unit = str(atr.get(key, '')).lower()
        if unit.startswith('meter') or unit.startswith('metre'):
            return True
        if unit.startswith('deg'):
            return False
    if atr.get('UTM_ZONE'):
        return True
    # no raster has a pixel a degree across - checked BEFORE the EPSG code,
    # because MintPy's GDAL reader can pick up the GEOGCS code (4326) of a
    # projected GeoTIFF whose WKT carries no code of its own
    try:
        if abs(float(atr['X_STEP'])) > 1.0:
            return True
    except (KeyError, TypeError, ValueError):
        pass
    code = atr.get('EPSG')
    if code not in (None, '', 'None'):
        try:
            from pyproj import CRS
            return bool(CRS.from_epsg(int(code)).is_projected)
        except Exception:        # noqa: BLE001
            pass
    return False


def get_pixel_size(atr, lat_c=None):
    """Ground pixel size of a MintPy file, in meters.

    Parameters: atr    - dict, file metadata
                lat_c  - float, latitude of the area of interest, used for the
                         degree -> meter conversion of a lat/lon geocoded file
    Returns:    x_step - float, ground size of one column, in meter (or pixel)
                y_step - float, ground size of one row,    in meter (or pixel)
                unit   - str, 'm' or 'pixel'
    """
    if 'Y_FIRST' in atr:
        y_step, x_step = abs(float(atr['Y_STEP'])), abs(float(atr['X_STEP']))
        if not is_projected(atr):
            # geographic degrees -> local equirectangular meters, the same
            # approximation MintPy uses in utils.utils.transect_yx()
            if lat_c is None:
                lat_c = float(atr['Y_FIRST']) + float(atr['Y_STEP']) * int(atr['LENGTH']) / 2.
            y_step *= np.pi / 180. * EARTH_RADIUS
            x_step *= np.pi / 180. * EARTH_RADIUS * np.cos(np.deg2rad(lat_c))
        return x_step, y_step, 'm'

    # radar coordinates - the formulas of MintPy's range / azimuth_ground_resolution
    try:
        inc = np.deg2rad(center_incidence(atr))
        x_step = float(atr['RANGE_PIXEL_SIZE']) / np.sin(inc)
        y_step = float(atr['AZIMUTH_PIXEL_SIZE'])
        if atr.get('PROCESSOR', 'isce') in ('roipac', 'isce'):
            radius = float(atr.get('EARTH_RADIUS', EARTH_RADIUS))
            y_step *= radius / (radius + float(atr['HEIGHT']))
        return x_step, y_step, 'm'
    except (KeyError, ValueError, ZeroDivisionError):
        return 1., 1., 'pixel'


def center_incidence(atr):
    """The incidence angle at mid-range, degrees: from the metadata, or from the
    slant range and the orbit height as MintPy's utils0.incidence_angle()."""
    for key in ('CENTER_INCIDENCE_ANGLE', 'INCIDENCE_ANGLE'):
        if key in atr:
            return float(atr[key])
    near = float(atr['STARTING_RANGE'])
    far = near + float(atr['RANGE_PIXEL_SIZE']) * int(atr['WIDTH'])
    radius = float(atr.get('EARTH_RADIUS', EARTH_RADIUS))
    height = float(atr['HEIGHT'])

    def inc(rng):
        return 180.0 - np.degrees(np.arccos((radius ** 2 + rng ** 2 - (radius + height) ** 2)
                                            / (2 * radius * rng)))
    return float((inc(near) + inc(far)) / 2.0)


def line_geometry(atr, start_yx, end_yx, half_width):
    """Set up the metric frame of the swath corridor.

    A local Cartesian frame is used with +X to the east (increasing column) and
    +Y to the north (decreasing row), so the sign of the off-profile distance
    follows the usual convention: positive = left of the start -> end direction.

    Parameters: atr        - dict, file metadata
                start_yx   - list of 2 int, [row, col] of the start point
                end_yx     - list of 2 int, [row, col] of the end   point
                half_width - float, half width of the corridor, in km (or pixel)
    Returns:    geom       - dict describing the corridor
    """
    if 'Y_FIRST' in atr and not is_projected(atr):
        lat_c = float(atr['Y_FIRST']) + float(atr['Y_STEP']) * (start_yx[0] + end_yx[0]) / 2.
    else:
        lat_c = None
    x_step, y_step, unit = get_pixel_size(atr, lat_c=lat_c)

    p0 = np.array([start_yx[1] * x_step, -start_yx[0] * y_step], dtype=np.float64)
    p1 = np.array([end_yx[1] * x_step, -end_yx[0] * y_step], dtype=np.float64)
    vec = p1 - p0
    length = float(np.hypot(vec[0], vec[1]))
    if length == 0:
        raise ValueError('the start and end point of the profile are identical!')

    uvec = vec / length                        # along-profile unit vector
    nvec = np.array([-uvec[1], uvec[0]])       # off-profile unit vector, to the left
    hwid = half_width * 1000. if unit == 'm' else half_width
    corners = np.array([p0 + hwid * nvec, p1 + hwid * nvec,
                        p1 - hwid * nvec, p0 - hwid * nvec])

    return dict(
        x_step=x_step, y_step=y_step, unit=unit,
        p0=p0, p1=p1, uvec=uvec, nvec=nvec,
        length=length, half_width=hwid, corners=corners,
        # the corners back in (row, col), to draw the corridor on the map
        corners_yx=np.stack([-corners[:, 1] / y_step, corners[:, 0] / x_step], axis=1),
        dist_scale=1e-3 if unit == 'm' else 1.,
        dist_unit='km' if unit == 'm' else 'pixel',
    )


def project_points(geom, rows, cols):
    """Along- and off-profile distance of points at fractional (row, col).

    The same frame as the samples of swath_extract (see line_geometry), so a
    GNSS station or a fault crossing lands exactly where a pixel would.
    Returns: along, perp - 1D arrays in km (or pixel); perp + = left
    """
    dx = np.asarray(cols, dtype=np.float64) * geom['x_step'] - geom['p0'][0]
    dy = -np.asarray(rows, dtype=np.float64) * geom['y_step'] - geom['p0'][1]
    along = dx * geom['uvec'][0] + dy * geom['uvec'][1]
    perp = dx * geom['nvec'][0] + dy * geom['nvec'][1]
    return along * geom['dist_scale'], perp * geom['dist_scale']


def bilinear_sample(data, rows, cols, mask=None):
    """Sample `data` at fractional (row, col), NaN-aware.

    A plain bilinear interpolation returns NaN as soon as any of the four
    neighbours is NaN, which would erode the valid area of a speckled InSAR
    product at every void edge. Instead the weights are renormalised over the
    neighbours that ARE finite, so a sample survives as long as one of its four
    neighbours is. Only where all four are invalid does it return NaN.

    Parameters: data  - 2D np.ndarray
                rows  - 1D np.ndarray of float row coordinates
                cols  - 1D np.ndarray of float column coordinates
                mask  - 2D np.ndarray of bool, optional; False counts as invalid
    Returns:    2D np.ndarray of shape (rows.size, cols.size)
    """
    nrow, ncol = data.shape[-2:]
    r0 = np.clip(np.floor(rows).astype(np.intp), 0, nrow - 1)
    c0 = np.clip(np.floor(cols).astype(np.intp), 0, ncol - 1)
    r1 = np.minimum(r0 + 1, nrow - 1)
    c1 = np.minimum(c0 + 1, ncol - 1)
    fr = (rows - r0)[:, None]
    fc = (cols - c0)[None, :]

    val = np.ma.getdata(data).astype(np.float64)
    good = np.isfinite(val)
    if np.ma.isMaskedArray(data):
        good &= ~np.ma.getmaskarray(data)
    if mask is not None:
        good &= mask.astype(bool)

    num = np.zeros((rows.size, cols.size), dtype=np.float64)
    den = np.zeros_like(num)
    for ri, wr in ((r0, 1.0 - fr), (r1, fr)):
        for ci, wc in ((c0, 1.0 - fc), (c1, fc)):
            wgt = wr * wc * good[np.ix_(ri, ci)]
            num += wgt * np.nan_to_num(val[np.ix_(ri, ci)])
            den += wgt
    with np.errstate(invalid='ignore', divide='ignore'):
        out = np.where(den > 0, num / den, np.nan)
    return out


def swath_extract(data, atr, box, start_yx, end_yx, half_width,
                  mask=None, keep_zero=False, decimate=1, nlook=1,
                  max_samples=20_000_000):
    """Buffer a line and linearly reference every valid pixel inside the buffer.

    Parameters: data       - 2D np.ndarray, the data read within `box`
                atr        - dict, file metadata
                box        - tuple of 4 int, (x0, y0, x1, y1) of `data` in the full file
                start_yx   - list of 2 int, [row, col] of the start point, full file
                end_yx     - list of 2 int, [row, col] of the end   point, full file
                half_width - float, half width of the swath, in km
                mask       - 2D np.ndarray of bool, same shape as data, True = keep
                keep_zero  - bool, keep the pixels with an exactly-zero value
                decimate   - int, use every Nth pixel of the corridor only
                nlook      - int, OVERSAMPLE the corridor by this factor in each
                             direction, i.e. sample it on a grid nlook times finer
                             than the native pixel (nlook**2 more samples), with
                             NaN-aware bilinear interpolation. 1 = native pixels.
                             Use it to recover a dense profile from a multilooked
                             raster; it is the inverse of `decimate`.
                max_samples - int, refuse to build a sample grid larger than this
    Returns:    swath      - dict of 1D arrays 'dist' (km along the profile from
                             the start point), 'perp' (km from the line, + = left),
                             'value', 'row', 'col' and, for geocoded files,
                             'lat' / 'lon'; plus the scalars 'length',
                             'half_width', 'dist_unit', 'clipped', 'nlook' and 'geom'.
    """
    geom = line_geometry(atr, start_yx, end_yx, half_width)
    x_step, y_step = geom['x_step'], geom['y_step']
    p0, uvec, nvec = geom['p0'], geom['uvec'], geom['nvec']
    length, hwid = geom['length'], geom['half_width']

    # 1. bounding box of the corridor, in full-file pixel coordinates
    rows_c, cols_c = geom['corners_yx'][:, 0], geom['corners_yx'][:, 1]
    r0, r1 = int(np.floor(rows_c.min())) - 1, int(np.ceil(rows_c.max())) + 1
    c0, c1 = int(np.floor(cols_c.min())) - 1, int(np.ceil(cols_c.max())) + 1

    # 2. clip it to the window that was actually read
    len_d, wid_d = data.shape[-2:]
    rr0, rr1 = max(r0, box[1]), min(r1, box[1] + len_d - 1)
    cc0, cc1 = max(c0, box[0]), min(c1, box[0] + wid_d - 1)
    clipped = bool(rr0 > r0 or rr1 < r1 or cc0 > c0 or cc1 < c1)
    if rr1 < rr0 or cc1 < cc0:
        raise ValueError('the swath corridor falls entirely outside of the data coverage!')

    decimate = max(1, int(decimate))
    nlook = max(1, int(nlook))
    if nlook > 1 and decimate > 1:
        raise ValueError('--nlook and --decimate are opposites; use one or the other.')

    if nlook > 1:
        # sub-pixel grid: nlook samples per native pixel, per axis
        step = 1.0 / nlook
        n_r = int(round((rr1 - rr0) * nlook)) + 1
        n_c = int(round((cc1 - cc0) * nlook)) + 1
        if n_r * n_c > max_samples:
            raise ValueError(
                f'--nlook {nlook} would need {n_r * n_c:,} samples for this corridor '
                f'(limit {max_samples:,}). Use a smaller --nlook, a shorter line '
                f'or a narrower --width.')
        rows = rr0 + np.arange(n_r) * step
        cols = cc0 + np.arange(n_c) * step
    else:
        rows = np.arange(rr0, rr1 + 1, decimate).astype(np.float64)
        cols = np.arange(cc0, cc1 + 1, decimate).astype(np.float64)

    # 3. along- / off-profile distance of every sample of the sub-window
    #    (outer sums, so only 1D arrays are broadcast into 2D)
    dx = cols * x_step - p0[0]
    dy = -rows * y_step - p0[1]
    d_along = dy[:, None] * uvec[1] + dx[None, :] * uvec[0]
    d_perp = dy[:, None] * nvec[1] + dx[None, :] * nvec[0]

    # 4. the data values at those samples
    if nlook > 1:
        val = bilinear_sample(data, rows - box[1], cols - box[0], mask=mask)
        flag = np.isfinite(val)
    else:
        idx = np.ix_(rows.astype(np.intp) - box[1], cols.astype(np.intp) - box[0])
        sub = data[idx]
        val = np.ma.getdata(sub).astype(np.float64)
        flag = np.isfinite(val)
        if np.ma.isMaskedArray(sub):
            flag &= ~np.ma.getmaskarray(sub)
        if mask is not None:
            flag &= mask[idx].astype(bool)

    # 5. the flat-capped buffer
    flag &= (np.abs(d_perp) <= hwid) & (d_along >= 0) & (d_along <= length)
    if not keep_zero:
        flag &= val != 0.

    ridx, cidx = np.where(flag)
    scale = geom['dist_scale']
    swath = dict(
        dist=(d_along[ridx, cidx] * scale).astype(np.float32),
        perp=(d_perp[ridx, cidx] * scale).astype(np.float32),
        value=val[ridx, cidx].astype(np.float32),
        row=rows[ridx].astype(np.float32 if nlook > 1 else np.int32),
        col=cols[cidx].astype(np.float32 if nlook > 1 else np.int32),
        length=length * scale, half_width=hwid * scale,
        dist_unit=geom['dist_unit'], clipped=clipped, nlook=nlook, geom=geom,
    )
    if 'Y_FIRST' in atr:
        # vectorized on purpose: coordinate.yx2lalo() loops in python
        swath['lat'] = float(atr['Y_FIRST']) + (swath['row'] + 0.5) * float(atr['Y_STEP'])
        swath['lon'] = float(atr['X_FIRST']) + (swath['col'] + 0.5) * float(atr['X_STEP'])
    return swath


##############################################################################
#                       density grid + statistics                            #
##############################################################################
def weighted_quantile(values, weights, q):
    """Weighted quantile of `values`, as used on one column of a density grid."""
    flag = weights > 0
    if not np.any(flag):
        return np.nan
    vals, wgts = np.asarray(values)[flag], np.asarray(weights)[flag]
    order = np.argsort(vals)
    vals, wgts = vals[order], wgts[order]
    cum = np.cumsum(wgts)
    idx = int(np.searchsorted(cum, q * cum[-1]))
    return float(vals[min(idx, vals.size - 1)])


def nice_step(vrange, num_bin=54):
    """A round-ish bin size that splits `vrange` into about `num_bin` bins."""
    raw = float(vrange) / max(num_bin, 1)
    if not np.isfinite(raw) or raw <= 0:
        return 1.
    mag = 10 ** np.floor(np.log10(raw))
    for mult in (1., 1.5, 2., 2.5, 3., 5., 7.5):
        if raw <= mult * mag:
            return mult * mag
    return 10. * mag


def swath_density(dist, value, bin_x, bin_y=None, ylim=None, x_offset=0.):
    """2-D point-density histogram of a swath profile, with per-bin statistics.

    The X bin edges are built outward from x = 0 in both directions, so x = 0 is
    ALWAYS a cell boundary - which matters once the profile is referenced to a
    fault sitting at x = 0.

    Parameters: dist     - 1D np.ndarray, along-profile distance, in km
                value    - 1D np.ndarray, the data value
                bin_x    - float, X cell size, in km
                bin_y    - float, Y cell size, in the data unit (auto if None)
                ylim     - list of 2 float, Y range of the grid (auto if None)
                x_offset - float, shift applied to `dist`, in km
    Returns:    den      - dict with the grid, the counts and the statistics.
                           `median_dw` / `q25_dw` / `q75_dw` are the
                           density-weighted quantiles that get drawn (so they
                           are clipped to the Y range), while `median`, `mean`,
                           `std`, `q25`, `q75` are exact and un-clipped.
    """
    x = np.asarray(dist, dtype=np.float64) - float(x_offset)
    y = np.asarray(value, dtype=np.float64)
    if x.size == 0:
        raise ValueError('the swath holds no valid pixel to bin.')

    if ylim is None:
        lo, hi = np.nanpercentile(y, [0.5, 99.5])
        pad = 0.05 * (hi - lo) if hi > lo else 1.
        ylim = [lo - pad, hi + pad]
    else:
        ylim = [float(ylim[0]), float(ylim[1])]
    if bin_y is None:
        bin_y = nice_step(ylim[1] - ylim[0])

    x_min, x_max = float(np.nanmin(x)), float(np.nanmax(x))
    left = np.arange(0, np.floor(x_min) - bin_x, -bin_x)[::-1]
    right = np.arange(0, np.ceil(x_max) + bin_x, bin_x)
    x_edges = np.unique(np.concatenate([left, right]))
    y_edges = np.arange(ylim[0], ylim[1] + bin_y, bin_y)

    counts = np.histogram2d(x, y, bins=[x_edges, y_edges])[0]
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2.
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2.
    num_x = x_centers.size

    # density-weighted quantiles - what the figure shows
    med_dw = np.full(num_x, np.nan)
    q25_dw = np.full(num_x, np.nan)
    q75_dw = np.full(num_x, np.nan)
    for i in range(num_x):
        col = counts[i, :]
        if col.sum() > 0:
            med_dw[i] = weighted_quantile(y_centers, col, 0.50)
            q25_dw[i] = weighted_quantile(y_centers, col, 0.25)
            q75_dw[i] = weighted_quantile(y_centers, col, 0.75)

    # exact, un-binned, un-clipped statistics - what goes into the .csv
    num = np.zeros(num_x, dtype=np.int64)
    med = np.full(num_x, np.nan)
    avg = np.full(num_x, np.nan)
    std = np.full(num_x, np.nan)
    q25 = np.full(num_x, np.nan)
    q75 = np.full(num_x, np.nan)
    ibin = np.digitize(x, x_edges) - 1
    # np.histogram puts a point sitting exactly on the rightmost edge in the last
    # bin, while np.digitize sends it one past the end - match the histogram, so
    # `count` and `counts.sum()` agree.
    ibin[x == x_edges[-1]] = num_x - 1
    inside = (ibin >= 0) & (ibin < num_x)
    if np.any(inside):
        ibin, yin = ibin[inside], y[inside]
        order = np.argsort(ibin, kind='stable')
        ibin, yin = ibin[order], yin[order]
        bounds = np.searchsorted(ibin, np.arange(num_x + 1))
        for i in range(num_x):
            seg = yin[bounds[i]:bounds[i + 1]]
            if seg.size:
                num[i] = seg.size
                med[i], avg[i], std[i] = np.median(seg), np.mean(seg), np.std(seg)
                q25[i], q75[i] = np.percentile(seg, [25, 75])

    return dict(
        x_edges=x_edges, y_edges=y_edges, counts=counts,
        x_centers=x_centers, y_centers=y_centers,
        bin_x=bin_x, bin_y=bin_y, ylim=ylim,
        median_dw=med_dw, q25_dw=q25_dw, q75_dw=q75_dw,
        count=num, median=med, mean=avg, std=std, q25=q25, q75=q75,
        num_outside=int(np.sum((y < y_edges[0]) | (y > y_edges[-1]))),
        x_offset=float(x_offset),
    )

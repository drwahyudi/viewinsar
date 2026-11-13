"""Reading rasters in whatever format, and moving a line between two grids.

`--topoxsect` has to take a DEM in any format and on any grid - a MintPy
geometry `.h5`, a GeoTIFF, an SRTM `.dem` with a `.rsc` - and it generally will
not share the projection of the velocity file (UTM data against a lat/lon DEM
is the common case). These helpers hide that.
"""

import os
import re

import numpy as np

__all__ = ['read_attr_any', 'read_any', 'crs_of', 'epsg_of', 'utm_epsg', 'same_crs',
           'crs_name', 'line_to_other_grid', 'lonlat_line_to_grid']


##############################################################################
def _rasterio_attrs(src):
    tfm = src.transform
    atr = {
        'LENGTH': str(src.height), 'WIDTH': str(src.width),
        'X_FIRST': str(tfm.c), 'Y_FIRST': str(tfm.f),
        'X_STEP': str(tfm.a), 'Y_STEP': str(tfm.e),
        'FILE_TYPE': os.path.splitext(src.name)[1].replace('.', '') or 'raster',
        'UNIT': '1', 'DATA_TYPE': str(src.dtypes[0]),
    }
    crs = src.crs
    if crs is not None:
        code = crs.to_epsg()
        if code:
            atr['EPSG'] = str(code)
        if crs.is_projected:
            atr['Y_UNIT'] = atr['X_UNIT'] = 'meters'
    return atr


def read_attr_any(fname):
    """Metadata of a raster - MintPy's reader (or its stand-in), rasterio failing that."""
    from .compat import read_attribute
    try:
        return dict(read_attribute(fname))
    except Exception:            # noqa: BLE001 - any reader failure falls through
        pass
    try:
        import rasterio
    except ImportError as err:
        raise ValueError(
            f'cannot read {fname}: it is not a MintPy / ROI_PAC file and rasterio is '
            'not installed. `pip install rasterio` adds GeoTIFF / GDAL support.') from err
    with rasterio.open(fname) as src:
        return _rasterio_attrs(src)


def read_any(fname, box=None, dset=None):
    """Read a 2D raster window, via MintPy first and rasterio as a fallback.

    Parameters: fname - str, path to the raster
                box   - tuple of 4 int, (x0, y0, x1, y1) window to read
                dset  - str, dataset name, for multi-dataset HDF5
    Returns:    data  - 2D np.ndarray of float32, no-data as NaN
                atr   - dict, metadata
    """
    from .compat import read, read_attribute
    try:
        if dset is None and str(fname).lower().endswith(('.h5', '.he5')):
            # a MintPy geometry file holds azimuthAngle, height, incidenceAngle ...:
            # its DEM is `height` (as in MintPy's own pp.read_dem), never the
            # first dataset, which would be the azimuth angle
            if read_attribute(fname).get('FILE_TYPE') == 'geometry':
                dset = 'height'
        data, atr = read(fname, datasetName=dset, box=box, print_msg=False)
        return np.array(data, dtype=np.float32), dict(atr)
    except Exception:            # noqa: BLE001
        pass

    import rasterio
    from rasterio.windows import Window
    with rasterio.open(fname) as src:
        atr = _rasterio_attrs(src)
        if box is None:
            data = src.read(1)
        else:
            x0, y0, x1, y1 = box
            data = src.read(1, window=Window(x0, y0, x1 - x0, y1 - y0))
        data = np.array(data, dtype=np.float32)
        if src.nodata is not None and np.isfinite(src.nodata):
            data[data == np.float32(src.nodata)] = np.nan
    return data, atr


##############################################################################
#                     coordinate reference systems                           #
##############################################################################
_UTM_RE = re.compile(r'(?<!\d)(\d{1,2})\s*([NS])(?![A-Za-z])', re.IGNORECASE)


def utm_epsg(zone):
    """EPSG code of a WGS84 UTM zone, however it is written.

    '51S', '51 S', 'UTM zone 51S' (GDAL) and 'WGS_1984_UTM_Zone_51S' (ESRI) all
    give 32751. The letter is the hemisphere, as in MintPy. None if unreadable.
    """
    match = _UTM_RE.search(str(zone or ''))
    if not match:
        return None
    number = int(match.group(1))
    if not 1 <= number <= 60:
        return None
    return (32700 if match.group(2).upper() == 'S' else 32600) + number


def _metric_grid(atr):
    """Do the grid's own numbers say metres? (units, a UTM zone, or a big step)"""
    for key in ('X_UNIT', 'Y_UNIT'):
        unit = str(atr.get(key, '')).lower()
        if unit.startswith(('meter', 'metre')):
            return True
        if unit.startswith('deg'):
            return False
    if utm_epsg(atr.get('UTM_ZONE')):
        return True
    try:
        return abs(float(atr['X_STEP'])) > 1.0       # no pixel is a degree across
    except (KeyError, TypeError, ValueError):
        return False


def _is_geographic(code):
    try:
        from pyproj import CRS
        return bool(CRS.from_epsg(int(code)).is_geographic)
    except Exception:            # noqa: BLE001
        return False


def _file_crs(fname):
    """The CRS written in a raster file itself (GeoTIFF & co.), as WKT, or None."""
    if not fname or str(fname).lower().endswith(('.h5', '.he5', '.hdf5')):
        return None              # MintPy HDF5 carries attributes, not a WKT
    try:
        from osgeo import gdal
        dset = gdal.Open(os.fspath(fname))
        wkt = dset.GetProjection() if dset is not None else ''
        if wkt:
            return wkt
    except Exception:            # noqa: BLE001 - GDAL missing or cannot open it
        pass
    try:
        import rasterio
        with rasterio.open(fname) as src:
            return src.crs.to_wkt() if src.crs else None
    except Exception:            # noqa: BLE001
        return None


def crs_of(atr, fname=None):
    """CRS of a geocoded file, for pyproj: an EPSG int, or WKT when it has no code.

    Sources, most direct first:
      1. the CRS written in the file itself (GeoTIFF, any GDAL raster)
      2. an EPSG attribute - unless it is geographic on a grid in metres:
         MintPy's GDAL reader can pick up the GEOGCS code 4326 of a projected
         GeoTIFF whose WKT carries no code of its own
      3. a UTM zone attribute, however it is spelled
      4. a grid in degrees -> WGS84 (EPSG:4326)
    Returns None for radar coordinates, or when nothing settles it.
    """
    if 'Y_FIRST' not in atr:
        return None
    wkt = _file_crs(fname) or atr.get('CRS_WKT')     # CRS_WKT: written by terrain.py
    if wkt:
        try:
            from pyproj import CRS
            crs = CRS.from_user_input(wkt)
            return crs.to_epsg(min_confidence=70) or crs.to_wkt()
        except Exception:        # noqa: BLE001 - fall back to the attributes
            pass
    try:
        code = int(str(atr.get('EPSG')).strip())
    except (TypeError, ValueError):
        code = None
    if code and code > 0 and not (_is_geographic(code) and _metric_grid(atr)):
        return code
    code = utm_epsg(atr.get('UTM_ZONE'))
    if code:
        return code
    return None if _metric_grid(atr) else 4326


def epsg_of(atr):
    """EPSG code of a geocoded file, from its attributes; None if it has none."""
    crs = crs_of(atr)
    return crs if isinstance(crs, int) else None


def same_crs(a, b):
    """Do two pyproj-ready CRS definitions describe the same system?"""
    if a == b:
        return True
    try:
        from pyproj import CRS
        return CRS.from_user_input(a).equals(CRS.from_user_input(b), ignore_axis_order=True)
    except Exception:            # noqa: BLE001
        return False


def crs_name(crs):
    """Short, printable name of a CRS: 'EPSG:32751' or its WKT name."""
    if crs is None:
        return 'unknown CRS'
    if isinstance(crs, int):
        return f'EPSG:{crs}'
    try:
        from pyproj import CRS
        return CRS.from_user_input(crs).name
    except Exception:            # noqa: BLE001
        return 'custom CRS'


def line_to_other_grid(atr_src, atr_dst, start_yx, end_yx, verbose=True,
                       crs_src=None, crs_dst=None, check_inside=True):
    """Map the two end points of a line from one geocoded grid onto another.

    Re-projects through pyproj when the two files use different coordinate
    reference systems, which is the usual case for UTM data against a lat/lon
    DEM. Pass crs_src / crs_dst (from crs_of with the file names) when known:
    the attributes alone can be incomplete.

    check_inside - bool, raise when an end point falls off the destination
                   grid. False lets a DEM that covers only part of the line
                   serve that part; the pixel numbers may then be negative or
                   past the grid.

    Returns: (start_yx, end_yx) as [row, col] on the destination grid.
    """
    for name, atr in (('data', atr_src), ('DEM', atr_dst)):
        if 'Y_FIRST' not in atr:
            raise ValueError(f'the {name} file is in radar coordinates; '
                             'a topography cross-section needs both to be geocoded.')

    x0s, y0s = float(atr_src['X_FIRST']), float(atr_src['Y_FIRST'])
    dxs, dys = float(atr_src['X_STEP']), float(atr_src['Y_STEP'])
    xs = [x0s + (c + 0.5) * dxs for c in (start_yx[1], end_yx[1])]
    ys = [y0s + (r + 0.5) * dys for r in (start_yx[0], end_yx[0])]

    c_src = crs_src if crs_src is not None else crs_of(atr_src)
    c_dst = crs_dst if crs_dst is not None else crs_of(atr_dst)
    if c_src is not None and c_dst is not None and not same_crs(c_src, c_dst):
        from pyproj import Transformer
        tfm = Transformer.from_crs(c_src, c_dst, always_xy=True)
        xs, ys = tfm.transform(xs, ys)
        if verbose:
            print(f'   re-projected the profile {crs_name(c_src)} -> {crs_name(c_dst)} '
                  'for the DEM')
    elif verbose and (c_src is None or c_dst is None):
        print(f'   WARNING: could not determine both CRSs (data: {crs_name(c_src)}, '
              f'DEM: {crs_name(c_dst)}); assuming they match.')

    x0d, y0d = float(atr_dst['X_FIRST']), float(atr_dst['Y_FIRST'])
    dxd, dyd = float(atr_dst['X_STEP']), float(atr_dst['Y_STEP'])
    cols = [int(round((x - x0d) / dxd - 0.5)) for x in xs]
    rows = [int(round((y - y0d) / dyd - 0.5)) for y in ys]

    length, width = int(atr_dst['LENGTH']), int(atr_dst['WIDTH'])
    inside = all(0 <= r < length for r in rows) and all(0 <= c < width for c in cols)
    if check_inside and not inside:
        raise ValueError(
            f'the profile falls outside the DEM: it maps to rows {rows}, cols {cols} '
            f'on a {length} x {width} grid ({crs_name(c_dst)}). Use a DEM that covers '
            'the profile.')
    return [rows[0], cols[0]], [rows[1], cols[1]]


def lonlat_line_to_grid(lonlat0, lonlat1, atr_dst, crs_dst):
    """Two lon/lat end points -> [row, col] on a geocoded grid (any CRS).

    Used for a radar-coded section against a geocoded DEM: the section's end
    points go to lon/lat through the lookup, then onto the DEM here. No bounds
    check: a DEM covering part of the line serves that part.
    """
    if 'Y_FIRST' not in atr_dst:
        raise ValueError('the DEM is in radar coordinates, but not on the data grid')
    xs = [float(lonlat0[0]), float(lonlat1[0])]
    ys = [float(lonlat0[1]), float(lonlat1[1])]
    if crs_dst is not None and not same_crs(4326, crs_dst):
        from pyproj import Transformer
        xs, ys = Transformer.from_crs(4326, crs_dst, always_xy=True).transform(xs, ys)
    x0d, y0d = float(atr_dst['X_FIRST']), float(atr_dst['Y_FIRST'])
    dxd, dyd = float(atr_dst['X_STEP']), float(atr_dst['Y_STEP'])
    cols = [int(round((x - x0d) / dxd - 0.5)) for x in xs]
    rows = [int(round((y - y0d) / dyd - 0.5)) for y in ys]
    return [rows[0], cols[0]], [rows[1], cols[1]]

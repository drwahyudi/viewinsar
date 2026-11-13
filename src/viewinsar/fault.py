"""Fault traces: read them from any vector format, bring them onto the data's
grid, and find where a cross-section crosses them.

Readers, picked by the file extension:
  .geojson .json        built in (RFC 7946 WGS84, or a legacy "crs" member)
  .shp                  pyshp, with the CRS from the .prj side-car file
  .gpkg                 built in (sqlite3), every feature table, its own CRS
  .kml .kmz             built in, WGS84 by definition
  .gmt .lonlat .xy .ll  built in: GMT multi-segment text, '>' between segments,
  .txt .dat             the CRS from an OGR-GMT "# @Je" / "# @Jp" header if any
  anything else         pyogrio -> fiona -> GDAL/OGR, whichever is installed
                        (FileGDB, GML, FlatGeobuf, MapInfo, DXF, ...)

Lines, multi-lines and polygon rings are all used as traces; points are
ignored. Everything is re-projected to the data's CRS with pyproj and cut to
the data extent. A file without any CRS information is taken as WGS84 lon/lat
when its numbers fit (|x| <= 360, |y| <= 90), and as the data's own CRS
otherwise (e.g. a UTM shapefile that lost its .prj).
"""

import json
import os
import re

import numpy as np

from .raster import crs_of, same_crs
from .swath import is_projected
from .wkb import as_geojson, wkb_to_geojson

__all__ = ['read_fault_file', 'load_faults', 'load_faults_radar', 'profile_crossings',
           'to_display']

GMT_EXT = ('.gmt', '.lonlat', '.xy', '.ll', '.txt', '.dat')


##############################################################################
#                               readers                                      #
#   every reader returns a list of layers: [(lines, crs), ...]               #
#   lines: list of (N, 2) float arrays of x/y; crs: EPSG int, WKT/PROJ str,  #
#   or None when the file does not say                                       #
##############################################################################
def _as_line(coords):
    arr = np.asarray(coords, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[0] < 2 or arr.shape[1] < 2:
        return None
    return arr[:, :2]


def _geojson_lines(geom, out):
    """Collect the line work of a GeoJSON object into `out`."""
    if not geom:
        return
    kind = geom.get('type')
    coords = geom.get('coordinates')
    if kind == 'LineString':
        parts = [coords]
    elif kind in ('MultiLineString', 'Polygon'):
        parts = coords
    elif kind == 'MultiPolygon':
        parts = [ring for poly in coords for ring in poly]
    elif kind == 'GeometryCollection':
        for sub in geom.get('geometries', []):
            _geojson_lines(sub, out)
        return
    elif kind == 'Feature':
        _geojson_lines(geom.get('geometry'), out)
        return
    elif kind == 'FeatureCollection':
        for feat in geom.get('features', []):
            _geojson_lines(feat.get('geometry'), out)
        return
    else:
        return                          # Point / MultiPoint: not a trace
    for part in parts or []:
        line = _as_line(part)
        if line is not None:
            out.append(line)


def _read_geojson(fname):
    with open(fname, encoding='utf-8') as fid:
        obj = json.load(fid)
    lines = []
    _geojson_lines(obj, lines)
    # RFC 7946: always WGS84 lon/lat. Files from before it may name a CRS.
    name = str(((obj.get('crs') or {}).get('properties') or {}).get('name', ''))
    crs = 4326
    if name and 'CRS84' not in name.upper():
        crs = name
        match = re.search(r'EPSG:+(\d+)', name, flags=re.IGNORECASE)
        if match:
            crs = int(match.group(1))
    return [(lines, crs)]


def _read_gmt(fname):
    """GMT multi-segment text: "x y" per line, segments split by '>'."""
    lines, seg, found = [], [], {}
    with open(fname, encoding='utf-8', errors='replace') as fid:
        for raw in fid:
            text = raw.strip()
            if not text:
                continue
            if text.startswith('#'):
                # OGR-GMT headers: # @Je4326   # @Jp"+proj=..."   # @Jw"WKT"
                match = re.search(r'@Je(\d+)', text)
                if match:
                    found.setdefault('e', int(match.group(1)))
                match = re.search(r'@Jp"([^"]+)"', text)
                if match:
                    found.setdefault('p', match.group(1))
                match = re.search(r'@Jw"((?:[^"\\]|\\.)+)"', text)
                if match:
                    found.setdefault('w', match.group(1).replace('\\"', '"'))
                continue
            if text.startswith('>'):
                if len(seg) >= 2:
                    lines.append(np.asarray(seg, dtype=np.float64))
                seg = []
                continue
            parts = text.replace(',', ' ').split()
            try:
                seg.append([float(parts[0]), float(parts[1])])
            except (IndexError, ValueError):
                continue                # a header / comment row
    if len(seg) >= 2:
        lines.append(np.asarray(seg, dtype=np.float64))
    return [(lines, found.get('e') or found.get('p') or found.get('w'))]


def _kml_lines(text):
    import xml.etree.ElementTree as ET

    lines = []
    root = ET.fromstring(text)
    for elem in root.iter():
        if elem.tag.rsplit('}', 1)[-1] not in ('LineString', 'LinearRing'):
            continue
        for sub in elem.iter():
            if sub.tag.rsplit('}', 1)[-1] != 'coordinates' or not sub.text:
                continue
            pts = []
            for tok in sub.text.split():
                parts = tok.split(',')
                try:
                    pts.append([float(parts[0]), float(parts[1])])
                except (IndexError, ValueError):
                    continue
            if len(pts) >= 2:
                lines.append(np.asarray(pts, dtype=np.float64))
    return lines


def _read_kml(fname):
    """KML / KMZ - always WGS84 lon/lat."""
    if fname.lower().endswith('.kmz'):
        import zipfile
        lines = []
        with zipfile.ZipFile(fname) as zfile:
            for name in zfile.namelist():
                if name.lower().endswith('.kml'):
                    lines += _kml_lines(zfile.read(name))
    else:
        with open(fname, 'rb') as fid:
            lines = _kml_lines(fid.read())
    return [(lines, 4326)]


def _read_pyshp(fname):
    import shapefile  # pyshp

    lines = []
    with shapefile.Reader(fname) as reader:
        for shp in reader.iterShapes():
            pts = np.asarray(shp.points, dtype=np.float64)
            if pts.ndim != 2 or len(pts) < 2:
                continue                       # a point, or an empty shape
            bounds = list(shp.parts) + [len(pts)]
            for i0, i1 in zip(bounds[:-1], bounds[1:]):
                if i1 - i0 >= 2:
                    lines.append(pts[i0:i1, :2])

    crs = None
    stem = os.path.splitext(fname)[0]
    for prj in (stem + '.prj', stem + '.PRJ'):
        if os.path.isfile(prj):
            with open(prj, encoding='utf-8', errors='replace') as fid:
                crs = fid.read().strip() or None
            break
    return [(lines, crs)]


def _gpkg_wkb(blob):
    """Strip the GeoPackage header off a geometry blob -> plain WKB."""
    blob = bytes(blob)
    if blob[:2] != b'GP':
        return blob                     # some writers store plain WKB
    flags = blob[3]
    if flags & 0x10:                    # the empty-geometry flag
        return None
    envelope = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}.get((flags >> 1) & 0x07)
    return None if envelope is None else blob[8 + envelope:]


def _read_gpkg(fname):
    """GeoPackage, read directly with sqlite3 - every feature table."""
    import pathlib
    import sqlite3

    uri = pathlib.Path(fname).resolve().as_uri() + '?mode=ro'
    try:
        con = sqlite3.connect(uri, uri=True)
    except sqlite3.Error:
        con = sqlite3.connect(fname)
    try:
        tables = con.execute('SELECT table_name, column_name, srs_id '
                             'FROM gpkg_geometry_columns').fetchall()
        layers = []
        for table, column, srs_id in tables:
            row = con.execute('SELECT organization, organization_coordsys_id, definition '
                              'FROM gpkg_spatial_ref_sys WHERE srs_id = ?',
                              (srs_id,)).fetchone()
            crs = None
            if row is not None:
                org, code, definition = row
                if org and str(org).upper() == 'EPSG' and code and int(code) > 0:
                    crs = int(code)
                elif definition and str(definition).strip().lower() not in ('', 'undefined'):
                    crs = str(definition)
            lines = []
            col, tab = column.replace('"', '""'), table.replace('"', '""')
            for (blob,) in con.execute(f'SELECT "{col}" FROM "{tab}"'):
                if blob is None:
                    continue
                wkb = _gpkg_wkb(blob)
                if wkb is None:
                    continue
                try:
                    _geojson_lines(wkb_to_geojson(wkb), lines)
                except Exception:        # noqa: BLE001 - skip an odd geometry
                    continue
            layers.append((lines, crs))
        return layers
    finally:
        con.close()


def _read_ogr(fname):
    """Any OGR-readable format, through whichever binding is installed."""
    try:
        import pyogrio
        from pyogrio import raw
    except ImportError:
        pass
    else:
        layers = []
        for name, gtype in pyogrio.list_layers(fname):
            if gtype is None:
                continue                # an attribute-only table
            meta, _, wkbs = raw.read(fname, layer=name)[:3]
            lines = []
            for wkb in wkbs:
                if wkb is not None:
                    _geojson_lines(wkb_to_geojson(wkb), lines)
            layers.append((lines, meta.get('crs')))
        return layers

    try:
        import fiona
    except ImportError:
        pass
    else:
        layers = []
        for name in fiona.listlayers(fname):
            lines = []
            with fiona.open(fname, layer=name) as src:
                crs = src.crs_wkt or None
                for feat in src:
                    if feat['geometry'] is not None:
                        _geojson_lines(as_geojson(feat['geometry']), lines)
            layers.append((lines, crs))
        return layers

    try:
        from osgeo import ogr
    except ImportError as err:
        raise ValueError(
            f'cannot read {os.path.basename(fname)}: .geojson, .shp, .gpkg, .kml/.kmz '
            'and GMT text are read directly; for any other vector format install '
            'pyogrio (`pip install pyogrio`), fiona or the GDAL python bindings.') from err
    dsrc = ogr.Open(fname)
    if dsrc is None:
        raise ValueError(f'GDAL/OGR could not open {fname}')
    layers = []
    for i in range(dsrc.GetLayerCount()):
        layer = dsrc.GetLayer(i)
        srs = layer.GetSpatialRef()
        lines = []
        for feat in layer:
            geo = feat.GetGeometryRef()
            if geo is not None:
                _geojson_lines(wkb_to_geojson(bytes(geo.ExportToWkb())), lines)
        layers.append((lines, srs.ExportToWkt() if srs is not None else None))
    return layers


def read_fault_file(fname):
    """Read the line work of a vector file.

    Returns: list of (lines, crs) - one per layer; lines is a list of (N, 2)
             float arrays of x/y in the file's own CRS, crs anything pyproj
             accepts (EPSG int, WKT, PROJ string) or None when the file does
             not say.
    """
    if not os.path.exists(fname):
        raise ValueError(f'no such fault file: {fname}')
    ext = os.path.splitext(fname)[1].lower()
    if ext in ('.geojson', '.json'):
        return _read_geojson(fname)
    if ext in ('.kml', '.kmz'):
        return _read_kml(fname)
    if ext in GMT_EXT:
        return _read_gmt(fname)
    if ext == '.shp':
        try:
            return _read_pyshp(fname)
        except ImportError:
            pass
    if ext == '.gpkg':
        try:
            return _read_gpkg(fname)
        except Exception:            # noqa: BLE001 - let OGR have a go
            pass
    return _read_ogr(fname)


##############################################################################
#                        onto the data's grid                                #
##############################################################################
def _split_finite(line):
    """Split a line wherever a vertex failed to project (inf / nan)."""
    good = np.all(np.isfinite(line), axis=1)
    if good.all():
        return [line]
    out, start = [], None
    for i, ok in enumerate(np.r_[good, False]):
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            if i - start >= 2:
                out.append(line[start:i])
            start = None
    return out


def _data_crs(atr, data_crs=None):
    if 'Y_FIRST' not in atr:
        raise ValueError('--fault needs a geocoded data file (this one is in radar '
                         'coordinates).')
    crs = data_crs if data_crs is not None else crs_of(atr)
    if crs is None:
        raise ValueError('could not work out the coordinate reference system of the '
                         'data file (no EPSG / UTM_ZONE attribute).')
    return crs


def _data_box(atr, margin):
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    x0, y0 = float(atr['X_FIRST']), float(atr['Y_FIRST'])
    x1, y1 = x0 + width * float(atr['X_STEP']), y0 + length * float(atr['Y_STEP'])
    xlo, xhi, ylo, yhi = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
    pad_x, pad_y = (xhi - xlo) * margin, (yhi - ylo) * margin
    return xlo - pad_x, xhi + pad_x, ylo - pad_y, yhi + pad_y


def _keep_near(lines, box, wrap):
    """The lines whose bounding box meets `box`; on a lat/lon grid also try +/-360."""
    xlo, xhi, ylo, yhi = box
    kept = []
    for line in lines:
        if line[:, 1].max() < ylo or line[:, 1].min() > yhi:
            continue
        lx0, lx1 = line[:, 0].min(), line[:, 0].max()
        for shift in ((0.0, 360.0, -360.0) if wrap else (0.0,)):
            if lx1 + shift >= xlo and lx0 + shift <= xhi:
                kept.append(line + np.array([shift, 0.0]) if shift else line)
                break
    return kept


def _traces_in(files, dst_crs, box, wrap, verbose, what='the data extent', target=None):
    """Every trace of `files`, in dst_crs, whose bounding box meets `box`.

    target - str, how to call dst_crs in the log (default: the data's CRS)
    """
    out = []
    for fname in files:
        n_in, kept, notes = 0, [], set()
        for lines, src_crs in read_fault_file(fname):
            lines = [ln for ln in lines if np.isfinite(ln).all(axis=1).sum() >= 2]
            n_in += len(lines)
            if not lines:
                continue
            if src_crs is None:
                xy = np.vstack(lines)
                if np.all(np.abs(xy[:, 0]) <= 360) and np.all(np.abs(xy[:, 1]) <= 90):
                    src_crs = 4326
                    notes.add('no CRS in the file, the numbers look like lon/lat')
                else:
                    src_crs = dst_crs
                    notes.add("no CRS in the file, taken as the data's own")

            if same_crs(src_crs, dst_crs):
                notes.add(f'already {target}' if target else 'same CRS as the data')
            else:
                from pyproj import Transformer
                tfm = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
                xy = np.vstack(lines)
                tx, ty = tfm.transform(xy[:, 0], xy[:, 1])
                pieces = np.split(np.c_[tx, ty], np.cumsum([len(ln) for ln in lines])[:-1])
                lines = [part for piece in pieces for part in _split_finite(piece)]
                notes.add(f're-projected to {target}' if target else 're-projected to the data CRS')
            kept += _keep_near(lines, box, wrap=wrap)

        out += kept
        if verbose:
            print(f'fault: {os.path.basename(fname)}: {len(kept)} of {n_in} trace(s) '
                  f'within {what} ({"; ".join(sorted(notes)) or "empty file"})')
    return out


def load_faults(files, atr, margin=0.25, verbose=True, data_crs=None):
    """Read fault files and return their traces in the DATA's map coordinates.

    Parameters: files    - list of str, vector files, any mix of formats and CRSs
                atr      - dict, metadata of the (geocoded) data file
                margin   - float, keep traces this fraction beyond the data extent
                data_crs - the data's CRS (EPSG int or WKT), e.g. from
                           raster.crs_of(atr, fname); worked out from atr if None
    Returns:    list of (N, 2) arrays of x/y: lon/lat on a lat/lon grid,
                easting/northing on a projected one.
    """
    data_crs = _data_crs(atr, data_crs)
    box = _data_box(atr, margin)
    return _traces_in(files, data_crs, box, wrap=not is_projected(atr), verbose=verbose)


def _densify(line, max_step):
    """Extra vertices so no segment is longer than max_step (in its own units)."""
    seg = np.hypot(*np.diff(line, axis=0).T)
    num = np.maximum(1, np.ceil(seg / max_step).astype(int))
    if num.max() == 1:
        return line
    pieces = [line[i] + np.linspace(0, 1, n, endpoint=False)[:, None] * (line[i + 1] - line[i])
              for i, n in enumerate(num)]
    return np.vstack(pieces + [line[-1:]])


def load_faults_radar(files, lookup, margin=0.1, verbose=True, max_step_deg=0.001):
    """Fault traces on a RADAR grid, through the MintPy lookup (geometryRadar.h5).

    The traces are taken to lon/lat, cut to the radar frame, densified (a line
    straight in lon/lat bends a little on the radar grid) and mapped to pixels.

    Parameters: files  - list of str, vector files, any format and CRS
                lookup - radar.RadarLookup of the data
    Returns:    list of (N, 2) arrays of (col, row) - fractional pixel numbers,
                the "map coordinates" of radar.PIXEL_ATR
    """
    west, south, east, north = lookup.extent(margin)
    lines = _traces_in(files, 4326, (west, east, south, north), wrap=True, verbose=verbose,
                       what='the radar frame', target='lon/lat')
    out = []
    for line in lines:
        dense = _densify(line, max_step_deg)
        rows, cols = lookup.to_rowcol(dense[:, 0], dense[:, 1])
        out.extend(_split_finite(np.c_[cols, rows]))
    return out


def to_display(faults, atr, fig_coord='geo'):
    """Fault traces in the coordinates of the map axes, as ONE NaN-separated line.

    view.py draws a geocoded file in its map coordinates ('geo') or, with
    --coord yx, in pixel numbers where pixel centres sit on integers.
    """
    if not faults:
        return np.array([]), np.array([])
    xy = np.vstack([np.vstack([ln, [np.nan, np.nan]]) for ln in faults])[:-1]
    xs, ys = xy[:, 0], xy[:, 1]
    if fig_coord != 'geo':
        xs = (xs - float(atr['X_FIRST'])) / float(atr['X_STEP']) - 0.5
        ys = (ys - float(atr['Y_FIRST'])) / float(atr['Y_STEP']) - 0.5
    return xs, ys


##############################################################################
#                       crossings with a profile                             #
##############################################################################
def _cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def profile_crossings(faults, atr, start_yx, end_yx, length_km, merge_km=0.02):
    """Where a cross-section's centre line crosses the fault traces.

    The profile's distance axis is an affine image of the map coordinates along
    the straight line - exactly so on a projected grid, and in the local
    equirectangular frame swath_extract() uses on a lat/lon grid - so a crossing
    at parameter t along the line lies t * length_km from the start.

    Parameters: faults    - list of (N, 2) arrays, from load_faults()
                atr       - dict, metadata of the data file
                start_yx  - [row, col] of the section start, full file
                end_yx    - [row, col] of the section end,   full file
                length_km - float, the section length, from swath_extract()
                merge_km  - float, crossings closer than this count once (a
                            fault vertex sitting exactly on the line, or two
                            mapped strands a few metres apart)
    Returns:    dists  - sorted 1D array, km from the section start
                points - (n, 2) array, the crossings in map coordinates
    """
    x0, y0 = float(atr['X_FIRST']), float(atr['Y_FIRST'])
    dx, dy = float(atr['X_STEP']), float(atr['Y_STEP'])
    pnt0 = np.array([x0 + (start_yx[1] + 0.5) * dx, y0 + (start_yx[0] + 0.5) * dy])
    pnt1 = np.array([x0 + (end_yx[1] + 0.5) * dx, y0 + (end_yx[0] + 0.5) * dy])
    vec = pnt1 - pnt0
    if not np.any(vec):
        return np.array([]), np.empty((0, 2))

    ts = []
    for line in faults:
        if len(line) < 2:
            continue
        seg_a, edge = line[:-1], np.diff(line, axis=0)
        denom = _cross(vec[None, :], edge)
        rel = seg_a - pnt0
        with np.errstate(divide='ignore', invalid='ignore'):
            t = _cross(rel, edge) / denom          # along the profile
            s = _cross(rel, vec[None, :]) / denom  # along the fault segment
        hit = (denom != 0) & (t >= 0) & (t <= 1) & (s >= 0) & (s <= 1)
        ts.extend(t[hit].tolist())

    if not ts:
        return np.array([]), np.empty((0, 2))
    ts = np.sort(np.asarray(ts, dtype=np.float64))
    dists = ts * float(length_km)
    keep = np.r_[True, np.diff(dists) > merge_km]
    ts, dists = ts[keep], dists[keep]
    points = pnt0[None, :] + ts[:, None] * vec[None, :]
    return dists, points

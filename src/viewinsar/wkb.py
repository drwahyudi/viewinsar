"""Well-known binary (WKB) geometries -> GeoJSON-like dicts, without shapely.

GeoPackage, pyogrio and GDAL/OGR hand geometries over as WKB: ISO (with the
Z / M / ZM codes 1000 / 2000 / 3000) or extended (PostGIS) WKB, either byte
order. Only x / y are kept.
"""

import struct

import numpy as np

__all__ = ['wkb_to_geojson', 'as_geojson']

_KINDS = {1: 'Point', 2: 'LineString', 3: 'Polygon', 4: 'MultiPoint',
          5: 'MultiLineString', 6: 'MultiPolygon', 7: 'GeometryCollection'}


def wkb_to_geojson(blob):
    """{'type': 'LineString', 'coordinates': [[x, y], ...]} ... from WKB bytes."""
    geom, _ = _read(memoryview(bytes(blob)), 0)
    return geom


def as_geojson(geom):
    """A GeoJSON-like dict from WKB, a __geo_interface__ object or a dict."""
    if geom is None:
        return None
    if isinstance(geom, (bytes, bytearray, memoryview)):
        return wkb_to_geojson(geom)
    geom = getattr(geom, '__geo_interface__', geom)
    return dict(geom)


def _read(buf, pos):
    order = '<' if buf[pos] == 1 else '>'
    pos += 1
    (code,) = struct.unpack_from(order + 'I', buf, pos)
    pos += 4
    has_z, has_m = bool(code & 0x80000000), bool(code & 0x40000000)
    if code & 0x20000000:                       # EWKB: an SRID follows
        pos += 4
    code &= 0x0FFFFFFF
    if code >= 1000:
        flag, code = divmod(code, 1000)
        has_z, has_m = has_z or flag in (1, 3), has_m or flag in (2, 3)
    ndim = 2 + int(has_z) + int(has_m)
    kind = _KINDS.get(code)
    if kind is None:
        raise ValueError(f'unsupported WKB geometry type {code}')

    def count(pos):
        return struct.unpack_from(order + 'I', buf, pos)[0], pos + 4

    def points(num, pos):
        arr = np.frombuffer(buf, dtype=np.dtype(order + 'f8'), count=num * ndim, offset=pos)
        return arr.reshape(num, ndim)[:, :2].tolist(), pos + 8 * num * ndim

    if kind == 'Point':
        pts, pos = points(1, pos)
        return {'type': kind, 'coordinates': pts[0]}, pos
    if kind == 'LineString':
        num, pos = count(pos)
        pts, pos = points(num, pos)
        return {'type': kind, 'coordinates': pts}, pos
    if kind == 'Polygon':
        rings, (num, pos) = [], count(pos)
        for _ in range(num):
            npt, pos = count(pos)
            pts, pos = points(npt, pos)
            rings.append(pts)
        return {'type': kind, 'coordinates': rings}, pos
    parts, (num, pos) = [], count(pos)
    for _ in range(num):
        part, pos = _read(buf, pos)
        parts.append(part)
    if kind == 'GeometryCollection':
        return {'type': kind, 'geometries': parts}, pos
    return {'type': kind, 'coordinates': [p['coordinates'] for p in parts]}, pos

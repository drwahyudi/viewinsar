"""GNSS velocities on InSAR profiles.

Read (almost) any table or point file of GNSS velocities, take the component
the InSAR shows - east, north, up, or the satellite line of sight - and place
the stations on a swath profile.

Columns are found by their names - lon / lat, ve / vn / vu, their sigmas and
the site - in any order and in the spellings common tables use ("Ve(mm/yr)",
"vel_e", "sig_e", "dE/dt" ...). A file without a header row is read in the
order of GMT psvelo: lon lat ve vn se sn [corr] [site]. Read directly: .csv,
.tsv, .txt and other text tables, .xlsx (first sheet), .shp, .geojson, .gpkg,
.kml / .kmz; any other vector format through pyogrio / fiona / GDAL.

    from viewinsar.gnss import read_gnss, los_unit_vector
    gnss = read_gnss('velocities.csv')                 # mm/yr
    v, s = gnss.component('los', los_unit_vector(39.0, -102.0))

The LOS follows MintPy: the unit vector from the ground to the satellite,
from the incidence angle and the azimuth angle of that vector (measured from
north, anti-clockwise positive) - positive LOS = motion towards the satellite.
"""

import json
import os
import re

import numpy as np

__all__ = ['GNSS', 'read_gnss', 'los_unit_vector', 'heading_to_azimuth', 'LOSGeometry',
           'COMPONENTS']

COMPONENTS = ('ve', 'vn', 'vu', 'los')
TO_MM = {'mm': 1.0, 'cm': 10.0, 'm': 1000.0}

# normalised column names (lower case, letters and digits only, unit stripped)
_NAMES = {
    'lon': {'lon', 'long', 'longitude', 'lng', 'londeg', 'longitudedeg', 'elon', 'lonwgs84'},
    'lat': {'lat', 'latitude', 'latdeg', 'latitudedeg', 'nlat', 'latwgs84'},
    've': {'ve', 'vele', 'vee', 'evel', 'veast', 'velocityeast', 'velocitye', 'vx', 'dedt',
           'east', 'eastvel', 'veloe', 'velx', 'vlon', 'veastmm'},
    'vn': {'vn', 'veln', 'vnn', 'nvel', 'vnorth', 'velocitynorth', 'velocityn', 'vy', 'dndt',
           'north', 'northvel', 'velon', 'vely', 'vlat'},
    'vu': {'vu', 'velu', 'vuu', 'uvel', 'vup', 'velocityup', 'velocityu', 'vz', 'dudt', 'up',
           'upvel', 'velou', 'velz', 'vv', 'vvert', 'vertical', 'vertvel', 'vh'},
    'se': {'se', 'sige', 'sigmae', 'sigve', 'sigmave', 'sdve', 'stdve', 'stde', 'sde', 'dve',
           'errve', 'veerr', 'vesig', 'vesigma', 'vestd', 'eerr', 'errore', 'esig', 'esigma',
           'uncertaintye', 'sigeast', 'stdeast', 'erre', 'sve', 'vesd', 'sigvx', 'sigmavx',
           'sigx', 'vestderr', 'sigmaeast', 'eunc', 'veunc'},
    'sn': {'sn', 'sign', 'sigman', 'sigvn', 'sigmavn', 'sdvn', 'stdvn', 'stdn', 'sdn', 'dvn',
           'errvn', 'vnerr', 'vnsig', 'vnsigma', 'vnstd', 'nerr', 'errorn', 'nsig', 'nsigma',
           'uncertaintyn', 'signorth', 'stdnorth', 'errn', 'svn', 'vnsd', 'sigvy', 'sigmavy',
           'sigy', 'vnstderr', 'sigmanorth', 'nunc', 'vnunc'},
    'su': {'su', 'sigu', 'sigmau', 'sigvu', 'sigmavu', 'sdvu', 'stdvu', 'stdu', 'sdu', 'dvu',
           'errvu', 'vuerr', 'vusig', 'vusigma', 'vustd', 'uerr', 'erroru', 'usig', 'usigma',
           'uncertaintyu', 'sigup', 'stdup', 'erru', 'svu', 'vusd', 'sigvz', 'sigmavz', 'sigz',
           'vustderr', 'sigmaup', 'uunc', 'vuunc', 'sigvert'},
    'site': {'site', 'name', 'station', 'sta', 'stat', 'id', 'code', 'stationid', 'sitename',
             'stationname', 'siteid', 'marker', 'fourcharid', 'stn', 'sitecode'},
}
_UNIT_TAIL = re.compile(r'(mmyr|mmy|mmyear|mma|cmyr|cmy|cmyear|myr|myear)$')


##############################################################################
#                                 the stations                               #
##############################################################################
class GNSS:
    """GNSS stations: site, lon, lat, ve, vn, vu and their sigmas (NaN = none)."""

    FIELDS = ('lon', 'lat', 've', 'vn', 'vu', 'se', 'sn', 'su')

    def __init__(self, site, lon, lat, ve, vn, vu=None, se=None, sn=None, su=None,
                 unit='mm/yr', source=''):
        n = len(lon)
        def arr(a):                                            # noqa: E306
            return (np.full(n, np.nan) if a is None
                    else np.asarray(a, dtype=np.float64).reshape(-1))
        self.site = [str(s) for s in (site if site is not None else
                                      [f'S{i + 1}' for i in range(n)])]
        self.lon, self.lat = arr(lon), arr(lat)
        self.ve, self.vn, self.vu = arr(ve), arr(vn), arr(vu)
        self.se, self.sn, self.su = arr(se), arr(sn), arr(su)
        self.unit, self.source = unit, source

    def __len__(self):
        return self.lon.size

    def __repr__(self):
        return f'<GNSS {len(self)} stations, {self.unit}, {os.path.basename(self.source)}>'

    @property
    def has_vu(self):
        """Per station: does it have a vertical velocity (a 3D station)?"""
        return np.isfinite(self.vu)

    def subset(self, index):
        index = np.asarray(index)
        if index.dtype == bool:
            index = np.flatnonzero(index)
        return GNSS([self.site[i] for i in index],
                    *[getattr(self, f)[index] for f in self.FIELDS],
                    unit=self.unit, source=self.source)

    def scaled(self, factor, unit):
        """The same stations with every velocity and sigma times `factor`."""
        out = self.subset(np.arange(len(self)))
        for name in ('ve', 'vn', 'vu', 'se', 'sn', 'su'):
            setattr(out, name, getattr(out, name) * factor)
        out.unit = unit
        return out

    def component(self, comp, unit_vector=None):
        """Velocity and 1-sigma of one component at every station.

        comp        - 've', 'vn', 'vu' or 'los'
        unit_vector - for 'los': (e, n, u) of the LOS unit vector, scalars or
                      one per station (see los_unit_vector)
        A station without the component gives NaN; a missing sigma counts 0.
        For the LOS a 2D station (no vu) is taken as vu = 0.
        """
        sig = {'ve': self.se, 'vn': self.sn, 'vu': self.su}
        if comp in ('ve', 'vn', 'vu'):
            return getattr(self, comp).copy(), np.nan_to_num(sig[comp].copy())
        if comp != 'los':
            raise ValueError(f'unknown GNSS component {comp!r}: one of {", ".join(COMPONENTS)}')
        if unit_vector is None:
            raise ValueError('the LOS needs the unit vector of the line of sight')
        ue, un, uu = (np.broadcast_to(np.asarray(u, dtype=np.float64), self.lon.shape)
                      for u in unit_vector)
        vu = np.where(self.has_vu, self.vu, 0.0)
        value = self.ve * ue + self.vn * un + vu * uu
        sigma = np.sqrt((np.nan_to_num(self.se) * ue) ** 2 + (np.nan_to_num(self.sn) * un) ** 2
                        + (np.nan_to_num(self.su) * uu) ** 2)
        return value, sigma


##############################################################################
#                                  reading                                   #
##############################################################################
def _norm(name):
    key = re.sub(r'\(.*?\)|\[.*?\]', '', str(name).lower())
    key = re.sub(r'[^a-z0-9]', '', key)
    return _UNIT_TAIL.sub('', key) or key


def _unit_from_header(names):
    """mm/yr, cm/yr or m/yr from the header of the velocity columns, if said."""
    for name in names:
        low = str(name).lower().replace(' ', '')
        if re.match(r'^(v|vel|d)', low) or 'vel' in low:
            if 'mm' in low:
                return 'mm/yr'
            if 'cm' in low:
                return 'cm/yr'
            if re.search(r'\(m/|\[m/|m/yr|m/a\b|m/y\b|myr', low):
                return 'm/yr'
    return None


def _match_columns(header):
    """{field: column index} for the columns we know."""
    keys = [_norm(h) for h in header]
    found = {}
    for field, names in _NAMES.items():
        for i, key in enumerate(keys):
            if key in names and i not in found.values():
                found[field] = i
                break
    if 'lon' not in found or 'lat' not in found:                # x / y, if nothing better
        for field, alias in (('lon', 'x'), ('lat', 'y')):
            if field not in found and alias in keys:
                found[field] = keys.index(alias)
    return found


def _to_float(values):
    out = np.full(len(values), np.nan)
    for i, val in enumerate(values):
        try:
            out[i] = float(str(val).strip().replace('−', '-'))
        except (TypeError, ValueError):
            pass
    return out


def _is_number(text):
    try:
        float(str(text).strip())
        return True
    except ValueError:
        return False


def _split(line, delim):
    if delim is None:
        return line.split()
    return [cell.strip().strip('"').strip("'") for cell in line.split(delim)]


def _read_text(fname):
    """(header or None, rows) of a delimited or whitespace text table."""
    with open(fname, encoding='utf-8-sig', errors='replace') as fid:
        lines = [ln.rstrip('\n\r') for ln in fid]
    comments = [ln for ln in lines if ln.strip().startswith(('#', '%', '>'))]
    body = [ln for ln in lines if ln.strip() and not ln.strip().startswith(('#', '%', '>'))]
    if not body:
        raise ValueError(f'{fname}: no data')
    sample = body[0]
    delim = next((d for d in (',', '\t', ';', '|') if d in sample), None)
    rows = [_split(ln, delim) for ln in body]
    header = None
    first = rows[0]
    if sum(not _is_number(c) for c in first) > len(first) / 2:
        header, rows = first, rows[1:]
    elif comments:                          # "# lon lat ve vn ..." above the data
        last = comments[-1].strip().lstrip('#%> ').strip()
        cand = _split(last, delim)
        if len(cand) >= min(len(r) for r in rows[:5]) and _match_columns(cand).keys() >= {
                'lon', 'lat'}:
            header = cand
    return header, rows


def _read_xlsx(fname):
    """The first worksheet of an .xlsx as rows of text - no extra dependency."""
    import xml.etree.ElementTree as ET
    import zipfile

    main = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    rel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    ns = {'m': main}
    with zipfile.ZipFile(fname) as zfile:
        names = zfile.namelist()
        shared = []
        if 'xl/sharedStrings.xml' in names:
            root = ET.fromstring(zfile.read('xl/sharedStrings.xml'))
            for item in root.findall('m:si', ns):
                shared.append(''.join(t.text or '' for t in item.iter(f'{{{main}}}t')))
        sheet = None
        try:                                # the first sheet of the workbook, by its id
            book = ET.fromstring(zfile.read('xl/workbook.xml'))
            first = book.find('m:sheets/m:sheet', ns)
            rid = first.get(f'{{{rel}}}id')
            rels = ET.fromstring(zfile.read('xl/_rels/workbook.xml.rels'))
            for item in rels:
                if item.get('Id') == rid:
                    target = item.get('Target').lstrip('/')
                    sheet = target if target.startswith('xl/') else 'xl/' + target
        except (KeyError, AttributeError, ET.ParseError):
            sheet = None
        if sheet not in names:
            sheets = sorted(n for n in names if n.startswith('xl/worksheets/sheet'))
            if not sheets:
                raise ValueError(f'{fname}: no worksheet')
            sheet = sheets[0]
        root = ET.fromstring(zfile.read(sheet))

    rows = []
    for row in root.iter(f'{{{main}}}row'):
        cells = {}
        for cell in row.findall('m:c', ns):
            ref = cell.get('r') or ''
            letters = re.match(r'[A-Z]+', ref)
            col = 0
            if letters:
                for ch in letters.group(0):
                    col = col * 26 + ord(ch) - 64
                col -= 1
            else:
                col = len(cells)
            kind, val = cell.get('t'), cell.find('m:v', ns)
            if kind == 's' and val is not None:
                text = shared[int(val.text)]
            elif kind == 'inlineStr':
                text = ''.join(t.text or '' for t in cell.iter(f'{{{main}}}t'))
            else:
                text = val.text if val is not None and val.text is not None else ''
            cells[col] = text
        if cells:
            rows.append([cells.get(i, '') for i in range(max(cells) + 1)])
    rows = [r for r in rows if any(str(c).strip() for c in r)]
    if not rows:
        raise ValueError(f'{fname}: the first sheet is empty')
    first = rows[0]
    if sum(not _is_number(c) for c in first if str(c).strip()) > len(first) / 2:
        return first, rows[1:]
    return None, rows


def _from_records(records, coords, crs, fname):
    """Point features (attribute dicts + x/y) -> header + rows, in lon/lat."""
    if not records:
        raise ValueError(f'{fname}: no point feature')
    xy = np.asarray(coords, dtype=np.float64)
    if crs is not None:
        from pyproj import CRS, Transformer
        src = CRS.from_user_input(crs)
        if not src.is_geographic or src.to_epsg() not in (4326, None):
            lon, lat = Transformer.from_crs(src, 4326, always_xy=True).transform(xy[:, 0],
                                                                                xy[:, 1])
            xy = np.c_[lon, lat]
    header = list(records[0].keys())
    for rec in records[1:]:
        header += [k for k in rec if k not in header]
    rows = [[rec.get(k, '') for k in header] for rec in records]
    # the point geometry wins over any lon / lat attribute
    header = ['__lon__', '__lat__'] + header
    rows = [[x, y] + r for (x, y), r in zip(xy, rows)]
    return header, rows


def _read_shp(fname):
    import shapefile  # pyshp

    records, coords = [], []
    with shapefile.Reader(fname) as reader:
        fields = [f[0] for f in reader.fields[1:]]
        for shape_rec in reader.iterShapeRecords():
            pts = shape_rec.shape.points
            if not pts:
                continue
            records.append(dict(zip(fields, shape_rec.record)))
            coords.append(pts[0][:2])
    crs = None
    stem = os.path.splitext(fname)[0]
    for prj in (stem + '.prj', stem + '.PRJ'):
        if os.path.isfile(prj):
            with open(prj, encoding='utf-8', errors='replace') as fid:
                crs = fid.read().strip() or None
            break
    return _from_records(records, coords, crs, fname)


def _read_geojson(fname):
    with open(fname, encoding='utf-8') as fid:
        obj = json.load(fid)
    feats = obj.get('features', [obj] if obj.get('type') == 'Feature' else [])
    records, coords = [], []
    for feat in feats:
        geom = feat.get('geometry') or {}
        pts = geom.get('coordinates')
        if geom.get('type') == 'MultiPoint' and pts:
            pts = pts[0]
        if geom.get('type') not in ('Point', 'MultiPoint') or not pts:
            continue
        records.append(dict(feat.get('properties') or {}))
        coords.append(pts[:2])
    crs = None
    name = ((obj.get('crs') or {}).get('properties') or {}).get('name')
    if name:
        match = re.search(r'EPSG::?(\d+)', str(name))
        crs = int(match.group(1)) if match else None
        crs = None if crs in (4326, None) else crs
    return _from_records(records, coords, crs, fname)


def _read_kml(fname):
    import xml.etree.ElementTree as ET

    if fname.lower().endswith('.kmz'):
        import zipfile
        with zipfile.ZipFile(fname) as zfile:
            name = next(n for n in zfile.namelist() if n.lower().endswith('.kml'))
            text = zfile.read(name)
    else:
        with open(fname, 'rb') as fid:
            text = fid.read()
    root = ET.fromstring(text)
    records, coords = [], []
    for mark in root.iter():
        if not mark.tag.endswith('Placemark'):
            continue
        point = next((e for e in mark.iter() if e.tag.endswith('Point')), None)
        if point is None:
            continue
        text_c = next((e.text for e in point.iter() if e.tag.endswith('coordinates')), '')
        parts = (text_c or '').strip().split(',')
        if len(parts) < 2:
            continue
        rec = {}
        for elem in mark.iter():
            tag = elem.tag.split('}')[-1]
            if tag == 'name' and elem.text:
                rec['name'] = elem.text.strip()
            elif tag == 'Data' and elem.get('name'):
                val = next((e.text for e in elem if e.tag.endswith('value')), '')
                rec[elem.get('name')] = val
            elif tag == 'SimpleData' and elem.get('name'):
                rec[elem.get('name')] = elem.text
        records.append(rec)
        coords.append([float(parts[0]), float(parts[1])])
    return _from_records(records, coords, None, fname)


def _read_ogr_points(fname):
    """Any OGR-readable point layer: pyogrio, else fiona."""
    from .wkb import wkb_to_geojson

    try:
        from pyogrio import raw
    except ImportError:
        pass
    else:
        meta, _, wkbs, fields = raw.read(fname)[:4]
        names = list(meta.get('fields', []))
        records, coords = [], []
        for i, wkb in enumerate(wkbs):
            if wkb is None:
                continue
            geom = wkb_to_geojson(wkb)
            pts = geom.get('coordinates')
            if geom['type'] == 'MultiPoint' and pts:
                pts = pts[0]
            elif geom['type'] != 'Point':
                continue
            records.append({n: fields[j][i] for j, n in enumerate(names)})
            coords.append(pts[:2])
        return _from_records(records, coords, meta.get('crs'), fname)
    try:
        import fiona
    except ImportError as err:
        raise ValueError(f'cannot read {os.path.basename(fname)}: read directly are text '
                         'tables, .xlsx, .shp, .geojson and .kml; for other formats install '
                         'pyogrio (`pip install pyogrio`) or fiona.') from err
    records, coords = [], []
    with fiona.open(fname) as src:
        crs = src.crs_wkt or None
        for feat in src:
            geom = feat['geometry']
            if geom is None or geom['type'] not in ('Point', 'MultiPoint'):
                continue
            pts = geom['coordinates'] if geom['type'] == 'Point' else geom['coordinates'][0]
            records.append(dict(feat['properties']))
            coords.append(pts[:2])
    return _from_records(records, coords, crs, fname)


def read_gnss(fname, unit='auto', verbose=True):
    """GNSS velocities from a table or a point file.

    unit - 'mm/yr', 'cm/yr', 'm/yr', or 'auto': said by the column headers,
           else mm/yr
    Returns a GNSS object (velocities in that unit).
    """
    ext = os.path.splitext(fname)[1].lower()
    if ext in ('.xlsx', '.xlsm'):
        header, rows = _read_xlsx(fname)
    elif ext == '.shp':
        header, rows = _read_shp(fname)
    elif ext in ('.geojson', '.json'):
        header, rows = _read_geojson(fname)
    elif ext in ('.kml', '.kmz'):
        header, rows = _read_kml(fname)
    elif ext in ('.gpkg', '.gml', '.gdb', '.fgb', '.sqlite', '.tab', '.mif'):
        header, rows = _read_ogr_points(fname)
    else:
        header, rows = _read_text(fname)
    if not rows:
        raise ValueError(f'{fname}: no station')

    if header is None:                       # GMT psvelo: lon lat ve vn se sn [corr] [site]
        ncol = min(len(r) for r in rows)
        if ncol < 4:
            raise ValueError(f'{fname}: a file without a header needs at least lon lat ve vn')
        header = ['lon', 'lat', 've', 'vn', 'se', 'sn', 'corr', 'site'][:max(len(r)
                                                                             for r in rows)]
        if verbose:
            print(f'GNSS: {os.path.basename(fname)} has no header - read as GMT psvelo '
                  '(lon lat ve vn se sn corr site)')
    cols = _match_columns(header)
    if '__lon__' in header:                  # a point file: its geometry
        cols['lon'], cols['lat'] = header.index('__lon__'), header.index('__lat__')
    missing = [f for f in ('lon', 'lat') if f not in cols]
    if missing or not ({'ve', 'vn', 'vu'} & cols.keys()):
        raise ValueError(f'{fname}: cannot find the {" / ".join(missing) or "velocity"} '
                         f'column(s) in {header}. Name them lon, lat, ve, vn, vu (and se, '
                         'sn, su for the sigmas).')
    width = max(cols.values()) + 1
    rows = [r + [''] * (width - len(r)) if len(r) < width else r for r in rows]
    take = {f: [r[i] for r in rows] for f, i in cols.items()}
    num = {f: _to_float(v) for f, v in take.items() if f != 'site'}
    if unit in (None, 'auto'):
        unit = _unit_from_header([header[cols[f]] for f in ('ve', 'vn', 'vu') if f in cols])
        unit = unit or 'mm/yr'
    gnss = GNSS(take.get('site'), num['lon'], num['lat'], num.get('ve'), num.get('vn'),
                num.get('vu'), num.get('se'), num.get('sn'), num.get('su'),
                unit=unit, source=fname)
    good = np.isfinite(gnss.lon) & np.isfinite(gnss.lat)
    if not good.all():
        gnss = gnss.subset(good)
    if not len(gnss):
        raise ValueError(f'{fname}: no station with a position')
    if np.any(np.abs(gnss.lat) > 90) or np.any(np.abs(gnss.lon) > 360):
        raise ValueError(f'{fname}: the coordinates are not lon/lat in degrees')
    if verbose:
        comps = [f for f in ('ve', 'vn', 'vu') if f in cols]
        print(f'GNSS: {len(gnss)} station(s) from {os.path.basename(fname)} '
              f'({", ".join(comps)}; {gnss.unit}; {int(gnss.has_vu.sum())} with vu)')
    return gnss


##############################################################################
#                            the line of sight                               #
##############################################################################
def heading_to_azimuth(heading, look_direction='right'):
    """The azimuth angle of the LOS (MintPy convention) from the orbit heading."""
    heading = np.asarray(heading, dtype=np.float64)
    az = (heading - 90.0) * -1 if look_direction == 'right' else (heading + 90.0) * -1
    return az - np.round(az / 360.0) * 360.0


# nominal orbits: inclination (deg), ground-track speed (km/s), mid-swath
# incidence (deg) - Sentinel-1 IW spans 29.1-46.0 deg, a 350 km ALOS-2 ScanSAR
# swath from 628 km some 25-49 deg
NOMINAL_ORBITS = {
    'Sentinel-1': (98.18, 6.77, 37.5),
    'ALOS-2': (97.92, 6.85, 38.0),
    None: (98.0, 6.8, 38.0),
}


def nominal_heading(lat, ascending, inclination=98.0, ground_speed=6.8):
    """The heading (deg from north, clockwise) of a sun-synchronous ground
    track at latitude `lat`: the orbit's own direction, turned west by the
    rotation of the Earth (about -12 ascending, -168 descending, near the
    equator, for Sentinel-1)."""
    phi, inc = np.radians(lat), np.radians(inclination)
    alpha = np.arcsin(np.clip(np.cos(inc) / max(np.cos(phi), 1e-6), -1.0, 1.0))
    v_east = ground_speed * np.sin(alpha) - 0.4651 * np.cos(phi)
    v_north = ground_speed * np.cos(alpha) * (1.0 if ascending else -1.0)
    return float(np.degrees(np.arctan2(v_east, v_north)))


def _name_tokens(path, levels=2):
    """The words of a file name and of its folders' names (nearest first)."""
    path = os.path.abspath(str(path))
    names = [os.path.splitext(os.path.basename(path))[0]]
    folder = os.path.dirname(path)
    for _ in range(levels):
        names.append(os.path.basename(folder))
        folder = os.path.dirname(folder)
    return [[t for t in re.split(r'[^a-z0-9]+', n.lower()) if t] for n in names]


def _pass_of(tokens):
    """'asc' / 'desc' from the words of one name, or None (neither or both)."""
    found = set()
    for t in tokens:
        if (t in ('asc', 'ascend', 'ascending') or t.endswith('asc')
                or (t.startswith('asc') and not t.startswith('ascii'))):
            found.add('asc')
        if (t in ('desc', 'dsc', 'des', 'descend', 'descending') or t.endswith(('desc', 'dsc'))
                or (t.startswith(('desc', 'dsc')) and not t.startswith('descr'))):
            found.add('desc')
        frame = re.fullmatch(r'(?:\d{3}|t\d{1,3})([ad])', t)  # LiCSAR 105D, track T125A
        if frame:
            found.add('asc' if frame.group(1) == 'a' else 'desc')
    return found.pop() if len(found) == 1 else None


def _mission_of(tokens):
    for i, t in enumerate(tokens):
        if t.startswith(('alos2', 'palsar2', 'scansar')) or t in ('alos', 'palsar', 'wd1', 'wd2'):
            return 'ALOS-2'
        if t.startswith(('sentinel', 's1a', 's1b', 's1c', 'licsar')) or t == 's1':
            return 'Sentinel-1'
        if (re.fullmatch(r'\d{3}[ad]', t) and len(tokens) > i + 2      # 105D_09094_131313
                and re.fullmatch(r'\d{5}', tokens[i + 1]) and re.fullmatch(r'\d{6}', tokens[i + 2])):
            return 'Sentinel-1'
    return None


def nominal_geometry(data_file, atr):
    """(incidence, azimuth, description) of the LOS from what the file's name
    (or its folders', or its ORBIT_DIRECTION / PLATFORM) says - the orbit
    direction, and the mission - or None when it does not say the direction."""
    atr = atr or {}
    words = _name_tokens(data_file)
    direction = str(atr.get('ORBIT_DIRECTION', '')).strip().lower()
    orbit = ('asc' if direction.startswith('asc') else 'desc' if direction.startswith('desc')
             else None)
    said = 'ORBIT_DIRECTION' if orbit else None
    for tokens in words:                       # the file's own name first
        if orbit is None:
            orbit = _pass_of(tokens)
            said = 'its name' if tokens is words[0] else 'its folder'
    if orbit is None:
        return None
    platform = str(atr.get('PLATFORM', '')).lower().replace('-', '').replace('_', '')
    mission = _mission_of([platform]) if platform else None
    for tokens in words:
        mission = mission or _mission_of(tokens)
    incl, speed, inc = NOMINAL_ORBITS[mission]
    look = 'left' if (str(atr.get('ANTENNA_SIDE', '-1')).strip() == '1'
                      or any(t in ('left', 'leftlooking') for t in words[0])) else 'right'
    lat = 0.0
    try:
        lat = float(atr['Y_FIRST']) + float(atr['Y_STEP']) * int(atr['LENGTH']) / 2.0
        if abs(lat) > 90:                       # a projected grid: near enough the equator
            lat = 0.0
    except (KeyError, TypeError, ValueError):
        pass
    heading = nominal_heading(lat, orbit == 'asc', incl, speed)
    az = float(heading_to_azimuth(heading, look))
    what = (f'{mission or "a polar-orbiting SAR"}, '
            f'{"ascending" if orbit == "asc" else "descending"}, {look}-looking - from {said}')
    return inc, az, what


def los_unit_vector(inc_angle, az_angle):
    """(east, north, up) of the unit vector from the ground to the satellite.

    inc_angle - incidence angle from the vertical, degrees
    az_angle  - azimuth angle of that vector, from north, anti-clockwise
                positive, degrees (MintPy's azimuthAngle)
    """
    inc, az = np.deg2rad(inc_angle), np.deg2rad(az_angle)
    return -np.sin(inc) * np.sin(az), np.sin(inc) * np.cos(az), np.cos(inc)


GEOMETRY_KINDS = ("MintPy's geometryGeo.h5 / geometryRadar.h5; LiCSBAS's E.geo, N.geo, "
                  "U.geo with EQA.dem_par (a GEOCml* folder); E / N / U GeoTIFFs "
                  "({frame}.geo.E.tif, E.geo.tif, E.tif, *_east.tif ...); incidence + "
                  "azimuth GeoTIFFs (incidenceAngle.tif + azimuthAngle.tif, inc.tif + az.tif "
                  "...); or a folder that holds any of them")


def _is_one_geometry(folder):
    """A folder that is itself one geometry (LiCSBAS's GEOCml*/), not a place to search."""
    from .licsbas import has_enu_geo
    return has_enu_geo(folder)


class LOSGeometry:
    """The incidence and azimuth angles of a data file, at any lon/lat.

    From MintPy's geometry file (incidenceAngle, azimuthAngle) - geocoded, or
    radar-coded with its latitude / longitude - or from the east, north and up
    of the line of sight that LiCSBAS (GEOCml*/E.geo, N.geo, U.geo) and LiCSAR
    ({frame}.geo.E.tif ...) provide; or, failing those, one value from the
    data's own attributes (CENTER_INCIDENCE_ANGLE, HEADING).
    """

    NAMES = ('geometryGeo.h5', 'geometryRadar.h5')
    NO_GEOMETRY = ('the LOS needs the incidence and azimuth angles, and none were found: no '
                   "MintPy geometry file (geometryGeo.h5), no LiCSBAS GEOCml*/E.geo, N.geo, "
                   "U.geo, no LiCSAR *.geo.E.tif next to the data, no "
                   'CENTER_INCIDENCE_ANGLE / HEADING in it, and its name does not say '
                   'ascending or descending (asc, desc ...). Give --geom (any of those) or '
                   '--los-inc and --los-az')
    NOMINAL_NOTE = ('NOTE: no LOS geometry found for {name} - nominal angles instead, for '
                    '{what}: incidence {inc:.1f}, azimuth {az:.1f} deg, one value for the '
                    'whole scene.\n      Approximate: the incidence changes by 15-25 deg '
                    'across a swath, so the GNSS in the LOS can be off by a few mm/yr\n      '
                    'towards its edges. --geom (MintPy geometryGeo.h5, LiCSBAS GEOCml*/, '
                    'LiCSAR *.geo.E.tif) gives the real angles.')

    def __init__(self, inc=None, az=None, atr=None, lat=None, lon=None, source='',
                 constant=False, crs=None, kind=''):
        self.inc, self.az, self.atr, self.lat, self.lon = inc, az, atr, lat, lon
        self.source, self.constant, self.crs, self.kind = source, constant, crs, kind
        self._lookup = None

    @classmethod
    def for_file(cls, data_file, atr, geom_file=None, inc=None, az=None, verbose=True):
        if inc is not None and az is not None:
            return cls(float(inc), float(az), source=f'given: incidence {inc}, azimuth {az}',
                       constant=True, kind='given')
        if geom_file and os.path.isdir(geom_file) and not _is_one_geometry(geom_file):
            return cls._from_folder(geom_file, data_file, atr, verbose)
        for path in ([geom_file] if geom_file else cls._candidates(data_file)):
            try:
                geom = cls.from_path(path)
                if not geom_file and not geom.covers(atr, data_file):
                    continue
            except (OSError, KeyError, ValueError) as err:
                if geom_file:
                    raise ValueError(f'{geom_file}: {err}') from err
                continue
            if verbose:
                print(f'LOS geometry of {os.path.basename(data_file)}: {geom.source}')
            return geom
        try:
            inc0 = float(atr.get('CENTER_INCIDENCE_ANGLE', atr.get('INCIDENCE_ANGLE')))
            head = float(atr['HEADING'])
        except (KeyError, TypeError, ValueError) as err:
            nominal = nominal_geometry(data_file, atr)
            if nominal is None:
                raise ValueError(cls.NO_GEOMETRY) from err
            inc0, az0, what = nominal
            if verbose:
                print(cls.NOMINAL_NOTE.format(name=os.path.basename(data_file), what=what,
                                              inc=inc0, az=az0))
            return cls(inc0, az0, source=f'nominal angles: {what}', constant=True,
                       kind='nominal')
        look = 'left' if str(atr.get('ANTENNA_SIDE', '-1')).strip() == '1' else 'right'
        az0 = float(heading_to_azimuth(head, look))
        if verbose:
            print(f'LOS geometry: one value for the whole scene, from the attributes '
                  f'(incidence {inc0:.1f}, azimuth {az0:.1f}) - approximate')
        return cls(inc0, az0, source='the data attributes (one value)', constant=True,
                   kind='attributes')

    @classmethod
    def _candidates(cls, data_file):
        from .licsbas import geometry_candidates

        folder = os.path.dirname(os.path.abspath(data_file))
        prefix = 'geo_' if os.path.basename(data_file).startswith('geo_') else ''
        out = []
        up = os.path.join('..', '')
        for sub in ('', 'inputs', up + 'inputs', 'geo', up + 'geo', 'mintpy/inputs',
                    'mintpy/geo', up + 'mintpy/inputs', up + 'mintpy/geo', up + up + 'inputs',
                    up + up + 'geo'):
            for name in cls.NAMES:
                for cand in dict.fromkeys([prefix + name, 'geo_' + name, name]):
                    path = os.path.normpath(os.path.join(folder, sub, cand))
                    if os.path.isfile(path) and path not in out:
                        out.append(path)
        return out + [p for p in geometry_candidates(data_file) if p not in out]

    @classmethod
    def _from_folder(cls, folder, data_file, atr, verbose=True):
        """--geom FOLDER: the first geometry in it (or below it) that covers the data."""
        from .licsbas import geometry_in_folder

        found, tried = geometry_in_folder(folder), []
        for path in found:
            try:
                geom = cls.from_path(path)
            except (OSError, KeyError, ValueError) as err:
                tried.append(f'{os.path.relpath(path, folder)}: {err}')
                continue
            if geom.covers(atr, data_file):
                if verbose:
                    print(f'LOS geometry of {os.path.basename(data_file)}: {geom.source} '
                          f'(found in --geom {folder})')
                return geom
            tried.append(f'{os.path.relpath(path, folder)}: does not cover the data')
        if not found:
            raise ValueError(f'--geom {folder}: no LOS geometry in it or in its subfolders - '
                             + GEOMETRY_KINDS)
        raise ValueError(f'--geom {folder}: none of its geometries covers '
                         f'{os.path.basename(data_file)} - ' + '; '.join(tried))

    @classmethod
    def from_path(cls, path):
        """MintPy's geometry file; LiCSBAS / LiCSAR E, N, U (a folder or a
        file); an incidence / azimuth GeoTIFF pair; or a folder holding any."""
        from .licsbas import geometry_in_folder, is_angle_tif, is_enu_source

        if str(path).lower().endswith(('.h5', '.he5', '.hdf5')):
            return cls.from_geometry(path)
        if is_enu_source(path):
            return cls.from_enu(path)
        if os.path.isdir(path):
            found = geometry_in_folder(path)
            if found:
                return cls.from_path(found[0])
        else:
            pair = is_angle_tif(path)
            if pair:
                return cls.from_angles(pair)
        raise ValueError('not a geometry: give ' + GEOMETRY_KINDS)

    @classmethod
    def from_angles(cls, pair):
        from .licsbas import read_angle_tifs
        from .raster import crs_of

        inc, az, atr = read_angle_tifs(pair)
        return cls(inc, az, atr=atr, source=f'{pair[0]} + {os.path.basename(pair[1])}',
                   crs=crs_of(atr, pair[0]), kind='angles')

    @classmethod
    def from_enu(cls, path):
        from .licsbas import enu_to_angles, read_enu
        from .raster import crs_of

        east, north, up, atr, source = read_enu(path)
        inc, az = enu_to_angles(east, north, up)
        crs = crs_of(atr, atr.get('FILE_PATH'))
        return cls(inc, az, atr=atr, source=source, crs=crs, kind='enu')

    @classmethod
    def from_geometry(cls, path):
        import h5py

        from .raster import crs_of

        with h5py.File(path, 'r') as fid:
            if 'incidenceAngle' not in fid or 'azimuthAngle' not in fid:
                raise KeyError('no incidenceAngle / azimuthAngle')
            inc = np.asarray(fid['incidenceAngle'], dtype=np.float64)
            az = np.asarray(fid['azimuthAngle'], dtype=np.float64)
            atr = {k: (v.decode() if isinstance(v, bytes) else str(v))
                   for k, v in fid.attrs.items()}
            lat = lon = None
            if 'Y_FIRST' not in atr:
                if 'latitude' not in fid or 'longitude' not in fid:
                    raise KeyError('a radar-coded geometry without latitude / longitude')
                lat = np.asarray(fid['latitude'], dtype=np.float64)
                lon = np.asarray(fid['longitude'], dtype=np.float64)
        inc[(inc == 0) | ~np.isfinite(inc)] = np.nan
        crs = crs_of(atr, path) if 'Y_FIRST' in atr else None
        return cls(inc, az, atr=atr, lat=lat, lon=lon, source=path, crs=crs, kind='mintpy')

    def covers(self, atr, data_file=None):
        """Does the geometry reach the middle of a (geocoded) data file?"""
        if self.constant or self.atr is None or 'Y_FIRST' not in self.atr \
                or 'Y_FIRST' not in atr:
            return True
        from .raster import crs_of, same_crs

        x = float(atr['X_FIRST']) + float(atr['X_STEP']) * int(atr['WIDTH']) / 2.0
        y = float(atr['Y_FIRST']) + float(atr['Y_STEP']) * int(atr['LENGTH']) / 2.0
        crs_d = crs_of(atr, data_file)
        crs_g = self.crs if self.crs is not None else 4326
        if crs_d is not None and not same_crs(crs_d, crs_g):
            from pyproj import Transformer
            x, y = Transformer.from_crs(crs_d, crs_g, always_xy=True).transform(x, y)
        row = (y - float(self.atr['Y_FIRST'])) / float(self.atr['Y_STEP'])
        col = (x - float(self.atr['X_FIRST'])) / float(self.atr['X_STEP'])
        return bool(0 <= row < self.inc.shape[0] and 0 <= col < self.inc.shape[1])

    def at_rowcol(self, rows, cols):
        """The angles at pixels of the geometry's own grid (nearest)."""
        if self.constant:
            n = np.size(rows)
            return np.full(n, self.inc), np.full(n, self.az)
        rows = np.clip(np.round(np.asarray(rows, dtype=np.float64)).astype(int), 0,
                       self.inc.shape[0] - 1)
        cols = np.clip(np.round(np.asarray(cols, dtype=np.float64)).astype(int), 0,
                       self.inc.shape[1] - 1)
        return self.inc[rows, cols], self.az[rows, cols]

    def at_lonlat(self, lon, lat):
        """The angles at lon/lat points (nearest pixel; NaN off the grid)."""
        lon, lat = np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)
        if self.constant:
            return np.full(lon.shape, self.inc), np.full(lon.shape, self.az)
        if self.lat is not None:                         # radar-coded: through lat/lon
            if self._lookup is None:
                from .radar import RadarLookup
                self._lookup = RadarLookup.__new__(RadarLookup)
                self._lookup.lat, self._lookup.lon = self.lat, self.lon
                self._lookup.length, self._lookup.width = self.lat.shape
                self._lookup._interp = None
            rows, cols = self._lookup.to_rowcol(lon, lat)
        else:
            from .raster import crs_of, same_crs
            x, y = lon, lat
            crs = self.crs if self.crs is not None else crs_of(self.atr, self.source)
            if crs is not None and not same_crs(crs, 4326):
                from pyproj import Transformer
                x, y = Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat)
            rows = (np.asarray(y) - float(self.atr['Y_FIRST'])) / float(self.atr['Y_STEP']) - 0.5
            cols = (np.asarray(x) - float(self.atr['X_FIRST'])) / float(self.atr['X_STEP']) - 0.5
        inside = (np.isfinite(rows) & np.isfinite(cols) & (rows > -0.5) & (cols > -0.5)
                  & (rows < self.inc.shape[0] - 0.5) & (cols < self.inc.shape[1] - 0.5))
        inc, az = np.full(lon.shape, np.nan), np.full(lon.shape, np.nan)
        if inside.any():
            inc[inside], az[inside] = self.at_rowcol(rows[inside], cols[inside])
        return inc, az

    def unit_vector_at(self, lon, lat):
        return los_unit_vector(*self.at_lonlat(lon, lat))

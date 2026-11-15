"""LiCSBAS and LiCSAR products: their line-of-sight geometry and their unit.

LiCSAR (COMET) ships, with every frame, the east, north and up components of
the unit vector from the ground to the satellite, as GeoTIFFs
({frame}.geo.E.tif, .N.tif, .U.tif). LiCSBAS multilooks them into its
GEOCml*/ folder as E.geo, N.geo and U.geo - float32 rasters on the grid of
GEOCml*/EQA.dem_par, whose corner_lat / corner_lon is the centre of the first
pixel. Turned into MintPy's incidence and azimuth angles, they give the LOS of
every pixel, so `--comp los` and `--vel los` need nothing more.

A LiCSBAS velocity (vel, vel.filt.mskd ..., and their GeoTIFFs written by
LiCSBAS_flt2geotiff.py) is in mm/yr but says so nowhere.
"""

import glob
import os
import re

import numpy as np

ENU = ('E', 'N', 'U')
# a set of E, N, U GeoTIFFs: {frame}.geo.E.tif (LiCSAR), E.geo.tif
# (LiCSBAS_flt2geotiff.py of E.geo), E.tif, los_e.tif, los_east.tif ...
_TIF_RE = re.compile(r'^(?P<head>.*?)(?P<comp>east|north|up|[ENU])(?P<mid>\.geo)?'
                     r'(?P<tail>\.tiff?)$', re.IGNORECASE)
_WORDS = {'E': 'east', 'N': 'north', 'U': 'up'}
# incidence / azimuth GeoTIFFs: incidenceAngle.tif + azimuthAngle.tif (MintPy's
# save_gdal.py), inc.tif + az.tif, los_incidence.tif + los_azimuth.tif ...
_INC_WORDS = ('incidenceangle', 'incidence_angle', 'incidence', 'inc')
_AZ_WORDS = ('azimuthangle', 'azimuth_angle', 'azimuth', 'azi', 'az')
_RASTER_EXT = ('.tif', '.tiff')
_GEOM_H5 = ('geometryGeo.h5', 'geo_geometryRadar.h5', 'geo_geometryGeo.h5',
            'geometryRadar.h5')
_PAIR_DIR = re.compile(r'^\d{8}_\d{8}$')     # GEOC/20180101_20180113/: an interferogram
_SEARCH_DEPTH = 4            # how far up a data file may sit below the LiCSBAS folder


##############################################################################
#                                 the grid                                   #
##############################################################################
def read_dem_par(path):
    """The fields of a GAMMA / LiCSBAS EQA.dem_par file: {key: first word}."""
    par = {}
    with open(path, encoding='utf-8', errors='replace') as fid:
        for line in fid:
            if ':' not in line:
                continue
            key, _, value = line.partition(':')
            words = value.split()
            par[key.strip()] = words[0] if words else ''
    return par


def grid_attrs(par):
    """MintPy-style attributes (edges of the first pixel) of an EQA.dem_par grid."""
    try:
        width, length = int(par['width']), int(par['nlines'])
        lat0, lon0 = float(par['corner_lat']), float(par['corner_lon'])
        dlat, dlon = float(par['post_lat']), float(par['post_lon'])
    except (KeyError, ValueError) as err:
        raise ValueError(f'EQA.dem_par: no usable {err} field') from err
    return {'LENGTH': str(length), 'WIDTH': str(width),
            'Y_FIRST': repr(lat0 - dlat / 2.0), 'X_FIRST': repr(lon0 - dlon / 2.0),
            'Y_STEP': repr(dlat), 'X_STEP': repr(dlon),
            'Y_UNIT': 'degrees', 'X_UNIT': 'degrees', 'EPSG': '4326'}


##############################################################################
#                           E, N, U -> the angles                            #
##############################################################################
def enu_to_angles(east, north, up):
    """Incidence and azimuth angles (MintPy's convention: the azimuth of the
    ground-to-satellite vector, from north, anti-clockwise) from the east,
    north and up components of that unit vector. NaN where there is none."""
    e, n, u = (np.asarray(a, dtype=np.float64) for a in (east, north, up))
    bad = ~(np.isfinite(e) & np.isfinite(n) & np.isfinite(u)) | ((e == 0) & (n == 0) & (u == 0))
    inc = np.where(bad, np.nan, np.degrees(np.arccos(np.clip(u, -1.0, 1.0))))
    az = np.where(bad, np.nan, np.degrees(np.arctan2(-e, n)))
    return inc, az


def _plausible(arrays):
    """Do three arrays look like the components of unit vectors?"""
    e, n, u = arrays
    good = np.isfinite(e) & np.isfinite(n) & np.isfinite(u) & ((e != 0) | (n != 0) | (u != 0))
    if not good.any():
        return False
    norm = np.sqrt(e[good] ** 2 + n[good] ** 2 + u[good] ** 2)
    return bool(np.nanmedian(np.abs(norm - 1.0)) < 0.01)


def has_enu_geo(folder):
    return all(os.path.isfile(os.path.join(folder, f'{c}.geo')) for c in ENU) and \
        os.path.isfile(os.path.join(folder, 'EQA.dem_par'))


def read_enu_geo(folder):
    """(E, N, U, attributes) of LiCSBAS's GEOCml*/E.geo, N.geo, U.geo."""
    atr = grid_attrs(read_dem_par(os.path.join(folder, 'EQA.dem_par')))
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    raw = []
    for comp in ENU:
        path = os.path.join(folder, f'{comp}.geo')
        size = os.path.getsize(path)
        if size != 4 * length * width:
            raise ValueError(f'{path}: {size} bytes, not the {length} x {width} float32 of '
                             'its EQA.dem_par')
        raw.append(np.fromfile(path, dtype='<f4').reshape(length, width))
    arrays = [a.astype(np.float64) for a in raw]
    if not _plausible(arrays):            # GAMMA writes big-endian
        swapped = [a.byteswap().astype(np.float64) for a in raw]
        if _plausible(swapped):
            arrays = swapped
        else:
            raise ValueError(f'{folder}: E.geo, N.geo, U.geo are not unit vectors')
    return (*arrays, atr)


def enu_tif_triple(path):
    """The (E, N, U) GeoTIFFs of a set, given any one of them - LiCSAR's
    {frame}.geo.E.tif ..., E.geo.tif, E.tif, los_east.tif - or None when the
    set is incomplete."""
    match = _TIF_RE.match(os.path.basename(path))
    if not match:
        return None
    folder = os.path.dirname(path)
    head, comp = match.group('head'), match.group('comp')
    mid, tail = match.group('mid') or '', match.group('tail')
    triple = []
    for letter in ENU:
        if len(comp) > 1:                                    # east, North, UP ...
            word = _WORDS[letter]
            spellings = (word, word.capitalize(), word.upper())
        else:
            spellings = (letter, letter.lower())
        for name in spellings:
            cand = os.path.join(folder, f'{head}{name}{mid}{tail}')
            if os.path.isfile(cand):
                triple.append(cand)
                break
        else:
            return None
    return tuple(triple)


def _rasters(folder):
    try:
        return sorted(n for n in os.listdir(folder) if n.lower().endswith(_RASTER_EXT)
                      and os.path.isfile(os.path.join(folder, n)))
    except OSError:
        return []


def all_enu_tifs(folder):
    """Every (E, N, U) GeoTIFF set in a folder (one per frame, say)."""
    out = []
    for name in _rasters(folder):
        match = _TIF_RE.match(name)
        if match and match.group('comp').lower() in ('e', 'east'):
            triple = enu_tif_triple(os.path.join(folder, name))
            if triple and triple not in out:
                out.append(triple)
    return out


def find_enu_tifs(folder):
    """The (E, N, U) GeoTIFF set in a folder, or None."""
    found = all_enu_tifs(folder)
    return found[0] if found else None


def _without(name, words):
    """The name with the first of `words` standing on its own taken out (a
    key that pairs incidenceAngle.tif with azimuthAngle.tif), or None."""
    low = name.lower()
    for word in words:
        for m in re.finditer(re.escape(word), low):
            before = low[m.start() - 1] if m.start() else ''
            after = low[m.end()] if m.end() < len(low) else ''
            if not before.isalpha() and not after.isalpha():
                return low[:m.start()] + '\0' + low[m.end():]
    return None


def all_angle_tifs(folder):
    """Every (incidence, azimuth) GeoTIFF pair in a folder."""
    names = _rasters(folder)
    inc = {k: n for n in names if (k := _without(n, _INC_WORDS))}
    az = {k: n for n in names if (k := _without(n, _AZ_WORDS))}
    return [(os.path.join(folder, inc[k]), os.path.join(folder, az[k]))
            for k in sorted(inc) if k in az]


def is_angle_tif(path):
    """Is `path` one of an (incidence, azimuth) GeoTIFF pair? The pair, or None."""
    folder, name = os.path.split(os.path.abspath(str(path)))
    for pair in all_angle_tifs(folder):
        if os.path.basename(pair[0]) == name or os.path.basename(pair[1]) == name:
            return pair
    return None


def read_angle_tifs(pair):
    """(incidence, azimuth, attributes) of an incidence / azimuth GeoTIFF pair -
    the azimuth in MintPy's convention (ground to satellite, anticlockwise)."""
    from .raster import read_any

    inc, atr = read_any(pair[0])
    az, atr_az = read_any(pair[1])
    if (atr_az['LENGTH'], atr_az['WIDTH']) != (atr['LENGTH'], atr['WIDTH']):
        raise ValueError(f'{pair[1]} is not on the grid of {pair[0]}')
    inc, az = np.asarray(inc, dtype=np.float64), np.asarray(az, dtype=np.float64)
    good = np.isfinite(inc) & (inc != 0)
    if not good.any() or not (0 < np.nanmedian(inc[good]) < 90):
        raise ValueError(f'{pair[0]}: not incidence angles in degrees')
    if np.nanmax(np.abs(az[good])) > 360:
        raise ValueError(f'{pair[1]}: not azimuth angles in degrees')
    inc[~good] = np.nan
    atr = dict(atr)
    atr['FILE_PATH'] = os.path.abspath(pair[0])
    return inc, az, atr


def _geometry_here(folder):
    """The LOS geometries kept in one folder, best first."""
    out = []
    if has_enu_geo(folder):
        out.append(folder)
    out += [triple[0] for triple in all_enu_tifs(folder)]
    out += [pair[0] for pair in all_angle_tifs(folder)]
    out += [os.path.join(folder, n) for n in _GEOM_H5 if os.path.isfile(os.path.join(folder, n))]
    return out


def _subfolders(folder):
    """The folders below `folder` where a geometry may be, likeliest first."""
    try:
        names = [n for n in os.listdir(folder) if os.path.isdir(os.path.join(folder, n))
                 and not n.startswith('.') and not _PAIR_DIR.match(n)]
    except OSError:
        return []

    def rank(name):
        for i, pre in enumerate(('GEOCml', 'GEOC', 'metadata', 'inputs', 'geo', 'geom')):
            if name.startswith(pre):
                return (i, name)
        return (9, name)
    return [os.path.join(folder, n) for n in sorted(names, key=rank)]


def geometry_in_folder(folder, depth=2):
    """Every LOS geometry in a folder and up to `depth` folders below it,
    nearest (then likeliest) first: LiCSBAS's E.geo, N.geo, U.geo; E / N / U
    GeoTIFFs ({frame}.geo.E.tif, E.geo.tif, E.tif, *_east.tif ...);
    incidence / azimuth GeoTIFFs; MintPy's geometry file. The interferogram
    folders of LiCSBAS's GEOC/ are skipped."""
    out, level = [], [os.path.abspath(str(folder))]
    for _ in range(depth + 1):
        nxt = []
        for here in level:
            out += [p for p in _geometry_here(here) if p not in out]
            nxt += _subfolders(here)
        level = nxt
    return out


def read_enu_tifs(triple):
    """(E, N, U, attributes) of three GeoTIFFs on one grid."""
    from .raster import read_any

    arrays, atr = [], None
    for path in triple:
        data, atr_i = read_any(path)
        if atr is not None and (atr_i['LENGTH'], atr_i['WIDTH']) != (atr['LENGTH'], atr['WIDTH']):
            raise ValueError(f'{path} is not on the grid of {triple[0]}')
        atr = atr if atr is not None else dict(atr_i)
        arrays.append(np.asarray(data, dtype=np.float64))
    if not _plausible(arrays):
        raise ValueError(f'{triple[0]} ...: not the components of unit vectors')
    atr['FILE_PATH'] = os.path.abspath(triple[0])
    return (*arrays, atr)


##############################################################################
#                          where LiCSBAS keeps them                          #
##############################################################################
def _ancestors(path, depth=_SEARCH_DEPTH):
    folder = os.path.dirname(os.path.abspath(path))
    out = [folder]
    for _ in range(depth):
        parent = os.path.dirname(out[-1])
        if parent == out[-1]:
            break
        out.append(parent)
    return out


def in_licsbas_tree(path):
    """Is the file inside a LiCSBAS work folder (TS_GEOC*/, GEOC*/ ...)?"""
    for folder in _ancestors(path):
        name = os.path.basename(folder)
        if name.startswith(('TS_GEOC', 'GEOC')) or has_enu_geo(folder):
            return True
    return False


def geometry_candidates(data_file):
    """Folders and files that may hold the E, N, U of a data file, nearest
    first: the file's own folder; for a file under TS_GEOCmlX/ its sibling
    GEOCmlX/ (then the other GEOCml*/, then GEOC/ with LiCSAR's GeoTIFFs);
    LiCSAR's {frame}.geo.E.tif next to the file or in metadata/."""
    out = []

    def add(path):
        path = os.path.normpath(path)
        if path not in out:
            out.append(path)

    folders = _ancestors(data_file)
    if has_enu_geo(folders[0]):
        add(folders[0])
    for folder in folders:
        name = os.path.basename(folder)
        if name.startswith('TS_GEOC') or (name.startswith('GEOC') and has_enu_geo(folder)):
            work = os.path.dirname(folder)
            want = name[3:] if name.startswith('TS_') else name
            subs = sorted(glob.glob(os.path.join(work, 'GEOC*')))
            subs.sort(key=lambda p: (os.path.basename(p) != want,
                                     not os.path.basename(p).startswith('GEOCml'), p))
            for sub in subs:
                if not os.path.isdir(sub):
                    continue
                if has_enu_geo(sub):
                    add(sub)
                elif find_enu_tifs(sub):
                    add(sub)
            break
    for folder in folders[:2]:
        for sub in (folder, os.path.join(folder, 'metadata')):
            if os.path.isdir(sub) and find_enu_tifs(sub):
                add(sub)
    return out


def is_enu_source(path):
    """Is `path` LiCSBAS / LiCSAR line-of-sight geometry (a folder or a file)?"""
    if not path:
        return False
    if os.path.isdir(path):
        return has_enu_geo(path) or find_enu_tifs(path) is not None
    base = os.path.basename(path)
    if base in ('E.geo', 'N.geo', 'U.geo', 'EQA.dem_par'):
        return has_enu_geo(os.path.dirname(path) or '.')
    return enu_tif_triple(path) is not None


def read_enu(path):
    """(E, N, U, attributes, description) of the LOS geometry at `path`: a
    GEOCml*/ folder, any of its E.geo / N.geo / U.geo / EQA.dem_par, a folder
    with LiCSAR's GeoTIFFs, or any one of those GeoTIFFs."""
    if os.path.isdir(path):
        if has_enu_geo(path):
            return (*read_enu_geo(path), f'{path} (LiCSBAS E.geo, N.geo, U.geo)')
        triple = find_enu_tifs(path)
        if triple:
            return (*read_enu_tifs(triple), f'{triple[0]} (+ N, U; LiCSAR)')
        raise ValueError(f'{path}: no E.geo, N.geo, U.geo (LiCSBAS) nor *E.tif, *N.tif, *U.tif '
                         '(LiCSAR) in it')
    base = os.path.basename(path)
    if base in ('E.geo', 'N.geo', 'U.geo', 'EQA.dem_par'):
        folder = os.path.dirname(path) or '.'
        return (*read_enu_geo(folder), f'{folder} (LiCSBAS E.geo, N.geo, U.geo)')
    triple = enu_tif_triple(path)
    if triple:
        return (*read_enu_tifs(triple), f'{triple[0]} (+ N, U; LiCSAR)')
    raise ValueError(f'{path}: not LiCSBAS / LiCSAR line-of-sight geometry')


##############################################################################
#                                 the unit                                   #
##############################################################################
def licsbas_unit(data_file, geom_hint=None):
    """mm/year for a LiCSBAS velocity (vel*, vstd*) inside a LiCSBAS tree - or
    with LiCSBAS geometry given for it - else None."""
    base = os.path.basename(str(data_file)).lower()
    if not base.startswith(('vel', 'vstd')):
        return None
    hints = geom_hint if isinstance(geom_hint, (list, tuple)) else [geom_hint]
    if in_licsbas_tree(data_file) or any(is_enu_source(h) for h in hints if h):
        return 'mm/year'
    return None


##############################################################################
#                    LiCSBAS's own files, opened directly                    #
##############################################################################
H5_GRID = ('corner_lat', 'corner_lon', 'post_lat', 'post_lon')
# LiCSBAS's datasets and results/ files, and their units
UNITS = {'vel': 'mm/year', 'vstd': 'mm/year', 'vintercept': 'mm', 'resid_rms': 'mm',
         'cum': 'mm', 'cum_model': 'mm', 'stc': 'mm', 'hgt': 'm'}
_NOT_RAW = ('.h5', '.he5', '.hdf5', '.tif', '.tiff', '.vrt', '.xml', '.rsc', '.png', '.jpg',
            '.kml', '.kmz', '.csv', '.txt', '.json', '.geojson', '.gpkg', '.shp', '.nc',
            '.grd', '.img', '.par', '.py', '.log', '.cc')
_KIND_CACHE = {}


def _unit_of(name):
    """The unit of a LiCSBAS dataset or results/ file, from its name."""
    base = os.path.basename(str(name))
    if base.lower().endswith('.unw'):
        return 'radian'
    stem = base.split('.')[0].split('-')[0]
    return UNITS.get(stem, '1')


def _dem_par_for(fname):
    """The EQA.dem_par of the grid a LiCSBAS file without a header is on."""
    folders = _ancestors(fname, depth=3)
    for folder in folders[:3]:                            # GEOCml*/, GEOCml*/YYYYMMDD_YYYYMMDD/
        path = os.path.join(folder, 'EQA.dem_par')
        if os.path.isfile(path):
            return path
    for folder in folders:
        name = os.path.basename(folder)
        if name.startswith('TS_GEOC'):
            work, want = os.path.dirname(folder), name[3:]
            subs = sorted(glob.glob(os.path.join(work, 'GEOC*')))
            subs.sort(key=lambda p: (os.path.basename(p) != want, p))
            for sub in [os.path.join(folder, 'info')] + subs:
                path = os.path.join(sub, 'EQA.dem_par')
                if os.path.isfile(path):
                    return path
    return None


def kind_of(fname):
    """('h5', None) for LiCSBAS's cum.h5 / cum_filt.h5, ('raw', EQA.dem_par)
    for one of its float32 files without a header (results/vel, vel.filt.mskd,
    GEOCml*/hgt ...), else (None, None)."""
    try:
        path = os.path.abspath(str(fname))
        key = (path, os.path.getmtime(path))
    except OSError:
        return None, None
    if key in _KIND_CACHE:
        return _KIND_CACHE[key]
    out = (None, None)
    low = path.lower()
    if low.endswith(('.h5', '.hdf5')):
        try:
            import h5py
            with h5py.File(path, 'r') as h5:
                if all(k in h5 for k in H5_GRID) and ('vel' in h5 or 'cum' in h5):
                    out = ('h5', None)
        except Exception:            # noqa: BLE001 - not HDF5, or unreadable: not ours
            pass
    elif os.path.isfile(path) and not low.endswith(_NOT_RAW) \
            and not os.path.isfile(path + '.rsc'):
        par = _dem_par_for(path)
        if par:
            try:
                atr = grid_attrs(read_dem_par(par))
                if os.path.getsize(path) == 4 * int(atr['LENGTH']) * int(atr['WIDTH']):
                    out = ('raw', par)
            except ValueError:
                pass
    _KIND_CACHE[key] = out
    return out


def _h5_grid_attrs(h5):
    lat1, lon1 = float(h5['corner_lat'][()]), float(h5['corner_lon'][()])
    dlat, dlon = float(h5['post_lat'][()]), float(h5['post_lon'][()])
    shape = h5['vel'].shape if 'vel' in h5 else h5['cum'].shape[1:]
    return {'LENGTH': str(shape[0]), 'WIDTH': str(shape[1]),
            'Y_FIRST': repr(lat1 - dlat / 2.0), 'X_FIRST': repr(lon1 - dlon / 2.0),
            'Y_STEP': repr(dlat), 'X_STEP': repr(dlon),
            'Y_UNIT': 'degrees', 'X_UNIT': 'degrees', 'EPSG': '4326'}


def _h5_dates(h5):
    if 'imdates' not in h5:
        return []
    return [str(int(d)) for d in np.ravel(h5['imdates'][()])]


def attributes(fname, datasetName=None):
    """MintPy-style metadata of a LiCSBAS file (see kind_of)."""
    kind, par = kind_of(fname)
    if kind == 'h5':
        import h5py
        with h5py.File(fname, 'r') as h5:
            atr = _h5_grid_attrs(h5)
            if 'refarea' in h5:
                ref = h5['refarea'][()]
                atr['REF_AREA'] = ref.decode() if isinstance(ref, bytes) else str(ref)
        atr['UNIT'] = _unit_of(datasetName or 'vel')
    elif kind == 'raw':
        atr = grid_attrs(read_dem_par(par))
        atr['UNIT'] = _unit_of(fname)
        atr['DEM_PAR'] = par
    else:
        raise ValueError(f'{fname}: not a LiCSBAS file')
    atr.update(FILE_TYPE='LiCSBAS', PROCESSOR='LiCSBAS', DATA_TYPE='float32',
               FILE_PATH=os.path.abspath(fname))
    return atr


def slices(fname):
    """The 2D slices of a LiCSBAS file: vel first, the others, cum-YYYYMMDD ..."""
    kind, _ = kind_of(fname)
    if kind == 'raw':
        return [os.path.basename(fname)]
    import h5py
    with h5py.File(fname, 'r') as h5:
        shape = h5['vel'].shape if 'vel' in h5 else h5['cum'].shape[1:]
        names = [k for k, obj in h5.items() if isinstance(obj, h5py.Dataset)
                 and obj.shape == shape and k not in H5_GRID]
        names.sort(key=lambda k: (k != 'vel', k))
        cube = 'cum' if 'cum' in h5 else None
        dates = _h5_dates(h5) if cube else []
    return names + [f'{cube}-{d}' for d in dates]


def read(fname, datasetName=None, box=None):
    """(data, metadata) of a LiCSBAS file; box = (x0, y0, x1, y1)."""
    kind, _ = kind_of(fname)
    atr = attributes(fname, datasetName)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    x0, y0, x1, y1 = box if box else (0, 0, width, length)
    if kind == 'raw':
        data = np.memmap(fname, dtype='<f4', mode='r', shape=(length, width))
        return np.array(data[y0:y1, x0:x1], dtype=np.float32), atr
    import h5py
    name = datasetName or 'vel'
    with h5py.File(fname, 'r') as h5:
        base, _, date = name.partition('-')
        if date and base in h5 and h5[base].ndim == 3:
            dates = _h5_dates(h5)
            if date not in dates:
                raise ValueError(f'{fname}: no date {date} in {base}')
            data = h5[base][dates.index(date), y0:y1, x0:x1]
        elif name in h5 and h5[name].ndim == 2:
            data = h5[name][y0:y1, x0:x1]
        elif base in h5 and h5[base].ndim == 3:               # the last epoch
            data = h5[base][-1, y0:y1, x0:x1]
        else:
            raise ValueError(f'{fname}: no dataset {name!r}; there are '
                             f'{", ".join(slices(fname)[:8])} ...')
    return np.array(data, dtype=np.float32), atr


def mask_for(fname):
    """LiCSBAS's results/mask for one of its files, when it is there."""
    kind, _ = kind_of(fname)
    folder = os.path.dirname(os.path.abspath(fname))
    if kind == 'h5':
        path = os.path.join(folder, 'results', 'mask')
    elif kind == 'raw':
        path = os.path.join(folder, 'mask')
    else:
        return None
    if os.path.basename(path) == os.path.basename(fname) or not os.path.isfile(path):
        return None
    return path

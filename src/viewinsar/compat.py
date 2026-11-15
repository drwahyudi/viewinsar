"""MintPy where it is installed, light stand-ins where it is not.

viewinsar uses a little of MintPy: reading its files, the display unit and
the mask of view.py, its map, and tsview.py behind `--timeseries`. MintPy is
not a requirement: with it installed, all of these come from MintPy - and
every view.py option works on the map. Without it, the stand-ins here read
MintPy HDF5 products (velocity, time series, geometry, masks, HDF-EOS5),
GeoTIFF and every other GDAL raster (with rasterio), and ROI_PAC .dem / .rsc;
mapview.py draws the map, with the view.py options in most use.

Set VIEWINSAR_NO_MINTPY=1 to use the stand-ins even where MintPy is installed.
"""

import functools
import glob
import os
import re

import numpy as np

__all__ = ['have_mintpy', 'read_attribute', 'read', 'get_slice_list', 'get_file_list',
           'coordinate', 'Coordinate', 'auto_disp_unit', 'scale_data2disp_unit', 'read_mask',
           'auto_figure_size', 'auto_flip_direction', 'MINTPY_HINT']

MINTPY_HINT = 'pip install mintpy (or conda install -c conda-forge mintpy)'
_LENGTH_UNITS = {'mm': 0.001, 'cm': 0.01, 'dm': 0.1, 'm': 1.0, 'km': 1000.0,
                 'in': 0.0254, 'inch': 0.0254, 'ft': 0.3048, 'foot': 0.3048,
                 'yd': 0.9144, 'yard': 0.9144, 'mi': 1609.344, 'mile': 1609.344}
_UNIT_OF_DSET = {'height': 'm', 'latitude': 'degree', 'longitude': 'degree',
                 'incidenceAngle': 'degree', 'azimuthAngle': 'degree',
                 'slantRangeDistance': 'm', 'temporalCoherence': '1', 'coherence': '1',
                 'mask': '1', 'waterMask': '1', 'shadowMask': '1', 'velocity': 'm/year',
                 'velocityStd': 'm/year', 'intercept': 'm', 'timeseries': 'm',
                 'displacement': 'm', 'unwrapPhase': 'radian'}


@functools.cache
def _mintpy_importable():
    try:
        from mintpy.utils import readfile  # noqa: F401
        return True
    except Exception:            # noqa: BLE001 - missing, or broken in this environment
        return False


def have_mintpy():
    """Is MintPy installed (and not turned off with VIEWINSAR_NO_MINTPY)?"""
    if os.environ.get('VIEWINSAR_NO_MINTPY', '').strip().lower() in ('1', 'true', 'yes', 'on'):
        return False
    return _mintpy_importable()


@functools.cache
def _gdal_importable():
    try:
        from osgeo import gdal  # noqa: F401
        return True
    except Exception:            # noqa: BLE001
        return False


_GDAL_EXT = ('.tif', '.tiff', '.vrt', '.img', '.jp2', '.grd', '.nc', '.dem.wgs84')


def _licsbas(fname):
    """'h5' / 'raw' for LiCSBAS's own files (cum.h5, results/vel ...), else None."""
    from .licsbas import kind_of
    return kind_of(fname)[0]


def mintpy_reads(fname):
    """Read this file with MintPy? Not LiCSBAS's own files (MintPy knows
    neither their grid nor their unit), nor a GeoTIFF (or another GDAL raster)
    when MintPy has no GDAL to read it with: rasterio does it then."""
    if not have_mintpy() or _licsbas(fname):
        return False
    return not (str(fname).lower().endswith(_GDAL_EXT) and not _gdal_importable())


def is_unitless(unit):
    """No unit: missing, '1' (MintPy's for a GeoTIFF without one) or 'none'."""
    return str(unit or '').strip().lower() in ('', '1', 'none')


_UNIT_RE = re.compile(r'^(mm|cm|dm|m|km)(/(yr|year|y|a))?$')


_LENGTH_M = {'mm': 1e-3, 'cm': 1e-2, 'dm': 1e-1, 'm': 1.0, 'km': 1e3}


def unit_ratio(unit, other):
    """How many `other` one `unit` makes - mm/year in m/year: 0.001 - or None
    when the two do not convert (radian, another time base ...)."""
    def split(text):
        head, _, tail = str(text).strip().lower().replace(' ', '').partition('/')
        return head, ('year' if tail in ('yr', 'year', 'y', 'a') else tail)
    (h1, t1), (h2, t2) = split(unit), split(other)
    if t1 != t2 or h1 not in _LENGTH_M or h2 not in _LENGTH_M:
        return None
    return _LENGTH_M[h1] / _LENGTH_M[h2]


def guess_velocity_unit(values):
    """(unit, p99) of a velocity map that says no unit, from its numbers: in
    m/yr nearly every pixel of a deformation map is below 1 (1 m/yr is
    extreme), in mm/yr the noise alone exceeds 1 - so p99 of |v| below 1 is
    m/year, else mm/year. (None, None) without a finite, non-zero value."""
    v = np.abs(np.asarray(values, dtype=np.float64).ravel())
    v = v[np.isfinite(v) & (v > 0)]
    if not v.size:
        return None, None
    if v.size > 2_000_000:
        v = v[::v.size // 1_000_000]
    p99 = float(np.percentile(v, 99))
    return ('m/year' if p99 < 1.0 else 'mm/year'), p99


def velocity_name(fname, dset=None):
    """Does the file (or dataset) name say it is a velocity - vel*, *velocity*, rate?"""
    names = [os.path.basename(str(fname)).lower(), str(dset or '').lower()]
    tokens = {t for name in names for t in re.split(r'[^a-z0-9]+', name) if t}
    return any(t.startswith('vel') or 'velocity' in t or t in ('rate', 'rates', 'vlos')
               for t in tokens)


def normalise_unit(unit):
    """mm/yr, cm/a, m ... in MintPy's spelling (mm/year, cm/year, m); radian."""
    text = str(unit).strip().lower().replace(' ', '')
    if text in ('rad', 'radian', 'radians'):
        return 'radian'
    match = _UNIT_RE.match(text)
    if not match:
        raise ValueError(f'unknown unit {unit!r}: use mm/yr, cm/yr, m/yr, mm, cm, m or radian')
    return match.group(1) + ('/year' if match.group(2) else '')


##############################################################################
#                                  reading                                   #
##############################################################################
def _decode(value):
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace')
    if isinstance(value, np.ndarray) and value.size == 1:
        return _decode(value.ravel()[0])
    if isinstance(value, (np.generic,)):
        return str(value.item())
    return str(value)


def _is_hdf5(fname):
    return str(fname).lower().endswith(('.h5', '.he5', '.hdf5'))


def _h5_images(fid):
    """[(path, dataset)] of every 2D / 3D dataset, in the order h5py walks them."""
    import h5py

    out = []
    fid.visititems(lambda name, obj: out.append((name, obj))
                   if isinstance(obj, h5py.Dataset) and obj.ndim in (2, 3) else None)
    return out


def _h5_dates(fid, path, num):
    """The dates of the 3D dataset at `path`: the `date` beside it (YYYYMMDD,
    or pairs for an interferogram stack), else 1, 2, 3 ..."""
    group = os.path.dirname(path)
    for cand in (f'{group}/date' if group else 'date', 'date', 'imdates'):
        if cand in fid and fid[cand].shape[0] == num:
            dates = fid[cand][()]
            if dates.ndim == 2:                       # pairs
                return ['_'.join(_decode(d) for d in pair) for pair in dates]
            return [_decode(d) for d in dates]
    return [str(i + 1) for i in range(num)]


def _h5_attrs(fname):
    import h5py

    with h5py.File(fname, 'r') as fid:
        atr = {k: _decode(v) for k, v in fid.attrs.items()}
        images = _h5_images(fid)
        if images and ('LENGTH' not in atr or 'WIDTH' not in atr):
            atr['LENGTH'], atr['WIDTH'] = (str(n) for n in images[0][1].shape[-2:])
        names = [os.path.basename(p) for p, _ in images]
    if 'FILE_TYPE' not in atr:
        if fname.lower().endswith('.he5'):
            atr['FILE_TYPE'] = 'HDFEOS'
        elif 'timeseries' in names:
            atr['FILE_TYPE'] = 'timeseries'
        elif names:
            atr['FILE_TYPE'] = names[0]
    atr.setdefault('UNIT', _UNIT_OF_DSET.get(atr.get('FILE_TYPE'), '1'))
    return atr


def _rsc(fname):
    """The attributes of a ROI_PAC .rsc file."""
    atr = {}
    with open(fname) as fid:
        for line in fid:
            parts = line.split()
            if len(parts) >= 2:
                atr[parts[0]] = parts[1]
    if 'FILE_LENGTH' in atr and 'LENGTH' not in atr:
        atr['LENGTH'] = atr['FILE_LENGTH']
    return atr


def _roipac_attrs(fname):
    atr = _rsc(fname + '.rsc')
    ext = os.path.splitext(fname)[1].lower()
    atr.setdefault('FILE_TYPE', ext)
    atr.setdefault('UNIT', 'm' if ext in ('.dem', '.hgt') else '1')
    atr.setdefault('DATA_TYPE', 'int16' if ext == '.dem' else 'float32')
    return atr


def read_attribute(fname, datasetName=None):
    """The metadata of a file, in MintPy's terms (LENGTH, WIDTH, X_FIRST ...)."""
    if _licsbas(fname):
        from . import licsbas
        return licsbas.attributes(fname, datasetName)
    if mintpy_reads(fname):
        from mintpy.utils import readfile
        return dict(readfile.read_attribute(fname, datasetName=datasetName))
    if not os.path.isfile(fname) and not os.path.isdir(fname):
        raise FileNotFoundError(f'no such file: {fname}')
    if _is_hdf5(fname):
        atr = _h5_attrs(fname)
        if datasetName and atr.get('FILE_TYPE') == 'geometry':     # its own unit each
            atr['UNIT'] = _UNIT_OF_DSET.get(datasetName.split('-')[0], '1')
    elif os.path.isfile(fname + '.rsc'):
        atr = _roipac_attrs(fname)
    else:
        atr = _rasterio_read(fname, header_only=True)[1]
    atr['FILE_PATH'] = os.path.abspath(fname)
    return atr


def _rasterio_read(fname, box=None, header_only=False):
    try:
        import rasterio
        from rasterio.windows import Window
    except ImportError as err:
        raise ValueError(f'cannot read {fname}: GeoTIFF and the other GDAL rasters need '
                         f'rasterio (pip install rasterio), or MintPy ({MINTPY_HINT})') from err
    from .raster import _rasterio_attrs

    with rasterio.open(fname) as src:
        atr = _rasterio_attrs(src)
        if header_only:
            return None, atr
        if box is None:
            data = src.read(1)
        else:
            x0, y0, x1, y1 = box
            data = src.read(1, window=Window(x0, y0, x1 - x0, y1 - y0))
        nodata = src.nodata
    data = np.array(data)
    if nodata is not None and np.isfinite(nodata) and data.dtype.kind == 'f':
        data[data == np.asarray(nodata, dtype=data.dtype)] = np.nan
    return data, atr


def _pick_dataset(fid, datasetName, atr):
    """(dataset path, index into its first axis or None) for a slice name."""
    images = _h5_images(fid)
    if not images:
        raise ValueError('no 2D / 3D dataset in the file')
    if not datasetName:
        path, obj = images[0]
        return path, (obj.shape[0] - 1 if obj.ndim == 3 else None)
    base, _, date = datasetName.partition('-')
    for path, obj in images:
        leaf = os.path.basename(path)
        if leaf == datasetName and obj.ndim == 2:
            return path, None
        if leaf == base and obj.ndim == 3:
            if not date:
                return path, obj.shape[0] - 1
            dates = _h5_dates(fid, path, obj.shape[0])
            if date in dates:
                return path, dates.index(date)
            raise ValueError(f'{datasetName}: no date {date} in {leaf}')
    if re.fullmatch(r'\d{8}(_\d{8})?', datasetName):         # a bare date
        for path, obj in images:
            if obj.ndim == 3:
                dates = _h5_dates(fid, path, obj.shape[0])
                if datasetName in dates:
                    return path, dates.index(datasetName)
    raise ValueError(f'no dataset {datasetName!r}; there are '
                     f'{", ".join(os.path.basename(p) for p, _ in images)}')


def read(fname, datasetName=None, box=None, print_msg=False):
    """(data, metadata) of one 2D slice of a file; box = (x0, y0, x1, y1)."""
    if _licsbas(fname):
        from . import licsbas
        return licsbas.read(fname, datasetName=datasetName, box=box)
    if mintpy_reads(fname):
        from mintpy.utils import readfile
        data, atr = readfile.read(fname, datasetName=datasetName, box=box,
                                  print_msg=print_msg)
        return data, dict(atr)
    atr = read_attribute(fname, datasetName=datasetName)
    if _is_hdf5(fname):
        import h5py

        with h5py.File(fname, 'r') as fid:
            path, index = _pick_dataset(fid, datasetName, atr)
            dset = fid[path]
            rows = slice(box[1], box[3]) if box else slice(None)
            cols = slice(box[0], box[2]) if box else slice(None)
            data = dset[index, rows, cols] if index is not None else dset[rows, cols]
        leaf = os.path.basename(path)
        if atr.get('FILE_TYPE') in ('geometry', 'HDFEOS') or leaf in ('height',):
            atr['UNIT'] = _UNIT_OF_DSET.get(leaf, atr.get('UNIT', '1'))
    elif os.path.isfile(fname + '.rsc'):
        length, width = int(atr['LENGTH']), int(atr['WIDTH'])
        dtype = np.dtype(atr.get('DATA_TYPE', 'int16')).newbyteorder('<')
        data = np.memmap(fname, dtype=dtype, mode='r', shape=(length, width))
        data = np.array(data[box[1]:box[3], box[0]:box[2]] if box else data)
    else:
        data = _rasterio_read(fname, box=box)[0]
    return data, atr                     # the file's metadata, as MintPy's read()


def get_slice_list(fname):
    """The names of the 2D slices of a file: velocity, timeseries-20200101 ..."""
    if _licsbas(fname):
        from . import licsbas
        return licsbas.slices(fname)
    if mintpy_reads(fname):
        from mintpy.utils import readfile
        return list(readfile.get_slice_list(fname))
    if not _is_hdf5(fname):
        return [os.path.splitext(os.path.basename(fname))[0]]
    import h5py

    out = []
    with h5py.File(fname, 'r') as fid:
        for path, obj in _h5_images(fid):
            leaf = os.path.basename(path)
            if obj.ndim == 2:
                out.append(leaf)
            else:
                out += [f'{leaf}-{d}' for d in _h5_dates(fid, path, obj.shape[0])]
    return out


def get_file_list(patterns):
    """Existing files matching the names / wildcards, in the order given."""
    if have_mintpy():
        from mintpy.utils import utils as ut
        return ut.get_file_list(patterns)
    if isinstance(patterns, str):
        patterns = [patterns]
    out = []
    for pat in patterns or []:
        hits = sorted(glob.glob(os.path.expanduser(pat))) or (
            [pat] if os.path.isfile(pat) else [])
        out += [h for h in hits if h not in out]
    return out


##############################################################################
#                                coordinates                                 #
##############################################################################
class Coordinate:
    """The geocoded-grid part of MintPy's utils.coordinate: lat/lon <-> row/col."""

    def __init__(self, metadata, lookup_file=None):
        self.src_metadata = dict(metadata)
        self.lookup_file = lookup_file
        self.length, self.width = int(metadata['LENGTH']), int(metadata['WIDTH'])
        self.geocoded = 'Y_FIRST' in metadata
        if self.geocoded:
            self.lat0, self.lon0 = float(metadata['Y_FIRST']), float(metadata['X_FIRST'])
            self.lat_step = float(metadata['Y_STEP'])
            self.lon_step = float(metadata['X_STEP'])

    @staticmethod
    def _as_list(value):
        if value is None:
            return None
        return [float(v) for v in np.atleast_1d(np.asarray(value, dtype=np.float64))]

    def _degrees_to_grid(self, lat, lon):
        """Degrees given on a projected grid (UTM ...): into its own coordinates."""
        from .raster import crs_of
        from .swath import is_projected

        if (lat is None or lon is None or not is_projected(self.src_metadata)
                or np.max(np.abs(lat)) > 90 or np.max(np.abs(lon)) > 360):
            return lat, lon
        crs = crs_of(self.src_metadata)
        if crs is None:
            return lat, lon
        from pyproj import Transformer
        xs, ys = Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat)
        return list(np.atleast_1d(ys)), list(np.atleast_1d(xs))

    def lalo2yx(self, lat_in, lon_in):
        """(row, col) of lat / lon (or northing / easting), as MintPy: integers,
        a pair of lists for lists."""
        if not self.geocoded:
            raise ValueError('Input file is NOT geocoded.')
        scalar = np.ndim(lat_in) == 0 and np.ndim(lon_in) == 0
        lat, lon = self._as_list(lat_in), self._as_list(lon_in)
        lat, lon = self._degrees_to_grid(lat, lon)
        num = len(lat if lat is not None else lon)
        ys = [None if lat is None else int(np.floor((lat[i] - self.lat0) / self.lat_step + 0.01))
              for i in range(num)]
        xs = [None if lon is None else int(np.floor((lon[i] - self.lon0) / self.lon_step + 0.01))
              for i in range(num)]
        if scalar or num == 1 and not isinstance(lat_in, (list, tuple, np.ndarray)):
            return ys[0], xs[0]
        return ys, xs

    def yx2lalo(self, y_in, x_in):
        """lat / lon of pixel centres."""
        if not self.geocoded:
            raise ValueError('Input file is NOT geocoded.')
        lat = self.lat0 + (np.asarray(y_in, dtype=np.float64) + 0.5) * self.lat_step
        lon = self.lon0 + (np.asarray(x_in, dtype=np.float64) + 0.5) * self.lon_step
        return (float(lat), float(lon)) if np.ndim(lat) == 0 else (lat, lon)

    def box_pixel2geo(self, pixel_box):
        """(W, N, E, S) of a pixel box (x0, y0, x1, y1) - None on a radar grid."""
        if not self.geocoded or pixel_box is None:
            return None
        x0, y0, x1, y1 = pixel_box
        return (self.lon0 + x0 * self.lon_step, self.lat0 + y0 * self.lat_step,
                self.lon0 + x1 * self.lon_step, self.lat0 + y1 * self.lat_step)

    def box_geo2pixel(self, geo_box):
        """(x0, y0, x1, y1) of a geo box (W, N, E, S)."""
        y0, x0 = self.lalo2yx(geo_box[1], geo_box[0])
        y1, x1 = self.lalo2yx(geo_box[3], geo_box[2])
        return (x0, y0, x1, y1)

    def check_box_within_data_coverage(self, pixel_box, print_msg=True):
        """The box, clipped to the data."""
        x0, y0, x1, y1 = (int(v) for v in pixel_box)
        if y0 >= self.length or y1 <= 0 or x0 >= self.width or x1 <= 0:
            raise ValueError(f'the subset {pixel_box} is out of the data coverage '
                             f'(0, 0, {self.width}, {self.length})')
        out = (max(0, x0), max(0, y0), min(self.width, x1), min(self.length, y1))
        if out != (x0, y0, x1, y1) and print_msg:
            print(f'subset clipped to the data coverage: {out}')
        return out


def coordinate(metadata, lookup_file=None):
    """MintPy's utils.coordinate, or the stand-in."""
    if have_mintpy():
        from mintpy.utils import utils as ut
        return ut.coordinate(metadata, lookup_file=lookup_file)
    return Coordinate(metadata, lookup_file=lookup_file)


##############################################################################
#                           units, mask, figure                              #
##############################################################################
def auto_disp_unit(metadata):
    """The display unit view.py picks: cm (cm/year ...) for MintPy's length
    data, the file's unit otherwise."""
    ftype = str(metadata.get('FILE_TYPE', '')).replace('.', '')
    unit = str(metadata.get('UNIT', '1')).lower()
    if (ftype in ('timeseries', 'giantTimeseries', 'velocity', 'HDFEOS')
            and unit.split('/')[0].endswith('m')
            and not str(metadata.get('DATA_TYPE', 'float32')).startswith('complex')):
        return 'cm'
    return metadata.get('UNIT', '1')


def scale_data2disp_unit(data=None, metadata=None, disp_unit=None):
    """Scale data from the file's unit to the display unit, as MintPy does.
    Returns (data, display unit, scale)."""
    if have_mintpy():
        from mintpy.utils import plot as pp
        return pp.scale_data2disp_unit(data, metadata=metadata, disp_unit=disp_unit)
    metadata = dict(metadata or {'UNIT': 'm'})
    data_unit = str(metadata.get('UNIT', '1')).lower().split('/')
    disp = str(disp_unit or metadata.get('UNIT', '1')).split('/')
    if [d.lower() for d in disp] == data_unit:
        return data, metadata.get('UNIT', '1'), 1.0
    scale = 1.0
    head, want = data_unit[0], disp[0].lower()
    if head.endswith('m') and head in _LENGTH_UNITS:
        if want in _LENGTH_UNITS:
            scale = _LENGTH_UNITS[head] / _LENGTH_UNITS[want]
        elif want in ('radian', 'radians', 'rad', 'r'):
            scale = _LENGTH_UNITS[head] * -4.0 * np.pi / float(metadata['WAVELENGTH'])
            disp[0] = 'radian'
        else:
            print(f'Unrecognized display unit: {disp[0]}')
            disp = [metadata.get('UNIT', '1')]
    elif head == 'radian':
        if want in _LENGTH_UNITS:
            scale = -float(metadata['WAVELENGTH']) / (4 * np.pi) / _LENGTH_UNITS[want]
        elif want not in ('radian', 'radians', 'rad', 'r'):
            print(f'Unrecognized display unit: {disp[0]}')
    elif head != '1':
        disp = [metadata.get('UNIT', '1')]
    if len(data_unit) == 2:
        tail = disp[1].lower() if len(disp) > 1 else 'year'
        if tail in ('y', 'yr', 'year'):
            tail = 'year'
        elif tail in ('m', 'mon', 'month'):
            tail, scale = 'mon', scale / 12.0
        elif tail in ('d', 'day'):
            tail, scale = 'day', scale / 365.25
        unit = f'{disp[0]}/{tail}'
    else:
        unit = disp[0]
    if data is not None and scale != 1.0:
        data = np.asarray(data, dtype=np.float32 if np.asarray(data).dtype.kind != 'f'
                          else np.asarray(data).dtype) * scale
    return data, unit, scale


def read_mask(fname, mask_file=None, datasetName=None, box=None, print_msg=False):
    """The mask of a data file, as view.py finds it: maskTempCoh.h5 (or
    maskResInv.h5) beside a velocity / time series, -m FILE, or none (-m no).
    Returns (bool mask or None, mask file)."""
    if _licsbas(fname) and not mask_file and 'msk' not in os.path.basename(fname):
        from .licsbas import mask_for
        mask_file = mask_for(fname)              # LiCSBAS's results/mask, as LiCSBAS
        if not mask_file:
            return None, None
    if mintpy_reads(fname) and (not mask_file or mintpy_reads(mask_file)):
        from mintpy.utils import plot as pp
        return pp.read_mask(fname, mask_file=mask_file, datasetName=datasetName, box=box,
                            print_msg=print_msg)
    atr = read_attribute(fname)
    if (not mask_file and atr.get('FILE_TYPE') in ('velocity', 'timeseries')
            and 'msk' not in os.path.basename(fname)):
        prefix = 'geo_' if os.path.basename(fname).startswith('geo_') else ''
        for name in ('maskTempCoh.h5', 'maskResInv.h5'):
            cand = os.path.join(os.path.dirname(fname), prefix + name)
            if os.path.isfile(cand):
                mask_file = cand
                break
    if not mask_file or not os.path.isfile(str(mask_file)):
        if atr.get('FILE_TYPE') == 'HDFEOS':
            try:
                return read(fname, datasetName='mask', box=box)[0] != 0, fname
            except ValueError:
                pass
        return None, None
    try:
        atr_m = read_attribute(mask_file)
        if (atr_m['LENGTH'], atr_m['WIDTH']) != (atr['LENGTH'], atr['WIDTH']):
            if print_msg:
                print(f'WARNING: the mask {mask_file} has another size than the data '
                      '--> no mask')
            return None, None
        data = read(mask_file, box=box)[0]
    except Exception as err:     # noqa: BLE001 - the mask is a nicety
        if print_msg:
            print(f'WARNING: cannot read the mask {mask_file}: {err}')
        return None, None
    if print_msg:
        print(f'read mask from file: {os.path.basename(mask_file)}')
    data = np.asarray(data)
    return (data != 0) & (np.isfinite(data) if data.dtype.kind == 'f' else True), mask_file


def auto_figure_size(ds_shape, disp_cbar=False, cbar_ratio=0.25, print_msg=False):
    """A figure size for a 2D image [length, width], as MintPy picks it."""
    if have_mintpy():
        from mintpy.utils import plot as pp
        return pp.auto_figure_size(ds_shape, disp_cbar=disp_cbar, print_msg=print_msg)
    shape = [float(ds_shape[1]), float(ds_shape[0])]
    shape[0] *= 1 + cbar_ratio if disp_cbar else 1
    scale = min(6.0 / min(shape), 10.0 / max(shape), 8.0 / shape[1])
    return [float(f'{x * scale:.1f}') for x in shape]


def auto_flip_direction(metadata):
    """(flip left-right, flip up-down) of a radar grid, as MintPy shows it:
    ascending up-down, descending left-right."""
    if 'Y_FIRST' in metadata or 'ORBIT_DIRECTION' not in metadata:
        return False, False
    if str(metadata['ORBIT_DIRECTION']).lower().startswith('a'):
        return False, True
    return True, False

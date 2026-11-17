"""Time series for `viewinsar tseries`.

LiCSBAS cum.h5 / cum_filt.h5 and MintPy timeseries.h5 / HDF-EOS5 (.he5) are
read alike: displacement in mm, the dates, the grid's geocoding as MintPy
attributes, the reference area and date, the velocity and the extra layers
(LiCSBAS results/ noise indices; MintPy mask, temporal coherence, height).

Then what the viewer plots, as LiCSBAS_plot_ts.py does it: every pixel
relative to the reference area, date by date, and to its own value at the
reference date; for a circle of pixels the median of each date with the 25th
and 75th percentiles; and the four LiCSBAS fit models.

Numpy only, plus h5py when a file is opened.
"""

import datetime as dt
import os
import re
import warnings

import numpy as np

from .swath import get_pixel_size, is_projected

__all__ = ['MODELS', 'TO_UNIT', 'TimeSeries', 'parse_range', 'parse_geo_range',
           'circle_pixels_rc', 'relative_series', 'series_stats', 'model_matrix',
           'model_curve', 'linear_velocity', 'auto_limits', 'find_time_series']

MODELS = ('Linear', 'Annual+L', 'Quad', 'Annual+Q')
TO_UNIT = {'mm': 1.0, 'cm': 0.1, 'm': 0.001}            # from mm

# LiCSBAS results/: (name shown, file name, unit) - as LiCSBAS_plot_ts.py
LICSBAS_LAYERS = (('mask', 'mask', ''), ('coh_avg', 'coh_avg', ''), ('n_unw', 'n_unw', ''),
                  ('vstd', 'vstd', 'mm/yr'), ('maxTlen', 'maxTlen', 'yr'),
                  ('n_gap', 'n_gap', ''), ('stc', 'stc', 'mm'),
                  ('n_ifg_noloop', 'n_ifg_noloop', ''), ('n_loop_err', 'n_loop_err', ''),
                  ('resid', 'resid_rms', 'mm'), ('mli', 'slc.mli', 'log10'),
                  ('hgt', 'hgt', 'm'))
HDFEOS = 'HDFEOS/GRIDS/timeseries'
DEFAULT_FILES = ('cum_filt.h5', 'cum.h5', 'timeseries.h5')
_UNIT_TO_MM = {'m': 1000.0, 'cm': 10.0, 'mm': 1.0}


##############################################################################
#                               small helpers                                #
##############################################################################
def _text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode('utf-8', 'replace')
    if isinstance(value, np.ndarray) and value.size == 1:
        return _text(value.ravel()[0])
    return str(value)


def _attrs(h5obj):
    """HDF5 attributes as MintPy keeps them: every value a string."""
    out = {}
    for key, value in h5obj.attrs.items():
        out[key] = _text(value) if np.ndim(value) == 0 or isinstance(value, bytes) else value
    return out


def _to_mm(unit, default='m'):
    unit = (unit or default).lower().split('/')[0].strip()
    return _UNIT_TO_MM.get(unit, _UNIT_TO_MM[default])


def find_time_series(folder='.'):
    """The time-series file to open when none is named, or None."""
    for name in DEFAULT_FILES:
        path = name if folder in ('', '.') else os.path.join(folder, name)
        if os.path.isfile(path):
            return path
    return None


def parse_range(text, width, length):
    """LiCSBAS 'x1:x2/y1:y2' -> (x1, x2, y1, y2), python slices; 0 for x2 or y2
    means to the end (0:0/0:0 is the whole grid)."""
    match = re.fullmatch(r'\s*(\d+):(\d+)/(\d+):(\d+)\s*', str(text))
    if not match:
        raise ValueError(f'not a range x1:x2/y1:y2: {text!r}')
    x1, x2, y1, y2 = (int(v) for v in match.groups())
    x2, y2 = x2 or width, y2 or length
    if not (0 <= x1 < x2 <= width and 0 <= y1 < y2 <= length):
        raise ValueError(f'range {text} does not fit the {width} x {length} (width x length) grid')
    return x1, x2, y1, y2


def parse_geo_range(text):
    """LiCSBAS 'lon1/lon2/lat1/lat2' -> four floats."""
    parts = re.split(r'/', str(text).strip())
    if len(parts) != 4:
        raise ValueError(f'not a range lon1/lon2/lat1/lat2: {text!r}')
    return tuple(float(p) for p in parts)


def auto_limits(data, pct=99.0):
    """The (100 - pct, pct) percentiles of the finite values, or (-1, 1)."""
    good = np.asarray(data)[np.isfinite(data)]
    if not good.size:
        return -1.0, 1.0
    lo, hi = np.percentile(good, [100.0 - pct, pct])
    if hi <= lo:
        lo, hi = lo - 1.0, hi + 1.0
    return float(lo), float(hi)


##############################################################################
#                            the time-series file                            #
##############################################################################
class TimeSeries:
    """One LiCSBAS or MintPy time-series file.

    cube      - (n_date, length, width) HDF5 dataset; `to_mm` scales it to mm
    dates     - list of 'YYYYMMDD'; datetimes; years (since the first date)
    atr       - MintPy-style attributes (geocoding, pixel size, orbit)
    ref_area  - (x1, x2, y1, y2) of the file (LiCSBAS refarea / MintPy REF_X, REF_Y)
    ref_date  - the file's reference date (MintPy REF_DATE) or None
    vel       - 2D velocity in mm/yr, or None until velocity() is called
    keep      - 2D bool from the file's mask (True = use the pixel), or None
    layers    - {name: (2D array, unit)}: LiCSBAS noise indices, or MintPy mask,
                temporal coherence and height
    gap       - LiCSBAS gap dataset (n_date - 1, length, width), or None
    label     - the legend label (LiCSBAS: the filter of cum_filt.h5)
    """

    def __init__(self, fname, mask_file=None, vel_file=None, results_dir=None,
                 verbose=True, scale=None):
        import h5py

        self.fname, self.verbose = fname, verbose
        self.h5 = h5py.File(fname, 'r')
        self.keep, self.layers, self.gap, self.vel = None, {}, None, None
        self.ref_date, self.mask_source, self.vel_source = None, None, None
        keys = set(self.h5.keys())
        if 'imdates' in keys and ({'cum', 'cum_model'} & keys):
            self._open_licsbas(results_dir)
        elif 'timeseries' in keys and 'date' in keys:
            self._open_mintpy(self.h5['timeseries'], self.h5['date'], self.h5)
        elif 'HDFEOS' in keys and f'{HDFEOS}/observation/displacement' in self.h5:
            obs = self.h5[f'{HDFEOS}/observation']
            self._open_mintpy(obs['displacement'], obs['date'], self.h5, hdfeos=True)
        else:
            self.h5.close()
            raise ValueError(f'{fname}: neither a LiCSBAS cum*.h5 nor a MintPy time series')
        self.n_date, self.length, self.width = self.cube.shape
        self.atr.update(LENGTH=str(self.length), WIDTH=str(self.width))
        self.datetimes = [dt.datetime.strptime(d, '%Y%m%d') for d in self.dates]
        self.years = np.array([(d - self.datetimes[0]).days for d in self.datetimes]) / 365.25
        if mask_file:
            self.set_mask_file(mask_file)
        if vel_file:
            self.vel = self._read_velocity(vel_file)
            self.vel_source = vel_file
        self._scale(scale)

    def _scale(self, scale):
        """--data-scale (or a rule of viewinsar.ini): every value x the factor."""
        from .settings import config_scale, wavelength_note

        if scale is not None:
            self.scale, why = float(scale), '--data-scale'
        else:
            self.scale, why = config_scale(self.fname, self.atr)
        note = wavelength_note(self.fname, self.atr)
        if note:
            print(note)
        if self.scale != 1.0:
            self.to_mm *= self.scale
            if self.vel is not None:
                self.vel = self.vel * self.scale
            print(f'{os.path.basename(self.fname)}: its values x {self.scale:g} ({why}) - the '
                  'map and every time series')

    def vprint(self, *args):
        if self.verbose:
            print(*args)

    ##---------- LiCSBAS
    def _open_licsbas(self, results_dir):
        h5 = self.h5
        self.kind, self.to_mm = 'LiCSBAS', 1.0
        self.cube = h5['cum'] if 'cum' in h5 else h5['cum_model']
        self.dates = [str(int(d)) for d in np.ravel(h5['imdates'][()])]
        _, length, width = self.cube.shape
        atr = {'FILE_TYPE': 'timeseries', 'UNIT': 'mm', 'PROCESSOR': 'LiCSBAS'}
        if all(k in h5 for k in ('corner_lat', 'corner_lon', 'post_lat', 'post_lon')):
            # LiCSBAS gives the centre of the first pixel; MintPy its outer corner
            lat1, lon1 = float(h5['corner_lat'][()]), float(h5['corner_lon'][()])
            dlat, dlon = float(h5['post_lat'][()]), float(h5['post_lon'][()])
            atr.update(Y_FIRST=repr(lat1 - dlat / 2), X_FIRST=repr(lon1 - dlon / 2),
                       Y_STEP=repr(dlat), X_STEP=repr(dlon), Y_UNIT='degrees',
                       X_UNIT='degrees')
        self.atr = atr
        self.ref_area = (parse_range(_text(h5['refarea'][()]), width, length)
                         if 'refarea' in h5 else (0, width, 0, length))
        if 'vel' in h5:
            self.vel = np.asarray(h5['vel'], dtype=np.float64)
            self.vel_source = 'vel'
        self.gap = h5['gap'] if 'gap' in h5 else None
        self.label = '1: No filter'
        if 'filtwidth_km' in h5 and 'filtwidth_yr' in h5:
            deramp = ''
            if 'deramp_flag' in h5:
                flag = _text(h5['deramp_flag'][()])
                deramp = f', drmp={flag}' if flag and flag != '[]' else ''
            if 'hgt_linear_flag' in h5 and int(np.ravel(h5['hgt_linear_flag'][()])[0]) == 1:
                deramp += ', hgt-linear'
            km, yr = float(h5['filtwidth_km'][()]), float(h5['filtwidth_yr'][()])
            self.label = f'1: s={km:.1f}km, t={yr:.2f}yr ({int(round(yr * 365.25))}d){deramp}'

        folder = results_dir or os.path.join(os.path.dirname(os.path.abspath(self.fname)),
                                             'results')
        for name, fname, unit in LICSBAS_LAYERS:
            path = os.path.join(folder, fname)
            if not os.path.isfile(path) or os.path.getsize(path) != length * width * 4:
                continue
            data = np.fromfile(path, dtype='<f4').reshape(length, width).astype(np.float64)
            if name == 'mli':
                with np.errstate(divide='ignore', invalid='ignore'):
                    data = np.log10(data)
            if name == 'mask':
                self.keep = np.isfinite(data) & (data != 0)
                self.mask_source = path
            self.layers[name] = (data, unit)

    ##---------- MintPy (timeseries.h5, or the HDF-EOS5 of save_hdfeos5.py)
    def _open_mintpy(self, cube, dates, h5, hdfeos=False):
        self.kind = 'MintPy'
        self.cube = cube
        self.dates = [_text(d) for d in np.ravel(dates[()])]
        self.atr = _attrs(h5)
        self.atr.setdefault('FILE_TYPE', 'timeseries')
        self.to_mm = _to_mm(self.atr.get('UNIT', 'm'))
        _, length, width = cube.shape
        try:
            ref_y, ref_x = int(self.atr['REF_Y']), int(self.atr['REF_X'])
            self.ref_area = (ref_x, ref_x + 1, ref_y, ref_y + 1)
        except (KeyError, ValueError):
            self.ref_area = (0, width, 0, length)
        if self.atr.get('REF_DATE') in self.dates:
            self.ref_date = self.atr['REF_DATE']
        self.label = f'1: {os.path.basename(self.fname)}'

        if hdfeos:
            grid = h5[HDFEOS]
            if 'quality/mask' in grid:
                self.keep = np.asarray(grid['quality/mask']).astype(bool)
                self.mask_source = f'{os.path.basename(self.fname)}: quality/mask'
            extra = [('temporalCoherence', 'quality/temporalCoherence', ''),
                     ('hgt', 'geometry/height', 'm')]
            for name, path, unit in extra:
                if path in grid:
                    self.layers[name] = (np.asarray(grid[path], dtype=np.float64), unit)
        else:
            found = self._companion(['maskTempCoh.h5'], 'mask')
            if found:
                self.keep, self.mask_source = found[0].astype(bool), found[1]
            found = self._companion(['temporalCoherence.h5'], 'temporalCoherence')
            if found:
                self.layers['temporalCoherence'] = (found[0].astype(np.float64), '')
            found = self._companion(['geometryGeo.h5', 'geometryRadar.h5'], 'height')
            if found:
                self.layers['hgt'] = (found[0].astype(np.float64), 'm')
            found = self._companion(['velocity.h5'], 'velocity', with_attrs=True)
            if found:
                data, path, atr = found
                self.vel = data.astype(np.float64) * _to_mm(atr.get('UNIT', 'm/year'))
                self.vel_source = path
        if self.keep is not None:
            self.layers = {'mask': (self.keep.astype(np.float64), ''), **self.layers}

    def _companion(self, names, dset, with_attrs=False):
        """A MintPy file of the same grid next to the time series: in its
        folder, inputs/, ../inputs/ or geo/, with or without the geo_ prefix."""
        import h5py

        folder = os.path.dirname(os.path.abspath(self.fname))
        prefix = 'geo_' if os.path.basename(self.fname).startswith('geo_') else ''
        for sub in ('', 'inputs', os.path.join('..', 'inputs'), 'geo'):
            for name in names:
                for cand in dict.fromkeys([prefix + name, name]):
                    path = os.path.normpath(os.path.join(folder, sub, cand))
                    if not os.path.isfile(path):
                        continue
                    try:
                        with h5py.File(path, 'r') as fid:
                            if dset not in fid or fid[dset].shape != self.cube.shape[1:]:
                                continue
                            data = np.asarray(fid[dset])
                            atr = _attrs(fid)
                    except OSError:
                        continue
                    self.vprint(f'{dset}: {path}')
                    return (data, path, atr) if with_attrs else (data, path)
        return None

    def _read_2d(self, fname, prefer):
        import h5py

        with h5py.File(fname, 'r') as fid:
            names = [prefer] if prefer in fid else [k for k in fid.keys()
                                                    if getattr(fid[k], 'ndim', 0) == 2]
            if not names:
                raise ValueError(f'{fname}: no 2D dataset')
            data = np.asarray(fid[names[0]])
            atr = _attrs(fid)
        if data.shape != (self.length, self.width):
            raise ValueError(f'{fname} is {data.shape[0]} x {data.shape[1]}, the time series '
                             f'{self.length} x {self.width}')
        return data, atr

    def set_mask_file(self, fname):
        """--mask: this mask instead of the file's own."""
        self.keep = self._read_mask(fname)
        self.mask_source = fname
        self.layers['mask'] = (self.keep.astype(np.float64), '')

    def _read_mask(self, fname):
        if not fname.endswith(('.h5', '.he5')):
            data = np.fromfile(fname, dtype='<f4')
            if data.size != self.length * self.width:
                raise ValueError(f'{fname}: not a {self.length} x {self.width} float32 mask')
            data = data.reshape(self.length, self.width)
        else:
            data = self._read_2d(fname, 'mask')[0]
        return np.isfinite(data) & (data != 0)

    def _read_velocity(self, fname):
        data, atr = self._read_2d(fname, 'velocity' if self.kind == 'MintPy' else 'vel')
        default = 'mm/year' if self.kind == 'LiCSBAS' else 'm/year'
        return data.astype(np.float64) * _to_mm(atr.get('UNIT', default),
                                                default.split('/')[0])

    ##---------- data
    def velocity(self):
        """The velocity in mm/yr: the file's (LiCSBAS vel, MintPy velocity.h5),
        or a linear fit of the time series, worked out once."""
        if self.vel is None:
            self.vprint(f'velocity: a linear fit of {os.path.basename(self.fname)} '
                        '(no velocity file found)')
            self.vel = linear_velocity(self.cube, self.years, self.to_mm)
            self.vel_source = 'linear fit of the time series'
        return self.vel

    def slice(self, index):
        """One date, in mm."""
        return np.asarray(self.cube[index], dtype=np.float64) * self.to_mm

    def pixels(self, rows, cols):
        """(n_date, n_pixel) mm at the given pixels, read through their bounding box."""
        rows, cols = np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)
        if not rows.size:
            return np.empty((self.n_date, 0))
        r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
        sub = np.asarray(self.cube[:, r0:r1, c0:c1], dtype=np.float64) * self.to_mm
        return sub[:, rows - r0, cols - c0]

    def area_mean(self, area, keep=None):
        """Mean of the area (x1, x2, y1, y2) at every date, in mm (NaN-aware)."""
        x1, x2, y1, y2 = area
        sub = np.asarray(self.cube[:, y1:y2, x1:x2], dtype=np.float64) * self.to_mm
        if keep is not None:
            sub = np.where(keep[y1:y2, x1:x2][None], sub, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return np.nanmean(sub, axis=(1, 2))

    def gap_fraction(self, rows, cols):
        """For each interval between dates, the share of the pixels with a gap
        in their interferogram network (LiCSBAS), or None."""
        if self.gap is None or not np.size(rows):
            return None
        rows, cols = np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)
        r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
        sub = np.asarray(self.gap[:, r0:r1, c0:c1])[:, rows - r0, cols - c0]
        return (sub == 1).mean(axis=1)

    ##---------- coordinates
    @property
    def geocoded(self):
        return 'Y_FIRST' in self.atr

    @property
    def projected(self):
        return self.geocoded and is_projected(self.atr)

    def yx2lalo(self, row, col):
        """(lat, lon) - northing, easting on a projected grid - of a pixel centre."""
        lat = float(self.atr['Y_FIRST']) + (row + 0.5) * float(self.atr['Y_STEP'])
        lon = float(self.atr['X_FIRST']) + (col + 0.5) * float(self.atr['X_STEP'])
        return lat, lon

    def lalo2yx(self, lat, lon):
        """The pixel (row, col) holding a lat/lon (northing/easting) point."""
        row = int(round((lat - float(self.atr['Y_FIRST'])) / float(self.atr['Y_STEP']) - 0.5))
        col = int(round((lon - float(self.atr['X_FIRST'])) / float(self.atr['X_STEP']) - 0.5))
        return row, col

    def pixel_size_km(self, row):
        """(x_step, y_step) of the grid in km at a given row."""
        lat_c = None
        if self.geocoded and not self.projected:
            lat_c = self.yx2lalo(row, 0)[0]
        x_step, y_step, unit = get_pixel_size(self.atr, lat_c=lat_c)
        if unit != 'm':
            raise ValueError('the pixel size of this grid is unknown, so a radius in km '
                             'cannot be turned into pixels')
        return x_step / 1000.0, y_step / 1000.0

    def close(self):
        try:
            self.h5.close()
        except Exception:        # noqa: BLE001
            pass


##############################################################################
#                       what the viewer works out                            #
##############################################################################
def circle_pixels_rc(ts, row, col, radius_km):
    """Rows and cols of the pixels whose centre is within radius_km of the
    centre of pixel (row, col) - on the ground, from the grid's pixel size."""
    x_km, y_km = ts.pixel_size_km(row)
    d_row, d_col = int(np.floor(radius_km / y_km)), int(np.floor(radius_km / x_km))
    rows = np.arange(max(row - d_row, 0), min(row + d_row, ts.length - 1) + 1)
    cols = np.arange(max(col - d_col, 0), min(col + d_col, ts.width - 1) + 1)
    rr, cc = np.meshgrid(rows, cols, indexing='ij')
    inside = ((rr - row) * y_km) ** 2 + ((cc - col) * x_km) ** 2 <= radius_km ** 2 * (1 + 1e-12)
    return rr[inside], cc[inside]


def relative_series(values, ref, ix_ref):
    """LiCSBAS's displacement of each pixel: minus the reference area at every
    date, then minus its own value at the reference date (so it starts at 0).

    values - (n_date, n_pixel); ref - (n_date,) the reference area's mean
    """
    out = np.asarray(values, dtype=np.float64) - np.asarray(ref, dtype=np.float64)[:, None]
    return out - out[ix_ref][None, :]


def series_stats(values):
    """Per date (row): median, 25th and 75th percentiles and count of the finite
    values of (n_date, n_pixel). One pixel: all three are its value."""
    values = np.asarray(values, dtype=np.float64)
    count = np.isfinite(values).sum(axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)          # an all-NaN date
        median = np.nanmedian(values, axis=1) if values.shape[1] else np.full(len(values), np.nan)
        if values.shape[1]:
            q25, q75 = np.nanpercentile(values, [25, 75], axis=1)
        else:
            q25 = q75 = np.full(len(values), np.nan)
    return dict(median=median, q25=q25, q75=q75, count=count)


def model_matrix(years, model):
    """LiCSBAS's fit models: 0 Linear, 1 Annual+L, 2 Quad, 3 Annual+Q."""
    t = np.asarray(years, dtype=np.float64)
    cols = [np.ones_like(t), t]
    if model in (2, 3):
        cols.append(t ** 2)
    if model in (1, 3):
        cols += [np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)]
    return np.column_stack(cols)


def model_curve(years, values, years_out, model):
    """Least-squares fit of one model to a series (NaN dropped), evaluated at
    years_out; NaN when there are too few dates."""
    years, values = np.asarray(years, dtype=np.float64), np.asarray(values, dtype=np.float64)
    good = np.isfinite(values)
    design = model_matrix(years[good], model)
    if good.sum() < design.shape[1] + 1:
        return np.full(len(years_out), np.nan)
    coef = np.linalg.lstsq(design, values[good], rcond=None)[0]
    return model_matrix(years_out, model) @ coef


def linear_velocity(cube, years, to_mm=1.0, max_values=20_000_000):
    """NaN-aware least-squares slope of every pixel of a (n_date, length, width)
    cube, in mm/yr - read a band of rows at a time."""
    n_date, length, width = cube.shape
    years = np.asarray(years, dtype=np.float64)
    out = np.full((length, width), np.nan, dtype=np.float64)
    step = max(1, int(max_values // max(n_date * width, 1)))
    for r0 in range(0, length, step):
        data = np.asarray(cube[:, r0:r0 + step], dtype=np.float64) * to_mm
        good = np.isfinite(data)
        t = np.where(good, years[:, None, None], 0.0)
        y = np.where(good, data, 0.0)
        n = good.sum(axis=0)
        st, sy = t.sum(axis=0), y.sum(axis=0)
        den = n * (t * t).sum(axis=0) - st * st
        with np.errstate(divide='ignore', invalid='ignore'):
            slope = (n * (t * y).sum(axis=0) - st * sy) / den
        slope[(n < 2) | ~(den > 0)] = np.nan
        out[r0:r0 + step] = slope
    return out

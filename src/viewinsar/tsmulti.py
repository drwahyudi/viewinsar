"""`viewinsar tseries` with several files: one map, and at a double click the
time series of every file that covers the point.

The files - adjacent frames of a track, overlapping tracks, or versions of one
time series - are drawn on one map, each where it lies, its outline and its
series in a colour of its own. They must be geocoded in one CRS (LiCSBAS
GEOC, MintPy geo_*), or radar-coded on one grid (the MintPy time series of
one stack).

Each file keeps its own dates and its own reference area - the file's, until
a right drag (or --ref-geo, --ref-lalo) gives one that it covers: drag one in
the overlap of two frames and both are relative to the same ground. The
reference date is common: --ref-date, else the latest of the files' first
dates (MintPy: REF_DATE), so that every series starts together; a file
without that date takes its nearest one.

    left double click   the time series of every file there (circles with --radius)
    right drag          a new reference area, for every file with data in it
    "e"                 export the time series (.csv), its figure and the map
"""

import datetime as dt
import os
import warnings

import numpy as np
from matplotlib import dates as mdates
from matplotlib import patches, patheffects
from matplotlib import pyplot as plt
from matplotlib.widgets import CheckButtons, RadioButtons, RectangleSelector, Slider
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .tsdata import (
    MODELS,
    TO_UNIT,
    TimeSeries,
    auto_limits,
    circle_pixels_rc,
    model_curve,
    parse_geo_range,
    parse_range,
    relative_series,
    series_stats,
)
from .tsviewer import (
    FIT_DEFAULT,
    GAP_SHARE,
    HELP_TEXT,
    _fmt,
    get_cmap,
    map_window_size,
    place_map_window,
)
from .viewer import _new_figure, file_labels, per_file

__all__ = ['FRAME_COLORS', 'Frame', 'MultiTimeSeriesViewer']

# one colour per file, on the map (outline, reference area) and in the plot;
# the first two are LiCSBAS's blue and red
FRAME_COLORS = ('tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple',
                'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan', 'tab:gray')
HALO = [patheffects.withStroke(linewidth=2.8, foreground='w')]
MAX_VERTICES = 400                 # of a file's outline, per side


def _short(text, size=18):
    return text if len(text) <= size else text[:size - 1] + '…'


def _finite(arrays, most=2_000_000):
    """The finite values of some arrays, for colour limits: at most about
    `most` of them, evenly thinned - percentiles of millions are slow."""
    values = np.concatenate([np.asarray(a)[np.isfinite(a)] for a in arrays])
    return values[::max(1, values.size // most)]


def _decimal_year(day):
    start = dt.datetime(day.year, 1, 1)
    return day.year + (day - start).days / (dt.datetime(day.year + 1, 1, 1) - start).days


##############################################################################
class Frame:
    """One input file on the map: its time series, where its grid lies, its
    reference area and date."""

    def __init__(self, viewer, ts, index, label, fname):
        self.v, self.ts, self.index, self.label, self.file = viewer, ts, index, label, fname
        self.tag = str(index + 1)
        self.color = FRAME_COLORS[index % len(FRAME_COLORS)]
        atr = ts.atr
        if ts.geocoded:
            self.x0, self.dx = float(atr['X_FIRST']), float(atr['X_STEP'])
            self.y0, self.dy = float(atr['Y_FIRST']), float(atr['Y_STEP'])
        else:                            # radar-coded: the map is in its pixels
            self.x0, self.dx, self.y0, self.dy = -0.5, 1.0, -0.5, 1.0
        self.ref = ts.ref_area           # (x1, x2, y1, y2), pixels of this file
        self.ix_ref = 0
        self.ref_slice = self.base_keep = self.vel = self.ref1 = None
        self.image = self.outline = self.ref_line = None
        self._last_cum = None
        self.ords = np.array([d.toordinal() for d in ts.datetimes])
        steps = np.diff(self.ords)
        self.half_step = 0.5 * float(np.median(steps)) if steps.size else 6.0

    ##---------- the grid on the map
    def xy(self, row, col):
        """Map coordinates of a pixel centre (or of a fractional pixel)."""
        return self.x0 + (col + 0.5) * self.dx, self.y0 + (row + 0.5) * self.dy

    def rc_float(self, x, y):
        """The pixel coordinates of a map point: pixel centres at whole numbers."""
        return (y - self.y0) / self.dy - 0.5, (x - self.x0) / self.dx - 0.5

    def rc(self, x, y):
        """The pixel (row, col) holding map point (x, y), or None off the grid."""
        r, c = self.rc_float(x, y)
        row, col = int(np.floor(r + 0.5)), int(np.floor(c + 0.5))
        if 0 <= row < self.ts.length and 0 <= col < self.ts.width:
            return row, col
        return None

    @property
    def extent(self):
        """(left, right, bottom, top) of the grid, for imshow (origin upper)."""
        if not self.ts.geocoded:         # rows top-down, as any image
            return -0.5, self.ts.width - 0.5, self.ts.length - 0.5, -0.5
        xa, xb = self.x0, self.x0 + self.ts.width * self.dx
        ya, yb = self.y0, self.y0 + self.ts.length * self.dy
        return min(xa, xb), max(xa, xb), min(ya, yb), max(ya, yb)

    def display(self, data):
        """A grid as imshow draws it in the extent: north up, east right."""
        if self.ts.geocoded:
            if self.dy > 0:
                data = data[::-1]
            if self.dx < 0:
                data = data[:, ::-1]
        return data

    def box_xy(self, box):
        """The outline of a pixel box (x1, x2, y1, y2) on the map."""
        x1, x2, y1, y2 = box
        xa, ya = self.x0 + x1 * self.dx, self.y0 + y1 * self.dy
        xb, yb = self.x0 + x2 * self.dx, self.y0 + y2 * self.dy
        return [xa, xb, xb, xa, xa], [ya, ya, yb, yb, ya]

    def box_of(self, xs, ys):
        """The pixel box (x1, x2, y1, y2) of the pixels whose centres lie in a
        map rectangle - the one pixel at its centre if it holds none - or
        None when it misses the grid."""
        (ra, ca), (rb, cb) = self.rc_float(xs[0], ys[0]), self.rc_float(xs[1], ys[1])
        (r_lo, r_hi), (c_lo, c_hi) = sorted((ra, rb)), sorted((ca, cb))
        x1 = max(int(np.ceil(c_lo - 1e-9)), 0)
        x2 = min(int(np.floor(c_hi + 1e-9)), self.ts.width - 1) + 1
        y1 = max(int(np.ceil(r_lo - 1e-9)), 0)
        y2 = min(int(np.floor(r_hi + 1e-9)), self.ts.length - 1) + 1
        if x2 > x1 and y2 > y1:
            return x1, x2, y1, y2
        rc = self.rc(0.5 * (xs[0] + xs[1]), 0.5 * (ys[0] + ys[1]))
        return None if rc is None else (rc[1], rc[1] + 1, rc[0], rc[0] + 1)

    def outline_xy(self):
        """The outline of the pixels with data: the first and the last of
        each row - a frame's footprint, without its holes."""
        valid = self.base_keep
        rows = np.flatnonzero(valid.any(axis=1))
        if not rows.size:
            return [], []
        step = max(1, int(np.ceil(rows.size / MAX_VERTICES)))
        rows = np.unique(np.append(rows[::step], rows[-1]))
        sub = valid[rows]
        first = sub.argmax(axis=1).astype(float)
        last = (self.ts.width - 1 - sub[:, ::-1].argmax(axis=1)).astype(float)
        r = np.concatenate([[rows[0] - 0.5], rows, [rows[-1] + 0.5]])
        left = np.concatenate([[first[0]], first, [first[-1]]]) - 0.5
        right = np.concatenate([[last[0]], last, [last[-1]]]) + 0.5
        cols = np.concatenate([left, right[::-1], left[:1]])
        rws = np.concatenate([r, r[::-1], r[:1]])
        return self.x0 + (cols + 0.5) * self.dx, self.y0 + (rws + 0.5) * self.dy

    ##---------- dates
    def nearest_date(self, day):
        return int(np.argmin(np.abs(self.ords - day.toordinal())))

    def date_at(self, day):
        """The date nearest `day`, or None outside the file's time span (by more
        than half its usual interval)."""
        t = day.toordinal()
        if t < self.ords[0] - self.half_step or t > self.ords[-1] + self.half_step:
            return None
        return self.nearest_date(day)

    ##---------- the reference
    @property
    def keep(self):
        """The pixels used: with data at the reference date, and in the mask if on."""
        if self.v.use_mask and self.ts.keep is not None:
            return self.base_keep & self.ts.keep
        return self.base_keep

    def area_mean_2d(self, data):
        x1, x2, y1, y2 = self.ref
        sub = np.where(self.keep[y1:y2, x1:x2], data[y1:y2, x1:x2], np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return float(np.nanmean(sub))

    def has_data(self, box):
        x1, x2, y1, y2 = box
        return bool(self.keep[y1:y2, x1:x2].any())

    def update_ref_values(self):
        self.ref1 = self.ts.area_mean(self.ref, self.keep)

    def last_cum(self):
        if self._last_cum is None:
            self._last_cum = self.ts.slice(self.ts.n_date - 1) - self.ref_slice
        return self._last_cum

    ##---------- the time series at a pixel
    def compute(self, row, col):
        """The series of a pixel - or the date-by-date statistics of the circle
        around it - relative to this file's reference area and date."""
        ts, v = self.ts, self.v
        if v.radius:
            rows, cols = circle_pixels_rc(ts, row, col, v.radius)
        else:
            rows, cols = np.array([row]), np.array([col])
        good = self.keep[rows, cols]
        res = dict(frame=self, row=row, col=col, n_total=int(rows.size),
                   n_kept=int(good.sum()), stats=None, vel=None, gap=None, indices={})
        kr, kc = rows[good], cols[good]
        pick = (kr, kc) if kr.size else (rows, cols)
        for name, (data, unit) in ts.layers.items():
            if name != 'mask':
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    res['indices'][name] = (float(np.nanmedian(data[pick])), unit)
        if not kr.size:
            return res
        rel = relative_series(ts.pixels(kr, kc), self.ref1, self.ix_ref) * v.f
        res['stats'] = series_stats(rel)
        res['vel'] = float(np.nanmedian(self.vel[kr, kc] - self.area_mean_2d(self.vel)) * v.f)
        res['gap'] = ts.gap_fraction(kr, kc)
        return res


##############################################################################
class MultiTimeSeriesViewer:
    """See the module docstring. Built from the parsed `viewinsar tseries`
    options, with inps.files holding several files."""

    def __init__(self, inps):
        for key, value in vars(inps).items():
            setattr(self, key, value)
        self.frames, self.lookup, self.crs = [], None, None
        self.geocoded = self.projected = self.one_grid = False
        self.fig_map = self.ax_map = self.fig_ts = self.ax_ts = None
        self.cbar = self.slider = self.radio = self.mask_check = None
        self.point = None             # (x, y) of the selection, map coordinates
        self.result = None            # what the time-series window shows
        self.cum_mode, self.date_index, self.layer = False, None, None
        self.fit_visible = list(FIT_DEFAULT)
        self.fit_lines, self.marker_artists = [], []

    def vprint(self, *args):
        if self.print_msg:
            print(*args)

    ##---------- setup
    def open(self):
        if getattr(self, 'argv', None) is not None:
            print('viewinsar tseries ' + ' '.join(self.argv))
        if not self.disp_fig:
            plt.switch_backend('Agg')
        labels = self.labels or file_labels(self.files)
        for k, fname in enumerate(self.files):
            ts = TimeSeries(fname, verbose=self.print_msg,
                            scale=per_file(getattr(self, 'data_scale', None), k))
            frame = Frame(self, ts, k, labels[k], fname)
            self.frames.append(frame)
            self.vprint(f'{frame.tag}: {labels[k]}: {ts.kind} time series {fname}: '
                        f'{ts.n_date} dates {ts.dates[0]} - {ts.dates[-1]}, '
                        f'{ts.length} x {ts.width} pixels')
        self.check_grids()
        if self.mask_file:
            self.need_one_grid('--mask')
            for frame in self.frames:
                frame.ts.set_mask_file(self.mask_file)
        if not self.geocoded:
            self.init_lookup()
        self.f = TO_UNIT[self.unit]
        self.cmap = get_cmap(self.colormap)
        self.has_mask = any(fr.ts.keep is not None for fr in self.frames)
        self.use_mask = self.has_mask and not self.nomask

        self.set_ref_date()
        for fr in self.frames:
            fr.ref_slice = fr.ts.slice(fr.ix_ref)
            fr.base_keep = np.isfinite(fr.ref_slice)
            fr.vel = fr.ts.velocity()
            self.vprint(f'{fr.tag}: velocity: {fr.ts.vel_source}'
                        + (f'; mask: {fr.ts.mask_source}' if fr.ts.keep is not None else ''))
        if self.has_mask and self.nomask:
            self.vprint('mask: off (--nomask)')
        self.aspect = self.map_aspect()
        self.set_initial_ref()
        for fr in self.frames:
            fr.update_ref_values()
            if not fr.has_data(fr.ref):
                print(f'WARNING: {fr.tag}: the reference area {fr.ref[0]}:{fr.ref[1]}/'
                      f'{fr.ref[2]}:{fr.ref[3]} of {fr.label} has no valid pixel - '
                      'right-drag a new one')
        self.point = self.initial_point()
        if not self.fig_size:
            self.fig_size = [8.0, 4.8]
        return

    def check_grids(self):
        """All geocoded in one CRS, or all radar-coded on one grid."""
        from .raster import crs_name, crs_of, same_crs

        frames = self.frames
        geo = [fr.ts.geocoded for fr in frames]
        if all(geo):
            self.geocoded = True
            proj = {fr.ts.projected for fr in frames}
            crs = [crs_of(fr.ts.atr, fr.file) for fr in frames]
            if len(proj) > 1 or (None in crs and len(set(crs)) > 1):
                raise ValueError('several files must share one CRS: '
                                 + ', '.join(f'{fr.label} is {crs_name(c)}'
                                             for fr, c in zip(frames, crs)))
            for fr, code in zip(frames[1:], crs[1:]):
                if code is not None and not same_crs(crs[0], code):
                    raise ValueError(f'several files must share one CRS: {frames[0].label} '
                                     f'is {crs_name(crs[0])}, {fr.label} {crs_name(code)}')
            if crs[0] is None:
                print('NOTE: the files do not say their CRS: taken to be the same one')
            self.crs, self.projected = crs[0], proj.pop()
        elif not any(geo):
            shape = (frames[0].ts.length, frames[0].ts.width)
            for fr in frames[1:]:
                if (fr.ts.length, fr.ts.width) != shape:
                    raise ValueError('radar-coded files must share one grid: '
                                     f'{frames[0].label} is {shape[1]} x {shape[0]}, '
                                     f'{fr.label} {fr.ts.width} x {fr.ts.length} (width x '
                                     'length) - geocode them to see them together')
        else:
            raise ValueError('several files must all be geocoded (in one CRS), or all '
                             'radar-coded on one grid: '
                             + ', '.join(f'{fr.label} is '
                                         + ('geocoded' if fr.ts.geocoded else 'radar-coded')
                                         for fr in frames))
        self.one_grid = all(self.same_grid(frames[0], fr) for fr in frames[1:])

    @staticmethod
    def same_grid(a, b):
        """The same pixels: same size, and the same corner and steps to 1/1000 pixel."""
        if (a.ts.length, a.ts.width) != (b.ts.length, b.ts.width):
            return False
        tol_x, tol_y = 1e-3 * abs(a.dx), 1e-3 * abs(a.dy)
        return (abs(a.x0 - b.x0) <= tol_x and abs(a.dx - b.dx) <= tol_x
                and abs(a.y0 - b.y0) <= tol_y and abs(a.dy - b.dy) <= tol_y)

    def need_one_grid(self, option):
        if not self.one_grid:
            raise ValueError(f'{option} is in the pixels of one grid, and these files are '
                             'on different grids: give --ref-geo / --ref-lalo, --lalo / '
                             '--p-geo instead' if option != '--mask' else
                             '--mask: one mask cannot fit files on different grids - each '
                             'file brings its own')

    def init_lookup(self):
        """Radar-coded MintPy data: lat/lon through geometryRadar.h5, when found."""
        from .radar import RadarLookup, find_lookup

        try:
            path = find_lookup(self.files[0], self.frames[0].ts.atr, given=self.lookup_file,
                               verbose=self.print_msg)
            self.lookup = RadarLookup(path) if path else None
        except Exception as err:     # noqa: BLE001 - lat/lon is a nicety here
            print(f'WARNING: {err}')
            self.lookup = None

    def set_ref_date(self):
        """--ref-date, else the latest of the files' own (their first date, or
        MintPy's REF_DATE): each file takes its date nearest to it."""
        own = [fr.ts.ref_date or fr.ts.dates[0] for fr in self.frames]
        target = self.ref_date or max(own)
        self.ref_day = dt.datetime.strptime(target, '%Y%m%d')
        for fr in self.frames:
            fr.ix_ref = fr.nearest_date(self.ref_day)
            if fr.ts.dates[fr.ix_ref] != target:
                print(f'NOTE: {fr.tag}: {fr.label} has no {target}: its reference date is '
                      f'{fr.ts.dates[fr.ix_ref]}, the nearest')
        self.vprint(f'reference date: {target}')

    def map_aspect(self):
        """y / x of the map: equal distances on the ground."""
        if self.geocoded and not self.projected:
            _, _, bottom, top = self.map_extent()
            return 1.0 / max(np.cos(np.deg2rad(0.5 * (bottom + top))), 0.05)
        if self.geocoded:
            return 1.0
        ts = self.frames[0].ts
        try:
            x_km, y_km = ts.pixel_size_km(ts.length // 2)
            return y_km / x_km
        except ValueError:
            if self.radius:
                raise
            return 1.0

    def map_extent(self):
        ext = np.array([fr.extent for fr in self.frames])
        if not self.geocoded:
            return tuple(ext[0])
        return ext[:, 0].min(), ext[:, 1].max(), ext[:, 2].min(), ext[:, 3].max()

    ##---------- coordinates
    def to_map(self, lat, lon):
        """lat/lon - or northing/easting on a projected grid - to map x, y."""
        if self.geocoded:
            if self.projected and abs(lat) <= 90 and abs(lon) <= 360:
                if self.crs is None:
                    raise ValueError('the CRS of the files is unknown: give northing and '
                                     'easting instead of lat/lon')
                from pyproj import Transformer
                east, north = Transformer.from_crs(4326, self.crs,
                                                   always_xy=True).transform(lon, lat)
                return float(east), float(north)
            return float(lon), float(lat)
        if self.lookup is None:
            raise ValueError('radar-coded data without geometryRadar.h5 (--lookup): '
                             'give the pixel instead (--yx Y X, --point X/Y)')
        rows, cols = self.lookup.to_rowcol(np.array([lon]), np.array([lat]))
        if not np.isfinite(rows[0] + cols[0]):
            raise ValueError(f'{lat}, {lon} is outside the radar frame')
        return float(cols[0]), float(rows[0])

    def lalo_of(self, x, y):
        """(lat, lon) - northing, easting on a projected grid - of a map point."""
        if self.geocoded:
            return float(y), float(x)
        if self.lookup is not None:
            row, col = int(np.floor(y + 0.5)), int(np.floor(x + 0.5))
            lon, lat = self.lookup.to_lonlat(np.array([row]), np.array([col]))
            return float(lat[0]), float(lon[0])
        return None

    def covering(self, x, y):
        """The files whose grid holds map point (x, y)."""
        return [fr for fr in self.frames if fr.rc(x, y) is not None]

    def set_initial_ref(self):
        first = self.frames[0]
        if self.ref_area or self.ref_yx:
            self.need_one_grid('--ref' if self.ref_area else '--ref-yx')
            if self.ref_area:
                box = parse_range(self.ref_area, first.ts.width, first.ts.length)
            else:
                row, col = self.ref_yx
                if not (0 <= row < first.ts.length and 0 <= col < first.ts.width):
                    raise ValueError(f'--ref-yx {row} {col} falls outside the {first.ts.width}'
                                     f' x {first.ts.length} (width x length) grid')
                box = (col, col + 1, row, row + 1)
            for fr in self.frames:
                fr.ref = box
        elif self.ref_geo or self.ref_lalo:
            if self.ref_geo:
                lon1, lon2, lat1, lat2 = parse_geo_range(self.ref_geo)
                (xa, ya), (xb, yb) = self.to_map(lat1, lon1), self.to_map(lat2, lon2)
                what = f'--ref-geo {self.ref_geo}'
            else:
                xa, ya = xb, yb = self.to_map(*self.ref_lalo)
                what = f'--ref-lalo {self.ref_lalo[0]} {self.ref_lalo[1]}'
            if not self.move_ref((xa, xb), (ya, yb)):
                raise ValueError(f'{what}: no file has a valid pixel there')
        elif not self.one_grid:
            print('NOTE: each file is relative to its own reference area: right-drag one '
                  'where they overlap (or give --ref-geo) to put them on the same ground')

    def move_ref(self, xs, ys):
        """A map rectangle as the new reference area of every file with data
        in it; the others keep theirs. Returns the files moved."""
        moved = []
        for fr in self.frames:
            box = fr.box_of(xs, ys)
            if box is not None and fr.has_data(box):
                fr.ref = box
                moved.append(fr)
        if moved:
            for fr in self.frames:
                x1, x2, y1, y2 = fr.ref
                if fr in moved:
                    self.vprint(f'{fr.tag}: reference area {x1}:{x2}/{y1}:{y2}')
                else:            # its series are not on the same ground: always say so
                    print(f'NOTE: {fr.tag}: {fr.label} keeps its reference area {x1}:{x2}/'
                          f'{y1}:{y2} - it has no valid pixel in the new one')
        return moved

    def initial_point(self):
        first = self.frames[0]
        if self.yx or self.point_xy:
            self.need_one_grid('--yx' if self.yx else '--point')
            row, col = self.yx if self.yx else (self.point_xy[1], self.point_xy[0])
            if not (0 <= row < first.ts.length and 0 <= col < first.ts.width):
                raise ValueError(f'the point {col}/{row} (x/y) falls outside the '
                                 f'{first.ts.width} x {first.ts.length} (width x length) grid')
            return first.xy(row, col)
        if self.lalo or self.point_geo:
            lat, lon = self.lalo if self.lalo else (self.point_geo[1], self.point_geo[0])
            x, y = self.to_map(lat, lon)
            if not self.covering(x, y):
                raise ValueError(f'the point {lat}, {lon} is outside every file')
            return x, y
        x1, x2, y1, y2 = first.ref                     # LiCSBAS: the reference area's centre
        return first.xy((y1 + y2) // 2, (x1 + x2) // 2)

    ##---------- main
    def plot(self):
        self.make_map_figure()
        self.make_ts_figure()
        self.select(*self.point)
        if self.save_fig:
            self.save_outputs()
        if self.disp_fig:
            print(self.help_msg())
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                plt.show()
        return

    def help_msg(self):
        msg = '\n' + '-' * 74
        msg += '\n  left double click : the time series of every file there'
        if self.radius:
            msg += (f'\n                      - the median of a circle of {self.radius:g} km, '
                    'with its 25-75 % range')
        msg += '\n  right drag        : a new reference area, for every file with data in it'
        msg += '\n  slider            : the cumulative displacement at that date'
        msg += f'\n  key "e"           : export the time series (.csv) and the figures ({self.fig_ext}),'
        msg += f'\n                      into {os.path.abspath(self.outdir or ".")}'
        msg += '\n' + '-' * 74
        return msg

    def make_map_figure(self):
        # what the map can show: the velocity, then the layers of any file
        self.layers = {'vel': 'vel'}
        for fr in self.frames:
            for name, (_data, unit) in fr.ts.layers.items():
                self.layers.setdefault(name, 'mask' if name == 'mask' else unit)

        left, right, bottom, top = self.map_extent()
        ratio = abs(top - bottom) * self.aspect / max(abs(right - left), 1e-12)
        size = map_window_size(ratio, len(self.layers), self.ref_note().count('\n') + 1,
                               self.fig_size_map)
        fig = self.fig_map = _new_figure('Velocity / Cumulative Displacement', size)
        place = place_map_window(fig, len(self.layers), right=0.85)   # 0.85: colour bar
        ax = self.ax_map = fig.add_axes(place['ax'])
        self.ref_text = fig.text(0.06 / size[0], 1 - 0.06 / size[1], '', fontsize=8, va='top')
        for fr in self.frames:
            fr.image = ax.imshow(np.full((fr.ts.length, fr.ts.width), np.nan), cmap=self.cmap,
                                 extent=fr.extent, origin='upper', interpolation='nearest',
                                 zorder=1)
            fr.outline, = ax.plot(*fr.outline_xy(), '-', color=fr.color, lw=1.2, zorder=4,
                                  path_effects=HALO, label=f'{fr.tag}: {fr.label}')
            fr.ref_line, = ax.plot([], [], '--', color=fr.color, lw=1.4, zorder=6,
                                   path_effects=HALO)
        # a margin, so that the outlines on the edges show; radar: (length - 0.5,
        # -0.5), rows top-down
        pad_x, pad_y = 0.015 * (right - left), 0.015 * (top - bottom)
        ax.set_xlim(left - pad_x, right + pad_x)
        ax.set_ylim(bottom - pad_y, top + pad_y)
        if not self.geocoded and not self.noflip:
            orbit = str(self.frames[0].ts.atr.get('ORBIT_DIRECTION', '')).lower()
            if orbit.startswith('a'):
                ax.invert_yaxis()             # as MintPy: north up
            elif orbit.startswith('d'):
                ax.invert_xaxis()             # as MintPy: east right
        ax.set_aspect(self.aspect)
        # the colour bar as tall as the map, whatever its shape
        cax = make_axes_locatable(ax).append_axes('right', size='3%', pad=0.1)
        self.cbar = fig.colorbar(self.frames[0].image, cax=cax)
        ax.legend(handles=[fr.outline for fr in self.frames], fontsize=7, loc='best',
                  framealpha=0.7)

        rad_ax = fig.add_axes(place['radio'])
        self.radio = RadioButtons(rad_ax, tuple(self.layers))
        for label in self.radio.labels:
            label.set_fontsize(8)
        self.radio.on_clicked(self.show_layer)
        if self.has_mask:
            box = fig.add_axes(place['mask'])
            self.mask_check = CheckButtons(box, ['mask'], [self.use_mask])
            self.mask_check.on_clicked(self.toggle_mask)

        # the date slider: every date of every file
        self.days = sorted({d for fr in self.frames for d in fr.ts.datetimes})
        self.date_nums = mdates.date2num(self.days)
        ref_num = mdates.date2num(self.ref_day)
        sax = fig.add_axes(place['slider'], yticks=[])
        self.slider = Slider(sax, 'yr', self.date_nums[0] - 3, self.date_nums[-1] + 3,
                             valinit=ref_num, valfmt='')
        sax.bar(self.date_nums, np.ones(len(self.days)), facecolor='black', width=4)
        sax.bar(ref_num, 1, facecolor='red', width=8)
        sax.set_axis_on()
        for side in ('left', 'right', 'top'):
            sax.spines[side].set_visible(False)
        sax.set_yticks([])
        locator = mdates.AutoDateLocator()
        sax.xaxis.set_major_locator(locator)
        sax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        sax.tick_params(labelsize=8)
        self.slider.on_changed(self.show_date)

        self.selector = RectangleSelector(ax, self.on_ref_select, useblit=True, button=[3],
                                          spancoords='pixels', interactive=False)
        fig.canvas.mpl_connect('button_press_event', self.on_click)
        fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.draw_refs()
        self.show_layer('vel')
        return

    def make_ts_figure(self):
        fig = self.fig_ts = _new_figure('Time-series', self.fig_size)
        self.ax_ts = fig.add_axes([0.12, 0.14, 0.7, 0.8])
        self.side_text = fig.text(0.83, 0.95, '', fontsize=8, va='top')
        box = fig.add_axes([0.83, 0.10, 0.16, 0.25])
        self.fit_check = CheckButtons(box, list(MODELS), list(self.fit_visible))
        self.fit_check.on_clicked(self.toggle_fit)
        fig.canvas.mpl_connect('key_press_event', self.on_key)
        return

    ##---------- the map
    def _show(self, arrays, cmap, limits):
        for fr, data in zip(self.frames, arrays):
            fr.image.set_data(fr.display(data))
            fr.image.set_cmap(cmap)
            fr.image.set_clim(*limits)

    def show_layer(self, name):
        """The velocity (relative to each file's reference area) or another layer."""
        kind = self.layers[name]
        arrays = []
        for fr in self.frames:
            keep = fr.keep
            if kind == 'vel':
                data = np.where(keep, (fr.vel - fr.area_mean_2d(fr.vel)) * self.f, np.nan)
            elif name in fr.ts.layers:
                data = fr.ts.layers[name][0]
                data = data if kind == 'mask' else np.where(keep, data, np.nan)
            else:
                data = np.full((fr.ts.length, fr.ts.width), np.nan)
            arrays.append(data)
        values = _finite(arrays)
        cmap, unit = self.cmap, kind
        if kind == 'vel':
            limits = self.vlim or auto_limits(values, self.auto_crange)
            unit = f'{self.unit}/yr'
        elif kind == 'mask':
            limits, cmap, unit = (0, 1), 'viridis', ''
        else:
            limits = auto_limits(values, self.auto_crange)
            cmap = 'viridis_r'
            if name in ('coh_avg', 'n_unw', 'maxTlen', 'temporalCoherence'):
                cmap = 'viridis'
            elif name == 'mli':
                cmap = 'gray'
            elif name == 'hgt':
                cmap, limits = 'terrain', (-limits[1] / 3, limits[1])
        self._show(arrays, cmap, limits)
        self.cbar.set_label(unit)
        self.ax_map.set_title(name)
        self.cum_mode, self.layer = False, name
        self.fig_map.canvas.draw_idle()
        return

    def cum_limits(self):
        """LiCSBAS: the colour range of the cumulative displacement at the last date."""
        values = []
        for fr in self.frames:
            last = np.where(fr.keep, fr.last_cum(), np.nan)
            values.append((last - fr.area_mean_2d(last)) * self.f)
        return auto_limits(_finite(values), self.auto_crange)

    def show_date(self, value):
        """The slider: each file's cumulative displacement at its date nearest
        the slider's - none before its first date or after its last."""
        idx = int(np.argmin(np.abs(self.date_nums - value)))
        day = self.days[idx]
        arrays = []
        for fr in self.frames:
            j = fr.date_at(day)
            if j is None:
                arrays.append(np.full((fr.ts.length, fr.ts.width), np.nan))
                continue
            disp = np.where(fr.keep, fr.ts.slice(j) - fr.ref_slice, np.nan)
            arrays.append((disp - fr.area_mean_2d(disp)) * self.f)
        self._show(arrays, self.cmap, self.dlim or self.cum_limits())
        self.cbar.set_label(self.unit)
        self.ax_map.set_title(f'{day:%Y/%m/%d} (Ref: {self.ref_day:%Y/%m/%d})')
        self.cum_mode, self.date_index = True, idx
        self.fig_map.canvas.draw_idle()
        return

    def redraw_map(self):
        if self.cum_mode:
            self.show_date(self.date_nums[self.date_index])
        else:
            self.show_layer(self.layer)

    def ref_note(self):
        """The text at the top left of the map window."""
        note = f'{HELP_TEXT}\n\nRef area:\n dashed, in the\n colour of its file'
        if self.radius:
            note += f'\n\nCircle:\n r = {self.radius:g} km'
        return note

    def draw_refs(self):
        for fr in self.frames:
            fr.ref_line.set_data(*fr.box_xy(fr.ref))
        self.ref_text.set_text(self.ref_note())

    def draw_marker(self):
        """The point - its pixel, or the circle around it - on the map."""
        for art in self.marker_artists:
            art.remove()
        self.marker_artists = []
        x, y = self.point
        ax = self.ax_map
        here = self.covering(x, y)
        if here:
            fr = here[0]
            row, col = fr.rc(x, y)
            if self.radius:
                x_km, y_km = fr.ts.pixel_size_km(row)
                ring = patches.Ellipse((x, y), 2 * self.radius / x_km * abs(fr.dx),
                                       2 * self.radius / y_km * abs(fr.dy), fill=False,
                                       ec='k', lw=1.6, zorder=7,
                                       path_effects=[patheffects.withStroke(linewidth=3.4,
                                                                            foreground='w')])
                self.marker_artists.append(ax.add_patch(ring))
            else:
                self.marker_artists += ax.plot(*fr.box_xy((col, col + 1, row, row + 1)), 'k',
                                               linewidth=3, zorder=7)
        self.marker_artists += ax.plot([x], [y], 'Pk', zorder=8)
        self.fig_map.canvas.draw_idle()

    ##---------- events
    def _toolbar_busy(self):
        toolbar = getattr(self.fig_map.canvas, 'toolbar', None)
        return bool(toolbar is not None and getattr(toolbar, 'mode', ''))

    def on_click(self, event):
        """Left double click on the map: the time series of every file there."""
        if (event.inaxes is not self.ax_map or event.button != 1 or not event.dblclick
                or event.xdata is None or self._toolbar_busy()):
            return
        x, y = float(event.xdata), float(event.ydata)
        if self.covering(x, y):
            self.select(x, y)

    def on_ref_select(self, eclick, erelease):
        """Right drag: a new reference area, for every file with data in it."""
        corners = (eclick.xdata, erelease.xdata, eclick.ydata, erelease.ydata)
        if None in corners:
            return
        if not self.move_ref(corners[:2], corners[2:]):
            print('Selected ref area has all nan!! Reselect different ref area.')
            return
        self.refresh()

    def toggle_mask(self, _label=None):
        self.use_mask = not self.use_mask
        self.refresh()

    def toggle_fit(self, label):
        model = MODELS.index(label)
        self.fit_visible[model] = not self.fit_visible[model]
        for lines in self.fit_lines:
            lines[model].set_visible(self.fit_visible[model])
        self.fig_ts.canvas.draw_idle()

    def on_key(self, event):
        if event.key == 'e':
            self.save_outputs()

    def refresh(self):
        """After a new reference area or mask: the map and the time series again."""
        for fr in self.frames:
            fr.update_ref_values()
        self.draw_refs()
        self.redraw_map()
        if self.point is not None:
            self.select(*self.point)

    ##---------- the time series
    def select(self, x, y):
        self.point = (float(x), float(y))
        self.draw_marker()
        self.result = self.compute(*self.point)
        self.draw_ts()

    def compute(self, x, y):
        """The series of every file that covers the map point (x, y)."""
        res = dict(x=x, y=y, lalo=self.lalo_of(x, y), frames=[])
        for fr in self.frames:
            rc = fr.rc(x, y)
            if rc is not None:
                res['frames'].append(fr.compute(*rc))
        return res

    def where(self, res):
        if self.geocoded:
            return (f'@(N {res["y"]:.0f}, E {res["x"]:.0f})' if self.projected
                    else f'@({res["y"]:.4f}, {res["x"]:.4f})')
        here = res['frames'][0] if res['frames'] else None
        if here is None:
            return f'@({res["x"]:.0f}, {res["y"]:.0f})'
        return f'@({here["col"]}, {here["row"]})'

    def info_text(self, res):
        """The text beside the time series, as LiCSBAS - a block per file."""
        text = ''
        if res['lalo'] is not None:
            lat, lon = res['lalo']
            text = (f'N:{lat:.1f}\nE:{lon:.1f}\n\n' if self.projected
                    else f'Lat:{lat:.5f}\nLon:{lon:.5f}\n\n')
        dates = sorted({fr.ts.dates[fr.ix_ref] for fr in self.frames})
        if len(dates) == 1:
            text += f'Ref date:\n {dates[0]}\n'
        if self.radius:
            text += f'Circle: r = {self.radius:g} km\n'
        for r in res['frames']:
            fr = r['frame']
            x1, x2, y1, y2 = fr.ref
            text += f'\n{fr.tag}: {_short(fr.label)}\n'
            if self.radius:
                text += f' {r["n_total"]} px, {r["n_kept"]} used\n'
            elif fr.ts.keep is not None:
                text += f' mask: {int(bool(fr.ts.keep[r["row"], r["col"]]))}\n'
            text += f' Ref X {x1}:{x2} Y {y1}:{y2}\n'
            if len(dates) > 1:
                text += f' Ref date {fr.ts.dates[fr.ix_ref]}\n'
        here = [r['frame'] for r in res['frames']]
        others = [fr.tag for fr in self.frames if fr not in here]
        if others:
            text += f'\nnot here: {", ".join(others)}\n'
        return text

    def draw_ts(self):
        res, ax = self.result, self.ax_ts
        ax.cla()
        ax.grid(zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlabel('Time')
        ax.set_ylabel(f'Displacement ({self.unit})')
        self.side_text.set_text(self.info_text(res))
        self.fit_lines = []
        where = self.where(res)
        shown = [r for r in res['frames'] if r['stats'] is not None]
        if not shown:
            ax.set_title(f'NaN {where}' + (f' - no valid pixel within {self.radius:g} km'
                                           if self.radius else ''), fontsize=10)
            self.fig_ts.canvas.draw_idle()
            return

        locator = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        for k, r in enumerate(shown):
            fr, stats = r['frame'], r['stats']
            ts, color = fr.ts, fr.color
            dates = np.array(ts.datetimes, dtype=object)
            n_fit = max(int((ts.datetimes[-1] - ts.datetimes[0]).days // 10), 1)
            fit_dates = [ts.datetimes[0] + dt.timedelta(days=10 * i) for i in range(n_fit + 1)]
            fit_years = np.arange(n_fit + 1) * 10 / 365.25
            lines = []
            for model in range(len(MODELS)):
                curve = model_curve(ts.years, stats['median'], fit_years, model)
                lines += ax.plot(fit_dates, curve, '-', color=color, alpha=0.6, zorder=3,
                                 visible=self.fit_visible[model])
            self.fit_lines.append(lines)
            label = f'{fr.tag}: {fr.label}'
            if self.radius:
                med = stats['median']
                yerr = np.vstack([med - stats['q25'], stats['q75'] - med])
                ax.errorbar(dates, med, yerr=yerr, fmt='o', ms=4, color=color, ecolor=color,
                            elinewidth=0.8, capsize=2, alpha=0.6, zorder=5 + k, label=label)
            else:
                ax.scatter(dates, stats['median'], color=color, alpha=0.6, zorder=5 + k,
                           label=label)
            if r['gap'] is not None:
                gaps = r['gap'] >= (GAP_SHARE if self.radius else 1.0)
                if gaps.any():
                    mid = dates[:-1][gaps] + (dates[1:][gaps] - dates[:-1][gaps]) / 2
                    ax.vlines(list(mid), 0, 1, transform=ax.get_xaxis_transform(), zorder=1,
                              alpha=0.6, colors=color, linestyles='--',
                              label=f'Gap of ifg network ({fr.tag})')
        unit = f'{self.unit}/yr'
        title = ', '.join(f'vel({r["frame"].tag}) = {r["vel"]:.1f}' for r in shown)
        title += f' {unit} {where}'
        if self.radius:
            counts = ', '.join(f'{r["n_kept"]} ({r["frame"].tag})' for r in shown)
            title += f'\nmedian of {counts} pixels within {self.radius:g} km'
        ax.set_title(title, fontsize=10)
        if self.ylim:
            ax.set_ylim(self.ylim)
        ax.legend(fontsize=8, title='median, 25-75 %' if self.radius else None,
                  title_fontsize=8)
        self.fig_ts.canvas.draw_idle()
        return

    ##---------- output
    @property
    def outfile_base(self):
        x, y = self.point
        if self.geocoded:
            name = (f'tseries_N{y:.0f}_E{x:.0f}' if self.projected
                    else f'tseries_lat{y:.4f}_lon{x:.4f}')
        else:
            row, col = self.frames[0].rc(x, y) or (int(round(y)), int(round(x)))
            name = f'tseries_Y{row}X{col}'
        name += f'_r{self.radius:g}km' if self.radius else ''
        return os.path.join(self.outdir, name)

    def save_outputs(self):
        """The time series as .csv, how they were made (_info.txt), and both figures."""
        res = self.result
        shown = [r for r in res['frames'] if r['stats'] is not None] if res else []
        if not shown:
            print('nothing to save - double-click a point with data on the map.')
            return []
        base = self.outfile_base
        os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
        unit, written = self.unit, []

        # 1. the time series: every date of the files shown, a column set per file
        days = sorted({d for r in shown for d in r['frame'].ts.datetimes})
        header = ['date', 'year']
        for r in shown:
            tag = f'_{r["frame"].tag}'
            header += ([f'median_{unit}{tag}', f'q25_{unit}{tag}', f'q75_{unit}{tag}',
                        f'n_px{tag}'] if self.radius else [f'disp_{unit}{tag}'])
        where = [{d: i for i, d in enumerate(r['frame'].ts.datetimes)} for r in shown]
        out = f'{base}.csv'
        with open(out, 'w') as fid:
            fid.write(','.join(header) + '\n')
            for day in days:
                cells = [f'{day:%Y%m%d}', f'{_decimal_year(day):.4f}']
                for r, index in zip(shown, where):
                    i, stats = index.get(day), r['stats']
                    if self.radius:
                        cells += (['nan', 'nan', 'nan', '0'] if i is None else
                                  [f'{stats["median"][i]:.4f}', f'{stats["q25"][i]:.4f}',
                                   f'{stats["q75"][i]:.4f}', f'{int(stats["count"][i])}'])
                    else:
                        cells += ['nan' if i is None else f'{stats["median"][i]:.4f}']
                fid.write(','.join(cells) + '\n')
        written.append(out)

        # 2. how they were made
        out = f'{base}_info.txt'
        with open(out, 'w') as fid:
            fid.write(f'files          : {len(self.frames)}, on one map; the series of '
                      f'{len(shown)} of them here\n')
            if res['lalo'] is not None:
                tag = 'northing/easting' if self.projected else 'lat/lon'
                fid.write(f'{tag:<15}: {res["lalo"][0]:.6f}, {res["lalo"][1]:.6f}\n')
            if self.radius:
                fid.write(f'circle         : r = {self.radius:g} km\n')
                fid.write('statistics     : per date, the median of the pixels and their 25th '
                          'and 75th percentiles (the error bars)\n')
            fid.write(f'reference date : {self.ref_day:%Y%m%d} - each file at its own '
                      'nearest date (every pixel 0 there)\n')
            fid.write(f'mask           : {"on" if self.use_mask else "off"}\n')
            fid.write(f'unit           : {unit}\n')
            fid.write('columns        : _N of the .csv = file N below; "nan" where a file '
                      'has no such date\n')
            for r in res['frames']:
                fr = r['frame']
                x1, x2, y1, y2 = fr.ref
                xs, ys = fr.box_xy(fr.ref)
                fid.write(f'\n[{fr.tag}] {fr.label}\n')
                fid.write(f'  file           : {os.path.abspath(fr.file)} ({fr.ts.kind})\n')
                fid.write(f'  pixel          : x/y {r["col"]}/{r["row"]} (start from 0)\n')
                if self.radius:
                    fid.write(f'  circle         : {r["n_total"]} pixels, {r["n_kept"]} used '
                              '(with data, in the mask)\n')
                fid.write(f'  reference area : {x1}:{x2}/{y1}:{y2} (x/y, start from 0)')
                if self.geocoded:
                    xy = 'E/N' if self.projected else 'lon/lat'
                    fid.write(f' = {xy} {min(xs):.6f}/{max(xs):.6f}/{min(ys):.6f}/'
                              f'{max(ys):.6f}')
                fid.write(', its mean subtracted at every date\n')
                fid.write(f'  reference date : {fr.ts.dates[fr.ix_ref]}\n')
                if fr.ts.keep is not None:
                    fid.write(f'  mask from      : {fr.ts.mask_source}\n')
                if r['vel'] is not None:
                    fid.write(f'  velocity       : {r["vel"]:.3f} {unit}/yr'
                              + (' (median of the circle)' if self.radius else '') + '\n')
                else:
                    fid.write('  velocity       : - (no valid pixel here)\n')
                fid.write(f'  velocity from  : {fr.ts.vel_source}\n')
                for name, (value, unit_i) in r['indices'].items():
                    fid.write(f'  {name:<15}: {_fmt(value, name.startswith("n_"))} '
                              f'{unit_i}'.rstrip() + '\n')
            for fr in self.frames:
                if fr not in [r['frame'] for r in res['frames']]:
                    fid.write(f'\n[{fr.tag}] {fr.label}: not at this point\n'
                              f'  file           : {os.path.abspath(fr.file)}\n')
        written.append(out)

        # 3. the figures
        out = self.outfile[0] if self.outfile else f'{base}{self.fig_ext}'
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        self.fig_ts.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
        written.append(out)
        out = f'{base}_map{self.fig_ext}'
        self.fig_map.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
        written.append(out)
        for path in written:
            print(f'saved: {path}')
        return written

    def close(self):
        for fig in (self.fig_map, self.fig_ts):
            if fig is not None:
                plt.close(fig)
        for fr in self.frames:
            fr.ts.close()
        return

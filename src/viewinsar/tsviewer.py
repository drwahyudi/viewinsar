"""`viewinsar tseries`: an interactive time-series viewer after LiCSBAS_plot_ts.py
(and MintPy's tsview.py), for LiCSBAS cum.h5 / cum_filt.h5 and MintPy time
series.

Two windows, as in LiCSBAS. The map: the velocity (or a noise index, or the
cumulative displacement at the date of the slider), relative to the
reference area. The time series of the selected point, with LiCSBAS's four
fit models. And circles: with -r KM a double click takes every pixel within
that radius, and each date shows their median, with the 25th - 75th
percentiles as an error bar.

    left double click   the time series there (of the circle, with -r)
    right drag          a new reference area
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
from .viewer import _new_figure

__all__ = ['TimeSeriesViewer', 'get_cmap']

COLORS = ('b', 'r')                          # file 1, file 2 - as LiCSBAS
FIT_DEFAULT = (True, True, False, False)     # Linear and Annual+L shown
GAP_SHARE = 0.5                              # a circle's gap: half its pixels or more
HELP_TEXT = 'Left-doubleclick:\n Plot time series\nRight-drag:\n Change ref area'


def get_cmap(name):
    """A colormap by name: MintPy's (roma, vik, ... as in LiCSBAS), else matplotlib's."""
    try:
        from mintpy.objects.colors import ColormapExt
        return ColormapExt(name).colormap
    except Exception:            # noqa: BLE001 - MintPy missing, or a name it lacks
        pass
    try:
        return plt.get_cmap(name)
    except ValueError:
        print(f'WARNING: unknown colormap {name!r} --> RdBu_r')
        return plt.get_cmap('RdBu_r')


def _fmt(value, integer=False):
    if value is None or not np.isfinite(value):
        return 'nan'
    if integer and float(value).is_integer():
        return f'{int(value)}'
    return f'{value:.1f}' if integer else f'{value:.2f}'


class TimeSeriesViewer:
    """See the module docstring. Built from the parsed `viewinsar tseries` options."""

    def __init__(self, inps):
        for key, value in vars(inps).items():
            setattr(self, key, value)
        self.ts = self.ts2 = self.lookup = None
        self.fig_map = self.ax_map = self.fig_ts = self.ax_ts = None
        self.image = self.cbar = self.slider = self.radio = self.mask_check = None
        self.point = None             # (row, col) of the selection
        self.result = None            # what the time-series window shows
        self.cum_mode, self.date_index, self.layer = False, None, None
        self.fit_visible = list(FIT_DEFAULT)
        self.fit_lines, self.marker_artists = [], []
        self._last_cum = None

    def vprint(self, *args):
        if self.print_msg:
            print(*args)

    ##---------- setup
    def open(self):
        if getattr(self, 'argv', None) is not None:
            print('viewinsar tseries ' + ' '.join(self.argv))
        if not self.disp_fig:
            plt.switch_backend('Agg')
        ts = self.ts = TimeSeries(self.file, mask_file=self.mask_file, vel_file=self.vel_file,
                                  results_dir=self.results_dir, verbose=self.print_msg)
        self.vprint(f'{ts.kind} time series {self.file}: {ts.n_date} dates '
                    f'{ts.dates[0]} - {ts.dates[-1]}, {ts.length} x {ts.width} pixels')
        if self.file2:
            self.open_second()

        # the reference date: -m, else the file's (MintPy REF_DATE), else the first
        self.ix_ref = ts.dates.index(ts.ref_date) if ts.ref_date else 0
        if self.ref_date:
            if self.ref_date in ts.dates:
                self.ix_ref = ts.dates.index(self.ref_date)
            else:
                print(f'No {self.ref_date} in the dates: the reference date is '
                      f'{ts.dates[self.ix_ref]}')
        self.vprint(f'reference date: {ts.dates[self.ix_ref]}')

        if not ts.geocoded:
            self.init_lookup()
        self.f = TO_UNIT[self.unit]
        self.cmap = get_cmap(self.colormap)
        self.ref_slice = ts.slice(self.ix_ref)                  # mm
        self.base_keep = np.isfinite(self.ref_slice)
        self.has_mask = ts.keep is not None
        self.use_mask = self.has_mask and not self.nomask
        if self.has_mask:
            self.vprint(f'mask: {ts.mask_source}' + (' (off: --nomask)' if self.nomask else ''))

        self.vel = ts.velocity()
        self.vel2 = self.ts2.velocity() if self.ts2 is not None else None
        self.vprint(f'velocity: {ts.vel_source}')
        try:
            x_km, y_km = ts.pixel_size_km(ts.length // 2)
            self.aspect = y_km / x_km
        except ValueError:
            if self.radius:
                raise
            self.aspect = 1.0

        self.ref = self.initial_ref()
        if not self.keep[self.ref[2]:self.ref[3], self.ref[0]:self.ref[1]].any():
            print('WARNING: the reference area {}:{}/{}:{} has no valid pixel - '
                  'right-drag a new one'.format(*self.ref))
        self.update_ref_values()
        self.point = self.initial_point()
        if not self.fig_size:
            self.fig_size = [8.0, 4.8]
        return

    def open_second(self):
        """--i2: a second time series on the same grid (LiCSBAS: cum.h5 beside
        cum_filt.h5), drawn in red, matched date by date."""
        ts2 = TimeSeries(self.file2, results_dir=self.results_dir, verbose=False)
        if (ts2.length, ts2.width) != (self.ts.length, self.ts.width):
            print(f'WARNING: --i2 {self.file2} is {ts2.length} x {ts2.width}, '
                  f'{self.file} {self.ts.length} x {self.ts.width} --> not used')
            ts2.close()
            return
        ts2.label = '2: ' + ts2.label.split(': ', 1)[-1]
        self.ts2 = ts2
        self.index2 = np.array([ts2.dates.index(d) if d in ts2.dates else -1
                                for d in self.ts.dates])
        if (self.index2 < 0).any():
            print(f'NOTE: {int((self.index2 < 0).sum())} date(s) of {self.file} are not in '
                  f'{self.file2}')
        self.vprint(f'second time series: {self.file2} ({ts2.label})')

    def init_lookup(self):
        """Radar-coded MintPy data: lat/lon through geometryRadar.h5, when found."""
        from .radar import RadarLookup, find_lookup

        try:
            path = find_lookup(self.file, self.ts.atr, given=self.lookup_file,
                               verbose=self.print_msg)
            self.lookup = RadarLookup(path) if path else None
        except Exception as err:     # noqa: BLE001 - lat/lon is a nicety here
            print(f'WARNING: {err}')
            self.lookup = None

    @property
    def keep(self):
        """The pixels used: with data at the reference date, and in the mask if on."""
        if self.use_mask:
            return self.base_keep & self.ts.keep
        return self.base_keep

    ##---------- coordinates
    def yx2lalo(self, row, col):
        if self.ts.geocoded:
            return self.ts.yx2lalo(row, col)
        if self.lookup is not None:
            lon, lat = self.lookup.to_lonlat(np.array([row]), np.array([col]))
            return float(lat[0]), float(lon[0])
        return None

    def lalo2yx(self, lat, lon):
        ts = self.ts
        if ts.geocoded:
            if ts.projected and abs(lat) <= 90 and abs(lon) <= 360:
                lat, lon = self._lonlat_to_grid(lat, lon)
            row, col = ts.lalo2yx(lat, lon)
        elif self.lookup is not None:
            rows, cols = self.lookup.to_rowcol(np.array([lon]), np.array([lat]))
            if not np.isfinite(rows[0] + cols[0]):
                raise ValueError(f'{lat}, {lon} is outside the radar frame')
            row, col = int(round(rows[0])), int(round(cols[0]))
        else:
            raise ValueError('radar-coded data without geometryRadar.h5 (--lookup): '
                             'give the pixel instead (--yx Y X, -p X/Y)')
        return self._inside(row, col, f'{lat}, {lon}')

    def _lonlat_to_grid(self, lat, lon):
        """lat/lon -> northing/easting of a projected grid (as MintPy's lalo2yx does)."""
        from pyproj import Transformer

        from .raster import crs_of

        crs = crs_of(self.ts.atr, self.file)
        east, north = Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat)
        return north, east

    def _inside(self, row, col, what):
        if not (0 <= row < self.ts.length and 0 <= col < self.ts.width):
            raise ValueError(f'{what} falls outside the {self.ts.width} x {self.ts.length} '
                             '(width x length) grid')
        return int(row), int(col)

    def initial_ref(self):
        ts = self.ts
        if self.ref_area:
            return parse_range(self.ref_area, ts.width, ts.length)
        if self.ref_geo:
            lon1, lon2, lat1, lat2 = parse_geo_range(self.ref_geo)
            (r1, c1), (r2, c2) = self.lalo2yx(lat1, lon1), self.lalo2yx(lat2, lon2)
            return min(c1, c2), max(c1, c2) + 1, min(r1, r2), max(r1, r2) + 1
        if self.ref_yx or self.ref_lalo:
            row, col = (self._inside(*self.ref_yx, 'the reference pixel') if self.ref_yx
                        else self.lalo2yx(*self.ref_lalo))
            return col, col + 1, row, row + 1
        return ts.ref_area

    def initial_point(self):
        if self.yx:
            return self._inside(*self.yx, f'--yx {self.yx[0]} {self.yx[1]}')
        if self.lalo:
            return self.lalo2yx(*self.lalo)
        if self.point_xy:
            return self._inside(self.point_xy[1], self.point_xy[0],
                                f'-p {self.point_xy[0]}/{self.point_xy[1]}')
        if self.point_geo:
            return self.lalo2yx(self.point_geo[1], self.point_geo[0])
        x1, x2, y1, y2 = self.ref                      # LiCSBAS: the reference area's centre
        return (y1 + y2) // 2, (x1 + x2) // 2

    ##---------- the reference
    def area_mean_2d(self, data):
        x1, x2, y1, y2 = self.ref
        sub = np.where(self.keep[y1:y2, x1:x2], data[y1:y2, x1:x2], np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return float(np.nanmean(sub))

    def update_ref_values(self):
        """The reference area's mean at every date (both files)."""
        self.ref1 = self.ts.area_mean(self.ref, self.keep)
        self.ref2 = None
        if self.ts2 is not None:
            self.ref2 = self.align2(self.ts2.area_mean(self.ref, self.keep))

    def align2(self, values):
        """The second file's dates -> the first file's (NaN where it has none)."""
        values = np.asarray(values, dtype=np.float64)
        out = np.full((len(self.index2),) + values.shape[1:], np.nan)
        good = self.index2 >= 0
        out[good] = values[self.index2[good]]
        return out

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
        msg += '\n  left double click : the time series there'
        if self.radius:
            msg += f' - the median of a circle of {self.radius:g} km, with its 25-75 % range'
        msg += '\n  right drag        : a new reference area'
        msg += '\n  slider            : the cumulative displacement at that date'
        msg += f'\n  key "e"           : export the time series (.csv) and the figures ({self.fig_ext})'
        msg += '\n' + '-' * 74
        return msg

    def make_map_figure(self):
        ts = self.ts
        fig_x = 6 if ts.length > ts.width else 9
        fig_y = max((fig_x - 2) * ts.length * self.aspect / ts.width + 1, 5)
        fig = self.fig_map = _new_figure('Velocity / Cumulative Displacement',
                                         self.fig_size_map or [fig_x, fig_y])
        ax = self.ax_map = fig.add_axes([0.15, 0.15, 0.83, 0.83])
        self.ref_text = fig.text(0.01, 0.99, '', fontsize=8, va='top')
        self.ref_box, = ax.plot([], [], '--k', alpha=0.8, zorder=6)
        self.image = ax.imshow(np.full((ts.length, ts.width), np.nan), cmap=self.cmap,
                               aspect=self.aspect, interpolation='nearest')
        self.cbar = fig.colorbar(self.image, ax=ax, orientation='vertical')
        if not ts.geocoded and not self.noflip:
            orbit = str(ts.atr.get('ORBIT_DIRECTION', '')).lower()
            if orbit.startswith('a'):
                ax.invert_yaxis()             # as MintPy: north up
            elif orbit.startswith('d'):
                ax.invert_xaxis()             # as MintPy: east right

        # what the map can show: the velocity, then the noise indices / layers
        self.layers = {}
        if self.ts2 is not None:
            self.layers['vel(1)'] = (self.vel, 'vel')
            self.layers['vel(2)'] = (self.vel2, 'vel')
        else:
            self.layers['vel'] = (self.vel, 'vel')
        for name, (data, unit) in ts.layers.items():
            self.layers[name] = (data, 'mask' if name == 'mask' else unit)
        rad_ax = fig.add_axes([0.01, 0.3, 0.13, len(self.layers) * 0.025 + 0.04])
        self.radio = RadioButtons(rad_ax, tuple(self.layers))
        for label in self.radio.labels:
            label.set_fontsize(8)
        self.radio.on_clicked(self.show_layer)
        if self.has_mask:
            box = fig.add_axes([0.01, 0.2, 0.1, 0.08])
            self.mask_check = CheckButtons(box, ['mask'], [self.use_mask])
            self.mask_check.on_clicked(self.toggle_mask)

        # the date slider: the cumulative displacement at that date
        self.date_nums = mdates.date2num(ts.datetimes)
        sax = fig.add_axes([0.1, 0.08, 0.8, 0.05], yticks=[])
        self.slider = Slider(sax, 'yr', self.date_nums[0] - 3, self.date_nums[-1] + 3,
                             valinit=self.date_nums[self.ix_ref], valfmt='')
        sax.bar(self.date_nums, np.ones(ts.n_date), facecolor='black', width=4)
        sax.bar(self.date_nums[self.ix_ref], 1, facecolor='red', width=8)
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
        self.draw_ref()
        self.show_layer(next(iter(self.layers)))
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
    def show_layer(self, name):
        """The velocity (relative to the reference area) or another layer."""
        data, kind = self.layers[name]
        keep = self.keep
        cmap, unit = self.cmap, kind
        if kind == 'vel':
            shown = np.where(keep, (data - self.area_mean_2d(data)) * self.f, np.nan)
            limits = self.vlim or auto_limits(shown, self.auto_crange)
            unit = f'{self.unit}/yr'
        elif kind == 'mask':
            shown, limits, cmap, unit = data, (0, 1), 'viridis', ''
        else:
            shown = np.where(keep, data, np.nan)
            limits = auto_limits(shown, self.auto_crange)
            cmap = 'viridis_r'
            if name in ('coh_avg', 'n_unw', 'maxTlen', 'temporalCoherence'):
                cmap = 'viridis'
            elif name == 'mli':
                cmap = 'gray'
            elif name == 'hgt':
                cmap, limits = 'terrain', (-limits[1] / 3, limits[1])
        self.image.set_data(shown)
        self.image.set_cmap(cmap)
        self.image.set_clim(*limits)
        self.cbar.set_label(unit)
        self.ax_map.set_title(name)
        self.cum_mode, self.layer = False, name
        self.fig_map.canvas.draw_idle()
        return

    def cum_limits(self):
        """LiCSBAS: the colour range of the cumulative displacement at the last date."""
        if self._last_cum is None:
            self._last_cum = self.ts.slice(self.ts.n_date - 1) - self.ref_slice
        last = np.where(self.keep, self._last_cum, np.nan)
        last = (last - self.area_mean_2d(last)) * self.f
        return auto_limits(last, self.auto_crange)

    def show_date(self, value):
        """The slider: the cumulative displacement at the nearest date."""
        ts = self.ts
        idx = int(np.argmin(np.abs(self.date_nums - value)))
        disp = np.where(self.keep, ts.slice(idx) - self.ref_slice, np.nan)
        disp = (disp - self.area_mean_2d(disp)) * self.f
        self.image.set_data(disp)
        self.image.set_cmap(self.cmap)
        self.image.set_clim(*(self.dlim or self.cum_limits()))
        self.cbar.set_label(self.unit)
        day = ts.datetimes[idx].strftime('%Y/%m/%d')
        ref = ts.datetimes[self.ix_ref].strftime('%Y/%m/%d')
        self.ax_map.set_title(f'{day} (Ref: {ref})')
        self.cum_mode, self.date_index = True, idx
        self.fig_map.canvas.draw_idle()
        return

    def redraw_map(self):
        if self.cum_mode:
            self.show_date(self.date_nums[self.date_index])
        else:
            self.show_layer(self.layer)

    def draw_ref(self):
        x1, x2, y1, y2 = self.ref
        xs = [x1 - 0.5, x2 - 0.5, x2 - 0.5, x1 - 0.5, x1 - 0.5]
        ys = [y1 - 0.5, y1 - 0.5, y2 - 0.5, y2 - 0.5, y1 - 0.5]
        self.ref_box.set_data(xs, ys)
        note = f'{HELP_TEXT}\n\nRef area:\n X {x1}:{x2}\n Y {y1}:{y2}\n (start from 0)'
        if self.radius:
            note += f'\n\nCircle:\n r = {self.radius:g} km'
        self.ref_text.set_text(note)

    def draw_marker(self):
        """The pixel (a square, as LiCSBAS) or the circle, on the map."""
        for art in self.marker_artists:
            art.remove()
        self.marker_artists = []
        row, col = self.point
        ax = self.ax_map
        if self.radius:
            x_km, y_km = self.ts.pixel_size_km(row)
            ring = patches.Ellipse((col, row), 2 * self.radius / x_km, 2 * self.radius / y_km,
                                   fill=False, ec='k', lw=1.6, zorder=7,
                                   path_effects=[patheffects.withStroke(linewidth=3.4,
                                                                        foreground='w')])
            self.marker_artists.append(ax.add_patch(ring))
        else:
            xs = [col - 0.5, col + 0.5, col + 0.5, col - 0.5, col - 0.5]
            ys = [row - 0.5, row - 0.5, row + 0.5, row + 0.5, row - 0.5]
            self.marker_artists += ax.plot(xs, ys, 'k', linewidth=3, zorder=7)
        self.marker_artists += ax.plot([col], [row], 'Pk', zorder=8)
        self.fig_map.canvas.draw_idle()

    ##---------- events
    def _toolbar_busy(self):
        toolbar = getattr(self.fig_map.canvas, 'toolbar', None)
        return bool(toolbar is not None and getattr(toolbar, 'mode', ''))

    def on_click(self, event):
        """Left double click on the map: the time series there."""
        if (event.inaxes is not self.ax_map or event.button != 1 or not event.dblclick
                or event.xdata is None or self._toolbar_busy()):
            return
        row, col = int(round(event.ydata)), int(round(event.xdata))
        if 0 <= row < self.ts.length and 0 <= col < self.ts.width:
            self.select(row, col)

    def on_ref_select(self, eclick, erelease):
        """Right drag: a new reference area (refused if it has no valid pixel)."""
        ts = self.ts
        xs = sorted([eclick.xdata, erelease.xdata])
        ys = sorted([eclick.ydata, erelease.ydata])
        if None in xs or None in ys:
            return
        x1, x2 = int(np.clip(round(xs[0]), 0, ts.width - 1)), int(np.clip(round(xs[1]), 0, ts.width))
        y1, y2 = int(np.clip(round(ys[0]), 0, ts.length - 1)), int(np.clip(round(ys[1]), 0, ts.length))
        x2, y2 = max(x2, x1 + 1), max(y2, y1 + 1)
        if not self.keep[y1:y2, x1:x2].any():
            print(f'Selected ref {x1}:{x2}/{y1}:{y2} has all nan!! Reselect different ref area.')
            return
        self.ref = (x1, x2, y1, y2)
        self.vprint(f'reference area: {x1}:{x2}/{y1}:{y2}')
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
        self.update_ref_values()
        self.draw_ref()
        self.redraw_map()
        if self.point is not None:
            self.select(*self.point)

    ##---------- the time series
    def select(self, row, col):
        self.point = (int(row), int(col))
        self.draw_marker()
        self.result = self.compute(*self.point)
        self.draw_ts()

    def compute(self, row, col):
        """The time series of a pixel, or the date-by-date statistics of a circle."""
        ts, f = self.ts, self.f
        if self.radius:
            rows, cols = circle_pixels_rc(ts, row, col, self.radius)
        else:
            rows, cols = np.array([row]), np.array([col])
        good = self.keep[rows, cols]
        res = dict(row=row, col=col, n_total=int(rows.size), n_kept=int(good.sum()),
                   series=[], vel=[], gap=None, lalo=self.yx2lalo(row, col), indices={})
        kr, kc = rows[good], cols[good]
        pick = (kr, kc) if kr.size else (rows, cols)
        for name, (data, unit) in ts.layers.items():
            if name != 'mask':
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    res['indices'][name] = (float(np.nanmedian(data[pick])), unit)
        if not kr.size:
            return res
        rel = relative_series(ts.pixels(kr, kc), self.ref1, self.ix_ref) * f
        res['series'].append(dict(series_stats(rel), label=ts.label))
        res['vel'].append(float(np.nanmedian(self.vel[kr, kc] - self.area_mean_2d(self.vel)) * f))
        if self.ts2 is not None:
            rel2 = relative_series(self.align2(self.ts2.pixels(kr, kc)), self.ref2,
                                   self.ix_ref) * f
            res['series'].append(dict(series_stats(rel2), label=self.ts2.label))
            res['vel'].append(float(np.nanmedian(self.vel2[kr, kc]
                                                 - self.area_mean_2d(self.vel2)) * f))
        res['gap'] = ts.gap_fraction(kr, kc)
        return res

    def info_text(self, res):
        """The text beside the time series, as LiCSBAS."""
        x1, x2, y1, y2 = self.ref
        text = ''
        if res['lalo'] is not None:
            lat, lon = res['lalo']
            text = (f'N:{lat:.1f}\nE:{lon:.1f}\n\n' if self.ts.projected
                    else f'Lat:{lat:.5f}\nLon:{lon:.5f}\n\n')
        text += (f'Ref area:\n X {x1}:{x2}\n Y {y1}:{y2}\n (start from 0)\n'
                 f'Ref date:\n {self.ts.dates[self.ix_ref]}\n\n')
        if self.radius:
            text += f'Circle: r = {self.radius:g} km\n {res["n_total"]} px, {res["n_kept"]} used\n'
            text += '(medians below)\n'
        elif self.has_mask:
            text += f'mask: {int(bool(self.ts.keep[res["row"], res["col"]]))}\n'
        for name, (value, unit) in res['indices'].items():
            integer = name.startswith('n_')
            text += f'{name}: {_fmt(value, integer)} {unit}'.rstrip() + '\n'
        return text

    def draw_ts(self):
        res, ax, ts = self.result, self.ax_ts, self.ts
        ax.cla()
        ax.grid(zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlabel('Time')
        ax.set_ylabel(f'Displacement ({self.unit})')
        self.side_text.set_text(self.info_text(res))
        self.fit_lines = []
        where = f'@({res["col"]}, {res["row"]})'
        if not res['series']:
            ax.set_title(f'NaN {where}' + (f' - no valid pixel within {self.radius:g} km'
                                           if self.radius else ''), fontsize=10)
            self.fig_ts.canvas.draw_idle()
            return

        locator = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        dates = np.array(ts.datetimes, dtype=object)
        n_fit = max(int((ts.datetimes[-1] - ts.datetimes[0]).days // 10), 1)
        fit_dates = [ts.datetimes[0] + dt.timedelta(days=10 * i) for i in range(n_fit + 1)]
        fit_years = np.arange(n_fit + 1) * 10 / 365.25
        for k, stats in enumerate(res['series']):
            color = COLORS[k]
            lines = []
            for model in range(len(MODELS)):
                curve = model_curve(ts.years, stats['median'], fit_years, model)
                lines += ax.plot(fit_dates, curve, f'{color}-', alpha=0.6, zorder=3 - k,
                                 visible=self.fit_visible[model])
            self.fit_lines.append(lines)
            if self.radius:
                med = stats['median']
                yerr = np.vstack([med - stats['q25'], stats['q75'] - med])
                ax.errorbar(dates, med, yerr=yerr, fmt='o', ms=4, color=color, ecolor=color,
                            elinewidth=0.8, capsize=2, alpha=0.6, zorder=5 - k,
                            label=f'{stats["label"]}\n(median, 25-75 %)')
            else:
                ax.scatter(dates, stats['median'], c=color, alpha=0.6, zorder=5 - k,
                           label=stats['label'])
        if res['gap'] is not None:
            gaps = res['gap'] >= (GAP_SHARE if self.radius else 1.0)
            if gaps.any():
                mid = dates[:-1][gaps] + (dates[1:][gaps] - dates[:-1][gaps]) / 2
                ax.vlines(list(mid), 0, 1, transform=ax.get_xaxis_transform(), zorder=1,
                          alpha=0.6, colors='k', label='Gap of ifg network')
        unit = f'{self.unit}/yr'
        if len(res['vel']) == 2:
            title = (f'vel(1) = {res["vel"][0]:.1f} {unit}, vel(2) = {res["vel"][1]:.1f} '
                     f'{unit} {where}')
        else:
            title = f'vel = {res["vel"][0]:.1f} {unit} {where}'
        if self.radius:
            title += f'\nmedian of {res["n_kept"]} pixels within {self.radius:g} km'
        ax.set_title(title, fontsize=10)
        if self.ylim:
            ax.set_ylim(self.ylim)
        ax.legend(fontsize=8)
        self.fig_ts.canvas.draw_idle()
        return

    ##---------- output
    @property
    def outfile_base(self):
        row, col = self.point
        name = f'tseries_Y{row}X{col}' + (f'_r{self.radius:g}km' if self.radius else '')
        return os.path.join(self.outdir, name)

    def save_outputs(self):
        """The time series as .csv, how it was made (_info.txt), and both figures."""
        res = self.result
        if res is None or not res['series']:
            print('nothing to save - double-click a point with data on the map.')
            return []
        base = self.outfile_base
        os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
        ts, unit, written = self.ts, self.unit, []

        # 1. the time series
        header = ['date', 'year']
        many = len(res['series']) > 1
        for k in range(len(res['series'])):
            tag = f'_{k + 1}' if many else ''
            if self.radius:
                header += [f'median_{unit}{tag}', f'q25_{unit}{tag}', f'q75_{unit}{tag}',
                           f'n_px{tag}']
            else:
                header += [f'disp_{unit}{tag}']
        out = f'{base}.csv'
        with open(out, 'w') as fid:
            fid.write(','.join(header) + '\n')
            for i, (date, day) in enumerate(zip(ts.dates, ts.datetimes)):
                start = dt.datetime(day.year, 1, 1)
                length = (dt.datetime(day.year + 1, 1, 1) - start).days
                cells = [date, f'{day.year + (day - start).days / length:.4f}']
                for stats in res['series']:
                    if self.radius:
                        cells += [f'{stats["median"][i]:.4f}', f'{stats["q25"][i]:.4f}',
                                  f'{stats["q75"][i]:.4f}', f'{int(stats["count"][i])}']
                    else:
                        cells += [f'{stats["median"][i]:.4f}']
                fid.write(','.join(cells) + '\n')
        written.append(out)

        # 2. how it was made
        out = f'{base}_info.txt'
        x1, x2, y1, y2 = self.ref
        with open(out, 'w') as fid:
            fid.write(f'file           : {os.path.abspath(self.file)} ({ts.kind})\n')
            if self.ts2 is not None:
                fid.write(f'second file    : {os.path.abspath(self.file2)}\n')
            fid.write(f'point          : x/y {res["col"]}/{res["row"]} (start from 0)\n')
            if res['lalo'] is not None:
                tag = 'northing/easting' if ts.projected else 'lat/lon'
                fid.write(f'{tag:<15}: {res["lalo"][0]:.6f}, {res["lalo"][1]:.6f}\n')
            if self.radius:
                fid.write(f'circle         : r = {self.radius:g} km, {res["n_total"]} pixels, '
                          f'{res["n_kept"]} used (with data, in the mask)\n')
                fid.write('statistics     : per date, the median of the pixels and their 25th '
                          'and 75th percentiles (the error bars)\n')
            fid.write(f'reference area : {x1}:{x2}/{y1}:{y2} (x/y, start from 0), '
                      'its mean subtracted at every date\n')
            fid.write(f'reference date : {ts.dates[self.ix_ref]} (every pixel 0 there)\n')
            fid.write(f'mask           : {"on" if self.use_mask else "off"}'
                      + (f' ({ts.mask_source})' if self.has_mask else ' (none)') + '\n')
            for k, vel in enumerate(res['vel']):
                fid.write(f'velocity{"" if len(res["vel"]) == 1 else f"({k + 1})":<7}: '
                          f'{vel:.3f} {unit}/yr'
                          + (' (median of the circle)' if self.radius else '') + '\n')
            fid.write(f'velocity from  : {ts.vel_source}\n')
            fid.write(f'unit           : {unit}\n')
        written.append(out)

        # 3. the figures
        out = self.outfile[0] if self.outfile else f'{base}{self.fig_ext}'
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
        for ts in (self.ts, self.ts2):
            if ts is not None:
                ts.close()
        return

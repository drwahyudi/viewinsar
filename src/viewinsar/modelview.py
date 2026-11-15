"""`viewinsar model`: a swath across a fault, fitted with an elastic dislocation model.

The map and the profile line work as in `viewinsar swath`. Every section gets
a window with the pixels of its swath as points - fault-parallel velocity
against the distance from the fault - their binned medians, the GNSS, and the
maximum a posteriori (MAP) model with its 95 % range and the parameters
written on the figure, after Hussain et al. (2016, 2018) and Ou et al. (2022);
and a second window with the posterior distributions of the parameters.

    left click x2  the profile, across the fault     "="  another section
    "n"            redo the current section          "e"  export everything
"""

import os
import re
import textwrap

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText

from .model import NAMES, UNITS, describe, fit_profile, model_band, strain_rate
from .radar import PIXEL_ATR
from .raster import same_crs
from .swath import EARTH_RADIUS, is_projected, line_geometry, project_points
from .viewer import (
    GNSS_EDGE,
    NO_LOOKUP,
    Section,
    _new_figure,
    export_dir,
    per_file,
    swathProfileViewer,
    velocity_factor,
    white_halo,
)

__all__ = ['ModelViewer', 'ModelSection', 'azimuth_of', 'local_en', 'strike_of_points',
           'orient_strike', 'slip_sense', 'fault_normal']

POINT_COLORS = ('#d62728', '#9467bd', '#2ca02c', '#ff7f0e', '#8c564b', '#e377c2',
                '#bcbd22', '#17becf')            # the pixels of each file, red first
GNSS_COLOR = '#1f77b4'                           # blue, as Hussain et al.
MAP_COLOR = 'k'                                  # the MAP model: bold black dashes
BAND_COLOR = '0.62'                              # the 95 % range: light grey
POST_COLOR = '#d62728'                           # the MAP in the posterior plots
HIST_COLOR, HIST_EDGE = '#1f4fd8', '#0a1f73'      # the marginal histograms, blue
# names that mark a component of a decomposition, not a line of sight
DECOMPOSED = {'hz', 'horz', 'horizontal', 'up', 'vert', 'vertical', 'east', 'north', 'ew',
              'ns', 've', 'vn', 'vu', 'decomp', 'decomposed'}
# the GNSS component that matches a file named as one (horizontal: unknown direction)
DECOMPOSED_GNSS = {'up': 'vu', 'vert': 'vu', 'vertical': 'vu', 'vu': 'vu', 'east': 've',
                   'ew': 've', 've': 've', 'north': 'vn', 'ns': 'vn', 'vn': 'vn'}
POSTERIOR_ORDER = ('S', 'd1', 'a', 'C', 'd2', 'x0')   # as Hussain et al. (2016)
POSTERIOR_LABEL = {'S': 'Slip rate', 'd1': 'Locking depth', 'a': 'Offset',
                   'C': 'Creep rate', 'd2': 'Creep depth', 'x0': 'Fault position'}
KM_PER_DEG = np.pi / 180.0 * EARTH_RADIUS / 1000.0
PHYSICAL = ('S', 'd1', 'C', 'd2', 'x0')          # what the posterior figure shows
SYMBOL = {'S': '$S$', 'd1': '$d_1$', 'C': '$C$', 'd2': '$d_2$', 'x0': '$x_0$',
          'a': '$a$', 'ag': '$a_{GNSS}$'}
COMP_TEXT = {'as-is': 'the data as given', 'par': 'fault-parallel data',
             've': 'from the east velocity', 'vn': 'from the north velocity',
             'los': 'from the LOS velocity'}
COMP_LONG = {'as-is': 'the data as given, no conversion',
             'par': 'fault-parallel velocity, as given',
             've': 'east velocity / sin(strike)', 'vn': 'north velocity / cos(strike)',
             'los': 'LOS velocity / (LOS unit vector . strike direction), pixel by pixel'}


##############################################################################
#                               geometry                                     #
##############################################################################
def azimuth_of(d_east, d_north):
    """Azimuth of a vector, degrees clockwise from north, in [0, 360)."""
    return float(np.degrees(np.arctan2(d_east, d_north)) % 360.0)


def local_en(lon, lat, lon0, lat0):
    """km east and north of (lon0, lat0) - local equirectangular, as the swath."""
    lon, lat = np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)
    dlon = (lon - lon0 + 180.0) % 360.0 - 180.0
    return dlon * np.cos(np.radians(lat0)) * KM_PER_DEG, (lat - lat0) * KM_PER_DEG


def strike_of_points(east, north):
    """The azimuth (0-180) of the principal direction of points (a fault trace)."""
    pts = np.c_[np.asarray(east, dtype=np.float64), np.asarray(north, dtype=np.float64)]
    pts = pts[np.all(np.isfinite(pts), axis=1)]
    if pts.shape[0] < 2:
        return None
    pts = pts - pts.mean(axis=0)
    vec = np.linalg.svd(pts, full_matrices=False)[2][0]
    return azimuth_of(vec[0], vec[1]) % 180.0


def orient_strike(strike, profile_az):
    """Of the two directions of a strike, the one for which S > 0 is
    right-lateral: the strike points to the right of the profile direction
    (sin(profile_az - strike) < 0)."""
    strike = float(strike) % 360.0
    if np.sin(np.radians(profile_az - strike)) > 0:
        strike = (strike + 180.0) % 360.0
    return strike


def slip_sense(S, strike, profile_az):
    """'right-lateral' or 'left-lateral' for a slip rate S: the velocity
    along `strike` rises by S across the fault, from the start of the profile
    to its end."""
    side = S * np.sin(np.radians(profile_az - strike))
    return 'left-lateral' if side > 0 else 'right-lateral'


def fault_normal(strike, profile_az):
    """The horizontal unit normal of the fault (east, north), pointing to the
    end of the profile, and the angle (deg) between it and the profile."""
    nrm = np.array([np.cos(np.radians(strike)), -np.sin(np.radians(strike))])
    prof = np.array([np.sin(np.radians(profile_az)), np.cos(np.radians(profile_az))])
    if nrm @ prof < 0:
        nrm = -nrm
    return nrm, float(np.degrees(np.arccos(np.clip(nrm @ prof, -1.0, 1.0))))


def _fmt(value):
    if not np.isfinite(value):
        return 'nan'
    return f'{value:.2f}' if abs(value) < 1 else f'{value:.1f}'


def _num(value):
    """_fmt with a true minus sign, so a range reads (−10.7 – −2.7)."""
    return _fmt(value).replace('-', '\u2212')


def _unit_text(unit):
    """' mm/yr', or '' for data without a unit."""
    unit = str(unit or '').strip()
    return '' if unit.lower() in ('', '1', 'none', 'value') else f' {unit}'


def strain_unit(model_unit):
    """The unit of dv/dx: µstrain/yr for mm/yr over km."""
    if model_unit == 'mm/yr':
        return 'µstrain/yr'
    unit = _unit_text(model_unit).strip()
    return f'{unit}/km' if unit else '1/km'


class _Placer:
    """Puts boxes (the parameters, the legend) inside an axes where they hide
    the fewest weighted points, clear of the A / A' labels and of each other -
    making room above or below the data when every spot would hide some."""

    SPOTS = (('upper left', 0.012, None), ('upper right', 0.988, None),
             ('lower left', 0.012, 0.015), ('lower right', 0.988, 0.015),
             ('upper center', 0.5, 0.985), ('lower center', 0.5, 0.015),
             ('center left', 0.012, 0.5), ('center right', 0.988, 0.5))

    def __init__(self, ax, xy, weights, labels=()):
        self.ax = ax
        self.data_xy = np.asarray(xy, dtype=np.float64)
        self.data_w = np.asarray(weights, dtype=np.float64)
        self.renderer = None
        try:
            self.renderer = ax.figure.canvas.get_renderer()
        except Exception:                         # noqa: BLE001 - a non-Agg canvas
            pass
        self.refresh()
        self.taken = [r for r in (self._rect_of(a) for a in labels) if r is not None]
        tops = [r[1] for r in self.taken]
        self.top = (min(tops) - 0.03) if tops else 0.86

    def refresh(self):
        """The points in axes fraction, after a change of the axis limits."""
        trans = self.ax.transData + self.ax.transAxes.inverted()
        xy = trans.transform(self.data_xy)
        good = np.all(np.isfinite(xy), axis=1)
        self.xy, self.w = xy[good], self.data_w[good]

    def _rect_of(self, artist):
        """(x0, y0, x1, y1) of an artist, in axes fraction."""
        if self.renderer is None:
            return None
        try:
            ext = artist.get_window_extent(self.renderer)
            box = self.ax.get_window_extent(self.renderer)
        except Exception:                         # noqa: BLE001
            return None
        return ((ext.x0 - box.x0) / box.width, (ext.y0 - box.y0) / box.height,
                (ext.x1 - box.x0) / box.width, (ext.y1 - box.y0) / box.height)

    def _cost(self, rect):
        x0, y0, x1, y1 = rect
        if x0 < -1e-3 or y0 < -1e-3 or x1 > 1.001 or y1 > 1.001:
            return np.inf
        for t in self.taken:
            if x0 < t[2] and x1 > t[0] and y0 < t[3] and y1 > t[1]:
                return np.inf
        inside = ((self.xy[:, 0] >= x0) & (self.xy[:, 0] <= x1)
                  & (self.xy[:, 1] >= y0) & (self.xy[:, 1] <= y1))
        return float(self.w[inside].sum())

    def _best(self, make):
        artist = make('upper left', (0.012, self.top))
        rect = self._rect_of(artist)
        if rect is None:
            return artist, None
        width, height = rect[2] - rect[0], rect[3] - rect[1]
        best = None
        for loc, ax_x, ax_y in self.SPOTS:
            ax_y = self.top if ax_y is None else ax_y
            x0 = {'left': ax_x, 'right': ax_x - width}.get(loc.split()[-1], ax_x - width / 2)
            y0 = {'upper': ax_y - height, 'lower': ax_y}.get(loc.split()[0], ax_y - height / 2)
            cand = (x0, y0, x0 + width, y0 + height)
            cost = self._cost(cand)
            if best is None or cost < best[0]:
                best = (cost, loc, (ax_x, ax_y), cand)
        artist.remove()
        return make(best[1], best[2]), best

    def place(self, make, room=True, budget=0.005, max_grow=0.4):
        """Draw a box - make(loc, anchor) - at the best spot. With room, a box
        that would hide more than `budget` of the points makes room for itself
        above (or below) the data: at most max_grow of the y-range, so the
        data keep most of the panel."""
        artist, best = self._best(make)
        if best is None:
            self.taken.append((0.0, self.top - 0.3, 0.5, self.top))
            return artist
        vertical = best[1].split()[0]
        if (room and np.isfinite(best[0]) and best[0] > budget * self.w.sum()
                and vertical in ('upper', 'lower')):
            cand = best[3]
            height = cand[3] - cand[1] + ((1 - self.top) if vertical == 'upper' else 0.015)
            lo, hi = self.ax.get_ylim()
            grow = min((hi - lo) * height / max(1.0 - height, 0.4) * 1.05, max_grow * (hi - lo))
            self.ax.set_ylim((lo, hi + grow) if vertical == 'upper' else (lo - grow, hi))
            artist.remove()
            self.refresh()
            artist, best = self._best(make)
        self.taken.append(best[3])
        return artist


##############################################################################
class ModelSection(Section):
    """One profile across a fault, its fit and its two windows."""

    def __init__(self, viewer, index):
        super().__init__(viewer, index)
        self.ax_strain = None
        self.fig_post = None
        self.map_marks = []
        self.clear_model()

    def clear_model(self):
        self.result = None
        self.fp = []                  # per file: x, v (fault-parallel), ...
        self.frame = None             # the local east/north frame, at the fault
        self.x_fault = None           # km along the profile
        self.fault_how = ''
        self.prof_az = self.strike = None
        self.strike_how = ''
        self.normal = None
        self.obliquity = 0.0
        self.sense_known = False
        self.along_scale = 1.0        # km from the fault per km along the profile
        self.gnss_fit_shift = 0.0

    def clear_data(self):
        super().clear_data()
        self.clear_model()

    @property
    def has_data(self):
        return self.result is not None

    ##---------- the windows
    def create_figure(self):
        v = self.v
        if self.fig is not None and plt.fignum_exists(self.fig.number):
            return self.fig
        panels = [3.2] + ([1.2] if v.show_strain else []) + ([1.1] if v.topoxsect else [])
        height = v.fig_size[1] + 1.5 * (len(panels) - 1)
        self.fig = _new_figure(f'Model {self.name}', [v.fig_size[0], height])
        gsp = self.fig.add_gridspec(len(panels), 1, height_ratios=panels, hspace=0.08,
                                    left=0.12, right=0.97, bottom=0.1, top=0.9)
        self.ax_swath = self.fig.add_subplot(gsp[0])
        k = 1
        self.ax_strain = self.ax_topo = None
        if v.show_strain:
            self.ax_strain = self.fig.add_subplot(gsp[k], sharex=self.ax_swath)
            k += 1
        if v.topoxsect:
            self.ax_topo = self.fig.add_subplot(gsp[k], sharex=self.ax_swath)
        self.fig.canvas.mpl_connect('key_press_event', v.on_key)
        self.reset_axes()
        self.ax_swath.text(0.5, 0.5, f'section {self.name}\nclick 2 points on the map,\n'
                           'across the fault', ha='center', va='center', color='0.55',
                           fontsize=v.font_size, transform=self.ax_swath.transAxes)
        return self.fig

    def reset_axes(self):
        v = self.v
        axes = [ax for ax in (self.ax_swath, self.ax_strain, self.ax_topo) if ax is not None]
        for ax in axes:
            ax.cla()
            ax.tick_params(which='both', direction='in', labelsize=v.font_size,
                           bottom=True, top=True, left=True, right=True)
        for ax in axes[:-1]:
            ax.tick_params(labelbottom=False)
        axes[-1].set_xlabel('Distance from the fault [km]', fontsize=v.font_size + 1)
        what = 'Velocity' if v.comp in ('as-is', 'auto') else 'Fault-parallel velocity'
        unit = _unit_text(v.model_unit).strip()
        self.ax_swath.set_ylabel(what + (f' [{unit}]' if unit else ''), fontsize=v.font_size + 1)
        if self.ax_strain is not None:
            self.ax_strain.set_ylabel(f'dv/dx\n[{strain_unit(v.model_unit)}]',
                                      fontsize=v.font_size)
        if self.ax_topo is not None:
            self.ax_topo.set_ylabel('Elevation\n[m]', fontsize=v.font_size)
        return

    def clear_map_artists(self):
        super().clear_map_artists()
        for art in self.map_marks:
            try:
                art.remove()
            except (ValueError, NotImplementedError):
                pass
        self.map_marks = []

    ##---------- the data
    def extract(self):
        if not super().extract():
            return False
        try:
            self.prepare()
        except ValueError as err:
            print(f'WARNING: section {self.name}: {err}')
            return False
        return True

    def extract_gnss(self):
        """The GNSS inside the swath (--gnss-width, flat caps); their values
        follow once the strike is known (prepare)."""
        v = self.v
        self.gnss_in = None
        if v.gnss is None:
            return
        geom = line_geometry(v.atr, self.start_yx, self.end_yx, v.width)
        along, perp = project_points(geom, v.gnss_rows, v.gnss_cols)
        length = geom['length'] * geom['dist_scale']
        width = v.gnss_width or v.width
        inside = (np.isfinite(along) & np.isfinite(perp) & (np.abs(perp) <= width)
                  & (along >= 0) & (along <= length))
        idx = np.flatnonzero(inside)
        self.gnss_in = dict(index=idx, dist=along[idx], perp=perp[idx],
                            site=[v.gnss.site[i] for i in idx], is3d=v.gnss.has_vu[idx])
        return

    def prepare(self):
        """Where the fault is, its strike, and every pixel as (distance from
        the fault, fault-parallel velocity)."""
        v = self.v
        geom = line_geometry(v.atr, self.start_yx, self.end_yx, v.width)
        length = geom['length'] * geom['dist_scale']
        if geom['dist_unit'] != 'km':
            raise ValueError('the pixel size is unknown, so are the distances in km')

        # 1. the fault along the profile
        if v.x0_km is not None:
            self.x_fault, self.fault_how = float(v.x0_km), '--x0'
        elif self.crossings.size:
            k = int(np.argmin(np.abs(self.crossings - length / 2)))
            self.x_fault = float(self.crossings[k])
            self.fault_how = 'where the profile crosses the fault trace'
            if self.crossings.size > 1:
                print(f'   section {self.name} crosses {self.crossings.size} fault traces, at '
                      + ', '.join(f'{c:.2f}' for c in self.crossings) + ' km: modelled at '
                      f'{self.x_fault:.2f} km (the nearest the middle; --x0 picks another)')
        else:
            self.x_fault = length / 2
            self.fault_how = 'the middle of the profile (no --x0, no fault crossing)'
        if not 0 <= self.x_fault <= length:
            print(f'WARNING: section {self.name}: the fault, at {self.x_fault:.2f} km, is '
                  f'beyond the profile (0 - {length:.2f} km)')

        # 2. the local east / north frame at the fault, the profile azimuth
        t = self.x_fault / length
        row_c = self.start_yx[0] + t * (self.end_yx[0] - self.start_yx[0])
        col_c = self.start_yx[1] + t * (self.end_yx[1] - self.start_yx[1])
        self.frame = v.make_frame(row_c, col_c)
        e01, n01, _, _ = v.en_km(self.frame, 0, [self.start_yx[0], self.end_yx[0]],
                                 [self.start_yx[1], self.end_yx[1]])
        if not np.all(np.isfinite(np.r_[e01, n01])):
            raise ValueError('the ends of the profile have no lon/lat')
        self.prof_az = azimuth_of(e01[1] - e01[0], n01[1] - n01[0])

        # 3. the strike: given, from the trace, or perpendicular to the profile
        strike, self.strike_how = None, ''
        if v.strike is not None:
            strike, self.strike_how = float(v.strike) % 360.0, '--strike'
        elif v.faults:
            strike = v.trace_strike(self.frame)
            if strike is not None:
                self.strike_how = 'the fault trace at the fault'
        if strike is None:
            strike = (self.prof_az + 90.0) % 360.0
            self.strike_how = 'perpendicular to the profile'
        given_sign = v.comp in ('as-is', 'par') and v.strike is not None
        self.strike = strike if given_sign else orient_strike(strike, self.prof_az)
        self.sense_known = v.comp in ('los', 've', 'vn') or (v.comp == 'par'
                                                             and v.strike is not None)
        self.normal, self.obliquity = fault_normal(self.strike, self.prof_az)
        prof = np.array([np.sin(np.radians(self.prof_az)), np.cos(np.radians(self.prof_az))])
        self.along_scale = float(self.normal @ prof)
        if self.obliquity > 30:
            print(f'NOTE: section {self.name} runs {self.obliquity:.0f} deg off the fault '
                  'normal: distances are taken perpendicular to the fault (a straight line '
                  'through the fault position, along the strike)')

        # 4. every pixel: distance from the fault, fault-parallel velocity
        p_vec = (np.sin(np.radians(self.strike)), np.cos(np.radians(self.strike)))
        self.fp = []
        for sw in self.swath_list:
            i = sw['file_idx']
            east, north, lon, lat = v.en_km(self.frame, i, sw['row'], sw['col'])
            x = east * self.normal[0] + north * self.normal[1]
            value = np.asarray(sw['value'], dtype=np.float64)
            if v.comp == 'los':
                ue, un, _ = v.los_vectors(i, sw['row'], sw['col'], lon, lat)
                proj = ue * p_vec[0] + un * p_vec[1]
            elif v.comp == 've':
                proj = np.full(value.size, p_vec[0])
            elif v.comp == 'vn':
                proj = np.full(value.size, p_vec[1])
            else:
                proj = np.ones(value.size)
            weak = ~np.isfinite(proj) | (np.abs(proj) < 0.05)
            with np.errstate(divide='ignore', invalid='ignore'):
                vpar = np.where(weak, np.nan, value * v.to_model / proj)
            good = np.isfinite(x) & np.isfinite(vpar)
            med_proj = float(np.nanmedian(np.abs(proj))) if np.isfinite(proj).any() else np.nan
            if v.comp in ('los', 've', 'vn') and med_proj < 0.25:
                print(f'WARNING: {v.labels[i]}: the {v.comp.upper()} sees only '
                      f'{med_proj:.2f} of the fault-parallel motion (strike '
                      f'{self.strike:.0f} deg): its noise is amplified {1 / med_proj:.1f}x')
            self.fp.append(dict(file_idx=i, x=x[good], v=vpar[good], dist=sw['dist'][good],
                                perp=sw['perp'][good], value=value[good],
                                row=np.asarray(sw['row'])[good], col=np.asarray(sw['col'])[good],
                                lon=np.asarray(lon)[good], lat=np.asarray(lat)[good],
                                proj=med_proj))
        self.fp = [f for f in self.fp if f['x'].size]
        if not self.fp:
            raise ValueError('no pixel left after the conversion to fault-parallel velocity')

        # 5. the GNSS in the swath: the same component, the same distance
        gin = self.gnss_in
        if gin is not None and gin['index'].size:
            idx = gin['index']
            lon_g, lat_g = v.gnss.lon[idx], v.gnss.lat[idx]
            if self.frame['kind'] == 'lonlat':
                east, north = local_en(lon_g, lat_g, self.frame['lon0'], self.frame['lat0'])
            else:
                east = north = np.full(idx.size, np.nan)
            gin['x'] = east * self.normal[0] + north * self.normal[1]
            value, sigma = v.gnss_values(v.gnss.subset(idx), p_vec)
            gin['value'], gin['sigma'] = value, sigma
        return

    ##---------- the fit
    def fit(self):
        v = self.v
        gin = self.gnss_in
        gnss = None
        if gin is not None and gin['index'].size and v.gnss_fit:
            ok = np.isfinite(gin['x']) & np.isfinite(gin['value'])
            if ok.any():
                gnss = dict(x=gin['x'][ok], v=gin['value'][ok], sigma=gin['sigma'][ok])
        bounds = {'S': (-v.s_max_used, v.s_max_used), 'd1': (0.1, v.d1_max),
                  'C': (0.0, v.c_max), 'd2': (0.01, v.d2_max)}
        free_x0 = (-v.free_x0, v.free_x0) if v.free_x0 else None
        print(f'section {self.name}: fault at {self.x_fault:.2f} km ({self.fault_how}), strike '
              f'{self.strike:.0f} deg ({self.strike_how}), profile azimuth {self.prof_az:.0f}'
              f' deg; {COMP_LONG[v.comp]}')
        self.result = fit_profile([(f['x'], f['v']) for f in self.fp], bin_km=v.bin_km,
                                  gnss=gnss, noise=v.noise, walkers=v.walkers, steps=v.steps,
                                  burn=v.burn, seed=v.seed, creep=v.creep, free_x0=free_x0,
                                  rot=v.rot, rot_side=1 if v.rot_side == 'end' else -1,
                                  bounds=bounds, d1_prior=v.d1_prior,
                                  verbose=v.print_msg)
        res = self.result
        if self.sense_known:
            S = res['summary']['S']['map']
            print(f'   S {S:+.2f} {v.model_unit}: {slip_sense(S, self.strike, self.prof_az)}')
        # the GNSS left out of the fit still sit on the model: their median misfit
        self.gnss_fit_shift = 0.0
        if gin is not None and gin['index'].size and gnss is None:
            ok = np.isfinite(gin['x']) & np.isfinite(gin['value'])
            if ok.any():
                pred = res['model'].predict(res['best'], gin['x'][ok])
                self.gnss_fit_shift = float(np.median(pred - gin['value'][ok]))
        return res

    def offset_shift(self, set_index):
        """What brings file `set_index` (of the fit) onto the offset of the first."""
        res = self.result
        model, best = res['model'], res['best']
        return best[model.index('a')] - best[model.index(model.offsets[set_index])]

    def gnss_shift(self):
        res = self.result
        model = res['model']
        if model.gnss:
            return res['best'][model.index('a')] - res['best'][model.index('ag')]
        return self.gnss_fit_shift

    ##---------- draw
    def draw(self):
        v = self.v
        self.create_figure()
        self.result = None
        if not self.extract():
            return
        try:
            self.fit()
        except ValueError as err:
            print(f'WARNING: section {self.name}: no model - {err}')
            self.result = None
            self.reset_axes()
            self.ax_swath.text(0.5, 0.5, '\n'.join(textwrap.wrap(f'no model: {err}', 60)),
                               ha='center', va='center', color='0.4', fontsize=v.font_size,
                               transform=self.ax_swath.transAxes)
            self.fig.canvas.draw_idle()
            return
        self.reset_axes()
        self.draw_fit()
        if self.ax_strain is not None:
            self.draw_strain()
        if self.ax_topo is not None:
            self.draw_topo()
        self.draw_posterior()
        self.mark_fault_on_map()
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return

    def plot_points(self):
        """(x, v) of the points to draw, per file of the fit - thinned to
        --max-points in all."""
        res = self.result
        total = sum(self.fp[k]['x'].size for k in res['sets'])
        keep = 1.0
        if self.v.max_points and total > self.v.max_points:
            keep = self.v.max_points / total
        rng = np.random.default_rng(0)
        out = []
        for j, k in enumerate(res['sets']):
            f = self.fp[k]
            sel = (rng.random(f['x'].size) < keep) if keep < 1 else slice(None)
            out.append((k, f['x'][sel], f['v'][sel] + self.offset_shift(j)))
        return out, total

    def draw_fit(self):
        v, ax, res = self.v, self.ax_swath, self.result
        model, summ = res['model'], res['summary']
        handles = []

        # the pixels, every file brought onto the offset of the first
        pts, total = self.plot_points()
        n_drawn = sum(p[1].size for p in pts)
        alpha = float(np.clip(1.5e4 / max(n_drawn, 1), 0.05, 0.6))
        multi = len(pts) > 1
        for k, x, y in pts:
            color = POINT_COLORS[k % len(POINT_COLORS)]
            ax.scatter(x, y, s=2.5, c=color, alpha=alpha, lw=0, zorder=1, rasterized=True)
            name = v.labels[self.fp[k]['file_idx']] if multi else 'InSAR'
            handles.append(Line2D([], [], ls='', marker='o', ms=4, color=color, label=name))

        # the model: its 95 % range and the MAP
        lo_x = min(float(x.min()) for _, x, _ in pts)
        hi_x = max(float(x.max()) for _, x, _ in pts)
        gin = self.gnss_in
        has_gnss = (gin is not None and gin['index'].size > 0
                    and bool(np.isfinite(gin['value']).any()))
        if has_gnss:
            lo_x = min(lo_x, float(np.nanmin(gin['x'])))
            hi_x = max(hi_x, float(np.nanmax(gin['x'])))
        xx = np.linspace(lo_x, hi_x, 801)
        x0 = summ['x0']['map'] if 'x0' in summ else 0.0
        xx = np.unique(np.r_[xx, x0 - 1e-6, x0 + 1e-6]) if model.creep else xx
        best, band_lo, band_hi = model_band(res, xx)

        # the binned medians, +/- one robust standard deviation
        bins_xy = []
        for j, k in enumerate(res['sets']):
            b = res['bins'][j]
            face = 'w' if not multi else POINT_COLORS[k % len(POINT_COLORS)]
            ax.errorbar(b['x'], b['median'] + self.offset_shift(j), yerr=b['std'], fmt='o',
                        ms=3.2, mfc=face, mec='0.1', mew=0.7, ecolor='0.25', elinewidth=0.5,
                        capsize=0, zorder=4)
            bins_xy.append(np.c_[b['x'], b['median'] + self.offset_shift(j)])
        handles.append(Line2D([], [], ls='', marker='o', ms=4.5, mfc='w', mec='0.1', mew=0.7,
                              label=f'Median per {res["bin_km"]:g} km (±1σ)'))

        ax.fill_between(xx, band_lo, band_hi, color=BAND_COLOR, alpha=0.55, lw=0, zorder=2)
        ax.plot(xx, best, ls='--', color=MAP_COLOR, lw=2.2, zorder=6)
        handles.append(Line2D([], [], ls='--', color=MAP_COLOR, lw=2.2, label='MAP model'))
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=BAND_COLOR, alpha=0.55,
                                     label='95% range'))

        # the GNSS
        if has_gnss:
            shift = self.gnss_shift()
            white_halo(ax.errorbar(gin['x'], gin['value'] + shift, yerr=gin['sigma'], fmt='o',
                                   ms=6.5, mfc=GNSS_COLOR, ecolor=GNSS_COLOR, elinewidth=1.4,
                                   capsize=3, capthick=1.2, zorder=7, **GNSS_EDGE))
            handles.append(Line2D([], [], ls='', marker='o', ms=6, mfc=GNSS_COLOR,
                                  label='GNSS' + ('' if model.gnss else ' (not fitted)'),
                                  **GNSS_EDGE))

        # the fault (dotted), and the traces the profile crosses (dash-dot)
        free = 'x0' in summ
        ax.axvline(x0, color='k', ls=':', lw=1.2, zorder=3)
        handles.append(Line2D([], [], color='k', ls=':', lw=1.2,
                              label='Fault, fitted position' if free else 'Fault'))
        traces = [(c - self.x_fault) * self.along_scale for c in self.crossings]
        traces = [xc for xc in traces if lo_x <= xc <= hi_x and (free or abs(xc) > 0.05)]
        for xc in traces:
            ax.axvline(xc, color=v.fault_color, ls='-.', lw=0.9, alpha=0.75, zorder=3)
        if traces:
            handles.append(Line2D([], [], color=v.fault_color, ls='-.', lw=0.9,
                                  label='Fault trace' if free else 'Other fault trace'))

        # limits: the bulk of the points, the model and the GNSS
        ally = np.concatenate([y for _, _, y in pts])
        lo_y, hi_y = np.percentile(ally, [0.5, 99.5])
        lo_y, hi_y = min(lo_y, band_lo.min()), max(hi_y, band_hi.max())
        if has_gnss:
            g = gin['value'] + self.gnss_shift()
            lo_y = min(lo_y, np.nanmin(g - gin['sigma']))
            hi_y = max(hi_y, np.nanmax(g + gin['sigma']))
        pad = 0.06 * (hi_y - lo_y or 1.0)
        ax.set_ylim(v.model_ylim if v.model_ylim else (lo_y - pad, hi_y + 0.25 * (hi_y - lo_y)
                                                        + pad))
        span = hi_x - lo_x
        ax.set_xlim(lo_x - 0.01 * span, hi_x + 0.01 * span)

        # the A / A' ends, inside the top corners
        kw = dict(xycoords='axes fraction', fontsize=v.font_size + 3, fontweight='bold',
                  va='top', zorder=12, bbox=dict(boxstyle='round,pad=0.2', fc='w', ec='0.4',
                                                 lw=0.6, alpha=0.85))
        ends = [ax.annotate(self.label, xy=(0.012, 0.975), ha='left', **kw),
                ax.annotate(f"{self.label}'", xy=(0.988, 0.975), ha='right', **kw)]

        # the parameters, as Hussain et al. (2016): "slip rate = 20 mm/yr (11 – 29, 95% CI)"
        lines = self.parameter_lines()

        # ... and the legend: each where it hides the least (points, bins,
        # the model, the GNSS), never on the A / A' labels nor on each other
        occupied, weights = [np.c_[x, y] for _, x, y in pts], [np.ones(x.size) for _, x, _ in pts]
        occupied += bins_xy + [np.c_[xx[::4], best[::4]]]
        weights += [np.full(b.shape[0], 25.0) for b in bins_xy] + [np.full(best[::4].size, 25.0)]
        if has_gnss:
            occupied.append(np.c_[gin['x'], gin['value'] + self.gnss_shift()])
            weights.append(np.full(gin['x'].size, 200.0))
        placer = _Placer(ax, np.vstack(occupied), np.concatenate(weights), ends)

        def make_box(loc, anchor):
            box = AnchoredText('\n'.join(lines), loc=loc, bbox_to_anchor=anchor,
                               bbox_transform=ax.transAxes, borderpad=0.0, frameon=True,
                               prop=dict(fontsize=max(v.font_size - 4, 6)))
            box.patch.set(boxstyle='round,pad=0.35', fc='w', ec='0.5', alpha=0.88)
            box.set_zorder(10)
            ax.add_artist(box)
            return box

        def make_legend(loc, anchor):
            return ax.legend(handles=handles, loc=loc, bbox_to_anchor=anchor,
                             bbox_transform=ax.transAxes, fontsize=max(v.font_size - 5, 6),
                             framealpha=0.88, handlelength=1.5, borderpad=0.35,
                             labelspacing=0.3, borderaxespad=0.0)

        room = not v.model_ylim
        placer.place(make_box, room=room)
        placer.place(make_legend, room=room, max_grow=0.25)

        n_files = len(res['sets'])
        title = (f'{self.name}: swath ±{v.width:g} km, {total:,} pixels'
                 + (f' in {n_files} files' if n_files > 1 else '')
                 + (f', {res["gnss"]["x"].size} GNSS' if res['gnss'] is not None else '')
                 + f'\nstrike {self.strike:.0f}°, {COMP_TEXT[v.comp]}')
        ax.set_title(title, fontsize=v.font_size - 1, pad=6)
        return

    def parameter_lines(self):
        """The parameters as Hussain et al. (2016) write them: the MAP value, its
        unit, and the 95 % range - the slip rate as a rate, with its sense."""
        v, res = self.v, self.result
        summ = res['summary']
        lines = []
        for name in [n for n in PHYSICAL if n in summ]:
            s = summ[name]
            unit = v.model_unit if UNITS[name] == 'mm/yr' else UNITS[name]
            val, lo, hi, extra = s['map'], s['lo'], s['hi'], ''
            if name == 'S' and self.sense_known:
                sense = slip_sense(s['map'], self.strike, self.prof_az)
                if lo * hi > 0 and val * lo > 0:      # one sense all over the range
                    val, (lo, hi) = abs(val), sorted((abs(lo), abs(hi)))
                    extra = f', {sense}'
                else:
                    extra = f', {sense} (the range crosses 0: S > 0 is right-lateral)'
            # LOS data fitted as they are: the rates are what the LOS sees of them
            seen = ' in the LOS' if v.data_los and name in ('S', 'C') else ''
            lines.append(f'{NAMES[name]}{seen} = {_num(val)}{_unit_text(unit)} '
                         f'({_num(lo)} – {_num(hi)}, 95% CI){extra}')
        d = res['derived']['strain']
        lines.append(f'strain rate{" in the LOS" if v.data_los else ""} = {d["map"]:.2f} '
                     f'{strain_unit(v.model_unit)} ({d["lo"]:.2f} – {d["hi"]:.2f}, 95% CI)')
        return lines

    def draw_strain(self):
        ax, res = self.ax_strain, self.result
        lo, hi = self.ax_swath.get_xlim()
        xx = np.linspace(lo, hi, 801)
        best, band_lo, band_hi = strain_rate(res, xx)
        ax.fill_between(xx, band_lo, band_hi, color=BAND_COLOR, alpha=0.55, lw=0)
        ax.plot(xx, best, ls='--', color=MAP_COLOR, lw=1.8)
        ax.axhline(0, color='k', lw=0.4)
        x0 = res['summary']['x0']['map'] if 'x0' in res['summary'] else 0.0
        ax.axvline(x0, color='k', ls=':', lw=1.1)
        return

    def draw_topo(self):
        v, ax = self.v, self.ax_topo
        tden = self.topo_den
        if tden is None:
            ax.text(0.5, 0.5, 'no topography for this profile' + (
                f'\n{textwrap.shorten(self.topo_error, 90)}' if self.topo_error else ''),
                ha='center', va='center', color='0.45', fontsize=max(v.font_size - 3, 6),
                transform=ax.transAxes)
            return
        x = (tden['x_centers'] - self.x_fault) * self.along_scale
        ax.plot(x, tden['median'], color='k', lw=1.3)
        finite = np.isfinite(tden['median'])
        if finite.any():
            lo, hi = np.nanmin(tden['median'][finite]), np.nanmax(tden['median'][finite])
            pad = max((hi - lo) * 0.12, 1.0)
            ax.set_ylim(lo - pad, hi + pad)
        x0 = self.result['summary']['x0']['map'] if 'x0' in self.result['summary'] else 0.0
        ax.axvline(x0, color='k', ls=':', lw=1.1)
        return

    def draw_posterior(self):
        """The posterior, as Hussain et al. (2016, Fig. 9): the joint distribution
        of every pair of parameters as a 2D histogram (jet; each panel scaled to
        its own peak, as their bin areas differ and a strong correlation such as
        slip rate - locking depth would otherwise leave every other panel dark
        blue), the marginal histograms along the bottom row (blue, one frequency
        scale), the MAP as a red line and dot, the 95 % range dashed."""
        from matplotlib.ticker import MaxNLocator

        v, res = self.v, self.result
        names = [n for n in POSTERIOR_ORDER if n in res['names']]
        idx = [res['names'].index(n) for n in names]
        samples, best = res['samples'][:, idx], res['best'][idx]
        num = len(names)
        fsize = max(v.font_size - 4, 6)
        gap = 0.3                                      # inches: room for the MAP titles
        head = 0.22 + 2 * 1.3 * fsize / 72.0           # inches: the two caption lines
        size = (1.75 * num + 1.0, 1.6 * num + 0.65 + head + gap)
        if self.fig_post is None or not plt.fignum_exists(self.fig_post.number):
            self.fig_post = _new_figure(f'Posterior {self.name}', size)
            self.fig_post.canvas.mpl_connect('key_press_event', v.on_key)
            if v.disp_fig:               # made after plt.show(): show it ourselves
                try:
                    self.fig_post.show()
                except Exception:        # noqa: BLE001 - a backend without windows
                    pass
        else:
            self.fig_post.clf()
        fig = self.fig_post
        # rows 0 .. num-2 the joint distributions, row num-1 a spacer, row num the
        # marginals
        cell = (size[1] - 0.7 - head - gap) / num
        gsp = fig.add_gridspec(num + 1, num, height_ratios=[1] * (num - 1) + [gap / cell, 1],
                               hspace=0.12, wspace=0.12, left=1.0 / size[0],
                               right=1 - 0.15 / size[0], bottom=0.7 / size[1],
                               top=1 - head / size[1])
        lims = []
        for k in range(num):
            lo, hi = np.percentile(samples[:, k], [0.25, 99.75])
            lo, hi = min(lo, best[k]), max(hi, best[k])
            pad = 0.04 * (hi - lo) if hi > lo else max(abs(lo) * 0.05, 1e-3)
            lims.append((lo - pad, hi + pad))

        def unit(name):
            return v.model_unit if UNITS.get(name, 'mm/yr') == 'mm/yr' else UNITS[name]

        def title(name):
            text = POSTERIOR_LABEL.get(name, name)
            if name == 'a' and len(res['sets']) > 1:       # the offset of the first file
                text += f' {v.labels[self.fp[res["sets"][0]]["file_idx"]]}'
            return text

        # the joint distributions: row i, column j <= i: x = name j, y = name i + 1
        hists = {}
        for i in range(num - 1):
            for j in range(i + 1):
                hist, xe, ye = np.histogram2d(samples[:, j], samples[:, i + 1], bins=35,
                                              range=[lims[j], lims[i + 1]])
                hists[(i, j)] = (hist, xe, ye)
        mesh = None
        for (i, j), (hist, xe, ye) in hists.items():
            ax = fig.add_subplot(gsp[i, j])
            rel = hist / hist.max() if hist.max() > 0 else hist
            mesh = ax.pcolormesh(xe, ye, np.ma.masked_equal(rel.T, 0), cmap='jet', vmin=0,
                                 vmax=1, shading='flat', rasterized=True)
            ax.plot(best[j], best[i + 1], 'o', mfc=POST_COLOR, mec='k', ms=5.5, mew=0.9,
                    zorder=4)
            ax.set_xlim(lims[j])
            ax.set_ylim(lims[i + 1])
            ax.tick_params(direction='in', labelsize=fsize, top=True, right=True,
                           labelbottom=False, labelleft=j == 0)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
            if j == 0:
                ax.set_ylabel(f'{title(names[i + 1])}\n[{unit(names[i + 1])}]',
                              fontsize=fsize + 1)

        # the marginals, along the bottom row, on one frequency scale
        for j in range(num):
            ax = fig.add_subplot(gsp[num, j])
            counts = ax.hist(samples[:, j], bins=40, range=lims[j], color=HIST_COLOR,
                             ec=HIST_EDGE, lw=0.3)[0]
            if j == 0:
                first = ax
                top = 0.0
            top = max(top, float(counts.max()))
            s = res['summary'][names[j]]
            ax.axvline(s['map'], color=POST_COLOR, lw=1.6, zorder=4)
            for edge in (s['lo'], s['hi']):
                ax.axvline(edge, color='0.25', ls='--', lw=0.7, zorder=3)
            ax.set_xlim(lims[j])
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.set_xlabel(f'{title(names[j])} [{unit(names[j])}]', fontsize=fsize + 1)
            ax.tick_params(direction='in', labelsize=fsize, top=True, right=True,
                           labelleft=j == 0)
            if j == 0:
                ax.set_ylabel('Frequency', fontsize=fsize + 1)
            else:
                ax.sharey(first)
            ax.set_title(f'{_fmt(s["map"])} [{_fmt(s["lo"])}, {_fmt(s["hi"])}]',
                         fontsize=fsize, pad=3)
        first.set_ylim(0, top * 1.06)

        # the section and the colour scale, in the empty corner (top right)
        corner = gsp[0, num - 1].get_position(fig)
        fig.text(corner.x1, corner.y1, f" {self.label} - {self.label}' ", ha='right',
                 va='top', fontsize=fsize + 7, fontweight='bold', color='w',
                 bbox=dict(boxstyle='square,pad=0.25', fc='k', ec='k'))
        if mesh is not None and num > 2:
            width = min(1.5, num - 1.2) * corner.width
            cax = fig.add_axes([corner.x1 - width, corner.y0 + 0.12 * corner.height,
                                width, 0.16 / size[1]])
            cbar = fig.colorbar(mesh, cax=cax, orientation='horizontal', ticks=[0, 0.5, 1])
            cbar.set_label('frequency', fontsize=fsize + 1)
            cax.xaxis.set_label_position('top')
            cbar.ax.tick_params(labelsize=fsize)
        fig.text(1.0 / size[0], 1 - 0.08 / size[1],
                 f'{res["samples"].shape[0]:,} samples ({res["walkers"]} walkers x '
                 f'{res["steps"]} steps)\nMAP red, 95 % range dashed', fontsize=fsize,
                 ha='left', va='top', color='0.3', linespacing=1.3)
        fig.canvas.draw_idle()
        return

    def mark_fault_on_map(self):
        """A small black diamond where the modelled fault meets the profile."""
        v = self.v
        for art in self.map_marks:
            try:
                art.remove()
            except (ValueError, NotImplementedError):
                pass
        self.map_marks = []
        geom = line_geometry(v.atr, self.start_yx, self.end_yx, v.width)
        t = self.x_fault / (geom['length'] * geom['dist_scale'])
        row = self.start_yx[0] + t * (self.end_yx[0] - self.start_yx[0])
        col = self.start_yx[1] + t * (self.end_yx[1] - self.start_yx[1])
        x, y = v._yx2disp([row], [col])
        self.map_marks.append(v.ax_img.plot(x, y, 'D', ms=7, mfc='w', mec='k', mew=1.2,
                                            zorder=12)[0])
        v.fig_map.canvas.draw_idle()
        return

    ##---------- output
    def save(self):
        if not self.has_data:
            print(f'section {self.name}: nothing to save yet.')
            return []
        v, res = self.v, self.result
        model, summ = res['model'], res['summary']
        base = self.outfile_base
        os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
        written = []

        # 1. the figures
        out = f'{base}{v.fig_ext}'
        self.fig.savefig(out, bbox_inches='tight', dpi=v.fig_dpi)
        written.append(out)
        if self.fig_post is not None:
            out = f'{base}_posterior{v.fig_ext}'
            self.fig_post.savefig(out, bbox_inches='tight', dpi=v.fig_dpi)
            written.append(out)

        # 2. the parameters
        unit_of = {n: (v.model_unit if UNITS.get(n, 'mm/yr') == 'mm/yr' else UNITS[n])
                   for n in res['names']}
        strain_unit = 'microstrain/yr' if v.model_unit == 'mm/yr' else f'{v.model_unit}/km'
        out = f'{base}_params.csv'
        with open(out, 'w') as fid:
            fid.write('name,description,unit,map,median,lower_95,upper_95\n')
            for name in res['names']:
                s = summ[name]
                fid.write(f'{name},{describe(name)},{unit_of[name]},{s["map"]:.6g},'
                          f'{s["median"]:.6g},{s["lo"]:.6g},{s["hi"]:.6g}\n')
            d = res['derived']['strain']
            fid.write(f'strain,shear strain rate at the fault |S|/(pi d1),{strain_unit},'
                      f'{d["map"]:.6g},{d["median"]:.6g},{d["lo"]:.6g},{d["hi"]:.6g}\n')
        written.append(out)

        out = f'{base}_params.txt'
        with open(out, 'w') as fid:
            fid.write(self.report(unit_of, strain_unit))
        written.append(out)

        # 3. the posterior samples (thinned)
        out = f'{base}_samples.csv'
        samples = res['samples']
        if samples.shape[0] > 20000:
            samples = samples[np.random.default_rng(2).choice(samples.shape[0], 20000,
                                                              replace=False)]
        i_s, i_d = model.index('S'), model.index('d1')
        strain = np.abs(samples[:, i_s]) / (np.pi * samples[:, i_d])
        np.savetxt(out, np.column_stack([samples, strain]), fmt='%.6g', delimiter=',',
                   header=','.join(res['names'] + ['strain']), comments='')
        written.append(out)

        # 4. the bins, with the model
        out = f'{base}_bins.csv'
        rows = []
        for j, k in enumerate(res['sets']):
            b = res['bins'][j]
            pred = model.predict(res['best'], b['x'], offset=model.offsets[j])
            shift = self.offset_shift(j)
            rows.append(np.column_stack([np.full(b['x'].size, self.fp[k]['file_idx']), b['x'],
                                         b['median'], b['std'], b['q25'], b['q75'],
                                         b['count'], pred, b['median'] - pred,
                                         b['median'] + shift]))
        np.savetxt(out, np.vstack(rows), fmt=['%d', '%.4f', '%.6g', '%.6g', '%.6g', '%.6g',
                                              '%d', '%.6g', '%.6g', '%.6g'], delimiter=',',
                   header='file_idx,x_km,median,std,q25,q75,count,model_map,residual,'
                          'median_shifted', comments='')
        written.append(out)

        # 5. every pixel
        out = f'{base}_points.csv'
        blocks = []
        for j, k in enumerate(res['sets']):
            f = self.fp[k]
            blocks.append(np.column_stack([
                f['dist'], f['x'], f['perp'], f['value'], f['v'], f['v'] + self.offset_shift(j),
                f['row'], f['col'], f['lon'], f['lat'], np.full(f['x'].size, f['file_idx'])]))
        rc_fmt = '%.4f' if v.nlook > 1 else '%d'
        np.savetxt(out, np.vstack(blocks), delimiter=',', comments='',
                   fmt=['%.5f', '%.5f', '%.5f', '%.6g', '%.6g', '%.6g', rc_fmt, rc_fmt,
                        '%.7f', '%.7f', '%d'],
                   header='dist_km,x_km,perp_km,value,v_model,v_shifted,row,col,lon,lat,'
                          'file_idx')
        written.append(out)

        # 6. the GNSS
        gin = self.gnss_in
        if gin is not None and gin['index'].size:
            out = f'{base}_gnss.csv'
            shift = self.gnss_shift()
            lon, lat = v.gnss.lon[gin['index']], v.gnss.lat[gin['index']]
            with open(out, 'w') as fid:
                fid.write('site,lon,lat,dist_km,x_km,perp_km,value,sigma,value_shifted,'
                          'fitted\n')
                for k in range(gin['index'].size):
                    fid.write(f"{gin['site'][k]},{lon[k]:.6f},{lat[k]:.6f},"
                              f"{gin['dist'][k]:.4f},{gin['x'][k]:.4f},{gin['perp'][k]:.4f},"
                              f"{gin['value'][k]:.6g},{gin['sigma'][k]:.6g},"
                              f"{gin['value'][k] + shift:.6g},{int(model.gnss)}\n")
            written.append(out)
        for path in written:
            print(f'saved: {path}')
        return written

    def report(self, unit_of, strain_unit):
        """The readable _params.txt."""
        v, res = self.v, self.result
        model, summ = res['model'], res['summary']
        lines = [f'viewinsar model - section {self.name}', '']
        lines.append('parameter   MAP          95% range                 median       '
                     'description')
        for name in res['names']:
            s = summ[name]
            rng = f'[{s["lo"]:.4g}, {s["hi"]:.4g}]'
            lines.append(f'{name:<10s}  {s["map"]:<11.4g}  {rng:<24s}  {s["median"]:<11.4g}  '
                         f'{describe(name)} [{unit_of[name]}]')
        d = res['derived']['strain']
        rng = f'[{d["lo"]:.4g}, {d["hi"]:.4g}]'
        lines.append(f'{"strain":<10s}  {d["map"]:<11.4g}  {rng:<24s}  {d["median"]:<11.4g}  '
                     f'shear strain rate at the fault, |S| / (pi d1) [{strain_unit}]')
        lines.append('')
        if self.sense_known:
            lines.append(f'sense of slip  : {slip_sense(summ["S"]["map"], self.strike, self.prof_az)}'
                         f' (S > 0 = right-lateral with this strike)')
        lines += [
            'model          : v(x) = S/pi atan((x - x0)/d1) + a'
            + (' + C [H(x - x0) - 1/2 - 1/pi atan((x - x0)/d2)]' if model.creep else '')
            + (f' + {model.rot:g} (x - x0) on the {"x > 0" if model.rot_side > 0 else "x < 0"}'
               ' side' if model.rot else ''),
            '                 x: km from the fault, perpendicular to it; + towards '
            f"{self.label}'",
            f'data           : {COMP_LONG[v.comp]} (--comp {v.comp}), in {v.model_unit}'
            + (f' (read in {v.disp_unit})' if v.disp_unit != v.model_unit else ''),
        ]
        for j, k in enumerate(res['sets']):
            f = self.fp[k]
            noise = res['noise'][j]
            lines.append(f'file {j + 1:<10d}: {os.path.abspath(v.file[f["file_idx"]])}')
            lines.append(f'                 {f["x"].size:,} pixels, {res["bins"][j]["x"].size} '
                         f'bins, offset {model.offsets[j]}'
                         + (f', LOS sees {f["proj"]:.2f} of the fault-parallel motion'
                            if v.comp == 'los' else ''))
            lines.append('                 noise: ' + (
                'white, the scatter of each bin' if res['noise_mode'] == 'white' else
                f'exponential covariance, sill {noise[0]:.4g} ({v.model_unit})^2, '
                f'length {noise[1]:.4g} km' + (' (given)' if res['noise_mode'] == 'given'
                                                else ' (from the residuals)')))
        dropped = [self.fp[k]['file_idx'] for k in range(len(self.fp)) if k not in res['sets']]
        if dropped:
            lines.append(f'left out       : {", ".join(v.labels[i] for i in dropped)} '
                         '(too few bins)')
        gin = self.gnss_in
        if v.gnss is not None:
            n_in = 0 if gin is None else int(gin['index'].size)
            lines.append(f'GNSS           : {os.path.abspath(v.gnss_file)}, {n_in} station(s) '
                         f'within {v.gnss_width or v.width:g} km of the line, component '
                         f'{v.gnss_comp_used}, ' + ('fitted with their own offset ag'
                                                    if model.gnss else 'shown, not fitted'))
        lines += [
            f'profile        : start y/x {list(self.start_yx)}, end y/x {list(self.end_yx)}'
            + (f', lat/lon {self.start_lalo} -> {self.end_lalo}' if self.start_lalo else ''),
            f'                 azimuth {self.prof_az:.2f} deg, swath half width {v.width:g} km',
            f'fault          : at {self.x_fault:.3f} km along the profile - {self.fault_how}',
            f'strike         : {self.strike:.2f} deg - {self.strike_how}; the profile is '
            f'{self.obliquity:.1f} deg off the fault normal',
            f'bins           : {res["bin_km"]:g} km, the median of each (3+ pixels)',
            'priors         : uniform, ' + ', '.join(
                f'{n} {model.bounds[n][0]:.4g} to {model.bounds[n][1]:.4g}'
                for n in res['names'] if n in ('S', 'd1', 'C', 'd2', 'x0'))
            + (f'; d1 ~ N({model.d1_prior[0]:g}, {model.d1_prior[1]:g})' if model.d1_prior
               else '') + ('; d2 <= d1' if model.creep else ''),
            f'sampler        : affine-invariant ensemble MCMC (Goodman & Weare 2010), '
            f'{res["walkers"]} walkers x {res["steps"]} steps, burn-in {res["burn"]}, seed '
            f'{v.seed}',
            f'                 acceptance {res["acceptance"]:.2f}, autocorrelation time '
            f'{np.max(res["tau"]):.0f} steps, {res["samples"].shape[0]:,} samples kept',
            'MAP            : the maximum a posteriori model (the best sample, polished)',
            '95% range      : the 2.5 and 97.5 percentiles of the posterior samples',
        ]
        return '\n'.join(lines) + '\n'


##############################################################################
class ModelViewer(swathProfileViewer):
    """Interactive dislocation-model fits of profiles across a fault.

    Example:
        from viewinsar.cli import get_view_cmd, model_cmd_line_parse
        from viewinsar.modelview import ModelViewer

        inps = model_cmd_line_parse('geo_velocity.h5 --comp los --fault faults.shp'.split())
        obj = ModelViewer(inps, get_view_cmd(inps))
        obj.open()
        obj.plot()
    """

    COMMAND = 'model'
    SECTION_CLASS = ModelSection

    def __init__(self, inps, view_cmd):
        super().__init__(inps, view_cmd)
        self._los = {}
        self.gnss_comp_used = None
        self.data_los = False            # LOS data, fitted as they are (--comp auto)

    ##---------- units
    @property
    def to_model(self):
        """The factor from the data's unit to that of the model (mm/yr)."""
        factor = velocity_factor(self.disp_unit, 'mm/yr')
        return 1.0 if factor is None else factor

    @property
    def model_unit(self):
        return 'mm/yr' if velocity_factor(self.disp_unit, 'mm/yr') is not None else (
            str(self.disp_unit or '').strip() or 'value')

    @property
    def s_max_used(self):
        if self.s_max is not None:
            return self.s_max
        if self.model_unit == 'mm/yr':
            return 60.0
        spans = [np.nanpercentile(d, 99) - np.nanpercentile(d, 1) for d in self.data_list
                 if np.isfinite(d).any()]
        return float(max(spans + [1.0]) * 4)

    ##---------- setup
    def open(self):
        super().open()
        if 'Y_FIRST' not in self.atr and self.lookup is None:
            self.init_lookup()
        self.resolve_comp()
        if self.model_unit != 'mm/yr':
            print(f'NOTE: the data are in {self.disp_unit!r}, not a velocity: the model is '
                  f'fitted in that unit (S bound +/-{self.s_max_used:.3g}). A file without a '
                  'unit (a GeoTIFF) takes one with --data-unit, e.g. --data-unit mm/yr')
        elif self.print_msg and velocity_factor(self.disp_unit, 'mm/yr') != 1.0:
            print(f'model: in mm/yr (the data are read in {self.disp_unit})')
        if self.comp == 'as-is' and self.print_msg:
            print('NOTE: --comp as-is: the data are fitted as they are; for LOS data '
                  '--comp los gives the fault-parallel velocity (and the sense of slip)')
        return

    def resolve_comp(self):
        """--comp auto: los when the LOS geometry of every file is known and no
        file is named as a component of a decomposition; as-is otherwise."""
        from .cli import check_model_gnss

        was_auto = self.comp == 'auto'
        if was_auto:
            why, named = None, []
            for i, fname in enumerate(self.file):
                tokens = set(re.split(r'[^a-z0-9]+', os.path.basename(fname).lower()))
                named = sorted(tokens & DECOMPOSED)
                if named:
                    why = f'{os.path.basename(fname)} is named as a component ({named[0]})'
                    break
                try:
                    geom = self._los_geometry(i)
                except ValueError:
                    why = f'no LOS geometry for {os.path.basename(fname)}'
                    break
                if geom.kind == 'nominal':
                    why = (f'only nominal LOS angles for {os.path.basename(fname)}, too rough '
                           'to turn its LOS into fault-parallel velocity')
                    break
            self.comp = 'as-is' if why else 'los'
            self.data_los = bool(why) and not named
            if self.print_msg:
                print(f'--comp auto -> {self.comp}: ' + (
                    f'{why}; the data are fitted as they are' if why else
                    'the LOS velocity becomes fault-parallel velocity, pixel by pixel '
                    '(--comp as-is fits the LOS itself)'))
            if self.comp == 'as-is' and self.gnss_file and self.gnss_comp is None:
                # the GNSS component that matches the data, as they are
                self.gnss_comp = (DECOMPOSED_GNSS.get(named[0]) if named else 'los')
                if self.gnss_comp and self.print_msg:
                    print(f'GNSS: --vel {self.gnss_comp}, to match the data as they are')
            check_model_gnss(self)
        elif self.comp == 'los':
            for i, fname in enumerate(self.file):
                try:
                    nominal = self._los_geometry(i).kind == 'nominal'
                except ValueError:
                    continue                        # said where the LOS is needed
                if nominal:
                    raise ValueError(
                        f'--comp los: no LOS geometry for {os.path.basename(fname)} - only '
                        'nominal angles from its name, too rough to turn the LOS into '
                        'fault-parallel velocity. Give --geom (MintPy geometryGeo.h5, LiCSBAS '
                        'GEOCml*/, LiCSAR *.geo.E.tif), or --los-inc and --los-az; or fit the '
                        'LOS as it is with --comp as-is')
        if self.gnss is not None:
            self.gnss_comp_used = self.gnss_comp
            self.vprint(f'GNSS: the {self.gnss_comp} component, '
                        f'{"fitted" if self.gnss_fit else "shown"} [{self.model_unit}]')
        return

    def init_gnss(self):
        """Read the GNSS and place them on the map; their component is taken
        per section, once the strike is known."""
        from .gnss import read_gnss

        try:
            gnss = read_gnss(self.gnss_file, unit=self.gnss_unit, verbose=self.print_msg)
            factor = velocity_factor(gnss.unit, self.model_unit)
            if factor is None:
                print(f'WARNING: --gnss: the data are in {self.disp_unit!r}, not a velocity;'
                      f' the GNSS stay in {gnss.unit}')
                factor = 1.0
            self.gnss = gnss.scaled(factor, self.model_unit)
            self.gnss_rows, self.gnss_cols = self.lonlat_to_rowcol(gnss.lon, gnss.lat)
            self.gnss_comp_used = self.gnss_comp
        except Exception as err:      # noqa: BLE001 - the GNSS are a nicety
            print(f'WARNING: cannot use --gnss: {err}')
            print('         --> continuing without GNSS.')
            self.gnss = None
        return

    def gnss_values(self, gnss, p_vec):
        """The GNSS component that matches the data: fault-parallel (par), east,
        north, up, or the LOS of the first file."""
        comp = self.gnss_comp
        if comp == 'par':
            return gnss.component('los', (p_vec[0], p_vec[1], 0.0))
        if comp == 'los':
            geom = self._los_geometry(0)
            return gnss.component('los', geom.unit_vector_at(gnss.lon, gnss.lat))
        return gnss.component(comp)

    ##---------- geometry
    def pixel_lonlat(self, index, rows, cols):
        """lon/lat of pixels of file `index` (its own grid); None for a
        projected grid of unknown CRS."""
        rows = np.asarray(rows, dtype=np.float64)
        cols = np.asarray(cols, dtype=np.float64)
        atr_i = self.atr_list[index]
        if 'Y_FIRST' not in atr_i:
            if self.lookup is None:
                raise ValueError(NO_LOOKUP)
            return self.lookup.to_lonlat(rows, cols)
        xs = float(atr_i['X_FIRST']) + (cols + 0.5) * float(atr_i['X_STEP'])
        ys = float(atr_i['Y_FIRST']) + (rows + 0.5) * float(atr_i['Y_STEP'])
        crs = self.crs_list[index] if index < len(self.crs_list) else self.data_crs
        if crs is None:
            return None if is_projected(atr_i) else (xs, ys)
        if same_crs(crs, 4326):
            return xs, ys
        from pyproj import Transformer
        lon, lat = Transformer.from_crs(crs, 4326, always_xy=True).transform(xs, ys)
        return np.asarray(lon), np.asarray(lat)

    def make_frame(self, row, col):
        """The local east / north frame (km) at a pixel of the first file."""
        lonlat = self.pixel_lonlat(0, [row], [col])
        if lonlat is None:             # a projected grid of unknown CRS: its metres
            return dict(kind='metric',
                        x0=float(self.atr['X_FIRST']) + (col + 0.5) * float(self.atr['X_STEP']),
                        y0=float(self.atr['Y_FIRST']) + (row + 0.5) * float(self.atr['Y_STEP']))
        lon, lat = float(np.ravel(lonlat[0])[0]), float(np.ravel(lonlat[1])[0])
        if not (np.isfinite(lon) and np.isfinite(lat)):
            raise ValueError('the fault position has no lon/lat (outside geometryRadar.h5?)')
        return dict(kind='lonlat', lon0=lon, lat0=lat)

    def en_km(self, frame, index, rows, cols):
        """(east, north, lon, lat) of pixels of file `index`, in km from the frame."""
        rows = np.asarray(rows, dtype=np.float64)
        cols = np.asarray(cols, dtype=np.float64)
        nan = np.full(rows.shape, np.nan)
        if frame['kind'] == 'metric':
            atr_i = self.atr_list[index]
            xs = float(atr_i['X_FIRST']) + (cols + 0.5) * float(atr_i['X_STEP'])
            ys = float(atr_i['Y_FIRST']) + (rows + 0.5) * float(atr_i['Y_STEP'])
            return (xs - frame['x0']) / 1000.0, (ys - frame['y0']) / 1000.0, nan, nan
        lonlat = self.pixel_lonlat(index, rows, cols)
        if lonlat is None:
            raise ValueError(f'{self.labels[index]}: a projected grid of unknown CRS')
        east, north = local_en(lonlat[0], lonlat[1], frame['lon0'], frame['lat0'])
        return east, north, np.asarray(lonlat[0]), np.asarray(lonlat[1])

    def trace_en(self, frame, xy):
        """A fault trace (in the map coordinates of the data) in km from the frame."""
        xy = np.asarray(xy, dtype=np.float64)
        if self.fault_atr is PIXEL_ATR:                   # radar: (col, row)
            lon, lat = self.lookup.to_lonlat(xy[:, 1], xy[:, 0])
        elif frame['kind'] == 'metric':
            return (xy[:, 0] - frame['x0']) / 1000.0, (xy[:, 1] - frame['y0']) / 1000.0
        else:
            lon, lat = xy[:, 0], xy[:, 1]
            if self.data_crs is not None and not same_crs(self.data_crs, 4326):
                from pyproj import Transformer
                lon, lat = Transformer.from_crs(self.data_crs, 4326,
                                                always_xy=True).transform(lon, lat)
        return local_en(lon, lat, frame['lon0'], frame['lat0'])

    def trace_strike(self, frame, radius_km=5.0, max_dist_km=1.0):
        """The strike of the fault trace passing the frame's origin - its
        vertices within radius_km, or the segment through it."""
        best = None
        for line in self.faults:
            if len(line) < 2:
                continue
            east, north = self.trace_en(frame, line)
            pts = np.c_[east, north]
            seg_a, seg_d = pts[:-1], np.diff(pts, axis=0)
            dd = np.einsum('ij,ij->i', seg_d, seg_d)
            with np.errstate(divide='ignore', invalid='ignore'):
                t = np.clip(-np.einsum('ij,ij->i', seg_a, seg_d) / dd, 0, 1)
            dist = np.hypot(*(seg_a + t[:, None] * seg_d).T)
            dist[~np.isfinite(dist)] = np.inf
            k = int(np.argmin(dist))
            if best is None or dist[k] < best[0]:
                best = (dist[k], pts, k)
        if best is None or best[0] > max_dist_km:
            return None
        _, pts, k = best
        near = np.all(np.isfinite(pts), axis=1) & (np.hypot(*pts.T) <= radius_km)
        sel = pts[near] if near.sum() >= 2 else pts[k:k + 2]
        return strike_of_points(sel[:, 0], sel[:, 1])

    def _los_geometry(self, index):
        from .gnss import LOSGeometry
        if index not in self._los:
            files = self.geom_file or []
            given = files[index] if index < len(files) else None
            self._los[index] = LOSGeometry.for_file(self.file[index], self.atr_list[index],
                                                    geom_file=given,
                                                    inc=per_file(self.los_inc, index),
                                                    az=per_file(self.los_az, index),
                                                    verbose=self.print_msg)
        return self._los[index]

    def los_vectors(self, index, rows, cols, lon=None, lat=None):
        """(east, north, up) of the LOS unit vector at pixels of file `index`."""
        from .gnss import los_unit_vector
        geom = self._los_geometry(index)
        atr_i = self.atr_list[index]
        same = geom.constant or (
            geom.atr is not None and geom.inc.shape == (int(atr_i['LENGTH']), int(atr_i['WIDTH']))
            and (('Y_FIRST' in atr_i) == ('Y_FIRST' in geom.atr))
            and all(abs(float(geom.atr.get(k, 0)) - float(atr_i.get(k, 0)))
                    <= 1e-6 * max(1.0, abs(float(atr_i.get(k, 0))))
                    for k in ('X_FIRST', 'Y_FIRST', 'X_STEP', 'Y_STEP') if k in atr_i))
        if same:
            inc, az = geom.at_rowcol(rows, cols)
        else:
            if lon is None or not np.isfinite(lon).any():
                raise ValueError(f'{self.labels[index]}: the LOS geometry is on another grid '
                                 'and the pixels have no lon/lat')
            inc, az = geom.at_lonlat(lon, lat)
        return los_unit_vector(inc, az)

    ##---------- windows
    def help_msg(self):
        msg = '\n' + '-' * 74
        msg += '\n  left click x2 : draw a profile ACROSS the fault -> the model is fitted'
        msg += '\n  key "="  or  "+" : start ANOTHER section in a new window'
        msg += '\n  key "n"       : redraw the CURRENT section from scratch'
        msg += '\n  key "e"       : export the map, every fit and posterior, and the tables,'
        msg += f'\n                  into {export_dir(self.outdir)}'
        msg += '\n' + '-' * 74
        return msg

    def close(self):
        super().close()
        for sec in self.sections:
            if sec.fig_post is not None:
                plt.close(sec.fig_post)
        return

    ##---------- output
    def save_outputs(self):
        """Every section's files, the map, and - for two sections or more - one
        table of the parameters of each, A, B, C ... (slip rate along strike)."""
        written = super().save_outputs()
        fitted = [s for s in self.sections if s.has_data]
        if len(fitted) < 2:
            return written
        out = os.path.join(self.outdir,
                           f'{self.COMMAND}_sections_{"".join(s.label for s in fitted)}.csv')
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        names = ('S', 'd1', 'C', 'd2', 'x0')
        head = ['section', 'start_lat', 'start_lon', 'end_lat', 'end_lon', 'fault_km',
                'fault_lat', 'fault_lon', 'strike_deg', 'sense', 'pixels', 'rate_unit']
        for name in names + ('strain',):
            head += [f'{name}_map', f'{name}_lower_95', f'{name}_upper_95']
        with open(out, 'w') as fid:
            fid.write(','.join(head) + '\n')
            for sec in fitted:
                fid.write(','.join(self._section_row(sec, names)) + '\n')
        print(f'saved: {out}  (one row per section; S, C in {self.model_unit}, d1, d2, x0 in '
              'km, strain in ' + ('microstrain/yr)' if self.model_unit == 'mm/yr' else
                                  f'{self.model_unit}/km)'))
        written.append(out)
        return written

    def _section_row(self, sec, names):
        """One row of the _sections.csv table."""
        res, summ = sec.result, sec.result['summary']

        def lonlat(row, col):
            try:
                ll = self.pixel_lonlat(0, [row], [col])
            except ValueError:
                return '', ''
            if ll is None:
                return '', ''
            lon, lat = float(np.ravel(ll[0])[0]), float(np.ravel(ll[1])[0])
            if not (np.isfinite(lon) and np.isfinite(lat)):
                return '', ''
            return f'{lat:.6f}', f'{lon:.6f}'

        (y0, x0), (y1, x1) = sec.start_yx, sec.end_yx
        row = [sec.label, *lonlat(y0, x0), *lonlat(y1, x1)]
        geom = line_geometry(self.atr, sec.start_yx, sec.end_yx, self.width)
        t = sec.x_fault / (geom['length'] * geom['dist_scale'])
        row += [f'{sec.x_fault:.4f}', *lonlat(y0 + t * (y1 - y0), x0 + t * (x1 - x0)),
                f'{sec.strike:.2f}',
                slip_sense(summ['S']['map'], sec.strike, sec.prof_az) if sec.sense_known
                else '', str(res['n_pixels']), self.model_unit]
        for name in names:
            s = summ.get(name)
            row += ([f'{s["map"]:.6g}', f'{s["lo"]:.6g}', f'{s["hi"]:.6g}'] if s
                    else ['', '', ''])
        d = res['derived']['strain']
        row += [f'{d["map"]:.6g}', f'{d["lo"]:.6g}', f'{d["hi"]:.6g}']
        return row

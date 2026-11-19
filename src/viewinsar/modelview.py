"""`viewinsar model`: a swath across a fault, fitted with an elastic dislocation model.

The map and the profile line work as in `viewinsar swath`. Every section gets
a window with the pixels of its swath as points - fault-parallel velocity
against the distance from the fault - their binned medians, the GNSS, and the
maximum a posteriori (MAP) model with its 95 % and 68 % ranges and the
parameters written on the figure, after Hussain et al. (2016, 2018) and Ou et
al. (2022); and a second window with the posterior distributions of the
parameters. With --cmap-prof the pixels are a pixel density instead, as in
swath, under red ranges and bigger medians. A reverse or normal fault
(`viewinsar model reverse|normal`) gets a panel per component of its data.

The fit always takes the pixels as they are: --nlook oversamples what is
drawn only - interpolated samples are no new data, and counted as such they
would narrow the 95 % ranges.

    left click x2  the profile, across the fault     "="  another section
    "n"            redo the current section          "e"  export everything
"""

import os
import re
import textwrap

import numpy as np
from matplotlib import colors as mcolors
from matplotlib import patheffects
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText

from .lineview import LineTools
from .model import (
    NAMES,
    UNITS,
    describe,
    fault_offset,
    fit_profile,
    forward,
    model_bands,
    split_offsets,
    strain_rate,
)
from .radar import PIXEL_ATR
from .raster import same_crs
from .swath import (
    EARTH_RADIUS,
    line_geometry,
    nice_step,
    project_points,
    swath_density,
    swath_extract,
)
from .viewer import (
    GNSS_EDGE,
    Section,
    _new_figure,
    export_dir,
    swathProfileViewer,
    velocity_factor,
    white_halo,
)

__all__ = ['ModelViewer', 'ModelSection', 'azimuth_of', 'local_en', 'strike_of_points',
           'orient_strike', 'slip_sense', 'fault_normal']

POINT_COLORS = ('#d62728', '#9467bd', '#2ca02c', '#ff7f0e', '#8c564b', '#e377c2',
                '#bcbd22', '#17becf')            # the pixels of each file, red first
GNSS_COLOR = '#1f77b4'                           # blue, as Hussain et al.
GNSS_ON_DENSITY = '#d62728'                      # --cmap-prof: red, outlined in white
MAP_COLOR = 'k'                                  # the MAP model: bold black dashes
BAND_COLOR = '0.62'                              # the 95 % range: grey ...
BAND_68 = '0.86'                                 # ... the 68 % range: light grey
BAND_ALPHA = (0.75, 0.6)                         # ... 25 % and 40 % transparent
# the y axis says what is plotted: fault-parallel velocity once the data are
# turned into it; else what they are, as they are (a component by its name);
# for reverse and normal faults: across the fault (h), vertical (z), a LOS
VEL_LABEL = {'par': 'Fault-parallel velocity', 'vu': 'Fault vertical velocity',
             'los': 'LOS velocity', 'east': 'East velocity', 'north': 'North velocity',
             'horizontal': 'Horizontal velocity', 'h': 'Fault-normal velocity',
             'z': 'Fault vertical velocity'}
# the far-field rates of a dip-slip fault: the hanging wall against the footwall
RATE_WORDS = {'reverse': ('shortening rate', 'uplift rate'),
              'normal': ('extension rate', 'subsidence rate')}
DIP_COMPS = ('vu', 've', 'vn', 'los')            # the data of a reverse or normal fault
NAMED_AS = {'up': 'vu', 'vert': 'vu', 'vertical': 'vu', 'vu': 'vu', 'east': 'east',
            'ew': 'east', 've': 'east', 'north': 'north', 'ns': 'north', 'vn': 'north',
            'hz': 'horizontal', 'horz': 'horizontal', 'horizontal': 'horizontal'}
LEGEND_CYCLE = ('in', 'out', 'off')              # the "l" key, in turn
# --cmap-prof: the pixels as a density; the 95 % range red, as the points were
DENSITY_BAND = '#d62728'
OFFSET_COLOR = '#ff7f0e'                         # --creep-offset: the two lines, orange
POST_COLOR = '#d62728'                           # the MAP in the posterior plots
HIST_COLOR, HIST_EDGE = '#1f4fd8', '#0a1f73'      # the marginal histograms, blue
# names that mark a component of a decomposition, not a line of sight
DECOMPOSED = {'hz', 'horz', 'horizontal', 'up', 'vert', 'vertical', 'east', 'north', 'ew',
              'ns', 've', 'vn', 'vu', 'decomp', 'decomposed'}
# the GNSS component that matches a file named as one (horizontal: unknown direction)
DECOMPOSED_GNSS = {'up': 'vu', 'vert': 'vu', 'vertical': 'vu', 'vu': 'vu', 'east': 've',
                   'ew': 've', 've': 've', 'north': 'vn', 'ns': 'vn', 'vn': 'vn'}
POSTERIOR_ORDER = ('S', 'd1', 'dip', 'a', 'C', 'd2', 'x0')   # as Hussain et al. (2016)
POSTERIOR_LABEL = {'S': 'Slip rate', 'd1': 'Locking depth', 'a': 'Offset', 'dip': 'Dip',
                   'C': 'Creep rate', 'd2': 'Creep depth', 'x0': 'Fault position'}
KM_PER_DEG = np.pi / 180.0 * EARTH_RADIUS / 1000.0
PHYSICAL = ('S', 'd1', 'dip', 'C', 'd2', 'x0')   # what the posterior figure shows
SYMBOL = {'S': '$S$', 'd1': '$d_1$', 'C': '$C$', 'd2': '$d_2$', 'x0': '$x_0$', 'dip': r'$\delta$',
          'a': '$a$', 'ag': '$a_{GNSS}$'}
COMP_TEXT = {'as-is': 'the data as given', 'par': 'fault-parallel data',
             've': 'from the east velocity', 'vn': 'from the north velocity',
             'los': 'from the LOS velocity', 'vu': 'the vertical velocity',
             'mixed': 'several components'}
DIP_COMP_TEXT = {'vu': 'the vertical velocity', 've': 'the east velocity, across the fault',
                 'vn': 'the north velocity, across the fault',
                 'los': 'the LOS velocity, as it is'}
COMP_LONG = {'as-is': 'the data as given, no conversion',
             'par': 'fault-parallel velocity, as given',
             've': 'east velocity / sin(strike)', 'vn': 'north velocity / cos(strike)',
             'los': 'LOS velocity / (LOS unit vector . strike direction), pixel by pixel',
             'vu': 'vertical velocity, as given'}


def profile_cmap(name):
    """The colormap of --cmap-prof: any MintPy or matplotlib name, its top end
    white beyond the colour scale - as swath draws its density."""
    from .tsviewer import get_cmap

    cmap = get_cmap(name)
    try:
        return cmap.with_extremes(over='white')          # matplotlib >= 3.4
    except AttributeError:
        cmap = cmap.copy()
        cmap.set_over('white')
        return cmap


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


def ci_text(value, lo, hi, unit=''):
    """'20.1 mm/yr (19.4 – 22.2, 95% CI)' - a value and its 95 % range; an
    angle in degrees: '29.9° (27.4 – 33.0, 95% CI)'."""
    unit = '°' if unit in ('deg', '°') else (_unit_text(unit) if unit else '')
    return f'{_num(value)}{unit} ({_num(lo)} – {_num(hi)}, 95% CI)'


def keep_default_keys(fig, keys):
    """matplotlib's own shortcuts on a figure (s to save, o to zoom ...), but
    for `keys`, which are ours - "l" would otherwise set a log y axis."""
    from matplotlib.backend_bases import key_press_handler

    manager = getattr(fig.canvas, 'manager', None)
    handler_id = getattr(manager, 'key_press_handler_id', None)
    if handler_id is None:
        return

    def handler(event):
        if event.key not in keys:
            key_press_handler(event)

    fig.canvas.mpl_disconnect(handler_id)
    manager.key_press_handler_id = fig.canvas.mpl_connect('key_press_event', handler)


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
    making room above or below the data when every spot would hide some. An
    upper spot is tried at the top of the axes and under each label there."""

    SPOTS = (('upper left', 0.012), ('upper right', 0.988), ('lower left', 0.012),
             ('lower right', 0.988), ('upper center', 0.5), ('lower center', 0.5),
             ('center left', 0.012), ('center right', 0.988))

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
        tops = sorted({r[1] - 0.03 for r in self.taken if r[1] + r[3] > 1.0}, reverse=True)
        floors = sorted({r[3] + 0.02 for r in self.taken if r[1] + r[3] <= 1.0})
        self.top = tops[-1] if tops else 0.985
        # the heights a box may hang from (upper) or stand on (lower)
        self.anchors = {'upper': [0.985] + tops, 'lower': [0.015] + floors, 'center': [0.5]}

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
        for loc, ax_x in self.SPOTS:
            vertical, horizontal = loc.split()
            for ax_y in self.anchors[vertical]:
                x0 = {'left': ax_x, 'right': ax_x - width}.get(horizontal, ax_x - width / 2)
                y0 = {'upper': ax_y - height, 'lower': ax_y}.get(vertical, ax_y - height / 2)
                cand = (x0, y0, x0 + width, y0 + height)
                cost = self._cost(cand)
                if best is None or cost < best[0]:
                    best = (cost, loc, (ax_x, ax_y), cand)
        artist.remove()
        return make(best[1], best[2]), best

    budget = 0.005                            # a box may hide this much of the points

    def place(self, make, room=True, budget=None, max_grow=0.4):
        """Draw a box - make(loc, anchor) - at the best spot. With room, a box
        that would hide more than `budget` of the points makes room for itself
        above (or below) the data: at most max_grow of the y-range, so the
        data keep most of the panel."""
        budget = self.budget if budget is None else budget
        artist, best = self._best(make)
        self.last = best
        if best is None:
            self.taken.append((0.0, self.top - 0.3, 0.5, self.top))
            return artist
        if room and not np.isfinite(best[0]):
            # no spot clear of the boxes already drawn: room below the data - kept
            # only when it frees a spot
            height = best[3][3] - best[3][1] + 0.03
            lo, hi = self.ax.get_ylim()
            grow = min((hi - lo) * height / max(1.0 - height, 0.4) * 1.05,
                       max(max_grow, 0.5) * (hi - lo))
            self.ax.set_ylim(lo - grow, hi)
            artist.remove()
            self.refresh()
            artist, best = self._best(make)
            if not np.isfinite(best[0]):
                artist.remove()
                self.ax.set_ylim(lo, hi)
                self.refresh()
                artist, best = self._best(make)
        if room and np.isfinite(best[0]) and best[0] > budget * self.w.sum():
            # the best spot hides points: room above the data (an upper spot) or
            # below it (a lower or a centre one) - kept only when it hides fewer
            vertical = best[1].split()[0]
            up = vertical == 'upper'
            cand, anchor = best[3], best[2][1]
            height = cand[3] - cand[1] + ((1 - anchor) if up else
                                          anchor if vertical == 'lower' else 0.015)
            lo, hi = self.ax.get_ylim()
            grow = min((hi - lo) * height / max(1.0 - height, 0.4) * 1.05, max_grow * (hi - lo))
            self.ax.set_ylim((lo, hi + grow) if up else (lo - grow, hi))
            artist.remove()
            self.refresh()
            grown, then = self._best(make)
            if then[0] < best[0]:
                artist, best = grown, then
            else:
                grown.remove()
                self.ax.set_ylim(lo, hi)
                self.refresh()
                artist, best = self._best(make)
        self.taken.append(best[3])
        self.last = best
        return artist


##############################################################################
class ModelSection(Section):
    """One profile across a fault, its fit and its two windows."""

    def __init__(self, viewer, index):
        super().__init__(viewer, index)
        self.ax_strain = None
        self.fig_post = None
        self.map_marks = []
        self.legend_where = None      # the "l" key: in / out / off (else --legend)
        self.legend_at = None         # where the legend went: in / out / off
        self.on_line = False          # a profile of --spacing: the fault is the line
        self.line_point = self.line_dist = None      # its point, km along the line
        self.clear_model()

    def clear_model(self):
        self.result = None
        self.fp = []                  # per file: x, v (fault-parallel), ... - what is fitted
        self.info = ''                # the window's title: pixels, files, strike ...
        self.fp_draw = None           # --nlook: the oversampled samples, drawn only
        self.density = None           # --cmap-prof: the density drawn
        self.offset = None            # --creep-offset: the creep rate from the step
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
    def panel_specs(self):
        """The panels of the profile: one for a strike-slip fault; for a reverse
        or normal fault one per component - across the fault (ve, vn), vertical
        (vu), and each line of sight on its own - as dict(key, kind, files)."""
        v = self.v
        if getattr(v, 'fault_type', 'strike') == 'strike':
            return [dict(key='main', kind=None, files=list(range(len(v.file))))]
        comps = list(v.comps)
        specs = []
        horizontal = [i for i, c in enumerate(comps) if c in ('ve', 'vn')]
        vertical = [i for i, c in enumerate(comps) if c == 'vu']
        if horizontal:
            specs.append(dict(key='h', kind='h', files=horizontal))
        if vertical:
            specs.append(dict(key='z', kind='z', files=vertical))
        for i, c in enumerate(comps):
            if c == 'los':
                specs.append(dict(key=f'los{i}', kind='los', files=[i]))
        kind = 'z' if v.gnss_comp == 'vu' else 'h'
        if getattr(v, 'gnss', None) is not None and all(sp['kind'] != kind for sp in specs):
            specs.append(dict(key='gnss', kind=kind, files=[]))   # the GNSS, on their own
        return specs

    def create_figure(self):
        v = self.v
        if self.fig is not None and plt.fignum_exists(self.fig.number):
            return self.fig
        self.specs = self.panel_specs()
        n_data = len(self.specs)
        panels = ([3.2] * n_data + ([1.2] if v.show_strain else [])
                  + ([1.1] if v.topoxsect else []))
        height = v.fig_size[1] + 3.0 * (n_data - 1) + 1.5 * (len(panels) - n_data)
        self.fig = _new_figure(f'Model {self.name}', [v.fig_size[0], height])
        margin = 0.1 if n_data == 1 else 0.55 / height     # stacked panels: in inches
        gsp = self.fig.add_gridspec(len(panels), 1, height_ratios=panels, hspace=0.08,
                                    left=0.12, right=0.86 if v.cmap_prof else 0.97,
                                    bottom=margin, top=1.0 - margin)
        self.ax_panels = [self.fig.add_subplot(gsp[0])]
        self.ax_swath = self.ax_panels[0]
        for k in range(1, n_data):
            self.ax_panels.append(self.fig.add_subplot(gsp[k], sharex=self.ax_swath))
        self.ax_cbars = []
        if v.cmap_prof:               # a fixed slot for each density's colour bar
            for ax in self.ax_panels:
                pos = ax.get_position()
                cax = self.fig.add_axes([pos.x1 + 0.015, pos.y0, 0.016, pos.height])
                cax.set_axis_off()
                self.ax_cbars.append(cax)
        self.ax_cbar = self.ax_cbars[0] if self.ax_cbars else None
        k = n_data
        self.ax_strain = self.ax_topo = None
        if v.show_strain:
            self.ax_strain = self.fig.add_subplot(gsp[k], sharex=self.ax_swath)
            k += 1
        if v.topoxsect:
            self.ax_topo = self.fig.add_subplot(gsp[k], sharex=self.ax_swath)
        # where the panels sit: a legend beside the plot narrows them, each draw anew
        self.home = [(ax, ax.get_position().frozen()) for ax in
                     (*self.ax_panels, self.ax_strain, self.ax_topo, *self.ax_cbars)
                     if ax is not None]
        self.fig.canvas.mpl_connect('key_press_event', v.on_key)
        keep_default_keys(self.fig, ('l', 'L'))
        self.reset_axes()
        self.ax_swath.text(0.5, 0.5, f'section {self.name}\nclick 2 points on the map,\n'
                           'across the fault', ha='center', va='center', color='0.55',
                           fontsize=v.font_size, transform=self.ax_swath.transAxes)
        return self.fig

    def reset_axes(self):
        v = self.v
        panels = getattr(self, 'ax_panels', None) or [self.ax_swath]
        axes = [ax for ax in (*panels, self.ax_strain, self.ax_topo) if ax is not None]
        if getattr(self, '_narrow', None) is not None:     # widened for a legend beside
            self.fig.set_size_inches(*self._narrow, forward=True)
        self._narrow = None
        for ax, box in getattr(self, 'home', []):
            ax.set_position(box)
        for ax in axes:
            ax.cla()
            ax.tick_params(which='both', direction='in', labelsize=v.font_size,
                           bottom=True, top=True, left=True, right=True)
        for cax in getattr(self, 'ax_cbars', []):
            cax.cla()
            cax.set_axis_off()
        for ax in axes[:-1]:
            ax.tick_params(labelbottom=False)
        axes[-1].set_xlabel('Distance from the fault (km)', fontsize=v.font_size + 1)
        unit = _unit_text(v.model_unit).strip()
        for ax, spec in zip(panels, getattr(self, 'specs', None) or self.panel_specs()):
            if spec['kind'] is None:
                what = VEL_LABEL.get(getattr(v, 'plotted', None), 'Velocity')
            elif spec['kind'] == 'los':
                what = f"{VEL_LABEL['los']}, {v.labels[spec['files'][0]]}"
            else:
                what = VEL_LABEL[spec['kind']]
            gap = '\n' if len(panels) > 1 else ' '      # stacked panels: two lines
            ax.set_ylabel(what + (f'{gap}({unit})' if unit else ''), fontsize=v.font_size + 1)
        if self.ax_strain is not None:
            self.ax_strain.set_ylabel(f'dv/dx\n({strain_unit(v.model_unit)})',
                                      fontsize=v.font_size)
        if self.ax_topo is not None:
            self.ax_topo.set_ylabel('Elevation\n(m)', fontsize=v.font_size)
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
        # the fit takes the real pixels; --nlook oversamples what is drawn only
        if not super().extract(nlook=1):
            return False
        try:
            self.prepare()
        except ValueError as err:
            print(f'WARNING: section {self.name}: {err}')
            return False
        self.fp_draw = self.oversampled() if self.v.nlook > 1 else None
        return True

    def oversampled(self):
        """--nlook: the swath of every file sampled NUM x NUM times as densely,
        as fault-parallel velocity - for the points or the density drawn."""
        v = self.v
        dense = []
        for sw in self.swath_list:
            i = sw['file_idx']
            try:
                one = swath_extract(v.data_list[i], v.atr_list[i], v.box_list[i],
                                    sw['start_yx'], sw['end_yx'], v.width, mask=v.msk_list[i],
                                    keep_zero=v.keep_zero, nlook=v.nlook)
            except ValueError:
                continue
            one['file_idx'] = i
            dense.append(one)
        return {f['file_idx']: f for f in self.fault_parallel(dense, warn=False)}

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
        if self.on_line:              # a profile of --spacing: where it crosses the line
            self.x_fault, self.fault_how = length / 2, 'the line along the fault'
        elif v.x0_km is not None:
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
        elif self.on_line:            # at right angles to the profile, by construction
            strike, self.strike_how = (self.prof_az + 90.0) % 360.0, 'the line along the fault'
        elif v.faults:
            strike = v.trace_strike(self.frame)
            if strike is not None:
                self.strike_how = 'the fault trace at the fault'
        if strike is None:
            strike = (self.prof_az + 90.0) % 360.0
            self.strike_how = 'perpendicular to the profile'
        given_sign = v.comp in ('as-is', 'par', 'vu') and v.strike is not None
        self.strike = strike if given_sign else orient_strike(strike, self.prof_az)
        self.sense_known = (v.comp in ('los', 've', 'vn') or (v.comp == 'par'
                                                              and v.strike is not None))
        if self.dip_slip:             # the sense is the fault type's
            self.sense_known = False
        self.normal, self.obliquity = fault_normal(self.strike, self.prof_az)
        if self.dip_slip:             # the right-hand rule: the fault dips to its right
            dip_az = azimuth_of(*self.normal) + (0.0 if v.dip_toward == 'end' else 180.0)
            self.strike = (dip_az - 90.0) % 360.0
        prof = np.array([np.sin(np.radians(self.prof_az)), np.cos(np.radians(self.prof_az))])
        self.along_scale = float(self.normal @ prof)
        if self.obliquity > 30:
            print(f'NOTE: section {self.name} runs {self.obliquity:.0f} deg off the fault '
                  'normal: distances are taken perpendicular to the fault (a straight line '
                  'through the fault position, along the strike)')

        # 4. every pixel: distance from the fault, fault-parallel velocity - or,
        #    across a reverse or normal fault, its component of the motion
        self.fp = (self.dip_slip_data(self.swath_list) if self.dip_slip
                   else self.fault_parallel(self.swath_list))
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
            p_vec = (tuple(self.normal) if v.gnss_comp == 'perp' else
                     (np.sin(np.radians(self.strike)), np.cos(np.radians(self.strike))))
            value, sigma = v.gnss_values(v.gnss.subset(idx), p_vec)
            gin['value'], gin['sigma'] = value, sigma
        return

    @property
    def dip_slip(self):
        """A reverse or a normal fault (`viewinsar model reverse|normal`)."""
        return getattr(self.v, 'fault_type', 'strike') != 'strike'

    def dip_slip_data(self, swaths, warn=True):
        """Every pixel across a reverse or normal fault as (distance from the
        fault, velocity), one dict per file, each with its kind: vertical (vu:
        z); across the fault, + towards the end of the profile (ve, vn: the
        east / north velocity divided by the fault normal's east / north,
        assuming no motion along the fault: h); or the LOS as it is, the model
        seen through it (ph of h, pu of z: the medians over the swath)."""
        v = self.v
        n_e, n_n = (float(c) for c in self.normal)
        out = []
        for sw in swaths:
            i = sw['file_idx']
            comp = v.comps[i]
            east, north, lon, lat = v.en_km(self.frame, i, sw['row'], sw['col'])
            x = east * n_e + north * n_n
            value = np.asarray(sw['value'], dtype=np.float64)
            kind, proj, med_proj = 'z', np.ones(value.size), 1.0
            if comp in ('ve', 'vn'):
                kind = 'h'
                proj = np.full(value.size, n_e if comp == 've' else n_n)
                med_proj = float(proj[0]) if proj.size else np.nan
                if warn and abs(med_proj) < 0.25:
                    print(f'WARNING: {v.labels[i]}: the {"east" if comp == "ve" else "north"} '
                          f'velocity sees only {abs(med_proj):.2f} of the motion across the '
                          f'fault (strike {self.strike:.0f} deg): its noise is amplified '
                          f'{1 / max(abs(med_proj), 1e-3):.1f}x - the '
                          f'{"north" if comp == "ve" else "east"} velocity sees more of it')
            elif comp == 'los':
                ue, un, up = v.los_vectors(i, sw['row'], sw['col'], lon, lat)
                ph = np.broadcast_to(ue * n_e + un * n_n, value.shape)
                pu = np.broadcast_to(up, value.shape)
                ok = np.isfinite(ph) & np.isfinite(pu)
                kind = ((float(np.median(ph[ok])), float(np.median(pu[ok]))) if ok.any()
                        else (np.nan, np.nan))
                med_proj = kind[0]
            weak = ~np.isfinite(proj) | (np.abs(proj) < 0.05)
            with np.errstate(divide='ignore', invalid='ignore'):
                vel = np.where(weak, np.nan, value * v.to_model / proj)
            good = np.isfinite(x) & np.isfinite(vel)
            one = dict(file_idx=i, x=x[good], v=vel[good], dist=sw['dist'][good],
                       perp=sw['perp'][good], value=value[good],
                       row=np.asarray(sw['row'])[good], col=np.asarray(sw['col'])[good],
                       lon=np.asarray(lon)[good], lat=np.asarray(lat)[good], proj=med_proj,
                       kind=kind, comp=comp)
            out.append(one)
        return [f for f in out if f['x'].size and np.all(np.isfinite(
            f['kind'] if isinstance(f['kind'], tuple) else 0.0))]

    def fault_parallel(self, swaths, warn=True):
        """Every pixel of some swaths as (distance from the fault, fault-parallel
        velocity), one dict per file - the files with no pixel left out."""
        v = self.v
        p_vec = (np.sin(np.radians(self.strike)), np.cos(np.radians(self.strike)))
        out = []
        for sw in swaths:
            i = sw['file_idx']
            east, north, lon, lat = v.en_km(self.frame, i, sw['row'], sw['col'])
            x = east * self.normal[0] + north * self.normal[1]
            value = np.asarray(sw['value'], dtype=np.float64)
            up = None
            if v.comp == 'los':
                ue, un, up = v.los_vectors(i, sw['row'], sw['col'], lon, lat)
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
            if warn and v.comp in ('los', 've', 'vn') and med_proj < 0.25:
                print(f'WARNING: {v.labels[i]}: the {v.comp.upper()} sees only '
                      f'{med_proj:.2f} of the fault-parallel motion (strike '
                      f'{self.strike:.0f} deg): its noise is amplified {1 / med_proj:.1f}x')
            one = dict(file_idx=i, x=x[good], v=vpar[good], dist=sw['dist'][good],
                       perp=sw['perp'][good], value=value[good],
                       row=np.asarray(sw['row'])[good], col=np.asarray(sw['col'])[good],
                       lon=np.asarray(lon)[good], lat=np.asarray(lat)[good], proj=med_proj)
            if up is not None:        # the LOS of each pixel: along the strike, and up
                one['los_par'] = np.broadcast_to(proj, value.shape)[good]
                one['los_up'] = np.broadcast_to(up, value.shape)[good]
            out.append(one)
        return [f for f in out if f['x'].size]

    ##---------- the fit
    def fit(self):
        v = self.v
        gin = self.gnss_in
        gnss = None
        if gin is not None and gin['index'].size and v.gnss_fit:
            ok = np.isfinite(gin['x']) & np.isfinite(gin['value'])
            if ok.any():
                gnss = dict(x=gin['x'][ok], v=gin['value'][ok], sigma=gin['sigma'][ok])
        free_x0 = (-v.free_x0, v.free_x0) if v.free_x0 else None
        if self.dip_slip:
            bounds = {'S': (0.0, v.s_max_used), 'd1': (0.1, v.d1_max)}
            data = ', '.join(f'{v.labels[f["file_idx"]]}: {DIP_COMP_TEXT[f["comp"]]}'
                             for f in self.fp)
            what = f'{v.fault_type} fault, dipping towards {self.dip_side()}; {data}'
        else:
            bounds = {'S': (-v.s_max_used, v.s_max_used), 'd1': (0.1, v.d1_max),
                      'C': (0.0, v.c_max), 'd2': (0.01, v.d2_max)}
            what = COMP_LONG[v.comp]
        print(f'section {self.name}: fault at {self.x_fault:.2f} km ({self.fault_how}), strike '
              f'{self.strike:.0f} deg ({self.strike_how}), profile azimuth {self.prof_az:.0f}'
              f' deg; {what}')
        self.result = fit_profile([(f['x'], f['v']) for f in self.fp], bin_km=v.bin_km,
                                  gnss=gnss, noise=v.noise, walkers=v.walkers, steps=v.steps,
                                  burn=v.burn, seed=v.seed, creep=v.creep, free_x0=free_x0,
                                  rot=v.rot, rot_side=1 if v.rot_side == 'end' else -1,
                                  bounds=bounds, d1_prior=v.d1_prior,
                                  verbose=v.print_msg, fault=getattr(v, 'fault_type', 'strike'),
                                  kinds=[f.get('kind') for f in self.fp],
                                  gnss_kind=self.gnss_kind, dip=v.dip_range,
                                  toward=1 if v.dip_toward == 'end' else -1)
        res = self.result
        if self.sense_known:
            S = res['summary']['S']['map']
            print(f'   S {S:+.2f} {v.model_unit}: {slip_sense(S, self.strike, self.prof_az)}')
        elif self.dip_slip:
            par = res['model'].unpack(res['best'])
            h, z = (res['derived'][k]['map'] for k in ('h_rate', 'v_rate'))
            words = RATE_WORDS[v.fault_type]
            print(f'   S {float(par["S"]):.2f} {v.model_unit} along the dip: {v.fault_type}, dip '
                  f'{float(par["dip"]):.0f} deg towards {self.dip_side()}; {words[0]} '
                  f'{h:.2f}, {words[1]} {z:.2f} {v.model_unit}')
        # the GNSS left out of the fit still sit on the model: their median misfit
        self.gnss_fit_shift = 0.0
        if gin is not None and gin['index'].size and gnss is None:
            ok = np.isfinite(gin['x']) & np.isfinite(gin['value'])
            if ok.any():
                ref = self.gnss_ref()
                pred = res['model'].predict(res['best'], gin['x'][ok],
                                            offset=res['model'].offsets[ref],
                                            kind=res['model'].gnss_kind)
                self.gnss_fit_shift = float(np.median(pred - gin['value'][ok]))
        return res

    @property
    def gnss_kind(self):
        """What the GNSS of a dip-slip fault are: vertical (z) or across the
        fault (h)."""
        return 'z' if self.v.gnss_comp == 'vu' else 'h'

    def dip_side(self):
        """A' or A: the side the fault dips towards."""
        return f"{self.label}'" if self.v.dip_toward == 'end' else self.label

    def gnss_ref(self):
        """The data set the GNSS are drawn against: the first of their kind
        (dip-slip), else the first."""
        model = self.result['model']
        if model.gnss_kind is None:
            return 0
        same = [j for j, kind in enumerate(model.kinds) if kind == model.gnss_kind]
        return same[0] if same else 0

    def measure_offset(self):
        """--creep-offset: the creep rate straight from the data, as Hussain et
        al. (2016) measure it - per file of the fit, a straight line through the
        pixels within KM of the fault and WIDTH of the profile line on each
        side, and their offset at the fault (in the sense of S), its 95 % range
        from a bootstrap that resamples each side on its own. LOS files of two
        geometries or more (ascending and descending) also give it split into
        fault-parallel and vertical, as Hussain et al. split theirs."""
        v, res = self.v, self.result
        km, width = v.creep_offset
        summ = res['summary']
        x0 = summ['x0']['map'] if 'x0' in summ else 0.0
        sense = 1.0 if summ['S']['map'] >= 0 else -1.0
        files, geoms = [], []
        for j, k in enumerate(res['sets']):
            f = self.fp[k]
            u = f['x'] - x0
            near = (np.abs(u) <= km) & (np.abs(f['perp']) <= width)
            a, b = near & (u < 0), near & (u > 0)
            fit = fault_offset(u[a], f['v'][a], u[b], f['v'][b], seed=v.seed + j)
            one = dict(j=j, k=k, label=v.labels[f['file_idx']], n_a=int(a.sum()),
                       n_b=int(b.sum()), fit=fit)
            if fit is not None:
                boot = sense * fit['boot']
                one.update(rate=sense * fit['offset'], median=float(np.median(boot)),
                           lo=float(np.percentile(boot, 2.5)),
                           hi=float(np.percentile(boot, 97.5)))
                if 'los_par' in f:    # back to the LOS, for the split
                    h = float(np.median(f['los_par'][a | b]))
                    up = float(np.median(f['los_up'][a | b]))
                    geoms.append((fit['boot'] * h, fit['offset'] * h, h, up))
            files.append(one)
        split = None
        if len(geoms) >= 2:
            point = split_offsets([g[1] for g in geoms], [g[2] for g in geoms],
                                  [g[3] for g in geoms])
            draws = split_offsets(np.vstack([g[0] for g in geoms]), [g[2] for g in geoms],
                                  [g[3] for g in geoms])
            if point is not None:
                par, up = sense * draws[0], draws[1]
                split = dict(par=float(sense * point[0][0]), up=float(point[1][0]),
                             par_median=float(np.median(par)), up_median=float(np.median(up)),
                             par_lo=float(np.percentile(par, 2.5)),
                             par_hi=float(np.percentile(par, 97.5)),
                             up_lo=float(np.percentile(up, 2.5)),
                             up_hi=float(np.percentile(up, 97.5)), n=len(geoms))
        for one in files:
            if one['fit'] is None:
                print(f'WARNING: section {self.name}: {one["label"]}: too few pixels within '
                      f'{km:g} km of the fault and {width:g} km of the line for the creep '
                      f'offset ({one["n_a"]} | {one["n_b"]}, 3 a side at least)')
        return dict(km=km, width=width, x0=x0, files=files, split=split)

    def offset_lines(self):
        """The creep offset as text, as Hussain et al. write their numbers - no
        longer than the lines of the model, so the box stays narrow."""
        v, off = self.v, self.offset
        seen = ' in the LOS' if v.data_los else ''
        out = []
        many = len(off['files']) > 1
        sp = off['split']
        for one in off['files'] if sp is None else []:    # split: the files in _params
            who = f', {one["label"]}' if many else ''
            if one['fit'] is None:
                out.append(f'creep offset{who}: too few pixels near the fault')
                continue
            out.append(f'creep offset{seen}{who} = '
                       + ci_text(one['rate'], one['lo'], one['hi'], v.model_unit))
        if sp is not None:
            unit = _unit_text(v.model_unit)
            out.append('creep offset, fault-parallel = '
                       + ci_text(sp['par'], sp['par_lo'], sp['par_hi'], v.model_unit))
            out.append(f'vertical offset = {sp["up"]:+.2g}{unit} ({sp["up_lo"]:+.2g} – '
                       f"{sp['up_hi']:+.2g}, 95% CI), {self.label}' side up")
        return out

    def offset_rows(self):
        """The creep offsets as rows (name, description, value, median, lo, hi)
        of _params.csv - the value is the least-squares one, the rest from the
        bootstrap."""
        off = self.offset
        if off is None:
            return []
        where = f'within {off["km"]:g} km of the fault and {off["width"]:g} km of the line'
        rows, many = [], len(off['files']) > 1
        for one in off['files']:
            if one['fit'] is not None:
                rows.append((f'C_offset_{one["j"] + 1}' if many else 'C_offset',
                             f'creep rate from the offset of two lines {where}'
                             + (f' - {one["label"]}' if many else ''),
                             one['rate'], one['median'], one['lo'], one['hi']))
        sp = off['split']
        if sp is not None:
            rows.append(('C_offset_par', f'fault-parallel creep rate from the offsets of '
                         f'{sp["n"]} geometries', sp['par'], sp['par_median'], sp['par_lo'],
                         sp['par_hi']))
            rows.append(('U_offset', f"vertical offset ({self.label}' side up) from the "
                         f'offsets of {sp["n"]} geometries', sp['up'], sp['up_median'],
                         sp['up_lo'], sp['up_hi']))
        return rows

    def set_ref(self, j):
        """The data set (of the fit) whose offset set j is drawn on: the first
        of its panel - the first of all for a strike-slip fault."""
        if not self.dip_slip:
            return 0
        for spec in getattr(self, 'specs', None) or self.panel_specs():
            sets = self.panel_sets(spec)
            if j in sets:
                return sets[0]
        return j

    def offset_shift(self, set_index, ref=0):
        """What brings file `set_index` (of the fit) onto the offset of file
        `ref` - by default the first."""
        res = self.result
        model, best = res['model'], res['best']
        return (best[model.index(model.offsets[ref])]
                - best[model.index(model.offsets[set_index])])

    def gnss_shift(self, ref=-1):
        """What brings the GNSS onto the offset of file `ref` of the fit - by
        default the one they are drawn against; None: their own offset."""
        res = self.result
        model = res['model']
        ref = self.gnss_ref() if ref == -1 else ref
        if model.gnss:
            if ref is None:
                return 0.0
            return (res['best'][model.index(model.offsets[ref])]
                    - res['best'][model.index('ag')])
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
        self.offset = self.measure_offset() if v.creep_offset else None
        self.reset_axes()
        self.draw_fit()
        if self.ax_strain is not None:
            self.draw_strain()
        if self.ax_topo is not None:
            self.draw_topo()
        if v.line_mode:               # a window each would be too many: made for the export
            self.fig_post = None
        else:
            self.draw_posterior()
        self.mark_fault_on_map()
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return

    def drawn(self, k):
        """What is drawn of file k of the fit: its oversampled samples with
        --nlook, else its pixels."""
        f = self.fp[k]
        if self.fp_draw and f['file_idx'] in self.fp_draw:
            return self.fp_draw[f['file_idx']]
        return f

    def plot_points(self, thin=True, sets=None, ref=0):
        """(x, v) of the points to draw, per file of the fit (or of `sets`),
        each brought onto the offset of file `ref` - thinned to --max-points in
        all (thin=True)."""
        res = self.result
        sets = list(range(len(res['sets']))) if sets is None else list(sets)
        total = sum(self.drawn(res['sets'][j])['x'].size for j in sets)
        keep = 1.0
        if thin and self.v.max_points and total > self.v.max_points:
            keep = self.v.max_points / total
        rng = np.random.default_rng(0)
        out = []
        for j in sets:
            k = res['sets'][j]
            f = self.drawn(k)
            sel = (rng.random(f['x'].size) < keep) if keep < 1 else slice(None)
            out.append((k, f['x'][sel], f['v'][sel] + self.offset_shift(j, ref)))
        return out, total

    def draw_density(self, ax, sets=None, ref=0, cax=None):
        """--cmap-prof: every pixel drawn (oversampled with --nlook) as a point
        density, as swath draws it: a log colour scale up to the 98th
        percentile of the occupied cells, the densest beyond it white."""
        v = self.v
        cax = self.ax_cbar if cax is None else cax
        pts, _ = self.plot_points(thin=False, sets=sets, ref=ref)
        x = np.concatenate([p[1] for p in pts])
        y = np.concatenate([p[2] for p in pts])
        # cells as swath's: ~55 over the data (not over the room left for the
        # legend), half a --bin wide; the grid spans the whole plot
        lo, hi = np.percentile(y, [0.5, 99.5])
        den = self.density = swath_density(x, y, bin_x=self.result['bin_km'] / 2.0,
                                       bin_y=nice_step(1.1 * (hi - lo)), ylim=ax.get_ylim())
        cmap = profile_cmap(v.cmap_prof)
        counts = np.ma.masked_where(den['counts'] == 0, den['counts'])
        occupied = den['counts'][den['counts'] > 0]
        cmax = max(float(np.percentile(occupied, 98)) if occupied.size else 10.0, 10.0)
        pcm = ax.pcolormesh(den['x_edges'], den['y_edges'], counts.T, cmap=cmap,
                            shading='flat', norm=mcolors.LogNorm(vmin=1.0, vmax=cmax),
                            zorder=0.5, rasterized=True)
        cax.set_axis_on()
        cbar = self.fig.colorbar(pcm, cax=cax, extend='max')
        unit = v.model_unit.split('/')[0] if _unit_text(v.model_unit).strip() else ''
        cbar.set_label(f"Pixel density\n[count / {den['bin_x']:g} km x {den['bin_y']:.3g}"
                       f"{' ' + unit if unit else ''} cell]", fontsize=max(v.font_size - 2, 6))
        cbar.ax.tick_params(labelsize=max(v.font_size - 2, 6))
        return pcm

    def draw_fit(self):
        """The profile: the pixels as red points (--cmap-prof: a density), the
        median of each bin with its spread, the MAP model as bold black dashes
        in its 95 % and 68 % ranges, the GNSS (blue; red on a density), the
        fault dotted, A and A' in the top corners, the parameters in a box and
        the legend where they hide the fewest points - the legend beside the
        plot when it is too full for it (--legend; the "l" key) - and a title.
        A reverse or normal fault: a panel per component, the model in each."""
        v, res = self.v, self.result
        specs = getattr(self, 'specs', None) or self.panel_specs()
        axes = getattr(self, 'ax_panels', None) or [self.ax_swath]
        cbars = getattr(self, 'ax_cbars', None) or []
        handles = {}                  # the legend: one of each

        # one x range for every panel: the points drawn of every file, the GNSS
        pts, _ = self.plot_points()
        lo_x = min(float(x.min()) for _, x, _ in pts if x.size)
        hi_x = max(float(x.max()) for _, x, _ in pts if x.size)
        gin = self.gnss_in
        if (gin is not None and gin['index'].size > 0
                and bool(np.isfinite(gin['value']).any())):
            lo_x = min(lo_x, float(np.nanmin(gin['x'])))
            hi_x = max(hi_x, float(np.nanmax(gin['x'])))
        if len(res['sets']) > 1 and len(axes) > 1:    # the files first in the legend
            for k in res['sets']:
                color = POINT_COLORS[k % len(POINT_COLORS)]
                name = v.labels[self.fp[k]['file_idx']]
                handles[name] = (
                    Line2D([], [], ls='', marker='o', ms=6, mfc=color, mec='k', mew=1.0,
                           label=name) if v.cmap_prof else
                    Line2D([], [], ls='', marker='o', ms=4, color=color, label=name))
        panels = []
        for i, (ax, spec) in enumerate(zip(axes, specs)):
            panels.append(self.draw_panel(ax, spec, handles, (lo_x, hi_x), top=i == 0,
                                          cax=cbars[i] if i < len(cbars) else None))
        handles = list(handles.values())
        ax, first = axes[0], panels[0]

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
        placer = _Placer(ax, np.vstack(first['occupied']), np.concatenate(first['weights']),
                         ends)

        def make_box(loc, anchor):
            small = 5 if len(lines) > 6 else 4        # a long box: a size smaller
            box = AnchoredText('\n'.join(lines), loc=loc, bbox_to_anchor=anchor,
                               bbox_transform=ax.transAxes, borderpad=0.0, frameon=True,
                               prop=dict(fontsize=max(v.font_size - small, 6)))
            box.patch.set(boxstyle='round,pad=0.35', fc='w', ec='0.5', alpha=0.88)
            box.set_zorder(10)
            ax.add_artist(box)
            return box

        def make_legend(loc, anchor):
            # a long legend in two columns: half as tall, so it finds a free spot
            return ax.legend(handles=handles, loc=loc, bbox_to_anchor=anchor,
                             bbox_transform=ax.transAxes, fontsize=max(v.font_size - 5, 6),
                             framealpha=0.88, handlelength=1.5, borderpad=0.35,
                             labelspacing=0.3, borderaxespad=0.0,
                             ncol=2 if len(handles) > 6 else 1, columnspacing=1.0)

        room = not v.model_ylim
        placer.place(make_box, room=room)
        # the legend (--legend, the "l" key): inside where it hides the least
        # (in; auto, unless it would hide data even so), beside the plot on the
        # right (out; auto, the plot too full), or none (off)
        where = self.legend_where or getattr(v, 'model_legend', 'auto') or 'auto'
        self.legend_at = where if where in ('out', 'off') else 'in'
        if where == 'out':
            self.legend_beside(ax, handles)
        elif where in ('auto', 'in'):
            ylim = ax.get_ylim()
            legend = placer.place(make_legend, room=room, max_grow=0.25)
            best = placer.last
            hidden = np.inf if best is None else best[0]
            if where == 'auto' and not hidden <= placer.budget * placer.w.sum():
                legend.remove()       # too full: beside the plot, the plot as it was
                if best is not None:
                    placer.taken.pop()
                ax.set_ylim(ylim)
                self.legend_beside(ax, handles)
                self.legend_at = 'out'

        n_files = len(res['sets'])
        total = sum(self.fp[k]['x'].size for k in res['sets'])    # fitted: the real pixels
        what = (f'{v.fault_type} fault dipping towards {self.dip_side()}, from '
                + ' + '.join(dict.fromkeys(f['comp'] for f in self.fp))
                if self.dip_slip else COMP_TEXT.get(v.comp, v.comp))
        title = (f'{self.name}: swath ±{v.width:g} km, {total:,} pixels'
                 + (f' in {n_files} files' if n_files > 1 else '')
                 + (f' (drawn {v.nlook}x oversampled)' if self.fp_draw else '')
                 + (f', {res["gnss"]["x"].size} GNSS' if res['gnss'] is not None else '')
                 + f'\nstrike {self.strike:.0f}°, {what}')
        ax.set_title(title, fontsize=v.font_size - 1, pad=6)
        self.info = title.replace('\n', '; ')
        return

    def panel_sets(self, spec):
        """The data sets of the fit (their indices) that a panel shows."""
        return [j for j, k in enumerate(self.result['sets'])
                if self.fp[k]['file_idx'] in spec['files']]

    def draw_panel(self, ax, spec, handles, x_range, top=True, cax=None):
        """One panel of the profile: its pixels (brought onto the offset of its
        first file), their medians, the model of its component with its ranges,
        its GNSS, the fault; returns what the box and the legend must not
        hide (the top panel's)."""
        v, res = self.v, self.result
        model, summ = res['model'], res['summary']
        density = bool(v.cmap_prof)
        multi = len(res['sets']) > 1
        sets = self.panel_sets(spec)
        ref = sets[0] if sets else None
        kind = model.kinds[ref] if ref is not None else spec['kind']
        lo_x, hi_x = x_range
        span = hi_x - lo_x
        x0 = summ['x0']['map'] if 'x0' in summ else 0.0

        def add(handle):
            handles.setdefault(handle.get_label(), handle)

        def draw_fault():
            # the fault (dotted), and the traces the profile crosses (dash-dot)
            free = 'x0' in summ
            ax.axvline(x0, color='k', ls=':', lw=1.2, zorder=3)
            add(Line2D([], [], color='k', ls=':', lw=1.2,
                       label='Fault, fitted position' if free else 'Fault'))
            traces = [(c - self.x_fault) * self.along_scale for c in self.crossings]
            traces = [xc for xc in traces if lo_x <= xc <= hi_x and (free or abs(xc) > 0.05)]
            for xc in traces:
                ax.axvline(xc, color=v.fault_color, ls='-.', lw=0.9, alpha=0.75, zorder=3)
            if traces:
                add(Line2D([], [], color=v.fault_color, ls='-.', lw=0.9,
                           label='Fault trace' if free else 'Other fault trace'))

        if spec['files'] and not sets:  # its files: no pixel here, or too few bins
            ax.text(0.5, 0.5, ', '.join(v.labels[i] for i in spec['files'])
                    + ': too few pixels in this swath to fit', transform=ax.transAxes,
                    ha='center', va='center', color='0.45', fontsize=max(v.font_size - 2, 6))
            if kind == 'los':         # no line of sight to see the model through
                ax.set_xlim(lo_x - 0.01 * span, hi_x + 0.01 * span)
                draw_fault()
                return dict(occupied=[np.empty((0, 2))], weights=[np.empty(0)])

        # the pixels, every file brought onto the offset of the first: red
        # points (Hussain et al.), or - --cmap-prof - a density, drawn below
        # once the limits are known
        pts, _ = self.plot_points(sets=sets, ref=ref if ref is not None else 0)
        n_drawn = sum(p[1].size for p in pts)
        alpha = float(np.clip(1.5e4 / max(n_drawn, 1), 0.05, 0.6))
        for k, x, y in pts:
            color = POINT_COLORS[k % len(POINT_COLORS)]
            name = v.labels[self.fp[k]['file_idx']] if multi else 'InSAR'
            if density:
                if multi:             # the file of each median: its colour
                    add(Line2D([], [], ls='', marker='o', ms=6, mfc=color, mec='k', mew=1.0,
                               label=name))
                continue
            ax.scatter(x, y, s=2.5, c=color, alpha=alpha, lw=0, zorder=1, rasterized=True)
            add(Line2D([], [], ls='', marker='o', ms=4, color=color, label=name))

        # the model of this component: its 95 % and 68 % ranges and the MAP
        gin = self.gnss_in
        has_gnss = (gin is not None and gin['index'].size > 0
                    and bool(np.isfinite(gin['value']).any())
                    and (spec['kind'] is None or spec['kind'] == self.gnss_kind))
        offset = (model.offsets[ref] if ref is not None
                  else ('ag' if model.gnss else model.offsets[0]))
        xx = np.linspace(lo_x, hi_x, 801)
        xx = np.unique(np.r_[xx, x0 - 1e-6, x0 + 1e-6]) if model.creep else xx
        outer, inner = (float(b) for b in getattr(v, 'bands', None) or (95.0, 68.0))
        best, bands = model_bands(res, xx, levels=(outer, inner), offset=offset, kind=kind)
        band_lo, band_hi = bands[outer]

        # the binned medians, +/- one robust standard deviation - bigger on a
        # density, their error bars haloed in white
        bins_xy = []
        big = dict(ms=6.5, mec='k', mew=1.1, ecolor='k', elinewidth=1.0) if density else \
            dict(ms=3.2, mec='0.1', mew=0.7, ecolor='0.25', elinewidth=0.5)
        for j in sets:
            k = res['sets'][j]
            b = res['bins'][j]
            face = 'w' if not multi else POINT_COLORS[k % len(POINT_COLORS)]
            shift = self.offset_shift(j, ref)
            bars = ax.errorbar(b['x'], b['median'] + shift, yerr=b['std'],
                               fmt='o', mfc=face, capsize=0, zorder=4, **big)
            if density:
                white_halo(bars, width=2.6)
            bins_xy.append(np.c_[b['x'], b['median'] + shift])
        if sets:
            add(Line2D([], [], ls='', marker='o', mfc='w',
                       label=f'Median per {res["bin_km"]:g} km (±1σ)',
                       **(dict(ms=6.5, mec='k', mew=1.1) if density else
                          dict(ms=4.5, mec='0.1', mew=0.7))))

        if density:                   # red, as the points were; the inner darker -
            # light enough for the density under them to show
            shades = ((outer, DENSITY_BAND, 0.4), (inner, DENSITY_BAND, 0.3))
        else:                         # grey, 25 % and 40 % transparent
            shades = ((outer, BAND_COLOR, BAND_ALPHA[0]), (inner, BAND_68, BAND_ALPHA[1]))
        for i, (level, color, alpha_band) in enumerate(shades):
            ax.fill_between(xx, *bands[level], color=color, alpha=alpha_band, lw=0,
                            zorder=2 + 0.1 * i)
        ax.plot(xx, best, ls='--', color=MAP_COLOR, lw=2.2, zorder=6)
        add(Line2D([], [], ls='--', color=MAP_COLOR, lw=2.2, label='MAP model'))
        for level, color, alpha_band in shades:
            add(plt.Rectangle((0, 0), 1, 1, fc=color, alpha=alpha_band,
                              label=f'{level:g}% range'))

        # --creep: the MAP model without its creep - what the creep adds
        if model.creep:
            par = model.unpack(res['best'])
            plain = forward(xx, par['S'], par['d1'], par['a'], x0=par['x0'], rot=model.rot,
                            rot_side=model.rot_side)
            ax.plot(xx, plain, ls=':', color=MAP_COLOR, lw=1.2, zorder=5)
            add(Line2D([], [], ls=':', color=MAP_COLOR, lw=1.2, label='MAP without creep'))

        # --creep-offset: the two lines either side of the fault (Hussain et al. 2016)
        if self.offset is not None:
            off = self.offset
            halo = [patheffects.withStroke(linewidth=4.4, foreground='w')]
            for one in off['files']:
                if one['fit'] is None:
                    continue
                shift = self.offset_shift(one['j'])
                for side, (slope, icpt) in ((-1, one['fit']['line_a']),
                                            (1, one['fit']['line_b'])):
                    uu = np.array([side * off['km'], 0.0])
                    ax.plot(off['x0'] + uu, slope * uu + icpt + shift, '-', color=OFFSET_COLOR,
                            lw=2.4, zorder=6.5, path_effects=halo, solid_capstyle='butt')
            if any(one['fit'] is not None for one in off['files']):
                add(Line2D([], [], color=OFFSET_COLOR, lw=2.4,
                           label=f'Offset lines, ±{off["km"]:g} km'))

        # the GNSS: blue - red on a density - outlined in white
        if has_gnss:
            shift = self.gnss_shift(ref)
            color = GNSS_ON_DENSITY if density else GNSS_COLOR
            white_halo(ax.errorbar(gin['x'], gin['value'] + shift, yerr=gin['sigma'], fmt='o',
                                   ms=6.5, mfc=color, ecolor=color, elinewidth=1.4, capsize=3,
                                   capthick=1.2, zorder=7, **GNSS_EDGE))
            add(Line2D([], [], ls='', marker='o', ms=6, mfc=color, **GNSS_EDGE,
                       label='GNSS' + ('' if model.gnss else ' (not fitted)')))

        draw_fault()

        # limits: the bulk of the points, the model and the GNSS - room on top
        # for the box and the legend in the top panel
        ally = [y for _, _, y in pts if y.size]
        lo_y, hi_y = (np.percentile(np.concatenate(ally), [0.5, 99.5]) if ally
                      else (band_lo.min(), band_hi.max()))
        lo_y, hi_y = min(lo_y, band_lo.min()), max(hi_y, band_hi.max())
        occupied = [np.c_[x, y] for _, x, y in pts] + bins_xy + [np.c_[xx[::4], best[::4]]]
        weights = ([np.ones(x.size) for _, x, _ in pts]
                   + [np.full(b.shape[0], 25.0) for b in bins_xy]
                   + [np.full(best[::4].size, 25.0)])
        if has_gnss:
            g = gin['value'] + self.gnss_shift(ref)
            lo_y = min(lo_y, np.nanmin(g - gin['sigma']))
            hi_y = max(hi_y, np.nanmax(g + gin['sigma']))
            occupied.append(np.c_[gin['x'], g])
            weights.append(np.full(gin['x'].size, 200.0))
        pad = 0.06 * (hi_y - lo_y or 1.0)
        head = 0.25 * (hi_y - lo_y) if top else 0.0
        ax.set_ylim(v.model_ylim if v.model_ylim else (lo_y - pad, hi_y + head + pad))
        ax.set_xlim(lo_x - 0.01 * span, hi_x + 0.01 * span)
        if density and sets:
            self.draw_density(ax, sets=sets, ref=ref, cax=cax)
        return dict(occupied=occupied, weights=weights)

    def cycle_legend(self):
        """The "l" key: the legend inside the plot, outside it, or none - the
        next of these from where it is now."""
        if not self.has_data:
            return
        now = self.legend_at if self.legend_at in LEGEND_CYCLE else 'in'
        self.legend_where = LEGEND_CYCLE[(LEGEND_CYCLE.index(now) + 1) % len(LEGEND_CYCLE)]
        self.redraw()
        print(f'section {self.name}: legend {self.legend_where} '
              '("l" again: ' + ', '.join(LEGEND_CYCLE) + ' in turn)')

    def redraw(self):
        """The profile of a fitted section again, as it is - no new fit."""
        if not self.has_data:
            return
        self.reset_axes()
        self.draw_fit()
        if self.ax_strain is not None:
            self.draw_strain()
        if self.ax_topo is not None:
            self.draw_topo()
        self.fig.canvas.draw_idle()

    def legend_beside(self, ax, handles):
        """The legend outside the plot, on its right - beyond the colour bar of
        a density. The figure grows to the right to hold it: the panels keep
        their size in inches, so all that is placed in them stays as it is."""
        v, fig = self.v, self.fig
        kw = dict(handles=handles, loc='upper left', frameon=False, handlelength=1.5,
                  fontsize=max(v.font_size - 5, 6), labelspacing=0.45, borderaxespad=0.0)
        try:
            renderer = fig.canvas.get_renderer()
        except Exception:                         # noqa: BLE001 - a non-Agg canvas
            return ax.legend(bbox_to_anchor=(1.02, 1.0), bbox_transform=ax.transAxes, **kw)
        width, height = fig.get_size_inches()
        boxes = [(a, a.get_position()) for a, _ in self.home]
        probe = ax.legend(bbox_to_anchor=(1.0, 1.0), bbox_transform=fig.transFigure, **kw)
        need = probe.get_window_extent(renderer).width / fig.dpi          # inches
        probe.remove()
        right = max(box.x1 for _, box in boxes) * width                     # inches
        for cax in getattr(self, 'ax_cbars', None) or [self.ax_cbar]:
            if cax is not None and cax.axison:
                right = max(right, cax.get_tightbbox(renderer).x1 / fig.dpi)
        gap, margin = 0.18, 0.12
        grown = max(width, right + gap + need + margin)
        if grown > width + 1e-6:
            if getattr(self, '_narrow', None) is None:
                self._narrow = (width, height)
            fig.set_size_inches(grown, height, forward=True)
            for a, box in boxes:                  # the same inches, in the wider figure
                a.set_position([box.x0 * width / grown, box.y0, box.width * width / grown,
                                box.height])
        return ax.legend(bbox_to_anchor=((right + gap) / grown, ax.get_position().y1),
                         bbox_transform=fig.transFigure, **kw)

    def parameter_lines(self):
        """The parameters as Hussain et al. (2016) write them: the MAP value, its
        unit, and the 95 % range - the slip rate as a rate, with its sense."""
        v, res = self.v, self.result
        summ = res['summary']
        if self.dip_slip:
            return self.dip_slip_lines()
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
            lines.append(f'{NAMES[name]}{seen} = {ci_text(val, lo, hi, unit)}{extra}')
        d = res['derived']['strain']
        lines.append(f'strain rate{" in the LOS" if v.data_los else ""} = {d["map"]:.2f} '
                     f'{strain_unit(v.model_unit)} ({d["lo"]:.2f} – {d["hi"]:.2f}, 95% CI)')
        if self.offset is not None:
            lines += self.offset_lines()
        return lines

    def dip_slip_lines(self):
        """The parameters of a reverse or normal fault: the slip rate along the
        dip, with the fault type; the locking depth; the dip, towards A or A';
        the far-field rates of the hanging wall against the footwall -
        shortening and uplift (reverse), extension and subsidence (normal)."""
        v, res = self.v, self.result
        summ, model = res['summary'], res['model']
        lines = []
        for name in ('S', 'd1', 'dip', 'x0'):
            if name == 'dip':
                text = (ci_text(summ['dip']['map'], summ['dip']['lo'], summ['dip']['hi'], 'deg')
                        if 'dip' in summ else f'{model.dip_fixed:g}° (fixed)')
                lines.append(f'dip = {text}, towards {self.dip_side()}')
                continue
            if name not in summ:
                continue
            s = summ[name]
            unit = v.model_unit if UNITS[name] == 'mm/yr' else UNITS[name]
            extra = f', {v.fault_type}' if name == 'S' else ''
            lines.append(f'{NAMES[name]} = {ci_text(s["map"], s["lo"], s["hi"], unit)}{extra}')
        for key, word in zip(('h_rate', 'v_rate'), RATE_WORDS[v.fault_type]):
            d = res['derived'][key]
            lines.append(f'{word} = ' + ci_text(d['map'], d['lo'], d['hi'], v.model_unit))
        return lines

    def derived_rows(self):
        """What is derived from the parameters, as rows (name, description,
        unit, MAP, median, lo, hi) of _params.csv: the strain rate at the
        fault (strike-slip), the far-field rates (reverse, normal) - and a
        fixed dip."""
        v, res = self.v, self.result
        if not self.dip_slip:
            d = res['derived']['strain']
            return [('strain', 'shear strain rate at the fault |S|/(pi d1)',
                     'microstrain/yr' if v.model_unit == 'mm/yr' else f'{v.model_unit}/km',
                     d['map'], d['median'], d['lo'], d['hi'])]
        words = RATE_WORDS[v.fault_type]
        rows = []
        if res['model'].dip_fixed is not None:
            dip = res['model'].dip_fixed
            rows.append(('dip', 'dip (fixed)', 'deg', dip, dip, dip, dip))
        for key, word, what in zip(('h_rate', 'v_rate'), words,
                                   ('S cos(dip) - across the fault', 'S sin(dip) - vertical')):
            d = res['derived'][key]
            rows.append((key, f'{word} {what}; the hanging wall against the footwall far off',
                         v.model_unit, d['map'], d['median'], d['lo'], d['hi']))
        return rows

    def dip_slip_model_lines(self):
        """The model of a reverse or normal fault, for _params.txt."""
        v, model = self.v, self.result['model']
        hanging = self.dip_side()
        moves = ('up the dip, over the footwall' if v.fault_type == 'reverse'
                 else 'down the dip, away from the footwall')
        h_word, z_word = RATE_WORDS[v.fault_type]
        return [
            f'fault type     : {v.fault_type} - dipping towards {hanging}; the hanging wall '
            f'(the {hanging} side) moves {moves}',
            'model          : locked from the surface down to d1, slipping at S along its dip '
            'below (the back-slip of',
            '                 Savage 1983; Segall 2010, eqs. 3.71-3.72), the trace at x0:',
            '                   u = t (x - x0),  z = u / d1 - cot(dip)',
            '                   across the fault (+ towards '
            f"{self.label}'): h(x) = -t s S/pi [cos(dip) atan(z) + (sin(dip) - z cos(dip)) "
            '/ (1 + z^2)] + a',
            '                   vertical (+ up)                 : v(x) =    s S/pi [sin(dip) '
            'atan(z) + (cos(dip) + z sin(dip)) / (1 + z^2)] + a',
            '                   a line of sight                 : ph h(x) + pu v(x) + a',
            f'                 s = {model.sense:+d} ({v.fault_type}), t = {model.toward:+d} '
            f'(dipping towards {hanging}); far off, the hanging wall moves S cos(dip) '
            f'across ({h_word.split()[0]}) and S sin(dip) vertically ({z_word.split()[0]})',
            '                 x: km from the fault, perpendicular to it; + towards '
            f"{self.label}'",
            f'data           : per file, below, in {v.model_unit}'
            + (f' (read in {v.disp_unit})' if v.disp_unit != v.model_unit else '')
            + '; ve / vn assume no motion along the fault',
        ]

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

    def draw_posterior(self, window=True):
        """The posterior, as Hussain et al. (2016, central NAF, Fig. 9): the joint
        distribution of every pair of parameters as a 2D histogram (jet; each panel scaled to
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
        if not window:                   # for the export only: no window
            from matplotlib.figure import Figure
            self.fig_post = Figure(figsize=size)
        elif self.fig_post is None or not plt.fignum_exists(self.fig_post.number):
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
        if self.fig_post is None and v.line_mode:
            self.draw_posterior(window=False)
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
            for name, what, unit, val, med, lo, hi in self.derived_rows():
                fid.write(f'{name},{what},{unit},{val:.6g},{med:.6g},{lo:.6g},{hi:.6g}\n')
            for name, what, val, med, lo, hi in self.offset_rows():
                what = what.replace(',', ';')         # one column: a label may hold a comma
                fid.write(f'{name},{what},{v.model_unit},{val:.6g},{med:.6g},{lo:.6g},'
                          f'{hi:.6g}\n')
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
        if self.dip_slip:            # the far-field rates, sample by sample
            extra, extra_names = np.column_stack(model.rates(samples)), ['h_rate', 'v_rate']
        else:
            i_s, i_d = model.index('S'), model.index('d1')
            extra = np.abs(samples[:, i_s]) / (np.pi * samples[:, i_d])
            extra_names = ['strain']
        np.savetxt(out, np.column_stack([samples, extra]), fmt='%.6g', delimiter=',',
                   header=','.join(res['names'] + extra_names), comments='')
        written.append(out)

        # 4. the bins, with the model
        out = f'{base}_bins.csv'
        rows = []
        for j, k in enumerate(res['sets']):
            b = res['bins'][j]
            pred = model.predict(res['best'], b['x'], offset=model.offsets[j],
                                 kind=model.kinds[j])
            shift = self.offset_shift(j, self.set_ref(j))
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
                f['dist'], f['x'], f['perp'], f['value'], f['v'],
                f['v'] + self.offset_shift(j, self.set_ref(j)),
                f['row'], f['col'], f['lon'], f['lat'], np.full(f['x'].size, f['file_idx'])]))
        np.savetxt(out, np.vstack(blocks), delimiter=',', comments='',
                   fmt=['%.5f', '%.5f', '%.5f', '%.6g', '%.6g', '%.6g', '%d', '%d',
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
        for name, what, unit, val, med, lo, hi in self.derived_rows():
            rng = f'[{lo:.4g}, {hi:.4g}]'
            if name == 'strain':
                what = 'shear strain rate at the fault, |S| / (pi d1)'
            lines.append(f'{name:<10s}  {val:<11.4g}  {rng:<24s}  {med:<11.4g}  {what} [{unit}]')
        lines.append('')
        if self.sense_known:
            lines.append(f'sense of slip  : {slip_sense(summ["S"]["map"], self.strike, self.prof_az)}'
                         f' (S > 0 = right-lateral with this strike)')
        if self.dip_slip:
            lines += self.dip_slip_model_lines()
        else:
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
            if self.dip_slip:
                comp = f['comp']
                if comp == 'vu':
                    about = ', --comp vu: the vertical velocity, as it is'
                elif comp in ('ve', 'vn'):
                    word = 'east' if comp == 've' else 'north'
                    about = (f', --comp {comp}: the {word} velocity / {f["proj"]:.2f} (the {word} '
                             'component of the fault normal) - the motion across the fault')
                else:                         # a line of sight: what it sees of each
                    ph, pu = f['kind']
                    about = (f', --comp los: the LOS velocity, as it is - it sees {ph:+.2f} of '
                             f'the motion across the fault and {pu:+.2f} of the vertical (the '
                             'medians of the swath)')
            else:
                about = (f', LOS sees {f["proj"]:.2f} of the fault-parallel motion'
                         if v.comp == 'los' else '')
            lines.append(f'                 {f["x"].size:,} pixels, {res["bins"][j]["x"].size} '
                         f'bins, offset {model.offsets[j]}{about}')
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
            f'strike         : {self.strike:.2f} deg - {self.strike_how}'
            + (' (the right-hand rule: the fault dips to its right)' if self.dip_slip else '')
            + f'; the profile is {self.obliquity:.1f} deg off the fault normal',
            f'bins           : {res["bin_km"]:g} km, the median of each (3+ pixels)'
            + (f'; the real pixels - --nlook {v.nlook} oversampled what is drawn only'
               if v.nlook > 1 else ''),
            'priors         : uniform, ' + ', '.join(
                f'{n} {model.bounds[n][0]:.4g} to {model.bounds[n][1]:.4g}'
                for n in res['names'] if n in ('S', 'd1', 'dip', 'C', 'd2', 'x0'))
            + (f'; dip fixed at {model.dip_fixed:g} deg'
               if getattr(model, 'dip_fixed', None) is not None else '')
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
        off = self.offset
        if off is not None:
            lines += [
                f'creep offset   : as Hussain et al. (2016) - per file, a straight line '
                f'through the pixels within {off["km"]:g} km of the fault (at {off["x0"]:.3f} '
                f'km) and {off["width"]:g} km of the profile line on each side, and their '
                'offset at the fault, in the sense of S;',
                '                 95% range: 2000 bootstrap draws resampling each side on '
                'its own, the pixels taken as independent']
            for one in off['files']:
                if one['fit'] is None:
                    lines.append(f'                 {one["label"]}: too few pixels '
                                 f'({one["n_a"]} | {one["n_b"]})')
                    continue
                lines.append(f'                 {one["label"]}: {one["rate"]:.3f} '
                             f'({one["lo"]:.3f} to {one["hi"]:.3f}) {v.model_unit}, '
                             f'{one["n_a"]} | {one["n_b"]} pixels')
            sp = off['split']
            if sp is not None:
                lines.append(f'                 split of the {sp["n"]} LOS offsets: '
                             f'fault-parallel {sp["par"]:.3f} ({sp["par_lo"]:.3f} to '
                             f'{sp["par_hi"]:.3f}), vertical {sp["up"]:+.3f} '
                             f'({sp["up_lo"]:+.3f} to {sp["up_hi"]:+.3f}) {v.model_unit}, '
                             f"+ = the {self.label}' side up")
        return '\n'.join(lines) + '\n'


##############################################################################
class ModelViewer(LineTools, swathProfileViewer):
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

    def on_key(self, event):
        if event.key in ('l', 'L'):           # the legend of the profile under the mouse
            figure = getattr(event.canvas, 'figure', None)
            for sec in self.sections:
                if sec.fig is not None and sec.fig is figure:
                    sec.cycle_legend()
            return
        if self.line_mode:
            self.on_line_key(event)
            return
        super().on_key(event)

    def __init__(self, inps, view_cmd):
        super().__init__(inps, view_cmd)
        self._los = {}
        self.gnss_comp_used = None
        self.data_los = False            # LOS data, fitted as they are (--comp auto)
        # --spacing: a line along the fault, a profile across it every SPACING km
        self.vertices, self.line_artists = [], []
        self.line, self.line_done = None, False
        self.fig_along = None

    ##---------- a line of profiles (--spacing), as stat
    @property
    def line_mode(self):
        """--spacing: the profiles come from a line along the fault."""
        return getattr(self, 'spacing', None) is not None

    def plot(self):
        if not self.line_mode:
            return super().plot()
        self.make_map()
        if self.vertices:
            self.draw_line_on_map()
            self.run_line()
        if self.save_fig:
            if any(s.has_data for s in self.sections):
                self.save_outputs()
            elif self.disp_fig:
                print(f'--save: nothing fitted yet - "e" exports into {export_dir(self.outdir)}')
        if self.disp_fig:
            print(self.help_msg())
            plt.show()
        return

    def select_point(self, event):
        """--spacing: a left click adds a point to the line along the fault."""
        if not self.line_mode:
            return super().select_point(event)
        if (event.inaxes is not self.ax_img or event.xdata is None or event.button != 1
                or self._toolbar_busy()):
            return
        if self.line_done:                    # a click after Enter starts a new line
            self.clear_line()
        yx = self.click_rowcol(event)
        self.vertices.append(yx)
        self.vprint(f'line point {len(self.vertices)}: y/x = {yx[0]:.2f}, {yx[1]:.2f}')
        self.draw_line_on_map()
        return

    def on_line_key(self, event):
        if event.key in ('enter', 'return'):
            if self.line_done:
                print('the profiles of this line are fitted - "n" (or a click) for a new line')
            else:
                self.run_line()
        elif event.key in ('u', 'backspace'):
            if self.vertices and not self.line_done:
                self.vertices.pop()
                self.draw_line_on_map()
        elif event.key == 'n':
            self.clear_line()
            print('new line - left-click its points along the fault, then press Enter.')
        elif event.key == 'e':
            self.save_outputs()
        elif event.key in ('=', '+'):
            print('--spacing: the profiles come from the line along the fault (without '
                  '--spacing, "=" starts the next profile)')
        return

    def run_line(self):
        """A profile across the line every --spacing km, --length km either side
        of it, at right angles to its direction over +/- --width - each a
        section of its own (A-A', B-B' ...), fitted; then the model along it."""
        from .stat import grid_mirrored, polyline_frame, profile_axes, sample_points, to_yx

        if len(self.vertices) < 2:
            print('the line needs at least two points.')
            return
        try:
            frame = polyline_frame(self.atr, self.vertices,
                                   mirrored=grid_mirrored(self.atr, self.lookup))
            dist, points, normals = sample_points(frame, self.spacing, smooth_km=self.width)
        except ValueError as err:
            print(f'cannot use this line: {err}')
            return
        self.line = dict(frame=frame, dist_km=dist, points=points, normals=normals,
                         length_km=frame['cum'][-1] * frame['scale'])
        self.line_done = True
        half = self.length_km / frame['scale']
        num = len(dist)
        print(f"line: {self.line['length_km']:.1f} km - {num} profile(s) every "
              f'{self.spacing:g} km, {self.length_km:g} km either side of it, the pixels '
              f'within {self.width:g} km of each')
        if self.spacing < 2 * self.width:
            print(f'NOTE: spacing {self.spacing:g} km < 2 x width: neighbouring swaths '
                  'overlap and share pixels')
        for k in range(num):
            across = profile_axes(normals[k])[0]                # towards the right side
            start = to_yx(frame, points[k] - across * half)
            end = to_yx(frame, points[k] + across * half)
            sec = self.SECTION_CLASS(self, len(self.sections))
            self.sections.append(sec)
            sec.start_yx = [float(start[0]), float(start[1])]
            sec.end_yx = [float(end[0]), float(end[1])]
            sec.on_line, sec.line_point, sec.line_dist = True, k, float(dist[k])
            sec.draw_on_map()
            self.vprint(f'\n----- section {sec.name} ({k + 1} of {num}), '
                        f'{dist[k]:.1f} km along the line')
            try:
                sec.draw()
            except (ValueError, IndexError) as err:   # one bad profile: go on with the rest
                print(f'WARNING: section {sec.name}: {err}')
                sec.clear_data()
            if self.disp_fig and sec.fig is not None:  # made after plt.show(): show it
                try:
                    sec.fig.show()
                except Exception:        # noqa: BLE001 - a backend without windows
                    pass
        fitted = [s for s in self.sections if s.on_line and s.has_data]
        print(f'line: {len(fitted)} of {num} profile(s) fitted')
        self.fig_map.canvas.draw_idle()
        self.draw_along()
        return

    def clear_line(self):
        """Forget the line and its profiles - their windows close."""
        for sec in self.sections:
            sec.clear_map_artists()
            for art in sec.map_marks:
                try:
                    art.remove()
                except (ValueError, NotImplementedError):
                    pass
            for fig in (sec.fig, sec.fig_post):
                if fig is not None:
                    plt.close(fig)
        self.sections = []
        self.vertices, self.line, self.line_done = [], None, False
        self._remove(self.line_artists)
        if self.fig_along is not None:
            plt.close(self.fig_along)
            self.fig_along = None
        self.fig_map.canvas.draw_idle()
        return

    def along_names(self):
        """The parameters drawn along the line, in this order."""
        fitted = [s for s in self.sections if s.on_line and s.has_data]
        names = set().union(*(s.result['summary'] for s in fitted)) if fitted else set()
        order = ('S', 'd1', 'dip') if self.dip_slip else ('S', 'd1', 'C', 'd2')
        return [n for n in order if n in names]

    def draw_along(self):
        """The model along the line: a panel per parameter - the MAP of every
        profile with its 95 % range, the profiles named over the points."""
        names = self.along_names()
        secs = [s for s in self.sections if s.on_line]
        if not names or not secs:
            return
        size = (max(self.fig_size[0], 8.0), 1.0 + 2.1 * len(names))
        if self.fig_along is None or not plt.fignum_exists(self.fig_along.number):
            self.fig_along = _new_figure('Along the line', size)
            self.fig_along.canvas.mpl_connect('key_press_event', self.on_key)
            if self.disp_fig:
                try:
                    self.fig_along.show()
                except Exception:        # noqa: BLE001
                    pass
        fig = self.fig_along
        fig.clf()
        axes = fig.subplots(len(names), 1, sharex=True, squeeze=False)[:, 0]
        fig.subplots_adjust(left=0.12, right=0.97, top=1 - 0.5 / size[1],
                            bottom=0.55 / size[1], hspace=0.08)
        dist = np.array([s.line_dist for s in secs])
        unit = _unit_text(self.model_unit).strip()
        sense = (not self.dip_slip and all(s.sense_known for s in secs if s.has_data))
        hi_top = None
        for ax, name in zip(axes, names):
            val, lo, hi = (np.array([s.result['summary'][name][key] if s.has_data and name in
                                     s.result['summary'] else np.nan for s in secs])
                           for key in ('map', 'lo', 'hi'))
            hi_top = hi if hi_top is None else hi_top        # the top of the first panel
            yerr = np.vstack([val - lo, hi - val])
            cont = ax.errorbar(dist, val, yerr=yerr, fmt='o-', color='k', ms=5, lw=1.2,
                               elinewidth=1.0, capsize=3, mfc='w', mew=1.2, zorder=4)
            cont.lines[0].set_gid(f'along_{name}')
            u = {'mm/yr': unit, 'km': 'km', 'deg': '°'}.get(UNITS.get(name, ''), '')
            label = describe(name).capitalize() + (f' ({u})' if u else '')
            if name == 'S' and sense:
                label += '\n+ right-lateral'
            ax.set_ylabel(label, fontsize=self.font_size)
            if name in ('S', 'C'):
                ax.axhline(0, color='0.5', lw=0.6, zorder=1)
            ax.tick_params(which='both', direction='in', labelsize=max(self.font_size - 1, 6),
                           top=True, right=True)
        lo, hi = axes[0].get_ylim()                  # room for the names over the points
        axes[0].set_ylim(lo, hi + 0.18 * (hi - lo))
        for sec, x, y, top in zip(secs, dist, axes[0].lines[0].get_ydata(), hi_top):
            if np.isfinite(y):
                axes[0].annotate(sec.label, xy=(x, top), xytext=(0, 4),
                                 textcoords='offset points', ha='center', va='bottom',
                                 fontsize=max(self.font_size - 3, 6), fontweight='bold')
        step = self.spacing
        axes[-1].set_xlim(min(0.0, dist.min() - 0.6 * step),
                          max(self.line['length_km'], dist.max() + 0.6 * step))
        axes[-1].set_xlabel('Distance along the line (km)', fontsize=self.font_size)
        axes[0].set_title(f'{len(secs)} profiles every {step:g} km, \u00b1{self.length_km:g} km'
                          f' \u00d7 \u00b1{self.width:g} km | MAP and 95% range',
                          fontsize=max(self.font_size - 1, 6))
        fig.canvas.draw_idle()
        return

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
        if self.line_mode:
            self.vertices = self.cli_line()
        return

    @property
    def dip_slip(self):
        """A reverse or a normal fault (`viewinsar model reverse|normal`)."""
        return getattr(self, 'fault_type', 'strike') != 'strike'

    def resolve_comp(self):
        """--comp auto: los when the LOS geometry of every file is known and no
        file is named as a component of a decomposition; as-is otherwise. A
        reverse or normal fault: what each file is (resolve_dip_comps)."""
        from .cli import check_model_gnss

        if self.dip_slip:
            return self.resolve_dip_comps()
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
            if named and NAMED_AS.get(named[0]) == 'vu':
                raise ValueError(
                    f'{why}: the vertical velocity tells of a reverse or a normal fault, not of '
                    'a strike-slip one - viewinsar model reverse ... (or normal), or --comp '
                    'as-is to fit it as it is all the same')
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
            self.comps = [self.comp] * len(self.file)
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
        if self.comp == 'as-is' and not was_auto:
            # given as they are: LOS data, unless a file is named as a component
            self.data_los = not any(self.named_as(f) for f in self.file)
        self.plotted = self.what_is_plotted()
        return

    def resolve_dip_comps(self):
        """reverse / normal: what each file is - vu (vertical), ve / vn (east /
        north: across the fault) or los (as it is); auto: by the file's name,
        else los, its line of sight known (not the nominal angles of a name)."""
        from .cli import check_model_gnss

        comps, said = list(self.comps), []
        for i, fname in enumerate(self.file):
            base = os.path.basename(fname)
            if comps[i] == 'auto':
                named = self.named_as(fname)
                if named == 'horizontal':
                    raise ValueError(f'{base} is named as a horizontal velocity, of no known '
                                     'direction: say which with --comp ve or vn (one per file)')
                comps[i] = {'vu': 'vu', 'east': 've', 'north': 'vn'}.get(named, 'los')
                said.append(f'{self.labels[i] if self.labels else base}: {comps[i]}'
                            + (' (its name)' if named else ''))
            if comps[i] != 'los':
                continue
            try:
                geom = self._los_geometry(i)
            except ValueError as err:
                raise ValueError(f'{base}, a LOS velocity: {err} - give --geom (or --los-inc and '
                                 '--los-az), or say what it is with --comp vu, ve or vn (one per '
                                 'file)') from None
            if geom.kind == 'nominal':
                raise ValueError(f'{base}, a LOS velocity: only nominal LOS angles from its name, '
                                 'too rough to see the model through - give --geom (MintPy '
                                 'geometryGeo.h5, LiCSBAS GEOCml*/, LiCSAR *.geo.E.tif) or '
                                 '--los-inc and --los-az; or say what it is with --comp vu, ve '
                                 'or vn (one per file)')
        self.comps = comps
        self.comp = comps[0] if len(set(comps)) == 1 else 'mixed'
        self.data_los = False
        if said and self.print_msg:
            print('--comp auto -> ' + ', '.join(said))
        check_model_gnss(self)
        if self.gnss is not None:
            self.gnss_comp_used = self.gnss_comp
            if self.gnss_comp == 'vu' and not np.any(self.gnss.has_vu):
                print('WARNING: --vel vu: no GNSS station has a vertical velocity')
            self.vprint(f'GNSS: the {self.gnss_comp} component '
                        + ('(across the fault, + towards the end of the profile), '
                           if self.gnss_comp == 'perp' else '')
                        + f'{"fitted" if self.gnss_fit else "shown"} [{self.model_unit}]')
        self.plotted = None
        return

    @staticmethod
    def named_as(fname):
        """The component a file is named as - 'vu', 'east', 'north',
        'horizontal' - or None."""
        tokens = set(re.split(r'[^a-z0-9]+', os.path.basename(fname).lower()))
        for token in sorted(tokens):
            if token in NAMED_AS:
                return NAMED_AS[token]
        return None

    def what_is_plotted(self):
        """What the y axis says: fault-parallel velocity once the data are
        turned into it (par, ve, vn, los); the vertical (vu); the data as they
        are - the LOS, or a component by the name of its file; else None, just
        "Velocity"."""
        if self.comp in ('par', 've', 'vn', 'los'):
            return 'par'
        if self.comp == 'vu':
            return 'vu'
        if self.data_los:
            return 'los'
        named = {self.named_as(f) for f in self.file} - {None}
        return named.pop() if len(named) == 1 else None

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
        """The GNSS component that matches the data: fault-parallel (par),
        across the fault (perp: + towards the end of the profile), east, north,
        up, or the LOS of the first file."""
        comp = self.gnss_comp
        if comp in ('par', 'perp'):  # along the strike, or across the fault
            return gnss.component('los', (p_vec[0], p_vec[1], 0.0))
        if comp == 'los':
            geom = self._los_geometry(0)
            return gnss.component('los', geom.unit_vector_at(gnss.lon, gnss.lat))
        return gnss.component(comp)

    ##---------- geometry
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

    ##---------- windows
    def help_msg(self):
        if self.line_mode:
            msg = '\n' + '-' * 74
            msg += '\n  left click     : add a point to the line ALONG the fault'
            msg += (f'\n  key Enter      : a profile across it every {self.spacing:g} km '
                    f'(+/-{self.length_km:g} km,')
            msg += f'\n                   +/-{self.width:g} km wide) -> each fitted, and the model'
            msg += '\n                   along the line'
            msg += '\n  key "u" / Back : remove the last point'
            msg += '\n  key "n"        : start a new line (its profiles close)'
            msg += '\n  key "l"        : the legend of a profile: in / out / off'
            msg += '\n  key "e"        : export the map, every fit and posterior, the tables and'
            msg += f'\n                   the model along the line, into {export_dir(self.outdir)}'
            msg += '\n' + '-' * 74
            return msg
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
        if self.fig_along is not None:
            plt.close(self.fig_along)
        return

    ##---------- output
    def save_outputs(self):
        """Every section's files, the map, and - for two sections or more - one
        table of the parameters of each, A, B, C ... (slip rate along strike)."""
        written = super().save_outputs()
        fitted = [s for s in self.sections if s.has_data]
        if len(fitted) < 2:
            if self.line_mode and self.line is not None:
                written += self.save_line()
            return written
        out = os.path.join(self.outdir,
                           f'{self.COMMAND}_sections_{"".join(s.label for s in fitted)}.csv')
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        if self.dip_slip:
            names, derived = ('S', 'd1', 'dip', 'x0'), ('h_rate', 'v_rate')
        else:
            names, derived = ('S', 'd1', 'C', 'd2', 'x0'), ('strain',)
        head = (['point', 'dist_km'] if self.line_mode else []) + [
            'section', 'start_lat', 'start_lon', 'end_lat', 'end_lon', 'fault_km',
            'fault_lat', 'fault_lon', 'strike_deg',
            'fault_type' if self.dip_slip else 'sense', 'pixels', 'rate_unit']
        for name in names + derived + self.offset_names(fitted):
            head += [f'{name}_map', f'{name}_lower_95', f'{name}_upper_95']
        with open(out, 'w') as fid:
            fid.write(','.join(head) + '\n')
            for sec in fitted:
                fid.write(','.join(self._section_row(sec, names, derived)) + '\n')
        if self.dip_slip:
            print(f'saved: {out}  (one row per section; S, h_rate, v_rate in {self.model_unit}, '
                  'd1, x0 in km, the dip in degrees)')
        else:
            print(f'saved: {out}  (one row per section; S, C in {self.model_unit}, d1, d2, x0 in '
                  'km, strain in ' + ('microstrain/yr)' if self.model_unit == 'mm/yr' else
                                      f'{self.model_unit}/km)'))
        written.append(out)
        if self.line_mode and self.line is not None:
            written += self.save_line()
        return written

    @property
    def line_base(self):
        (y0, x0), (y1, x1) = self.vertices[0], self.vertices[-1]
        return os.path.join(self.outdir, f'{self.COMMAND}_line_Y{int(round(y0))}X'
                                         f'{int(round(x0))}_Y{int(round(y1))}X{int(round(x1))}_'
                                         f'sp{self.spacing:g}km_L{self.length_km:g}km_'
                                         f'w{self.width:g}km')

    def save_line(self):
        """--spacing: the model along the line, and where the line is."""
        base, written = self.line_base, []
        if self.fig_along is not None:
            out = f'{base}_along{self.fig_ext}'
            self.fig_along.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
            written.append(out)
        out = f'{base}_info.txt'
        with open(out, 'w') as fid:
            fid.write(f"line           : {self.line['length_km']:.4f} km, "
                      f'{len(self.vertices)} vertices\n')
            for i, (r, c) in enumerate(self.vertices):
                ll = self.pixel_lonlat(0, [r], [c]) if 'Y_FIRST' in self.atr or self.lookup \
                    else None
                where = (f', lat/lon {float(np.ravel(ll[1])[0]):.6f}, '
                         f'{float(np.ravel(ll[0])[0]):.6f}') if ll is not None else ''
                fid.write(f'vertex {i + 1:<8d}: y/x {r:.2f}, {c:.2f}{where}\n')
            fid.write(f'profiles       : every {self.spacing:g} km along the line (centred on '
                      f'it), {self.length_km:g} km either side of it, the pixels within '
                      f'{self.width:g} km of each; at right angles to the line over +/- '
                      f'{self.width:g} km; the fault of each = where it crosses the line\n')
            for sec in self.sections:
                if sec.on_line:
                    state = 'fitted' if sec.has_data else 'no model'
                    fid.write(f'section {sec.label:<7s}: point {sec.line_point + 1}, '
                              f'{sec.line_dist:.2f} km along the line - {state}\n')
        written.append(out)
        for path in written:
            print(f'saved: {path}')
        return written

    def _section_row(self, sec, names, derived=('strain',)):
        """One row of the _sections.csv table."""
        res, summ = sec.result, sec.result['summary']
        fixed = getattr(res['model'], 'dip_fixed', None)

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
        row = ([str(sec.line_point + 1), f'{sec.line_dist:.4f}'] if self.line_mode else [])
        row += [sec.label, *lonlat(y0, x0), *lonlat(y1, x1)]
        geom = line_geometry(self.atr, sec.start_yx, sec.end_yx, self.width)
        t = sec.x_fault / (geom['length'] * geom['dist_scale'])
        row += [f'{sec.x_fault:.4f}', *lonlat(y0 + t * (y1 - y0), x0 + t * (x1 - x0)),
                f'{sec.strike:.2f}',
                slip_sense(summ['S']['map'], sec.strike, sec.prof_az) if sec.sense_known
                else (self.fault_type if self.dip_slip else ''), str(res['n_pixels']),
                self.model_unit]
        for name in names:
            s = summ.get(name)
            if s is None and name == 'dip' and fixed is not None:
                s = dict(map=fixed, lo=fixed, hi=fixed)
            row += ([f'{s["map"]:.6g}', f'{s["lo"]:.6g}', f'{s["hi"]:.6g}'] if s
                    else ['', '', ''])
        for key in derived:
            d = res['derived'][key]
            row += [f'{d["map"]:.6g}', f'{d["lo"]:.6g}', f'{d["hi"]:.6g}']
        rows = {name: (val, lo, hi) for name, _, val, _, lo, hi in sec.offset_rows()}
        for name in self.offset_names(self.sections):
            val = rows.get(name)
            row += [f'{val[0]:.6g}', f'{val[1]:.6g}', f'{val[2]:.6g}'] if val else ['', '', '']
        return row

    def offset_names(self, sections):
        """--creep-offset: the offsets of every section, as columns."""
        names = []
        for sec in sections:
            for name, *_ in (sec.offset_rows() if sec.has_data else []):
                if name not in names:
                    names.append(name)
        return tuple(names)

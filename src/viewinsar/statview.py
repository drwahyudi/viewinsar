"""`viewinsar stat`: the creep rate along a line drawn on a fault.

As Hussain et al. (2016) measured the creep of the North Anatolian Fault:
every SPACING km along the line a short profile across it (stat.py), a
straight line through the pixels of each side, and their offset at the line
(model.fit_two_lines - their eq. 2) with a bootstrap that resamples each side
on its own. Three windows:

    the map         view.py's, as in `viewinsar swath`: the line is drawn on
                    it, and its points (black circles) and their profiles show
    the creep rate  along the line, a series per file with 1-sigma error bars
                    (descending solid, ascending dashed - their Fig. 5a); with
                    two or more LOS geometries also split into fault-parallel
                    and vertical (their Fig. 5c)
    a profile       a click on a point of the creep rate: its pixels, the
                    least-squares (BLUE, red) and bootstrap (green) lines, and
                    the PDFs of the creep rate (their Fig. S3)

    left click   add a point to the line          Enter  finish: profiles + creep
    "u" / Backspace  remove the last point        "n"    start a new line
    "e"          export the figures and the .csv files
    profile window: left / right arrow the point before / after, up / down the
                    other files
"""

import csv
import os
import warnings

import numpy as np
from matplotlib import patches, patheffects
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.lines import Line2D

from .compat import is_unitless
from .gnss import _name_tokens, _pass_of
from .lineview import LineTools
from .model import split_offsets
from .modelview import azimuth_of, keep_default_keys, local_en
from .stat import (
    creep_profiles,
    from_yx,
    grid_mirrored,
    polyline_frame,
    profile_axes,
    profile_outline,
    sample_points,
    to_yx,
)
from .swath import is_projected
from .viewer import (
    _new_figure,
    export_dir,
    swathProfileViewer,
    velocity_factor,
)

__all__ = ['StatViewer']

FINISH_KEYS = ('enter', 'return')
UNDO_KEYS = ('u', 'backspace')
PROFILE_KEYS = ('left', 'right', 'up', 'down', 'e')
# a colour per file (Okabe & Ito): black first, as the first track of Hussain et al.
SERIES_COLORS = ('k', '#E69F00', '#56B4E9', '#009E73', '#CC79A7', '#0072B2', '#D55E00',
                 '#999999')
PAR_COLOR, UP_COLOR = '#d62728', '#1f77b4'      # the split: red and blue, as their Fig. 5c
BLUE_COLOR, BOOT_COLOR = '#d62728', '#2ca02c'   # the profile: red BLUE, green bootstrap
PIXEL_COLOR = '0.45'
INK = '#222222'
FILL_ALPHA = 0.16           # a profile's half on the map: its block's colour, see-through
PICK_PIXELS = 25            # a click this close to a point of the creep rate picks it


def _gauss(x, mean, sigma):
    if not (np.isfinite(mean) and np.isfinite(sigma) and sigma > 0):
        return np.full(np.shape(x), np.nan)
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def _fmt(value, digits=2):
    """A number with a true minus sign; '-' for NaN."""
    if value is None or not np.isfinite(value):
        return '–'
    return f'{value:.{digits}f}'.replace('-', '−')


class StatViewer(LineTools, swathProfileViewer):
    """The creep rate along a line on a fault (see the module docstring).

    Example:
        from viewinsar.cli import get_view_cmd, stat_cmd_line_parse
        from viewinsar.statview import StatViewer

        inps = stat_cmd_line_parse('geo_velocity.h5 --spacing 5 --length 5'.split())
        obj = StatViewer(inps, get_view_cmd(inps))
        obj.open()
        obj.plot()
    """

    COMMAND = 'stat'
    MULTI_FILE_NOTE = 'the same profiles across the line on each, a series each'

    def __init__(self, inps, view_cmd):
        super().__init__(inps, view_cmd)
        self.vertices = []            # [row, col] of the line, full file
        self.mirrored = None          # the grid is a mirror image of the map (radar)
        self.finished = False
        self.line = None              # the points along the line, on the first file
        self.results = []             # per file: its points and profiles (or None)
        self.result = None            # the first file's
        self.split = None             # the creep split into fault-parallel and vertical
        self.split_files = []         # the files with an LOS geometry for the split
        self.geom_why = {}            # file -> why it has none (no geometry, nominal)
        self.orbits = []              # 'asc' / 'desc' / None, per file
        self.factor, self.rate_unit = 1.0, 'mm/yr'
        self.fig_stat = self.ax_stat = None
        self.axes_stat = []
        self.fig_prof = None
        self.selected = None          # (file, point) of the profile on screen
        self.line_artists, self.profile_artists, self.select_artists = [], [], []

    ##---------- setup
    def open(self):
        if not self.fig_size:
            self.fig_size = [9.0, 4.8]
        super().open()
        # radar-coded data: geometryRadar.h5 gives the lat/lon of the line and
        # of every pixel used (written to the .csv files), when it is around
        if 'Y_FIRST' not in self.atr and self.lookup is None:
            self.init_lookup()
        # left / right are those of the ground, also on a (mirrored) radar grid
        self.mirrored = grid_mirrored(self.atr, self.lookup)
        factor = velocity_factor(self.disp_unit, 'mm/yr')
        if factor is None:
            self.factor, self.rate_unit = 1.0, str(self.disp_unit or '').strip() or 'value'
            print(f'NOTE: the data are in {self.disp_unit!r}, not a velocity: the offsets are '
                  'in that unit (--data-unit mm/yr if they are velocities)')
        else:
            self.factor, self.rate_unit = factor, 'mm/yr'
            if factor != 1.0:
                self.vprint(f'creep rates in mm/yr (the data are read in {self.disp_unit})')
        self.find_geometries()
        self.vertices = self.cli_line()
        return

    def find_geometries(self):
        """The LOS geometry of every file, when there are several: two or more
        (ascending and descending) split the creep into fault-parallel and
        vertical. Nominal angles - from a name, asc or desc - are left out."""
        self.split_files, self.geom_why = [], {}
        if len(self.file) == 1 and (self.geom_file or self.los_inc):
            print('NOTE: --geom / --los-inc / --los-az split the creep of two or more files '
                  '(ascending and descending) - not used with one')
        if len(self.file) > 1:
            files, pool = self.geom_list(), self.geom_pool()
            if pool:
                print(f"--geom: {len(pool)} folder(s) for {len(self.file)} files - each file "
                      'takes its own track from them')
            for i in range(len(self.file)):
                try:
                    geom = self.find_los_geometry(i, verbose=False)
                except ValueError as err:
                    if i < len(files) and not pool:          # named, and wrong: say why
                        raise
                    self.geom_why[i] = 'no LOS geometry found'
                    print(f'   {self.labels[i]}: {str(err)[:300]}')
                    continue
                if geom.kind == 'nominal':
                    self.geom_why[i] = 'only nominal angles (from its name)'
                    continue
                self._los[i] = geom
                self.split_files.append(i)
        self.orbits = [self.orbit_of(i) for i in range(len(self.file))]
        if len(self.file) > 1:
            print('LOS geometry - for the split into fault-parallel and vertical:')
            for i in range(len(self.file)):
                pas = {'asc': 'ascending', 'desc': 'descending'}.get(self.orbits[i], '?')
                what = (self._los[i].source if i in self._los else
                        f'NONE - {self.geom_why[i]}: --geom or --los-inc / --los-az')
                print(f'   {self.labels[i]:<14s} {pas:<11s} {what}')
            if len(self.split_files) < 2:
                print('NOTE: fewer than two files with an LOS geometry: no split into '
                      'fault-parallel and vertical')
        return

    @property
    def unit_word(self):
        """' mm/yr' - or nothing for data without a unit."""
        return '' if is_unitless(self.rate_unit) or self.rate_unit == 'value' else (
            f' {self.rate_unit}')

    @property
    def unit_paren(self):
        """' (mm/yr)' - or nothing for data without a unit."""
        return f' ({self.unit_word.strip()})' if self.unit_word else ''

    @property
    def measured(self):
        """'creep rate' for velocities, 'offset' for anything else (an epoch)."""
        return 'creep rate' if self.rate_unit == 'mm/yr' else 'offset'

    @property
    def split_planned(self):
        return len(self.split_files) >= 2

    def orbit_of(self, index):
        """'asc' or 'desc': the ORBIT_DIRECTION of file `index`, else its name
        (asc, desc, T014A, 105D ...), else the look of its LOS geometry; None
        when nothing says."""
        atr = self.atr_list[index]
        direction = str(atr.get('ORBIT_DIRECTION', '')).strip().lower()
        if direction.startswith('asc'):
            return 'asc'
        if direction.startswith(('desc', 'dsc')):
            return 'desc'
        for tokens in _name_tokens(self.file[index]):
            found = _pass_of(tokens)
            if found:
                return found
        geom = self._los.get(index)
        if geom is not None:
            with np.errstate(invalid='ignore'):
                inc = geom.inc if geom.constant else np.nanmedian(geom.inc)
                az = geom.az if geom.constant else np.nanmedian(geom.az)
            east = -np.sin(np.radians(inc)) * np.sin(np.radians(az))
            right = str(atr.get('ANTENNA_SIDE', '-1')).strip() != '1'
            if np.isfinite(east) and east != 0:
                # a right-looking SAR sees from the west on an ascending pass
                return 'asc' if (east < 0) == right else 'desc'
        return None

    ##---------- main
    def plot(self):
        self.make_map()
        self.create_stat_figure()
        if self.vertices:
            self.draw_line_on_map()
            self.finish()
        if self.save_fig:
            if self.finished:
                self.save_outputs()
            elif self.disp_fig:
                print(f'--save: no line yet - "e" exports into {export_dir(self.outdir)}')
        if self.disp_fig:
            print(self.help_msg())
            plt.show()
        return

    def help_msg(self):
        msg = '\n' + '-' * 78
        msg += '\n  left click     : add a point to the line (follow the fault)'
        msg += '\n  key Enter      : finish the line -> its profiles and the creep rate'
        msg += '\n  key "u" / Back : remove the last point'
        msg += '\n  key "n"        : start a new line'
        msg += '\n  click a point of the creep rate: its profile (pixels, BLUE and bootstrap'
        msg += '\n                   lines, PDFs) - there left / right: the point before /'
        msg += '\n                   after, up / down: the other files'
        msg += f'\n  key "e"        : export the figures and the .csv files ({self.fig_ext}), into'
        msg += f'\n                   {export_dir(self.outdir)}'
        msg += f'\n  block 1 = {self.block_names[0]}: LEFT of the line, walking from its '
        msg += f'first point to its last;\n  block 2 = {self.block_names[1]}: RIGHT of it. '
        msg += 'The creep rate: block 2 - block 1 at the line'
        msg += '\n' + '-' * 78
        return msg

    def create_stat_figure(self):
        n_panel = 2 if self.split_planned else 1
        size = list(self.fig_size)
        if n_panel > 1 and self.fig_size == [9.0, 4.8]:
            size[1] = 7.4
        self.fig_stat = _new_figure('Creep rate', size)
        axes = self.fig_stat.subplots(n_panel, 1, sharex=True, squeeze=False)[:, 0]
        self.axes_stat, self.ax_stat = list(axes), axes[0]
        if n_panel == 1:
            self.fig_stat.subplots_adjust(left=0.1, right=0.97, bottom=0.13, top=0.9)
        else:
            self.fig_stat.subplots_adjust(left=0.1, right=0.97, bottom=0.08, top=0.93,
                                          hspace=0.08)
        self.fig_stat.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig_stat.canvas.mpl_connect('button_release_event', self.on_stat_click)
        self.reset_stat_axes('draw a line along the fault on the map:\n'
                             'left-click its points, then press Enter')
        return

    def reset_stat_axes(self, message=None):
        unit = self.unit_paren
        for k, ax in enumerate(self.axes_stat):
            ax.cla()
            if ax is self.axes_stat[-1]:
                ax.set_xlabel('Distance along the line (km)', fontsize=self.font_size + 1)
            ylabel = f'LOS {self.measured}' if k == 0 else self.measured.capitalize()
            ax.set_ylabel(ylabel + unit, fontsize=self.font_size + 1)
            ax.tick_params(which='both', direction='in', labelsize=self.font_size,
                           bottom=True, top=True, left=True, right=True)
        if message:
            self.ax_stat.text(0.5, 0.5, message, ha='center', va='center', color='0.55',
                              fontsize=self.font_size, transform=self.ax_stat.transAxes)
        return

    ##---------- events
    def select_point(self, event):
        """Left click: one more point of the line."""
        if (event.inaxes is not self.ax_img or event.xdata is None or event.button != 1
                or self._toolbar_busy()):
            return
        if self.finished:                       # a click after Enter starts a new line
            self.clear_line()
        yx = self.click_rowcol(event)
        self.vertices.append(yx)
        self.vprint(f'line point {len(self.vertices)}: y/x = {yx[0]:.2f}, {yx[1]:.2f}')
        self.draw_line_on_map()
        return

    def on_key(self, event):
        if event.key in FINISH_KEYS:
            self.finish()
        elif event.key in UNDO_KEYS:
            if self.vertices and not self.finished:
                self.vertices.pop()
                self.draw_line_on_map()
        elif event.key == 'n':
            self.clear_line()
            print('new line - left-click its points on the map, then press Enter.')
        elif event.key == 'e':
            self.save_outputs()
        return

    def on_stat_click(self, event):
        """A click on a point of the creep rate: its profile."""
        if (not self.finished or event.inaxes not in self.axes_stat or event.button != 1
                or event.xdata is None or self._toolbar_busy(self.fig_stat)):
            return
        picked = self.pick(event)
        if picked is not None:
            self.show_profile(*picked)
        return

    def pick(self, event):
        """(file, point) of the plotted creep rate nearest a click, or None."""
        ax = event.inaxes
        dist = self.line['dist_km']
        if ax is self.axes_stat[0]:
            series = [(i, self.creep_of(i)[0]) for i, res in enumerate(self.results)
                      if res is not None]
        elif self.split is not None:             # the split: the profile of its first file
            first = self.split['files'][0]
            series = [(first, self.split['par_mean']), (first, self.split['up_mean'])]
        else:
            return None
        best, where = None, np.array([event.x, event.y], dtype=np.float64)
        for i, ys in series:
            ok = np.isfinite(ys)
            if not ok.any():
                continue
            pts = ax.transData.transform(np.c_[dist[ok], ys[ok]])
            gap = np.hypot(*(pts - where).T)
            j = int(np.argmin(gap))
            if gap[j] <= PICK_PIXELS and (best is None or gap[j] < best[0]):
                best = (gap[j], i, int(np.flatnonzero(ok)[j]))
        return None if best is None else best[1:]

    def on_profile_key(self, event):
        """The profile window: left / right the point before / after, up / down
        the other files."""
        if self.selected is None or not self.finished:
            return
        i, k = self.selected
        n_pt = len(self.line['dist_km'])
        files = [j for j, res in enumerate(self.results) if res is not None]
        if event.key == 'right':
            self.show_profile(i, min(k + 1, n_pt - 1))
        elif event.key == 'left':
            self.show_profile(i, max(k - 1, 0))
        elif event.key in ('up', 'down') and len(files) > 1:
            j = files.index(i) + (1 if event.key == 'up' else -1)
            self.show_profile(files[j % len(files)], k)
        elif event.key == 'e':
            self.save_outputs()
        return

    def clear_line(self):
        self.vertices, self.finished = [], False
        self.line = self.result = self.split = self.selected = None
        self.results = []
        for artists in (self.line_artists, self.profile_artists, self.select_artists):
            self._remove(artists)
        self.reset_stat_axes('draw a line along the fault on the map:\n'
                             'left-click its points, then press Enter')
        if self.fig_prof is not None:
            self.fig_prof.clf()
            self.fig_prof.canvas.draw_idle()
        for fig in (self.fig_map, self.fig_stat):
            fig.canvas.draw_idle()
        return

    ##---------- compute
    def finish(self):
        """The points along the line, a profile at each, the creep - and draw them."""
        if len(self.vertices) < 2:
            print('the line needs at least two points.')
            return
        verts = np.asarray(self.vertices, dtype=np.float64)
        try:
            frame = polyline_frame(self.atr, verts, mirrored=self.mirrored)
            dist, points, normals = sample_points(frame, self.spacing, smooth_km=self.width_km)
        except ValueError as err:
            print(f'cannot use this line: {err}')
            return
        self.line = dict(frame=frame, dist_km=dist, points=points, normals=normals,
                         length_km=frame['cum'][-1] * frame['scale'])
        self.results = []
        for i in range(len(self.file)):
            try:
                frame_i, pts_i, nrm_i = self.points_on_file(i)
                profiles = creep_profiles(
                    self.data_list[i], self.box_list[i], frame_i, pts_i, nrm_i,
                    self.length_km, self.width_km, mask=self.msk_list[i],
                    keep_zero=self.keep_zero, factor=self.factor, n_boot=self.n_boot,
                    seed=(self.seed, i), min_points=self.min_pixels)
                res = dict(frame=frame_i, points=pts_i, normals=nrm_i, profiles=profiles)
            except ValueError as err:
                print(f'{self.labels[i]}: cannot use this line: {err}')
                res = None
            self.results.append(res)
        self.result = self.results[0]
        self.finished = True
        self.split = self.split_creep()
        self.report()
        self.draw_profiles_on_map()
        self.draw_stat()
        if self.fig_prof is not None and plt.fignum_exists(self.fig_prof.number):
            i, k = self.selected or (0, len(dist) // 2)     # the same point, if still there
            if i >= len(self.results) or self.results[i] is None or k >= len(dist):
                i, k = 0, len(dist) // 2
            self.show_profile(i, k)
        return

    def points_on_file(self, index):
        """The line's points on the grid of file `index`, in its own frame: the
        same places on the ground, the same directions, for every file."""
        line = self.line
        if index == 0:
            return line['frame'], line['points'], line['normals']
        verts = np.asarray(self.vertices, dtype=np.float64)
        rows, cols = self.rowcol_on_file(index, verts[:, 0], verts[:, 1])
        mirrored = grid_mirrored(self.atr_list[index]) if 'Y_FIRST' in self.atr else self.mirrored
        frame = polyline_frame(self.atr_list[index], np.c_[rows, cols], mirrored=mirrored)
        frame0 = line['frame']
        tips = line['points'] + line['normals'] / frame0['scale']        # 1 km to the left
        pts = from_yx(frame, *self.rowcol_on_file(index, *to_yx(frame0, line['points'])))
        tip = from_yx(frame, *self.rowcol_on_file(index, *to_yx(frame0, tips)))
        normals = tip - pts
        size = np.hypot(normals[:, 0], normals[:, 1])
        with np.errstate(invalid='ignore', divide='ignore'):
            normals = normals / size[:, None]
        return frame, pts, normals

    def creep_of(self, index):
        """mean, std (1 sigma), BLUE of the creep rate of file `index` at every
        point - NaN where it has none."""
        n_pt = len(self.line['dist_km'])
        out = np.full((3, n_pt), np.nan)
        res = self.results[index] if index < len(self.results) else None
        if res is not None:
            for k, prof in enumerate(res['profiles']):
                if prof['fit'] is not None:
                    out[:, k] = prof['fit']['mean'], prof['fit']['std'], prof['fit']['offset']
        return out

    def strike_vectors(self):
        """(n, 2): east and north of the unit vector along the line (first
        point -> last) at each of its points, and its azimuth (degrees from
        north, clockwise)."""
        line, frame = self.line, self.line['frame']
        along = np.c_[line['normals'][:, 1], -line['normals'][:, 0]]
        half = 0.5 / frame['scale']
        ra, ca = to_yx(frame, line['points'] - along * half)
        rb, cb = to_yx(frame, line['points'] + along * half)
        if 'Y_FIRST' in self.atr and is_projected(self.atr):
            # a projected grid: its easting / northing (grid north: within a degree or two)
            de = (cb - ca) * float(self.atr['X_STEP'])
            dn = (rb - ra) * float(self.atr['Y_STEP'])
        else:
            lla, llb = self.pixel_lonlat(0, ra, ca), self.pixel_lonlat(0, rb, cb)
            de, dn = local_en(llb[0], llb[1], lla[0], lla[1])
        vec = np.c_[de, dn]
        with np.errstate(invalid='ignore', divide='ignore'):
            vec = vec / np.hypot(vec[:, 0], vec[:, 1])[:, None]
        azim = np.array([azimuth_of(e, n) if np.isfinite(e) else np.nan for e, n in vec])
        return vec, azim

    def los_at(self, index, prof, max_pixels=4000):
        """(east, north, up) of the LOS unit vector of file `index`: the median
        over the pixels of a profile."""
        rows, cols = prof['rows'], prof['cols']
        if rows.size > max_pixels:
            step = int(np.ceil(rows.size / max_pixels))
            rows, cols = rows[::step], cols[::step]
        try:
            vec = self.los_vectors(index, rows, cols)
        except ValueError:            # the geometry is on another grid: through lon/lat
            lonlat = self.pixel_lonlat(index, rows, cols)
            if lonlat is None:
                raise
            vec = self.los_vectors(index, rows, cols, lon=lonlat[0], lat=lonlat[1])
        with warnings.catch_warnings():          # all NaN: the geometry does not reach
            warnings.simplefilter('ignore', RuntimeWarning)
            return np.array([np.nanmedian(np.asarray(c, dtype=np.float64)) for c in vec])

    def split_creep(self):
        """The LOS creep rates of two or more geometries as one fault-parallel
        and one vertical rate at every point - assuming no motion across the
        fault, as Hussain et al. (2016) split theirs where they had ascending
        AND descending tracks - draw by draw of the bootstraps, and from the
        least-squares rates.

        los_i = dh (L_i . s) + du L_i,up, s the line's direction there: the
        fault-parallel creep is -dh (+ right-lateral), the vertical du (+ the
        block 2 side up). Where there is no split, why[k] says why."""
        files = [i for i in self.split_files if self.results[i] is not None]
        if len(files) < 2:
            return None
        n_pt = len(self.line['dist_km'])
        strike, azim = self.strike_vectors()
        out = {key: np.full(n_pt, np.nan) for key in (
            'par', 'par_mean', 'par_std', 'par_median', 'par_lo', 'par_hi', 'up', 'up_mean',
            'up_std', 'up_median', 'up_lo', 'up_hi')}
        out.update(files=files, used=[[] for _ in range(n_pt)], strike=azim,
                   los={i: np.full((n_pt, 3), np.nan) for i in files}, why=[''] * n_pt)
        for k in range(n_pt):
            fitted = [i for i, res in enumerate(self.results)
                      if res is not None and res['profiles'][k]['fit'] is not None]
            draws, blue, horiz, up, seen, notes = [], [], [], [], [], []
            for i in fitted:
                if i not in files:
                    notes.append(f'{self.labels[i]} has no LOS geometry')
                    continue
                prof = self.results[i]['profiles'][k]
                try:
                    los = self.los_at(i, prof)
                except (ValueError, ImportError):      # no lon/lat for the geometry
                    los = np.full(3, np.nan)
                if not (np.all(np.isfinite(los)) and np.all(np.isfinite(strike[k]))):
                    notes.append(f'the LOS geometry of {self.labels[i]} does not reach here')
                    continue
                out['los'][i][k] = los
                draws.append(prof['fit']['boot'][:, 3])
                blue.append(prof['fit']['offset'])
                horiz.append(los[0] * strike[k, 0] + los[1] * strike[k, 1])
                up.append(los[2])
                seen.append(i)
            if len(fitted) < 2:
                out['why'][k] = ('no track has a creep rate here' if not fitted else
                                 f'only {self.labels[fitted[0]]} has a creep rate here')
                continue
            if len(seen) < 2:
                out['why'][k] = '; '.join(notes)
                continue
            east = np.array([out['los'][i][k][0] for i in seen])
            if not (east.min() < 0 < east.max()):         # all from one side
                names = ' + '.join(self.labels[i] for i in seen)
                out['why'][k] = (f'{names} all look from the {"west" if east.max() < 0 else "east"}'
                                 ' (all ascending or all descending) - the split needs both')
                continue
            point = split_offsets(blue, horiz, up)
            if point is None:
                out['why'][k] = ('the lines of sight are too alike here to tell the fault-'
                                 'parallel motion from the vertical (the fault runs along them)')
                continue
            out['used'][k] = seen
            dh, du = split_offsets(np.vstack(draws), horiz, up)
            par = -dh
            out['par'][k], out['up'][k] = -float(point[0][0]), float(point[1][0])
            for name, vals in (('par', par), ('up', du)):
                lo, med, hi = np.percentile(vals, [2.5, 50, 97.5])
                out[f'{name}_mean'][k], out[f'{name}_std'][k] = vals.mean(), vals.std(ddof=1)
                out[f'{name}_median'][k], out[f'{name}_lo'][k], out[f'{name}_hi'][k] = med, lo, hi
        return out

    def split_ranges(self):
        """[(first point, last point, what)] - the split, or why there is none,
        along the line: consecutive points alike in one range."""
        sp, out = self.split, []
        if sp is None:
            return out
        for k in range(len(self.line['dist_km'])):
            what = (sp['why'][k] if sp['why'][k] else
                    'split from ' + ' + '.join(self.labels[i] for i in sp['used'][k]))
            if out and out[-1][2] == what and out[-1][1] == k - 1:
                out[-1] = (out[-1][0], k, what)
            else:
                out.append((k, k, what))
        return out

    def report(self):
        line = self.line
        n_pt = len(line['dist_km'])
        self.vprint(f"line: {line['length_km']:.2f} km, {n_pt} profiles every "
                    f"{self.spacing:g} km, {self.length_km:g} km either side of the line, "
                    f"the pixels within {self.width_km:g} km of each")
        for i, res in enumerate(self.results):
            if res is None:
                continue
            profs = res['profiles']
            fitted = [p for p in profs if p['fit'] is not None]
            least = [min(p['n_1'], p['n_2']) for p in profs]
            self.vprint(f'   {self.labels[i]}: the creep rate at {len(fitted)} of {n_pt} '
                        f'points; pixels a side: fewest {min(least)}, most '
                        f"{max(max(p['n_1'], p['n_2']) for p in profs)}")
            if len(fitted) < n_pt:
                print(f'   {self.labels[i]}: {n_pt - len(fitted)} point(s) with fewer than '
                      f'{self.min_pixels} pixels on a side - no creep rate there')
            dropped = sum(p['dropped'] for p in profs)
            if dropped:
                print(f'   {self.labels[i]}: {dropped} pixel(s) beyond a bend of the line - '
                      'on its other side than their profile says - left out')
            if any(p['clipped'] for p in profs):
                print(f'   WARNING: {self.labels[i]}: some profiles reach past the data / '
                      'the displayed area; only the part inside is used')
        if self.result is not None:
            bends = [p['trace_km'] for p in self.result['profiles']
                     if p['trace_km'] > 0.1 * self.length_km]
            if bends:
                print(f'   NOTE: the line bends inside {len(bends)} profile(s), up to '
                      f'{max(bends):.2f} km from their centre - straighten it there, or '
                      'shorten --width')
        if self.spacing < 2 * self.width_km:
            print(f'   NOTE: spacing {self.spacing:g} km < 2 x width: neighbouring profiles '
                  'overlap and share pixels')
        sp = self.split
        if sp is not None:
            dist = line['dist_km']
            print(f"split into fault-parallel and vertical at {int(np.isfinite(sp['par']).sum())}"
                  f' of {n_pt} points:')
            for k0, k1, what in self.split_ranges():
                where = (f'{dist[k0]:.1f} km (point {k0 + 1})' if k0 == k1 else
                         f'{dist[k0]:.1f}-{dist[k1]:.1f} km (points {k0 + 1}-{k1 + 1})')
                print(f"   {where}: {what if what.startswith('split') else 'none - ' + what}")
        return

    ##---------- draw: the map
    def block_style(self, b, alpha=1.0):
        """Fill and edge of block b, on the map and in the profiles."""
        red, green, blue, _ = to_rgba(self.block_colors[b])
        return dict(facecolor=(red, green, blue, alpha), edgecolor=(red, green, blue, 1.0))

    def block_handle(self, b, lw=1.2):
        return patches.Patch(label=self.block_names[b], lw=lw,
                             **self.block_style(b, alpha=max(FILL_ALPHA, 0.35)))

    def draw_profiles_on_map(self):
        """Each profile as two halves, in the colours of their blocks, and its
        point on the line as a black circle (as Hussain et al.'s Fig. 5b)."""
        self._remove(self.profile_artists)
        self._remove(self.select_artists)
        line, frame = self.line, self.line['frame']
        halves = ((-self.length_km, 0.0), (0.0, self.length_km))
        for k in range(len(line['dist_km'])):
            for b, (x0, x1) in enumerate(halves):
                rows, cols = profile_outline(frame, line['points'][k], line['normals'][k],
                                             self.length_km, self.width_km, x0, x1)
                xs, ys = self._yx2disp(rows, cols)
                patch = patches.Polygon(np.c_[xs, ys], closed=True, lw=0.9, zorder=8.5,
                                        **self.block_style(b, alpha=FILL_ALPHA))
                patch.set_gid(f'profile_{k + 1}_b{b + 1}')
                self.profile_artists.append(self.ax_img.add_patch(patch))
        rows, cols = to_yx(frame, line['points'])
        xs, ys = self._yx2disp(rows, cols)
        halo = [patheffects.Stroke(linewidth=2.6, foreground='w'), patheffects.Normal()]
        self.profile_artists += self.ax_img.plot(xs, ys, 'o', ms=5.5, mfc='none', mec='k',
                                                 mew=1.3, zorder=11, path_effects=halo,
                                                 gid='points')
        self.label_blocks()
        self.map_legend()
        self.fig_map.canvas.draw_idle()
        return

    def label_blocks(self):
        """Name each block once, just beyond the end of the middle profile."""
        line, frame, ax = self.line, self.line['frame'], self.ax_img
        k = len(line['dist_km']) // 2
        across = profile_axes(line['normals'][k])[0]
        to_disp = lambda xy: ax.transData.transform(  # noqa: E731
            np.c_[self._yx2disp(*[np.atleast_1d(v) for v in to_yx(frame, xy)])])
        centre = to_disp(line['points'][k])[0]
        for b, sign in ((0, -1.0), (1, 1.0)):
            end = line['points'][k] + sign * across * self.length_km / frame['scale']
            out = to_disp(end)[0] - centre
            out = out / max(np.hypot(*out), 1e-9)
            ha = 'left' if out[0] > 0.35 else ('right' if out[0] < -0.35 else 'center')
            va = 'bottom' if out[1] > 0.35 else ('top' if out[1] < -0.35 else 'center')
            x, y = self._yx2disp(*[np.atleast_1d(v) for v in to_yx(frame, end)])
            self.profile_artists.append(ax.annotate(
                self.block_names[b], xy=(x[0], y[0]), textcoords='offset points',
                xytext=tuple(out * 6.0), ha=ha, va=va,
                fontsize=max(self.font_size - 2, 6), fontweight='bold', zorder=12, color=INK,
                bbox=dict(boxstyle='round,pad=0.25', lw=1.6, fc=(1, 1, 1, 0.9),
                          ec=self.block_style(b)['edgecolor'])))
        return

    def map_legend(self):
        """Which colour is which block: a key on the map, beside view.py's own."""
        handles = [self.block_handle(b, lw=1.4) for b in range(2)]
        handles.append(Line2D([], [], ls='none', marker='o', ms=5.5, mfc='none', mec='k',
                              mew=1.3, label='Profile point'))
        leg = Legend(self.ax_img, handles, [h.get_label() for h in handles], loc='best',
                     fontsize=max(self.font_size - 3, 6), framealpha=0.9, handlelength=1.6)
        leg.set_zorder(13)
        self.ax_img.add_artist(leg)
        self.profile_artists.append(leg)
        return

    def highlight_on_map(self, k):
        """The profile on screen, outlined on the map."""
        self._remove(self.select_artists)
        line, frame = self.line, self.line['frame']
        rows, cols = profile_outline(frame, line['points'][k], line['normals'][k],
                                     self.length_km, self.width_km)
        xs, ys = self._yx2disp(rows, cols)
        halo = [patheffects.Stroke(linewidth=4.2, foreground='w'), patheffects.Normal()]
        self.select_artists += self.ax_img.plot(xs, ys, '-', color='k', lw=2.0, zorder=12,
                                                path_effects=halo, gid='selected')
        self.fig_map.canvas.draw_idle()
        return

    ##---------- draw: the creep rate
    def series_style(self, index):
        return dict(color=SERIES_COLORS[index % len(SERIES_COLORS)],
                    ls='--' if self.orbits[index] == 'asc' else '-')

    def draw_stat(self):
        """The LOS creep rate along the line, a series per file (Hussain et al.
        2016, Fig. 5a); and its split into fault-parallel and vertical (Fig. 5c)."""
        self.reset_stat_axes()
        ax, dist = self.axes_stat[0], self.line['dist_km']
        handles = []
        for i, res in enumerate(self.results):
            if res is None:
                continue
            mean, std, _ = self.creep_of(i)
            style = self.series_style(i)
            cont = ax.errorbar(dist, mean, yerr=std, fmt='o', ls=style['ls'], lw=1.4, ms=4.5,
                               color=style['color'], ecolor=style['color'], elinewidth=1.0,
                               capsize=2.5, mfc='w' if style['ls'] == '--' else style['color'],
                               mec=style['color'], mew=1.1, zorder=4 + i, label=self.labels[i])
            cont.lines[0].set_gid(f'creep_{i + 1}')
            handles.append(cont)
        ax.axhline(0, color='k', lw=0.6, zorder=1)
        orbits = {self.orbits[i] for i, res in enumerate(self.results) if res is not None}
        if {'asc', 'desc'} <= orbits:          # say which line is which pass
            handles += [Line2D([], [], color='0.4', ls='-', lw=1.4, label='descending'),
                        Line2D([], [], color='0.4', ls='--', lw=1.4, label='ascending')]
        if handles:
            ax.legend(handles=handles, fontsize=max(self.font_size - 3, 6), loc='best',
                      framealpha=0.85, ncol=2 if len(handles) > 4 else 1)
        if self.split_planned:
            self.draw_split(self.axes_stat[1])
        lo = min(0.0, dist.min() - 0.6 * self.spacing)
        hi = max(self.line['length_km'], dist.max() + 0.6 * self.spacing)
        ax.set_xlim(lo, hi)
        if self.ylim:
            for one in self.axes_stat:
                one.set_ylim(self.ylim)
        names = self.block_names
        ax.set_title(f'creep = {names[1]} − {names[0]} | {len(dist)} profiles every '
                     f'{self.spacing:g} km, ±{self.length_km:g} × ±{self.width_km:g} km | '
                     f'1σ: bootstrap ({self.n_boot} draws)', fontsize=max(self.font_size - 1, 6))
        self.fig_stat.canvas.draw_idle()
        return

    def draw_split(self, ax):
        sp, dist = self.split, self.line['dist_km']
        note = self.split_note()
        if sp is None or not np.isfinite(sp['par']).any():
            ax.text(0.5, 0.5, 'no split into fault-parallel and vertical\n' + note,
                    ha='center', va='center', color='0.45', transform=ax.transAxes,
                    fontsize=max(self.font_size - 2, 6), wrap=True)
            return
        if note:
            ax.text(0.01, 0.02, note, ha='left', va='bottom', color='0.4', zorder=6,
                    fontsize=max(self.font_size - 4, 6), transform=ax.transAxes,
                    bbox=dict(boxstyle='square,pad=0.2', fc='w', ec='none', alpha=0.8))
        for name, color, label in (
                ('par', PAR_COLOR, 'Fault-parallel (+ right-lateral)'),
                ('up', UP_COLOR, f'Vertical (+ {self.block_names[1]} up)')):
            cont = ax.errorbar(dist, sp[f'{name}_mean'], yerr=sp[f'{name}_std'], fmt='o-',
                               lw=1.4, ms=4.5, color=color, ecolor=color, elinewidth=1.0,
                               capsize=2.5, zorder=4, label=label)
            cont.lines[0].set_gid(f'split_{name}')
        ax.axhline(0, color='k', lw=0.6, zorder=1)
        ax.legend(fontsize=max(self.font_size - 3, 6), loc='best', framealpha=0.85)
        return

    def split_note(self):
        """Why the split is missing at some points - a line or two for the panel."""
        sp, n_pt = self.split, len(self.line['dist_km'])
        if sp is None:
            missing = [self.labels[i] for i in sorted(self.geom_why)]
            return (f"no LOS geometry for {', '.join(missing)} (--geom or --los-inc / --los-az)"
                    if missing else '')
        gaps = [k for k in range(n_pt) if sp['why'][k]]
        if not gaps:
            return ''
        reasons = []
        nogeo = sorted({i for i in self.geom_why for k in gaps
                        if f'{self.labels[i]} has no LOS geometry' in sp['why'][k]})
        if nogeo:
            reasons.append('no LOS geometry for ' + ', '.join(self.labels[i] for i in nogeo)
                           + ' (--geom or --los-inc / --los-az)')
        for key, text in (('does not reach', 'an LOS geometry that does not reach the point'),
                          ('all look from', 'tracks of one pass only (ascending or '
                                            'descending)'),
                          ('has a creep rate here', 'one track only'),
                          ('too alike', 'lines of sight too alike')):
            if any(key in sp['why'][k] for k in gaps):
                reasons.append(text)
        return (f'no split at {len(gaps)} of {n_pt} points: ' + '; '.join(reasons)
                + ' - see the run or _split.csv')

    ##---------- draw: one profile (Hussain et al. 2016, Fig. S3)
    def show_profile(self, index, k):
        """The profile of point k of file `index`, in the profile window."""
        if not self.finished or self.results[index] is None:
            return
        if self.fig_prof is None or not plt.fignum_exists(self.fig_prof.number):
            self.fig_prof = _new_figure('Profile', [10.0, 4.0])
            self.fig_prof.canvas.mpl_connect('key_press_event', self.on_profile_key)
            keep_default_keys(self.fig_prof, PROFILE_KEYS)
            if self.disp_fig:            # made after plt.show(): show it ourselves
                try:
                    self.fig_prof.show()
                except Exception:        # noqa: BLE001 - a backend without windows
                    pass
        self.selected = (index, k)
        self.fig_prof.clf()
        self.draw_profile(self.fig_prof, index, k)
        self.highlight_on_map(k)
        self.fig_prof.canvas.draw_idle()
        return

    def _profile_lines(self, ax, fit, color, lw, label, zorder):
        """The two lines of a fit, (slope_1, slope_2, b1, creep), each on its side."""
        slope_1, slope_2, icpt, creep = fit
        xa, xb = np.array([-self.length_km, 0.0]), np.array([0.0, self.length_km])
        ax.plot(xa, slope_1 * xa + icpt, '-', color=color, lw=lw, zorder=zorder, label=label)
        ax.plot(xb, slope_2 * xb + icpt + creep, '-', color=color, lw=lw, zorder=zorder)
        return

    def draw_profile_axes(self, ax, index, k, small=False):
        """The pixels of a profile and the two fits: the BLUE (red) and the mean
        of the bootstrap draws (green), with the fault location."""
        prof = self.results[index]['profiles'][k]
        fit, fsize = prof['fit'], max(self.font_size - (4 if small else 1), 5)
        x, v = prof['x_km'], prof['values']
        ax.plot(x, v, '.', ms=2.2 if small else 2.6, alpha=0.6, zorder=2, color=PIXEL_COLOR,
                mec='none', rasterized=True, gid='pixels')
        if fit is not None:
            self._profile_lines(ax, fit['params'], BLUE_COLOR, 1.6 if small else 2.4, 'BLUE', 4)
            self._profile_lines(ax, fit['boot'].mean(axis=0), BOOT_COLOR, 1.0 if small else 1.4,
                                'Bootstrap', 5)
        ax.axvline(0, color='k', lw=0.6, ls=':', zorder=1)
        if not small:
            ax.set_xlim(-self.length_km, self.length_km)
        if v.size and not small:
            lo, hi = np.percentile(v, [0.5, 99.5])
            span = max(hi - lo, 1e-9)
            ax.set_ylim(lo - 0.08 * span, hi + 0.3 * span)
        ax.tick_params(which='both', direction='in', labelsize=fsize, top=True, right=True)
        if not small:
            ax.annotate('Fault location', xy=(0, 0.83), xycoords=('data', 'axes fraction'),
                        xytext=(0, 0.97), textcoords=('data', 'axes fraction'), ha='center',
                        va='top', fontsize=fsize,
                        arrowprops=dict(arrowstyle='-|>', color='k', lw=1.2, shrinkA=1,
                                        shrinkB=0))
            for b, (xpos, ha) in enumerate(((0.02, 'left'), (0.98, 'right'))):
                ax.text(xpos, 0.03, self.block_names[b], transform=ax.transAxes, ha=ha,
                        va='bottom', fontsize=fsize, fontweight='bold',
                        color=self.block_style(b)['edgecolor'])
        if fit is None:
            ax.text(0.5, 0.5, f'fewer than {self.min_pixels} pixels\non a side', ha='center',
                    va='center', color='0.4', fontsize=fsize, transform=ax.transAxes)
        return fit

    def draw_profile(self, fig, index, k):
        """Hussain et al. (2016), Fig. S3: (a) the profile with the BLUE and
        bootstrap lines, (b) the PDFs of the creep rate of both."""
        ax1, ax2 = fig.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.45, 1]))
        fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.84, wspace=0.25)
        unit = self.unit_paren
        prof = self.results[index]['profiles'][k]
        fit = self.draw_profile_axes(ax1, index, k)
        ax1.set_xlabel(f'Distance from the line (km), + towards {self.block_names[1]}',
                       fontsize=self.font_size)
        ax1.set_ylabel(f'LOS velocity{unit}', fontsize=self.font_size)
        if fit is not None:
            ax1.legend(fontsize=max(self.font_size - 3, 6), loc='best', framealpha=0.85)
            draws = fit['boot'][:, 3]
            ax2.hist(draws, bins=min(60, max(15, int(np.sqrt(draws.size)))), density=True,
                     color='0.82', edgecolor='0.6', lw=0.5, label='Bootstrap draws', zorder=1)
            sig = max(fit['std'], fit['offset_std'] if np.isfinite(fit['offset_std']) else 0)
            centre = np.mean([fit['mean'], fit['offset']])
            grid = np.linspace(min(draws.min(), centre - 4.5 * sig),
                               max(draws.max(), centre + 4.5 * sig), 400)
            # the numbers in the legend, as the PDFs they give
            ax2.plot(grid, _gauss(grid, fit['mean'], fit['std']), color=BOOT_COLOR, lw=2.0,
                     zorder=3, label=f"Bootstrap: {_fmt(fit['mean'])} ± {_fmt(fit['std'])}\n"
                                     f"(95%: {_fmt(fit['lo'])} – {_fmt(fit['hi'])})")
            ax2.plot(grid, _gauss(grid, fit['offset'], fit['offset_std']), color=BLUE_COLOR,
                     lw=2.0, zorder=2,
                     label=f"BLUE: {_fmt(fit['offset'])} ± {_fmt(fit['offset_std'])}")
            ax2.set_ylim(0, ax2.get_ylim()[1] * 1.55)
            ax2.legend(fontsize=max(self.font_size - 3, 6), loc='upper left', framealpha=0.9)
        else:
            ax2.text(0.5, 0.5, f'no {self.measured}', ha='center', va='center', color='0.4',
                     fontsize=self.font_size, transform=ax2.transAxes)
        ax2.set_xlabel(f'{self.measured.capitalize()}{unit}', fontsize=self.font_size)
        ax2.set_ylabel('Normalised PDF', fontsize=self.font_size)
        ax2.tick_params(which='both', direction='in', labelsize=max(self.font_size - 1, 5),
                        top=True, right=True)
        n_pt = len(self.line['dist_km'])
        fig.suptitle(f"{self.labels[index]} | point {k + 1} of {n_pt}, "
                     f"{self.line['dist_km'][k]:.1f} km along the line | "
                     f"{prof['n_1']:,} | {prof['n_2']:,} pixels ({self.block_names[0]} | "
                     f"{self.block_names[1]})", fontsize=self.font_size)
        return

    def profile_grid(self, index):
        """Every profile of file `index` in one figure, for the export."""
        from matplotlib.ticker import MaxNLocator

        profs = self.results[index]['profiles']
        num = len(profs)
        ncol = min(num, 5)
        nrow = int(np.ceil(num / ncol))
        fig = Figure(figsize=(2.6 * ncol + 0.5, 1.95 * nrow + 0.9), layout='constrained')
        axes = fig.subplots(nrow, ncol, sharex=True, sharey=True, squeeze=False).ravel()
        u = self.unit_word
        fsize = max(self.font_size - 4, 5)
        for k, ax in enumerate(axes):
            if k >= num:
                ax.set_visible(False)
                continue
            fit = self.draw_profile_axes(ax, index, k, small=True)
            ax.set_title(f"{k + 1}  ({self.line['dist_km'][k]:.1f} km)", fontsize=fsize, pad=2)
            if fit is not None:
                ax.text(0.03, 0.96, f"{_fmt(fit['mean'])} \u00b1 {_fmt(fit['std'])}{u}",
                        transform=ax.transAxes, ha='left', va='top', fontsize=fsize, zorder=6,
                        bbox=dict(boxstyle='square,pad=0.15', fc='w', ec='none', alpha=0.8))
        # the ticks inside the panels, so that neighbours do not run into each other
        ticks = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]).tick_values(-self.length_km,
                                                                            self.length_km)
        axes[0].set_xticks([t for t in ticks if abs(t) < 0.999 * self.length_km])
        axes[0].set_xlim(-self.length_km, self.length_km)
        # every panel the same y range: the 1 - 99 % of all the pixels
        allv = np.concatenate([p['values'] for p in profs if p['values'].size] or [np.zeros(1)])
        lo, hi = np.percentile(allv, [1, 99])
        span = max(hi - lo, 1e-9)
        axes[0].set_ylim(lo - 0.1 * span, hi + 0.3 * span)
        unit = self.unit_paren
        fig.supxlabel(f'Distance from the line (km), + towards {self.block_names[1]}',
                      fontsize=self.font_size - 1)
        fig.supylabel(f'LOS velocity{unit}', fontsize=self.font_size - 1)
        fig.suptitle(f'{self.labels[index]}: the profiles across the line, the creep rate '
                     '(bootstrap mean \u00b1 1\u03c3) - BLUE red, bootstrap green',
                     fontsize=self.font_size - 1)
        return fig

    def close(self):
        super().close()
        for fig in (self.fig_stat, self.fig_prof):
            if fig is not None:
                try:
                    plt.close(fig)
                except Exception:        # noqa: BLE001
                    pass
        return

    ##---------- output
    @property
    def outfile_base(self):
        (y0, x0), (y1, x1) = self.vertices[0], self.vertices[-1]
        return os.path.join(self.outdir, f'stat_Y{int(round(y0))}X{int(round(x0))}_'
                                         f'Y{int(round(y1))}X{int(round(x1))}_'
                                         f'sp{self.spacing:g}km_L{self.length_km:g}km_'
                                         f'w{self.width_km:g}km')

    def _lalo(self, rows, cols):
        """lat/lon (or northing/easting) of full-file pixels, when there are any."""
        rows, cols = np.asarray(rows, dtype=np.float64), np.asarray(cols, dtype=np.float64)
        if 'Y_FIRST' in self.atr:
            lat = float(self.atr['Y_FIRST']) + (rows + 0.5) * float(self.atr['Y_STEP'])
            lon = float(self.atr['X_FIRST']) + (cols + 0.5) * float(self.atr['X_STEP'])
            return lat, lon
        if self.lookup is not None:
            lon, lat = self.lookup.to_lonlat(rows, cols)
            return lat, lon
        return None

    def _lalo_file(self, index, rows, cols):
        """lat/lon (or northing/easting) of pixels of file `index`."""
        if index == 0 or 'Y_FIRST' not in self.atr:
            return self._lalo(rows, cols)
        atr = self.atr_list[index]
        rows, cols = np.asarray(rows, dtype=np.float64), np.asarray(cols, dtype=np.float64)
        lat = float(atr['Y_FIRST']) + (rows + 0.5) * float(atr['Y_STEP'])
        lon = float(atr['X_FIRST']) + (cols + 0.5) * float(atr['X_STEP'])
        return lat, lon

    @property
    def _tag(self):
        return ['northing', 'easting'] if is_projected(self.atr) else ['lat', 'lon']

    def save_outputs(self):
        """The creep rate, the map, the profiles, and the .csv + _info.txt."""
        if not self.finished or self.line is None:
            print('nothing to save yet - draw a line and press Enter.')
            return []
        base = self.outfile_base
        os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
        written = [self.save_creep(base)]
        if self.split is not None:
            written.append(self.save_split(base))
        written.append(self.save_values(base))
        written += self.save_info_and_figures(base)
        for path in written:
            print(f'saved: {path}')
        return written

    def save_creep(self, base):
        """A row per point and file: where, and the creep rate both ways."""
        out = f'{base}_creep.csv'
        _, azim = self.strike_vectors()
        n_pt = len(self.line['dist_km'])
        cols = (['point', 'dist_km', 'row', 'col', *self._tag, 'strike_deg', 'file', 'label',
                 'orbit', 'n_1', 'n_2', 'creep', 'creep_std', 'creep_median', 'creep_lo95',
                 'creep_hi95', 'creep_lo68', 'creep_hi68', 'creep_blue', 'creep_blue_std',
                 'slope_1', 'slope_2', 'intercept_1', 'sigma', 'trace_km', 'dropped',
                 'clipped'])
        los_cols = self.split is not None
        if los_cols:
            cols += ['los_east', 'los_north', 'los_up']
        with open(out, 'w', newline='') as fid:
            writer = csv.writer(fid)
            writer.writerow(cols)
            for i, res in enumerate(self.results):
                if res is None:
                    continue
                prow, pcol = to_yx(res['frame'], res['points'])
                lalo = self._lalo_file(i, prow, pcol)
                for k in range(n_pt):
                    prof, fit = res['profiles'][k], res['profiles'][k]['fit']
                    where = [lalo[0][k], lalo[1][k]] if lalo else [np.nan, np.nan]
                    row = [k + 1, self.line['dist_km'][k], prow[k], pcol[k], *where, azim[k],
                           i + 1, self.labels[i], self.orbits[i] or '', prof['n_1'],
                           prof['n_2']]
                    if fit is None:
                        row += [np.nan] * 13
                    else:
                        row += [fit['mean'], fit['std'], fit['median'], fit['lo'], fit['hi'],
                                fit['lo68'], fit['hi68'], fit['offset'], fit['offset_std'],
                                fit['params'][0], fit['params'][1], fit['params'][2],
                                fit['sigma']]
                    row += [prof['trace_km'], prof['dropped'], int(prof['clipped'])]
                    if los_cols:
                        los = self.split['los'].get(i)
                        row += list(los[k]) if los is not None else [np.nan] * 3
                    writer.writerow([_cell(v) for v in row])
        return out

    def save_split(self, base):
        """A row per point: the creep split into fault-parallel and vertical."""
        sp, out = self.split, f'{base}_split.csv'
        prow, pcol = to_yx(self.line['frame'], self.line['points'])
        lalo = self._lalo(prow, pcol)
        cols = ['point', 'dist_km', *self._tag, 'strike_deg', 'files', 'par', 'par_mean',
                'par_std', 'par_median', 'par_lo95', 'par_hi95', 'up', 'up_mean', 'up_std',
                'up_median', 'up_lo95', 'up_hi95', 'why_not']
        with open(out, 'w', newline='') as fid:
            writer = csv.writer(fid)
            writer.writerow(cols)
            for k in range(len(self.line['dist_km'])):
                where = [lalo[0][k], lalo[1][k]] if lalo else [np.nan, np.nan]
                used = '+'.join(self.labels[i] for i in sp['used'][k])
                row = [k + 1, self.line['dist_km'][k], *where, sp['strike'][k], used]
                for name in ('par', 'up'):
                    row += [sp[name][k], sp[f'{name}_mean'][k], sp[f'{name}_std'][k],
                            sp[f'{name}_median'][k], sp[f'{name}_lo'][k], sp[f'{name}_hi'][k]]
                row.append(sp['why'][k])
                writer.writerow([_cell(v) for v in row])
        return out

    def save_values(self, base):
        """Every pixel used: its file, point, block, place in the profile, value."""
        out = f'{base}_values.csv'
        cols = ['file', 'point', 'block', 'x_km', 'y_km', 'value', 'row', 'col', *self._tag]
        blocks = []
        for i, res in enumerate(self.results):
            if res is None:
                continue
            for k, prof in enumerate(res['profiles']):
                num = prof['values'].size
                if not num:
                    continue
                ll = self._lalo_file(i, prof['rows'], prof['cols'])
                ll = ll if ll else (np.full(num, np.nan), np.full(num, np.nan))
                blocks.append(np.column_stack([
                    np.full(num, i + 1), np.full(num, k + 1),
                    np.where(prof['x_km'] < 0, 1, 2), prof['x_km'], prof['y_km'],
                    prof['values'], prof['rows'], prof['cols'], ll[0], ll[1]]))
        table = np.vstack(blocks) if blocks else np.empty((0, len(cols)))
        fmt = ['%d', '%d', '%d', '%.5f', '%.5f', '%.6g', '%d', '%d', '%.7f', '%.7f']
        np.savetxt(out, table, delimiter=',', fmt=fmt, header=','.join(cols), comments='')
        return out

    def save_info_and_figures(self, base):
        line, written = self.line, []
        out = f'{base}_info.txt'
        vyx = np.asarray(self.vertices, dtype=np.float64)
        vll = self._lalo(vyx[:, 0], vyx[:, 1])
        tag = ','.join(self._tag)
        names = self.block_names
        with open(out, 'w') as fid:
            fid.write(f'file           : {os.path.abspath(self.file[0])}\n')
            for i, fname in enumerate(self.file[1:], start=2):
                fid.write(f'file {i:<10d}: {os.path.abspath(fname)}\n')
            fid.write(f'labels         : {self.labels}\n')
            fid.write(f"orbits         : {[o or '?' for o in self.orbits]}\n")
            fid.write(f'dataset        : {self.dset}\n')
            fid.write(f'unit           : {self.rate_unit} (the data read in {self.disp_unit})\n')
            fid.write(f"line length    : {line['length_km']:.4f} km, "
                      f'{len(self.vertices)} vertices\n')
            for i, (r, c) in enumerate(vyx):
                where = f', {tag} {vll[0][i]:.6f}, {vll[1][i]:.6f}' if vll else ''
                fid.write(f'vertex {i + 1:<8d}: y/x {r:.2f}, {c:.2f}{where}\n')
            fid.write(f"spacing        : {self.spacing} km  ({len(line['dist_km'])} points, "
                      'centred on the line)\n')
            fid.write(f'profiles       : {self.length_km} km either side of the line, the '
                      f'pixels within {self.width_km} km either side of the profile '
                      'projected onto it; perpendicular to the line over +/- '
                      f'{self.width_km} km around each point\n')
            fid.write(f'block 1        : {names[0]} - LEFT of the line, walking from its '
                      'first vertex to its last\n')
            fid.write(f'block 2        : {names[1]} - RIGHT of the line\n')
            fid.write('method         : Hussain et al. (2016), eq. 2 - v = slope_1 x + b1 '
                      '(block 1, x < 0), v = slope_2 x + b1 + creep (block 2, x > 0), x the '
                      'distance from the line along the profile, + towards block 2; the '
                      f'creep rate = {names[1]} - {names[0]} at the line\n')
            fid.write(f'uncertainty    : percentile bootstrap, {self.n_boot} draws (seed '
                      f'{self.seed}), each block resampled on its own (n1 of n1, n2 of n2 '
                      'pixels, Bekaert et al. 2015); creep = the mean of the draws, creep_std '
                      'their 1 sigma; creep_blue = least squares (BLUE), creep_blue_std its '
                      'sigma from the scatter - too narrow when the errors are correlated\n')
            fid.write(f'min pixels     : {self.min_pixels} a side\n')
            fid.write('\ncreep.csv : a row per point and file; strike_deg = the direction of '
                      'the line there, degrees clockwise from north; slope_1 / slope_2 in '
                      f'{self.rate_unit} per km; intercept_1 = block 1 at the line; sigma = '
                      'the scatter of the pixels about the lines; trace_km = how far the line '
                      'strays from the centre of the profile; dropped = pixels beyond a bend '
                      'of the line, left out\n')
            fid.write('values.csv: a row per pixel used; x_km along the profile from the '
                      'line (+ towards block 2), y_km along the line\n')
            sp = self.split
            if sp is not None:
                fid.write('split.csv : the LOS creep rates of the files '
                          + ', '.join(self.labels[i] for i in sp['files'])
                          + ' as los_i = dh (L_i . s) + du L_i,up at every point - L_i the '
                          "LOS unit vector (ground to satellite, the median over the "
                          "profile's pixels), s the line's direction - assuming no motion "
                          'across the fault (Hussain et al. 2016); par = -dh, the fault-'
                          f'parallel creep (+ right-lateral), up = du (+ {names[1]} up); par '
                          '/ up from the least-squares rates, _mean / _std / ... from the '
                          'bootstrap draws; only where ascending AND descending tracks with an '
                          'LOS geometry have a creep rate - why_not says why elsewhere\n')
                dist = self.line['dist_km']
                for k0, k1, what in self.split_ranges():
                    fid.write(f'split          : {dist[k0]:.1f}-{dist[k1]:.1f} km (points '
                              f'{k0 + 1}-{k1 + 1}): {what}\n')
            for i in range(len(self.file) if len(self.file) > 1 else 0):
                geom = self._los.get(i)
                fid.write(f'LOS geometry   : {self.labels[i]}: '
                          + (geom.source if geom is not None else
                             f"none - {self.geom_why.get(i, '?')}") + '\n')
        written.append(out)

        # the figures
        out = f'{base}{self.fig_ext}'
        self.fig_stat.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
        written.append(out)
        out = self.outfile[0] if self.outfile else f'{base}_map{self.fig_ext}'
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        self.fig_map.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
        written.append(out)
        for i, res in enumerate(self.results):
            if res is None:
                continue
            fig = self.profile_grid(i)
            out = f"{base}_profiles{'' if i == 0 else f'_{i + 1}'}{self.fig_ext}"
            fig.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
            written.append(out)
        if (self.selected is not None and self.fig_prof is not None
                and plt.fignum_exists(self.fig_prof.number)):
            i, k = self.selected
            out = f"{base}_profile{k + 1}{'' if i == 0 else f'_{i + 1}'}{self.fig_ext}"
            self.fig_prof.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
            written.append(out)
        return written


def _cell(value):
    """A .csv cell: numbers in full, NaN as an empty cell, text as it is."""
    if isinstance(value, str):
        return value
    if isinstance(value, (bool, np.bool_)):
        return int(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    value = float(value)
    return '' if not np.isfinite(value) else f'{value:.9g}'

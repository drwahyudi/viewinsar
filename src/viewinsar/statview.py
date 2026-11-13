"""`viewinsar stat`: block statistics on both sides of a line along a fault.

Two windows: the map (view.py, as in `viewinsar swath`), where the line is
drawn, and the statistics - a violin per block at every point along the line,
with strikes at its median, mean and mode, and the medians of each block
joined by a dashed line. Each block has one colour, the same fill and edge in
its circles on the map, its violins and both legends.

    left click   add a point to the line          Enter  finish: circles + violins
    "u" / Backspace  remove the last point        "n"    start a new line
    "e"          export the figures and the .csv files
"""

import os

import numpy as np
from matplotlib import patches, patheffects
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.legend import Legend
from matplotlib.lines import Line2D

from .fault import load_faults, load_faults_radar
from .stat import circle_outline, distance_to_polyline, from_yx, grid_mirrored, stat_line, to_yx
from .swath import is_projected
from .viewer import NO_LOOKUP, _new_figure, export_dir, swathProfileViewer, unit_label

__all__ = ['StatViewer']

FINISH_KEYS = ('enter', 'return')
UNDO_KEYS = ('u', 'backspace')
TINT = 0.5                  # the fill: the block's colour, half way to white
CIRCLE_ALPHA = 0.9          # the circles keep a hint of the data under them
INK = '#222222'             # text, and the mean / mode strikes
MODE_STYLE = (0, (1, 1.3))  # dotted
STRIKE = 0.45               # half length of a strike, as a fraction of the violin width
RATE_COLORS = ('#d62728', '#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b', '#e377c2',
               '#17becf')            # one per file in the -rate plot (red, blue ...)
ROBUST = 1.349                      # interquartile range / standard deviation (normal)


class StatViewer(swathProfileViewer):
    """Interactive block statistics along a line (see the module docstring).

    Example:
        from viewinsar.cli import get_view_cmd, stat_cmd_line_parse
        from viewinsar.statview import StatViewer

        inps = stat_cmd_line_parse('geo_velocity.h5 -sp 5 -r 2'.split())
        obj = StatViewer(inps, get_view_cmd(inps))
        obj.open()
        obj.plot()
    """

    COMMAND = 'stat'

    def __init__(self, inps, view_cmd):
        super().__init__(inps, view_cmd)
        self.results = []             # one stat_line result per input file
        self.vertices = []            # [row, col] of the line, full file
        self.mirrored = None          # the grid is a mirror image of the map (radar)
        self.finished = False
        self.result = None
        self.fig_stat = self.ax_stat = None
        self.line_artists, self.circle_artists = [], []

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
        self.vertices = self.cli_line()
        return

    def cli_line(self):
        """The line given on the command line, as [row, col] vertices (or [])."""
        if self.line_yx:
            return [list(pair) for pair in self.line_yx]
        if self.line_lalo:
            lats, lons = [p[0] for p in self.line_lalo], [p[1] for p in self.line_lalo]
            if 'Y_FIRST' in self.atr:
                rows, cols = self.coord.lalo2yx(lats, lons)
                return [[int(r), int(c)] for r, c in zip(np.atleast_1d(rows), np.atleast_1d(cols))]
            if self.lookup is None:
                raise ValueError(f'--line-lalo on radar-coded data: {NO_LOOKUP}')
            rows, cols = self.lookup.to_rowcol(np.array(lons), np.array(lats))
            if not np.all(np.isfinite(np.r_[rows, cols])):
                raise ValueError('--line-lalo: a point falls outside the radar frame')
            return [[float(r), float(c)] for r, c in zip(rows, cols)]
        if self.line_file:
            return self.line_from_file(self.line_file)
        return []

    def line_from_file(self, fname):
        """The longest line in a vector file, as [row, col] vertices."""
        if 'Y_FIRST' in self.atr:
            traces = load_faults([fname], self.atr, verbose=self.print_msg,
                                 data_crs=self.data_crs)
            x0, y0 = float(self.atr['X_FIRST']), float(self.atr['Y_FIRST'])
            dx, dy = float(self.atr['X_STEP']), float(self.atr['Y_STEP'])
            to_pix = lambda xy: np.c_[(xy[:, 1] - y0) / dy - 0.5, (xy[:, 0] - x0) / dx - 0.5]  # noqa: E731
        else:
            if self.lookup is None:
                raise ValueError(f'--line-file on radar-coded data: {NO_LOOKUP}')
            traces = load_faults_radar([fname], self.lookup, verbose=self.print_msg)
            to_pix = lambda xy: np.c_[xy[:, 1], xy[:, 0]]  # noqa: E731
        if not traces:
            raise ValueError(f'--line-file: no line of {os.path.basename(fname)} '
                             'falls on the data')
        longest = max(traces, key=lambda t: np.sum(np.hypot(*np.diff(t, axis=0).T)))
        if len(traces) > 1:
            print(f'--line-file: {len(traces)} lines on the data - using the longest')
        return to_pix(longest).tolist()

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
        msg = '\n' + '-' * 74
        msg += '\n  left click     : add a point to the line (follow the fault)'
        msg += '\n  key Enter      : finish the line -> circles on the map + the violins'
        msg += '\n  key "u" / Back : remove the last point'
        msg += '\n  key "n"        : start a new line'
        msg += f'\n  key "e"        : export the figures and the .csv files ({self.fig_ext}), into'
        msg += f'\n                   {export_dir(self.outdir)}'
        msg += f'\n  block 1 = {self.block_names[0]}: LEFT of the line, walking from its '
        msg += f'first point to its last;\n  block 2 = {self.block_names[1]}: RIGHT of it'
        msg += '\n' + '-' * 74
        return msg

    def create_stat_figure(self):
        n_panel = 1 if (self.rate or len(self.file) == 1) else len(self.file)
        size = list(self.fig_size)
        if n_panel > 1:
            size[1] = self.fig_size[1] * (0.45 + 0.42 * n_panel)
        self.fig_stat = _new_figure('Rate' if self.rate else 'Statistics', size)
        axes = self.fig_stat.subplots(n_panel, 1, sharex=True, squeeze=False)[:, 0]
        self.axes_stat, self.ax_stat = list(axes), axes[0]
        if n_panel == 1:
            self.fig_stat.subplots_adjust(left=0.1, right=0.97, bottom=0.25, top=0.9)
        else:
            self.fig_stat.subplots_adjust(left=0.1, right=0.97, bottom=0.25 / n_panel ** 0.8,
                                          top=1 - 0.1 / n_panel ** 0.8, hspace=0.1)
        self.fig_stat.canvas.mpl_connect('key_press_event', self.on_key)
        self.reset_stat_axes('draw a line along the fault on the map:\n'
                             'left-click its points, then press Enter')
        return

    def reset_stat_axes(self, message=None):
        unit = unit_label(self.disp_unit)
        for k, ax in enumerate(self.axes_stat):
            ax.cla()
            if ax is self.axes_stat[-1]:
                ax.set_xlabel('Distance along the line [km]', fontsize=self.font_size + 1)
            if self.rate:
                ylabel = f'{self.block_names[0]} \u2212 {self.block_names[1]}{unit}'
            elif len(self.axes_stat) > 1:
                ylabel = f'{self.labels[k]}{unit}'
            else:
                ylabel = f'LOS velocity{unit}'
            ax.set_ylabel(ylabel, fontsize=self.font_size + 1)
            ax.tick_params(which='both', direction='in', labelsize=self.font_size,
                           bottom=True, top=True, left=True, right=True)
        if message:
            self.ax_stat.text(0.5, 0.5, message, ha='center', va='center', color='0.55',
                              fontsize=self.font_size, transform=self.ax_stat.transAxes)
        return

    ##---------- events
    def _toolbar_busy(self):
        toolbar = getattr(self.fig_map.canvas, 'toolbar', None)
        return bool(toolbar is not None and getattr(toolbar, 'mode', ''))

    def select_point(self, event):
        """Left click: one more point of the line."""
        if (event.inaxes is not self.ax_img or event.xdata is None or event.button != 1
                or self._toolbar_busy()):
            return
        if self.finished:                       # a click after Enter starts a new line
            self.clear_line()
        yx = self.event_yx(event)[0]
        self.vertices.append(yx)
        self.vprint(f'line point {len(self.vertices)}: y/x = {yx}')
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

    def clear_line(self):
        self.vertices, self.finished, self.result = [], False, None
        self._remove(self.line_artists)
        self._remove(self.circle_artists)
        self.reset_stat_axes('draw a line along the fault on the map:\n'
                             'left-click its points, then press Enter')
        for fig in (self.fig_map, self.fig_stat):
            fig.canvas.draw_idle()
        return

    @staticmethod
    def _remove(artists):
        for art in artists:
            try:
                art.remove()
            except (ValueError, NotImplementedError):
                pass
        artists.clear()

    ##---------- compute + draw
    def finish(self):
        """Points, circles, statistics - and draw them."""
        if len(self.vertices) < 2:
            print('the line needs at least two points.')
            return
        verts = np.asarray(self.vertices, dtype=np.float64)
        self.results = []
        for i in range(len(self.file)):
            # the same line on this file's grid (frames, tracks, any CRS)
            rows, cols = self.rowcol_on_file(i, verts[:, 0], verts[:, 1])
            mirrored = self.mirrored
            if i and 'Y_FIRST' in self.atr:
                mirrored = grid_mirrored(self.atr_list[i])
            try:
                res = stat_line(self.data_list[i], self.atr_list[i], self.box_list[i],
                                np.c_[rows, cols], self.spacing, self.radius, self.gap,
                                mask=self.msk_list[i], keep_zero=self.keep_zero,
                                mirrored=mirrored)
            except ValueError as err:
                print(f'{self.labels[i]}: cannot use this line: {err}')
                res = None
            self.results.append(res)
        self.result = self.results[0]
        if self.result is None:
            return
        self.finished = True
        self.report()
        self.draw_circles_on_map()
        self.draw_stat()
        return

    def report(self):
        res = self.result
        n_pts = len(res['dist_km'])
        self.vprint(f"line: {res['length_km']:.2f} km, {n_pts} points every "
                    f"{self.spacing:g} km, circles r = {self.radius:g} km, "
                    f"{self.gap:g} km clear of the line")
        for b in range(2):
            counts = [s['count'] for s in res['stats'][b]]
            self.vprint(f"   {self.block_names[b]}: {sum(counts):,} pixels in {n_pts} circles "
                        f"(fewest {min(counts)}, most {max(counts)})")
        pushed = np.isfinite(res['offset_km']) & (res['offset_km'] > self.radius + self.gap + 1e-6)
        if pushed.any():
            print(f'   {pushed.sum()} circle(s) on the inside of a bend were pushed further out, '
                  f'to stay {self.gap:g} km clear of the line')
        dropped = ~np.isfinite(res['offset_km'])
        if dropped.any():
            print(f'   WARNING: {dropped.sum()} circle(s) found no spot {self.gap:g} km clear of '
                  'the line (a too tight bend) and were left out')
        if any(any(c) for c in res['clipped']):
            print('   WARNING: some circles reach past the data / the displayed area; '
                  'only the part inside is used')
        if self.spacing < 2 * self.radius:
            print(f'   NOTE: spacing {self.spacing:g} km < 2 x radius: neighbouring circles '
                  'overlap and share pixels')
        return

    def draw_line_on_map(self):
        self._remove(self.line_artists)
        if self.vertices:
            yx = np.asarray(self.vertices, dtype=np.float64)
            xs, ys = self._yx2disp(yx[:, 0], yx[:, 1])
            halo = [patheffects.Stroke(linewidth=3.6, foreground='w', alpha=0.7),
                    patheffects.Normal()]
            self.line_artists += self.ax_img.plot(xs, ys, '-', color='k', lw=1.8, zorder=9,
                                                  path_effects=halo)
            self.line_artists += self.ax_img.plot(xs, ys, 'o', ms=3.5, mfc='w', mec='k',
                                                  mew=0.8, zorder=10)
        self.fig_map.canvas.draw_idle()
        return

    def block_style(self, b, alpha=1.0):
        """Fill and edge of block b: the same in its circles, violins and legends.

        An opaque tint rather than a transparent colour, so a circle over red or
        blue data still shows the very colour of its violins.
        """
        red, green, blue, _ = to_rgba(self.block_colors[b])
        fill = tuple(1.0 - TINT * (1.0 - c) for c in (red, green, blue))
        return dict(facecolor=(*fill, alpha), edgecolor=(red, green, blue, 1.0))

    def block_handle(self, b, lw=1.2):
        return patches.Patch(label=self.block_names[b], lw=lw, **self.block_style(b))

    def draw_circles_on_map(self):
        self._remove(self.circle_artists)
        res, frame = self.result, self.result['frame']
        for b in range(2):
            for k in range(len(res['dist_km'])):
                cen = res['centres'][b, k]
                if not np.all(np.isfinite(cen)):
                    continue
                rows, cols = circle_outline(frame, cen, self.radius)
                xs, ys = self._yx2disp(rows, cols)
                # a thin white ring keeps the circle apart from the data under it
                self.circle_artists += self.ax_img.plot(xs, ys, '-', color='w', lw=3.4,
                                                        alpha=0.9, zorder=8)
                patch = patches.Polygon(np.c_[xs, ys], closed=True, lw=1.4, zorder=8.5,
                                        **self.block_style(b, alpha=CIRCLE_ALPHA))
                self.circle_artists.append(self.ax_img.add_patch(patch))
        # the points along the line, in red as in the sketch
        rows, cols = to_yx(frame, res['points'])
        xs, ys = self._yx2disp(rows, cols)
        self.circle_artists += self.ax_img.plot(xs, ys, 'o', ms=4.5, color='red', mec='k',
                                                mew=0.5, zorder=11)
        self.label_blocks()
        self.map_legend()
        self.fig_map.canvas.draw_idle()
        return

    def label_blocks(self):
        """Name each block once, just outside a circle near the middle of the line."""
        res, frame, ax = self.result, self.result['frame'], self.ax_img
        both = np.flatnonzero(np.all(np.isfinite(res['centres'][:, :, 0]), axis=0))
        if not both.size:
            return
        k = int(both[np.argmin(np.abs(both - (len(res['dist_km']) - 1) / 2))])
        to_disp = lambda rows, cols: ax.transData.transform(  # noqa: E731
            np.c_[self._yx2disp(np.atleast_1d(rows), np.atleast_1d(cols))])
        point = to_disp(*to_yx(frame, res['points'][k]))[0]
        points_per_pixel = 72.0 / self.fig_map.dpi
        for b in range(2):
            cen = res['centres'][b, k]
            centre = to_disp(*to_yx(frame, cen))[0]
            edge = to_disp(*circle_outline(frame, cen, self.radius, num=5))
            reach = np.max(np.hypot(*(edge - centre).T))          # the circle, on screen
            out = centre - point
            out = out / max(np.hypot(*out), 1e-9)
            ha = 'left' if out[0] > 0.35 else ('right' if out[0] < -0.35 else 'center')
            va = 'bottom' if out[1] > 0.35 else ('top' if out[1] < -0.35 else 'center')
            x, y = self._yx2disp(*[np.atleast_1d(v) for v in to_yx(frame, cen)])
            # dark text on white; the block's colour is the frame, as on the circles
            self.circle_artists.append(ax.annotate(
                self.block_names[b], xy=(x[0], y[0]), textcoords='offset points',
                xytext=tuple(out * (reach + 4) * points_per_pixel), ha=ha, va=va,
                fontsize=max(self.font_size - 2, 6), fontweight='bold', zorder=12, color=INK,
                bbox=dict(boxstyle='round,pad=0.25', lw=1.6, fc=(1, 1, 1, 0.9),
                          ec=self.block_style(b)['edgecolor'])))
        return

    def map_legend(self):
        """Which colour is which block: a key on the map, beside view.py's own."""
        handles = [self.block_handle(b, lw=1.4) for b in range(2)]
        leg = Legend(self.ax_img, handles, [h.get_label() for h in handles], loc='best',
                     fontsize=max(self.font_size - 3, 6), framealpha=0.9, handlelength=1.6)
        leg.set_zorder(13)
        self.ax_img.add_artist(leg)
        self.circle_artists.append(leg)
        return

    def draw_stat(self):
        """The violins (one panel per file), or with -rate the rates."""
        self.reset_stat_axes()
        if self.rate:
            self.draw_rate()
        else:
            for k, ax in enumerate(self.axes_stat):
                if self.results[k] is not None:
                    self.draw_violins(ax, self.results[k], first=k == 0,
                                      last=ax is self.axes_stat[-1])
        self.fig_stat.canvas.draw_idle()
        return

    def rate_of(self, res):
        """Block 1 - block 2 at every point, and its 1-sigma.

        --rate-err iqr (default): the robust spread of the pixels, sqrt(s1^2 +
        s2^2) with s = (q75 - q25) / 1.349 - conservative, as the pixels of a
        circle are correlated; sem: the standard error of the two medians,
        1.2533 sqrt(s1^2/n1 + s2^2/n2), as if they were independent.
        """
        st1, st2 = res['stats']
        get = lambda st, key: np.array([d[key] for d in st], dtype=np.float64)  # noqa: E731
        rate = get(st1, 'median') - get(st2, 'median')
        s1 = (get(st1, 'q75') - get(st1, 'q25')) / ROBUST
        s2 = (get(st2, 'q75') - get(st2, 'q25')) / ROBUST
        if self.rate_err == 'sem':
            n1 = np.maximum(get(st1, 'count'), 1)
            n2 = np.maximum(get(st2, 'count'), 1)
            err = 1.2533 * np.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2)
        else:
            err = np.hypot(s1, s2)
        return rate, err

    def draw_rate(self):
        """Block 1 - block 2 along the line, one series per file (as Hussain et al.
        2016 for the creep along the North Anatolian Fault)."""
        ax = self.ax_stat
        lengths, dists = [], []
        for k, res in enumerate(self.results):
            if res is None:
                continue
            rate, err = self.rate_of(res)
            good = np.isfinite(rate)
            color = RATE_COLORS[k % len(RATE_COLORS)]
            ax.errorbar(res['dist_km'][good], rate[good], yerr=err[good], fmt='o-', ms=5,
                        lw=1.3, color=color, ecolor=color, elinewidth=1.0, capsize=3,
                        alpha=0.95, zorder=4 + k, label=self.labels[k])
            lengths.append(res['length_km'])
            dists.append(res['dist_km'])
        ax.axhline(0, color='k', lw=0.7, zorder=1)
        step = self.spacing
        dist = np.concatenate(dists)
        ax.set_xlim(min(0.0, dist.min() - 0.6 * step), max(max(lengths), dist.max() + 0.6 * step))
        if self.ylim:
            ax.set_ylim(self.ylim)
        ax.legend(fontsize=max(self.font_size - 3, 6), loc='best', framealpha=0.85)
        res = self.result
        ax.set_title(f"rate = {self.block_names[0]} \u2212 {self.block_names[1]} | line "
                     f"{res['length_km']:.2f} km, {len(res['dist_km'])} points every {step:g} km"
                     f" | r = {self.radius:g} km, {self.gap:g} km clear | error: "
                     f"{'IQR / 1.349' if self.rate_err == 'iqr' else 'SE of the medians'}",
                     fontsize=max(self.font_size - 1, 6))
        return

    def draw_violins(self, ax, res, first=True, last=True):
        """A violin per block at every point, with strikes at its median (the
        block's colour), mean (solid) and mode (dotted); the medians of each
        block joined by a dashed line."""
        dist, step = res['dist_km'], self.spacing
        shift, width = 0.2 * step, 0.36 * step
        half = STRIKE * width
        for b in range(2):
            style = self.block_style(b)
            edge = style['edgecolor']
            xs = dist + (-shift if b == 0 else shift)
            stats = res['stats'][b]
            idx = [k for k, v in enumerate(res['values'][b]) if v.size >= 3 and np.ptp(v) > 0]
            if idx:
                pos = xs[idx]
                parts = ax.violinplot([res['values'][b][k] for k in idx], positions=pos,
                                      widths=width, points=200, showmeans=False,
                                      showmedians=False, showextrema=False)
                for body in parts['bodies']:
                    body.set_alpha(None)             # violinplot's own 0.3 would win
                    body.set_facecolor(style['facecolor'])
                    body.set_edgecolor(edge)
                    body.set_linewidth(1.0)
                # the strikes: median (the block's colour), mean (solid), mode (dotted,
                # the peak of the density - where the violin is widest)
                for key, kwargs in (('median', dict(colors=[edge], linewidths=2.4, zorder=3)),
                                    ('mean', dict(colors=INK, linewidths=1.3, zorder=4)),
                                    ('mode', dict(colors=INK, linewidths=1.8, zorder=4,
                                                  linestyles=[MODE_STYLE]))):
                    lines = ax.hlines([stats[k][key] for k in idx], pos - half, pos + half,
                                      **kwargs)
                    lines.set_gid(f'{key}_b{b + 1}')
                    lines.set_capstyle('round' if key == 'mode' else 'butt')
            # the dashed line joins the medians of this block (empty circles skipped)
            meds = np.array([st['median'] for st in stats])
            ok = np.isfinite(meds)
            ax.plot(xs[ok], meds[ok], ls='--', lw=1.3, color=edge, marker='o', ms=3.8,
                    mfc='w', mec=edge, mew=1.1, zorder=5)
        ax.axhline(0, color='k', lw=0.4)
        lo = min(0.0, dist.min() - 0.6 * step)
        hi = max(res['length_km'], dist.max() + 0.6 * step)
        ax.set_xlim(lo, hi)
        if self.ylim:
            ax.set_ylim(self.ylim)
        # the blocks by colour, the statistics by line - under the axes, in one row
        grey = '0.35'
        handles = [self.block_handle(b, lw=1.0) for b in range(2)]
        handles += [Line2D([], [], ls='--', lw=1.3, color=grey, marker='o', ms=3.8, mfc='w',
                           mec=grey, mew=1.1, label='Median'),
                    Line2D([], [], ls='-', lw=1.3, color=INK, label='Mean'),
                    Line2D([], [], ls=MODE_STYLE, lw=1.8, color=INK, dash_capstyle='round',
                           label='Mode')]
        if last:
            ax.legend(handles=handles, fontsize=max(self.font_size - 3, 6),
                      ncol=len(handles), loc='upper center', bbox_to_anchor=(0.5, -0.17),
                      frameon=False, handlelength=2.4, columnspacing=1.8)
        if first:
            ax.set_title(f"line {res['length_km']:.2f} km | {len(dist)} points every "
                         f"{step:g} km | circles r = {self.radius:g} km, {self.gap:g} km "
                         'clear of the line', fontsize=self.font_size)
        return

    def close(self):
        super().close()
        if self.fig_stat is not None:
            try:
                plt.close(self.fig_stat)
            except Exception:        # noqa: BLE001
                pass
        return

    ##---------- output
    @property
    def outfile_base(self):
        (y0, x0), (y1, x1) = self.vertices[0], self.vertices[-1]
        return os.path.join(self.outdir, f'stat_Y{int(round(y0))}X{int(round(x0))}_'
                                         f'Y{int(round(y1))}X{int(round(x1))}_'
                                         f'sp{self.spacing:g}km_r{self.radius:g}km')

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

    def save_outputs(self):
        """The violins (or rates), the map, and the statistics as .csv + _info.txt."""
        if not self.finished or self.result is None:
            print('nothing to save yet - draw a line and press Enter.')
            return []
        base = self.outfile_base
        os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
        written = []
        for i, res in enumerate(self.results):
            if res is not None:
                written += self.save_tables(i, res, base + ('' if i == 0 else f'_{i + 1}'))
        if self.rate:
            written.append(self.save_rates(base))
        written += self.save_info_and_figures(base)
        for path in written:
            print(f'saved: {path}')
        return written

    def save_rates(self, base):
        """One row per point of the first file: the rate of every file and its error."""
        res = self.result
        prow, pcol = to_yx(res['frame'], res['points'])
        lalo = self._lalo(prow, pcol)
        tag = 'northing,easting' if is_projected(self.atr) else 'lat,lon'
        cols = ['point', 'dist_km', 'row', 'col'] + (tag.split(',') if lalo else [])
        table = [np.arange(len(res['dist_km'])), res['dist_km'], prow, pcol]
        if lalo:
            table += [lalo[0], lalo[1]]
        many = len(self.results) > 1
        for k, other in enumerate(self.results):
            tail = f'_{k + 1}' if many else ''
            rate = np.full(len(res['dist_km']), np.nan)
            err = np.full(len(res['dist_km']), np.nan)
            if other is not None:
                r_k, e_k = self.rate_of(other)
                # match the points by their distance along the line (the grids differ)
                for j, d in enumerate(res['dist_km']):
                    near = int(np.argmin(np.abs(other['dist_km'] - d)))
                    if abs(other['dist_km'][near] - d) <= 0.25 * self.spacing:
                        rate[j], err[j] = r_k[near], e_k[near]
            cols += [f'rate{tail}', f'rate_err{tail}']
            table += [rate, err]
        out = f'{base}_rate.csv'
        np.savetxt(out, np.column_stack(table), delimiter=',', fmt='%.6f',
                   header=','.join(cols), comments='')
        return out

    def save_tables(self, index, res, base):
        """stats.csv and values.csv of one file."""
        tag = 'northing,easting' if is_projected(self.atr) else 'lat,lon'
        written = []

        # 1. one row per point: both blocks side by side, and their difference
        prow, pcol = to_yx(res['frame'], res['points'])
        lalo = self._lalo_file(index, prow, pcol)
        keys = ('count', 'median', 'mean', 'mode', 'std', 'q25', 'q75', 'min', 'max')
        cols = ['point', 'dist_km', 'row', 'col'] + (tag.split(',') if lalo else [])
        for b in (1, 2):
            cols += [f'{key}_b{b}' for key in keys] + [f'offset_km_b{b}']
        cols += ['median_b1_minus_b2']
        rows = []
        for k, dist in enumerate(res['dist_km']):
            row = [k, dist, prow[k], pcol[k]] + ([lalo[0][k], lalo[1][k]] if lalo else [])
            for b in range(2):
                row += [res['stats'][b][k][key] for key in keys] + [res['offset_km'][b, k]]
            row += [res['stats'][0][k]['median'] - res['stats'][1][k]['median']]
            rows.append(row)
        out = f'{base}_stats.csv'
        np.savetxt(out, np.array(rows, dtype=np.float64), delimiter=',', fmt='%.6f',
                   header=','.join(cols), comments='')
        written.append(out)

        # 2. every pixel used, with its distance from the line
        blocks = []
        for b in range(2):
            for k in range(len(res['dist_km'])):
                val = res['values'][b][k]
                if not val.size:
                    continue
                rr, cc = res['rows'][b][k], res['cols'][b][k]
                xy = from_yx(res['frame'], rr, cc)
                line_km = distance_to_polyline(res['frame'], xy) * res['frame']['scale']
                arrs = [np.full(val.size, k), np.full(val.size, b + 1),
                        np.full(val.size, res['dist_km'][k]), val, rr, cc, line_km]
                ll = self._lalo_file(index, rr, cc)
                if ll:
                    arrs += [ll[0], ll[1]]
                blocks.append(np.column_stack(arrs))
        vcols = ['point', 'block', 'dist_km', 'value', 'row', 'col', 'dist_to_line_km']
        if lalo:
            vcols += tag.split(',')
        out = f'{base}_values.csv'
        table = np.vstack(blocks) if blocks else np.empty((0, len(vcols)))
        np.savetxt(out, table, delimiter=',', fmt='%.6f', header=','.join(vcols), comments='')
        written.append(out)
        return written

    def save_info_and_figures(self, base):
        res, written = self.result, []
        tag = 'northing,easting' if is_projected(self.atr) else 'lat,lon'
        out = f'{base}_info.txt'
        vyx = np.asarray(self.vertices, dtype=np.float64)
        vll = self._lalo(vyx[:, 0], vyx[:, 1])
        with open(out, 'w') as fid:
            fid.write(f'file           : {os.path.abspath(self.file[0])}\n')
            for i, fname in enumerate(self.file[1:], start=2):
                fid.write(f'file {i:<10d}: {os.path.abspath(fname)}  (tables *_{i}.csv)\n')
            if len(self.file) > 1:
                fid.write(f'labels         : {self.labels}\n')
            fid.write(f'dataset        : {self.dset}\n')
            fid.write(f'unit           : {self.disp_unit}\n')
            fid.write(f"line length    : {res['length_km']:.4f} km, "
                      f'{len(self.vertices)} vertices\n')
            for i, (r, c) in enumerate(vyx):
                where = f', {tag} {vll[0][i]:.6f}, {vll[1][i]:.6f}' if vll else ''
                fid.write(f'vertex {i + 1:<8d}: y/x {r:.2f}, {c:.2f}{where}\n')
            fid.write(f'spacing        : {self.spacing} km  ({len(res["dist_km"])} points, '
                      'centred on the line)\n')
            fid.write(f'radius         : {self.radius} km\n')
            fid.write(f'gap            : {self.gap} km between the line and every circle\n')
            fid.write(f"normals        : the line's direction over +/- {res['smooth_km']:g} km "
                      'around each point\n')
            fid.write(f'block 1        : {self.block_names[0]} - LEFT of the line, walking '
                      'from its first vertex to its last\n')
            fid.write(f'block 2        : {self.block_names[1]} - RIGHT of the line\n')
            fid.write('\nstats.csv : one row per point; offset_km = distance from the point '
                      'to the circle centre (radius + gap, more on the inside of a bend; '
                      'NaN = no clear spot, circle left out); mode = the peak of the '
                      "values' Gaussian kernel density (Scott's rule), the widest point of "
                      'the violin\n')
            fid.write('values.csv: one row per pixel used; dist_to_line_km is never below '
                      'the gap\n')
            if self.rate:
                fid.write('rate.csv  : rate = median of block 1 - median of block 2 at every '
                          'point, per file; rate_err = '
                          + ('the robust spread sqrt(s1^2 + s2^2), s = (q75 - q25) / 1.349\n'
                             if self.rate_err == 'iqr' else
                             'the standard error of the medians, 1.2533 sqrt(s1^2/n1 + '
                             's2^2/n2)\n'))
        written.append(out)

        # 4. the figures
        out = f'{base}{self.fig_ext}'
        self.fig_stat.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
        written.append(out)
        out = self.outfile[0] if self.outfile else f'{base}_map{self.fig_ext}'
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        self.fig_map.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
        written.append(out)
        return written

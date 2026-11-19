"""The line along a fault that `viewinsar stat` and `viewinsar model --spacing`
take: clicked on the map, or given with --line-lalo, --line-yx or --line-file
(the longest line of a vector file) - as [row, col] of the first file, where
it was put, not snapped to the pixels (a pixel can be a km wide, and the
two sides of the fault meet at the line).
"""

import os

import numpy as np
from matplotlib import patheffects

from .fault import load_faults, load_faults_radar
from .viewer import NO_LOOKUP

__all__ = ['LineTools']


class LineTools:
    """Mixed into a swathProfileViewer that takes a line along a fault: the
    viewer holds the line as self.vertices ([row, col] of the first file)
    and its drawing on the map as self.line_artists."""

    def cli_line(self):
        """The line given on the command line, as [row, col] vertices (or [])."""
        if self.line_yx:
            return [[float(c) for c in pair] for pair in self.line_yx]
        if self.line_lalo:
            lats, lons = [p[0] for p in self.line_lalo], [p[1] for p in self.line_lalo]
            if 'Y_FIRST' in self.atr:
                return [list(self.geo_to_rowcol(la, lo)) for la, lo in zip(lats, lons)]
            if self.lookup is None:
                raise ValueError(f'--line-lalo on radar-coded data: {NO_LOOKUP}')
            rows, cols = self.lookup.to_rowcol(np.array(lons), np.array(lats))
            if not np.all(np.isfinite(np.r_[rows, cols])):
                raise ValueError('--line-lalo: a point falls outside the radar frame')
            return [[float(r), float(c)] for r, c in zip(rows, cols)]
        if self.line_file:
            return self.line_from_file(self.line_file)
        return []

    def geo_to_rowcol(self, lat, lon):
        """Fractional (row, col) of a lat/lon (or northing/easting) on the
        first file's grid - not snapped to a pixel."""
        row = (float(lat) - float(self.atr['Y_FIRST'])) / float(self.atr['Y_STEP']) - 0.5
        col = (float(lon) - float(self.atr['X_FIRST'])) / float(self.atr['X_STEP']) - 0.5
        return row, col

    def click_rowcol(self, event):
        """[row, col] of a click on the map, where it was."""
        if 'Y_FIRST' in self.atr and self.fig_coord == 'geo':
            return list(self.geo_to_rowcol(event.ydata, event.xdata))
        return [float(event.ydata), float(event.xdata)]           # the map is in pixels

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

    def _toolbar_busy(self, fig=None):
        toolbar = getattr((fig or self.fig_map).canvas, 'toolbar', None)
        return bool(toolbar is not None and getattr(toolbar, 'mode', ''))

    def draw_line_on_map(self):
        """The line and its vertices, black on a white halo."""
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

    @staticmethod
    def _remove(artists):
        for art in artists:
            try:
                art.remove()
            except (ValueError, NotImplementedError):
                pass
        artists.clear()

"""The interactive swath profile viewer.

The map lives in its own window, drawn by MintPy's view.py. Every cross-section
gets a window of its own: A-A', B-B', C-C' ... Press "=" on any window to start
another section; its line stays on the map, labelled, alongside the earlier
ones. The interaction follows plot_transection.py and the optional point
time-series window reuses tsview.py.
"""

import os
import string
import textwrap

import numpy as np
from matplotlib import colors as mcolors
from matplotlib import patches, patheffects
from matplotlib import pyplot as plt
from mintpy import view
from mintpy.utils import plot as pp
from mintpy.utils import readfile
from mintpy.utils import utils as ut

from .fault import load_faults, load_faults_radar, profile_crossings, to_display
from .radar import PIXEL_ATR, RadarLookup, find_lookup, has_height
from .raster import (
    crs_name,
    crs_of,
    line_to_other_grid,
    lonlat_line_to_grid,
    read_any,
    read_attr_any,
    same_crs,
)
from .swath import is_projected, line_geometry, swath_density, swath_extract
from .terrain import ATTRIBUTION, is_auto, terrain_dem, terrain_dem_radar

__all__ = ['swathProfileViewer', 'Section', 'topo_sampling']

SECTION_LABELS = string.ascii_uppercase      # A, B, C, ...
DEM_VALID = (-1000.0, 10000.0)               # plausible elevations, in meter
TOPO_MAX_SAMPLES = 4_000_000                 # DEM samples per section: plenty for a median
NO_LOOKUP = ("the data are in radar coordinates, and that needs MintPy's geometryRadar.h5 "
             '(latitude / longitude of every pixel) - none was found in inputs/, ./ or '
             '../inputs/ next to the data or here; name it with --lookup')


def topo_sampling(n_pixels, pixel_km, bin_x, budget=TOPO_MAX_SAMPLES):
    """How to sample the DEM under a section: (decimate, nlook).

    The DEM is sampled on its own terms - --nlook is for the data. A 30 m DEM
    under a 100 x 30 km corridor already holds millions of pixels, so taking
    --nlook 10 over to it asked for hundreds of millions of samples. Instead:
    thin a DEM finer than the median line needs down to `budget` samples, and
    oversample only a DEM coarser than half a distance bin, so every bin still
    gets samples.

    Parameters: n_pixels - int, native DEM pixels in the corridor's bounding box
                pixel_km - float, DEM pixel size, in km
                bin_x    - float, the along-profile bin size, in km
    Returns:    decimate - int, take every Nth DEM pixel (1 = all)
                nlook    - int, oversampling factor (1 = native)
    """
    n_pixels = max(int(n_pixels), 1)
    if n_pixels > budget:
        return int(np.ceil(np.sqrt(n_pixels / budget))), 1
    want = int(np.ceil(pixel_km / (bin_x / 2.0)))
    room = int(np.floor(np.sqrt(budget / n_pixels)))
    return 1, max(1, min(want, room))


def _ordinal(num):
    """2 -> '2nd', 3 -> '3rd', 11 -> '11th'."""
    suffix = 'th' if 10 <= num % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(num % 10, 'th')
    return f'{num}{suffix}'


def unit_label(unit):
    """' [unit]' for an axis label, or '' when the unit is unknown.

    MintPy gives a GeoTIFF without a unit the unit '1', which would otherwise
    read "LOS velocity [1]".
    """
    unit = str(unit or '').strip()
    return '' if unit.lower() in ('', '1', 'none') else f' [{unit}]'


##############################################################################
def _new_figure(name, figsize):
    """A FRESH figure with a window title.

    plt.figure('name') returns an EXISTING figure when one already carries that
    label, which would draw a second section on top of the first whenever two
    viewers live in the same session. Always make a new one, then label it.
    """
    fig = plt.figure(figsize=figsize)
    fig.set_label(name)
    try:
        fig.canvas.manager.set_window_title(name)
    except AttributeError:
        pass
    return fig


class Section:
    """One cross-section: its line on the map, its own window, its data."""

    def __init__(self, viewer, index):
        self.v = viewer
        self.index = index
        self.label = (SECTION_LABELS[index] if index < len(SECTION_LABELS)
                      else f'S{index + 1}')
        self.start_yx = self.end_yx = None
        self.start_lalo = self.end_lalo = None
        # its window
        self.fig = self.ax_swath = self.ax_cbar = self.ax_topo = None
        # its artists on the map
        self.line_obj = self.poly_obj = self.txt0 = self.txt1 = None
        # its data
        self.swath_list = None
        self.den = None
        self.topo_swath = None
        self.topo_den = None
        self.topo_error = None       # why there is no topography, shown in its panel
        self.topo_note = ''          # how the DEM was sampled
        self.clear_crossings()

    def clear_crossings(self):
        """No fault crossing known (yet): km along the line, and map x/y."""
        self.crossings = np.array([])
        self.crossing_pts = np.empty((0, 2))

    def __repr__(self):
        return f'<Section {self.name} {self.start_yx} -> {self.end_yx}>'

    @property
    def name(self):
        return f"{self.label}-{self.label}'"

    @property
    def has_data(self):
        return self.swath_list is not None and self.den is not None

    ##---------- its window
    def create_figure(self):
        """Open this section's window, with the axes it needs."""
        v = self.v
        if self.fig is not None and plt.fignum_exists(self.fig.number):
            return self.fig

        height = v.fig_size[1] + (1.7 if v.topoxsect else 0)
        self.fig = _new_figure(f'Section {self.name}', [v.fig_size[0], height])
        gs_kw = dict(left=0.13, right=0.86, bottom=0.12, top=0.88)
        if v.topoxsect:
            gsp = self.fig.add_gridspec(2, 1, height_ratios=[3.0, 1.15],
                                        hspace=0.08, **gs_kw)
            self.ax_swath = self.fig.add_subplot(gsp[0])
            self.ax_topo = self.fig.add_subplot(gsp[1], sharex=self.ax_swath)
        else:
            gsp = self.fig.add_gridspec(1, 1, **gs_kw)
            self.ax_swath = self.fig.add_subplot(gsp[0])

        # a fixed slot for the density colorbar, so redraws do not stack them
        pos = self.ax_swath.get_position()
        self.ax_cbar = self.fig.add_axes([pos.x1 + 0.015, pos.y0, 0.016, pos.height])
        self.ax_cbar.set_axis_off()

        self.fig.canvas.mpl_connect('key_press_event', v.on_key)
        self.reset_axes()
        self.ax_swath.text(0.5, 0.5, f'section {self.name}\nclick 2 points on the map',
                           ha='center', va='center', color='0.55',
                           fontsize=v.font_size, transform=self.ax_swath.transAxes)
        return self.fig

    def reset_axes(self):
        v = self.v
        self.ax_swath.cla()
        self.ax_cbar.cla()
        self.ax_cbar.set_axis_off()
        self.ax_swath.set_ylabel(f'LOS velocity{unit_label(v.disp_unit)}',
                                 fontsize=v.font_size + 1, labelpad=8)
        self.ax_swath.tick_params(which='both', direction='in', labelsize=v.font_size,
                                  bottom=True, top=True, left=True, right=True)
        if self.ax_topo is not None:
            self.ax_topo.cla()
            self.ax_topo.set_ylabel('Elevation [m]', fontsize=v.font_size + 1, labelpad=8)
            self.ax_topo.set_xlabel('Distance along profile [km]',
                                    fontsize=v.font_size + 1, labelpad=6)
            self.ax_topo.tick_params(which='both', direction='in', labelsize=v.font_size,
                                     bottom=True, top=True, left=True, right=True)
            self.ax_swath.tick_params(labelbottom=False)     # shared x-axis
        else:
            self.ax_swath.set_xlabel('Distance along profile [km]',
                                     fontsize=v.font_size + 1, labelpad=6)
        return

    ##---------- its artists on the map
    def clear_map_artists(self):
        for attr in ('line_obj', 'poly_obj', 'txt0', 'txt1'):
            obj = getattr(self, attr, None)
            if obj is not None:
                try:
                    obj.remove()
                except (ValueError, NotImplementedError):
                    pass
                setattr(self, attr, None)

    def draw_on_map(self):
        """The centre line, the buffer polygon and the A / A' labels."""
        v = self.v
        self.clear_map_artists()
        geom = line_geometry(v.atr, self.start_yx, self.end_yx, v.width)

        pxs, pys = v._yx2disp(geom['corners_yx'][:, 0], geom['corners_yx'][:, 1])
        self.poly_obj = patches.Polygon(np.stack([pxs, pys], axis=1), closed=True,
                                        facecolor='w', edgecolor='k', alpha=0.22,
                                        lw=1.0, zorder=8)
        v.ax_img.add_patch(self.poly_obj)

        lxs, lys = v._yx2disp([self.start_yx[0], self.end_yx[0]],
                              [self.start_yx[1], self.end_yx[1]])
        self.line_obj = v.ax_img.plot(lxs, lys, '-', color='k', lw=1.4, zorder=9)[0]

        kw = dict(fontsize=v.font_size + 1, fontweight='bold', zorder=10,
                  ha='center', va='center',
                  bbox=dict(boxstyle='round,pad=0.18', fc='w', ec='k', lw=0.7, alpha=0.85))
        self.txt0 = v.ax_img.annotate(self.label, xy=(lxs[0], lys[0]),
                                      xytext=(-11, -11), textcoords='offset points', **kw)
        self.txt1 = v.ax_img.annotate(f"{self.label}'", xy=(lxs[1], lys[1]),
                                      xytext=(11, 11), textcoords='offset points', **kw)
        return

    ##---------- its data
    def extract(self):
        """Pull the swath (and the topography swath) for this line."""
        v = self.v
        self.swath_list = []
        for i, fname in enumerate(v.file):
            atr_i = v.atr_list[i]
            if self.start_lalo is not None and 'Y_FIRST' in atr_i:
                coord_i = ut.coordinate(atr_i)
                (y0, y1), (x0, x1) = coord_i.lalo2yx(
                    [self.start_lalo[0], self.end_lalo[0]],
                    [self.start_lalo[1], self.end_lalo[1]])
                syx, eyx = [y0, x0], [y1, x1]
            else:
                syx, eyx = list(self.start_yx), list(self.end_yx)
            try:
                swath = swath_extract(v.data_list[i], atr_i, v.box_list[i], syx, eyx,
                                      v.width, mask=v.msk_list[i],
                                      keep_zero=v.keep_zero, decimate=v.decimate,
                                      nlook=v.nlook)
            except ValueError as err:
                print(f'WARNING: {os.path.basename(fname)}: {err}')
                continue
            swath['file'] = fname
            swath['start_yx'], swath['end_yx'] = syx, eyx
            self.swath_list.append(swath)

        self.swath_list = [s for s in self.swath_list if s['value'].size > 0]
        if not self.swath_list:
            print(f'WARNING: section {self.name}: no valid pixel inside the corridor.')
            self.swath_list = None
            return False

        swath = self.swath_list[0]
        if swath['clipped']:
            print(f'WARNING: section {self.name}: the corridor reaches beyond the data'
                  ' coverage (or the displayed subset) and was clipped there.')

        self.den = swath_density(swath['dist'], swath['value'], bin_x=v.bin_x,
                                 bin_y=v.bin_y, ylim=v.ylim, x_offset=v.x_offset)

        # where the centre line crosses the fault traces, in km from the start
        self.clear_crossings()
        if v.faults:
            self.crossings, self.crossing_pts = profile_crossings(
                v.faults, v.fault_atr, self.start_yx, self.end_yx, swath['length'])

        if v.topoxsect and v.topo_atr is not None:
            self.extract_topo(swath['length'])
        return True

    def extract_topo(self, main_length):
        """The same corridor, sampled from the DEM on its own grid."""
        v = self.v
        self.topo_swath = self.topo_den = self.topo_error = None
        self.topo_note = ''
        try:
            # the section on the DEM's grid (see init_topo for the three modes); the
            # end points may lie off the DEM: whatever part it covers is used
            if v.topo_mode == 'grid':
                syx, eyx, grid_atr = list(self.start_yx), list(self.end_yx), v.atr
            elif v.topo_mode == 'lookup':
                lon, lat = v.lookup.to_lonlat([self.start_yx[0], self.end_yx[0]],
                                              [self.start_yx[1], self.end_yx[1]])
                if not np.all(np.isfinite(np.r_[lon, lat])):
                    raise ValueError('the section end points have no lon/lat in the lookup')
                syx, eyx = lonlat_line_to_grid((lon[0], lat[0]), (lon[1], lat[1]),
                                               v.topo_atr, v.topo_crs)
                grid_atr = v.topo_atr
            else:
                syx, eyx = line_to_other_grid(v.atr, v.topo_atr, self.start_yx, self.end_yx,
                                              verbose=False, crs_src=v.data_crs,
                                              crs_dst=v.topo_crs, check_inside=False)
                grid_atr = v.topo_atr
            geom = line_geometry(grid_atr, syx, eyx, v.width)
            rows_c, cols_c = geom['corners_yx'][:, 0], geom['corners_yx'][:, 1]
            length, width = int(v.topo_atr['LENGTH']), int(v.topo_atr['WIDTH'])
            box = (max(0, int(np.floor(cols_c.min())) - 1),
                   max(0, int(np.floor(rows_c.min())) - 1),
                   min(width, int(np.ceil(cols_c.max())) + 2),
                   min(length, int(np.ceil(rows_c.max())) + 2))
            if box[2] <= box[0] or box[3] <= box[1]:
                raise ValueError('the profile does not overlap the DEM')

            # its own sampling, not --nlook (see topo_sampling)
            pixel_km = max(geom['x_step'], geom['y_step']) * geom['dist_scale']
            n_pixels = (box[2] - box[0]) * (box[3] - box[1])
            decimate, nlook = topo_sampling(n_pixels, pixel_km, v.bin_x)
            how = (f'every {_ordinal(decimate)} pixel' if decimate > 1 else
                   f'{nlook}x oversampled' if nlook > 1 else 'every pixel')
            self.topo_note = f'{how} of a {pixel_km * 1000:.0f} m DEM'

            dem = read_any(v.topoxsect, box=box)[0].astype(np.float64)
            dem[(dem < DEM_VALID[0]) | (dem > DEM_VALID[1])] = np.nan   # voids, no-data

            # keep_zero: 0 m is sea level, a real elevation
            sw = swath_extract(dem, grid_atr, box, syx, eyx, v.width,
                               keep_zero=True, decimate=decimate, nlook=nlook)
            if sw['value'].size == 0:
                raise ValueError('no valid DEM pixel inside the corridor')

            # the two grids may differ slightly in projected length; rescale so the
            # topography panel registers exactly with the velocity panel above it
            if sw['length'] > 0:
                sw['dist'] = sw['dist'] * (main_length / sw['length'])
                sw['length'] = main_length

            self.topo_swath = sw
            self.topo_den = swath_density(sw['dist'], sw['value'], bin_x=v.bin_x,
                                          x_offset=v.x_offset)
            v.vprint(f"   topography: {sw['value'].size:,} DEM samples ({self.topo_note})")
        except Exception as err:      # noqa: BLE001 - topography is a nicety
            self.topo_error = str(err)
            print(f'WARNING: section {self.name}: no topography profile ({err})')
        return

    ##---------- draw
    def draw(self):
        """Extract and (re)draw this section's window."""
        v = self.v
        self.create_figure()
        if not self.extract():
            return
        self.reset_axes()
        ax, den = self.ax_swath, self.den

        cmap = plt.get_cmap(v.pcmap)
        try:
            cmap = cmap.with_extremes(over='white')      # matplotlib >= 3.4
        except AttributeError:
            cmap = cmap.copy()
            cmap.set_over('white')
        counts = np.ma.masked_where(den['counts'] == 0, den['counts'])
        cmax = v.cmax
        if cmax is None:
            # auto: the 98th percentile of the OCCUPIED cells, so the colour
            # scale adapts to how densely the corridor was sampled instead of
            # leaving coarse data stuck on the darkest step of the log scale
            occupied = den['counts'][den['counts'] > 0]
            cmax = float(np.percentile(occupied, 98)) if occupied.size else 10.0
            cmax = max(cmax, v.cmin * 10)
        self.cmax_used = cmax
        pcm = ax.pcolormesh(den['x_edges'], den['y_edges'], counts.T, cmap=cmap,
                            shading='flat',
                            norm=mcolors.LogNorm(vmin=v.cmin, vmax=cmax))

        self.ax_cbar.set_axis_on()
        cbar = self.fig.colorbar(pcm, cax=self.ax_cbar, extend='max')
        cell_unit = v.disp_unit.split('/')[0] + ' ' if unit_label(v.disp_unit) else ''
        cbar.set_label(f"Point density\n[count / {den['bin_x']:g} km"
                       f" x {den['bin_y']:g} {cell_unit}cell]",
                       fontsize=max(v.font_size - 2, 6))
        cbar.ax.tick_params(labelsize=max(v.font_size - 2, 6))

        ax.plot(den['x_centers'], den['median_dw'], color='red', lw=1.3,
                label='Median (density-weighted)')
        ax.plot(den['x_centers'], den['q25_dw'], color='0.55', ls='--', lw=0.8,
                label='25th-75th percentile')
        ax.plot(den['x_centers'], den['q75_dw'], color='0.55', ls='--', lw=0.8)

        for i, sw in enumerate(self.swath_list[1:], start=1):
            den_i = swath_density(sw['dist'], sw['value'], bin_x=v.bin_x,
                                  bin_y=den['bin_y'], ylim=den['ylim'],
                                  x_offset=v.x_offset)
            ax.plot(den_i['x_centers'], den_i['median_dw'], color=f'C{i}', lw=1.2,
                    label=f"Median - {os.path.basename(sw['file'])}")

        ax.axhline(0, color='k', ls='-', lw=0.4)
        if v.x_offset != 0:
            ax.axvline(0, color='blue', ls='--', lw=1.0)
        ax.set_xlim(den['x_edges'][0], den['x_edges'][-1])
        ax.set_ylim(den['ylim'])
        ax.legend(loc='lower right', fontsize=max(v.font_size - 3, 6))

        # the A / A' end labels, echoed inside the panel corners
        lbl_kw = dict(xycoords='axes fraction', fontsize=v.font_size + 3,
                      fontweight='bold', va='top', zorder=12,
                      bbox=dict(boxstyle='round,pad=0.2', fc='w', ec='0.4',
                                lw=0.6, alpha=0.85))
        ax.annotate(self.label, xy=(0.015, 0.975), ha='left', **lbl_kw)
        ax.annotate(f"{self.label}'", xy=(0.985, 0.975), ha='right', **lbl_kw)

        swath = self.swath_list[0]
        msg = (f"swath +/-{v.width:g} km | length {swath['length']:.2f} km"
               f" | N = {swath['value'].size:,}")
        if v.nlook > 1:
            msg += f" (nlook {v.nlook})"
        if self.start_lalo is not None:
            # X/Y_FIRST are northing/easting on a projected grid, not lat/lon
            tag, dp = ('N/E', 1) if is_projected(v.atr) else ('lat/lon', 4)
            msg += (f"\n{tag}: ({self.start_lalo[0]:.{dp}f}, {self.start_lalo[1]:.{dp}f})"
                    f" --> ({self.end_lalo[0]:.{dp}f}, {self.end_lalo[1]:.{dp}f})")
        else:
            msg += (f"\ny/x: ({self.start_yx[0]}, {self.start_yx[1]})"
                    f" --> ({self.end_yx[0]}, {self.end_yx[1]})")
        ax.set_title(msg, fontsize=v.font_size, pad=8)

        # the topography panel - median line only, in black. The EXACT median
        # per bin: the density-weighted one snaps to the height cells of the
        # histogram, 50-100 m steps over a mountain range.
        if self.ax_topo is not None:
            tden = self.topo_den
            if tden is not None:
                self.ax_topo.plot(tden['x_centers'], tden['median'],
                                  color='k', lw=1.4, label='Topography (median)')
                self.ax_topo.set_xlim(den['x_edges'][0], den['x_edges'][-1])
                finite = np.isfinite(tden['median'])
                if finite.any():
                    lo = float(np.nanmin(tden['median'][finite]))
                    hi = float(np.nanmax(tden['median'][finite]))
                    pad = max((hi - lo) * 0.12, 1.0)
                    self.ax_topo.set_ylim(lo - pad, hi + pad)
            else:
                # say WHY, here, not only in the terminal
                kw = dict(ha='center', va='center', color='0.45',
                          transform=self.ax_topo.transAxes)
                self.ax_topo.text(0.5, 0.68, 'no topography for this profile',
                                  fontsize=v.font_size, **kw)
                if self.topo_error:
                    reason = textwrap.shorten(self.topo_error, width=170, placeholder=' ...')
                    self.ax_topo.text(0.5, 0.28, '\n'.join(textwrap.wrap(reason, 75)),
                                      fontsize=max(v.font_size - 4, 6), **kw)

        self.draw_fault_marks()

        v.vprint(f"section {self.name}: {swath['value'].size:,} samples, "
                 f"length {swath['length']:.2f} km, half width {swath['half_width']:.2f} km")
        if self.crossings.size:
            v.vprint(f"   crosses a fault trace at "
                     f"{', '.join(f'{x:.2f}' for x in self.crossings - v.x_offset)} km")
        if den['num_outside']:
            v.vprint(f"   {den['num_outside']:,} samples fall outside the Y range "
                     f"[{den['y_edges'][0]:g}, {den['y_edges'][-1]:g}] {v.disp_unit}")

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return

    def draw_fault_marks(self):
        """A dashed vertical line at every fault crossing, on every panel."""
        v = self.v
        if not self.crossings.size:
            return
        xs = self.crossings - v.x_offset
        ax = self.ax_swath

        # a crossing over a no-data stretch near either end may fall outside the
        # density grid: widen the (shared) x-axis rather than hide it
        lo, hi = ax.get_xlim()
        pad = 0.01 * (hi - lo)
        lo, hi = min(lo, xs.min() - pad), max(hi, xs.max() + pad)

        kw = dict(color=v.fault_color, ls='--', lw=1.2, zorder=6)
        for x in xs:
            ax.axvline(x, **kw)
            if self.ax_topo is not None:
                self.ax_topo.axvline(x, **kw)
        ax.set_xlim(lo, hi)              # the topography panel shares it

        if v.fault_label:
            self._label_crossings(xs, lo, hi)
        return

    def _label_crossings(self, xs, lo, hi):
        """Annotate the crossings without letting the labels collide.

        A label is vertical text hugging its dashed line, on the right of it
        or, failing that, on the left. Where neither side is free - crossings
        a few km apart - the line stays unlabelled: its neighbour already says
        what it is. Near the A / A' boxes in the top corners a label drops below them.
        """
        v, ax = self.v, self.ax_swath
        fsize = max(v.font_size - 1, 6)
        width_pt = max(ax.get_position().width * self.fig.get_figwidth() * 72, 1.0)
        thick = 1.5 * fsize / width_pt               # rotated text + its box, in axes fraction
        gap = 3.0 / width_pt
        corner = (1.4 * (v.font_size + 3) + 6) / width_pt + 0.015   # the A / A' boxes

        taken = []
        for x in xs:
            frac = (x - lo) / (hi - lo)
            for side in ('right', 'left'):
                if side == 'right':
                    span = (frac + gap, frac + gap + thick)
                else:
                    span = (frac - gap - thick, frac - gap)
                if span[0] < 0 or span[1] > 1:
                    continue                          # would leave the panel
                if any(span[0] < b and span[1] > a for a, b in taken):
                    continue                          # would hit another label
                taken.append(span)
                in_corner = span[0] < corner or span[1] > 1 - corner
                ax.annotate(v.fault_label, xy=(x, 0.80 if in_corner else 0.97),
                            xycoords=('data', 'axes fraction'),
                            xytext=(3 if side == 'right' else -3, 0),
                            textcoords='offset points', rotation=90, va='top',
                            ha='left' if side == 'right' else 'right',
                            fontsize=fsize, color=v.fault_color, zorder=12,
                            bbox=dict(boxstyle='round,pad=0.15', fc='w', ec='none',
                                      alpha=0.75))
                break
        return

    ##---------- output
    @property
    def outfile_base(self):
        sw = self.swath_list[0]
        y0, x0 = sw['start_yx']
        y1, x1 = sw['end_yx']
        return os.path.join(
            self.v.outdir,
            f"swath_{self.label}_Y{y0}X{x0}_Y{y1}X{x1}_w{self.v.width:g}km")

    def save(self):
        """Write this section's points, statistics, metadata and figure."""
        if not self.has_data:
            print(f'section {self.name}: nothing to save yet.')
            return []
        v = self.v
        base = self.outfile_base
        os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
        written = []

        # 1. the located points - one row per sample (the LRS / swath table)
        cols = ['dist_km', 'perp_km', 'value', 'row', 'col']
        has_lalo = any('lat' in s for s in self.swath_list)
        if has_lalo:
            cols += ['lat', 'lon']
        cols += ['file_idx']
        blocks = []
        for i, sw in enumerate(self.swath_list):
            # dist_km is from the start of the line, minus --x-offset, so the
            # .csv, the statistics and the figure share one X axis.
            arrs = [sw['dist'] - v.x_offset, sw['perp'], sw['value'], sw['row'], sw['col']]
            if has_lalo:
                nan = np.full(sw['value'].size, np.nan)
                arrs += [sw.get('lat', nan), sw.get('lon', nan)]
            arrs += [np.full(sw['value'].size, i)]
            blocks.append(np.column_stack([np.asarray(a, dtype=np.float64) for a in arrs]))
        table = np.vstack(blocks)
        rc_fmt = '%.4f' if v.nlook > 1 else '%d'
        fmt = (['%.5f', '%.5f', '%.6f', rc_fmt, rc_fmt]
               + (['%.8f', '%.8f'] if has_lalo else []) + ['%d'])
        out = f'{base}_points.csv'
        np.savetxt(out, table, fmt=fmt, delimiter=',', header=','.join(cols), comments='')
        print(f'saved {table.shape[0]:,} swath samples to: {out}')
        written.append(out)

        # 2. the binned statistics, with the topography median alongside
        den = self.den
        cols = ['dist_km', 'count', 'median', 'mean', 'std', 'q25', 'q75',
                'median_density_weighted', 'q25_density_weighted', 'q75_density_weighted']
        stat = [den['x_centers'], den['count'], den['median'], den['mean'],
                den['std'], den['q25'], den['q75'], den['median_dw'],
                den['q25_dw'], den['q75_dw']]
        if self.topo_den is not None:
            # the exact median, as drawn; NaN where the DEM had nothing
            topo_on_grid = np.interp(den['x_centers'], self.topo_den['x_centers'],
                                     self.topo_den['median'],
                                     left=np.nan, right=np.nan)
            stat.append(topo_on_grid)
            cols.append('elevation_median_m')
        out = f'{base}_stats.csv'
        np.savetxt(out, np.column_stack(stat), fmt='%.6f', delimiter=',',
                   header=','.join(cols), comments='')
        print(f'saved the binned statistics to: {out}')
        written.append(out)

        # 3. everything needed to reproduce the profile
        sw = self.swath_list[0]
        out = f'{base}_info.txt'
        with open(out, 'w') as fid:
            fid.write(f'section        : {self.name}\n')
            fid.write(f'file           : {[os.path.abspath(f) for f in v.file]}\n')
            fid.write(f'dataset        : {v.dset}\n')
            fid.write(f'unit           : {v.disp_unit}\n')
            if v.data_crs is not None:
                fid.write(f'CRS            : {crs_name(v.data_crs)}\n')
            fid.write(f"start y/x      : {sw['start_yx']}\n")
            fid.write(f"end   y/x      : {sw['end_yx']}\n")
            if self.start_lalo:
                fid.write(f'start lat/lon  : {self.start_lalo}\n')
                fid.write(f'end   lat/lon  : {self.end_lalo}\n')
            fid.write(f'swath halfwidth: {v.width} km  (flat-capped buffer)\n')
            fid.write(f"profile length : {sw['length']:.4f} km\n")
            fid.write(f"num of samples : {sw['value'].size}\n")
            fid.write(f'nlook          : {v.nlook}  '
                      f'({"native pixels" if v.nlook == 1 else f"{v.nlook}x oversampled"})\n')
            fid.write(f"bin size       : {den['bin_x']:g} km x {den['bin_y']:.6g} {v.disp_unit}\n")
            fid.write(f'x offset       : {v.x_offset} km  (already applied to dist_km)\n')
            if v.topoxsect:
                fid.write(f'DEM            : {os.path.abspath(v.topoxsect)}'
                          f'  ({crs_name(v.topo_crs)})\n')
                if v.dem_download and v.topoxsect == v.dem_download:
                    fid.write(f'DEM source     : {ATTRIBUTION}\n')
                if self.topo_error:
                    fid.write(f'DEM sampling   : none - {self.topo_error}\n')
                else:
                    fid.write(f'DEM sampling   : {self.topo_note}; elevation_median_m '
                              'is the exact median per bin\n')
            if v.fault_file:
                fid.write(f'fault file     : {[os.path.abspath(f) for f in v.fault_file]}\n')
                radar = v.fault_atr is PIXEL_ATR
                tag, dp = ('N/E', 1) if is_projected(v.atr) else ('lat/lon', 6)
                for dist, (px, py) in zip(self.crossings, self.crossing_pts):
                    where = f'{tag} {py:.{dp}f}, {px:.{dp}f}'
                    if radar:              # (px, py) are (col, row) on the radar grid
                        lon, lat = v.lookup.to_lonlat([py], [px])
                        where = f'y/x {py:.1f}, {px:.1f}, lat/lon {lat[0]:.6f}, {lon[0]:.6f}'
                    fid.write(f'fault crossing : {dist - v.x_offset:.4f} km  (dist_km), '
                              f'{where}\n')
                if not self.crossings.size:
                    fid.write('fault crossing : none - the centre line crosses no trace\n')
            fid.write('\ndist_km = distance along the profile from the start point,'
                      ' minus the x offset\n')
            fid.write('perp_km = distance from the centre line,'
                      " + = left of the start->end direction\n")
        print(f'saved the profile metadata to: {out}')
        written.append(out)

        # 4. the figure
        out = f'{base}{v.fig_ext}'
        self.fig.savefig(out, bbox_inches='tight', dpi=v.fig_dpi)
        print(f'saved section {self.name} to: {out}')
        written.append(out)
        return written


##############################################################################
class swathProfileViewer:
    """Interactive swath (buffered) profile viewer.

    Example:
        from viewinsar.cli import cmd_line_parse, get_view_cmd
        from viewinsar.viewer import swathProfileViewer

        inps = cmd_line_parse('geo_velocity.h5 -w 15'.split())
        obj = swathProfileViewer(inps, get_view_cmd(inps))
        obj.open()
        obj.plot()
    """

    COMMAND = 'swath'                # the sub-command, for the echoed command line

    def __init__(self, inps, view_cmd):
        self.fig_map = self.ax_img = None
        self.fig_ts = self.ax_ts = None
        self.sections = []
        self.pts_idx = 0
        self.ts_obj = None
        self.topo_atr = None
        self.data_crs = self.topo_crs = None     # EPSG int or WKT, read from the files
        self.dem_file = self.dem_zoom = None
        self.dem_download = None                 # path of a downloaded (Terrain Tiles) DEM
        self.lookup_file = self.lookup = None    # radar-coded data: geometryRadar.h5
        self.topo_mode = None                    # 'geo', 'grid' (radar, same grid), 'lookup'
        self.fault_atr = None                    # the "map" the fault traces are in
        self.faults = []             # fault traces, in the data's map coordinates
        self.fault_obj = None
        self.fault_file = None
        self.fault_label, self.fault_color, self.fault_lw = 'Fault line', 'k', 1.2

        self.view_cmd = view_cmd
        # view.py owns names like inps.width / inps.length (the file dimensions),
        # so anything that came from OUR parser must win over its namesakes.
        self._own_keys = set(inps.__dict__.keys()) - {'print_msg', 'atr'}
        for key, value in inps.__dict__.items():
            setattr(self, key, value)

    def vprint(self, *args, **kwargs):
        if self.print_msg:
            print(*args, **kwargs)

    @property
    def section(self):
        """The section currently being drawn."""
        return self.sections[-1] if self.sections else None

    ##---------- setup
    def open(self):
        if getattr(self, 'argv', None) is not None:
            print(f'viewinsar {self.COMMAND} ' + ' '.join(self.argv))
        if not self.disp_fig:
            plt.switch_backend('Agg')

        # the map: let view.py do all of the work
        self.vprint(f'run: {self.view_cmd}')
        self.data_img, self.atr, view_inps = view.prep_slice(self.view_cmd)
        for key, value in view_inps.__dict__.items():
            if key not in self._own_keys:
                setattr(self, key, value)

        # read the data of every input file, within the displayed window
        self.data_list, self.atr_list, self.box_list, self.msk_list = [], [], [], []
        kept = []
        base_geo = 'Y_FIRST' in self.atr
        for fname in self.file:
            atr_i = readfile.read_attribute(fname)
            if ('Y_FIRST' in atr_i) != base_geo:
                print(f'WARNING: {os.path.basename(fname)} is in '
                      f'{"geo" if "Y_FIRST" in atr_i else "radar"} coordinates while '
                      f'{os.path.basename(self.file[0])} is not --> skipped.')
                continue
            kept.append(fname)

            coord_i = ut.coordinate(atr_i)
            if self.geo_box and 'Y_FIRST' in atr_i:
                box_i = coord_i.check_box_within_data_coverage(
                    coord_i.box_geo2pixel(self.geo_box), print_msg=False)
            else:
                box_i = tuple(self.pix_box)

            data_i = readfile.read(fname, datasetName=self.dset, box=box_i,
                                   print_msg=False)[0]
            data_i = pp.scale_data2disp_unit(data_i, metadata=atr_i,
                                             disp_unit=self.disp_unit)[0]
            msk_i = pp.read_mask(fname, mask_file=self.mask_file, datasetName=self.dset,
                                 box=box_i, print_msg=False)[0]

            self.data_list.append(data_i)
            self.atr_list.append(atr_i)
            self.box_list.append(box_i)
            self.msk_list.append(msk_i)
            self.vprint(f'read {os.path.basename(fname)} within box {box_i}'
                        f' -> {data_i.shape} [{self.disp_unit}]')

        if not kept:
            raise ValueError('no usable input file left.')
        self.file = kept

        # the data's CRS, from the file itself where it carries one (GeoTIFF):
        # MintPy's attributes alone can misreport a projected GeoTIFF as EPSG:4326
        self.data_crs = crs_of(self.atr, self.file[0])

        # radar-coded data: MintPy's geometry file is the way to lon/lat
        if 'Y_FIRST' not in self.atr and (self.fault_file or self.topoxsect or self.dem_file
                                          or self.lookup_file):
            self.init_lookup()

        # -d / --topoxsect without a file: one download, for the map extent
        if is_auto(self.dem_file) or is_auto(self.topoxsect):
            self.init_terrain()

        if self.nlook > 1:
            self.vprint(f'--nlook {self.nlook}: the corridor is sampled {self.nlook}x '
                        f'finer than the native pixel ({self.nlook**2}x more samples, '
                        'NaN-aware bilinear)')

        if self.topoxsect:
            self.init_topo()
        if self.fault_file:
            self.init_faults()
        if self.timeseries:
            self.init_timeseries()

        # figure sizes: the map and the sections live in separate windows
        length = self.pix_box[3] - self.pix_box[1]
        width = self.pix_box[2] - self.pix_box[0]
        if not self.fig_size_map:
            self.fig_size_map = pp.auto_figure_size((length, width), disp_cbar=True)
        if not self.fig_size:
            self.fig_size = [7.6, 4.8]
        return

    def init_lookup(self):
        """Find and read the radar-coded geometry file (--lookup, or MintPy's places)."""
        try:
            path = find_lookup(self.file[0], self.atr, given=self.lookup_file,
                               verbose=self.print_msg)
            self.lookup_file = path
            self.lookup = RadarLookup(path) if path else None
        except Exception as err:      # noqa: BLE001
            print(f'WARNING: {err}')
            self.lookup_file, self.lookup = None, None
        return

    def init_terrain(self):
        """A DEM for -d and/or --topoxsect given without a file.

        Geocoded data: the Terrain Tiles over the map, downloaded. Radar-coded
        data: the `height` of the geometry file - the DEM of the processing, on
        the very grid of the data - or, if it has none, the Terrain Tiles
        resampled onto the radar grid through its latitude / longitude.
        """
        path = None
        try:
            if 'Y_FIRST' not in self.atr:
                if self.lookup is None:
                    raise ValueError(NO_LOOKUP)
                if has_height(self.lookup_file):
                    path = self.lookup_file
                    print(f'DEM: the height of {self.lookup_file} '
                          '(the DEM of the processing, on the radar grid)')
                else:
                    path = terrain_dem_radar(self.lookup, self.atr, zoom=self.dem_zoom,
                                             verbose=self.print_msg)
                    self.dem_download = path
                    print(f'DEM: {ATTRIBUTION}')
            else:
                if not self.geo_box:
                    raise ValueError('the map extent is unknown')
                if self.data_crs is None:
                    raise ValueError("the data's coordinate reference system is unknown")
                path = terrain_dem(self.geo_box, self.data_crs, is_projected(self.atr),
                                   zoom=self.dem_zoom, verbose=self.print_msg)
                self.dem_download = path
                print(f'DEM: {ATTRIBUTION}')
        except Exception as err:      # noqa: BLE001 - the DEM is a nicety
            print(f'WARNING: no DEM for the map: {err}')
            print('         --> continuing without it. A DEM file works too: -d dem.tif')

        if is_auto(self.dem_file):
            self.dem_file = path
            # what view.py does for a DEM it is given: let the relief show through
            if path and self.disp_dem_shade and '--alpha' not in self.view_extra:
                self.transparency = 0.8
        if is_auto(self.topoxsect):
            self.topoxsect = path
        return

    def init_topo(self):
        """Validate the DEM for --topoxsect, and work out how sections reach it.

        topo_mode: 'geo'    - data and DEM geocoded; sections re-projected if the
                              CRSs differ
                   'grid'   - radar-coded DEM on the data's own radar grid (the
                              height of geometryRadar.h5): sections as they are
                   'lookup' - radar-coded data, geocoded DEM: the section end
                              points go to lon/lat through geometryRadar.h5
        """
        try:
            self.topo_atr = read_attr_any(self.topoxsect)
            data_geo, dem_geo = 'Y_FIRST' in self.atr, 'Y_FIRST' in self.topo_atr
            size = f"{self.topo_atr['LENGTH']} x {self.topo_atr['WIDTH']}"
            if not dem_geo:
                if data_geo:
                    raise ValueError('the DEM is in radar coordinates, the data are geocoded')
                if (self.topo_atr['LENGTH'], self.topo_atr['WIDTH']) != (self.atr['LENGTH'],
                                                                          self.atr['WIDTH']):
                    raise ValueError(f"a radar-coded DEM must be on the data's grid: it is "
                                     f"{size}, the data {self.atr['LENGTH']} x "
                                     f"{self.atr['WIDTH']}")
                self.topo_mode = 'grid'
                msg = f'topography: {self.topoxsect} ({size}, the radar grid of the data)'
            elif not data_geo:
                if self.lookup is None:
                    raise ValueError(NO_LOOKUP)
                self.topo_mode = 'lookup'
                self.topo_crs = crs_of(self.topo_atr, self.topoxsect)
                msg = (f'topography: {self.topoxsect} ({size}, {crs_name(self.topo_crs)}); '
                       f'sections reach it through {os.path.basename(self.lookup_file)}')
            else:
                self.topo_mode = 'geo'
                self.topo_crs = crs_of(self.topo_atr, self.topoxsect)
                msg = f'topography: {self.topoxsect} ({size}, {crs_name(self.topo_crs)})'
                if self.topo_crs is None or self.data_crs is None:
                    msg += (' - WARNING: the CRS of the '
                            f'{"DEM" if self.topo_crs is None else "data"} is unknown, '
                            'taken to be the same as the other')
                elif not same_crs(self.topo_crs, self.data_crs):
                    msg += f'; sections are re-projected from {crs_name(self.data_crs)}'
            self.vprint(msg)
        except Exception as err:      # noqa: BLE001
            print(f'WARNING: cannot use the DEM for --topoxsect: {err}')
            print('         --> continuing without the topography panel.')
            self.topoxsect, self.topo_atr = None, None
        return

    def init_faults(self):
        """Read the --fault traces and bring them onto the data's grid.

        Geocoded data: into the data's CRS. Radar-coded data: to radar pixels
        through geometryRadar.h5 - fault_atr is then radar.PIXEL_ATR, the "map"
        whose coordinates are (col, row).
        """
        try:
            if 'Y_FIRST' in self.atr:
                self.faults = load_faults(self.fault_file, self.atr, verbose=self.print_msg,
                                          data_crs=self.data_crs)
                self.fault_atr = self.atr
            else:
                if self.lookup is None:
                    raise ValueError(NO_LOOKUP)
                self.faults = load_faults_radar(self.fault_file, self.lookup,
                                                verbose=self.print_msg)
                self.fault_atr = PIXEL_ATR
        except Exception as err:      # noqa: BLE001 - faults are a nicety
            print(f'WARNING: cannot use --fault: {err}')
            print('         --> continuing without fault lines.')
            self.faults = []
            return
        if not self.faults:
            print('WARNING: --fault: no fault trace falls within the data extent. '
                  'Does the file cover this area, with the right CRS (.prj)?')
        return

    def init_timeseries(self):
        """Set up tsview.py to serve the time-series window."""
        try:
            from mintpy import tsview
            from mintpy.cli.tsview import cmd_line_parse as ts_cmd_line_parse
        except ImportError as err:
            print(f'WARNING: cannot import tsview.py ({err}) --> no time-series window.')
            self.timeseries = None
            return

        ts_args = list(self.timeseries)
        if not self.print_msg:
            ts_args += ['--noverbose']
        # keep the time-series cube within the displayed area, to save memory
        if self.subset_lat and self.subset_lon:
            ts_args += ['--sub-lat', *[str(i) for i in self.subset_lat],
                        '--sub-lon', *[str(i) for i in self.subset_lon]]
        elif self.subset_y and self.subset_x:
            ts_args += ['--sub-y', *[str(i) for i in self.subset_y],
                        '--sub-x', *[str(i) for i in self.subset_x]]

        try:
            ts_inps = ts_cmd_line_parse(ts_args)
            self.ts_obj = tsview.timeseriesViewer(ts_inps)
            self.ts_obj.open()
            self.ts_obj.ts_data, self.ts_obj.mask = tsview.read_timeseries_data(self.ts_obj)[0:2]
            self.vprint(f'time-series window ready: {self.timeseries}')
        except Exception as err:      # noqa: BLE001 - never kill the main tool
            print(f'WARNING: failed to initiate the time-series window: {err}')
            self.timeseries, self.ts_obj = None, None
        return

    ##---------- main
    def make_map(self):
        """The map window: view.py's figure, the fault traces, and the mouse / keys."""
        subplot_kw = dict(projection=self.map_proj_obj) if self.map_proj_obj is not None else {}
        self.fig_map = _new_figure('Map', self.fig_size_map)
        self.ax_img = self.fig_map.add_subplot(1, 1, 1, **subplot_kw)
        self.ax_img = view.plot_slice(self.ax_img, self.data_img, self.atr, self)[0]
        self.draw_faults_on_map()
        self.fig_map.tight_layout()
        self.cid = self.fig_map.canvas.mpl_connect('button_release_event', self.select_point)
        self.fig_map.canvas.mpl_connect('key_press_event', self.on_key)
        return

    def event_yx(self, event):
        """(row, col) of the full file under a mouse event, and lat/lon when shown."""
        if 'Y_FIRST' in self.atr and self.fig_coord == 'geo':
            return list(self.coord.lalo2yx(event.ydata, event.xdata)), [event.ydata, event.xdata]
        return [int(event.ydata + 0.5), int(event.xdata + 0.5)], None

    def plot(self):
        # window 1 - the map
        self.make_map()

        # window 2 - section A
        self.sections.append(Section(self, 0))
        self.section.create_figure()

        # window 3 - the point time series
        if self.ts_obj is not None:
            self.fig_ts = _new_figure('Point time series', [7.6, 3.8])
            self.ax_ts = self.fig_ts.add_subplot(1, 1, 1)
            self.ts_obj.ax_pts = self.ax_ts
            self.ts_obj.fig_pts = self.fig_ts
            self.ax_ts.text(0.5, 0.5, 'right-click a pixel on the map',
                            ha='center', va='center', color='0.55',
                            fontsize=self.font_size, transform=self.ax_ts.transAxes)
            self.ax_ts.set_xticks([])
            self.ax_ts.set_yticks([])
            self.fig_ts.canvas.mpl_connect('key_press_event', self.on_key)

        # an initial line given on the command line
        if self.start_yx and self.end_yx:
            sec = self.section
            sec.start_yx, sec.end_yx = list(self.start_yx), list(self.end_yx)
            sec.start_lalo, sec.end_lalo = self.start_lalo, self.end_lalo
            sec.draw_on_map()
            sec.draw()
        if self.ts_obj is not None and self.yx:
            self.draw_timeseries(self.yx)

        if self.save_fig:
            self.save_outputs()
        if self.disp_fig:
            print(self.help_msg())
            plt.show()
        return

    def help_msg(self):
        msg = '\n' + '-' * 74
        msg += '\n  left click x2 : draw the centre line of the current section on the map'
        if self.ts_obj is not None:
            msg += '\n  right click   : plot the time series of that pixel'
        msg += '\n  key "="  or  "+" : start ANOTHER section in a new window'
        msg += "\n                     (the current line stays on the map as A-A', B-B', ...)"
        msg += '\n  key "n"       : redraw the CURRENT section from scratch'
        msg += '\n  key "e"       : export the map and every section'
        msg += f'  (figures as {self.fig_ext})'
        msg += '\n' + '-' * 74
        return msg

    ##---------- events
    def select_point(self, event):
        """Pick the profile end points (left) / the time-series pixel (right)."""
        if event.inaxes is not self.ax_img or event.xdata is None:
            return
        yx, lalo = self.event_yx(event)

        if event.button == 3:
            self.draw_timeseries(yx)
            return
        if event.button != 1:
            return

        sec = self.section
        if self.pts_idx == 0:
            sec.start_lalo, sec.start_yx = lalo, yx
            sec.end_lalo, sec.end_yx = None, None
            self.vprint(f'section {sec.name} start: y/x = {yx}'
                        + (f', lat/lon = {lalo}' if lalo else ''))
        else:
            sec.end_lalo, sec.end_yx = lalo, yx
            self.vprint(f'section {sec.name} end  : y/x = {yx}'
                        + (f', lat/lon = {lalo}' if lalo else ''))

        self.pts_idx = (self.pts_idx + 1) % 2
        if self.pts_idx == 0:
            sec.draw_on_map()
            self.fig_map.canvas.draw_idle()
            sec.draw()
        return

    def on_key(self, event):
        if event.key == 'e':
            self.save_outputs()

        elif event.key in ('=', '+'):
            sec = self.section
            if sec is not None and not sec.has_data:
                print(f'section {sec.name} has no profile yet - '
                      'draw it first (2 left clicks on the map).')
                return
            new = Section(self, len(self.sections))
            self.sections.append(new)
            self.pts_idx = 0
            new.create_figure()
            print(f'started section {new.name} - click 2 points on the map.'
                  f'  (sections so far: '
                  f"{', '.join(s.name for s in self.sections[:-1])})")
            if self.disp_fig:
                new.fig.show()

        elif event.key == 'n':
            sec = self.section
            if sec is None:
                return
            sec.clear_map_artists()
            sec.start_yx = sec.end_yx = sec.start_lalo = sec.end_lalo = None
            sec.swath_list = sec.den = sec.topo_swath = sec.topo_den = None
            sec.clear_crossings()
            self.pts_idx = 0
            sec.reset_axes()
            sec.fig.canvas.draw_idle()
            self.fig_map.canvas.draw_idle()
            print(f'section {sec.name} cleared - click 2 new points on the map.')
        return

    ##---------- drawing helpers
    def draw_faults_on_map(self):
        """The --fault traces on the map: one line, a light halo, under the sections."""
        if not self.faults:
            return
        ax = self.ax_img
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        fig_coord = 'geo' if ('Y_FIRST' in self.atr and self.fig_coord == 'geo') else 'yx'
        xs, ys = to_display(self.faults, self.fault_atr, fig_coord)
        halo = [patheffects.Stroke(linewidth=self.fault_lw + 1.8, foreground='w', alpha=0.6),
                patheffects.Normal()]
        self.fault_obj = ax.plot(xs, ys, '-', color=self.fault_color, lw=self.fault_lw,
                                 zorder=7, path_effects=halo)[0]
        # traces run past the data: keep the map extent view.py chose
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        return

    def _yx2disp(self, ys, xs):
        """(row, col) of the full file -> the coordinates of the map axes."""
        if 'Y_FIRST' in self.atr and self.fig_coord == 'geo':
            lat = float(self.atr['Y_FIRST']) + (np.asarray(ys, dtype=float) + 0.5) * float(self.atr['Y_STEP'])
            lon = float(self.atr['X_FIRST']) + (np.asarray(xs, dtype=float) + 0.5) * float(self.atr['X_STEP'])
            return lon, lat
        return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)

    def draw_timeseries(self, yx):
        """The point time series, via tsview.py, in its own window."""
        if self.ts_obj is None:
            return
        x0, y0, x1, y1 = self.ts_obj.pix_box
        if not (x0 <= yx[1] < x1 and y0 <= yx[0] < y1):
            print(f'pixel {yx} is outside of the time-series coverage {self.ts_obj.pix_box}')
            return
        try:
            self.ts_obj.plot_point_timeseries(tuple(yx))
        except Exception as err:      # noqa: BLE001
            print(f'WARNING: cannot plot the time series at {yx}: {err}')
        return

    def close(self):
        """Close every window this viewer opened (headless runs accumulate them)."""
        for fig in [self.fig_map, self.fig_ts] + [s.fig for s in self.sections]:
            if fig is not None:
                try:
                    plt.close(fig)
                except Exception:        # noqa: BLE001
                    pass
        return

    ##---------- output
    def save_outputs(self):
        """Write the map (with every section line) and each section."""
        drawn = [s for s in self.sections if s.has_data]
        if not drawn:
            print('nothing to save yet - draw a section first (2 left clicks).')
            return []

        os.makedirs(self.outdir or '.', exist_ok=True)
        written = []
        for sec in drawn:
            written += sec.save()

        out = self.outfile[0] if self.outfile else os.path.join(
            self.outdir, f'swath_map_{"".join(s.label for s in drawn)}{self.fig_ext}')
        self.fig_map.savefig(out, bbox_inches='tight', dpi=self.fig_dpi)
        print(f'saved the map with {len(drawn)} section line(s) to: {out}')
        written.append(out)
        return written

"""The interactive swath profile viewer.

Panel 1 is drawn by MintPy's view.py, the interaction follows
plot_transection.py, and the optional 3rd panel reuses tsview.py.
"""

import os

import numpy as np
from matplotlib import colors as mcolors
from matplotlib import patches
from matplotlib import pyplot as plt
from mintpy import view
from mintpy.utils import plot as pp
from mintpy.utils import readfile
from mintpy.utils import utils as ut

from .swath import line_geometry, swath_density, swath_extract

__all__ = ['swathProfileViewer']


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

    def __init__(self, inps, view_cmd):
        self.figname = 'Swath Profile'
        self.fig = self.ax_img = self.ax_swath = self.ax_cbar = self.ax_ts = None
        self.line_ann = self.line_obj = self.poly_obj = None
        self.pts_idx = 0
        self.swath_list = None
        self.den = None
        self.ts_obj = None
        self.outfile_base = 'swath_profile'

        self.view_cmd = view_cmd
        # view.py owns names like inps.width / inps.length (the file dimensions),
        # so anything that came from OUR parser must win over its namesakes.
        self._own_keys = set(inps.__dict__.keys()) - {'print_msg', 'atr'}
        for key, value in inps.__dict__.items():
            setattr(self, key, value)

    ##---------- setup
    def open(self):
        global vprint
        vprint = print if self.print_msg else (lambda *a, **k: None)
        if getattr(self, 'argv', None) is not None:
            print('viewinsar swath ' + ' '.join(self.argv))
        if not self.disp_fig:
            plt.switch_backend('Agg')

        # panel 1: let view.py do all of the map work
        vprint(f'run: {self.view_cmd}')
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

            data_i = readfile.read(fname, datasetName=self.dset, box=box_i, print_msg=False)[0]
            data_i = pp.scale_data2disp_unit(data_i, metadata=atr_i, disp_unit=self.disp_unit)[0]
            msk_i = pp.read_mask(fname, mask_file=self.mask_file, datasetName=self.dset,
                                 box=box_i, print_msg=False)[0]

            self.data_list.append(data_i)
            self.atr_list.append(atr_i)
            self.box_list.append(box_i)
            self.msk_list.append(msk_i)
            vprint(f'read {os.path.basename(fname)} within box {box_i}'
                   f' -> {data_i.shape} [{self.disp_unit}]')

        if not kept:
            raise ValueError('no usable input file left.')
        self.file = kept

        # optional 3rd panel
        if self.timeseries:
            self.init_timeseries()

        # figure size
        if not self.fig_size:
            length = self.pix_box[3] - self.pix_box[1]
            width = self.pix_box[2] - self.pix_box[0]
            fsize = pp.auto_figure_size((length, width), disp_cbar=True)
            self.fig_size = [fsize[0] + 1.35 * fsize[1], max(fsize[1], 4.6)]
        return

    def init_timeseries(self):
        """Set up tsview.py to serve the 3rd panel."""
        try:
            from mintpy import tsview
            from mintpy.cli.tsview import cmd_line_parse as ts_cmd_line_parse
        except ImportError as err:
            print(f'WARNING: cannot import tsview.py ({err}) --> no time-series panel.')
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
            vprint(f'time-series panel ready: {self.timeseries}')
        except Exception as err:      # noqa: BLE001 - never let this kill the main tool
            print(f'WARNING: failed to initiate the time-series panel: {err}')
            print('         --> continuing without it.')
            self.timeseries, self.ts_obj = None, None
        return

    ##---------- main
    def plot(self):
        self.fig = plt.figure(self.figname, figsize=self.fig_size)
        subplot_kw = dict(projection=self.map_proj_obj) if self.map_proj_obj is not None else {}

        # view.py hangs the map colorbar off the right side of ax_img and the
        # density colorbar lives in its own fixed slot, hence the wide wspace
        # and the right margin reserved here.
        gs_kw = dict(width_ratios=[1, 1.1], left=0.055, right=0.855, bottom=0.13, top=0.90)
        if self.ts_obj is not None:
            gs = self.fig.add_gridspec(2, 2, hspace=0.55, wspace=0.42, **gs_kw)
            self.ax_img = self.fig.add_subplot(gs[:, 0], **subplot_kw)
            self.ax_swath = self.fig.add_subplot(gs[0, 1])
            self.ax_ts = self.fig.add_subplot(gs[1, 1])
            self.ts_obj.ax_pts = self.ax_ts
            self.ts_obj.fig_pts = self.fig
            self.ax_ts.text(0.5, 0.5, 'right-click a pixel on the map\nto plot its time series',
                            ha='center', va='center', color='0.5', fontsize=self.font_size,
                            transform=self.ax_ts.transAxes)
            self.ax_ts.set_xticks([])
            self.ax_ts.set_yticks([])
        else:
            gs = self.fig.add_gridspec(1, 2, wspace=0.42, **gs_kw)
            self.ax_img = self.fig.add_subplot(gs[0, 0], **subplot_kw)
            self.ax_swath = self.fig.add_subplot(gs[0, 1])

        # a permanent slot for the density colorbar, so redraws do not stack them
        pos = self.ax_swath.get_position()
        self.ax_cbar = self.fig.add_axes([pos.x1 + 0.012, pos.y0, 0.011, pos.height])
        self.ax_cbar.set_axis_off()

        # panel 1 - the map, drawn by view.py
        self.ax_img = view.plot_slice(self.ax_img, self.data_img, self.atr, self)[0]

        # panel 2 - the swath profile
        self.reset_swath_axes()
        if self.start_yx and self.end_yx:
            self.draw_line(self.start_yx, self.end_yx)
            self.draw_swath(self.start_yx, self.end_yx, self.start_lalo, self.end_lalo)
        if self.ts_obj is not None and self.yx:
            self.draw_timeseries(self.yx)

        self.cid = self.fig.canvas.mpl_connect('button_release_event', self.select_point)
        self.kid = self.fig.canvas.mpl_connect('key_press_event', self.on_key)

        if self.save_fig:
            self.save_outputs()
        if self.disp_fig:
            print(self.help_msg())
            plt.show()
        return

    def help_msg(self):
        msg = '\n' + '-' * 74
        msg += '\n  left click x2 : draw the swath centre line on the map'
        if self.ts_obj is not None:
            msg += '\n  right click   : plot the time series of that pixel'
        msg += '\n  key "e"       : export the swath points + stats (.csv) and the figure'
        msg += '\n  key "n"       : clear the line and start a new one'
        msg += '\n' + '-' * 74
        return msg

    ##---------- events
    def select_point(self, event):
        """Pick the profile end points (left) / the time-series pixel (right)."""
        if event.inaxes is not self.ax_img or event.xdata is None:
            return
        if 'Y_FIRST' in self.atr and self.fig_coord == 'geo':
            lalo = [event.ydata, event.xdata]
            yx = list(self.coord.lalo2yx(event.ydata, event.xdata))
        else:
            lalo = None
            yx = [int(event.ydata + 0.5), int(event.xdata + 0.5)]

        if event.button == 3:
            self.draw_timeseries(yx)
            return
        if event.button != 1:
            return

        if self.pts_idx == 0:
            self.start_lalo, self.start_yx = lalo, yx
            self.end_lalo, self.end_yx = None, None
            vprint(f'start point: y/x = {yx}' + (f', lat/lon = {lalo}' if lalo else ''))
        else:
            self.end_lalo, self.end_yx = lalo, yx
            vprint(f'end   point: y/x = {yx}' + (f', lat/lon = {lalo}' if lalo else ''))

        self.pts_idx = (self.pts_idx + 1) % 2
        if self.pts_idx == 0:
            self.draw_line(self.start_yx, self.end_yx)
            self.draw_swath(self.start_yx, self.end_yx, self.start_lalo, self.end_lalo)
        return

    def on_key(self, event):
        if event.key == 'e':
            if self.swath_list is None:
                print('no swath profile yet - draw a line first (2 left clicks).')
            else:
                self.save_outputs()
        elif event.key == 'n':
            self.pts_idx = 0
            self.start_yx = self.end_yx = self.start_lalo = self.end_lalo = None
            self.swath_list = None
            self.clear_line()
            self.reset_swath_axes()
            self.fig.canvas.draw_idle()
            print('cleared - click 2 new points on the map.')
        return

    ##---------- drawing
    def clear_line(self):
        for attr in ('line_ann', 'line_obj', 'poly_obj'):
            obj = getattr(self, attr, None)
            if obj is not None:
                try:
                    obj.remove()
                except (ValueError, NotImplementedError):
                    pass
                setattr(self, attr, None)

    def _yx2disp(self, ys, xs):
        """(row, col) of the full file -> the coordinates of the map axes."""
        if 'Y_FIRST' in self.atr and self.fig_coord == 'geo':
            lat = float(self.atr['Y_FIRST']) + (np.asarray(ys, dtype=float) + 0.5) * float(self.atr['Y_STEP'])
            lon = float(self.atr['X_FIRST']) + (np.asarray(xs, dtype=float) + 0.5) * float(self.atr['X_STEP'])
            return lon, lat
        return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)

    def draw_line(self, start_yx, end_yx):
        """Draw the centre line + the buffer polygon on the map."""
        self.clear_line()
        geom = line_geometry(self.atr, start_yx, end_yx, self.width)

        pxs, pys = self._yx2disp(geom['corners_yx'][:, 0], geom['corners_yx'][:, 1])
        self.poly_obj = patches.Polygon(np.stack([pxs, pys], axis=1), closed=True,
                                        facecolor='w', edgecolor='k', alpha=0.25,
                                        lw=1.2, zorder=8)
        self.ax_img.add_patch(self.poly_obj)

        lxs, lys = self._yx2disp([start_yx[0], end_yx[0]], [start_yx[1], end_yx[1]])
        self.line_obj = self.ax_img.plot(lxs, lys, 'k--', lw=1.2, zorder=9)[0]
        try:
            self.line_ann = pp.add_arrow(self.line_obj, position=lxs[1])
        except Exception:      # noqa: BLE001 - the arrow is cosmetic only
            self.line_ann = self.ax_img.annotate(
                '', xytext=(lxs[0], lys[0]), xy=(lxs[1], lys[1]), zorder=9,
                arrowprops=dict(arrowstyle='-|>', color='k', lw=1.2, shrinkA=0, shrinkB=0))

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return

    def reset_swath_axes(self):
        self.ax_swath.cla()
        self.ax_cbar.cla()
        self.ax_cbar.set_axis_off()
        self.ax_swath.set_xlabel('Distance along profile [km]',
                                 fontsize=self.font_size + 2, labelpad=6)
        self.ax_swath.set_ylabel(f'LOS velocity [{self.disp_unit}]',
                                 fontsize=self.font_size + 2, labelpad=10)
        self.ax_swath.tick_params(which='both', direction='in', labelsize=self.font_size,
                                  bottom=True, top=True, left=True, right=True)
        return

    def draw_swath(self, start_yx, end_yx, start_lalo=None, end_lalo=None):
        """Extract the swath and draw the point-density profile."""
        self.swath_list = []
        for i, fname in enumerate(self.file):
            atr_i = self.atr_list[i]
            if start_lalo is not None and 'Y_FIRST' in atr_i:
                # go through lat/lon, so files on different grids still line up
                coord_i = ut.coordinate(atr_i)
                (y0, y1), (x0, x1) = coord_i.lalo2yx([start_lalo[0], end_lalo[0]],
                                                     [start_lalo[1], end_lalo[1]])
                syx, eyx = [y0, x0], [y1, x1]
            else:
                syx, eyx = list(start_yx), list(end_yx)
            try:
                swath = swath_extract(self.data_list[i], atr_i, self.box_list[i],
                                      syx, eyx, self.width, mask=self.msk_list[i],
                                      keep_zero=self.keep_zero, decimate=self.decimate)
            except ValueError as err:
                print(f'WARNING: {os.path.basename(fname)}: {err}')
                continue
            swath['file'] = fname
            swath['start_yx'], swath['end_yx'] = syx, eyx
            self.swath_list.append(swath)

        self.swath_list = [s for s in self.swath_list if s['value'].size > 0]
        if not self.swath_list:
            print('WARNING: no valid pixel found inside the swath corridor.')
            self.swath_list = None
            return

        swath = self.swath_list[0]
        if swath['clipped']:
            print('WARNING: the swath corridor reaches beyond the data coverage (or beyond'
                  ' the displayed subset) --> it was clipped there. Shorten the line, reduce'
                  ' --width, or widen --sub-lat/--sub-lon / --sub-y/--sub-x if unintended.')

        self.den = swath_density(swath['dist'], swath['value'], bin_x=self.bin_x,
                                 bin_y=self.bin_y, ylim=self.ylim, x_offset=self.x_offset)
        den = self.den

        # ---- plot
        self.reset_swath_axes()
        ax = self.ax_swath
        cmap = plt.get_cmap(self.pcmap)
        try:
            cmap = cmap.with_extremes(over='white')     # matplotlib >= 3.4
        except AttributeError:
            cmap = cmap.copy()
            cmap.set_over('white')
        counts = np.ma.masked_where(den['counts'] == 0, den['counts'])
        pcm = ax.pcolormesh(den['x_edges'], den['y_edges'], counts.T, cmap=cmap,
                            shading='flat',
                            norm=mcolors.LogNorm(vmin=self.cmin, vmax=self.cmax))

        self.ax_cbar.set_axis_on()
        cbar = self.fig.colorbar(pcm, cax=self.ax_cbar, extend='max')
        cbar.set_label(f"Point density\n[count / {den['bin_x']:g} km"
                       f" x {den['bin_y']:g} {self.disp_unit.split('/')[0]} cell]",
                       fontsize=max(self.font_size - 2, 6))
        cbar.ax.tick_params(labelsize=max(self.font_size - 2, 6))

        ax.plot(den['x_centers'], den['median_dw'], color='red', lw=1.2,
                label='Median (density-weighted)')
        ax.plot(den['x_centers'], den['q25_dw'], color='0.55', ls='--', lw=0.8,
                label='25th-75th percentile')
        ax.plot(den['x_centers'], den['q75_dw'], color='0.55', ls='--', lw=0.8)

        for i, sw in enumerate(self.swath_list[1:], start=1):
            den_i = swath_density(sw['dist'], sw['value'], bin_x=self.bin_x,
                                  bin_y=den['bin_y'], ylim=den['ylim'],
                                  x_offset=self.x_offset)
            ax.plot(den_i['x_centers'], den_i['median_dw'], color=f'C{i}', lw=1.2,
                    label=f"Median - {os.path.basename(sw['file'])}")

        ax.axhline(0, color='k', ls='-', lw=0.4)
        if self.x_offset != 0:
            ax.axvline(0, color='blue', ls='--', lw=1.0)
        ax.set_xlim(den['x_edges'][0], den['x_edges'][-1])
        ax.set_ylim(den['ylim'])
        ax.legend(loc='lower right', fontsize=max(self.font_size - 3, 6))

        y0, x0 = swath['start_yx']
        y1, x1 = swath['end_yx']
        msg = (f"swath +/-{self.width:g} km | length {swath['length']:.2f} km"
               f" | N = {swath['value'].size:,}")
        if start_lalo is not None:
            msg += (f"\nlat/lon: ({start_lalo[0]:.4f}, {start_lalo[1]:.4f})"
                    f" --> ({end_lalo[0]:.4f}, {end_lalo[1]:.4f})")
        else:
            msg += f'\ny/x: ({y0}, {x0}) --> ({y1}, {x1})'
        ax.set_title(msg, fontsize=self.font_size)

        self.outfile_base = os.path.join(
            self.outdir, f'swath_Y{y0}X{x0}_Y{y1}X{x1}_w{self.width:g}km')
        vprint(f"swath: {swath['value'].size:,} pixels, length {swath['length']:.2f} km,"
               f" half width {swath['half_width']:.2f} km")
        if den['num_outside']:
            vprint(f"note: {den['num_outside']:,} pixels fall outside the Y range "
                   f"[{den['y_edges'][0]:g}, {den['y_edges'][-1]:g}] {self.disp_unit}"
                   f" and are not counted in the density / the weighted median.")

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return

    def draw_timeseries(self, yx):
        """Panel 3 - the point time series, via tsview.py."""
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

    ##---------- output
    def save_outputs(self):
        """Write the located swath points, the binned statistics and the figure."""
        if self.swath_list is None or self.den is None:
            print('nothing to save yet - draw a line first.')
            return
        base = self.outfile_base
        os.makedirs(os.path.dirname(base) or '.', exist_ok=True)

        # 1. the located points - one row per pixel (the LRS / swath table)
        cols = ['dist_km', 'perp_km', 'value', 'row', 'col']
        has_lalo = any('lat' in s for s in self.swath_list)
        if has_lalo:
            cols += ['lat', 'lon']
        cols += ['file_idx']

        blocks = []
        for i, sw in enumerate(self.swath_list):
            # dist_km is measured from the start of the line, minus --x-offset,
            # so the .csv, the stats and the figure share one X axis.
            arrs = [sw['dist'] - self.x_offset, sw['perp'], sw['value'], sw['row'], sw['col']]
            if has_lalo:
                nan = np.full(sw['value'].size, np.nan)
                arrs += [sw.get('lat', nan), sw.get('lon', nan)]
            arrs += [np.full(sw['value'].size, i)]
            blocks.append(np.column_stack([np.asarray(a, dtype=np.float64) for a in arrs]))
        table = np.vstack(blocks)
        fmt = (['%.5f', '%.5f', '%.6f', '%d', '%d']
               + (['%.8f', '%.8f'] if has_lalo else []) + ['%d'])
        out_pts = f'{base}_points.csv'
        np.savetxt(out_pts, table, fmt=fmt, delimiter=',', header=','.join(cols), comments='')
        print(f'saved {table.shape[0]:,} swath points to: {out_pts}')

        # 2. the binned statistics
        den = self.den
        stat = np.column_stack([den['x_centers'], den['count'], den['median'], den['mean'],
                                den['std'], den['q25'], den['q75'], den['median_dw'],
                                den['q25_dw'], den['q75_dw']])
        hdr = ('dist_km,count,median,mean,std,q25,q75,'
               'median_density_weighted,q25_density_weighted,q75_density_weighted')
        out_stat = f'{base}_stats.csv'
        np.savetxt(out_stat, stat, fmt='%.6f', delimiter=',', header=hdr, comments='')
        print(f'saved the binned statistics to: {out_stat}')

        # 3. everything needed to reproduce the profile
        sw = self.swath_list[0]
        with open(f'{base}_info.txt', 'w') as fid:
            fid.write(f'file           : {[os.path.abspath(f) for f in self.file]}\n')
            fid.write(f'dataset        : {self.dset}\n')
            fid.write(f'unit           : {self.disp_unit}\n')
            fid.write(f"start y/x      : {sw['start_yx']}\n")
            fid.write(f"end   y/x      : {sw['end_yx']}\n")
            if self.start_lalo:
                fid.write(f'start lat/lon  : {self.start_lalo}\n')
                fid.write(f'end   lat/lon  : {self.end_lalo}\n')
            fid.write(f'swath halfwidth: {self.width} km  (flat-capped buffer)\n')
            fid.write(f"profile length : {sw['length']:.4f} km\n")
            fid.write(f"num of pixels  : {sw['value'].size}\n")
            fid.write(f"bin size       : {den['bin_x']} km x {den['bin_y']} {self.disp_unit}\n")
            fid.write(f'x offset       : {self.x_offset} km  (already applied to dist_km)\n')
            fid.write('\ndist_km = distance along the profile from the start point,'
                      ' minus the x offset\n')
            fid.write('perp_km = distance from the centre line,'
                      ' + = left of the start->end direction\n')
        print(f'saved the profile metadata to: {base}_info.txt')

        # 4. the figure
        out_fig = self.outfile[0] if self.outfile else f'{base}{self.fig_ext}'
        self.fig.savefig(out_fig, bbox_inches='tight', dpi=self.fig_dpi)
        print(f'saved the figure to: {out_fig}')
        return

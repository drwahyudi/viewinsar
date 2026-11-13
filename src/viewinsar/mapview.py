"""The map, without MintPy: the stand-in for view.py's prep_slice / plot_slice.

With MintPy installed the map comes from view.py itself, and every view.py
option works. Without it this module draws it, with the options in most use:

    -c/--colormap  -v/--vlim  -u/--unit  -m/--mask (FILE or no)  --zero-mask
    --sub-lat  --sub-lon  --sub-y  --sub-x  -d/--dem  --dem-noshade
    --alpha  --ref-lalo  --ref-yx  --noreference  --fontsize  --title
    --notitle  --nocbar  --cbar-label  --flip-lr  --flip-ud  --noflip
    --noverbose

Other view.py options are reported and left out. As in view.py: the display
unit is cm (cm/year) for MintPy's velocity and time series, the mask is
maskTempCoh.h5 beside them, radar grids are flipped as the satellite sees
the ground (ascending up-down, descending left-right), and -d shades the
relief under the data.
"""

import argparse
import os
import shlex

import numpy as np

from .compat import (
    MINTPY_HINT,
    Coordinate,
    auto_disp_unit,
    auto_flip_direction,
    read,
    read_attribute,
    read_mask,
    scale_data2disp_unit,
)

__all__ = ['prep_slice', 'plot_slice', 'view_functions']

GRAY_DATASETS = ('coherence', 'temporalCoherence', 'waterMask', 'shadowMask', 'mask',
                 'height', 'avgSpatialCoh')


def view_functions():
    """(prep_slice, plot_slice): view.py's with MintPy, else these."""
    from .compat import have_mintpy
    if have_mintpy():
        from mintpy import view
        return view.prep_slice, view.plot_slice
    return prep_slice, plot_slice


def _parser():
    p = argparse.ArgumentParser(prog='view.py (viewinsar, without MintPy)', add_help=False,
                                allow_abbrev=False)
    p.add_argument('file')
    p.add_argument('dset', nargs='?')
    p.add_argument('-c', '--colormap', '--cmap', dest='colormap')
    p.add_argument('-v', '--vlim', dest='vlim', type=float, nargs=2)
    p.add_argument('-u', '--unit', dest='disp_unit')
    p.add_argument('-m', '--mask', dest='mask_file')
    p.add_argument('--zero-mask', dest='zero_mask', action='store_true')
    p.add_argument('--sub-lat', '--subset-lat', dest='subset_lat', type=float, nargs=2)
    p.add_argument('--sub-lon', '--subset-lon', dest='subset_lon', type=float, nargs=2)
    p.add_argument('--sub-y', '--subset-y', dest='subset_y', type=int, nargs=2)
    p.add_argument('--sub-x', '--subset-x', dest='subset_x', type=int, nargs=2)
    p.add_argument('-d', '--dem', dest='dem_file')
    p.add_argument('--dem-noshade', dest='disp_dem_shade', action='store_false')
    p.add_argument('--dem-nocontour', dest='dem_nocontour', action='store_true')
    p.add_argument('--alpha', '--transparency', dest='transparency', type=float)
    p.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2)
    p.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2)
    p.add_argument('--noreference', dest='disp_ref_pixel', action='store_false')
    p.add_argument('--fontsize', '--font-size', dest='font_size', type=float, default=12.0)
    p.add_argument('--title', '--figtitle', dest='fig_title')
    p.add_argument('--notitle', dest='disp_title', action='store_false')
    p.add_argument('--nocbar', '--nocolorbar', dest='disp_cbar', action='store_false')
    p.add_argument('--cbar-label', dest='cbar_label')
    p.add_argument('--flip-lr', dest='flip_lr', action='store_true')
    p.add_argument('--flip-ud', dest='flip_ud', action='store_true')
    p.add_argument('--noflip', dest='auto_flip', action='store_false')
    p.add_argument('--lalo-label', dest='lalo_label', action='store_true')
    p.add_argument('--coord', dest='fig_coord', choices=['radar', 'geo', 'yx'])
    p.add_argument('--noverbose', dest='print_msg', action='store_false')
    return p


def prep_slice(cmd, print_msg=None):
    """Read and prepare the map, as view.prep_slice: (data, metadata, inps).

    cmd - 'view.py FILE [DSET] [options]', or the same as a list without
          'view.py' (safer: file names may hold spaces)
    """
    args = shlex.split(cmd)[1:] if isinstance(cmd, str) else list(cmd)
    inps, unknown = _parser().parse_known_args(args)
    if print_msg is not None:
        inps.print_msg = print_msg
    vprint = print if inps.print_msg else (lambda *a, **k: None)
    if unknown:
        print('NOTE: without MintPy the map leaves out these view.py options: '
              f'{" ".join(unknown)} ({MINTPY_HINT} for all of them)')

    fname, dset = inps.file, inps.dset
    atr = read_attribute(fname, datasetName=dset)
    atr.setdefault('FILE_PATH', os.path.abspath(fname))
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    coord = Coordinate(atr)
    geocoded = 'Y_FIRST' in atr

    # the displayed window
    pix_box = (0, 0, width, length)
    if (inps.subset_lat or inps.subset_lon) and not geocoded:
        print('WARNING: --sub-lat / --sub-lon need geocoded data --> left out')
    elif inps.subset_lat or inps.subset_lon:
        west, north, east, south = coord.box_pixel2geo(pix_box)
        if inps.subset_lat:
            south, north = min(inps.subset_lat), max(inps.subset_lat)
        if inps.subset_lon:
            west, east = min(inps.subset_lon), max(inps.subset_lon)
        if coord.lat_step > 0:
            north, south = south, north
        pix_box = coord.box_geo2pixel((west, north, east, south))
    if inps.subset_y:
        pix_box = (pix_box[0], min(inps.subset_y), pix_box[2], max(inps.subset_y))
    if inps.subset_x:
        pix_box = (min(inps.subset_x), pix_box[1], max(inps.subset_x), pix_box[3])
    pix_box = coord.check_box_within_data_coverage(pix_box, print_msg=inps.print_msg)
    geo_box = coord.box_pixel2geo(pix_box)
    vprint(f'data   coverage in y/x: {(0, 0, width, length)}')
    vprint(f'subset coverage in y/x: {pix_box}')

    data = np.array(read(fname, datasetName=dset, box=pix_box)[0], dtype=np.float32)

    # a new reference point
    ref_yx = list(inps.ref_yx) if inps.ref_yx else None
    if inps.ref_lalo:
        if geocoded:
            ref_yx = list(coord.lalo2yx(*inps.ref_lalo))
        else:
            print('WARNING: --ref-lalo on radar-coded data needs MintPy --> left out')
    if ref_yx:
        y, x = ref_yx
        if not (0 <= y < length and 0 <= x < width):
            print(f'WARNING: the reference point {ref_yx} is off the data --> left out')
            ref_yx = None
        else:
            ref = float(np.asarray(read(fname, datasetName=dset, box=(x, y, x + 1, y + 1))[0],
                                   dtype=np.float64).ravel()[0])
            if np.isfinite(ref):
                data -= np.float32(ref)
                atr['REF_Y'], atr['REF_X'] = str(y), str(x)
                vprint(f'referenced to y/x {ref_yx}')
            else:
                print(f'WARNING: no data at the reference point {ref_yx} --> left out')
                ref_yx = None

    # the mask
    mask_file = None if str(inps.mask_file).lower() == 'no' else inps.mask_file
    msk = None
    if str(inps.mask_file).lower() != 'no':
        msk, mask_file = read_mask(fname, mask_file=mask_file, datasetName=dset, box=pix_box,
                                   print_msg=inps.print_msg)
    if msk is not None:
        data[~msk] = np.nan
    if inps.zero_mask:
        data[data == 0] = np.nan

    # the unit
    disp_unit = inps.disp_unit or auto_disp_unit(atr)
    data, disp_unit, _ = scale_data2disp_unit(data, metadata=atr, disp_unit=disp_unit)
    data = np.asarray(data, dtype=np.float32)

    # colours, range, orientation, title
    name = str(dset or atr.get('FILE_TYPE', '')).split('-')[0]
    if not inps.colormap:
        inps.colormap = 'gray' if name in GRAY_DATASETS else 'jet'
    if not inps.vlim:
        if name in ('coherence', 'temporalCoherence', 'avgSpatialCoh'):
            inps.vlim = [0.0, 1.0]
        elif np.isfinite(data).any():
            inps.vlim = [float(np.nanmin(data)), float(np.nanmax(data))]
        else:
            inps.vlim = [0.0, 1.0]
    if inps.auto_flip and not (inps.flip_lr or inps.flip_ud):
        inps.flip_lr, inps.flip_ud = auto_flip_direction(atr)
    if not inps.fig_title:
        base = os.path.splitext(os.path.basename(fname))[0]
        inps.fig_title = dset if dset and dset.split('-')[0] != atr.get('FILE_TYPE') \
            or (dset and '-' in dset) else base
    if inps.transparency is None:
        inps.transparency = 0.8 if inps.dem_file and inps.disp_dem_shade else 1.0
    vprint(f'colormap: {inps.colormap}')
    vprint(f'display range: {inps.vlim} {disp_unit}')

    inps.disp_unit, inps.mask_file = disp_unit, mask_file
    inps.pix_box, inps.geo_box, inps.coord = pix_box, geo_box, coord
    inps.data_box = (0, 0, width, length)
    inps.fig_coord = 'geo' if geocoded and inps.fig_coord in (None, 'geo') else 'radar'
    inps.map_proj_obj = None
    inps.ref_yx, inps.key, inps.dset = ref_yx, atr.get('FILE_TYPE'), dset
    return data, atr, inps


def _relief(ax, dem_file, atr, pix_box, extent, geo):
    """Shaded relief from a DEM, under the data (zorder 0)."""
    from matplotlib.colors import LightSource

    from .raster import read_any, read_attr_any
    from .swath import get_pixel_size

    datr = read_attr_any(dem_file)
    same = ((datr['LENGTH'], datr['WIDTH']) == (atr['LENGTH'], atr['WIDTH'])
            and ('Y_FIRST' in datr) == ('Y_FIRST' in atr)
            and all(abs(float(datr.get(k, 0)) - float(atr.get(k, 0))) < 1e-9 * max(
                1.0, abs(float(atr.get(k, 0)))) for k in ('X_FIRST', 'Y_FIRST', 'X_STEP', 'Y_STEP')
                if k in atr))
    if same:
        box, dem_extent = pix_box, extent
    elif geo and 'Y_FIRST' in datr:
        xf, yf = float(datr['X_FIRST']), float(datr['Y_FIRST'])
        dx, dy = float(datr['X_STEP']), float(datr['Y_STEP'])
        west, east = sorted(extent[:2])
        south, north = sorted(extent[2:])
        cols = sorted([(west - xf) / dx, (east - xf) / dx])
        rows = sorted([(north - yf) / dy, (south - yf) / dy])
        box = (max(0, int(np.floor(cols[0]))), max(0, int(np.floor(rows[0]))),
               min(int(datr['WIDTH']), int(np.ceil(cols[1]))),
               min(int(datr['LENGTH']), int(np.ceil(rows[1]))))
        if box[2] <= box[0] or box[3] <= box[1]:
            raise ValueError('the DEM does not cover the map')
        dem_extent = (xf + box[0] * dx, xf + box[2] * dx, yf + box[3] * dy, yf + box[1] * dy)
    else:
        raise ValueError('the DEM is on another grid than the radar-coded data')
    dem = np.array(read_any(dem_file, box=box)[0], dtype=np.float64)
    step = max(1, int(np.ceil(max(dem.shape) / 2000)))
    dem = dem[::step, ::step]
    dem[(dem < -1000) | (dem > 10000)] = np.nan
    if not np.isfinite(dem).any():
        raise ValueError('the DEM has no valid height there')
    x_step, y_step, unit = get_pixel_size(datr)
    shade = LightSource(azdeg=315, altdeg=45).hillshade(
        np.where(np.isfinite(dem), dem, np.nanmedian(dem)), vert_exag=0.7,
        dx=x_step * step if unit == 'm' else 30.0, dy=y_step * step if unit == 'm' else 30.0)
    ax.imshow(np.ma.masked_where(~np.isfinite(dem), shade), cmap='gray', vmin=0, vmax=1,
              extent=dem_extent, origin='upper', interpolation='nearest', zorder=0)


def plot_slice(ax, data, atr, inps):
    """Draw the map on ax, as view.plot_slice: (ax, inps, image, colorbar)."""
    from matplotlib import pyplot as plt

    from .swath import get_pixel_size, is_projected

    geo = 'Y_FIRST' in atr and getattr(inps, 'fig_coord', 'geo') == 'geo'
    x0, y0, x1, y1 = inps.pix_box
    if geo:
        xf, yf = float(atr['X_FIRST']), float(atr['Y_FIRST'])
        dx, dy = float(atr['X_STEP']), float(atr['Y_STEP'])
        extent = (xf + x0 * dx, xf + x1 * dx, yf + y1 * dy, yf + y0 * dy)
    else:
        extent = (x0 - 0.5, x1 - 0.5, y1 - 0.5, y0 - 0.5)
    font = getattr(inps, 'font_size', 12) or 12

    dem_file = getattr(inps, 'dem_file', None)
    if dem_file and getattr(inps, 'disp_dem_shade', True):
        try:
            _relief(ax, dem_file, atr, inps.pix_box, extent, geo)
        except Exception as err:          # noqa: BLE001 - the relief is a nicety
            print(f'WARNING: no shaded relief from {dem_file}: {err}')

    try:
        cmap = plt.get_cmap(inps.colormap)
    except ValueError:
        print(f'WARNING: colormap {inps.colormap!r} needs MintPy --> jet')
        cmap = plt.get_cmap('jet')
    vmin, vmax = inps.vlim
    image = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent, origin='upper',
                      interpolation='nearest', alpha=getattr(inps, 'transparency', 1.0),
                      zorder=1)

    # the ground's proportions
    if geo and not is_projected(atr):
        lat_c = np.radians((extent[2] + extent[3]) / 2)
        ax.set_aspect(1.0 / max(np.cos(lat_c), 0.05))
    elif geo:
        ax.set_aspect('equal')
    else:
        x_step, y_step, unit = get_pixel_size(atr)
        ax.set_aspect(y_step / x_step if unit == 'm' and x_step else 'equal')

    # the reference point, a black square
    if getattr(inps, 'disp_ref_pixel', True) and 'REF_Y' in atr and 'REF_X' in atr:
        try:
            ry, rx = int(float(atr['REF_Y'])), int(float(atr['REF_X']))
            if y0 <= ry < y1 and x0 <= rx < x1:
                px, py = ((xf + (rx + 0.5) * dx, yf + (ry + 0.5) * dy) if geo else (rx, ry))
                ax.plot(px, py, 'ks', ms=6, zorder=3)
        except (TypeError, ValueError):
            pass

    cbar = None
    if getattr(inps, 'disp_cbar', True):
        cbar = ax.figure.colorbar(image, ax=ax, pad=0.02, fraction=0.04)
        cbar.set_label(getattr(inps, 'cbar_label', None) or inps.disp_unit, fontsize=font)
        cbar.ax.tick_params(labelsize=font)
    if getattr(inps, 'disp_title', True):
        ax.set_title(getattr(inps, 'fig_title', ''), fontsize=font)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    if getattr(inps, 'flip_lr', False):
        ax.invert_xaxis()
    if getattr(inps, 'flip_ud', False):
        ax.invert_yaxis()
    ax.tick_params(labelsize=font)
    return ax, inps, image, cbar

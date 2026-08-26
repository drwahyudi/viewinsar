"""Command line interface of viewinsar.

    viewinsar swath velocity.h5 -w 15
    viewinsar list  velocity.h5
    viewinsar --version

The first argument may also be a file, in which case `swath` is implied:

    viewinsar velocity.h5 -w 15
"""

import argparse
import sys

from . import __version__

PROG = 'viewinsar'

GMT_FILE = """GMT xy file, i.e. transect_lonlat.xy:
>
131.1663    33.1157
131.2621    33.0860
"""

EXAMPLE = """example:
  # interactive: click 2 points on the map to draw the swath centre line
  viewinsar swath velocity.h5
  viewinsar swath geo_velocity_msk.h5 -w 15 --bin-x 0.5 --bin-y 1.5 --ylim -40 40

  # non-interactive: give the line up front, save the figure + the .csv files
  viewinsar swath geo_velocity.h5 --start-lalo 30.125 129.988 \\
      --end-lalo 30.250 130.116 -w 15 --save --nodisplay
  viewinsar swath velocity.h5 --start-yx 5290 5579 --end-yx 12177 482 -w 20 --save
  viewinsar swath velocity.h5 --line-file transect_lonlat.xy --save

  # a 3rd panel with the point time series - right-click a pixel on the map
  viewinsar swath geo_velocity.h5 --timeseries geo_timeseries_ERA5_demErr.h5

  # every view.py option is forwarded to the map panel
  viewinsar swath geo_velocity.h5 -c jet -v -3 3 -u cm -d srtm1.dem \\
      --mask geo_maskTempCoh.h5 --sub-lat 30.0 30.5 --sub-lon 129.9 130.3

mouse / keyboard:
  left  click x2 : set the start / end point of the swath centre line
  right click    : plot the time series of that pixel  (needs --timeseries)
  key "e"        : export the swath points (.csv), the binned stats (.csv)
                   and the figure
  key "n"        : clear the current line and start a new one

Any option not listed below is forwarded verbatim to MintPy's view.py for the
map panel: -c/--colormap  -v/--vlim  -u/--unit  -d/--dem  -m/--mask
--sub-lat/--sub-lon  --ref-lalo  --coastline  --lalo-label  --fontsize ...
Run `view.py -h` for the full list.
"""


##############################################################################
def create_swath_parser(prog=f'{PROG} swath'):
    """The argument parser of the `swath` sub-command."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description='Interactive SWATH (buffered) profile viewer, built on MintPy - '
                    'view.py + plot_transection.py + tsview.py',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='input file(s) to show the swath profile of')
    parser.add_argument('--dset', dest='dset',
                        help='dataset name to read, e.g. velocity, timeseries-20200101')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='disable the verbose message printing.')

    line = parser.add_argument_group('Profile location',
                                     'start/end point of the swath centre line')
    line.add_argument('--start-yx', '--yx0', dest='start_yx', metavar=('Y0', 'X0'),
                      type=int, nargs=2,
                      help='start point of the profile, in pixel number [y, x]')
    line.add_argument('--end-yx', '--yx1', dest='end_yx', metavar=('Y1', 'X1'),
                      type=int, nargs=2,
                      help='end   point of the profile, in pixel number [y, x]')
    line.add_argument('--start-lalo', '--lalo0', dest='start_lalo', metavar=('LAT0', 'LON0'),
                      type=float, nargs=2, help='start point of the profile, in [lat, lon]')
    line.add_argument('--end-lalo', '--lalo1', dest='end_lalo', metavar=('LAT1', 'LON1'),
                      type=float, nargs=2, help='end   point of the profile, in [lat, lon]')
    line.add_argument('--line-file', dest='lola_file',
                      help='file with the start/end point in lon lat, GMT format.\n' + GMT_FILE)

    swath = parser.add_argument_group('Swath', 'buffer + linear referencing settings')
    swath.add_argument('-w', '--width', '--buffer', dest='width', type=float, default=15.0,
                       metavar='KM',
                       help='HALF width of the swath buffer, in km (default: %(default)s).\n'
                            'The corridor is flat-capped, i.e. it covers\n'
                            '  |off-profile distance| <= WIDTH  AND  0 <= along-profile <= length.')
    swath.add_argument('--bin-x', dest='bin_x', type=float, default=0.5, metavar='KM',
                       help='along-profile bin size of the density grid, in km '
                            '(default: %(default)s).')
    swath.add_argument('--bin-y', dest='bin_y', type=float, default=None, metavar='NUM',
                       help='value bin size of the density grid, in the display unit\n'
                            '(default: auto, ~55 bins over the Y range).')
    swath.add_argument('--ylim', dest='ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'),
                       help='value range of the profile plot '
                            '(default: the 0.5-99.5%% percentile).')
    swath.add_argument('--cmin', dest='cmin', type=float, default=1.0, metavar='NUM',
                       help='min point density of the log colour scale (default: %(default)s).')
    swath.add_argument('--cmax', dest='cmax', type=float, default=100.0, metavar='NUM',
                       help='max point density of the log colour scale (default: %(default)s).')
    swath.add_argument('--pcmap', dest='pcmap', default='viridis', metavar='CMAP',
                       help='colormap of the point density (default: %(default)s).')
    swath.add_argument('--x-offset', dest='x_offset', type=float, default=0.0, metavar='KM',
                       help='shift the along-profile distance by this much, in km, so a known\n'
                            'feature (e.g. a fault) sits at x = 0 (default: %(default)s = off).\n'
                            'x = 0 is always a cell boundary of the density grid.')
    swath.add_argument('--keep-zero', dest='keep_zero', action='store_true',
                       help='keep the pixels with an exactly-zero value (dropped by default,\n'
                            'to match the behaviour of plot_transection.py).')
    swath.add_argument('--decimate', dest='decimate', type=int, default=1, metavar='NUM',
                       help='use every Nth pixel of the corridor only (default: %(default)s).')

    ts = parser.add_argument_group('Time series', 'optional 3rd panel, via tsview.py')
    ts.add_argument('--timeseries', '--ts', dest='timeseries', nargs='+', metavar='FILE',
                    help='time-series file(s); right-click a pixel on the map to plot it.')
    ts.add_argument('--yx', dest='yx', type=int, nargs=2, metavar=('Y', 'X'),
                    help='initial pixel of the time-series panel.')

    out = parser.add_argument_group('Save/Output', 'save the figure and the swath data')
    out.add_argument('-o', '--outfile', type=str, nargs='*',
                     help='file name of the saved figure.')
    out.add_argument('--outdir', dest='outdir', default='.',
                     help='output directory (default: %(default)s).')
    out.add_argument('--save', dest='save_fig', action='store_true',
                     help='save the figure + the swath .csv files.')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not display the figure.')
    out.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='size of the whole figure, in inches.')
    return parser


def read_lonlat_file(fname):
    """Read the start/end lat/lon from a GMT-format lonlat text file."""
    with open(fname) as fid:
        lines = [ln for ln in fid.read().splitlines()
                 if ln.strip() and not ln.lstrip().startswith('>')]
    lon0, lat0 = (float(i) for i in lines[0].split()[:2])
    lon1, lat1 = (float(i) for i in lines[1].split()[:2])
    return [lat0, lon0], [lat1, lon1]


def cmd_line_parse(iargs=None):
    """Parse the `swath` arguments; unknown ones are forwarded to view.py."""
    from mintpy.utils import readfile
    from mintpy.utils import utils as ut

    parser = create_swath_parser()
    inps, view_extra = parser.parse_known_args(args=iargs)
    inps.argv = list(iargs) if iargs is not None else sys.argv[1:]
    inps.view_extra = view_extra

    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    if inps.width <= 0:
        raise ValueError(f'--width must be positive, got {inps.width}')
    if inps.bin_x <= 0:
        raise ValueError(f'--bin-x must be positive, got {inps.bin_x}')

    inps.file = ut.get_file_list(inps.file)
    inps.atr = readfile.read_attribute(inps.file[0])
    inps.coord = ut.coordinate(inps.atr)
    inps.num_file = len(inps.file)
    if not inps.dset:
        inps.dset = readfile.get_slice_list(inps.file[0])[0]

    if inps.lola_file:
        inps.start_lalo, inps.end_lalo = read_lonlat_file(inps.lola_file)
    if inps.start_lalo and inps.end_lalo:
        (y0, y1), (x0, x1) = inps.coord.lalo2yx([inps.start_lalo[0], inps.end_lalo[0]],
                                                [inps.start_lalo[1], inps.end_lalo[1]])
        inps.start_yx, inps.end_yx = [y0, x0], [y1, x1]

    if view_extra and inps.print_msg:
        print(f'options forwarded to view.py: {" ".join(view_extra)}')
    return inps


def get_view_cmd(inps):
    """Assemble the view.py command line for the map panel."""
    cmd = f'view.py {inps.file[0]} '
    cmd += f'{inps.dset} ' if inps.dset else ''
    cmd += ' '.join(inps.view_extra)
    if not inps.print_msg:
        cmd += ' --noverbose'
    return cmd


def swath_main(iargs=None):
    """Entry point of `viewinsar swath`."""
    from .viewer import swathProfileViewer

    inps = cmd_line_parse(iargs)
    obj = swathProfileViewer(inps, get_view_cmd(inps))
    obj.open()
    obj.plot()
    if obj.fig is not None:
        obj.fig.canvas.mpl_disconnect(obj.cid)
        obj.fig.canvas.mpl_disconnect(obj.kid)
    return obj


def list_main(iargs=None):
    """Entry point of `viewinsar list` - print the datasets of a file."""
    from mintpy.utils import readfile

    parser = argparse.ArgumentParser(
        prog=f'{PROG} list', description='List the datasets available in a file.')
    parser.add_argument('file', help='input file')
    parser.add_argument('-n', '--num', type=int, default=0, metavar='NUM',
                        help='print the first NUM names only (default: all).')
    inps = parser.parse_args(args=iargs)

    names = readfile.get_slice_list(inps.file)
    for name in (names[:inps.num] if inps.num > 0 else names):
        print(name)
    if inps.num > 0 and len(names) > inps.num:
        print(f'... and {len(names) - inps.num} more (total {len(names)})')
    return names


SUBCOMMANDS = {
    'swath': (swath_main, 'interactive swath (buffered) profile viewer'),
    'list': (list_main, 'list the datasets available in a file'),
}


def main(argv=None):
    """Top-level `viewinsar` entry point."""
    argv = list(sys.argv[1:] if argv is None else argv)

    if not argv or argv[0] in ('-h', '--help', 'help'):
        width = max(len(k) for k in SUBCOMMANDS)
        print(f'usage: {PROG} <command> [options]\n\ncommands:')
        for name, (_, helptxt) in SUBCOMMANDS.items():
            print(f'  {name:<{width}}  {helptxt}')
        print(f'\nRun `{PROG} <command> -h` for the options of a command.')
        print(f'`{PROG} FILE ...` is a shortcut for `{PROG} swath FILE ...`.')
        print(f'\nviewinsar {__version__}')
        return 0

    if argv[0] in ('-V', '--version'):
        print(f'viewinsar {__version__}')
        return 0

    # "viewinsar velocity.h5 ..." is a shortcut for "viewinsar swath velocity.h5 ..."
    if argv[0] not in SUBCOMMANDS:
        argv.insert(0, 'swath')

    SUBCOMMANDS[argv[0]][0](argv[1:])
    return 0


if __name__ == '__main__':
    sys.exit(main())

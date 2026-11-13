"""Command line interface of viewinsar.

    viewinsar swath velocity.h5 -w 15
    viewinsar list  velocity.h5
    viewinsar --version

The first argument may also be a file, in which case `swath` is implied:

    viewinsar velocity.h5 -w 15
"""

import argparse
import os
import re
import sys

from . import __version__
from .terrain import AUTO_DEM, is_auto

PROG = 'viewinsar'

GMT_FILE = """GMT xy file, e.g. transect_lonlat.xy:
  >
  131.1663    33.1157
  131.2621    33.0860"""

DESCRIPTION = """\
Interactive SWATH (buffered) profile viewer for InSAR, built on MintPy.

Draw a line on the map with two clicks. The line is buffered by +/-WIDTH km
(flat caps); every valid pixel inside that corridor is projected onto it to get
its distance along the profile - QGIS LRS / ArcGIS "Locate Features Along
Routes" - and the cloud is drawn as a 2-D point-density profile with the
density-weighted median (red) and the 25th-75th percentile envelope.

The map, every cross-section and the point time series each open in their OWN
window. Press "=" to add another section: they are labelled A-A', B-B', C-C' ...
on the map and on each profile. Optional: an elevation profile under the LOS
panel (--topoxsect), fault traces from any vector file on the map, marked on
every section that crosses them (--fault), oversampling for coarse data
(--nlook), vector output (--svg).

No DEM at hand? -d and --topoxsect without a file download one for the map
extent from the Terrain Tiles (the "Mapzen Global Terrain" of QGIS).

Radar-coded data (no geocoding) work too: MintPy's geometryRadar.h5, found
next to the data or given with --lookup, supplies the DEM (its height) and
places fault traces and geocoded DEMs on the radar grid.

-nlook, -topoxsect, -fault and -svg work with one dash or two.
"""

EXAMPLE = """\
examples:
  # interactive: click 2 points on the map to draw a section
  viewinsar swath velocity.h5
  viewinsar swath geo_velocity_msk.h5 -w 15 --bin-x 0.5 --bin-y 1.5 --ylim -40 40

  # non-interactive: give the line up front, save everything, show nothing
  viewinsar swath geo_velocity.h5 --start-lalo 30.125 129.988 \\
      --end-lalo 30.250 130.116 -w 15 --save --nodisplay
  viewinsar swath velocity.h5 --start-yx 5290 5579 --end-yx 12177 482 -w 20 --save
  viewinsar swath velocity.h5 --line-file transect_lonlat.xy --save

  # a UTM file takes lat/lon OR northing/easting - these two are the same line
  viewinsar swath geo_velocity_utm.h5 --start-lalo -2.2329 120.4382 \\
      --end-lalo -1.9263 120.9328
  viewinsar swath geo_velocity_utm.h5 --start-lalo 9752950 215050 \\
      --end-lalo 9786950 270050

  # a coarse / multilooked raster: oversample the corridor so it fills in
  viewinsar swath geo_velocity_ml10.h5 --nlook 10

  # an elevation profile under the LOS panel (DEM in any format and any CRS)
  viewinsar swath geo_velocity.h5 --topoxsect srtm1.dem
  viewinsar swath geo_velocity.h5 --topoxsect dem_utm.tif --nlook 4

  # no DEM file: download one for the map extent (Terrain Tiles, cached) and use
  # it for the shaded-relief map AND the elevation profile - or the profile only
  viewinsar swath geo_velocity.h5 -d --topoxsect
  viewinsar swath geo_velocity.h5 --topoxsect
  # finer tiles for a small area (zoom 13 ~ 19 m)
  viewinsar swath geo_velocity.h5 -d --topoxsect --dem-zoom 13 \\
      --sub-lat -1.6 -1.2 --sub-lon 119.9 120.4
  # your own DEM for both: -d names it once, --topoxsect alone re-uses it
  viewinsar swath geo_velocity.h5 -d srtm1.dem --topoxsect

  # RADAR-coded data: inputs/geometryRadar.h5 is found next to the data (or in
  # ../inputs/) and gives the DEM, the elevation profile and the fault positions
  viewinsar swath velocity.h5 -d --topoxsect -fault faults.shp
  viewinsar swath velocity.h5 --lookup ../inputs/geometryRadar.h5 \\
      --topoxsect dem_wgs84.tif -fault faults.shp

  # fault traces from a vector file (.shp .geojson .gpkg .kml .gmt ...), any CRS:
  # drawn on the map, and a dashed "Fault line" where a section crosses one
  viewinsar swath geo_velocity.h5 -fault faults.shp
  viewinsar swath geo_velocity.h5 -fault gem_active_faults.geojson -topoxsect dem.tif
  viewinsar swath geo_velocity_utm.h5 -fault faults_wgs84.gpkg mapped_2024.kml \\
      --fault-label "Main fault" --fault-color darkred

  # vector figures (.svg) for publication
  viewinsar swath geo_velocity.h5 --topoxsect dem.tif --svg --save

  # the point time series, in its own window: right-click a pixel on the map
  viewinsar swath geo_velocity.h5 --timeseries geo_timeseries_ERA5_demErr.h5

  # every view.py option is forwarded to the map
  viewinsar swath geo_velocity.h5 -c jet -v -3 3 -u cm -d srtm1.dem \\
      --mask geo_maskTempCoh.h5 --sub-lat 30.0 30.5 --sub-lon 129.9 130.3

mouse / keyboard (the keys work in ANY window):
  left  click x2 : set the start / end point of the current section
  right click    : plot the time series of that pixel   (needs --timeseries)
  "=" or "+"     : start ANOTHER section, in a new window
  "n"            : clear the CURRENT section and draw it again
  "e"            : export the map and every section (see "output" below)

several cross-sections, step by step:
  1. left-click twice on the map  -> section A-A' opens in its own window
  2. press "="                    -> an empty window for section B-B' opens
  3. left-click twice on the map  -> B-B' is drawn; A-A' stays on the map
  4. press "=" again for C-C', D-D', ... as many as you need
  5. press "e"                    -> saves the map with every line and label,
                                     plus each section

output (written by "e" or --save, into --outdir):
  swath_map_ABC.png         the map, with every section line + label
  and, for each section X (A, B, C ...), named after its line, e.g.
  swath_A_Y470X150_Y130X700_w15km  (start row/col, end row/col, half width):
    <name>.png              the section figure; .svg with --svg
    <name>_points.csv       one row per sample: dist_km, perp_km, value,
                            row, col, lat, lon, file_idx
    <name>_stats.csv        per distance bin: count, median, mean, std, q25,
                            q75, the density-weighted median / q25 / q75 the
                            figure draws, and elevation_median_m (--topoxsect)
    <name>_info.txt         everything needed to reproduce the section, where
                            it crosses each fault trace (--fault), and the
                            DEM, its CRS, sampling and source (--topoxsect)
  dist_km runs from the section start (X) to its end (X'), minus --x-offset;
  perp_km is the distance from the line, + = left of the start -> end direction.

map options:
  Any option not listed above is forwarded verbatim to MintPy's view.py for the
  map: -c/--colormap  -v/--vlim  -u/--unit  -m/--mask  --sub-lat  --sub-lon
  --ref-lalo  --coastline  --lalo-label  --alpha  --dem-noshade  --fontsize ...
  -d/--dem FILE goes to view.py as before; -d alone downloads (see above).
  Run `view.py -h` for the full list.

downloaded DEM:
  Terrain Tiles on AWS Open Data, terrarium PNG tiles in Web Mercator:
  https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png
  SRTM, GMTED2010 and ETOPO1 (bathymetry) worldwide, finer national DEMs where
  available; please cite them - the attribution is printed and written to
  _info.txt. The tiles over the map are mosaicked, re-gridded onto the data's
  CRS and saved as a MintPy geometry .h5 (height, m) in ~/.cache/viewinsar,
  re-used by later runs. Needs internet access (HTTPS_PROXY is honoured).

radar-coded data:
  view.py draws them in pixel rows / columns (flipped left-right for a
  descending track). Everything lon/lat goes through MintPy's geometry file
  geometryRadar.h5 - the latitude / longitude / height of every pixel, on the
  same grid - found in inputs/, ./ or ../inputs/ next to the data or the
  current folder, or named with --lookup:
    -d, --topoxsect (no file)   its height: the DEM of your own processing,
                                already on the radar grid - nothing downloaded
                                (only if it has none, Terrain Tiles resampled
                                onto the radar grid through its lat/lon)
    --topoxsect dem.tif         a geocoded DEM, reached through its lat/lon
    --fault faults.shp          traces mapped from lon/lat to radar pixels;
                                _info.txt gives each crossing's row/col and
                                lat/lon
"""


##############################################################################
def create_swath_parser(prog=f'{PROG} swath'):
    """The argument parser of the `swath` sub-command."""
    # allow_abbrev=False: an unknown option is forwarded to view.py, so it must
    # never be taken for an abbreviation of ours (view.py's -n is not -nlook)
    parser = argparse.ArgumentParser(
        prog=prog, description=DESCRIPTION, epilog=EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)

    parser.add_argument('file', nargs='+',
                        help='input file(s): velocity.h5, a time-series epoch, a\n'
                             'GeoTIFF ... With several files the density comes from\n'
                             'the first and the others are overlaid as median lines.')
    parser.add_argument('--dset', dest='dset',
                        help='dataset to read, e.g. velocity, timeseries-20200101.\n'
                             'Default: the one matching FILE_TYPE (so "velocity",\n'
                             'never "intercept"); for a time series, the LAST epoch.\n'
                             '`viewinsar list FILE` shows what is available.')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='disable the verbose message printing.')

    line = parser.add_argument_group(
        'Profile location',
        'the first section, given up front instead of by clicking')
    line.add_argument('--start-yx', '--yx0', dest='start_yx', metavar=('Y0', 'X0'),
                      type=int, nargs=2,
                      help='start point, in pixel number [row, col]')
    line.add_argument('--end-yx', '--yx1', dest='end_yx', metavar=('Y1', 'X1'),
                      type=int, nargs=2,
                      help='end   point, in pixel number [row, col]')
    line.add_argument('--start-lalo', '--lalo0', dest='start_lalo', metavar=('LAT0', 'LON0'),
                      type=float, nargs=2,
                      help='start point in [lat, lon]. On a UTM / projected file\n'
                           '[northing, easting] in metres works too - the two\n'
                           'are told apart by their size.')
    line.add_argument('--end-lalo', '--lalo1', dest='end_lalo', metavar=('LAT1', 'LON1'),
                      type=float, nargs=2,
                      help='end   point, same convention as --start-lalo.')
    line.add_argument('--line-file', dest='lola_file', metavar='FILE',
                      help='start/end point as "lon lat" lines, GMT format.\n' + GMT_FILE)

    swath = parser.add_argument_group('Swath', 'the buffer and the density grid')
    swath.add_argument('-w', '--width', '--buffer', dest='width', type=float, default=15.0,
                       metavar='KM',
                       help='HALF width of the buffer, in km (default: %(default)s).\n'
                            'Flat caps: a pixel is in when |off-profile| <= WIDTH\n'
                            'and 0 <= along-profile <= the line length.')
    swath.add_argument('--bin-x', dest='bin_x', type=float, default=0.5, metavar='KM',
                       help='along-profile cell size, in km (default: %(default)s).')
    swath.add_argument('--bin-y', dest='bin_y', type=float, default=None, metavar='NUM',
                       help='value cell size, in the display unit\n'
                            '(default: auto, ~55 cells over the Y range).')
    swath.add_argument('--ylim', dest='ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'),
                       help='value range of the profile\n'
                            '(default: the 0.5-99.5%% percentile of the data).')
    swath.add_argument('--cmin', dest='cmin', type=float, default=1.0, metavar='NUM',
                       help='lower limit of the log density colour scale\n'
                            '(default: %(default)s).')
    swath.add_argument('--cmax', dest='cmax', type=float, default=None, metavar='NUM',
                       help='upper limit of the log density colour scale.\n'
                            'Default: auto, the 98th percentile of the occupied\n'
                            'cells - so a sparse profile stays readable and\n'
                            '--nlook does not turn the panel white.')
    swath.add_argument('--pcmap', dest='pcmap', default='viridis', metavar='CMAP',
                       help='colormap of the point density (default: %(default)s).')
    swath.add_argument('--x-offset', dest='x_offset', type=float, default=0.0, metavar='KM',
                       help='shift the distances by KM so a known feature (a\n'
                            'fault) sits at x = 0, marked by a dashed line\n'
                            '(default: %(default)s = off). x = 0 is always a cell edge.\n'
                            'With --fault, _info.txt lists the crossing distance.')
    swath.add_argument('--keep-zero', dest='keep_zero', action='store_true',
                       help='keep exactly-zero pixels (dropped by default, as in\n'
                            'plot_transection.py). --topoxsect always keeps them:\n'
                            '0 m is sea level.')
    swath.add_argument('--decimate', dest='decimate', type=int, default=1, metavar='NUM',
                       help='use every NUMth pixel only, to speed up huge swaths\n'
                            '(default: %(default)s). The opposite of --nlook.')
    swath.add_argument('--nlook', '-nlook', dest='nlook', type=int, default=1, metavar='NUM',
                       help='OVERSAMPLE the corridor NUM times per axis, so NUM**2\n'
                            'more samples, by NaN-aware bilinear interpolation\n'
                            '(default: %(default)s = native pixels). Fills in the sparse\n'
                            'profile of a multilooked raster: data looked by 10\n'
                            'needs about --nlook 10. It restores sampling density,\n'
                            'not information. Cannot be combined with --decimate.')

    topo = parser.add_argument_group(
        'Topography', 'the shaded-relief map background, and an elevation profile under\n'
                      'the LOS panel - from a DEM file, or downloaded for the map extent')
    topo.add_argument('-d', '--dem', dest='dem_file', nargs='?', const=AUTO_DEM, metavar='DEM',
                      help='DEM for the shaded-relief background of the map, as in\n'
                           'view.py. WITHOUT a file, it is downloaded for the map\n'
                           'extent from the Terrain Tiles (AWS Open Data - the\n'
                           '"Mapzen Global Terrain" layer of QGIS) and cached in\n'
                           '~/.cache/viewinsar ($VIEWINSAR_CACHE), so a second run\n'
                           'downloads nothing. On RADAR-coded data it is the height\n'
                           'of geometryRadar.h5 instead (see --lookup). --topoxsect\n'
                           'without a file uses the same DEM.')
    topo.add_argument('--topoxsect', '-topoxsect', '--topo', dest='topoxsect', nargs='?',
                      const=AUTO_DEM, metavar='DEM',
                      help='DEM to profile along the SAME swath, drawn as a black\n'
                           'median line in a panel under the LOS profile, sharing\n'
                           'its distance axis. WITHOUT a file: the -d DEM, or the\n'
                           'Terrain Tiles download if -d names no file. Any format\n'
                           '(MintPy .h5, GeoTIFF, .dem + .rsc, any GDAL raster) and\n'
                           'any CRS - read from the file itself; the line is\n'
                           're-projected when the DEM and the data differ. On\n'
                           'radar-coded data: geometryRadar.h5 itself, or a geocoded\n'
                           'DEM reached through it (--lookup). A DEM covering part of\n'
                           'a section gives that part. Sampled at its own\n'
                           'resolution, not --nlook (thinned\n'
                           'past ~4 million samples per section). A section that\n'
                           'gets no topography says why in its panel. Written to\n'
                           '_stats.csv as elevation_median_m, the exact median.')
    topo.add_argument('--dem-zoom', dest='dem_zoom', type=int, metavar='Z',
                      help='zoom level of the downloaded tiles, 0-15 (default: auto,\n'
                           'the finest that keeps the download within ~300 tiles,\n'
                           'at most 12). Ground pixel ~ 156 km / 2**Z at the equator:\n'
                           '10 ~ 150 m, 12 ~ 38 m, 13 ~ 19 m.')
    topo.add_argument('--lookup', dest='lookup_file', metavar='FILE',
                      help='RADAR-coded data only: MintPy geometry file on the same\n'
                           'grid - geometryRadar.h5, with latitude / longitude /\n'
                           'height of every pixel. Its height is the DEM of -d and\n'
                           '--topoxsect given without a file (the DEM of your own\n'
                           'processing, so nothing is downloaded); its lat/lon put\n'
                           '--fault traces and a geocoded DEM file on the radar grid.\n'
                           'Default: found as MintPy does - inputs/geometryRadar.h5,\n'
                           './geometryRadar.h5, ../inputs/geometryRadar.h5 - next to\n'
                           'the data file, then in the current folder.')

    flt = parser.add_argument_group('Fault', 'fault traces on the map and on every section')
    flt.add_argument('--fault', '-fault', dest='fault_file', nargs='+', metavar='FILE',
                     help='fault trace(s) from any geocoded vector file, in any CRS:\n'
                          '.shp (+ .prj), .geojson, .gpkg, .kml/.kmz, GMT text\n'
                          '(.gmt .lonlat .xy) - and everything else GDAL reads\n'
                          'once pyogrio, fiona or GDAL is installed. Lines and\n'
                          'polygon outlines are drawn on the map, re-projected to\n'
                          'the data CRS (a WGS84 file on a UTM map is fine) - or,\n'
                          'on radar-coded data, mapped onto the radar grid through\n'
                          'geometryRadar.h5 (see --lookup).\n'
                          'Wherever a section crosses a trace, a dashed vertical\n'
                          'line marks the spot on its profile, annotated "Fault\n'
                          'line", and the distance goes to _info.txt. Several\n'
                          'files may be given. (view.py\'s --faultline, forwarded as\n'
                          'usual, only draws .shp / .lonlat on a lat/lon map and\n'
                          'marks no profile.)')
    flt.add_argument('--fault-label', dest='fault_label', default='Fault line', metavar='TEXT',
                     help='annotation of a crossing on the profile\n'
                          '(default: "%(default)s"; "" for none).')
    flt.add_argument('--fault-color', dest='fault_color', default='k', metavar='COLOR',
                     help='colour of the traces and of the crossing lines, any\n'
                          'matplotlib colour (default: %(default)s = black).')
    flt.add_argument('--fault-lw', dest='fault_lw', type=float, default=1.2, metavar='NUM',
                     help='line width of the traces on the map (default: %(default)s).')

    ts = parser.add_argument_group('Time series', 'a point time-series window, via tsview.py')
    ts.add_argument('--timeseries', '--ts', dest='timeseries', nargs='+', metavar='FILE',
                    help='time-series file(s); right-click a pixel on the map to\n'
                         'plot its displacement history in its own window.')
    ts.add_argument('--yx', dest='yx', type=int, nargs=2, metavar=('Y', 'X'),
                    help='pixel to plot in the time-series window at start-up.')

    out = parser.add_argument_group('Save/Output', 'figures and data - see "output" below')
    out.add_argument('-o', '--outfile', type=str, nargs='*',
                     help='file name of the saved MAP figure. Sections are named\n'
                          'after their label: swath_A_..., swath_B_...')
    out.add_argument('--outdir', dest='outdir', default='.',
                     help='output directory (default: %(default)s).')
    out.add_argument('--save', dest='save_fig', action='store_true',
                     help='save the map + every section at start-up, as "e" does.')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not open any window (implies --save).')
    out.add_argument('--svg', '-svg', dest='svg', action='store_true',
                     help='write the figures as vector .svg instead of .png.\n'
                          'Shorthand for --figext .svg.')
    out.add_argument('--figext', dest='fig_ext', default='.png', metavar='EXT',
                     choices=['.png', '.pdf', '.svg', '.eps', '.jpg'],
                     help='figure format: .png .pdf .svg .eps .jpg\n'
                          '(default: %(default)s).')
    out.add_argument('--dpi', dest='fig_dpi', type=int, default=300, metavar='NUM',
                     help='resolution of raster figures (default: %(default)s).')
    out.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='size of each SECTION window, in inches\n'
                          '(default: 7.6 4.8, taller with --topoxsect).')
    out.add_argument('--figsize-map', dest='fig_size_map', metavar=('WID', 'LEN'),
                     type=float, nargs=2,
                     help='size of the MAP window, in inches\n'
                          '(default: from the data shape).')
    return parser


def read_lonlat_file(fname):
    """Read the start/end lat/lon from a GMT-format lonlat text file."""
    with open(fname) as fid:
        lines = [ln for ln in fid.read().splitlines()
                 if ln.strip() and not ln.lstrip().startswith('>')]
    lon0, lat0 = (float(i) for i in lines[0].split()[:2])
    lon1, lat1 = (float(i) for i in lines[1].split()[:2])
    return [lat0, lon0], [lat1, lon1]


def default_dset(fname, atr=None, print_msg=True):
    """The dataset to show when the user gave no --dset.

    A MintPy velocity.h5 holds `intercept`, `interceptStd`, `velocity` and
    `velocityStd`, and get_slice_list returns them in file order, so blindly
    taking the first can hand back `intercept` - the y-intercept of the linear
    fit, in metres, anti-correlated with the velocity. Match FILE_TYPE instead.
    For a time-series file the first slice is the all-zero reference date, so
    use the LAST epoch, the cumulative displacement, as tsview.py does.
    """
    from mintpy.utils import readfile

    atr = atr or readfile.read_attribute(fname)
    slices = readfile.get_slice_list(fname)
    if not slices:
        raise ValueError(f'no 2D dataset found in {fname}')
    names = [s.split('-')[0] for s in slices]
    ftype = atr.get('FILE_TYPE', '').replace('.', '')

    if ftype in ('timeseries', 'giantTimeseries', 'HDFEOS'):
        base = slices[0].split('-')[0]
        epochs = [s for s in slices if s.split('-')[0] == base]
        if len(epochs) > 1:
            if print_msg:
                print(f'no --dset given: showing the last epoch, {epochs[-1]} '
                      f'(the cumulative displacement). `viewinsar list {fname}` '
                      'lists the others.')
            return epochs[-1]

    for want in ([ftype] if ftype else []) + ['velocity', 'displacement', 'unwrapPhase',
                                              'temporalCoherence', 'coherence', 'height']:
        if want and want in names:
            pick = slices[names.index(want)]
            others = sorted({n for n in names if n != want})
            if print_msg and others:
                print(f'no --dset given: using "{pick}" (matches FILE_TYPE = '
                      f'{ftype or "?"}). Others in this file: {others}')
            return pick

    if print_msg and len(set(names)) > 1:
        print(f'no --dset given and nothing matches FILE_TYPE = "{ftype}"; '
              f'using "{slices[0]}". Available: {sorted(set(names))}')
    return slices[0]


def check_map_options(inps):
    """-d, --dem-zoom, --lookup and --fault - the same for swath and stat."""
    # -d FILE stays view.py's background, exactly as before; -d alone = download
    if is_auto(inps.dem_file):
        inps.dem_file = AUTO_DEM
    elif inps.dem_file:
        if not os.path.isfile(inps.dem_file):
            raise ValueError(f'-d/--dem: no such file: {inps.dem_file}')
        inps.view_extra = list(inps.view_extra) + ['-d', inps.dem_file]
    if inps.dem_zoom is not None and not 0 <= inps.dem_zoom <= 15:
        raise ValueError(f'--dem-zoom must be 0-15, got {inps.dem_zoom}')
    if inps.lookup_file and not os.path.isfile(inps.lookup_file):
        raise ValueError(f'--lookup: no such file: {inps.lookup_file}')
    if inps.fault_file:
        # os.path.exists, not isfile: a FileGDB (.gdb) is a directory
        missing = [f for f in inps.fault_file if not os.path.exists(f)]
        if missing:
            raise ValueError(f'--fault: no such file: {", ".join(missing)}')
    if inps.fault_lw <= 0:
        raise ValueError(f'--fault-lw must be positive, got {inps.fault_lw}')
    return inps


def read_file_info(inps):
    """The input file list, its metadata and the dataset to show."""
    from mintpy.utils import readfile
    from mintpy.utils import utils as ut

    inps.file = ut.get_file_list(inps.file)
    if not inps.file:
        raise ValueError('no such input file')
    inps.atr = readfile.read_attribute(inps.file[0])
    inps.coord = ut.coordinate(inps.atr)
    inps.num_file = len(inps.file)
    if not inps.dset:
        inps.dset = default_dset(inps.file[0], inps.atr, print_msg=inps.print_msg)
    return inps


def cmd_line_parse(iargs=None):
    """Parse the `swath` arguments; unknown ones are forwarded to view.py."""
    parser = create_swath_parser()
    inps, view_extra = parser.parse_known_args(args=iargs)
    inps.argv = list(iargs) if iargs is not None else sys.argv[1:]
    inps.view_extra = view_extra

    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    if inps.svg:
        inps.fig_ext = '.svg'
    if inps.width <= 0:
        raise ValueError(f'--width must be positive, got {inps.width}')
    if inps.bin_x <= 0:
        raise ValueError(f'--bin-x must be positive, got {inps.bin_x}')
    if inps.nlook < 1:
        raise ValueError(f'--nlook must be >= 1, got {inps.nlook}')
    if inps.nlook > 1 and inps.decimate > 1:
        raise ValueError('--nlook and --decimate are opposites; use one or the other.')
    if inps.cmax is not None and inps.cmax <= inps.cmin:
        raise ValueError(f'--cmax ({inps.cmax}) must exceed --cmin ({inps.cmin})')
    check_map_options(inps)
    if is_auto(inps.topoxsect):
        # --topoxsect alone: the -d file if there is one, else the download
        inps.topoxsect = inps.dem_file if inps.dem_file else AUTO_DEM
    elif inps.topoxsect and not os.path.isfile(inps.topoxsect):
        raise ValueError(f'--topoxsect: no such file: {inps.topoxsect}')

    read_file_info(inps)

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
    if obj.fig_map is not None and getattr(obj, 'cid', None) is not None:
        obj.fig_map.canvas.mpl_disconnect(obj.cid)
    if not obj.disp_fig:
        obj.close()          # nothing is on screen, so do not leak the figures
    return obj


##############################################################################
#                                 stat                                       #
##############################################################################
STAT_DESCRIPTION = """\
Block statistics on both sides of a line - typically drawn along a fault.

Draw the line on the map (left-click its points, Enter to finish), or give it
with --line-lalo / --line-yx / --line-file. Every SPACING km along it a point
is placed (centred on the line), and at each point two circles of RADIUS km,
one on each side, their edges GAP km (default 1) clear of the line. On the
inside of a bend a circle is pushed further out, so no circle ever reaches
across the line or within GAP of it. Block 1 is the LEFT side of the line,
walking from its first point to its last; block 2 the RIGHT.

The values in each circle are drawn as a violin - block 1 and block 2 side by
side at every point, each in the colour of its circles on the map - with
strikes at the median (in the block's colour), the mean (solid) and the mode
(dotted: the peak of the density, where the violin is widest). The medians of
each block are joined by a dashed line along the line.
"""

STAT_EXAMPLE = """\
examples:
  # draw the line along a fault shown on the map: circles every 5 km, 2 km across
  viewinsar stat velocity.h5 -sp 5 -r 2 -fault faults.shp
  viewinsar stat geo_velocity.h5 -sp 10 -r 3 --gap 1 -d -c RdBu_r -v -40 40

  # the line given up front: lat/lon or y/x pairs, or the longest line of a
  # vector file (.shp .geojson .gpkg .kml .gmt ..., e.g. drawn in QGIS)
  viewinsar stat geo_velocity.h5 -sp 5 -r 2 --save --nodisplay \\
      --line-lalo -1.20 119.95 -1.55 120.05 -1.80 120.20
  viewinsar stat velocity.h5 -sp 5 -r 2 --line-yx 3670 744 3300 720 2967 696
  viewinsar stat velocity.h5 -sp 5 -r 2 --line-file fault_trace.geojson --svg --save

  # your own names and colours for the two blocks
  viewinsar stat velocity.h5 -sp 5 -r 2 --names "West block" "East block" \\
      --colors tab:green tab:purple

mouse / keyboard (in either window):
  left click       : add a point to the line
  Enter            : finish the line -> points, circles and violins
  "u" or Backspace : remove the last point
  "n"              : start a new line (so does a click after Enter)
  "e"              : export (see "output")

output (written by "e" or --save, into --outdir), named after the line's first
and last point, e.g. stat_Y3670X744_Y2967X696_sp5km_r2km:
  <name>.png          the violins (.svg with --svg)
  <name>_map.png      the map, with the line, the points and the circles
  <name>_stats.csv    one row per point: dist_km along the line, its position,
                      count / median / mean / mode / std / q25 / q75 / min /
                      max of each block (_b1, _b2), the circle offsets, and
                      median_b1_minus_b2
  <name>_values.csv   every pixel used: point, block, value, row, col,
                      dist_to_line_km (never below the gap), and lat/lon
  <name>_info.txt     the line, the parameters, and which side is which block

map options: anything not listed above goes to view.py, as for `viewinsar swath`
(-c/--colormap -v/--vlim -u/--unit -m/--mask --sub-lat --sub-lon ...).
"""


def create_stat_parser(prog=f'{PROG} stat'):
    """The argument parser of the `stat` sub-command."""
    parser = argparse.ArgumentParser(
        prog=prog, description=STAT_DESCRIPTION, epilog=STAT_EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)
    parser.add_argument('file', help='input file: velocity.h5, a GeoTIFF, a time-series\n'
                                     'epoch ... geocoded or in radar coordinates.')
    parser.add_argument('--dset', dest='dset',
                        help='dataset to read (default: the one matching FILE_TYPE;\n'
                             'see `viewinsar list FILE`).')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='disable the verbose message printing.')

    line = parser.add_argument_group('Line', 'given up front instead of drawn on the map')
    line.add_argument('--line-yx', dest='line_yx', type=float, nargs='+', metavar='Y X',
                      help='its points as row col pairs: Y0 X0 Y1 X1 [Y2 X2 ...].')
    line.add_argument('--line-lalo', dest='line_lalo', type=float, nargs='+',
                      metavar='LAT LON',
                      help='its points as lat lon pairs (or northing easting on a\n'
                           'UTM file; through geometryRadar.h5 on radar-coded data).')
    line.add_argument('--line-file', dest='line_file', metavar='FILE',
                      help='the longest line in a vector file - any format and CRS\n'
                           'that --fault reads (.shp .geojson .gpkg .kml .gmt ...).')

    circ = parser.add_argument_group('Circles', 'where the blocks are sampled')
    circ.add_argument('-sp', '--spacing', dest='spacing', type=float, default=5.0,
                      metavar='KM',
                      help='distance between the points along the line, in km\n'
                           '(default: %(default)s). The points are centred on it.')
    circ.add_argument('-r', '--radius', dest='radius', type=float, default=2.0, metavar='KM',
                      help='radius of the circles, in km (default: %(default)s).')
    circ.add_argument('--gap', dest='gap', type=float, default=1.0, metavar='KM',
                      help='space between the line and every circle, in km\n'
                           '(default: %(default)s): a circle centre sits RADIUS + GAP\n'
                           'from its point, more on the inside of a bend.')
    circ.add_argument('--keep-zero', dest='keep_zero', action='store_true',
                      help='keep exactly-zero pixels (dropped by default).')

    look = parser.add_argument_group('Plot', 'the violins')
    look.add_argument('--names', dest='block_names', nargs=2, metavar=('NAME1', 'NAME2'),
                      default=['Block 1', 'Block 2'],
                      help='names of block 1 (left) and block 2 (right)\n'
                           '(default: "Block 1" "Block 2").')
    look.add_argument('--colors', dest='block_colors', nargs=2, metavar=('COLOR1', 'COLOR2'),
                      default=['#0072B2', '#D55E00'],
                      help='colours of block 1 and block 2, any matplotlib colour\n'
                           '(default: blue #0072B2 and vermilion #D55E00).')
    look.add_argument('--ylim', dest='ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'),
                      help='value range of the violin plot (default: auto).')

    mapg = parser.add_argument_group('Map', 'the guides drawn on the map')
    mapg.add_argument('-d', '--dem', dest='dem_file', nargs='?', const=AUTO_DEM, metavar='DEM',
                      help='DEM for the shaded-relief background, as in view.py.\n'
                           'Without a file: downloaded for the map (Terrain Tiles),\n'
                           'or the height of geometryRadar.h5 on radar-coded data.')
    mapg.add_argument('--dem-zoom', dest='dem_zoom', type=int, metavar='Z',
                      help='zoom of the downloaded tiles, 0-15 (default: auto).')
    mapg.add_argument('--lookup', dest='lookup_file', metavar='FILE',
                      help='radar-coded data: MintPy geometryRadar.h5 (default:\n'
                           'found in inputs/, ./ or ../inputs/ as MintPy does).')
    mapg.add_argument('--fault', '-fault', dest='fault_file', nargs='+', metavar='FILE',
                      help='fault traces to draw on the map, to follow with the line;\n'
                           'any vector format and CRS, as in `viewinsar swath`.')
    mapg.add_argument('--fault-color', dest='fault_color', default='k', metavar='COLOR',
                      help='colour of the fault traces (default: %(default)s).')
    mapg.add_argument('--fault-lw', dest='fault_lw', type=float, default=1.2, metavar='NUM',
                      help='line width of the fault traces (default: %(default)s).')

    out = parser.add_argument_group('Save/Output', 'figures and data - see "output" below')
    out.add_argument('-o', '--outfile', type=str, nargs='*',
                     help='file name of the saved MAP figure.')
    out.add_argument('--outdir', dest='outdir', default='.',
                     help='output directory (default: %(default)s).')
    out.add_argument('--save', dest='save_fig', action='store_true',
                     help='save everything at start-up (needs a --line-* option).')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not open any window (implies --save).')
    out.add_argument('--svg', '-svg', dest='svg', action='store_true',
                     help='write the figures as vector .svg instead of .png.')
    out.add_argument('--figext', dest='fig_ext', default='.png', metavar='EXT',
                     choices=['.png', '.pdf', '.svg', '.eps', '.jpg'],
                     help='figure format: .png .pdf .svg .eps .jpg (default: %(default)s).')
    out.add_argument('--dpi', dest='fig_dpi', type=int, default=300, metavar='NUM',
                     help='resolution of raster figures (default: %(default)s).')
    out.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='size of the violin window, in inches (default: 9 4.8).')
    out.add_argument('--figsize-map', dest='fig_size_map', metavar=('WID', 'LEN'),
                     type=float, nargs=2,
                     help='size of the map window, in inches (default: auto).')
    return parser


def _pairs(values, option):
    if values is None:
        return None
    if len(values) < 4 or len(values) % 2:
        raise ValueError(f'{option} needs at least two points, as pairs of numbers; '
                         f'got {len(values)} number(s)')
    return [list(values[i:i + 2]) for i in range(0, len(values), 2)]


def stat_cmd_line_parse(iargs=None):
    """Parse the `stat` arguments; unknown ones are forwarded to view.py."""
    from matplotlib.colors import is_color_like

    parser = create_stat_parser()
    inps, view_extra = parser.parse_known_args(args=iargs)
    inps.argv = list(iargs) if iargs is not None else sys.argv[1:]
    inps.view_extra = view_extra
    inps.file = [inps.file]

    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    if inps.svg:
        inps.fig_ext = '.svg'
    if inps.spacing <= 0:
        raise ValueError(f'-sp/--spacing must be positive, got {inps.spacing}')
    if inps.radius <= 0:
        raise ValueError(f'-r/--radius must be positive, got {inps.radius}')
    if inps.gap < 0:
        raise ValueError(f'--gap must not be negative, got {inps.gap}')
    bad = [c for c in inps.block_colors if not is_color_like(c)]
    if bad:
        raise ValueError(f'--colors: not a colour: {", ".join(bad)}')
    inps.line_yx = _pairs(inps.line_yx, '--line-yx')
    inps.line_lalo = _pairs(inps.line_lalo, '--line-lalo')
    given = [opt for opt in ('line_yx', 'line_lalo', 'line_file') if getattr(inps, opt)]
    if len(given) > 1:
        raise ValueError('give the line once: --line-yx, --line-lalo OR --line-file')
    if inps.line_file and not os.path.exists(inps.line_file):
        raise ValueError(f'--line-file: no such file: {inps.line_file}')
    if not inps.disp_fig and not given:
        raise ValueError('--nodisplay needs the line up front '
                         '(--line-yx, --line-lalo or --line-file)')
    check_map_options(inps)
    read_file_info(inps)

    # what the shared map machinery of the swath viewer expects to find
    for key, value in dict(nlook=1, topoxsect=None, timeseries=None, yx=None,
                           start_yx=None, end_yx=None, start_lalo=None, end_lalo=None,
                           fault_label='').items():
        setattr(inps, key, value)

    if view_extra and inps.print_msg:
        print(f'options forwarded to view.py: {" ".join(view_extra)}')
    return inps


def stat_main(iargs=None):
    """Entry point of `viewinsar stat`."""
    from .statview import StatViewer

    inps = stat_cmd_line_parse(iargs)
    obj = StatViewer(inps, get_view_cmd(inps))
    obj.open()
    obj.plot()
    if not obj.disp_fig:
        obj.close()
    return obj


##############################################################################
#                                tseries                                     #
##############################################################################
TSERIES_DESCRIPTION = """\
Interactive time-series viewer, after LiCSBAS_plot_ts.py (and MintPy's
tsview.py), for LiCSBAS cum_filt.h5 / cum.h5 and MintPy time series
(timeseries*.h5, geo_timeseries*.h5, HDF-EOS5 .he5).

Two windows, as in LiCSBAS. The map shows the velocity relative to the
reference area - or a noise index (LiCSBAS results/; MintPy mask, temporal
coherence, height), or, with the slider, the cumulative displacement at a
date. Double-click a pixel for its time series, with LiCSBAS's fit models
(Linear, Annual+L, Quad, Annual+Q). Right-drag a new reference area.

With -r KM a double click takes a CIRCLE of that radius: every pixel within
it, relative to the reference area and the reference date as LiCSBAS does a
single pixel, and each date shows their MEDIAN, with the 25th - 75th
percentiles as an error bar. The fits follow the medians.
"""

TSERIES_EXAMPLE = """\
examples:
  # LiCSBAS: in TS_GEOCml*/ (cum_filt.h5, and cum.h5 in red, as LiCSBAS does)
  viewinsar tseries
  viewinsar tseries cum_filt.h5 -r 2                 # circles of 2 km
  viewinsar tseries cum_filt.h5 -r 1.5 -m 20170105 --ref 140:145/125:130

  # MintPy (velocity.h5, maskTempCoh.h5 ... found next to it)
  viewinsar tseries timeseries_ERA5_ramp_demErr.h5 -r 2 -u cm
  viewinsar tseries geo/geo_timeseries_ERA5_demErr.h5 -r 3 --lalo -2.51 121.68

  # no window: save the time series of one circle (.csv + figures)
  viewinsar tseries cum_filt.h5 -r 2 --p-geo 121.68/-2.51 --nodisplay

mouse / keyboard:
  left double click : the time series of that pixel (of the circle, with -r)
  right drag        : a new reference area
  slider            : the cumulative displacement at that date
  "e"               : export, see below

output (written by "e", --save or --nodisplay, into --outdir), named after the
point, e.g. tseries_Y156X144_r2km:
  <name>.csv        date, year, then per file: median_mm, q25_mm, q75_mm,
                    n_px (with -r), or disp_mm (one pixel)
  <name>.png        the time series (.svg with --svg); <name>_map.png the map
  <name>_info.txt   the file(s), the point, the circle, the reference area
                    and date, the velocity

LiCSBAS users: -r here is the RADIUS of the circle; the reference area that
LiCSBAS_plot_ts.py takes with -r is --ref X1:X2/Y1:Y2. -p X/Y and
--p_geo / --ref_geo work as in LiCSBAS.
"""


def _radius_km(text):
    if re.search(r'[:/]', str(text)):
        raise argparse.ArgumentTypeError(
            f'{text!r}: -r is the circle radius in km here; give the reference area '
            '(LiCSBAS -r) with --ref X1:X2/Y1:Y2')
    try:
        value = float(text)
    except ValueError as err:
        raise argparse.ArgumentTypeError(f'{text!r} is not a number of km') from err
    if not value > 0:
        raise argparse.ArgumentTypeError(f'the radius must be positive, got {text}')
    return value


def _pair(text, cast, what):
    parts = str(text).split('/')
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(f'{text!r}: expected {what}')
    try:
        return [cast(v) for v in parts]
    except ValueError as err:
        raise argparse.ArgumentTypeError(f'{text!r}: expected {what}') from err


def create_tseries_parser(prog=f'{PROG} tseries'):
    """The argument parser of the `tseries` sub-command."""
    parser = argparse.ArgumentParser(
        prog=prog, description=TSERIES_DESCRIPTION, epilog=TSERIES_EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)
    parser.add_argument('file', nargs='?',
                        help='the time series: LiCSBAS cum_filt.h5 / cum.h5, or MintPy\n'
                             '(default: ./cum_filt.h5, ./cum.h5 or ./timeseries.h5).')
    parser.add_argument('-r', '--radius', dest='radius', type=_radius_km, metavar='KM',
                        help='a circle of that radius, in km: each date shows the\n'
                             'median of its pixels, with the 25th-75th percentiles\n'
                             'as an error bar (default: one pixel).')
    parser.add_argument('--i2', dest='file2', metavar='FILE',
                        help='a second time series, in red (default: cum.h5 when the\n'
                             'file is cum_filt.h5, as LiCSBAS).')
    parser.add_argument('--no-i2', dest='no_i2', action='store_true',
                        help='no second time series.')
    parser.add_argument('-m', '--ref-date', dest='ref_date', metavar='YYYYMMDD',
                        help='reference date, where every series is 0 (default: the\n'
                             'first date; MintPy: its REF_DATE).')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='disable the verbose message printing.')

    ref = parser.add_argument_group('Reference area', 'its mean is subtracted at every '
                                    'date (default: the file\'s)')
    ref.add_argument('--ref', '--ref-area', dest='ref_area', metavar='X1:X2/Y1:Y2',
                     help='in pixels, as LiCSBAS -r (0 for X2 / Y2 = to the end).')
    ref.add_argument('--ref-geo', '--ref_geo', dest='ref_geo', metavar='LON1/LON2/LAT1/LAT2',
                     help='in geographic coordinates, as LiCSBAS --ref_geo.')
    ref.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, metavar=('Y', 'X'),
                     help='one reference pixel, as MintPy.')
    ref.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2,
                     metavar=('LAT', 'LON'), help='one reference pixel, in lat/lon.')

    pnt = parser.add_argument_group('Point', 'the first time series shown (default: the '
                                    'centre of the reference area)')
    pnt.add_argument('--yx', dest='yx', type=int, nargs=2, metavar=('Y', 'X'),
                     help='row and column, as MintPy tsview.py.')
    pnt.add_argument('--lalo', dest='lalo', type=float, nargs=2, metavar=('LAT', 'LON'),
                     help='lat and lon (northing easting on a UTM grid), as tsview.py.')
    pnt.add_argument('-p', '--point', dest='point_xy', metavar='X/Y',
                     type=lambda t: _pair(t, int, 'X/Y in pixels'),
                     help='x/y in pixels, as LiCSBAS -p.')
    pnt.add_argument('--p-geo', '--p_geo', dest='point_geo', metavar='LON/LAT',
                     type=lambda t: _pair(t, float, 'LON/LAT'),
                     help='lon/lat, as LiCSBAS --p_geo.')

    data = parser.add_argument_group('Data')
    data.add_argument('--results', dest='results_dir', metavar='DIR',
                      help='LiCSBAS noise indices (mask, coh_avg, n_unw, vstd ...)\n'
                           '(default: results/ beside the file, as LiCSBAS -d).')
    data.add_argument('--mask', dest='mask_file', metavar='FILE',
                      help='mask file (default: LiCSBAS results/mask; MintPy\n'
                           'maskTempCoh.h5 beside the file).')
    data.add_argument('--nomask', dest='nomask', action='store_true',
                      help='start with the mask off (the box on the map turns it on).')
    data.add_argument('--vel', dest='vel_file', metavar='FILE',
                      help='velocity file (default: LiCSBAS vel; MintPy velocity.h5\n'
                           'beside the file, else a linear fit of the time series).')
    data.add_argument('--lookup', dest='lookup_file', metavar='FILE',
                      help='radar-coded MintPy data: geometryRadar.h5, for lat/lon\n'
                           '(default: found where MintPy keeps it).')

    disp = parser.add_argument_group('Display')
    disp.add_argument('-c', '--colormap', dest='colormap', default='roma_r', metavar='CMAP',
                      help='colormap of the velocity / displacement (default: roma_r,\n'
                           'as LiCSBAS; any MintPy or matplotlib name).')
    disp.add_argument('-u', '--unit', dest='unit', default='mm', choices=['mm', 'cm', 'm'],
                      help='displacement unit (default: %(default)s).')
    disp.add_argument('--vlim', dest='vlim', type=float, nargs=2, metavar=('VMIN', 'VMAX'),
                      help='colour range of the velocity (default: auto).')
    disp.add_argument('--dlim', dest='dlim', type=float, nargs=2, metavar=('DMIN', 'DMAX'),
                      help='colour range of the cumulative displacement (default: auto).')
    disp.add_argument('--auto-crange', dest='auto_crange', type=float, default=99.0,
                      metavar='PCT', help='percentile of the automatic colour ranges\n'
                                          '(default: %(default)s).')
    disp.add_argument('--ylim', dest='ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'),
                      help='range of the time-series plot (default: auto).')
    disp.add_argument('--noflip', dest='noflip', action='store_true',
                      help='radar-coded data: do not flip the map as MintPy does.')
    disp.add_argument('--figsize', dest='fig_size', type=float, nargs=2, metavar=('WID', 'LEN'),
                      help='size of the time-series window, in inches (default: 8 4.8).')
    disp.add_argument('--figsize-map', dest='fig_size_map', type=float, nargs=2,
                      metavar=('WID', 'LEN'), help='size of the map window (default: auto).')

    out = parser.add_argument_group('Save/Output', 'see "output" below')
    out.add_argument('-o', '--outfile', type=str, nargs='*',
                     help='file name of the saved time-series figure.')
    out.add_argument('--outdir', dest='outdir', default='.',
                     help='output directory (default: %(default)s).')
    out.add_argument('--save', dest='save_fig', action='store_true',
                     help='save the first time series at start-up.')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not open any window (implies --save).')
    out.add_argument('--svg', '-svg', dest='svg', action='store_true',
                     help='write the figures as vector .svg instead of .png.')
    out.add_argument('--figext', dest='fig_ext', default='.png', metavar='EXT',
                     choices=['.png', '.pdf', '.svg', '.eps', '.jpg'],
                     help='figure format: .png .pdf .svg .eps .jpg (default: %(default)s).')
    out.add_argument('--dpi', dest='fig_dpi', type=int, default=300, metavar='NUM',
                     help='resolution of raster figures (default: %(default)s).')
    return parser


def tseries_cmd_line_parse(iargs=None):
    """Parse the `tseries` arguments."""
    from .tsdata import find_time_series

    parser = create_tseries_parser()
    inps = parser.parse_args(args=iargs)
    inps.argv = list(iargs) if iargs is not None else sys.argv[1:]

    if not inps.file:
        inps.file = find_time_series('.')
        if not inps.file:
            raise ValueError('no time series given, and no cum_filt.h5, cum.h5 or '
                             'timeseries.h5 here')
    elif not os.path.isfile(inps.file):
        raise ValueError(f'no such file: {inps.file}')
    if inps.no_i2:
        inps.file2 = None
    elif inps.file2:
        if not os.path.isfile(inps.file2):
            raise ValueError(f'--i2: no such file: {inps.file2}')
    elif os.path.basename(inps.file) == 'cum_filt.h5':
        cand = os.path.join(os.path.dirname(os.path.abspath(inps.file)), 'cum.h5')
        inps.file2 = cand if os.path.isfile(cand) else None
    for opt, path in (('--mask', inps.mask_file), ('--vel', inps.vel_file),
                      ('--lookup', inps.lookup_file)):
        if path and not os.path.isfile(path):
            raise ValueError(f'{opt}: no such file: {path}')
    if inps.results_dir and not os.path.isdir(inps.results_dir):
        raise ValueError(f'--results: no such folder: {inps.results_dir}')
    if inps.ref_date and not re.fullmatch(r'\d{8}', inps.ref_date):
        raise ValueError(f'-m: a date is YYYYMMDD, got {inps.ref_date}')

    points = [o for o, v in (('--yx', inps.yx), ('--lalo', inps.lalo), ('-p', inps.point_xy),
                             ('--p-geo', inps.point_geo)) if v]
    if len(points) > 1:
        raise ValueError(f'give the point once: {" / ".join(points)}')
    refs = [o for o, v in (('--ref', inps.ref_area), ('--ref-geo', inps.ref_geo),
                           ('--ref-yx', inps.ref_yx), ('--ref-lalo', inps.ref_lalo)) if v]
    if len(refs) > 1:
        raise ValueError(f'give the reference area once: {" / ".join(refs)}')
    for opt, pair in (('--vlim', inps.vlim), ('--dlim', inps.dlim), ('--ylim', inps.ylim)):
        if pair and not pair[0] < pair[1]:
            raise ValueError(f'{opt}: the first value must be the smaller, got {pair}')
    if not 50 < inps.auto_crange <= 100:
        raise ValueError(f'--auto-crange must be within (50, 100], got {inps.auto_crange}')

    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    if inps.svg:
        inps.fig_ext = '.svg'
    return inps


def tseries_main(iargs=None):
    """Entry point of `viewinsar tseries`."""
    from .tsviewer import TimeSeriesViewer

    inps = tseries_cmd_line_parse(iargs)
    obj = TimeSeriesViewer(inps)
    obj.open()
    obj.plot()
    if not obj.disp_fig:
        obj.close()
    return obj


def list_main(iargs=None):
    """Entry point of `viewinsar list` - print the datasets of a file."""
    from mintpy.utils import readfile

    parser = argparse.ArgumentParser(
        prog=f'{PROG} list',
        description='List the datasets available in a file - the names that\n'
                    '`viewinsar swath --dset` accepts. The one swath uses by\n'
                    'default is marked with "<- default".',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('file', help='input file')
    parser.add_argument('-n', '--num', type=int, default=0, metavar='NUM',
                        help='print the first NUM names only (default: all).')
    inps = parser.parse_args(args=iargs)

    names = readfile.get_slice_list(inps.file)
    try:
        default = default_dset(inps.file, print_msg=False)
    except Exception:            # noqa: BLE001 - the marker is a nicety
        default = None
    for name in (names[:inps.num] if inps.num > 0 else names):
        print(f'{name}  <- default' if name == default else name)
    if inps.num > 0 and len(names) > inps.num:
        print(f'... and {len(names) - inps.num} more (total {len(names)})')
    return names


SUBCOMMANDS = {
    'swath': (swath_main, 'interactive swath (buffered) profile viewer'),
    'stat': (stat_main, 'block statistics on both sides of a line (violins)'),
    'tseries': (tseries_main, 'time series of a pixel or a circle (LiCSBAS / MintPy)'),
    'list': (list_main, 'list the datasets available in a file'),
}


def top_level_help():
    """What `viewinsar` / `viewinsar -h` prints."""
    width = max(len(k) for k in SUBCOMMANDS)
    cmds = '\n'.join(f'  {name:<{width}}  {helptxt}'
                     for name, (_, helptxt) in SUBCOMMANDS.items())
    return f"""\
viewinsar {__version__} - interactive swath (buffered) profiles of InSAR data,
built on MintPy.

usage: {PROG} <command> [options]
       {PROG} FILE [options]          shortcut for `{PROG} swath FILE ...`
       {PROG} --version

commands:
{cmds}

swath at a glance:
  {PROG} swath velocity.h5 -w 15                 draw sections with the mouse
  {PROG} swath vel_looked10.h5 --nlook 10        fill in a coarse raster
  {PROG} swath velocity.h5 --topoxsect dem.tif   add an elevation profile
  {PROG} swath velocity.h5 -d --topoxsect        ... from a DEM downloaded for you
  {PROG} swath velocity.h5 -fault faults.shp     faults: map + dashed on profiles
  {PROG} swath velocity.h5 --timeseries ts.h5    right-click: point history
  {PROG} swath velocity.h5 --svg --save          vector figures

  every window:  "="  new section (A-A', B-B', C-C' ...)
                 "n"  redo the current section
                 "e"  export the map + every section
                 right click on the map: the time series of that pixel

stat at a glance:
  {PROG} stat velocity.h5 -sp 5 -r 2 -fault faults.shp
        click a line along the fault, Enter: every 5 km a pair of circles
        (2 km radius, 1 km clear of the line), block 1 left / block 2 right,
        a violin per circle (median, mean and mode strikes), and the block
        medians joined by dashed lines
  {PROG} stat velocity.h5 -sp 5 -r 2 --line-file trace.geojson --save

tseries at a glance:
  {PROG} tseries cum_filt.h5 -r 2       LiCSBAS_plot_ts.py, with circles: double-click
                                        -> the median of a 2 km circle at every date,
                                        25th-75th percentiles as error bars
  {PROG} tseries timeseries.h5 -r 2     the same for a MintPy time series

Run `{PROG} <command> -h` for every option of a command."""


def main(argv=None):
    """Top-level `viewinsar` entry point."""
    argv = list(sys.argv[1:] if argv is None else argv)

    if not argv or argv[0] in ('-h', '--help', 'help'):
        print(top_level_help())
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

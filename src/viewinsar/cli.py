"""Command line interface of viewinsar.

    viewinsar swath velocity.h5 --width 15
    viewinsar list  velocity.h5
    viewinsar --version

The first argument may also be a file, in which case `swath` is implied:

    viewinsar velocity.h5 --width 15
"""

import argparse
import os
import re
import sys

from . import __version__
from .terrain import AUTO_DEM, DEM_SOURCES, is_auto

PROG = 'viewinsar'

GMT_FILE = """GMT xy file, e.g. transect_lonlat.xy - one section per segment:
  >  A-A'
  131.1663    33.1157
  131.2621    33.0860
  >  B-B'
  131.2105    33.1630
  131.3050    33.1342"""
VECTOR_EXT = ('.geojson', '.json', '.kml', '.kmz', '.shp', '.gpkg', '.gml', '.fgb', '.gpx')

DESCRIPTION = """\
Interactive SWATH (buffered) profile viewer for InSAR (MintPy, LiCSBAS ...).

Draw a line on the map with two clicks. The line is buffered by +/-WIDTH km
(flat caps); every valid pixel inside that corridor is projected onto it to get
its distance along the profile - QGIS LRS / ArcGIS "Locate Features Along
Routes" - and the cloud is drawn as a 2-D pixel-density profile with the
density-weighted median (red) and the 25th-75th percentile envelope.

The map, every cross-section and the point time series each open in their OWN
window. Press "=" to add another section: they are labelled A-A', B-B', C-C' ...
on the map and on each profile. Several sections can be given up front, too:
--profile (repeat it) or a --line-file with one line per section.

Optional: an elevation profile under the LOS panel (--topoxsect), fault traces
from any vector file on the map, marked on every section that crosses them
(--fault), oversampling for coarse data (--nlook), vector output (--svg).

No DEM at hand? --dem and --topoxsect without a file download the Copernicus
GLO-30 DEM (1", ~30 m) for the map extent.

Radar-coded data (no geocoding) work too: MintPy's geometryRadar.h5, found
next to the data or given with --lookup, supplies the DEM (its height) and
places fault traces and geocoded DEMs on the radar grid.

Several files (the frames of a track, several tracks) make one density, with a
median line for each. --gnss puts GNSS velocities on the map and, inside each
swath, on its profile - in the LOS (--vel los) or as one component.

Every option is written with two dashes: --fault, --gnss, --vel, --dem,
--width ... The one-dash spellings of before (-fault, -gnss, -vel, -d, -w ...)
still work, with a note.
"""

EXAMPLE = """\
examples:
  # interactive: click 2 points on the map to draw a section
  viewinsar swath velocity.h5
  viewinsar swath geo_velocity_msk.h5 --width 15 --bin-x 0.5 --bin-y 1.5 --ylim -40 40

  # non-interactive: give the line up front, save everything, show nothing
  viewinsar swath geo_velocity.h5 --start-lalo 30.125 129.988 \\
      --end-lalo 30.250 130.116 --width 15 --save --nodisplay
  viewinsar swath velocity.h5 --start-yx 5290 5579 --end-yx 12177 482 --width 20 --save
  viewinsar swath velocity.h5 --line-file transect_lonlat.xy --save

  # several sections up front, A-A', B-B', C-C', each in its own window
  viewinsar swath geo_velocity.h5 --profile 30.10 129.95 30.30 130.10 \\
      --profile 30.20 130.05 30.40 130.20 --profile 30.30 130.15 30.50 130.30 \\
      --width 10 --save sections --nodisplay
  # ... or one section per line of a file (GMT segments, .geojson, .shp, .kml)
  viewinsar swath velocity.h5 --line-file profiles.gmt --width 10 --save

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

  # no DEM file: download the Copernicus GLO-30 DEM for the map extent (cached)
  # and use it for the shaded-relief map AND the elevation profile - or the
  # profile only
  viewinsar swath geo_velocity.h5 --dem --topoxsect
  viewinsar swath geo_velocity.h5 --topoxsect
  # the Terrain Tiles ("Mapzen Global Terrain" of QGIS) instead, zoom 13 ~ 19 m
  viewinsar swath geo_velocity.h5 --dem --topoxsect --dem-source terrain \\
      --dem-zoom 13 --sub-lat -1.6 -1.2 --sub-lon 119.9 120.4
  # your own DEM for both: --dem names it once, --topoxsect alone re-uses it
  viewinsar swath geo_velocity.h5 --dem srtm1.dem --topoxsect

  # RADAR-coded data: inputs/geometryRadar.h5 is found next to the data (or in
  # ../inputs/) and gives the DEM, the elevation profile and the fault positions
  viewinsar swath velocity.h5 --dem --topoxsect --fault faults.shp
  viewinsar swath velocity.h5 --lookup ../inputs/geometryRadar.h5 \\
      --topoxsect dem_wgs84.tif --fault faults.shp

  # fault traces from a vector file (.shp .geojson .gpkg .kml .gmt ...), any CRS:
  # drawn on the map, and a dashed "Fault line" where a section crosses one
  viewinsar swath geo_velocity.h5 --fault faults.shp
  viewinsar swath geo_velocity.h5 --fault gem_active_faults.geojson --topoxsect dem.tif
  viewinsar swath geo_velocity_utm.h5 --fault faults_wgs84.gpkg mapped_2024.kml \\
      --fault-label "Main fault" --fault-color darkred

  # several files at once - frames of a track, or several tracks: one density
  # for all of them, a median line for each, every file on the map
  viewinsar swath T160A/geo_velocity.h5 T087A/geo_velocity.h5 --labels T160A T087A

  # GNSS on the profiles (.csv .xlsx .shp .geojson .gpkg .kml, GMT psvelo ...):
  # the stations on the map, those inside a swath on its profile - black
  # triangles 3D, circles 2D, 1-sigma error bars
  viewinsar swath geo_velocity.h5 --gnss gnss.csv --vel los      # into the LOS
  viewinsar swath geo_vel_east.h5 --gnss gnss.xlsx --vel ve      # the east component

  # vector figures (.svg) for publication
  viewinsar swath geo_velocity.h5 --topoxsect dem.tif --svg --save

  # the point time series, in its own window: right-click a pixel on the map
  viewinsar swath geo_velocity.h5 --timeseries geo_timeseries_ERA5_demErr.h5

  # every view.py option is forwarded to the map
  viewinsar swath geo_velocity.h5 --colormap jet --vlim -3 3 --unit cm --dem srtm1.dem \\
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

output (written by "e" or --save, into the --save folder - made if it does not
exist - or --outdir, else the current folder):
  swath_map_ABC.png         the map, with every section line + label
  and, for each section X (A, B, C ...), named after its line, e.g.
  swath_A_Y470X150_Y130X700_w15km  (start row/col, end row/col, half width):
    <name>.png              the section figure; .svg with --svg
    <name>_points.csv       one row per sample: dist_km, perp_km, value,
                            row, col, lat, lon, file_idx
    <name>_gnss.csv         the GNSS in the swath (--gnss): site, lon, lat,
                            dist_km, perp_km, value, sigma, value_shifted ...
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
  map: --colormap  --vlim  --unit  --mask  --sub-lat  --sub-lon
  --ref-lalo  --coastline  --lalo-label  --alpha  --dem-noshade  --fontsize ...
  --dem FILE goes to view.py as before; --dem alone downloads (see above).
  Run `view.py -h` for the full list.

downloaded DEM:
  --dem-source cop30 (the default): the Copernicus DEM GLO-30 on AWS Open Data,
  1" (~30 m), heights above the EGM2008 geoid, one Cloud Optimized GeoTIFF per
  1 x 1 deg cell: https://copernicus-dem-30m.s3.amazonaws.com/
  Only the parts over the map are read, by byte range - no GDAL needed. The
  elevation profile gets 1" pixels up to 8000 a side (a map of ~2.2 deg), the
  map background up to 3000 (~0.8 deg); a bigger map gets wider pixels, read
  from the tiles' own overviews. The sea has no tiles: 0 m there.
  --dem-source terrain: the Terrain Tiles ("Mapzen Global Terrain" of QGIS),
  terrarium PNG tiles in Web Mercator - SRTM, GMTED2010 and ETOPO1 (bathymetry);
  the fallback when GLO-30 cannot be had and no source was named:
  https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png
  Either way, please cite them: the attribution is printed and written to _info.txt.
  The download is re-gridded onto the data's CRS and saved as a MintPy geometry
  .h5 (height, m) in ~/.cache/viewinsar, re-used by later runs. Needs internet
  access (HTTPS_PROXY is honoured).

radar-coded data:
  view.py draws them in pixel rows / columns (flipped left-right for a
  descending track). Everything lon/lat goes through MintPy's geometry file
  geometryRadar.h5 - the latitude / longitude / height of every pixel, on the
  same grid - found in inputs/, ./ or ../inputs/ next to the data or the
  current folder, or named with --lookup:
    --dem, --topoxsect (no file)   its height: the DEM of your own processing,
                                already on the radar grid - nothing downloaded
                                (only if it has none, GLO-30 resampled onto
                                the radar grid through its lat/lon)
    --topoxsect dem.tif         a geocoded DEM, reached through its lat/lon
    --fault faults.shp          traces mapped from lon/lat to radar pixels;
                                _info.txt gives each crossing's row/col and
                                lat/lon
"""


class ViewinsarParser(argparse.ArgumentParser):
    """No abbreviations, of single-dash options either.

    Python before 3.13 abbreviates single-dash options even with
    allow_abbrev=False. Every option of ours is --long, and the one-dash
    spellings of before are rewritten first (modern_options): only an exact
    option may match here; anything else - view.py's -v, -u, -c ... - goes on
    to view.py.
    """

    def _get_option_tuples(self, option_string):
        tuples = super()._get_option_tuples(option_string)
        if self.allow_abbrev:
            return tuples
        return [t for t in tuples if t[1] == option_string
                or (len(t[1]) == 2 and option_string.startswith(t[1]))]


##############################################################################
#          the option spellings of before 0.5.0 - still understood           #
##############################################################################
# Every option of viewinsar is written --long since 0.5.0. The one-dash words
# (-fault, -gnss, -vel ...) and letters (-d, -w, -r ...) of before still work:
# they are rewritten here, with a note. view.py's own options (-v, -u, -c, -m
# ... for the map of swath, stat and model) are view.py's, and pass as they are.
_OLD_EVERY = {'-save': '--save', '-svg': '--svg', '-o': '--outfile'}
_OLD_MAP = {'-fault': '--fault', '-d': '--dem'}
_OLD_SECTION = {'-profile': '--profile', '-nlook': '--nlook', '-topoxsect': '--topoxsect',
                '--topo': '--topoxsect', '-gnss': '--gnss', '-vel': '--vel',
                '--gnss-comp': '--vel', '-w': '--width', '--buffer': '--width'}
OLD_SPELLINGS = {
    'swath': {**_OLD_EVERY, **_OLD_MAP, **_OLD_SECTION},
    'model': {**_OLD_EVERY, **_OLD_MAP, **_OLD_SECTION},
    'stat': {**_OLD_EVERY, **_OLD_MAP, '-sp': '--spacing'},
    'tseries': {**_OLD_EVERY, '-r': '--radius', '-m': '--ref-date', '-p': '--point',
                '-c': '--colormap', '-u': '--unit', '--vel': '--vel-file'},
    'list': {'-n': '--num'},
}
# options that are gone: left out, with a note - and the values they took
# (0.5.1: swath draws no GNSS line; 0.5.8: stat fits two lines across the fault)
_STAT_NOW = 'stat measures the creep on profiles across the line now (--length, --width)'
DROPPED = {
    'swath': {'--gnss-nofit': ('the profiles draw no line through the GNSS any more', 0)},
    'stat': {'--radius': (_STAT_NOW, 1), '-r': (_STAT_NOW, 1), '--gap': (_STAT_NOW, 1),
             '--rate': ('the creep rate is what stat draws now', 0),
             '-rate': ('the creep rate is what stat draws now', 0),
             '--rate-err': ('its error bars come from a bootstrap now (--boot)', 1)},
}


def _is_value(tok):
    """A token that is a value, not an option: a word, or a number (-2.5)."""
    tok = str(tok)
    if not tok.startswith('-'):
        return True
    try:
        float(tok)
    except ValueError:
        return False
    return True


def modern_options(argv, command):
    """`argv` with the spellings of before 0.5.0 (-fault, -gnss, -vel, -d, -w
    ...) written as today (--fault, --gnss, --vel, --dem, --width ...), and a
    note that says which, once. A one-letter option glued to its value (-w15)
    is split as argparse would."""
    table, gone = OLD_SPELLINGS.get(command, {}), DROPPED.get(command, {})
    out, said, skip = [], {}, 0
    for tok in argv:
        head, eq, tail = str(tok).partition('=')
        if skip:                        # the value of an option that is gone
            skip -= 1
            if _is_value(tok):
                continue
            skip = 0
        if tok in gone or (eq and head in gone):
            note, num = gone[head if eq else tok]
            print(f'NOTE: {head if eq else tok} is not needed any more - {note}')
            skip = 0 if eq else num
            continue
        if tok in table:
            new, old = [table[tok]], tok
        elif eq and head in table:
            new, old = [f'{table[head]}={tail}'], head
        elif len(tok) > 2 and tok[:2] in table and tok[2] not in '-=':
            new, old = [table[tok[:2]], tok[2:]], tok[:2]
        else:
            out.append(tok)
            continue
        said.setdefault(old, table[old])
        out.extend(new)
    if said:
        pairs = ', '.join(f'{new} for {old}' for old, new in said.items())
        print(f'NOTE: viewinsar\'s options are all written --long now: {pairs} '
              '(the old spellings still work)')
    return out


def _argv(iargs, command):
    """The arguments of a command, in today's spellings."""
    return modern_options(list(sys.argv[1:] if iargs is None else iargs), command)


##############################################################################
def create_swath_parser(prog=f'{PROG} swath'):
    """The argument parser of the `swath` sub-command."""
    # allow_abbrev=False: an unknown option is forwarded to view.py, so it must
    # never be taken for an abbreviation of ours (view.py's -n is not --nlook)
    parser = ViewinsarParser(
        prog=prog, description=DESCRIPTION, epilog=EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)

    parser.add_argument('file', nargs='+',
                        help='input file(s): velocity.h5, a time-series epoch, a\n'
                             'GeoTIFF ... Several files - frames of a track, several\n'
                             'tracks - make one density, with a median line for each.\n'
                             'The map is the first file, with the others drawn under it.')
    parser.add_argument('--labels', dest='labels', nargs='+', metavar='NAME',
                        help='names of the input files in the legends (default: the\n'
                             'file names).')
    parser.add_argument('--dset', dest='dset',
                        help='dataset to read, e.g. velocity, timeseries-20200101.\n'
                             'Default: the one matching FILE_TYPE (so "velocity",\n'
                             'never "intercept"); for a time series, the LAST epoch.\n'
                             '`viewinsar list FILE` shows what is available.')
    parser.add_argument('--data-unit', dest='data_unit', metavar='UNIT',
                        help='the unit the input files are really in: mm/yr, cm/yr,\n'
                             'm/yr, mm ... - for a file that does not say it (a GeoTIFF),\n'
                             'or says a wrong one (e.g. mm/yr labelled m/year). Default:\n'
                             'the file\'s own; mm/yr for a LiCSBAS velocity (vel*, vstd*)\n'
                             'without a unit, inside its LiCSBAS folder; for any other\n'
                             'velocity without one (model, --gnss, or a vel* name), m/yr\n'
                             'when 99%% of its values are below 1, else mm/yr.')
    parser.add_argument('--data-scale', dest='data_scale', type=float, nargs='+',
                        metavar='FACTOR',
                        help='multiply the data by FACTOR before anything - the map, the\n'
                             'profiles, the fit, the statistics - one per input file (or\n'
                             'one for all): 0.1 for a product 10x too large. For every\n'
                             'file of a mission, a rule in your viewinsar.ini does it -\n'
                             '~/.config/viewinsar/viewinsar.ini: [data-scale]\n'
                             'ALOS-2 = 0.1 (--data-scale 1: not this time).')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='disable the verbose message printing.')

    line = parser.add_argument_group(
        'Profile location',
        "sections given up front instead of by clicking: A-A' from\n"
        "--start-* / --end-*, then one per --line-file line, then one per\n"
        '--profile')
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
                      help='the sections as "lon lat" lines, GMT format, one per\n'
                           'segment, from its first to its last point - or any\n'
                           'vector file (.geojson .kml .shp .gpkg), one per line.\n'
                           + GMT_FILE)
    add_profile_options(line)

    swath = parser.add_argument_group('Swath', 'the buffer and the density grid')
    swath.add_argument('--width', dest='width', type=float, default=15.0,
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
                       help='colormap of the pixel density (default: %(default)s).')
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
    swath.add_argument('--nlook', dest='nlook', type=int, default=1, metavar='NUM',
                       help='OVERSAMPLE the corridor NUM times per axis, so NUM**2\n'
                            'more samples, by NaN-aware bilinear interpolation\n'
                            '(default: %(default)s = native pixels). Fills in the sparse\n'
                            'profile of a multilooked raster: data looked by 10\n'
                            'needs about --nlook 10. It restores sampling density,\n'
                            'not information. Cannot be combined with --decimate.')

    topo = parser.add_argument_group(
        'Topography', 'the shaded-relief map background, and an elevation profile under\n'
                      'the LOS panel - from a DEM file, or downloaded for the map extent')
    topo.add_argument('--dem', dest='dem_file', nargs='?', const=AUTO_DEM, metavar='DEM',
                      help='DEM for the shaded-relief background of the map, as in\n'
                           'view.py. WITHOUT a file, the Copernicus GLO-30 DEM (1",\n'
                           '~30 m, AWS Open Data) is downloaded for the map extent\n'
                           '(see --dem-source) and cached in ~/.cache/viewinsar\n'
                           '($VIEWINSAR_CACHE), so a second run downloads nothing.\n'
                           'On RADAR-coded data it is the height of geometryRadar.h5\n'
                           'instead (see --lookup). --topoxsect without a file uses\n'
                           'the same download.')
    topo.add_argument('--topoxsect', dest='topoxsect', nargs='?',
                      const=AUTO_DEM, metavar='DEM',
                      help='DEM to profile along the SAME swath, drawn as a black\n'
                           'median line in a panel under the LOS profile, sharing\n'
                           'its distance axis. WITHOUT a file: the --dem DEM, or the\n'
                           'GLO-30 download if --dem names no file. Any format\n'
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
    add_dem_source(topo)
    topo.add_argument('--dem-zoom', dest='dem_zoom', type=int, metavar='Z',
                      help='--dem-source terrain only: zoom level of the Terrain\n'
                           'Tiles, 0-15 (default: auto, the finest that keeps the\n'
                           'download within ~300 tiles, at most 12). Ground pixel\n'
                           '~ 156 km / 2**Z at the equator: 10 ~ 150 m, 12 ~ 38 m,\n'
                           '13 ~ 19 m.')
    topo.add_argument('--lookup', dest='lookup_file', metavar='FILE',
                      help='RADAR-coded data only: MintPy geometry file on the same\n'
                           'grid - geometryRadar.h5, with latitude / longitude /\n'
                           'height of every pixel. Its height is the DEM of --dem and\n'
                           '--topoxsect given without a file (the DEM of your own\n'
                           'processing, so nothing is downloaded); its lat/lon put\n'
                           '--fault traces and a geocoded DEM file on the radar grid.\n'
                           'Default: found as MintPy does - inputs/geometryRadar.h5,\n'
                           './geometryRadar.h5, ../inputs/geometryRadar.h5 - next to\n'
                           'the data file, then in the current folder.')

    flt = parser.add_argument_group('Fault', 'fault traces on the map and on every section')
    flt.add_argument('--fault', dest='fault_file', nargs='+', metavar='FILE',
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

    add_gnss_options(parser)

    ts = parser.add_argument_group('Time series', 'a point time-series window, via tsview.py')
    ts.add_argument('--timeseries', '--ts', dest='timeseries', nargs='+', metavar='FILE',
                    help='time-series file(s); right-click a pixel on the map to\n'
                         'plot its displacement history in its own window.')
    ts.add_argument('--yx', dest='yx', type=int, nargs=2, metavar=('Y', 'X'),
                    help='pixel to plot in the time-series window at start-up.')

    out = parser.add_argument_group('Save/Output', 'figures and data - see "output" below')
    out.add_argument('--outfile', type=str, nargs='*',
                     help='file name of the saved MAP figure. Sections are named\n'
                          'after their label: swath_A_..., swath_B_...')
    add_save_options(out, 'the map and every section given up front\n(--start-lalo ...)')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not open any window (implies --save).')
    out.add_argument('--svg', dest='svg', action='store_true',
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


def read_profile_file(fname):
    """Every profile of a --line-file: one per GMT segment, or per vector line.

    A GMT text file has "x y" rows - lon lat, or easting northing on a
    projected data file - with segments split by '>' lines; a vector file
    (.geojson .kml .shp .gpkg ...) is re-projected to lon/lat when it says
    its CRS. A line of more than two points runs from its first to its last.

    Returns: list of ([lat0, lon0], [lat1, lon1]), in file order.
    """
    import numpy as np

    from .fault import _read_gmt, read_fault_file
    from .raster import same_crs

    if not os.path.isfile(fname):
        raise ValueError(f'--line-file: no such file: {fname}')
    ext = os.path.splitext(fname)[1].lower()
    layers = read_fault_file(fname) if ext in VECTOR_EXT else _read_gmt(fname)
    out = []
    for lines, crs in layers:
        for line in lines:
            line = np.asarray(line, dtype=np.float64)[:, :2]
            line = line[np.isfinite(line).all(axis=1)]
            if len(line) < 2:
                continue
            ends = line[[0, -1]]
            if crs is not None and not same_crs(crs, 4326):
                from pyproj import Transformer
                lon, lat = Transformer.from_crs(crs, 4326, always_xy=True).transform(
                    ends[:, 0], ends[:, 1])
                ends = np.c_[lon, lat]
            out.append(([float(ends[0, 1]), float(ends[0, 0])],
                        [float(ends[1, 1]), float(ends[1, 0])]))
    if not out:
        raise ValueError(f'--line-file: no line of two or more points in {fname}')
    return out


def read_lonlat_file(fname):
    """The start / end [lat, lon] of the first line of a --line-file."""
    return read_profile_file(fname)[0]


def add_profile_options(group):
    """--profile / --profile-yx: more sections up front (swath, model)."""
    group.add_argument('--profile', dest='profile', action='append', nargs=4,
                       type=float, metavar=('LAT0', 'LON0', 'LAT1', 'LON1'),
                       help='one more section up front, from [LAT0, LON0] to\n'
                            '[LAT1, LON1] (northing / easting on a projected file).\n'
                            "Repeat it: each is the next of A-A', B-B', C-C' ...")
    group.add_argument('--profile-yx', dest='profile_yx', action='append', nargs=4, type=int,
                       metavar=('Y0', 'X0', 'Y1', 'X1'),
                       help='the same, in pixel numbers [row, col].')
    return group


def collect_profiles(inps):
    """Every section given up front, in order: --start-lalo / --end-lalo (or
    --start-yx / --end-yx), one per --line-file line, one per --profile and
    --profile-yx. Sets inps.profiles, a list of dicts (start_yx, end_yx,
    start_lalo, end_lalo), and inps.start_yx ... to the first of them."""
    for a, b in (('start_lalo', 'end_lalo'), ('start_yx', 'end_yx')):
        if bool(getattr(inps, a, None)) != bool(getattr(inps, b, None)):
            given, missing = (a, b) if getattr(inps, a, None) else (b, a)
            raise ValueError(f'--{given.replace("_", "-")} needs --{missing.replace("_", "-")}')
    profs = []

    def add(what, lalo=None, yx=None):
        if lalo is not None:
            (y0, y1), (x0, x1) = inps.coord.lalo2yx([lalo[0][0], lalo[1][0]],
                                                    [lalo[0][1], lalo[1][1]])
            yx = ([y0, x0], [y1, x1])
        start, end = list(yx[0]), list(yx[1])
        if [int(c) for c in start] == [int(c) for c in end]:
            raise ValueError(f'{what}: both ends fall on the same pixel {start}')
        profs.append(dict(start_yx=start, end_yx=end,
                          start_lalo=[float(c) for c in lalo[0]] if lalo else None,
                          end_lalo=[float(c) for c in lalo[1]] if lalo else None))

    if inps.start_lalo and inps.end_lalo:
        add('--start-lalo / --end-lalo', lalo=(inps.start_lalo, inps.end_lalo))
    elif inps.start_yx and inps.end_yx:
        add('--start-yx / --end-yx', yx=(inps.start_yx, inps.end_yx))
    if getattr(inps, 'lola_file', None):
        lines = read_profile_file(inps.lola_file)
        for k, (start, end) in enumerate(lines):
            add(f'--line-file line {k + 1}', lalo=(start, end))
        if len(lines) > 1 and inps.print_msg:
            print(f'--line-file: {len(lines)} lines -> {len(lines)} sections')
    for k, p in enumerate(getattr(inps, 'profile', None) or []):
        add(f'--profile #{k + 1}', lalo=(p[:2], p[2:]))
    for k, p in enumerate(getattr(inps, 'profile_yx', None) or []):
        add(f'--profile-yx #{k + 1}', yx=(p[:2], p[2:]))

    inps.profiles = profs
    if profs:
        first = profs[0]
        inps.start_yx, inps.end_yx = first['start_yx'], first['end_yx']
        inps.start_lalo, inps.end_lalo = first['start_lalo'], first['end_lalo']
    return profs


def default_dset(fname, atr=None, print_msg=True):
    """The dataset to show when the user gave no --dset.

    A MintPy velocity.h5 holds `intercept`, `interceptStd`, `velocity` and
    `velocityStd`, and get_slice_list returns them in file order, so blindly
    taking the first can hand back `intercept` - the y-intercept of the linear
    fit, in metres, anti-correlated with the velocity. Match FILE_TYPE instead.
    For a time-series file the first slice is the all-zero reference date, so
    use the LAST epoch, the cumulative displacement, as tsview.py does.
    """
    from .compat import get_slice_list, read_attribute

    atr = atr or read_attribute(fname)
    slices = get_slice_list(fname)
    if not slices:
        raise ValueError(f'no 2D dataset found in {fname}')
    names = [s.split('-')[0] for s in slices]
    ftype = atr.get('FILE_TYPE', '').replace('.', '')
    if atr.get('PROCESSOR') == 'LiCSBAS' and 'vel' in slices:
        if print_msg and len(slices) > 1:
            print('no --dset given: using "vel", the LiCSBAS velocity. `viewinsar list '
                  f'{fname}` lists the others (cum-YYYYMMDD ...)')
        return 'vel'

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


def add_gnss_options(parser):
    """--gnss and its companions - for swath and model."""
    grp = parser.add_argument_group('GNSS', 'GNSS velocities on the map and the profiles')
    grp.add_argument('--gnss', dest='gnss_file', metavar='FILE',
                     help='GNSS velocities: .csv .txt .xlsx .shp .geojson .gpkg .kml\n'
                          '... with columns lon, lat, ve, vn [vu] [se sn su] [site],\n'
                          'any spelling (Ve(mm/yr), sig_e ...); no header = GMT psvelo\n'
                          '(lon lat ve vn se sn corr site).')
    grp.add_argument('--vel', dest='gnss_comp', default=['los'], nargs='+',
                     choices=['ve', 'vn', 'vu', 'los'], metavar='COMP',
                     help='the GNSS component(s) to compare with the data, one or\n'
                          'several of ve vn vu los (default: los, projected into the\n'
                          'LOS of the first file). Several - --vel ve vn vu los - are\n'
                          'drawn together, each in a colour of its own and shifted\n'
                          'onto the InSAR on its own (--gnss-ref).')
    grp.add_argument('--gnss-unit', dest='gnss_unit', default='auto',
                     choices=['auto', 'mm/yr', 'cm/yr', 'm/yr'],
                     help='unit of the GNSS file (default: auto - said by the column\n'
                          'names, else mm/yr). Converted to the data\'s unit.')
    grp.add_argument('--gnss-ref', dest='gnss_ref', default='insar', choices=['insar', 'none'],
                     help='insar (default): shift the GNSS onto the InSAR, by the median\n'
                          'of InSAR - GNSS at the stations in the swath; none: as read.')
    grp.add_argument('--gnss-radius', dest='gnss_radius', type=float, default=1.0,
                     metavar='KM', help='the InSAR "at" a station: the pixels within this\n'
                                        'distance (default: %(default)s km).')
    grp.add_argument('--gnss-label', dest='gnss_label', action='store_true',
                     help='write the site names on the map.')
    grp.add_argument('--geom', dest='geom_file', metavar='PATH',
                     help='--vel los: the LOS geometry - MintPy\'s geometryGeo.h5,\n'
                          'LiCSBAS\'s E.geo, N.geo, U.geo (a GEOCml* folder), E / N / U\n'
                          'GeoTIFFs ({frame}.geo.E.tif, E.geo.tif, E.tif ...), an\n'
                          'incidence + azimuth GeoTIFF pair - or a FOLDER, searched\n'
                          'with its subfolders for one that covers the data (default:\n'
                          'found beside the data or in its LiCSBAS folder, else the\n'
                          'data\'s attributes, else nominal angles for the orbit its\n'
                          'name says: asc, desc ...).')
    grp.add_argument('--los-inc', dest='los_inc', type=float, nargs='+', metavar='DEG',
                     help='--vel los: one incidence angle for the scene instead.')
    grp.add_argument('--los-az', dest='los_az', type=float, nargs='+', metavar='DEG',
                     help='--vel los: and its azimuth angle (MintPy convention: of the\n'
                          'ground-to-satellite vector, from north, anti-clockwise;\n'
                          'Sentinel-1 about 100 ascending, -100 descending).')
    return grp


def add_save_options(group, what):
    """--save [DIR] and --outdir - the same for every command."""
    group.add_argument('--save', dest='save_dir', nargs='?', const='', metavar='DIR',
                       help='export into DIR: the folder "e" writes to, made if it\n'
                            'does not exist (default: --outdir, else the current\n'
                            f'folder). Also saves {what} at start-up.')
    group.add_argument('--outdir', dest='outdir', default='.', metavar='DIR',
                       help='the folder to export into, as --save DIR but without the\n'
                            'save at start-up; made if it does not exist (default:\n'
                            'the current folder).')
    return group


# options whose value may be left out: a file after them is an input instead
BARE_BEFORE_A_FILE = {'--save': '', '--cmap-prof': 'viridis', '--legend': 'auto'}
# ... and options of any number of numbers: moved to the end, where they take none
BARE_LIST_BEFORE_A_FILE = ('--creep-offset',)


def _bare_save(argv):
    """`--save velocity.h5`: a data file after --save is an input, not the folder
    to export into - keep --save bare then (and --cmap-prof: no colormap;
    --creep-offset: no numbers)."""
    if argv is None:
        argv = sys.argv[1:]
    argv = list(argv)
    for i, tok in enumerate(argv[:-1]):
        if tok in BARE_BEFORE_A_FILE and os.path.isfile(argv[i + 1]):
            argv[i] = f'{tok}={BARE_BEFORE_A_FILE[tok]}'
    moved = [i for i, tok in enumerate(argv[:-1])
             if tok in BARE_LIST_BEFORE_A_FILE and os.path.isfile(argv[i + 1])]
    return [tok for i, tok in enumerate(argv) if i not in moved] + [argv[i] for i in moved]


def check_colormap(name, option):
    """A colormap name matplotlib, MintPy or cmcrameri knows - else a clear error."""
    import matplotlib

    from .compat import have_mintpy

    if name in matplotlib.colormaps:
        return name
    if have_mintpy():
        try:
            from mintpy.objects.colors import ColormapExt
            ColormapExt(name)
            return name
        except Exception:        # noqa: BLE001 - a name MintPy does not know either
            pass
    try:
        import cmcrameri.cm as scm
        if hasattr(scm, name[:-2] if name.endswith('_r') else name):
            return name
    except ImportError:
        pass
    raise ValueError(f'{option}: no colormap called {name!r} (try viridis, magma, inferno, '
                     'plasma, cividis, Greys, hot_r ...)')


def check_save_options(inps):
    """What to save at start-up, and the folder every export goes to.

    --save [DIR], --outfile and --nodisplay save at start-up. --save DIR (or else
    --outdir) is the folder of every export - "e" too - and it is made now,
    if it does not exist. A bare --outfile file name goes into that folder.
    """
    save_dir = getattr(inps, 'save_dir', None)
    inps.save_fig = save_dir is not None or bool(inps.outfile) or not inps.disp_fig
    folder = inps.outdir or '.'
    if save_dir:
        if folder != '.' and os.path.normpath(folder) != os.path.normpath(save_dir):
            print(f'NOTE: --save {save_dir} and --outdir {folder}: exporting into {save_dir}')
        folder = save_dir
    folder = os.path.normpath(os.path.expandvars(os.path.expanduser(folder)))
    inps.outdir = make_folder(folder, '--save' if save_dir else '--outdir')
    if inps.outfile:
        names = []
        for name in inps.outfile:
            name = os.path.expandvars(os.path.expanduser(name))
            if not os.path.dirname(name):
                name = os.path.join(inps.outdir, name)
            make_folder(os.path.dirname(name) or '.', '--outfile')
            names.append(name)
        inps.outfile = names
    return inps


def make_folder(folder, option='--save'):
    """Make `folder` (and its parents) if it does not exist; return it."""
    if os.path.exists(folder) and not os.path.isdir(folder):
        raise ValueError(f'{option}: {folder} is a file, not a folder')
    if not os.path.isdir(folder):
        try:
            os.makedirs(folder, exist_ok=True)
        except OSError as err:
            raise ValueError(f'{option}: cannot make the folder {folder}: {err}') from err
        print(f'made the folder {os.path.abspath(folder)}')
    return folder


def add_dem_source(group):
    """--dem-source: where --dem / --topoxsect without a file download from."""
    group.add_argument('--dem-source', dest='dem_source', type=str.lower, choices=DEM_SOURCES,
                       metavar='SOURCE',
                       help='where --dem / --topoxsect without a file download from:\n'
                            'cop30 - the Copernicus GLO-30 DEM, 1" (~30 m), the\n'
                            'default; terrain - the Terrain Tiles ("Mapzen Global\n'
                            'Terrain" of QGIS), also the fallback when GLO-30 cannot\n'
                            'be had and no source is named.')


def check_map_options(inps):
    """--dem, --dem-source, --dem-zoom, --lookup and --fault - the same for swath,
    stat and model."""
    # --dem FILE stays view.py's background, exactly as before; --dem alone = download
    if is_auto(inps.dem_file):
        inps.dem_file = AUTO_DEM
    elif inps.dem_file:
        if not os.path.isfile(inps.dem_file):
            raise ValueError(f'--dem: no such file: {inps.dem_file}')
        inps.view_extra = list(inps.view_extra) + ['--dem', inps.dem_file]
    if inps.dem_zoom is not None and not 0 <= inps.dem_zoom <= 15:
        raise ValueError(f'--dem-zoom must be 0-15, got {inps.dem_zoom}')
    if inps.dem_zoom is not None and getattr(inps, 'dem_source', None) != 'terrain':
        print(f'NOTE: --dem-zoom {inps.dem_zoom} is for the Terrain Tiles (--dem-source '
              'terrain); the Copernicus GLO-30 DEM comes at its own 1" (~30 m)')
    if inps.lookup_file and not os.path.isfile(inps.lookup_file):
        raise ValueError(f'--lookup: no such file: {inps.lookup_file}')
    if inps.fault_file:
        # os.path.exists, not isfile: a FileGDB (.gdb) is a directory
        missing = [f for f in inps.fault_file if not os.path.exists(f)]
        if missing:
            raise ValueError(f'--fault: no such file: {", ".join(missing)}')
    if inps.fault_lw <= 0:
        raise ValueError(f'--fault-lw must be positive, got {inps.fault_lw}')
    geoms = getattr(inps, 'geom_file', None) or []
    for opt, path in ([('--gnss', getattr(inps, 'gnss_file', None))]
                      + [('--geom', g) for g in ([geoms] if isinstance(geoms, str) else geoms)]):
        if path and not os.path.exists(path):
            raise ValueError(f'{opt}: no such file: {path}')
    inc, az = getattr(inps, 'los_inc', None), getattr(inps, 'los_az', None)
    if (inc is None) != (az is None):
        raise ValueError('--los-inc and --los-az go together')
    if inc is not None and len(inc) != len(az):
        raise ValueError(f'--los-inc has {len(inc)} value(s), --los-az {len(az)}')
    if getattr(inps, 'data_unit', None):
        from .compat import normalise_unit
        inps.data_unit = normalise_unit(inps.data_unit)
    check_data_scale(inps)
    if getattr(inps, 'gnss_radius', 1.0) <= 0:
        raise ValueError(f'--gnss-radius must be positive, got {inps.gnss_radius}')
    return inps


def check_data_scale(inps, files=None):
    """--data-scale: finite, not 0, one per input file at most."""
    scale = getattr(inps, 'data_scale', None)
    if not scale:
        return inps
    files = files if files is not None else (getattr(inps, 'file', None) or [])
    for value in scale:
        if not value or value != value or abs(value) == float('inf'):
            raise ValueError(f'--data-scale: {value} is no factor')
    if len(scale) > 1 and files and len(scale) != len(files):
        raise ValueError(f'--data-scale has {len(scale)} values for {len(files)} file(s): '
                         'give one per file, or one for all')
    return inps


def read_file_info(inps):
    """The input file list, its metadata and the dataset to show."""
    from .compat import coordinate, get_file_list, read_attribute

    inps.file = get_file_list(inps.file)
    if not inps.file:
        raise ValueError('no such input file')
    labels = getattr(inps, 'labels', None)
    if labels and len(labels) != len(inps.file):
        raise ValueError(f'--labels: {len(labels)} name(s) for {len(inps.file)} file(s)')
    inps.atr = read_attribute(inps.file[0])
    inps.coord = coordinate(inps.atr)
    inps.num_file = len(inps.file)
    if not inps.dset:
        inps.dset = default_dset(inps.file[0], inps.atr, print_msg=inps.print_msg)
    return inps


def cmd_line_parse(iargs=None):
    """Parse the `swath` arguments; unknown ones are forwarded to view.py."""
    iargs = _argv(iargs, 'swath')
    parser = create_swath_parser()
    inps, view_extra = parser.parse_known_args(args=_bare_save(iargs))
    inps.argv = list(iargs)
    inps.view_extra = view_extra
    # --vel: one GNSS component or several, each once, in the order given
    comps = inps.gnss_comp if isinstance(inps.gnss_comp, list) else [inps.gnss_comp]
    inps.gnss_comps = list(dict.fromkeys(comps))
    inps.gnss_comp = inps.gnss_comps[0]

    check_save_options(inps)
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
        # --topoxsect alone: the --dem file if there is one, else the download
        inps.topoxsect = inps.dem_file if inps.dem_file else AUTO_DEM
    elif inps.topoxsect and not os.path.isfile(inps.topoxsect):
        raise ValueError(f'--topoxsect: no such file: {inps.topoxsect}')

    read_file_info(inps)
    collect_profiles(inps)

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
The creep rate along a fault, from short profiles across a line drawn on it -
as Hussain et al. (2016) measured the creep of the North Anatolian Fault.

Draw the line on the map (left-click its points, Enter to finish), or give it
with --line-lalo / --line-yx / --line-file - ON the fault trace (--fault shows
it, --line-file takes it): the two sides are compared at the line. Every
--spacing km along it (centred on the line) a profile crosses it at right
angles: --length km (default 5) either side of the line, the pixels within
--width km (default 2.5) either side of the profile projected onto it. Block 1
is the LEFT side of the line, walking from its first point to its last; block
2 the RIGHT.

A straight line through the pixels of each block, and the creep rate is their
offset at the line - block 2 minus block 1 - one linear system (their eq. 2):

    | v1 |   | x1  0  1  0 |  | slope 1 |
    |    | = |             |  | slope 2 |      x: the distance from the line
    | v2 |   | 0  x2  1  1 |  | b1      |      along the profile, + towards
                              | creep   |      block 2

Its uncertainty: a percentile bootstrap that resamples each block on its own,
n1 of its n1 pixels and n2 of n2 (--boot draws), as Hussain et al. - the least-
squares one (the BLUE) is too narrow when the errors are large or correlated.

The figure: the LOS creep rate along the line, a series per file, 1-sigma
error bars (descending tracks solid, ascending dashed). Two or more LOS
geometries - ascending and descending, found beside the data or given with
--geom or --los-inc / --los-az - also give it split into fault-parallel (+
right-lateral) and vertical (+ block 2 up), assuming no motion across the
fault. Click a point of it: its profile, the pixels with both fits (BLUE red,
bootstrap green) and the PDFs of the creep rate (their Fig. S3).
"""

STAT_EXAMPLE = """\
examples:
  # draw the line along a fault shown on the map: a profile every 5 km
  viewinsar stat velocity.h5 --fault faults.shp
  viewinsar stat geo_velocity.h5 --spacing 2 --length 4 --width 1.5 --dem \\
      --colormap RdBu_r --vlim -1 1

  # the line given up front: lat/lon or y/x pairs, or the longest line of a
  # vector file (.shp .geojson .gpkg .kml .gmt ..., e.g. drawn in QGIS)
  viewinsar stat geo_velocity.h5 --save --nodisplay \\
      --line-lalo 40.70 29.90 40.73 30.40 40.70 30.70
  viewinsar stat velocity.h5 --line-yx 3670 744 3300 720 2967 696
  viewinsar stat velocity.h5 --line-file fault_trace.geojson --svg --save

  # ascending and descending: the LOS creep of each, and split into
  # fault-parallel and vertical
  viewinsar stat T036D/geo_velocity.h5 T131A/geo_velocity.h5 --labels T036D T131A \\
      --line-file fault_trace.geojson
  viewinsar stat desc.tif asc.tif --los-inc 39 34 --los-az -100 100 --data-unit mm/yr \\
      --line-file fault_trace.geojson

  # your own names and colours for the two blocks, more bootstrap draws
  viewinsar stat velocity.h5 --names North South --colors tab:green tab:purple \\
      --boot 5000 --line-file fault_trace.geojson

mouse / keyboard:
  left click (map)  : add a point to the line
  Enter             : finish the line -> the profiles and the creep rate
  "u" or Backspace  : remove the last point
  "n"               : start a new line (so does a click after Enter)
  left click on a point of the creep rate: its profile, in a window of its own -
                      there left / right arrow: the point before / after,
                      up / down arrow: the other files
  "e"               : export (see "output")

output (written by "e" or --save, into the --save folder - made if it does not
exist - or --outdir, else the current folder), named after the line's first and
last point, e.g. stat_Y3670X744_Y2967X696_sp5km_L5km_w2.5km:
  <name>.png           the creep rate along the line (.svg with --svg)
  <name>_map.png       the map: the line, its points and their profiles
  <name>_creep.csv     a row per point and file: dist_km along the line, row /
                       col, lat / lon, strike, the creep rate (the bootstrap
                       mean, its 1 sigma, median, 68 and 95 % ranges), the
                       least-squares (BLUE) one and its sigma, both slopes,
                       the pixels of each block
  <name>_split.csv     two or more LOS geometries: the creep split into
                       fault-parallel (+ right-lateral) and vertical (+ block
                       2 up) at every point
  <name>_values.csv    every pixel used: file, point, block, x_km (from the
                       line, + towards block 2), y_km, value, row, col, lat/lon
  <name>_profiles.png  every profile of the first file (_profiles_2.png ...)
  <name>_profile7.png  the profile on screen (point 7), when its window is open
  <name>_info.txt      the line, the parameters, which side is which block

map options: anything not listed above goes to view.py, as for `viewinsar swath`
(--colormap --vlim --unit --mask --sub-lat --sub-lon ...).
"""


def create_stat_parser(prog=f'{PROG} stat'):
    """The argument parser of the `stat` sub-command."""
    parser = ViewinsarParser(
        prog=prog, description=STAT_DESCRIPTION, epilog=STAT_EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)
    parser.add_argument('file', nargs='+',
                        help='input file(s): velocity.h5, a GeoTIFF, a time-series epoch\n'
                             '... geocoded or in radar coordinates. Several files -\n'
                             'tracks, frames - each get the same profiles, a series each.')
    parser.add_argument('--labels', dest='labels', nargs='+', metavar='NAME',
                        help='names of the input files in the legends (default: the\n'
                             'file names).')
    parser.add_argument('--dset', dest='dset',
                        help='dataset to read (default: the one matching FILE_TYPE;\n'
                             'see `viewinsar list FILE`).')
    parser.add_argument('--data-unit', dest='data_unit', metavar='UNIT',
                        help='the unit the input files are really in: mm/yr, cm/yr,\n'
                             'm/yr, mm ... - for a file that does not say it (a GeoTIFF),\n'
                             'or says a wrong one (e.g. mm/yr labelled m/year). Default:\n'
                             'the file\'s own; mm/yr for a LiCSBAS velocity (vel*, vstd*)\n'
                             'without a unit, inside its LiCSBAS folder; for any other\n'
                             'velocity without one (model, --gnss, or a vel* name), m/yr\n'
                             'when 99%% of its values are below 1, else mm/yr.')
    parser.add_argument('--data-scale', dest='data_scale', type=float, nargs='+',
                        metavar='FACTOR',
                        help='multiply the data by FACTOR before anything - the map, the\n'
                             'profiles, the fit, the statistics - one per input file (or\n'
                             'one for all): 0.1 for a product 10x too large. For every\n'
                             'file of a mission, a rule in your viewinsar.ini does it -\n'
                             '~/.config/viewinsar/viewinsar.ini: [data-scale]\n'
                             'ALOS-2 = 0.1 (--data-scale 1: not this time).')
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

    prof = parser.add_argument_group('Profiles', 'across the line, where the creep is '
                                     'measured (Hussain et al. 2016)')
    prof.add_argument('--spacing', dest='spacing', type=float, default=5.0,
                      metavar='KM',
                      help='distance between the points along the line - a profile\n'
                           'at each - in km (default: %(default)s). Centred on the\n'
                           'line.')
    prof.add_argument('--length', dest='length_km', type=float, default=5.0, metavar='KM',
                      help='how far each profile reaches either side of the line, in\n'
                           'km (default: %(default)s).')
    prof.add_argument('--width', dest='width_km', type=float, default=2.5, metavar='KM',
                      help='the pixels within this distance either side of a profile\n'
                           'are projected onto it, in km (default: %(default)s).')
    prof.add_argument('--boot', dest='n_boot', type=int, default=2000, metavar='NUM',
                      help='bootstrap draws for the uncertainty of each creep rate\n'
                           '(default: %(default)s).')
    prof.add_argument('--seed', dest='seed', type=int, default=0, metavar='NUM',
                      help='seed of the bootstrap: the same NUM, the same draws\n'
                           '(default: %(default)s).')
    prof.add_argument('--min-pixels', dest='min_pixels', type=int, default=10, metavar='NUM',
                      help='fewer pixels than this on a side of a profile: no creep\n'
                           'rate there (default: %(default)s).')
    prof.add_argument('--keep-zero', dest='keep_zero', action='store_true',
                      help='keep exactly-zero pixels (dropped by default).')

    geo = parser.add_argument_group('LOS geometry', 'the split into fault-parallel and '
                                    'vertical (two or more files)')
    geo.add_argument('--geom', dest='geom_file', nargs='+', metavar='PATH',
                     help='the LOS geometry, one per input file in order: MintPy\'s\n'
                          'geometryGeo.h5, LiCSBAS\'s E.geo, N.geo, U.geo (a GEOCml*\n'
                          'folder), E / N / U GeoTIFFs ({frame}.geo.E.tif, E.geo.tif,\n'
                          'E.tif ...), an incidence + azimuth GeoTIFF pair - or a\n'
                          'FOLDER, searched with its subfolders for the one that\n'
                          'covers the file (default: found beside each file or in\n'
                          'its LiCSBAS folder, else its attributes; the nominal\n'
                          'angles of a name - asc, desc ... - are not used).')
    geo.add_argument('--los-inc', dest='los_inc', type=float, nargs='+', metavar='DEG',
                     help='one incidence angle for all the pixels of a file instead,\n'
                          'one per file (or one for all).')
    geo.add_argument('--los-az', dest='los_az', type=float, nargs='+', metavar='DEG',
                     help='and the azimuth angle (MintPy convention: of the ground-\n'
                          'to-satellite vector, from north, anti-clockwise; Sentinel-1\n'
                          'about 100 ascending, -100 descending).')

    look = parser.add_argument_group('Plot', 'the creep rate and the profiles')
    look.add_argument('--names', dest='block_names', nargs=2, metavar=('NAME1', 'NAME2'),
                      default=['Block 1', 'Block 2'],
                      help='names of block 1 (left) and block 2 (right)\n'
                           '(default: "Block 1" "Block 2").')
    look.add_argument('--colors', dest='block_colors', nargs=2, metavar=('COLOR1', 'COLOR2'),
                      default=['#0072B2', '#D55E00'],
                      help='colours of block 1 and block 2 on the map and in the\n'
                           'profiles, any matplotlib colour (default: blue #0072B2\n'
                           'and vermilion #D55E00).')
    look.add_argument('--ylim', dest='ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'),
                      help='range of the creep-rate axis (default: auto).')

    mapg = parser.add_argument_group('Map', 'the guides drawn on the map')
    mapg.add_argument('--dem', dest='dem_file', nargs='?', const=AUTO_DEM, metavar='DEM',
                      help='DEM for the shaded-relief background, as in view.py.\n'
                           'Without a file: downloaded for the map (Copernicus\n'
                           'GLO-30), or the height of geometryRadar.h5 on radar-coded\n'
                           'data.')
    add_dem_source(mapg)
    mapg.add_argument('--dem-zoom', dest='dem_zoom', type=int, metavar='Z',
                      help='zoom of the Terrain Tiles (--dem-source terrain), 0-15\n'
                           '(default: auto).')
    mapg.add_argument('--lookup', dest='lookup_file', metavar='FILE',
                      help='radar-coded data: MintPy geometryRadar.h5 (default:\n'
                           'found in inputs/, ./ or ../inputs/ as MintPy does).')
    mapg.add_argument('--fault', dest='fault_file', nargs='+', metavar='FILE',
                      help='fault traces to draw on the map, to follow with the line;\n'
                           'any vector format and CRS, as in `viewinsar swath`.')
    mapg.add_argument('--fault-color', dest='fault_color', default='k', metavar='COLOR',
                      help='colour of the fault traces (default: %(default)s).')
    mapg.add_argument('--fault-lw', dest='fault_lw', type=float, default=1.2, metavar='NUM',
                      help='line width of the fault traces (default: %(default)s).')

    out = parser.add_argument_group('Save/Output', 'figures and data - see "output" below')
    out.add_argument('--outfile', type=str, nargs='*',
                     help='file name of the saved MAP figure.')
    add_save_options(out, 'everything, when the line is given up front\n(--line-yx / --line-lalo / --line-file)')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not open any window (implies --save).')
    out.add_argument('--svg', dest='svg', action='store_true',
                     help='write the figures as vector .svg instead of .png.')
    out.add_argument('--figext', dest='fig_ext', default='.png', metavar='EXT',
                     choices=['.png', '.pdf', '.svg', '.eps', '.jpg'],
                     help='figure format: .png .pdf .svg .eps .jpg (default: %(default)s).')
    out.add_argument('--dpi', dest='fig_dpi', type=int, default=300, metavar='NUM',
                     help='resolution of raster figures (default: %(default)s).')
    out.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='size of the creep-rate window, in inches (default: 9 4.8;\n'
                          '9 7.4 with the split).')
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

    iargs = _argv(iargs, 'stat')
    parser = create_stat_parser()
    inps, view_extra = parser.parse_known_args(args=_bare_save(iargs))
    inps.argv = list(iargs)
    inps.view_extra = view_extra

    check_save_options(inps)
    if inps.svg:
        inps.fig_ext = '.svg'
    for opt, value in (('--spacing', inps.spacing), ('--length', inps.length_km),
                       ('--width', inps.width_km)):
        if not value > 0:
            raise ValueError(f'{opt} must be positive, got {value}')
    if inps.n_boot < 100:
        raise ValueError(f'--boot: at least 100 draws, got {inps.n_boot}')
    if inps.min_pixels < 3:
        raise ValueError(f'--min-pixels: at least 3 (a line and its scatter), got '
                         f'{inps.min_pixels}')
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
    for opt, values in (('--geom', inps.geom_file), ('--los-inc / --los-az', inps.los_inc)):
        if values and len(values) > len(inps.file):
            raise ValueError(f'{opt}: {len(values)} for {len(inps.file)} input file(s)')

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
#                                  model                                     #
##############################################################################
MODEL_DESCRIPTION = """\
Fit an elastic dislocation model to a profile across a fault - of the type
named by the first word: strike (the default), reverse (or thrust), normal.

  strike   a strike-slip fault: the slip rate S below the locking depth d1
           (Savage & Burford 1973) - optionally with shallow creep C from
           the surface down to d2, and a rotation term - as Hussain et al.
           (2016, 2018) did along the North Anatolian Fault, fitted to the
           fault-parallel velocity (--comp par, ve, vn, los - one for all):
             v(x) = S/pi atan(x/d1) + C [H(x) - 1/2 - 1/pi atan(x/d2)] + a
             viewinsar model strike geo_ve.h5 --comp ve ...
  reverse  a reverse (thrust) or a normal fault dipping at dip (--dip DEG,
  normal   or a range: --dip 10 45; default 5 85) towards A' or A
           (--dip-toward end|start), locked from the surface down to d1,
           slipping at S along its dip below it - up the dip (reverse) or
           down it (normal): the back-slip model of Savage (1983), as Segall
           (2010, eqs. 3.71-3.72) writes it. Fitted to the vertical velocity
           (vu), the motion across the fault (ve / vn: the east / north
           velocity over the east / north of the fault normal) or a line of
           sight as it is (los) - one --comp per file, a panel each, the
           model in each; the far-field rates follow - shortening / extension
           S cos(dip), uplift / subsidence S sin(dip):
             viewinsar model reverse geo_vu.h5 --comp vu ...
             viewinsar model normal geo_ve.h5 geo_vu.h5 --comp ve vu ...

Draw the profile ACROSS the fault on the map (two clicks), or give it - or
several, A-A', B-B', C-C' ... - up front with --start-lalo / --end-lalo,
--profile (repeat it) or --line-file (one section per line); each section is
fitted on its own, with windows of its own. Every pixel of the swath
(+/-WIDTH km) becomes a point: velocity against its distance from the
fault. The points are binned (--bin km, the median of each bin, as Ou et
al. 2022) and the bins fitted with the affine-invariant ensemble MCMC
sampler of Goodman & Weare (2010) - the emcee sampler Hussain et al. used -
under uniform priors, with a noise covariance that decays with distance
(estimated from the residuals, or --noise). Several files (tracks, frames)
are fitted together, each with an offset of its own; GNSS (--gnss) too, with
theirs.

One window shows the points (no density: the pixels themselves, red), the
binned medians, the GNSS (blue), the maximum a posteriori (MAP) model (bold
dashes) and its 95 % and 68 % ranges (grey, --bands), with the parameters
and their 95 % ranges written on it, and a legend - inside the plot, or
beside it when the plot is too full for it (--legend; "l" on the profile:
inside, outside, off, in turn); a second one the posterior distributions, as
Hussain et al. (2016, central NAF, Fig. 9): every pair of parameters as a 2D
histogram (jet), the marginals as blue histograms along the bottom, the MAP
in red.
With --cmap-prof CMAP the pixels are a pixel DENSITY in that colormap
instead, as in swath (--nlook fills it in), under red ranges, with bigger
medians and the GNSS red.

Creep, as Hussain et al. (2016) constrained it, two ways: --creep fits the
creep rate C and its depth d2 in the model (written on the profile, the MAP
model without the creep dotted); --creep-offset measures the rate straight
from the data - the offset at the fault of two straight lines, one through
the pixels on each side within a few km of it, its 95 % range by bootstrap.
The LOS offsets of an ascending and a descending track are split into
fault-parallel and vertical, so a vertical step at the fault (subsidence of
one side) is not taken for creep.

The fault: where the profile crosses a --fault trace, with the strike of the
trace there - or --x0 and --strike. --comp says what the data are (strike-
slip; for reverse and normal faults see above and --comp below): par
(fault-parallel velocity), ve / vn (east / north velocity), los (line-of-
sight velocity), vu (vertical velocity) or as-is (fitted as they are); ve,
vn and los are turned into fault-parallel velocity assuming purely
fault-parallel motion (no vertical, no fault-normal motion). Then S > 0 is
right-lateral, and the sense of slip is written out. By default (auto) a LOS
velocity whose line of sight is known becomes fault-parallel velocity by
itself. The y axis says what is plotted: "Fault-parallel velocity (mm/yr)",
or "LOS velocity (mm/yr)" (the LOS fitted as it is); for a reverse or normal
fault "Fault vertical velocity", "Fault-normal velocity" or "LOS velocity".

The line of sight of every pixel is found next to the data: MintPy's
geometryGeo.h5 / geometryRadar.h5, LiCSBAS's GEOCml*/E.geo, N.geo, U.geo (for
a file anywhere under TS_GEOCml*/), or LiCSAR's {frame}.geo.E.tif, .N.tif,
.U.tif; --geom names it, --los-inc / --los-az replace it. A LiCSBAS velocity
(vel*, also as a GeoTIFF) is read in mm/yr; --data-unit gives the unit of any
other file that does not say its own.
"""

MODEL_EXAMPLE = """\
examples:
  # click a profile across the fault: the trace gives the fault position and
  # its strike, the LOS velocity becomes fault-parallel velocity
  viewinsar model geo_velocity.h5 --comp los --fault faults.shp --width 25

  # up front, saved, nothing shown; with shallow creep (C, d2)
  viewinsar model geo_velocity.h5 --comp los --fault faults.shp --creep \\
      --start-lalo 40.20 33.10 --end-lalo 41.30 33.40 --width 25 --nodisplay

  # three profiles along the fault, A-A', B-B', C-C', each fitted, saved in fits/
  viewinsar model geo_velocity.h5 --comp los --fault faults.shp --width 25 \\
      --profile 40.20 33.10 41.30 33.40 --profile 40.30 33.90 41.40 34.20 \\
      --profile 40.40 34.70 41.50 35.00 --nodisplay --save fits
  viewinsar model geo_velocity.h5 --comp los --fault faults.shp \\
      --line-file profiles.gmt --nodisplay --save fits    # one per segment

  # LiCSBAS, as it is: cum_filt.h5 (its velocity, vel), results/vel.filt.mskd or
  # a GeoTIFF of it - the grid, the unit (mm/yr), results/mask and the LOS of
  # every pixel (GEOCml*/E.geo, N.geo, U.geo) are found by themselves
  viewinsar model TS_GEOCml10GACOS/cum_filt.h5 --fault faults.shp --gnss gnss.csv \\
      --width 20 --strain
  viewinsar model TS_GEOCml10GACOS/results/vel.filt.mskd --fault faults.shp
  # ... the same file copied elsewhere: say where its geometry and unit are
  viewinsar model vel_asc.tif --geom ~/licsbas/105A/GEOCml10GACOS \\
      --data-unit mm/yr --fault faults.shp

  # east velocity (asc + desc decomposed) and GNSS, fitted together
  viewinsar model geo_ve.h5 --comp ve --gnss gnss.csv --fault faults.shp --strain

  # two tracks at once, each with its own offset
  viewinsar model T014A/geo_velocity.h5 T123D/geo_velocity.h5 --comp los \\
      --fault faults.shp --labels T014A T123D

  # the pixels as a pixel density (as swath), filled in for ml10 data; the
  # ranges red, bigger medians
  viewinsar model TS_GEOCml10GACOS/cum_filt.h5 --fault faults.shp --width 20 \\
      --cmap-prof magma --nlook 10

  # as Hussain et al. (2016): rotation of the block south of the fault,
  # a Gaussian prior on the locking depth, a given noise covariance
  viewinsar model geo_vpar.h5 --comp par --strike 100 --rot 0.0215 \\
      --d1-prior 17 5 --noise 6.35 35.8

  # a reverse fault, from the vertical velocity (asc + desc decomposed): the
  # dip free between 10 and 45 degrees, dipping towards A' (the end of the line)
  viewinsar model reverse geo_vu.h5 --comp vu --fault faults.shp --dip 10 45

  # a normal fault from the east and the vertical velocity together, a panel
  # each, dipping towards A; the GNSS across the fault (--vel perp, the default)
  viewinsar model normal geo_ve.h5 geo_vu.h5 --comp ve vu --fault faults.shp \\
      --dip-toward start --gnss gnss.csv

  # a thrust (= reverse) seen by two tracks, each LOS as it is, the dip fixed
  viewinsar model thrust T014A/geo_velocity.h5 T123D/geo_velocity.h5 --comp los \\
      --fault faults.shp --labels T014A T123D --dip 25

  # the creep: its rate and depth in the model, and its rate from the offset of
  # two lines +/-5 km either side of the fault - ascending and descending
  # together, the offsets split into fault-parallel and vertical
  viewinsar model T014A/geo_velocity.h5 T123D/geo_velocity.h5 --comp los \\
      --fault faults.shp --labels T014A T123D --creep --creep-offset
  viewinsar model geo_velocity.h5 --comp los --fault faults.shp --creep \\
      --creep-offset 3 2      # lines 3 km long, pixels within 2 km of the line

mouse / keyboard: as `viewinsar swath` - two left clicks draw a section (and
fit it), "=" starts another, "n" redoes the current one, "e" exports; "l" on a
profile moves its legend: inside, outside, off, in turn.

output (written by "e" or --save, into the --save folder - made if it does not
exist - or --outdir, else the current folder), for each section X, named after
its line, e.g. model_A_Y470X150_Y130X700_w25km:
  <name>.png             the fit: points, bins, GNSS, MAP model, 95 / 68 % ranges
  <name>_posterior.png   the posterior distributions, the MAP in red
  <name>_params.csv      every parameter: MAP, median, the 95 % range (2.5 and
                         97.5 percentiles: lower_95, upper_95),
                         and the shear strain rate at the fault |S| / (pi d1)
                         - reverse / normal: h_rate and v_rate, S cos(dip) and
                         S sin(dip) (shortening / extension, uplift /
                         subsidence);
                         --creep-offset: C_offset (C_offset_1, _2 ... per file;
                         C_offset_par and U_offset, split) - the offset of the
                         lines, the median and 95 % range of the bootstrap
  <name>_params.txt      the same, readable, with the fault position, strike,
                         sense of slip, noise, priors and sampler diagnostics
  <name>_samples.csv     posterior samples, for plots of your own
  <name>_bins.csv        the binned medians and the model at each
  <name>_points.csv      every pixel: distance along the profile and from the
                         fault, its value and the velocity fitted
  <name>_gnss.csv        the GNSS in the swath (--gnss)
  model_map_AB.png       the map, with every section and the fault positions
  model_sections_AB.csv  two sections or more: one row each - its ends, where
                         it crosses the fault, the strike, the sense of slip,
                         S, d1 (C, d2, x0) and the strain with their 95 % ranges
                         (reverse / normal: S, d1, dip, x0, h_rate, v_rate)

map options: anything not listed above goes to view.py, as for `viewinsar swath`
(--colormap --vlim --unit --mask --sub-lat --sub-lon ...).
"""


# viewinsar model strike|reverse|normal: the fault type, the first word
FAULT_TYPES = {'strike': 'strike', 'strike-slip': 'strike', 'strikeslip': 'strike',
               'reverse': 'reverse', 'thrust': 'reverse', 'normal': 'normal'}
MODEL_COMPS = ('auto', 'as-is', 'par', 've', 'vn', 'los', 'vu')
DIP_SLIP_COMPS = ('auto', 'vu', 've', 'vn', 'los')     # reverse, normal
DIP_SLIP_GNSS = ('perp', 'vu')
# options of one kind of fault only: (option, what it is for)
STRIKE_ONLY = {'--creep': 'shallow creep', '--creep-offset': 'the creep offset',
               '--rot': 'the rotation term', '--rot-side': 'the rotation term',
               '--strain': 'the shear strain rate', '--c-max': 'the creep rate',
               '--d2-max': 'the creep depth'}
DIP_SLIP_ONLY = {'--dip': 'the dip', '--dip-toward': 'the dip direction'}


def split_fault_type(argv):
    """`viewinsar model reverse FILE ...`: the fault type (strike, reverse,
    normal - thrust is reverse), and the arguments after it; strike when the
    first word is none of them."""
    argv = list(argv)
    if argv and str(argv[0]).lower() in FAULT_TYPES:
        return FAULT_TYPES[str(argv[0]).lower()], argv[1:], True
    return 'strike', argv, False


def _glue_values(argv, option, accept, most=None):
    """`--comp ve vu FILE ...`: the values after `option` that `accept` takes
    (at most `most`), glued to it - `--comp=ve,vu` - so the files after them
    stay inputs."""
    out, i = [], 0
    while i < len(argv):
        tok = argv[i]
        if tok == option:
            j, values = i + 1, []
            while j < len(argv) and (most is None or len(values) < most) and accept(argv[j]):
                values.append(argv[j])
                j += 1
            if values:
                out.append(f'{option}=' + ','.join(values))
                i = j
                continue
        out.append(tok)
        i += 1
    return out


def _is_number(text):
    try:
        float(text)
    except (TypeError, ValueError):
        return False
    return True


def model_argv(argv):
    """The arguments of `model` after its fault type, as its parser takes
    them: the values of --comp and --dip glued to them (--comp=ve,vu), so the
    files after them stay inputs; --save, --cmap-prof ... bare before a file."""
    argv = _glue_values(argv, '--comp', lambda tok: str(tok).lower() in MODEL_COMPS)
    argv = _glue_values(argv, '--dip', _is_number, most=2)
    return _bare_save(argv)


def given_options(argv, options):
    """Those of `options` that the command line gives (--opt or --opt=...)."""
    return [opt for opt in options
            if any(tok == opt or str(tok).startswith(opt + '=') for tok in argv)]


def create_model_parser(prog=f'{PROG} model [strike|reverse|normal]'):
    """The argument parser of the `model` sub-command."""
    parser = ViewinsarParser(
        prog=prog, description=MODEL_DESCRIPTION, epilog=MODEL_EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)
    parser.add_argument('file', nargs='+',
                        help='input file(s): velocity.h5, a GeoTIFF ... Several files\n'
                             '(tracks, frames) are fitted together, each with its own\n'
                             'offset.')
    parser.add_argument('--labels', dest='labels', nargs='+', metavar='NAME',
                        help='names of the input files in the legends (default: the\n'
                             'file names).')
    parser.add_argument('--dset', dest='dset',
                        help='dataset to read (default: the one matching FILE_TYPE;\n'
                             'see `viewinsar list FILE`).')
    parser.add_argument('--data-unit', dest='data_unit', metavar='UNIT',
                        help='the unit the input files are really in: mm/yr, cm/yr,\n'
                             'm/yr, mm ... - for a file that does not say it (a GeoTIFF),\n'
                             'or says a wrong one (e.g. mm/yr labelled m/year). Default:\n'
                             'the file\'s own; mm/yr for a LiCSBAS velocity (vel*, vstd*)\n'
                             'without a unit, inside its LiCSBAS folder; for any other\n'
                             'velocity without one (model, --gnss, or a vel* name), m/yr\n'
                             'when 99%% of its values are below 1, else mm/yr.')
    parser.add_argument('--data-scale', dest='data_scale', type=float, nargs='+',
                        metavar='FACTOR',
                        help='multiply the data by FACTOR before anything - the map, the\n'
                             'profiles, the fit, the statistics - one per input file (or\n'
                             'one for all): 0.1 for a product 10x too large. For every\n'
                             'file of a mission, a rule in your viewinsar.ini does it -\n'
                             '~/.config/viewinsar/viewinsar.ini: [data-scale]\n'
                             'ALOS-2 = 0.1 (--data-scale 1: not this time).')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='disable the verbose message printing.')

    line = parser.add_argument_group(
        'Profile location',
        "sections given up front, each fitted on its own: A-A' from\n"
        "--start-* / --end-*, then one per --line-file line, then one per\n"
        '--profile')
    line.add_argument('--start-yx', '--yx0', dest='start_yx', metavar=('Y0', 'X0'), type=int,
                      nargs=2, help='start point, in pixel number [row, col]')
    line.add_argument('--end-yx', '--yx1', dest='end_yx', metavar=('Y1', 'X1'), type=int,
                      nargs=2, help='end   point, in pixel number [row, col]')
    line.add_argument('--start-lalo', '--lalo0', dest='start_lalo', metavar=('LAT0', 'LON0'),
                      type=float, nargs=2,
                      help='start point in [lat, lon] (or [northing, easting] on a\n'
                           'projected file).')
    line.add_argument('--end-lalo', '--lalo1', dest='end_lalo', metavar=('LAT1', 'LON1'),
                      type=float, nargs=2, help='end   point, as --start-lalo.')
    line.add_argument('--line-file', dest='lola_file', metavar='FILE',
                      help='the sections as "lon lat" lines, GMT format, one per\n'
                           'segment (split by ">"), from its first to its last point\n'
                           '- or any vector file (.geojson .kml .shp .gpkg), one\n'
                           'per line.')
    add_profile_options(line)

    swath = parser.add_argument_group('Swath', 'the pixels that go into the fit')
    swath.add_argument('--width', dest='width', type=float, default=20.0,
                       metavar='KM',
                       help='HALF width of the swath, in km (default: %(default)s;\n'
                            'Hussain et al. took +/-20-25 km).')
    swath.add_argument('--bin', dest='bin_km', type=float, default=1.0, metavar='KM',
                       help='bin size along the profile, in km (default: %(default)s):\n'
                            'the median of each bin is fitted.')
    swath.add_argument('--keep-zero', dest='keep_zero', action='store_true',
                       help='keep exactly-zero pixels (dropped by default).')
    swath.add_argument('--decimate', dest='decimate', type=int, default=1, metavar='NUM',
                       help='use every NUMth pixel only (default: %(default)s).')
    swath.add_argument('--nlook', dest='nlook', type=int, default=1, metavar='NUM',
                       help='oversample a coarse raster NUM times per axis for the\n'
                            'picture - the points, or the --cmap-prof density, as in\n'
                            'swath (default: %(default)s). The fit always takes the real\n'
                            'pixels: interpolated samples are no new data.')

    flt = parser.add_argument_group('Fault', 'where the fault is, and its strike')
    flt.add_argument('--fault', dest='fault_file', nargs='+', metavar='FILE',
                     help='fault trace(s), any vector format and CRS (as swath). The\n'
                          'crossing nearest the middle of the profile is the fault\n'
                          'position, the trace there gives the strike.')
    flt.add_argument('--x0', dest='x0_km', type=float, metavar='KM',
                     help='the fault position: KM along the profile from its start\n'
                          '(default: the --fault crossing, else the middle).')
    flt.add_argument('--strike', dest='strike', type=float, metavar='DEG',
                     help='fault strike, degrees clockwise from north (default:\n'
                          'from the trace, else perpendicular to the profile). With\n'
                          '--comp par, the direction the data are positive along.')
    flt.add_argument('--fault-color', dest='fault_color', default='k', metavar='COLOR',
                     help='colour of the traces (default: %(default)s).')
    flt.add_argument('--fault-lw', dest='fault_lw', type=float, default=1.2, metavar='NUM',
                     help='line width of the traces on the map (default: %(default)s).')

    data = parser.add_argument_group('Data', 'what the data are')
    data.add_argument('--comp', dest='comp', default='auto', nargs='+', metavar='COMP',
                      help='what the data are - one for all the files, or one per\n'
                           'file (reverse / normal: --comp ve vu ve.h5 vu.h5).\n'
                           'strike (one for all): auto (default) - los when the LOS\n'
                           'geometry of every file is known (found, --geom,\n'
                           '--los-inc / --los-az; not the nominal angles of a name)\n'
                           'and no file is named as a component (hz, east, ve ...),\n'
                           'else as-is, and the GNSS --vel los (or the component of\n'
                           'the name); as-is: fitted as they are (LOS data: rates "in\n'
                           'the LOS"); par: fault-parallel velocity; ve / vn: east /\n'
                           'north velocity, divided by sin / cos(strike); los: LOS\n'
                           'velocity, divided by the LOS . strike direction, pixel\n'
                           'by pixel. ve, vn and los assume fault-parallel motion\n'
                           'only; the y axis says "Fault-parallel velocity (mm/yr)",\n'
                           '"LOS velocity (mm/yr)" (the LOS as it is).\n'
                           'reverse / normal: vu - the vertical velocity ("Fault\n'
                           'vertical velocity"); ve / vn - the east / north velocity,\n'
                           'divided by the east / north of the fault normal: the\n'
                           'motion across the fault, assuming none along it ("Fault-\n'
                           'normal velocity", + towards the end of the profile); los\n'
                           '- the LOS as it is, the model seen through the line of\n'
                           'sight of the file; auto (default) - by the file name (up,\n'
                           'vu, vertical; east, ve; north, vn), else los. A panel per\n'
                           'component, the model in each.')
    data.add_argument('--geom', dest='geom_file', nargs='+', metavar='PATH',
                      help='the LOS geometry, one per input file in order: MintPy\'s\n'
                           'geometryGeo.h5, LiCSBAS\'s E.geo, N.geo, U.geo (a GEOCml*\n'
                           'folder), E / N / U GeoTIFFs ({frame}.geo.E.tif, E.geo.tif,\n'
                           'E.tif ...), an incidence + azimuth GeoTIFF pair - or a\n'
                           'FOLDER, searched with its subfolders for the one that\n'
                           'covers the file (default: found beside each file or in\n'
                           'its LiCSBAS folder, else its attributes, else nominal\n'
                           'angles for the orbit its name says - asc, desc ... - for\n'
                           'the GNSS only: --comp auto then fits the LOS as it is).')
    data.add_argument('--los-inc', dest='los_inc', type=float, nargs='+', metavar='DEG',
                      help='one incidence angle for all the pixels of a file instead,\n'
                           'one per file (or one for all).')
    data.add_argument('--los-az', dest='los_az', type=float, nargs='+', metavar='DEG',
                      help='and the azimuth angle (MintPy convention: of the ground-\n'
                           'to-satellite vector, from north, anti-clockwise; Sentinel-1\n'
                           'about 100 ascending, -100 descending).')

    mod = parser.add_argument_group('Model', 'the parameters and their uniform priors')
    mod.add_argument('--creep', dest='creep', action='store_true',
                     help='add shallow creep: C (rate) from the surface down to d2\n'
                          '(depth), d2 <= d1, in the same sense as S - as Hussain et\n'
                          'al. (2016): written on the profile with their 95 %% ranges,\n'
                          'and the MAP model without the creep dotted.')
    mod.add_argument('--creep-offset', dest='creep_offset', type=float, nargs='*',
                     metavar='KM',
                     help='the creep rate straight from the data, as Hussain et al.\n'
                          '(2016): a straight line through the pixels within KM of\n'
                          'the fault (default 5) on each side, both within WIDTH of\n'
                          'the profile line (--creep-offset KM WIDTH, default 2.5 km),\n'
                          'and their offset at the fault, its 95 %% range by bootstrap\n'
                          '(the pixels taken as independent). LOS of an ascending and\n'
                          'a descending track: split into fault-parallel and vertical.\n'
                          'Drawn on the profile; with --creep or without.')
    mod.add_argument('--dip', dest='dip_range', nargs='+', metavar='DEG',
                     help='reverse / normal: the dip of the fault, in degrees from\n'
                          'the horizontal - one value: fixed (--dip 30); two: the\n'
                          'range of its uniform prior (--dip 10 45; default: 5 85).')
    mod.add_argument('--dip-toward', dest='dip_toward', default='end',
                     choices=['end', 'start'],
                     help="reverse / normal: the side the fault dips towards - end\n"
                          "(default: towards A', x > 0; the hanging wall on that\n"
                          'side) or start (towards A).')
    mod.add_argument('--free-x0', dest='free_x0', type=float, nargs='?', const=5.0,
                     metavar='KM',
                     help='fit the fault position too, within +/-KM of it\n'
                          '(KM default: 5).')
    mod.add_argument('--rot', dest='rot', type=float, default=0.0, metavar='RATE',
                     help='a fixed linear term, RATE x, on one side of the fault, in\n'
                          'mm/yr per km - the rotation of a block about a nearby\n'
                          'pole (Hussain et al. 2016: 0.0215). Default: none.')
    mod.add_argument('--rot-side', dest='rot_side', default='end', choices=['end', 'start'],
                     help="the side of --rot: end (default, towards A', x > 0) or\n"
                          'start (towards A, x < 0).')
    mod.add_argument('--s-max', dest='s_max', type=float, metavar='NUM',
                     help='|S| <= NUM (default: 60 mm/yr; reverse / normal: the slip\n'
                          'rate along the dip, 0 <= S <= NUM).')
    mod.add_argument('--d1-max', dest='d1_max', type=float, default=60.0, metavar='KM',
                     help='0.1 <= d1 <= KM (default: %(default)s).')
    mod.add_argument('--c-max', dest='c_max', type=float, default=30.0, metavar='NUM',
                     help='0 <= C <= NUM (default: %(default)s mm/yr).')
    mod.add_argument('--d2-max', dest='d2_max', type=float, default=40.0, metavar='KM',
                     help='0 < d2 <= KM, and d2 <= d1 (default: %(default)s).')
    mod.add_argument('--d1-prior', dest='d1_prior', type=float, nargs=2,
                     metavar=('MEAN', 'SIGMA'),
                     help='a Gaussian prior on d1, in km, on top of its bounds\n'
                          '(Hussain et al. 2016, Izmit: 17 5).')
    mod.add_argument('--noise', dest='noise', nargs='+', default=['auto'], metavar='NOISE',
                     help='auto (default): an exponential covariance, sill and\n'
                          'length fitted to the residuals of a least-squares fit,\n'
                          'for each file; white: the scatter of each bin, no\n'
                          'correlation; SILL LENGTH: a given exponential covariance,\n'
                          '(mm/yr)^2 and km (e.g. from a far-field area).')

    mc = parser.add_argument_group('Sampler', 'the ensemble MCMC')
    mc.add_argument('--walkers', dest='walkers', type=int, default=64, metavar='NUM',
                    help='number of walkers (default: %(default)s).')
    mc.add_argument('--steps', dest='steps', type=int, default=5000, metavar='NUM',
                    help='steps of every walker (default: %(default)s).')
    mc.add_argument('--burn', dest='burn', type=int, metavar='NUM',
                    help='steps dropped as burn-in (default: a third).')
    mc.add_argument('--seed', dest='seed', type=int, default=0, metavar='NUM',
                    help='random seed, for repeatable runs (default: %(default)s).')

    gns = parser.add_argument_group('GNSS', 'GNSS velocities, fitted with the InSAR')
    gns.add_argument('--gnss', dest='gnss_file', metavar='FILE',
                     help='GNSS velocities, any format swath reads (.csv .xlsx .shp\n'
                          '.geojson .gpkg .kml, GMT psvelo ...). Those within\n'
                          '--gnss-width of the line are fitted with an offset of\n'
                          'their own, and drawn in blue with their 1-sigma.')
    gns.add_argument('--vel', dest='gnss_comp', nargs='+',
                     choices=['par', 'perp', 've', 'vn', 'vu', 'los'], metavar='COMP',
                     help='the GNSS component that matches the data, one of par perp\n'
                          've vn vu los. strike: default par, fault-parallel, from\n'
                          've and vn (los when --comp auto fits LOS data as they\n'
                          'are); --comp as-is needs it. reverse / normal: perp -\n'
                          'across the fault, from ve and vn (default with ve / vn\n'
                          'data, or LOS only: then a panel of their own) - or vu\n'
                          '(default with vu data only). One only: it is fitted with\n'
                          'the data.')
    gns.add_argument('--gnss-unit', dest='gnss_unit', default='auto',
                     choices=['auto', 'mm/yr', 'cm/yr', 'm/yr'],
                     help='unit of the GNSS file (default: auto).')
    gns.add_argument('--gnss-width', dest='gnss_width', type=float, metavar='KM',
                     help='take the GNSS within KM of the line (default: --width).')
    gns.add_argument('--gnss-nofit', dest='gnss_fit', action='store_false',
                     help='show the GNSS, but leave them out of the fit.')
    gns.add_argument('--gnss-label', dest='gnss_label', action='store_true',
                     help='write the site names on the map.')

    look = parser.add_argument_group('Plot', 'the fit window')
    look.add_argument('--ylim', dest='model_ylim', type=float, nargs=2,
                      metavar=('YMIN', 'YMAX'), help='velocity range of the plot.')
    look.add_argument('--strain', dest='show_strain', action='store_true',
                      help='a panel with the shear strain rate dv/dx of the model.')
    look.add_argument('--cmap-prof', dest='cmap_prof', nargs='?', const='viridis',
                      metavar='CMAP',
                      help='draw the pixels as a pixel DENSITY in this colormap, as\n'
                           'swath does (alone: viridis) - with the 95 / 68 %% ranges in\n'
                           'red, bigger medians and the GNSS red, and no points.\n'
                           'Default: the pixels as red points, the ranges grey, the\n'
                           'GNSS blue.')
    look.add_argument('--bands', dest='bands', type=float, nargs=2, default=[95.0, 68.0],
                      metavar=('OUTER', 'INNER'),
                      help='the two ranges of the model shaded, in %% (default: 95 68 -\n'
                           'the 68 %% range is +/- 1 sigma, as Hussain et al.).')
    look.add_argument('--legend', dest='model_legend', nargs='?', const='auto',
                      default='auto', choices=['auto', 'in', 'out', 'off'], metavar='WHERE',
                      help='the legend: auto (default) inside the plot where it hides\n'
                           'the least, room made for it - outside, on the right, when\n'
                           'the plot is too full for it; in: inside anyway; out:\n'
                           'outside, on the right; off: none. On the profile, "l"\n'
                           'moves it: in, out, off, in turn.')
    look.add_argument('--max-points', dest='max_points', type=int, default=200000,
                      metavar='NUM',
                      help='draw at most NUM pixels as points, a random subset\n'
                           '(default: %(default)s); all of them are fitted.')

    topo = parser.add_argument_group('Topography', 'the map background, an elevation panel')
    topo.add_argument('--dem', dest='dem_file', nargs='?', const=AUTO_DEM, metavar='DEM',
                      help='DEM for the shaded-relief map, as in view.py; without a\n'
                           'file it is downloaded (Copernicus GLO-30), as in swath.')
    topo.add_argument('--topoxsect', dest='topoxsect', nargs='?',
                      const=AUTO_DEM, metavar='DEM',
                      help='an elevation panel under the fit (as in swath).')
    add_dem_source(topo)
    topo.add_argument('--dem-zoom', dest='dem_zoom', type=int, metavar='Z',
                      help='zoom of the Terrain Tiles (--dem-source terrain), 0-15\n'
                           '(default: auto).')
    topo.add_argument('--lookup', dest='lookup_file', metavar='FILE',
                      help='radar-coded data: MintPy geometryRadar.h5 (default:\n'
                           'found in inputs/, ./ or ../inputs/ as MintPy does).')

    out = parser.add_argument_group('Save/Output', 'figures and data - see "output" below')
    out.add_argument('--outfile', type=str, nargs='*',
                     help='file name of the saved MAP figure.')
    add_save_options(out, 'every section given up front\n(--start-lalo ...)')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not open any window (implies --save).')
    out.add_argument('--svg', dest='svg', action='store_true',
                     help='write the figures as vector .svg instead of .png.')
    out.add_argument('--figext', dest='fig_ext', default='.png', metavar='EXT',
                     choices=['.png', '.pdf', '.svg', '.eps', '.jpg'],
                     help='figure format: .png .pdf .svg .eps .jpg (default: %(default)s).')
    out.add_argument('--dpi', dest='fig_dpi', type=int, default=300, metavar='NUM',
                     help='resolution of raster figures (default: %(default)s).')
    out.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='size of each fit window, in inches (default: 7.6 4.8,\n'
                          'taller with --strain / --topoxsect).')
    out.add_argument('--figsize-map', dest='fig_size_map', metavar=('WID', 'LEN'),
                     type=float, nargs=2, help='size of the map window, in inches.')
    return parser


def _parse_noise(values):
    """--noise: 'auto', 'white' or (sill, length)."""
    if len(values) == 1 and values[0].lower() in ('auto', 'white'):
        return values[0].lower()
    if len(values) == 2:
        try:
            sill, length = (float(v) for v in values)
        except ValueError:
            sill = length = None
        if sill is not None and sill >= 0 and length > 0:
            return (sill, length)
    raise ValueError(f'--noise: auto, white, or SILL LENGTH (two positive numbers); '
                     f'got {" ".join(values)}')


def check_model_gnss(inps):
    """The GNSS component that matches --comp (inps: the options, or the viewer)."""
    fault = getattr(inps, 'fault_type', 'strike')
    if fault != 'strike':
        comps = list(getattr(inps, 'comps', None) or [inps.comp])
        if inps.gnss_comp is None:    # across the fault: GNSS horizontals are the good ones
            inps.gnss_comp = ('vu' if 'vu' in comps and not {'ve', 'vn'} & set(comps)
                              else 'perp')
        elif inps.gnss_comp not in DIP_SLIP_GNSS:
            raise ValueError(f'viewinsar model {fault}: --vel {inps.gnss_comp} - the GNSS of a '
                             'reverse or normal fault are perp (across the fault, from ve and '
                             'vn) or vu (vertical)')
        return inps
    if inps.gnss_comp == 'perp':
        raise ValueError('--vel perp: the GNSS across the fault are for a reverse or a normal '
                         'fault (viewinsar model reverse|normal ...); a strike-slip one moves '
                         'along it: --vel par')
    if inps.gnss_comp is None:
        if inps.gnss_file and inps.comp == 'as-is':
            raise ValueError('--comp as-is: say which GNSS component matches the data with '
                             '--vel (par, ve, vn, vu or los) - or give --comp')
        inps.gnss_comp = 'par'
    if (inps.gnss_file and inps.gnss_comp == 'par' and inps.comp in ('as-is', 'par')
            and inps.strike is None):
        raise ValueError(f'--comp {inps.comp} with fault-parallel GNSS: give --strike, the '
                         'azimuth the data are positive along')
    if inps.comp in ('par', 've', 'vn', 'los') and inps.gnss_comp != 'par':
        print(f'NOTE: --comp {inps.comp} fits fault-parallel velocity; the GNSS are taken as '
              f'--vel {inps.gnss_comp} all the same')
    return inps


def check_model_comps(inps, given):
    """--comp of `viewinsar model strike|reverse|normal`: one for all the
    files or one per file; inps.comps, one per file, and inps.comp - the one,
    or 'mixed'. `given`: the options on the command line."""
    fault, n_file = inps.fault_type, len(inps.file)
    given_comp = inps.comp if isinstance(inps.comp, str) else ','.join(inps.comp)
    comps = [c for c in re.split(r'[,\s]+', given_comp.strip().lower()) if c]
    for c in comps:
        if c not in MODEL_COMPS:
            raise ValueError(f'--comp: {c} is not a component - one of '
                             + ', '.join(MODEL_COMPS))
    if len(comps) not in (1, n_file):
        raise ValueError(f'--comp: {len(comps)} components for {n_file} input file(s) - one '
                         'for all, or one per file')
    for word in inps.file:            # the fault type, put after the options
        if word.lower() in FAULT_TYPES and not os.path.exists(word):
            raise ValueError(f'{word}: the fault type comes first - viewinsar model '
                             f'{FAULT_TYPES[word.lower()]} FILE ... [options]')
    if fault == 'strike':
        if len(set(comps)) > 1:
            raise ValueError(f'--comp {" ".join(comps)}: a strike-slip fault takes ONE --comp '
                             'for all the files - they all become fault-parallel velocity. One '
                             'per file is for reverse and normal faults')
        if comps[0] == 'vu':
            raise ValueError('--comp vu: the vertical velocity tells of a reverse or a normal '
                             'fault, not of a strike-slip one - viewinsar model reverse --comp vu '
                             '... (or normal); --comp as-is fits it as it is anyway')
        inps.comp, inps.comps = comps[0], [comps[0]] * n_file
        wrong = given_options(given, DIP_SLIP_ONLY)
        if wrong:
            raise ValueError(f'{wrong[0]}: {DIP_SLIP_ONLY[wrong[0]]} of a reverse or a normal '
                             f'fault - viewinsar model reverse|normal ... {wrong[0]}')
        inps.dip_range = (5.0, 85.0)
        return inps
    for c in comps:
        if c not in DIP_SLIP_COMPS:
            raise ValueError(f'viewinsar model {fault}: --comp {c} - a {fault} fault is fitted '
                             'to vu (the vertical velocity), ve / vn (the east / north velocity: '
                             'the motion across the fault) or los (the LOS as it is), one per '
                             'file - or auto (by the file name, else los)')
    wrong = given_options(given, STRIKE_ONLY)
    if wrong:
        raise ValueError(f'{wrong[0]}: {STRIKE_ONLY[wrong[0]]} is for a strike-slip fault '
                         f'(viewinsar model strike ...), not a {fault} one')
    inps.comps = comps * n_file if len(comps) == 1 else comps
    inps.comp = inps.comps[0] if len(set(inps.comps)) == 1 else 'mixed'
    # the dip: fixed, or the range of its prior
    given_dip = (inps.dip_range if isinstance(inps.dip_range, (str, type(None)))
                 else ','.join(inps.dip_range))
    values = ([5.0, 85.0] if given_dip is None
              else [v for v in re.split(r'[,\s]+', given_dip.strip()) if v])
    try:
        values = [float(v) for v in values]
    except ValueError:
        raise ValueError(f'--dip: one angle (fixed) or two (a range), in degrees; got '
                         f'{given_dip}') from None
    if len(values) not in (1, 2) or not all(0 < v <= 90 for v in values):
        raise ValueError(f'--dip: one angle (fixed) or two (a range), each above 0 and at most '
                         f'90 degrees; got {" ".join(f"{v:g}" for v in values)}')
    inps.dip_range = (min(values), max(values))
    return inps


def model_cmd_line_parse(iargs=None):
    """Parse the `model` arguments; unknown ones are forwarded to view.py."""
    iargs = _argv(iargs, 'model')
    fault_type, rest, named = split_fault_type(iargs)
    parser = create_model_parser()
    inps, view_extra = parser.parse_known_args(args=model_argv(rest))
    inps.argv = list(iargs)
    inps.fault_type = fault_type
    if named and inps.print_msg:
        print(f'model: a {fault_type} fault' + {
            'strike': ' - the arctan of Savage & Burford (1973)',
            'reverse': ' - a dislocation dipping into the crust, slipping up-dip below d1',
            'normal': ' - a dislocation dipping into the crust, slipping down-dip below '
                      'd1'}[fault_type])
    inps.view_extra = view_extra
    check_model_comps(inps, rest)
    if isinstance(inps.gnss_comp, list):
        comps = list(dict.fromkeys(inps.gnss_comp))
        if len(comps) > 1:
            raise ValueError(f'--vel {" ".join(comps)}: model fits ONE GNSS component with '
                             'the data - the one that matches them (par, ve, vn, vu or '
                             'los; perp or vu for a reverse or normal fault). Several at once '
                             'are for swath.')
        inps.gnss_comp = comps[0]

    check_save_options(inps)
    if inps.svg:
        inps.fig_ext = '.svg'
    positive = dict(width=('--width', inps.width), bin_km=('--bin', inps.bin_km),
                    d1_max=('--d1-max', inps.d1_max), c_max=('--c-max', inps.c_max),
                    d2_max=('--d2-max', inps.d2_max))
    for opt, value in positive.values():
        if value <= 0:
            raise ValueError(f'{opt} must be positive, got {value}')
    if inps.d1_max <= 0.1:
        raise ValueError(f'--d1-max must exceed 0.1 km, got {inps.d1_max}')
    for opt, value in (('--s-max', inps.s_max), ('--free-x0', inps.free_x0),
                       ('--gnss-width', inps.gnss_width)):
        if value is not None and value <= 0:
            raise ValueError(f'{opt} must be positive, got {value}')
    if inps.d1_prior is not None and inps.d1_prior[1] <= 0:
        raise ValueError(f'--d1-prior: SIGMA must be positive, got {inps.d1_prior[1]}')
    if inps.nlook < 1:
        raise ValueError(f'--nlook must be >= 1, got {inps.nlook}')
    if inps.nlook > 1 and inps.decimate > 1:
        raise ValueError('--nlook and --decimate are opposites; use one or the other.')
    if inps.walkers < 2:
        raise ValueError(f'--walkers must be at least 2, got {inps.walkers}')
    if inps.steps < 20:
        raise ValueError(f'--steps must be at least 20, got {inps.steps}')
    if inps.burn is not None and not 0 <= inps.burn < inps.steps:
        raise ValueError(f'--burn must be 0 or more and below --steps ({inps.steps}), '
                         f'got {inps.burn}')
    if inps.max_points is not None and inps.max_points < 1:
        raise ValueError(f'--max-points must be positive, got {inps.max_points}')
    if inps.cmap_prof:
        check_colormap(inps.cmap_prof, '--cmap-prof')
    if inps.creep_offset is not None:
        values = list(inps.creep_offset)
        if len(values) > 2:
            raise ValueError(f'--creep-offset takes KM and WIDTH at most, got {values}')
        km, width = (values + [5.0, 2.5][len(values):])[:2]
        if km <= 0 or width <= 0:
            raise ValueError(f'--creep-offset: KM and WIDTH must be positive, got {km} {width}')
        inps.creep_offset = (float(km), float(width))
    inps.noise = _parse_noise(inps.noise)
    outer, inner = inps.bands
    if not 0 < inner < outer < 100:
        raise ValueError(f'--bands OUTER INNER: 0 < INNER < OUTER < 100, got {outer:g} {inner:g}')

    # the GNSS component must match the data (--comp auto: once it is settled)
    if 'auto' not in inps.comps:
        check_model_gnss(inps)
    elif inps.gnss_comp is not None and (
            inps.gnss_comp not in DIP_SLIP_GNSS if inps.fault_type != 'strike'
            else inps.gnss_comp == 'perp'):
        check_model_gnss(inps)        # wrong whatever the data are: said now

    # what the shared machinery of the swath viewer expects to find
    for key, value in dict(bin_x=inps.bin_km, bin_y=None, ylim=None, x_offset=0.0, cmin=1.0,
                           cmax=None, pcmap='viridis', timeseries=None, yx=None,
                           gnss_ref='none', gnss_radius=1.0, fault_label='Fault').items():
        setattr(inps, key, value)
    check_map_options(inps)
    if is_auto(inps.topoxsect):
        inps.topoxsect = inps.dem_file if inps.dem_file else AUTO_DEM
    elif inps.topoxsect and not os.path.isfile(inps.topoxsect):
        raise ValueError(f'--topoxsect: no such file: {inps.topoxsect}')

    read_file_info(inps)
    if inps.geom_file and len(inps.geom_file) > len(inps.file):
        raise ValueError(f'--geom: {len(inps.geom_file)} file(s) for {len(inps.file)} '
                         'input file(s)')
    if inps.los_inc and len(inps.los_inc) > len(inps.file):
        raise ValueError(f'--los-inc / --los-az: {len(inps.los_inc)} values for '
                         f'{len(inps.file)} input file(s)')

    collect_profiles(inps)
    if not inps.disp_fig and not inps.profiles:
        raise ValueError('--nodisplay needs the profiles up front (--start-lalo / --end-lalo, '
                         '--start-yx / --end-yx, --line-file or --profile)')

    if view_extra and inps.print_msg:
        print(f'options forwarded to view.py: {" ".join(view_extra)}')
    return inps


def model_main(iargs=None):
    """Entry point of `viewinsar model`."""
    from .modelview import ModelViewer

    inps = model_cmd_line_parse(iargs)
    obj = ModelViewer(inps, get_view_cmd(inps))
    obj.open()
    obj.plot()
    if obj.fig_map is not None and getattr(obj, 'cid', None) is not None:
        obj.fig_map.canvas.mpl_disconnect(obj.cid)
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

With --radius KM a double click takes a CIRCLE of that radius: every pixel within
it, relative to the reference area and the reference date as LiCSBAS does a
single pixel, and each date shows their MEDIAN, with the 25th - 75th
percentiles as an error bar. The fits follow the medians.

Several files - adjacent frames of a track, overlapping tracks, or versions of
one time series - share ONE MAP, each drawn where it lies with its outline in a
colour of its own. A double click shows the series of every file that covers
the point, each relative to its own reference area: right-drag one where they
overlap (or give --ref-geo) to put them on the same ground.
"""

TSERIES_EXAMPLE = """\
examples:
  # LiCSBAS: in TS_GEOCml*/ (cum_filt.h5, and cum.h5 in red, as LiCSBAS does)
  viewinsar tseries
  viewinsar tseries cum_filt.h5 --radius 2                 # circles of 2 km
  viewinsar tseries cum_filt.h5 --radius 1.5 --ref-date 20170105 \\
      --ref 140:145/125:130

  # MintPy (velocity.h5, maskTempCoh.h5 ... found next to it)
  viewinsar tseries timeseries_ERA5_ramp_demErr.h5 --radius 2 --unit cm
  viewinsar tseries geo/geo_timeseries_ERA5_demErr.h5 --radius 3 --lalo -2.51 121.68

  # no window: save the time series of one circle (.csv + figures)
  viewinsar tseries cum_filt.h5 --radius 2 --p-geo 121.68/-2.51 --nodisplay

  # several files, one map: adjacent LiCSBAS frames of a track
  viewinsar tseries 061D_*/TS_GEOCml10GACOS/cum.h5 --radius 15
  viewinsar tseries frameA/cum.h5 frameB/cum.h5 --labels A B --radius 15 \\
      --ref-geo 121.60/121.70/-2.60/-2.50

mouse / keyboard:
  left double click : the time series of that pixel (of the circle, with --radius)
  right drag        : a new reference area
  slider            : the cumulative displacement at that date
  "e"               : export, see below

output (written by "e", --save or --nodisplay, into the --save folder - made if
it does not exist - or --outdir, else the current folder), named after the
point, e.g. tseries_Y156X144_r2km:
  <name>.csv        date, year, then per file: median_mm, q25_mm, q75_mm,
                    n_px (with --radius), or disp_mm (one pixel)
  <name>.png        the time series (.svg with --svg); <name>_map.png the map
  <name>_info.txt   the file(s), the point, the circle, the reference area
                    and date, the velocity

several files (one map):
  the files       geocoded in one CRS (LiCSBAS GEOC, MintPy geo_*), or radar-coded
                  on one grid (the MintPy time series of one stack)
  reference area  each file's own, until a right drag (or --ref-geo, --ref-lalo)
                  gives one: every file with data in it takes it, the others keep
                  theirs. Drag one into the overlap of two frames: their series
                  become comparable (dashed on the map, in the file's colour)
  reference date  --ref-date, else the latest of the files' first dates (MintPy:
                  REF_DATE), so that all start together; a file without that date
                  takes its nearest one
  dates           each file keeps its own; the slider shows each file at its
                  date nearest the slider's, and none before its first or after
                  its last
  pixel options   --ref X1:X2/Y1:Y2, --ref-yx, --yx, --point and --mask: files on
                  one grid only. --i2, --results and --vel-file: one file only
  output          named after the point: tseries_lat-2.5100_lon121.6800_r15km
                  (projected: _N.._E..); the .csv has every date of the files and
                  a set of columns per file: median_mm_1 ... n_px_2 ("nan" where a
                  file has no such date)

LiCSBAS users: --radius here is the RADIUS of the circle; the reference area that
LiCSBAS_plot_ts.py takes with -r is --ref X1:X2/Y1:Y2. LiCSBAS -p X/Y (here --point),
--p_geo / --ref_geo work as in LiCSBAS.
"""


def _radius_km(text):
    if re.search(r'[:/]', str(text)):
        raise argparse.ArgumentTypeError(
            f'{text!r}: --radius is the circle radius in km here; give the reference area '
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
    parser = ViewinsarParser(
        prog=prog, description=TSERIES_DESCRIPTION, epilog=TSERIES_EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter, allow_abbrev=False)
    parser.add_argument('file', nargs='*',
                        help='the time series: LiCSBAS cum_filt.h5 / cum.h5, or MintPy\n'
                             '(default: ./cum_filt.h5, ./cum.h5 or ./timeseries.h5).\n'
                             'Several files - frames of a track, tracks, versions of\n'
                             'one time series - share one map: see "several files".')
    parser.add_argument('--labels', dest='labels', nargs='+', metavar='NAME',
                        help='names of the input files in the legends (default: the\n'
                             'file names, or the folders that tell them apart).')
    parser.add_argument('--radius', dest='radius', type=_radius_km, metavar='KM',
                        help='a circle of that radius, in km: each date shows the\n'
                             'median of its pixels, with the 25th-75th percentiles\n'
                             'as an error bar (default: one pixel).')
    parser.add_argument('--i2', dest='file2', metavar='FILE',
                        help='a second time series, in red (default: cum.h5 when the\n'
                             'file is cum_filt.h5, as LiCSBAS).')
    parser.add_argument('--no-i2', dest='no_i2', action='store_true',
                        help='no second time series.')
    parser.add_argument('--ref-date', dest='ref_date', metavar='YYYYMMDD',
                        help='reference date, where every series is 0 (default: the\n'
                             'first date; MintPy: its REF_DATE; several files: the\n'
                             'latest of theirs).')
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
    pnt.add_argument('--point', dest='point_xy', metavar='X/Y',
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
    data.add_argument('--vel-file', dest='vel_file', metavar='FILE',
                      help='velocity file (default: LiCSBAS vel; MintPy velocity.h5\n'
                           'beside the file, else a linear fit of the time series).')
    data.add_argument('--data-scale', dest='data_scale', type=float, nargs='+',
                      metavar='FACTOR',
                      help='multiply the time series (and its velocity) by FACTOR\n'
                           'before anything: one per file (and one for --i2), or\n'
                           'one for all: 0.1 for a product 10x too large. A rule in\n'
                           '~/.config/viewinsar/viewinsar.ini does it for every file\n'
                           'of a mission: [data-scale] ALOS-2 = 0.1.')
    data.add_argument('--lookup', dest='lookup_file', metavar='FILE',
                      help='radar-coded MintPy data: geometryRadar.h5, for lat/lon\n'
                           '(default: found where MintPy keeps it).')

    disp = parser.add_argument_group('Display')
    disp.add_argument('--colormap', dest='colormap', default='roma_r', metavar='CMAP',
                      help='colormap of the velocity / displacement (default: roma_r,\n'
                           'as LiCSBAS; any MintPy or matplotlib name).')
    disp.add_argument('--unit', dest='unit', default='mm', choices=['mm', 'cm', 'm'],
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
    out.add_argument('--outfile', type=str, nargs='*',
                     help='file name of the saved time-series figure.')
    add_save_options(out, 'the first time series')
    out.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                     help='save and do not open any window (implies --save).')
    out.add_argument('--svg', dest='svg', action='store_true',
                     help='write the figures as vector .svg instead of .png.')
    out.add_argument('--figext', dest='fig_ext', default='.png', metavar='EXT',
                     choices=['.png', '.pdf', '.svg', '.eps', '.jpg'],
                     help='figure format: .png .pdf .svg .eps .jpg (default: %(default)s).')
    out.add_argument('--dpi', dest='fig_dpi', type=int, default=300, metavar='NUM',
                     help='resolution of raster figures (default: %(default)s).')
    return parser


def _time_series_files(names):
    """The input files, in the order given: names, or wildcards (quoted, or on
    Windows) - each name once."""
    import glob

    files = []
    for name in names:
        hits = sorted(glob.glob(name)) if glob.has_magic(name) else (
            [name] if os.path.isfile(name) else [])
        if not hits:
            raise ValueError(f'no such file: {name}')
        files += [f for f in hits
                  if os.path.abspath(f) not in {os.path.abspath(g) for g in files}]
    return files


def tseries_cmd_line_parse(iargs=None):
    """Parse the `tseries` arguments.

    inps.files is the list of input files; inps.file the first of them.
    """
    from .tsdata import find_time_series

    iargs = _argv(iargs, 'tseries')
    parser = create_tseries_parser()
    inps = parser.parse_args(args=_bare_save(iargs))
    inps.argv = list(iargs)

    if not inps.file:
        found = find_time_series('.')
        if not found:
            raise ValueError('no time series given, and no cum_filt.h5, cum.h5 or '
                             'timeseries.h5 here')
        inps.files = [found]
    else:
        inps.files = _time_series_files(inps.file)
    inps.file = inps.files[0]
    several = len(inps.files) > 1
    if inps.labels and len(inps.labels) != len(inps.files):
        raise ValueError(f'--labels: {len(inps.labels)} name(s) for {len(inps.files)} file(s)')
    if several:
        # what belongs to one file only; several files bring their own
        for opt, value in (('--i2', inps.file2), ('--results', inps.results_dir),
                           ('--vel-file', inps.vel_file)):
            if value:
                raise ValueError(f'{opt} is for one time series: with several files, each '
                                 'brings its own (LiCSBAS results/, MintPy velocity.h5 ...)'
                                 + (' - give every series as an input' if opt == '--i2'
                                    else ''))
        inps.file2 = None
    elif inps.no_i2:
        inps.file2 = None
    elif inps.file2:
        if not os.path.isfile(inps.file2):
            raise ValueError(f'--i2: no such file: {inps.file2}')
    elif os.path.basename(inps.file) == 'cum_filt.h5':
        cand = os.path.join(os.path.dirname(os.path.abspath(inps.file)), 'cum.h5')
        inps.file2 = cand if os.path.isfile(cand) else None
    check_data_scale(inps, files=inps.files + ([inps.file2] if inps.file2 else []))
    for opt, path in (('--mask', inps.mask_file), ('--vel-file', inps.vel_file),
                      ('--lookup', inps.lookup_file)):
        if path and not os.path.isfile(path):
            raise ValueError(f'{opt}: no such file: {path}')
    if inps.results_dir and not os.path.isdir(inps.results_dir):
        raise ValueError(f'--results: no such folder: {inps.results_dir}')
    if inps.ref_date and not re.fullmatch(r'\d{8}', inps.ref_date):
        raise ValueError(f'--ref-date: a date is YYYYMMDD, got {inps.ref_date}')

    points = [o for o, v in (('--yx', inps.yx), ('--lalo', inps.lalo), ('--point', inps.point_xy),
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

    check_save_options(inps)
    if inps.svg:
        inps.fig_ext = '.svg'
    return inps


def tseries_main(iargs=None):
    """Entry point of `viewinsar tseries`: one file, or several on one map."""
    inps = tseries_cmd_line_parse(iargs)
    if len(inps.files) > 1:
        from .tsmulti import MultiTimeSeriesViewer
        obj = MultiTimeSeriesViewer(inps)
    else:
        from .tsviewer import TimeSeriesViewer
        obj = TimeSeriesViewer(inps)
    obj.open()
    obj.plot()
    if not obj.disp_fig:
        obj.close()
    return obj


def create_list_parser(prog=f'{PROG} list'):
    """The argument parser of the `list` sub-command."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description='List the datasets available in a file - the names that\n'
                    '`viewinsar swath --dset` accepts. The one swath uses by\n'
                    'default is marked with "<- default".',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('file', help='input file')
    parser.add_argument('--num', type=int, default=0, metavar='NUM',
                        help='print the first NUM names only (default: all).')
    return parser


def list_main(iargs=None):
    """Entry point of `viewinsar list` - print the datasets of a file."""
    from .compat import get_slice_list

    inps = create_list_parser().parse_args(args=_argv(iargs, 'list'))

    names = get_slice_list(inps.file)
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
    'stat': (stat_main, 'fault creep along a line, from short profiles across it'),
    'model': (model_main, 'fit slip rate and locking depth to a profile (MCMC)'),
    'tseries': (tseries_main, 'time series of a pixel or a circle (LiCSBAS / MintPy)'),
    'list': (list_main, 'list the datasets available in a file'),
}


RULE = '-' * 88

TOP_HELP = """\
viewinsar {version} - swath profiles, fault creep along a line, slip-rate models and
time series of InSAR data (MintPy, LiCSBAS); MintPy itself is optional.

usage: {prog} <command> [options]
       {prog} FILE [options]          shortcut for `{prog} swath FILE ...`
       {prog} <command> -h            every option of one command
       {prog} help <command>          the same
       {prog} --help-all              this, then every option of every command
       {prog} --version

commands:
{cmds}

data: MintPy HDF5 (velocity.h5, timeseries.h5, geo_*.h5, .he5), LiCSBAS cum*.h5 (vel,
cum-YYYYMMDD) and results/ (vel, vel.filt.mskd ...), GeoTIFF (with MintPy or
rasterio) and ROI_PAC .dem/.rsc - geocoded (lat/lon or UTM)
or radar-coded (lat/lon from MintPy's geometryRadar.h5). Several files at once:
frames, tracks, components. The line of sight of every pixel comes from MintPy's
geometry file, LiCSBAS's GEOCml*/E.geo, N.geo, U.geo or LiCSAR's *.geo.E.tif, found
next to the data; a LiCSBAS velocity is read in mm/yr (--data-unit for any other
file without a unit).

{rule}
swath - swath (buffered) profiles
{rule}
  Two clicks on the map draw a section: every pixel within +/-WIDTH km of the line
  is projected onto it and drawn as a pixel-density profile with the median and the
  25th-75th percentiles. Each section opens in its own window: A-A', B-B', C-C' ...
  the line      two left clicks; or up front --start-lalo / --end-lalo (--start-yx /
                --end-yx), --profile LAT0 LON0 LAT1 LON1 (repeat it: B-B', C-C' ...),
                --line-file (GMT segments or a vector file: one section per line)
  the swath     --width half width in km, --bin-x / --bin-y density cells, --ylim,
                --cmin / --cmax / --pcmap, --x-offset, --keep-zero, --decimate,
                --nlook N (oversample a coarse raster so the swath fills in)
  elevation     --topoxsect [DEM] a profile under the data, --dem [DEM] shaded relief
                on the map - either without a file downloads Copernicus GLO-30
                (--dem-source); radar-coded data take the height of geometryRadar.h5
  faults        --fault FILE ... (.shp .geojson .gpkg .kml .gmt, any CRS): on the map,
                and a dashed line where a section crosses one (--fault-label,
                --fault-color, --fault-lw)
  GNSS          --gnss FILE (.csv .xlsx .shp .geojson .gpkg .kml, GMT psvelo): on the
                map and on every profile; --vel los ve vn vu (one or several),
                --gnss-ref, --gnss-radius, --gnss-label, --geom / --los-inc
  several files f1.h5 f2.h5 ... --labels: one density, a median line for each file
  time series   --timeseries ts.h5: right-click a pixel for its history (MintPy)
  keys          left click x2: draw | "=" new section | "n" redo the current one |
                "e" export the map + every section | right click: time series
  output        swath_<A>_..._points.csv (every pixel: distance along / off the
                line, value, row, col, lat, lon), _stats.csv (the bins), _info.txt,
                the figure; swath_map_ABC.png
  examples      {prog} swath velocity.h5 --width 15
                {prog} swath velocity.h5 --fault faults.shp --dem --topoxsect --nlook 4
                {prog} swath f1.h5 f2.h5 --gnss gnss.csv --vel los --svg --save
                {prog} swath velocity.h5 --timeseries timeseries.h5
                {prog} swath velocity.h5 --line-file profiles.gmt --nodisplay --save

{rule}
stat - fault creep along a line, from short profiles across it
{rule}
  Click points along the fault (Enter to finish), or give the line with --line-lalo,
  --line-yx or --line-file. Every --spacing km along it a profile crosses it: --length
  km either side (default 5), the pixels within --width km of it (2.5). A straight line
  through each side, and the creep rate is their offset at the line - block 2 (RIGHT,
  walking from the first point to the last) minus block 1 (LEFT) - with a bootstrap
  that resamples each side on its own (--boot), as Hussain et al. (2016).
  the figure    the LOS creep rate along the line, a series per file, 1-sigma error
                bars (descending solid, ascending dashed); ascending + descending
                (--geom, --los-inc / --los-az): split into fault-parallel and vertical
  a profile     click a point of the creep rate: its pixels, the least-squares
                (BLUE) and bootstrap lines, the PDFs (arrows: other points / files)
  options       --names / --colors of the blocks, --ylim, --keep-zero, --min-pixels,
                --dem [DEM] relief, --fault FILE on the map, --lookup (radar-coded)
  keys          left click: add a point | Enter: finish | "u" / Backspace: undo |
                "n": new line | "e": export
  output        stat_..._creep.csv (each point and file: the creep rate, its 1 sigma,
                68 / 95 % ranges, the BLUE, the slopes), _split.csv, _values.csv
                (every pixel), _profiles.png, _info.txt, the figure and the map
  examples      {prog} stat velocity.h5 --fault faults.shp
                {prog} stat velocity.h5 --spacing 2 --length 4 \\
                    --line-file trace.geojson --save
                {prog} stat desc.h5 asc.h5 --labels T036D T131A --line-file trace.geojson

{rule}
model - slip rate and locking depth from profiles across a fault
{rule}
  The pixels of a swath against their distance from the fault, binned (--bin km),
  and an elastic dislocation fitted by MCMC (the Goodman & Weare sampler)
  fault type    {prog} model [strike|reverse|normal] FILE ... - strike (the default):
                v(x) = S/pi atan(x/d1) + a (+ creep, rotation), as Hussain et al.
                (2016, 2018); reverse (thrust) / normal: a fault dipping at dip,
                locked down to d1, slipping at S along its dip below (Savage 1983)
                - --dip DEG [DEG], --dip-toward end|start; shortening / extension
                S cos(dip), uplift / subsidence S sin(dip)
  profiles      as swath: two clicks, or up front --start-lalo / --end-lalo,
                --profile (repeat it: A-A', B-B' ... each fitted on its own) and
                --line-file (one per line)
  the fault     where the profile crosses --fault FILE, with its strike there; or
                --x0 KM / --strike DEG
  the data      strike: --comp auto|los|par|ve|vn|as-is - auto (default) turns a
                LOS velocity into fault-parallel velocity, pixel by pixel, when its
                line of sight is found (MintPy geometry, LiCSBAS E/N/U, LiCSAR
                *.geo.E.tif), named by --geom or given by --los-inc / --los-az;
                reverse / normal: --comp vu|ve|vn|los, one per file (--comp ve vu),
                a panel each; several files fitted together, each with an offset
                of its own
  the model     S and d1; --creep (C, d2), --free-x0 [KM], --rot / --rot-side,
                --d1-prior MEAN SIGMA, --s-max --d1-max --c-max --d2-max,
                --noise auto|white|SILL LENGTH
  the sampler   --walkers, --steps, --burn, --seed
  keys          as swath: left click x2 draws a profile and fits it | "=" another |
                "n" redo the current one | "e" export
  GNSS          --gnss FILE --vel par|ve|vn|vu|los (reverse / normal: perp|vu),
                fitted with an offset of their own (--gnss-nofit: shown only),
                --gnss-width
  figures       the fit: the pixels as points, binned medians, GNSS, the MAP model
                and its 95 / 68 % ranges (--bands), the parameters written on it as
                Hussain et al. ("slip rate = 20 mm/yr (11 - 29, 95% CI), right-
                lateral", locking depth in km, strain rate in µstrain/yr), the
                legend inside or beside it (--legend; "l" moves it); --cmap-prof
                [CMAP]: the pixels as a density instead, as swath (--nlook fills
                it in), the ranges red, bigger medians, GNSS red; --strain: dv/dx
                below, --topoxsect: elevation; the posteriors as Hussain et al.
                (2016, central NAF, Fig. 9): every pair a jet 2D histogram, the
                marginals in blue, the MAP in red - one of each per section
  output        model_<A>_...: .png, _posterior.png, _params.csv / .txt, _samples.csv,
                _bins.csv, _points.csv, _gnss.csv; model_map_ABC.png, and
                model_sections_ABC.csv (two sections or more: one row each)
  examples      {prog} model velocity.h5 --comp los --fault faults.shp --width 25
                {prog} model asc.h5 desc.h5 --comp los --fault faults.shp --gnss gnss.csv
                {prog} model velocity.h5 --comp los --fault faults.shp --creep \\
                    --profile 40.2 33.1 41.3 33.4 --profile 40.3 33.9 41.4 34.2 \\
                    --nodisplay --save fits
                {prog} model reverse geo_vu.h5 --comp vu --fault faults.shp --dip 10 45
                {prog} model normal ve.h5 vu.h5 --comp ve vu --fault faults.shp

{rule}
tseries - time series, after LiCSBAS_plot_ts.py
{rule}
  The velocity map (or a noise index, or the displacement at a date on the slider);
  double-click a pixel for its time series with LiCSBAS's fits (Linear, Annual+L,
  Quad, Annual+Q); right-drag a new reference area.
  --radius KM   a circle instead of one pixel: the median of its pixels at every
                date, the 25th-75th percentiles as error bars
  the data      LiCSBAS cum_filt.h5 / cum.h5 (results/: mask, noise indices; --i2 a
                second series), MintPy timeseries*.h5 / .he5 (its mask and velocity
                found beside it); --mask, --nomask, --vel-file, --results, --lookup
  several files frames of a track, tracks, versions of one series: ONE MAP, a colour
                each; a double click gives the series of every file there, each on
                its own reference area until a right drag puts them on one;
                --labels NAME ...
  reference     --ref X1:X2/Y1:Y2, --ref-geo, --ref-yx, --ref-lalo; --ref-date YYYYMMDD
  the point     --yx, --lalo, --point X/Y, --p-geo LON/LAT (as LiCSBAS)
  display       --colormap CMAP, --unit mm|cm|m, --vlim, --dlim, --auto-crange,
                --ylim, --noflip
  keys          double click: the time series | right drag: reference | "e": export
  output        tseries_..._r2km.csv (date, median, q25, q75, n_px per file), the
                figure, _map.png, _info.txt
  examples      {prog} tseries cum_filt.h5 --radius 2
                {prog} tseries timeseries_ERA5_demErr.h5 --radius 2 --unit cm \\
                    --lalo -2.51 121.68
                {prog} tseries 061D_*/TS_GEOCml10GACOS/cum.h5 --radius 15

{rule}
common options (swath, stat, model, tseries)
{rule}
  --save [DIR]  write what is given up front at start-up, and make "e" export into
                DIR (made if it does not exist)
  --outdir DIR, --outfile FILE, --nodisplay, --svg / --figext EXT, --dpi, --figsize,
  --figsize-map, --noverbose; --dset for swath, stat and model (`{prog} list FILE`
  shows the datasets of a file)
  map           swath, stat, model: other options go to the map - all of view.py with
                MintPy (--colormap, --vlim, --unit, --mask, --sub-lat, --sub-lon ...),
                the common ones without it
  the values    --data-unit UNIT the unit they are really in; --data-scale FACTOR
                multiplies them first (0.1: a product 10x too large) - for every
                file of a mission, a rule in ~/.config/viewinsar/viewinsar.ini:
                    [data-scale]
                    ALOS-2 = 0.1
  spelling      every option is --long (-h for help); the one-dash spellings of
                before 0.5.0 (-fault, -gnss, -vel, -d, -w, -sp, -r ...) still work

install:
  pip install viewinsar              numpy scipy matplotlib h5py pyproj pyshp pillow
  pip install "viewinsar[mintpy]"    + MintPy: the full view.py map, swath --timeseries
  pip install "viewinsar[raster]"    + rasterio: GeoTIFF / GDAL rasters without MintPy
  pip install "viewinsar[vector]"    + pyogrio: every GDAL vector format
  pip install "viewinsar[all]"       all of them

Run `{prog} <command> -h` for every option of a command."""


def top_level_help():
    """What `viewinsar` / `viewinsar -h` prints: every command, what it does,
    its options in groups, its keys, its output and examples."""
    width = max(len(k) for k in SUBCOMMANDS)
    cmds = '\n'.join(f'  {name:<{width}}  {helptxt}'
                     for name, (_, helptxt) in SUBCOMMANDS.items())
    return TOP_HELP.format(version=__version__, prog=PROG, cmds=cmds, rule=RULE)


def command_help(name):
    """The full -h of one command."""
    makers = {'swath': create_swath_parser, 'stat': create_stat_parser,
              'model': create_model_parser, 'tseries': create_tseries_parser}
    if name == 'list':
        return create_list_parser().format_help()
    return makers[name]().format_help()


def full_help():
    """`viewinsar --help-all`: the overview, then every option of every command."""
    parts = [top_level_help()]
    for name in SUBCOMMANDS:
        parts.append('\n'.join(['', '#' * 88, f'#  {PROG} {name} -h', '#' * 88, '',
                                command_help(name)]))
    return '\n'.join(parts)


def main(argv=None):
    """Top-level `viewinsar` entry point."""
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv and argv[0] in ('--help-all', '-H') or argv[:2] == ['help', 'all']:
        print(full_help())
        return 0
    if argv[:1] == ['help'] and len(argv) > 1 and argv[1] in SUBCOMMANDS:
        print(command_help(argv[1]))
        return 0
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

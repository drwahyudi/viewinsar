# Changelog

All notable changes to viewinsar are recorded here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project uses [semantic versioning](https://semver.org/).

## [Unreleased]

## [0.5.5] - 2026-11-17

### Added
- **`model --cmap-prof [CMAP]`: the pixels as a point density** instead of
  points, drawn as `swath` draws its density (log colour scale up to the
  98th percentile of the occupied cells, the densest white, a colour bar;
  cells half a `--bin` wide, ~55 over the data). Under it the 95 % range is
  red, the medians per bin are twice the size (white, or a colour per file
  with several - each named in the legend), and no point is drawn. Alone,
  `--cmap-prof` takes viridis; any matplotlib, MintPy or cmcrameri name works,
  and an unknown one is refused. The default view - red points, a grey
  range - is as before.

### Fixed
- **`model --nlook` no longer goes into the fit.** It oversamples what is drawn
  (the points, or the density) and the fit takes the real pixels. Before, the
  interpolated samples were fitted as data: the medians hardly moved, but each
  bin seemed to hold N x N times more independent pixels, so the 95 % ranges
  came out too narrow (half as wide with `--nlook 3` in a test). Fits made
  with `--nlook` before 0.5.5 should be redone; without it nothing changes.
  `_points.csv` lists the fitted pixels, `_params.txt` says when the drawing
  was oversampled.

## [0.5.4] - 2026-11-17

### Added
- **`viewinsar tseries` takes several files** - adjacent frames of a track,
  overlapping tracks, or versions of one time series
  (`viewinsar tseries 061D_*/TS_GEOCml10GACOS/cum.h5 --radius 15`). They share
  one map, each drawn where it lies, with its outline, reference area and
  series in a colour of its own; a double click shows the series of every
  file that covers the point (in an overlap, all of them). Each file keeps
  its own dates and, until a right drag (or `--ref-geo`, `--ref-lalo`) gives
  one where it has data, its own reference area: drag one into the overlap
  of two frames and their series become comparable. The reference date is
  common: `--ref-date`, else the latest of the files' first dates (a file
  without it takes its nearest). The slider shows each file at its date
  nearest the slider's. The files must be geocoded in one CRS, or
  radar-coded on one grid (MintPy's `timeseries.h5` and
  `timeseries_ERA5.h5`, as `tsview.py`); the options in pixels (`--ref`,
  `--ref-yx`, `--yx`, `--point`) and `--mask` need one grid. The export is
  named after the lat/lon of the point, and the `.csv` has every date of the
  files, with columns `_1`, `_2` ... per file.
- **`--labels`** for `tseries`: the names of the files in the legends.
- Wildcards in quotes work for `tseries` too (`"061D_*/TS_GEOC*/cum.h5"`,
  for a shell that does not expand them).

### Changed
- Files with the same name in folders of the same name
  (`061D_09234_131313/TS_GEOCml10GACOS/cum.h5`, `061D_09012_161516/...`) are
  named in legends by the folders that tell them apart (`061D_09234_131313`),
  in `swath`, `stat` and `model` as in `tseries`; before, both were
  `TS_GEOCml10GACOS/cum`.

### Fixed
- `tseries`: the map window's tick labels no longer hide behind the date
  slider (a tall map) or the help text and buttons on the left (a wide one).
  Its parts are placed in inches now, and the window is an inch wider, and
  taller when the left column needs it (LiCSBAS's dozen noise indices).
- Docs: the GLO-30 tiles have overviews at 1/2, 1/4 and 1/8 resolution (not
  only 1/2 and 1/4). The reader always used whichever exist; checked bit for
  bit against GDAL on a real tile.

## [0.5.3] - 2026-11-15

### Changed
- **`--dem` / `--topoxsect` without a file now download the Copernicus DEM
  GLO-30** (1", ~30 m, heights above the EGM2008 geoid; AWS Open Data,
  `copernicus-dem-30m`) instead of the Terrain Tiles, in `swath`, `stat`
  and `model`. Only the parts of the 1 x 1 deg Cloud Optimized GeoTIFFs over
  the map are read, by HTTP byte range - no GDAL or rasterio needed - from
  the overview nearest the grid being made, and the tiles are joined
  seamlessly. The elevation profile gets 1" pixels up to 8000 a side (a map
  of ~2.2 deg), the map background up to 3000 (view.py shades every pixel
  it is given); a bigger map gets wider pixels. The sea, which has no tiles,
  is 0 m; a tile that cannot be downloaded leaves its area empty. Radar-coded
  data whose geometry file has no `height` get GLO-30 resampled onto their
  grid. The Copernicus licence notice is printed and written into the DEM
  file and `_info.txt`. Everything read is cached in `~/.cache/viewinsar`.
- **GNSS points on the profiles are outlined in white** - marker edges and
  error bars - in `swath` and `model`, so they stand out on the data behind
  them. With several components (`--vel ve vn vu los`), `vu` is now blue
  instead of white.

### Added
- **`--dem-source cop30|terrain`** (`swath`, `stat`, `model`): where a
  download comes from. `terrain` takes the Terrain Tiles as before
  (`--dem-zoom` applies to them only, and says so otherwise). Without
  `--dem-source`, a GLO-30 that cannot be had (no access to its server, say)
  falls back to the Terrain Tiles, with a warning.

## [0.5.2] - 2026-11-15

### Added
- **`swath --vel` takes several components at once**: `--vel ve vn vu los`
  draws all of them on the same profile, each in a colour of its own (ve
  magenta, vn orange, vu white, los black; triangles for 3D stations, circles
  for 2D), each shifted onto the InSAR on its own (`--gnss-ref insar`). The
  console and `_info.txt` give each shift; `_gnss.csv` has a row per station
  and component, with a new `component` column. A component that cannot be
  had - `vu` with no 3D station, `los` with no look angles - is left out
  with a warning, the others stay. One component is drawn as before.
  `model --vel` takes one, the component it fits, and says so if given more.

## [0.5.1] - 2026-11-15

### Removed
- **The smoothed line through the GNSS on the swath profiles.** The stations
  are drawn with their 1-sigma error bars, and nothing else. `swath
  --gnss-nofit`, which only turned that line off, is gone; given anyway, it
  is left out with a note. (`model` never drew such a line; its
  `--gnss-nofit` - the GNSS shown, but not fitted - stays.)

## [0.5.0] - 2026-11-14

0.4.3 was built but never released; its changes are part of this one.

### Changed
- **Every option is written `--long`**, one spelling each: `--fault`, `--gnss`,
  `--vel`, `--dem`, `--width`, `--nlook`, `--topoxsect`, `--profile`, `--svg`,
  `--save`, `--outfile`; `stat`: `--spacing`, `--radius`, `--rate`; `tseries`:
  `--radius`, `--ref-date`, `--point`, `--colormap`, `--unit`; `list`: `--num`.
  The spellings of before (`-fault`, `-gnss`, `-vel`, `--gnss-comp`, `-d`,
  `-w`, `--buffer`, `-nlook`, `-topoxsect`, `--topo`, `-profile`, `-svg`,
  `-save`, `-o`, `-sp`, `-r`, `-rate`, `-m`, `-p`, `-c`, `-u`, `-n`) still
  work, with a note that names the new ones. view.py's own options (`-v`,
  `-u`, `-c`, `-m` ... for the map) are view.py's and pass as they are; the
  help and the README write them long (`--vlim`, `--unit`, `--colormap`,
  `--mask`).
- `tseries --vel FILE` is `--vel-file FILE`, so that `--vel` means one thing
  (the GNSS component of `swath` and `model`); `--vel` still works there.
- `model --comp auto` fitting LOS data as they are takes the GNSS in the LOS
  (`--vel los`) without being told, or the component a file is named as
  (`up` -> `vu`, `east` -> `ve` ...), and writes the rates as "slip rate in
  the LOS", "strain rate in the LOS" - what the LOS sees of them.

### Added
- **`--geom FOLDER`**: a folder is searched, with its subfolders, for the line
  of sight - LiCSBAS's `GEOCml*/` (E.geo, N.geo, U.geo), E / N / U GeoTIFFs
  under any of their usual names (`{frame}.geo.E.tif`, `E.geo.tif`, `E.tif`,
  `*_e.tif`, `*_east.tif`), incidence + azimuth GeoTIFF pairs
  (`incidenceAngle.tif` + `azimuthAngle.tif`, `inc.tif` + `az.tif` ...) and
  MintPy's geometry file - and the first that covers the data is taken. The
  interferogram folders of LiCSBAS's `GEOC/` are skipped. `--geom` the LiCSBAS
  work folder, or a folder of geometry GeoTIFFs, and `model --comp auto` turns
  the LOS into fault-parallel velocity.
- **`--data-scale FACTOR`** (swath, stat, model, tseries): multiplies the data
  before anything - the map, the profiles, the fit, the statistics, the time
  series - one factor per file, or one for all; `0.1` for a product 10x too
  large.
- **Your own defaults, `~/.config/viewinsar/viewinsar.ini`**: `[data-scale]
  ALOS-2 = 0.1` scales every file of a mission (told by its `PLATFORM`, or by
  its name and folders), `*/path/pattern*.h5 = 0.1` every file whose path
  matches. Every run says when a rule applies; `--data-scale` comes first.
- A note when a file's `WAVELENGTH` is not its mission's (ALOS-2 0.238 m,
  Sentinel-1 0.0555 m): a displacement computed with it is wrong by the same
  ratio.
- **A velocity that says no unit gets one from its numbers.** A GeoTIFF
  (or any file) without a unit, used as a velocity - in `model`, with
  `--gnss`, or named `vel*` - is read in m/yr when 99% of its values are
  below 1 and in mm/yr otherwise, and is then shown in mm/yr, like the
  GNSS. The console says what it chose and why; `--data-unit` overrides it.
  Before, such a file stayed unitless ('1'): the GNSS were compared with it
  unconverted and the model was fitted "not a velocity", with a slip-rate
  bound of +/-4.
- **Nominal line-of-sight angles from the file's name**, as a last resort
  when no geometry is found: a name (or folder, or `ORBIT_DIRECTION`) that
  says the orbit - `asc`, `desc`, `scansardesc`, LiCSAR's `105D_...`,
  `T125A` - gives the heading of a sun-synchronous orbit at the scene's
  latitude and a mid-swath incidence (Sentinel-1 IW 37.5 deg, ALOS-2 38
  deg), right-looking. They put the GNSS in the LOS (swath and model), with
  a note that they are approximate; they never turn the InSAR into
  fault-parallel velocity: `--comp auto` then fits the LOS as it is, and
  `--comp los` asks for `--geom` or `--los-inc` / `--los-az`.
- MintPy's geometry file is also looked for in `mintpy/inputs/`,
  `mintpy/geo/` and one or two folders up.

### Fixed
- `-u` / `--unit` with a file whose unit viewinsar supplies (`--data-unit`,
  LiCSBAS, the guess): the map stayed in the file's numbers while the
  profiles were converted.

## [0.4.2] - 2026-11-14

### Added
- **LiCSBAS and LiCSAR without typing the look angles.** The line of sight of
  every pixel comes from LiCSBAS's `GEOCml*/E.geo`, `N.geo`, `U.geo` (with
  `EQA.dem_par`) - found by itself for a file anywhere under `TS_GEOCml*/`,
  e.g. `TS_GEOCml10GACOS/results/vel.filt.mskd.geo.tif` → `GEOCml10GACOS/` -
  or from LiCSAR's `{frame}.geo.E.tif`, `.N.tif`, `.U.tif`, as well as from
  MintPy's geometry file. `--geom` takes any of them (or the `GEOCml*`
  folder), one per input file, for `model --comp los` and `swath -vel los`.
- **LiCSBAS's own files open directly** in swath, stat and model, as MintPy
  files and GeoTIFFs do: `TS_GEOCml*/cum_filt.h5` / `cum.h5` (the velocity
  `vel` by default, `--dset cum-YYYYMMDD` an epoch, `viewinsar list` the
  rest) and the float32 files of `results/` (`vel`, `vel.filt.mskd`, `vstd`
  ...) on the grid of `GEOCml*/EQA.dem_par` - in mm/yr (mm for `cum`), with
  LiCSBAS's `results/mask`, and the line of sight of every pixel from
  `GEOCml*/E.geo`, `N.geo`, `U.geo`.
- **The unit of a LiCSBAS velocity**: `vel*` and `vstd*` in a LiCSBAS folder
  are read in mm/yr - the GeoTIFFs say nothing - so the slip rate is in
  mm/yr and the strain rate in µstrain/yr.
- **`--data-unit UNIT`** (swath, stat, model): the unit the input files are
  really in - for a file that does not say it, or says a wrong one (values
  in mm/yr labelled m/year or cm/year). The map, its colour range, the
  profiles, the fit and the comparison with the GNSS all follow; the console
  says which unit each file is read in.
- **`model --comp auto`, the new default**: a LOS velocity whose line of sight
  is known becomes fault-parallel velocity, pixel by pixel; without one (or
  for a file named as a component - hz, up, east, ve ...) the data are
  fitted as they are. The run says which.
- `--los-inc` / `--los-az` take one value per input file (ascending and
  descending in one run).

### Changed
- **The parameters on the fit are written as Hussain et al. (2016)**:
  `slip rate = 20.3 mm/yr (18.4 – 34.0, 95% CI), right-lateral`, `locking depth
  = 12.7 km (…)`, creep rate, creep depth, and `strain rate` in µstrain/yr -
  the slip rate as a rate with its sense. A velocity axis without a unit no
  longer reads "[1]".
- With MintPy but without GDAL, GeoTIFFs are read with rasterio (and drawn on
  viewinsar's own map) instead of failing.

### Fixed
- The map now follows the UNIT a file says, as the profiles always did:
  view.py takes any `velocity` dataset as m/year whatever its UNIT, so a
  velocity labelled mm/year or cm/year was drawn 1000x or 100x too large on
  the map.

## [0.4.1] - 2026-11-13

### Added
- **`viewinsar model`: slip rate and locking depth from a profile across a
  fault**, after Hussain et al. (2016, 2018) and Ou et al. (2022). An elastic
  screw dislocation (Savage & Burford 1973) - optionally with shallow creep
  (`--creep`: rate C down to d2 <= d1), a free fault position (`--free-x0`)
  and a fixed rotation term (`--rot`) - fitted to the swath with the
  affine-invariant ensemble MCMC sampler of Goodman & Weare (2010), under
  uniform priors (`--s-max`, `--d1-max`, `--c-max`, `--d2-max`, and a Gaussian
  `--d1-prior`), with an exponential noise covariance fitted to the residuals
  (or `--noise white` / `SILL LENGTH`). The pixels are binned (`--bin`, 1 km)
  and drawn as points, not as a density: the binned medians, the GNSS, the MAP
  model (dashed) with its 95 % range (grey), the parameters with their 95 %
  ranges - placed where they hide the fewest points - and a second window with
  the posterior of every parameter and pair. The fault position and strike
  come from a `--fault` trace (or `--x0`, `--strike`); distances are taken
  perpendicular to the fault. `--comp los / ve / vn` turn LOS (pixel by pixel,
  through MintPy's geometry file), east or north velocities into
  fault-parallel velocity, so that S > 0 is right-lateral and the sense of
  slip is written out. Several files (tracks, frames) are fitted together,
  each with its own offset; GNSS within `--gnss-width` with theirs.
  `--strain` adds the shear strain rate of the model. `e` / `--save` writes
  the fit, the posterior, `_params.csv` / `.txt` (MAP, median, 95 % range,
  the strain rate at the fault |S| / (pi d1)), `_samples.csv`, `_bins.csv`,
  `_points.csv` and `_gnss.csv`. In Python: `fit_profile()`.
- **Several files in one swath**: `viewinsar swath FILE1 FILE2 ...` - the
  frames of a track, or tracks side by side, in any CRS - samples the line in
  every file, pools them into one density and draws a median line per file
  (`--labels`); the other files are drawn under the first on the map.
- **GNSS on the swath** (`-gnss FILE`): .csv / .txt / .xlsx tables with any
  column spelling, GMT psvelo, .shp, .geojson, .gpkg, .kml / .kmz, and
  anything GDAL reads (`read_gnss()`). `-vel los` (default) projects them into
  the LOS of the data, through MintPy's geometry file or the data's
  attributes (`--geom`, `--los-inc`, `--los-az`); `-vel ve / vn / vu` takes a
  component. On the map every station, on each profile those inside its
  swath: black triangles (3D) and circles (2D) with 1-sigma error bars and a
  smoothed black line, shifted onto the InSAR (`--gnss-ref`). `_gnss.csv`.
- **`viewinsar stat -rate`**: the median of block 1 minus that of block 2 at
  every point, with error bars (`--rate-err iqr` or `sem`) - the slip or
  creep rate along the fault (Hussain et al. 2016). `stat` takes several
  files (e.g. the east and up velocities): a row of violins each, or a rate
  series each; `_rate.csv`, and `_2_stats.csv` ... for the other files.
- **Several sections up front** (`swath`, `model`): `--profile LAT0 LON0 LAT1
  LON1`, repeated, and `--profile-yx`; `--line-file` takes several lines - GMT
  segments split by `>`, or the lines of a vector file (.geojson .kml .shp
  .gpkg ..., re-projected from its CRS). They become `A-A'`, `B-B'`, `C-C'`
  ..., each in its own window and with its own files; `model` fits each on
  its own (its own fit and posterior figure) and writes
  `model_sections_ABC.csv`, one row per section: its ends, where it crosses
  the fault, the strike, the sense of slip and every parameter with its 95 %
  range - the slip rate along strike.
- **The posterior figure as Hussain et al. (2016, Fig. 9)**: every pair of
  parameters (S, d1, the offset, and C, d2, x0 when fitted) a 2D histogram in
  jet with a frequency colour bar, the marginals as blue histograms along the
  bottom with the MAP as a red line and the 95 % range dashed and written
  above, the MAP a red dot on every pair, and the section (`A - A'`) white on
  black.
- **`--save [DIR]` for every command**: what is given up front is written at
  start-up and `e` exports into DIR, made if it does not exist (`--outdir`
  and `-o` make their folders too).
- **`viewinsar -h` documents every command** - what swath, stat, model and
  tseries do, their options in groups, their keys, their output and examples;
  `viewinsar --help-all` adds every option of every command, `viewinsar help
  <command>` those of one.

### Changed
- **MintPy is optional, and the dependencies are what viewinsar imports**:
  numpy, scipy, matplotlib, h5py, pyproj, pyshp and pillow (shapely is gone:
  a small WKB reader of its own). Without MintPy, viewinsar reads MintPy HDF5,
  LiCSBAS HDF5, ROI_PAC .dem/.rsc and - with rasterio - GeoTIFF itself, and
  draws the map with the common view.py options; with it, everything is as
  before. Extras: `[mintpy]`, `[raster]` (rasterio), `[vector]` (pyogrio),
  `[all]`. `VIEWINSAR_NO_MINTPY=1` ignores an installed MintPy.
- `--line-file` with more than two points per line takes its first and last
  point (it took the first two).

### Fixed
- Single-dash options are no longer abbreviated on Python < 3.13 (`-v` was
  read as `-vel`, `-n` as `-nlook`), so view.py's `-v VMIN VMAX` reaches the
  map.

## [0.4.0] - 2026-11-13

### Added
- **`viewinsar tseries`: LiCSBAS_plot_ts.py, with circles.** An interactive
  time-series viewer for LiCSBAS cum_filt.h5 / cum.h5 and MintPy time series
  (timeseries*.h5, geo_timeseries*.h5, HDF-EOS5 .he5), laid out as LiCSBAS's:
  the velocity, the noise indices or the cumulative displacement of a date
  on the map, double-click for the time series, right-drag a new reference
  area, the four fit models, the second file (cum.h5) in red, the gaps of the
  interferogram network. With `-r KM` a double click takes every pixel within
  that radius, and each date shows their median with the 25th-75th
  percentiles as error bars. `e` / `--save` writes the series as .csv, both
  figures and `_info.txt`.
- **`viewinsar stat`: block statistics on both sides of a fault.** Draw a line
  along a fault (left clicks, Enter), or give it with `--line-lalo`,
  `--line-yx` or `--line-file`. Every `-sp` km along it a point is placed, and
  at each point two circles of radius `-r` km, one on each side, their edges
  `--gap` km (default 1) clear of the line. Block 1 is the left of the line,
  walking from its first point to its last, block 2 the right - those of the
  ground, also on (mirrored) radar grids. The values in each circle are drawn
  as violins, block 1 and block 2 side by side, each in the colour of its
  circles on the map, with strikes at the median, the mean and the mode (the
  peak of the kernel density), and the medians of each block joined by a
  dashed line. No circle ever reaches
  across, or within the gap of, the line: on the inside of a bend a circle is
  pushed further out, or left out where there is no room. `e` / `--save`
  writes the violins, the map, `_stats.csv` (one row per point: count,
  median, mean, mode, std, quartiles and range of both blocks, and their
  median difference), `_values.csv` (every pixel used, with its
  distance from the line) and `_info.txt`. Also from Python: `stat_line()`.
- **Radar-coded data** (MintPy products that are not geocoded) now get a DEM,
  an elevation profile and fault traces, through MintPy's geometry file
  `geometryRadar.h5` (the latitude, longitude and height of every pixel, on the
  same grid). It is found the way MintPy finds it - `inputs/`, `./` or
  `../inputs/`, next to the data file or in the current folder - or named with
  the new `--lookup FILE`.
  - `-d` / `--topoxsect` without a file use its `height`: the DEM of your own
    processing, already on the radar grid, so nothing is downloaded. Only a
    geometry file without height falls back to the Terrain Tiles, resampled
    onto the radar grid through its latitude / longitude.
  - `--topoxsect dem.tif` with a geocoded DEM reaches the radar section
    through the lookup.
  - `--fault` traces are mapped from lon/lat onto the radar pixels; `_info.txt`
    gives each crossing's row/col and lat/lon.
  Before, all three stopped with "needs geocoded data".

### Fixed
- **`--topoxsect geometryGeo.h5` (any MintPy geometry file) profiled the azimuth
  angle**, not the height: with no dataset named, MintPy returns the first one,
  `azimuthAngle`. A geometry file now always gives its `height`, as in
  MintPy's own DEM reader.

## [0.3.0] - 2026-11-11

### Added
- **A DEM without a DEM file.** `-d` given without a file downloads one for the
  map extent from the Terrain Tiles on AWS Open Data - the "Mapzen Global
  Terrain" layer of QGIS, terrarium PNG tiles in Web Mercator - and draws it as
  view.py's shaded-relief background. `--topoxsect` without a file profiles
  that same DEM, downloading it if `-d` was not given; `-d FILE --topoxsect`
  profiles your own file, named once. The tiles are mosaicked, re-gridded onto
  the data's CRS (lat/lon or UTM) and saved as a MintPy geometry `.h5` in
  `~/.cache/viewinsar` (`$VIEWINSAR_CACHE`), so a second run downloads nothing.
  The zoom is the finest that keeps the download within ~300 tiles (at most 12,
  ~38 m); `--dem-zoom` overrides it. The data attribution is printed and
  written to `_info.txt`.

### Fixed
- **`--topoxsect` gave "no topography for this profile" with `--nlook`.** The
  DEM was oversampled by the data's `--nlook`, so a 30-90 m DEM under a
  15 km-wide section asked for hundreds of millions of samples and was refused.
  The DEM is now sampled at its own resolution: thinned to at most ~4 million
  samples per section when it is finer than the median line needs, and
  oversampled only when its pixels are coarser than half a `--bin-x` bin.
- **A projected DEM could be taken for lat/lon.** MintPy's GDAL reader gives a
  UTM GeoTIFF whose WKT carries no EPSG code of its own (common from ArcGIS)
  the attribute `EPSG = 4326`, so the profile landed off the DEM. The CRS of a
  GeoTIFF is now read from the file itself, an EPSG:4326 label on a grid in
  metres is no longer trusted, and UTM zones are understood however they are
  written (`51S`, `UTM zone 51S`, `WGS_1984_UTM_Zone_51S`). `--fault` uses the
  same CRS detection.
- The elevation line and `elevation_median_m` are now the exact median per
  bin; the density-weighted median snapped to the histogram's 50-100 m height
  cells.
- A DEM covering only part of a section gives that part, instead of nothing.
- A section without topography says why in its panel, and in `_info.txt`,
  instead of only in the terminal. `_info.txt` also records the CRS of the data
  and of the DEM, and how the DEM was sampled.
- A file without a unit - a GeoTIFF, which MintPy gives the unit "1" - is no
  longer labelled "LOS velocity [1]".

## [0.2.0] - 2026-11-11

### Added
- **`--nlook N`** - oversample the corridor N times finer than the native pixel
  (N**2 more samples, NaN-aware bilinear), so a multilooked or otherwise coarse
  raster gives a filled density profile instead of a handful of scattered
  cells. The inverse of `--decimate`; the two cannot be combined.
- **`--topoxsect DEM`** - an elevation panel under the LOS profile, sharing its
  distance axis. The same swath is taken from the DEM on its own grid, in any
  format (MintPy `.h5`, GeoTIFF, `.dem` + `.rsc`, any GDAL raster) and any CRS -
  the profile is re-projected with pyproj when the two differ. Only the
  density-weighted median is drawn, as a black line, and it is written to the
  `_stats.csv` as `elevation_median_m`.
- **Several cross-sections.** Press `=` (or `+`) to start another one: it opens
  in a new window and the line you already drew stays on the map, labelled
  `A-A'`, then `B-B'`, `C-C'` ... `e` exports the map with every line plus each
  section's figure, points, statistics and metadata. `n` still clears just the
  current section.
- **Separate windows.** The map, each section and the point time series are now
  their own windows instead of panels crammed into one figure. `--figsize-map`
  sizes the map; `--figsize` sizes each section.
- **`--svg`** - write the figures as vector `.svg`. `--figext` and `--dpi` are
  now proper options of the swath command.
- **`--fault FILE ...`** - fault traces from any geocoded vector file, in any
  CRS: `.shp` (+ `.prj`), `.geojson`, `.gpkg`, `.kml`/`.kmz` and GMT text are
  read directly (no GDAL needed); with `pyogrio`, `fiona` or GDAL installed,
  every other OGR format too. The traces are re-projected to the data's CRS and
  drawn on the map. Wherever a section's centre line crosses one, its profile
  gets a dashed vertical line annotated "Fault line" (on the elevation panel
  too), and the crossing distance and coordinates go to `_info.txt`.
  `--fault-label`, `--fault-color` and `--fault-lw` style them. Unlike MintPy's
  own `--faultline`, it works on UTM maps and marks the profiles.
- **Single-dash spellings**: `-nlook`, `-topoxsect`, `-svg` and `-fault` work as
  well as their `--` forms. Option abbreviation is now off, so view.py's own
  short options (e.g. `-n`) can never be mistaken for one of ours.
- `viewinsar.is_projected()` and `viewinsar.bilinear_sample()` are public.
- **The help documents everything.** `viewinsar swath -h` now explains what the
  tool does, walks through the A-A' / B-B' workflow step by step, lists every
  output file and its columns, and has examples for each new option.
  `viewinsar -h` gives a one-screen overview; `viewinsar list` marks the
  dataset swath uses by default. A test parses every example in the help
  against the real parsers, so a mistyped flag in the docs fails the build.
- `--start-lalo` / `--end-lalo` on a UTM file accept either lat/lon or
  northing/easting (they are told apart by size); the help now says so.

### Fixed
- **Projected rasters without a `UTM_ZONE` tag had every distance inflated by
  ~111,320x.** `get_pixel_size()` decided degrees-vs-metres from `UTM_ZONE`
  alone, so a GeoTIFF carrying only `EPSG` (as GDAL writes) had its metre-sized
  pixels treated as degrees - a 58 km profile came out as 4,974,152 km. It now
  uses `Y_UNIT`/`X_UNIT`, `UTM_ZONE`, the EPSG code and finally the pixel size
  itself.
- **A `velocity.h5` with no `--dset` could export `intercept` instead of
  `velocity`** - the wrong band, anti-correlated with the velocity and larger by
  the stack length in years. The default dataset now matches `FILE_TYPE`, and a
  time-series file defaults to the LAST epoch (the cumulative displacement)
  rather than the all-zero reference date.
- The point density colour limit is now auto-scaled from the occupied cells, so
  a sparsely sampled profile is as readable as a dense one and `--nlook` does
  not saturate the panel. `--cmax` still overrides it.
- A profile on a projected grid is labelled `N/E`, not `lat/lon`.
- Figures are closed after a `--nodisplay` run instead of accumulating.

## [0.1.0] - 2026-08-25

First release.

### Added
- `viewinsar swath` - interactive swath (buffered) profile viewer.
  - Two left clicks draw the centre line on a `view.py` map; the line is
    buffered by `--width` km with flat caps, every valid pixel inside is
    linearly referenced onto it (QGIS LRS / ArcGIS *Locate Features Along
    Routes*), and the cloud is drawn as a log-scaled 2-D point-density profile
    with the density-weighted median and the 25th-75th percentile envelope.
  - Optional third panel with the point displacement time series, via
    `tsview.py` (`--timeseries`, right-click a pixel).
  - Key `e` exports the located points and the per-bin statistics to `.csv`,
    plus the profile metadata and the figure. Key `n` starts a new line.
  - Non-interactive use with `--start-lalo/--end-lalo`, `--start-yx/--end-yx` or
    `--line-file`, together with `--save --nodisplay`.
  - Multiple input files: the density comes from the first, the others are
    overlaid as median lines and included in the exported `.csv`.
  - Every unrecognised option is forwarded verbatim to MintPy's `view.py`.
- `viewinsar list` - print the datasets available in a file.
- `viewinsar.swath` - the numerical core (`swath_extract`, `swath_density`,
  `line_geometry`, `weighted_quantile`), importable with numpy alone.

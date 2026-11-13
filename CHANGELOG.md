# Changelog

All notable changes to viewinsar are recorded here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project uses [semantic versioning](https://semver.org/).

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

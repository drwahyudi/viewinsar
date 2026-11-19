# viewinsar

**Interactive swath (buffered) profile viewer for InSAR velocity and time-series
products, built on [MintPy](https://github.com/insarlab/MintPy).**

`plot_transection.py` samples a single line of pixels. `viewinsar` samples a
**corridor**: you draw the same line with two clicks, it gets buffered by
±*W* km, every valid pixel inside the buffer is linearly referenced onto the
line, and the resulting cloud is shown as a 2-D pixel-density profile with the
density-weighted median and the 25th–75th percentile envelope.

In GIS terms, one click-click does this:

| step | GIS equivalent |
|------|----------------|
| 1. buffer the line by ±*W* km, **flat caps** | QGIS `buffer()` with square cap style |
| 2. clip the raster with that polygon | *Clip raster by mask layer* |
| 3. project each pixel onto the line → along-profile *measure* + off-profile distance | QGIS **LRS** / ArcGIS **Locate Features Along Routes** |
| 4. bin (distance, value) into a log-scaled density grid | — |

The map, every cross-section and the point time series each get **their own
window**. Press `=` to add another section: the line you already drew stays on
the map as `A-A'`, the new one becomes `B-B'`, and so on. Give it fault traces
(`--fault faults.shp`, any vector format, any CRS) and they are drawn on the map
and marked with a dashed "Fault line" on every profile that crosses them.

The pieces come from three MintPy tools: the map is `view.py`, the
click-to-draw interaction is `plot_transection.py`, and the optional
displacement time series is `tsview.py`.

`viewinsar stat` measures **fault creep along a line**, as Hussain et al.
(2016) did on the North Anatolian Fault: draw the line on the fault, and every
*SP* km a short profile crosses it; a straight line through the pixels of each
side, and their offset at the fault is the creep rate, with a bootstrap for its
uncertainty - along the line, per track, and split into fault-parallel and
vertical when there are ascending and descending tracks (see
[Fault creep along a line](#fault-creep-along-a-line-viewinsar-stat)).

`viewinsar tseries` is LiCSBAS_plot_ts.py for LiCSBAS and MintPy time series,
with circles: `--radius 2` and a double click gives, at every date, the median of
all pixels within 2 km, with the 25th-75th percentiles as error bars. Several
files - the frames of a track, say - share one map, and a double click gives the
series of each file there (see [Time series](#time-series-viewinsar-tseries)).

Several files at once - the frames of a track, or several tracks - make one
density with a median line for each, and are drawn together on the map; GNSS
velocities (`--gnss`, any table or point format) appear on the map and, inside
a swath, on its profile with their error bars (see
[Several files and GNSS](#several-files-and-gnss)).

`viewinsar model` fits the **slip rate and locking depth** of an elastic
dislocation - and the rate and depth of shallow creep - to a profile across a
fault, by ensemble MCMC, as Hussain et al. (2016, 2018) did for the North
Anatolian Fault: the pixels of the swath as points, the maximum a posteriori
model with its 95 % range, and the posterior distributions; `--creep-offset`
measures the creep rate straight from the offset at the fault too. `viewinsar
model reverse` and `viewinsar model normal` fit a fault that dips - its slip
rate along the dip, locking depth and dip - to the vertical velocity, the
motion across the fault or a line of sight, a panel each (see
[Fault slip models](#fault-slip-models-viewinsar-model)).

---

## Install

```bash
pip install viewinsar
```

or, once the conda-forge feedstock is live:

```bash
conda install -c conda-forge viewinsar
```

From a checkout, for development:

```bash
git clone https://github.com/drwahyudi/viewinsar
cd viewinsar
pip install -e ".[test]"
pytest
```

**Requirements:** Python ≥ 3.9 and only what viewinsar imports itself:
`numpy` ≥ 1.20, `scipy`, `matplotlib` ≥ 3.5, `h5py`, `pyproj`, `pyshp` and
`pillow`. **MintPy is optional.** Without it viewinsar reads MintPy's HDF5
products, LiCSBAS `cum*.h5`, ROI_PAC `.dem/.rsc` and (with `rasterio`)
GeoTIFF itself, and draws the map with the common `view.py` options
(`--colormap --vlim --unit --mask --sub-lat --sub-lon --dem --ref-lalo --fontsize --title
--flip-lr --flip-ud --coord ...`). The extras add the rest:

| install | adds |
|---|---|
| `pip install viewinsar` | the core, above |
| `pip install "viewinsar[mintpy]"` | MintPy: the map from `view.py` itself (every option), `swath --timeseries`, every format MintPy reads |
| `pip install "viewinsar[raster]"` | `rasterio`: GeoTIFF and every GDAL raster (DEMs, LiCSBAS `.tif`) without MintPy |
| `pip install "viewinsar[vector]"` | `pyogrio`: `--fault`, `--gnss` and `--line-file` read every GDAL vector format, not only the common ones viewinsar reads by itself |
| `pip install "viewinsar[all]"` | all three |

If you already have a MintPy environment, install viewinsar into it and it
uses MintPy. `VIEWINSAR_NO_MINTPY=1` makes it ignore an installed MintPy.

---

## Quick start

```bash
# interactive - click 2 points on the map to draw the swath centre line
viewinsar swath geo_velocity_msk.h5 --width 15 --bin-x 0.5 --bin-y 1.5 --ylim -40 40

# non-interactive - give the line up front, write the figure + the .csv files
viewinsar swath geo_velocity.h5 \
    --start-lalo 30.125 129.988 --end-lalo 30.250 130.116 \
    --width 15 --save --nodisplay

# several sections up front - A-A', B-B', C-C' - saved into sections/
viewinsar swath geo_velocity.h5 --width 10 --nodisplay --save sections \
    --profile 30.10 129.95 30.30 130.10 --profile 30.20 130.05 30.40 130.20 \
    --profile 30.30 130.15 30.50 130.30
viewinsar swath geo_velocity.h5 --line-file profiles.gmt --save   # one per segment

# a coarse / multilooked raster: oversample the corridor so the profile fills in
viewinsar swath geo_velocity_ml10.h5 --nlook 10

# an elevation profile under the LOS panel, from a DEM in any format or CRS
viewinsar swath geo_velocity.h5 --topoxsect srtm1.dem

# no DEM file? the Copernicus GLO-30 DEM is downloaded for the map: shaded
# relief + elevation profile
viewinsar swath geo_velocity.h5 --dem --topoxsect

# fault traces on the map, and a dashed "Fault line" where a section crosses one
viewinsar swath geo_velocity.h5 --fault faults.shp
viewinsar swath geo_velocity.h5 --fault gem_active_faults.geojson --topoxsect dem.tif

# the point time series, in its own window (right-click a pixel on the map)
viewinsar swath geo_velocity.h5 --timeseries geo_timeseries_ERA5_demErr.h5

# vector output
viewinsar swath geo_velocity.h5 --topoxsect dem.tif --svg --save

# several frames or tracks in one swath, and GNSS on the profiles
viewinsar swath T014A/geo_velocity.h5 T123D/geo_velocity.h5 --labels T014A T123D
viewinsar swath geo_velocity.h5 --gnss gnss.csv --vel los

# fault creep along a line: a profile across it every 5 km, two lines and their offset
viewinsar stat geo_velocity.h5 --spacing 5 --fault faults.shp
# ... of an ascending and a descending track: also split into fault-parallel and vertical
viewinsar stat T036D/geo_velocity.h5 T131A/geo_velocity.h5 --line-file fault_trace.geojson

# slip rate and locking depth of a fault, by MCMC, from a profile across it
viewinsar model geo_velocity.h5 --comp los --fault faults.shp --width 25
# ... of a reverse fault, from the vertical velocity: slip rate, depth and dip
viewinsar model reverse geo_vu.h5 --comp vu --fault faults.shp --width 25

# time series, LiCSBAS-style: the median of a 2 km circle at every date
viewinsar tseries cum_filt.h5 --radius 2
viewinsar tseries geo_timeseries_ERA5_demErr.h5 --radius 2
# ... of every frame of a track, on one map
viewinsar tseries 061D_*/TS_GEOCml10GACOS/cum.h5 --radius 15

# list what is inside a file
viewinsar list timeseries.h5
```

`viewinsar FILE ...` is a shortcut for `viewinsar swath FILE ...`.
`viewinsar -h` prints what every command does - its options in groups, its
keys, its output and examples; `viewinsar --help-all` adds every option of
every command, and `viewinsar help model` (or `viewinsar model -h`) those of
one.

### Mouse and keyboard

| input | action |
|-------|--------|
| left click ×2 | set the start / end point of the current section's line |
| right click | plot the time series of that pixel (needs `--timeseries`) |
| <kbd>=</kbd> or <kbd>+</kbd> | start **another** section in a new window — `A-A'`, `B-B'`, `C-C'` … |
| <kbd>n</kbd> | redraw the **current** section from scratch |
| <kbd>e</kbd> | export the map (with every line) and each section |

Keys work from any window.

---

## Options

Anything `viewinsar swath` does not recognise is **forwarded verbatim to
`view.py`**, so the whole MintPy display vocabulary keeps working:

```bash
viewinsar swath geo_velocity.h5 --colormap jet --vlim -3 3 --unit cm --dem srtm1.dem \
    --mask geo_maskTempCoh.h5 --sub-lat 30.0 30.5 --sub-lon 129.9 130.3 \
    --coastline --lalo-label --fontsize 12
```

The swath-specific options:

| option | default | meaning |
|--------|---------|---------|
| `--width` | `15` | **half** width of the buffer, in km |
| `--bin-x` | `0.5` | along-profile cell size of the density grid, in km |
| `--bin-y` | auto | value cell size, in the display unit (~55 bins over the Y range) |
| `--ylim` | 0.5–99.5 % | value range of the profile plot |
| `--cmin` / `--cmax` | `1` / auto | limits of the log density colour scale; `--cmax` defaults to the 98th percentile of the occupied cells |
| `--pcmap` | `viridis` | colormap of the density |
| `--x-offset` | `0` | shift the distances so a known feature sits at *x* = 0 |
| `--keep-zero` | off | keep exactly-zero pixels (dropped by default, as in `plot_transection.py`) |
| `--decimate` | `1` | use every *N*th pixel of the corridor |
| `--nlook` | `1` | **oversample** the corridor *N*× per axis (*N*² more samples). The cure for a sparse profile from multilooked data; inverse of `--decimate` |
| `--dem` | — | DEM for the shaded-relief map background, as in view.py; **without a file, the Copernicus GLO-30 DEM is downloaded** for the map extent (see below) |
| `--topoxsect` | — | DEM to profile under the LOS panel, any format / any CRS; **without a file**, the `--dem` DEM or a download |
| `--dem-source` | `cop30` | where a download comes from: `cop30` (Copernicus GLO-30, 1″ ≈ 30 m) or `terrain` (the Terrain Tiles, also the fallback when GLO-30 cannot be had) |
| `--dem-zoom` | auto | `--dem-source terrain` only: zoom of the tiles, 0–15 (12 ≈ 38 m, 13 ≈ 19 m) |
| `--lookup` | auto | radar-coded data only: MintPy's `geometryRadar.h5` (found in `inputs/`, `./` or `../inputs/` when not given) |
| `--fault` | — | fault trace file(s), any vector format / any CRS: drawn on the map, and a dashed vertical line + "Fault line" on every profile at each crossing |
| `--fault-label` / `--fault-color` / `--fault-lw` | `Fault line` / `k` / `1.2` | the annotation text (`""` for none), the colour, and the map line width |
| `--svg` | off | write figures as vector `.svg` |
| `--figsize-map` | auto | size of the map window (`--figsize` sizes each section) |
| `--timeseries` | — | time-series file(s) for the 3rd panel |
| `--start-yx` / `--end-yx` | — | line end points in pixel numbers |
| `--start-lalo` / `--end-lalo` | — | line end points in lat/lon |
| `--profile LAT0 LON0 LAT1 LON1` | — | **one more section** up front; repeat it for `B-B'`, `C-C'` … (`--profile-yx` in pixels) |
| `--line-file` | — | sections from a file: GMT `lon lat` text, one section per segment (split by `>`), or a vector file (`.geojson .kml .shp .gpkg` …), one per line - each from its first to its last point |
| `--data-unit` | the file's own | the unit the data are really in, when the file says none (a GeoTIFF) or a wrong one - e.g. `mm/yr` for mm/yr values labelled m/year; a LiCSBAS velocity in its folder is read in mm/yr by itself, and a velocity without any unit (with `--gnss`, or a `vel*` name) in m/yr or mm/yr as its numbers say |
| `--data-scale` | `1` | multiply the data by a factor before anything - the map, the profiles, the fit, the statistics - one per file, or one for all (`0.1` for a product 10x too large). A rule in your `viewinsar.ini` does it for every file of a mission (*Your own defaults*, below) |
| `--save [DIR]` | — | write the sections given up front at start-up, and make <kbd>e</kbd> export into `DIR` (made if it does not exist) |
| `--nodisplay` / `--outfile` / `--outdir` | — | output control |

**Option names.** Since 0.5.0 every option of viewinsar is written with two
dashes. The spellings of before still work - with a note - so older commands
and scripts keep running:

| before 0.5.0 | now |
|---|---|
| `-fault`, `-gnss`, `-vel` / `--gnss-comp`, `-nlook`, `-topoxsect` / `--topo`, `-profile`, `-svg`, `-save` | `--fault`, `--gnss`, `--vel`, `--nlook`, `--topoxsect`, `--profile`, `--svg`, `--save` |
| `-d`, `-w` / `--buffer`, `-o` | `--dem`, `--width`, `--outfile` |
| `stat`: `-sp` | `--spacing` (0.5.8: `--radius`, `--gap`, `--rate`, `--rate-err` are gone - left out with a note) |
| `tseries`: `-r`, `-m`, `-p`, `-c`, `-u`, `--vel` | `--radius`, `--ref-date`, `--point`, `--colormap`, `--unit`, `--vel-file` |
| `list`: `-n` | `--num` |

view.py's own options for the map (`-v`, `-u`, `-c`, `-m` …) belong to view.py
and pass as they are; their long forms (`--vlim`, `--unit`, `--colormap`,
`--mask`) work with and without MintPy. Run `viewinsar swath -h` for the full
list.

### Topography: your DEM, or one downloaded for you

`--topoxsect dem.tif` profiles a DEM in any format and any CRS. Its CRS is read
from the file itself, the section is re-projected onto it, and it is sampled
at its own resolution, independently of `--nlook`. A DEM that covers part of
a section gives that part, and a section that gets no topography says why in
its panel.

Without a DEM file, `--dem` and `--topoxsect` fetch one:

| command | map background | elevation profile |
|---------|----------------|-------------------|
| `--dem` | downloaded | — |
| `--topoxsect` | — | downloaded |
| `--dem --topoxsect` | downloaded | the same download (finer than the background on a map over ~0.8°) |
| `--dem dem.tif --topoxsect` | `dem.tif` | `dem.tif` |

The download is the [Copernicus DEM GLO-30](https://registry.opendata.aws/copernicus-dem/)
on AWS Open Data: the 1-arc-second (~30 m) DEM of the TanDEM-X mission, with
heights above the EGM2008 geoid, as one Cloud Optimized GeoTIFF per 1 × 1°
cell, e.g. `Copernicus_DSM_COG_10_S01_00_E119_00_DEM` for lat −1…0,
lon 119…120. Only the parts over the map are read, by HTTP byte range, so
neither GDAL nor rasterio is needed. The tiles are joined seamlessly,
re-gridded onto the data's CRS (lat/lon or UTM) and saved as a MintPy geometry
`.h5` (`height`, metres) in `~/.cache/viewinsar` (or `$VIEWINSAR_CACHE`),
along with the pieces read. A second run over the same area downloads
nothing.

| what | grid |
|------|------|
| `--topoxsect` (the elevation profile) | 1″ pixels (30.9 m on a UTM grid), up to 8000 a side: full resolution over a map of up to ~2.2° |
| `--dem` (the map background) | the same, up to 3000 a side (~0.8°): view.py shades every pixel it is given |
| a bigger map | proportionally wider pixels, read from the tiles' own overviews (1/2, 1/4 and 1/8 resolution), so the download shrinks with them |

The sea has no tiles; it is 0 m, as the product specifies. So are the few
countries not released (Armenia, Azerbaijan). A tile that cannot be
downloaded leaves its area empty, and the run says so. The licence is free,
with the notice printed and written into the DEM file and `_info.txt`:
"produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus
Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European
Union and ESA; all rights reserved".

`--dem-source terrain` takes the [Terrain Tiles](https://registry.opendata.aws/terrain-tiles/)
instead, the "Mapzen Global Terrain" layer of QGIS:
`https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png`.
These are 256-px PNG tiles in Web Mercator, with the height coded in the colour
as R·256 + G + B/256 − 32768 m. The zoom is the finest that keeps the download
within about 300 tiles (at most 12, about 38 m); `--dem-zoom 13` with
`--sub-lat/--sub-lon` gives about 19 m over a small area. The sources are SRTM,
GMTED2010 and ETOPO1 (with bathymetry) worldwide, plus finer national DEMs
where available. They are also the fallback: when GLO-30 cannot be had (no
access to its server, say) and `--dem-source` was not given, the run says so
and takes the Terrain Tiles. Either way, please cite them: the attribution is
printed and written to `_info.txt`.

### Radar-coded data

A MintPy product that was never geocoded (view.py draws it in pixel rows and
columns, flipped left-right for a descending track) has no lon/lat of its
own. Its companion geometry file `geometryRadar.h5` has the latitude,
longitude and height of every pixel, on the same grid. viewinsar finds it the
way MintPy does: in `inputs/`, `./` or `../inputs/`, next to the data file or
in the current folder. You can also name it with `--lookup`. Then:

| command | what it uses |
|---------|--------------|
| `--dem`, `--topoxsect` (no file) | its `height`, the DEM of your own processing, already on the radar grid; nothing is downloaded |
| `--topoxsect dem.tif` | a geocoded DEM, reached through its latitude/longitude |
| `--fault faults.shp` | traces mapped from lon/lat onto the radar pixels; `_info.txt` gives each crossing's row/col and lat/lon |

Only a geometry file without `height` makes `--dem` download GLO-30 (or the
Terrain Tiles, `--dem-source terrain`), resampled onto the radar grid at about
its own pixel spacing.

### Fault traces

`--fault` takes one or more vector files and reads `.shp` (with its `.prj`),
`.geojson`, `.gpkg`, `.kml`/`.kmz` and GMT multi-segment text (`.gmt`,
`.lonlat`, `.xy`) directly, without GDAL. Anything else - FileGDB, GML,
FlatGeobuf, ... - goes through `pyogrio`, `fiona` or the GDAL bindings, whichever
is installed. Lines and polygon outlines count as traces; points are ignored.

Every file is re-projected to the data's CRS, so a WGS84 GeoJSON on a UTM
velocity map, or a UTM shapefile on a lat/lon one, just works. A file with no
CRS at all is taken as lon/lat when its numbers fit, and as the data's own CRS
otherwise (a UTM shapefile that lost its `.prj`).

A crossing is where the section's **centre line** meets a trace. Its distance
along the profile, and its coordinates, are written to `_info.txt`. Every
crossing gets its dashed line; when two are so close that their labels would
overlap, only the first is labelled.

### Several files and GNSS

**Several files.** `viewinsar swath FILE1 FILE2 ...` - the frames of a track,
or tracks side by side - samples the same line in every file, each on its own
grid and in its own CRS, and pools the samples into **one density**, with a
median line per file (`--labels` names them). The map is the first file, with
the others drawn under it in the same colours. Geocoded files may be in any
CRS; radar-coded files must share one grid. `_points.csv` says which file
each sample came from (`file_idx`).

**GNSS.** `--gnss FILE` reads GNSS velocities from a table (`.csv`, `.txt`,
`.xlsx`) or a point file (`.shp`, `.geojson`, `.gpkg`, `.kml`/`.kmz`, and
anything GDAL reads). The columns are found by name, in any spelling - `lon`,
`Longitude`, `Ve(mm/yr)`, `sig_e`, `Station` ...; a table without a header is
read as GMT `psvelo` (`lon lat ve vn se sn corr site`). The unit comes from the
headers (mm/yr, cm/yr, m/yr) or `--gnss-unit`, and is converted to the data's.

| `--vel` | the GNSS value compared with the data |
|--------|----------------------------------------|
| `los` (default) | projected into the line of sight of the first file, `ve·uₑ + vn·uₙ + vu·uᵤ`: the LOS unit vector from MintPy's geometry file (`geometryGeo.h5` / `geometryRadar.h5` beside the data), LiCSBAS's `GEOCml*/E.geo …` or LiCSAR's `*.geo.E.tif`, or `--geom`; else the data's `CENTER_INCIDENCE_ANGLE` and `HEADING`, or `--los-inc` / `--los-az`; else nominal angles for the orbit the file's name says (`asc`, `desc` …, see *model*). A station without `vu` counts `vu = 0` |
| `ve`, `vn`, `vu` | that component, for decomposed velocity maps |
| several, e.g. `--vel ve vn vu los` | all of them on the same profile, each in a colour of its own (ve magenta, vn orange, vu blue, los black), each shifted onto the InSAR on its own; `_gnss.csv` gets a row per station and component (`component` column). `vu` exists for 3D stations only. `model` takes one: the component it fits |

Every station is drawn on the map - a **triangle** for a 3D station (one with
`vu`), a **circle** for a 2D one, in black. The stations inside a swath (the
same flat-capped buffer as the pixels) are drawn on its profile, in black with
their 1-σ as error bars, as a comparison with the InSAR - outlined in white,
error bars too, so they stand out on the density behind them. By default
they are shifted onto the InSAR (`--gnss-ref insar`): by the median, over the
stations, of the InSAR within `--gnss-radius` km of a station minus the
station; `--gnss-ref none` leaves them as read. `--gnss-label` writes the site
names on the map.

---

## Output

Pressing <kbd>e</kbd> (or passing `--save`) writes four files per section,
named after its label (`swath_A_...`, `swath_B_...`), plus one
`swath_map_ABC...` figure of the map with every section line and fault trace.

**`swath_..._points.csv`** — one row per sample, the linear-referencing table:

| column | meaning |
|--------|---------|
| `dist_km` | distance along the profile from the start point, minus `--x-offset` |
| `perp_km` | distance from the centre line; **+ = left** of the start→end direction |
| `value` | the pixel value, in the display unit |
| `row`, `col` | pixel coordinates in the **full** file |
| `lat`, `lon` | pixel centre (geocoded files only) |
| `file_idx` | which input file the row came from |

**`swath_..._stats.csv`** — per along-profile bin: `count`, `median`, `mean`,
`std`, `q25`, `q75`, plus the density-weighted `median`/`q25`/`q75` that the
figure draws, and `elevation_median_m` (the exact median) when `--topoxsect` is used. The plain statistics are exact and un-clipped; the
density-weighted ones are clipped to the plotted Y range, exactly like the
histogram.

**`swath_..._gnss.csv`** — with `--gnss`, the stations in the swath: `site`,
`lon`, `lat`, `dist_km`, `perp_km`, `value`, `sigma`, `value_shifted` (as
drawn), `insar_near` and `is3d`.

**`swath_..._info.txt`** — everything needed to reproduce the profile, and with
`--fault` one `fault crossing : <km> km (dist_km), <lat/lon or N/E>` line per
crossing.
**`swath_....png`** — the figure.

### Referencing the profile to a fault, afterwards

By default *x* = 0 is the start of the line. To re-zero on a fault you already
know the position of:

```python
import pandas as pd

df = pd.read_csv('swath_A_Y300X150_Y300X750_w15km_points.csv')
df['dist_from_fault'] = df['dist_km'] - 27.5     # your fault position, in km
df['vel'] = df['value']                          # already in the display unit
```

With `--fault`, the position comes from `_info.txt`:

```python
import re

info = open('swath_A_Y300X150_Y300X750_w15km_info.txt').read()
fault_km = float(re.search(r'fault crossing : (\S+) km', info).group(1))  # the 1st
df['dist_from_fault'] = df['dist_km'] - fault_km
```

or do it up front with `--x-offset 27.5`, which also puts a dashed marker at
*x* = 0. The density grid always places a cell **boundary** at *x* = 0, so the
bins never straddle the fault.

---

## Fault creep along a line: `viewinsar stat`

How fast does a fault creep, and where along it? Draw a line **on the fault
trace**, and - as Hussain et al. (2016) measured the creep of the Izmit
rupture of the North Anatolian Fault -

1. every `--spacing` km along the line (default 5) a **point** is placed (a
   black circle on the map), centred on the line so the first and last sit
   equally far from its ends;
2. at every point a short **profile** crosses the line at right angles,
   `--length` km (5) either side of it, and the pixels within `--width` km
   (2.5) either side of the profile are projected onto it - on the map each
   profile is a rectangle in two halves, one per block;
3. **block 1** is the **left** of the line, walking from its first point to its
   last, **block 2** the right; *x* is the distance from the line along the
   profile, + towards block 2;
4. a **straight line through each side**, and the **creep rate is their offset
   at the line** - block 2 minus block 1 - solved as one linear system, their
   eq. 2:

   ```
   | v1 |   | x1  0  1  0 |   | slope 1 |
   |    | = |             | . | slope 2 |
   | v2 |   | 0  x2  1  1 |   | b1      |
                                | creep   |
   ```

5. its uncertainty from a **percentile bootstrap** that resamples each side on
   its own, *n₁* of its *n₁* pixels and *n₂* of *n₂* (`--boot` draws, 2000;
   after Bekaert et al. 2015): the creep rate is the mean of the draws, ±1σ
   their standard deviation. The least-squares solution (the BLUE) is written
   too, with its σ - as Hussain et al. show (their Fig. S3), it is too narrow
   when the errors are large or correlated.

```bash
# draw the line on the map (over the fault traces): left-click along it, Enter
viewinsar stat geo_velocity.h5 --fault faults.shp --colormap RdBu_r --vlim -1 1

# or give it: lat/lon or y/x pairs, or the longest line of a vector file
viewinsar stat geo_velocity.h5 --save --nodisplay \
    --line-lalo 40.70 29.90 40.73 30.40 40.70 30.70
viewinsar stat geo_velocity.h5 --line-file fault_trace.geojson --svg --save

# shorter, narrower profiles every 2 km
viewinsar stat geo_velocity.h5 --spacing 2 --length 4 --width 1.5 --line-file fault_trace.geojson

# an ascending and a descending track: the LOS creep of each, and the split
viewinsar stat T036D/geo_velocity.h5 T131A/geo_velocity.h5 --labels T036D T131A \
    --line-file fault_trace.geojson
viewinsar stat desc.tif asc.tif --los-inc 39 34 --los-az -100 100 --data-unit mm/yr \
    --line-file fault_trace.geojson

# your own names and colours for the blocks
viewinsar stat geo_velocity.h5 --names North South --colors tab:green tab:purple \
    --line-file fault_trace.geojson
```

**The figure** is their Fig. 5: the **LOS creep rate** along the line in
mm/yr, 1σ error bars, a series per file - descending tracks solid, ascending
dashed (the file's `ORBIT_DIRECTION`, else its name - `asc`, `T014A`, `105D`
... - else its line of sight). A file without a unit (a GeoTIFF) is taken as a
velocity: mm/yr, or m/yr when 99 % of its values are below 1 (the run says
which; `--data-unit` if not). Motion along an E-W fault shows with opposite
signs in ascending and descending tracks, vertical motion with the same sign;
so **where an ascending and a descending track both have a creep rate and an
LOS geometry**, the creep is also **split**, point by point and draw by draw,
into a **fault-parallel** creep rate (red, + right-lateral) and a **vertical**
one (blue, + block 2 up), assuming no motion across the fault, as Hussain et
al. did:

  `los_i = Δh (L_i · s) + Δu L_i,up`, fault-parallel = −Δh, vertical = Δu

with `L_i` the LOS unit vector of file *i* (ground to satellite, the median
over the profile's pixels) and `s` the line's direction there. Elsewhere the
run says why not, range by range - one track only, a track without an LOS
geometry, a geometry that does not reach the point, tracks of one pass only,
lines of sight too alike - and so do a note in the panel and the `why_not`
column of `_split.csv`. Along a fault that runs N-S the fault-parallel rate is
less certain (the satellites look E-W): its error bars show it.

**Each track its own LOS geometry.** MintPy's `geometryGeo.h5` /
`geometryRadar.h5`, LiCSBAS's `GEOCml*/E.geo …` or LiCSAR's / your own E, N, U
GeoTIFFs beside the data, or `--geom`, or `--los-inc` / `--los-az`. When the
names say the track - `039A.los.tif` and `039A.E.tif`, `039A_*.geo.E.tif` - a
file takes its own track's set only, never another track's from the same
folder, and the frames of a track mosaic (`039A_09001_…`, `039A_09202_…`) are
joined into one; a geometry of the other pass (ascending for a descending
file) is never taken. `--geom` gives one geometry per file, in order - or
FOLDERS, fewer than the files (`--geom geo_asc/ geo_dsc/` for four tracks),
all searched for every file. Nominal angles from a name are not used for the
split. The run lists what each file got:

```
LOS geometry - for the split into fault-parallel and vertical:
   039A.los       ascending   …/039A.E.tif (+ N, U; LiCSAR)
   112A.los       ascending   …/112A.E.tif (+ N, U; LiCSAR)
   061D.los       descending  …/061D.E.tif (+ N, U; LiCSAR)
   134D.los       descending  …/134D.E.tif (+ N, U; LiCSAR)
```

**Click a point of the creep rate**: its profile opens in a window of its own,
as their Fig. S3 - the pixels, the BLUE lines (red) and the bootstrap lines
(green), the fault location, and the PDFs of the creep rate from both, the
numbers in the legend; the profile is outlined on the map.
<kbd>←</kbd> / <kbd>→</kbd> step to the point before / after,
<kbd>↑</kbd> / <kbd>↓</kbd> to the other files.

**Mind the line.** It must lie on the fault trace: *x* = 0 is where the two
sides are compared, and a line 1 km off the trace puts that km of one side
into the other - the creep comes out low. Draw it over the `--fault` traces,
or take the trace itself with `--line-file`. The line is used where you click
or give it, not moved to the centre of a pixel (a pixel can be a km wide). A profile is perpendicular to the
line's direction over ±`--width` around its point, so a hand-digitised trace
with small wiggles still gives steady profiles; where the line bends inside a
profile, a pixel on its other side than its *x* says is left out (the run says
how many; `trace_km` in `_creep.csv` is how far the line strays from the
centre of the profile). Straight lines through ±5 km also see the curvature of
the interseismic signal - small while the locking depth is well beyond 5 km; a
shorter `--length` lessens it, with fewer pixels. Neighbouring profiles share
pixels when `--spacing` is less than twice `--width`; the run notes it.

"Left" and "right" are those of the ground, as the map shows it, also on
radar-coded data, whose grid is a mirror image of the ground (that is why
view.py flips it): there, `geometryRadar.h5` gives the lat/lon of the line
(`--line-lalo`, `--line-file`) and of every pixel in the `.csv` files.

| input | action |
|-------|--------|
| left click (map) | add a point to the line |
| <kbd>Enter</kbd> | finish the line: the profiles on the map, the creep rate |
| <kbd>u</kbd> or <kbd>Backspace</kbd> | remove the last point |
| <kbd>n</kbd> | start a new line (so does a click after <kbd>Enter</kbd>) |
| left click (creep rate) | the profile of that point, in its own window |
| <kbd>←</kbd> <kbd>→</kbd> / <kbd>↑</kbd> <kbd>↓</kbd> (profile window) | the point before / after; the other files |
| <kbd>e</kbd> | export |

| option | default | meaning |
|--------|---------|---------|
| `--spacing` | `5` | distance between the points along the line, km - a profile at each |
| `--length` | `5` | how far each profile reaches either side of the line, km |
| `--width` | `2.5` | the pixels within this distance either side of a profile are projected onto it, km |
| `--boot` / `--seed` | `2000` / `0` | bootstrap draws per profile, and their seed |
| `--min-pixels` | `10` | fewer pixels on a side: no creep rate at that point |
| `--geom` / `--los-inc` / `--los-az` | found beside the data | the LOS geometry of each file, for the split: one per file, or folders searched for every file (above) |
| `--line-yx` / `--line-lalo` | — | the line as row/col or lat/lon pairs, 2 points or more |
| `--line-file` | — | the longest line of a vector file (any format / CRS `--fault` reads) |
| `--names` / `--colors` | `Block 1` `Block 2` / blue, vermilion | the two blocks |
| `--ylim` | auto | range of the creep-rate axes |
| `--keep-zero` | off | keep exactly-zero pixels |
| `--dem`, `--dem-source`, `--dem-zoom`, `--lookup`, `--fault` | — | as in `viewinsar swath` |
| `--save` / `--nodisplay` / `--svg` / `--outfile` / `--outdir` | — | output control |

The rates are in mm/yr whatever unit the data are read in (a file that is not
a velocity keeps its own). Everything else goes to view.py, as for `swath`.
The export (<kbd>e</kbd> or `--save`) is named after the line's first and last
point and the parameters, e.g. `stat_Y3670X744_Y2967X696_sp5km_L5km_w2.5km`:

| file | contents |
|------|----------|
| `….png` | the creep rate along the line (and the split) |
| `…_map.png` | the map with the line, its points and their profiles |
| `…_creep.csv` | a row per point and file: `dist_km`, row/col, lat/lon, `strike_deg` (the line's direction), `label`, `orbit`, `n_1` `n_2` (pixels a side), `creep` (the bootstrap mean) `creep_std` `creep_median` `creep_lo95` `creep_hi95` `creep_lo68` `creep_hi68`, `creep_blue` `creep_blue_std`, `slope_1` `slope_2` (per km) `intercept_1`, `sigma` (the scatter), `trace_km`, `dropped`, `clipped`; with the split, the LOS unit vector `los_east` `los_north` `los_up` |
| `…_split.csv` | two LOS geometries or more: a row per point, `par` / `up` (from the least-squares rates) and `_mean` `_std` `_median` `_lo95` `_hi95` of the draws, the files used, and `why_not` where there is no split |
| `…_values.csv` | every pixel used: `file`, `point`, `block`, `x_km` (from the line, + towards block 2), `y_km` (along it), `value`, `row`, `col`, lat/lon |
| `…_profiles.png` | every profile of the first file - pixels and both lines, its creep rate (`_profiles_2.png` … for the others) |
| `…_profile7.png` | the profile on screen (point 7), when its window is open |
| `…_info.txt` | the line's vertices, the parameters, the method, which side is which block |

`model --creep-offset` measures the same offset, with the same code, on the
profile of a model fit.

---

## Fault slip models: `viewinsar model`

Three kinds of fault, named by the first word:

```bash
viewinsar model strike geo_ve.h5 --comp ve ...                 # strike-slip (the default)
viewinsar model reverse geo_vu.h5 --comp vu ...                # reverse (thrust)
viewinsar model normal geo_ve.h5 geo_vu.h5 --comp ve vu ...    # normal
```

A **strike-slip** fault (`strike`, or no word) is fitted to the velocity along
it; a **reverse** or **normal** one to the vertical velocity, the motion
across it or a line of sight - see [Reverse and normal
faults](#reverse-and-normal-faults) below. The rest of this section holds for
all three.

The interseismic velocity across a strike-slip fault, after Savage & Burford
(1973) and Hussain et al. (2016, 2018):

```
v(x) = S/π · atan(x/d₁) + a
     + C · [H(x) − ½ − atan(x/d₂)/π]     shallow creep (--creep), d₂ ≤ d₁
     + θ · x, on one side of the fault    rotation (--rot), fixed
```

*x* is the distance from the fault, perpendicular to it; *S* the slip rate
below the locking depth *d₁*; *C* the creep rate from the surface down to *d₂*,
in the same sense as *S*; *a* an offset.

Draw a profile **across** the fault (two clicks), or give it - or several,
`A-A'`, `B-B'`, `C-C'` … - up front with `--start-lalo` / `--end-lalo`,
`--profile` (repeat it) or `--line-file` (one per line); each section is
fitted on its own, with its own windows and files, and
`model_sections_ABC.csv` gathers them. Every pixel of the swath becomes a point:
velocity against distance from the fault. The points are binned (`--bin`,
1 km: the median of each bin, as Ou et al. 2022) and the bins fitted with the
affine-invariant ensemble MCMC sampler of Goodman & Weare (2010) - the
sampler of `emcee`, which Hussain et al. used - under uniform priors, with a
noise covariance that decays exponentially with distance, `σ² exp(−r/λ)`. The
result: the **maximum a posteriori (MAP)** model and the **95 % range** of
every parameter and of the curve.

```bash
# click a profile across the fault: the trace gives the fault position and its
# strike, and the LOS velocity becomes fault-parallel velocity, pixel by pixel
viewinsar model geo_velocity.h5 --comp los --fault faults.shp --width 25

# up front, saved, nothing shown; with shallow creep
viewinsar model geo_velocity.h5 --comp los --fault faults.shp --creep \
    --start-lalo 40.20 33.10 --end-lalo 41.30 33.40 --width 25 --nodisplay

# three profiles along the fault, each fitted: the slip rate along strike
viewinsar model geo_velocity.h5 --comp los --fault faults.shp --width 25 \
    --profile 40.20 33.10 41.30 33.40 --profile 40.30 33.90 41.40 34.20 \
    --profile 40.40 34.70 41.50 35.00 --nodisplay --save fits

# LiCSBAS as it is - cum_filt.h5 (its velocity), results/vel.filt.mskd or a
# GeoTIFF of it: the grid, the unit (mm/yr), results/mask and the line of sight
# of every pixel (GEOCml10GACOS/E.geo, N.geo, U.geo) are found
viewinsar model TS_GEOCml10GACOS/cum_filt.h5 --fault faults.shp --gnss gnss.csv \
    --width 20 --strain

# ascending and descending together (each with its own offset), and GNSS
viewinsar model T014A/geo_velocity.h5 T123D/geo_velocity.h5 --comp los \
    --fault faults.shp --gnss gnss.csv --labels T014A T123D --strain

# as Hussain et al. (2016): the rotation south of the fault, a prior on d1,
# the noise covariance of a far-field area
viewinsar model geo_vpar.h5 --comp par --strike 100 --rot 0.0215 \
    --d1-prior 17 5 --noise 6.35 35.8
```

**A line of profiles, as `viewinsar stat`.** With `--spacing KM` the
profiles are not drawn one by one: draw a line **along** the fault instead -
left-click its points, <kbd>Enter</kbd> - or give it with `--line-lalo`,
`--line-yx` or `--line-file` (with `--spacing`, the longest line of the file:
the trace itself). Every `--spacing` km along it (centred on it) a profile
crosses it at right angles - `--length` km either side (default 50), its swath
`--width` km either side - and the fault of each profile is where it crosses
the line, its strike the line's. Each profile is a section of its own, A-A',
B-B' … with its window, its fit and its files, as above; their posteriors are
written by the export rather than opened, one window each would be too many.
Then a window **along the line**: the slip rate and locking depth of every
profile (and *C*, *d₂* with `--creep`; the dip of a reverse or normal fault),
the MAP and its 95 % range - the slip rate along strike.

```bash
# click the line along the fault, Enter: a profile every 20 km, +/-60 km long, +/-15 km wide
viewinsar model geo_velocity.h5 --comp los --fault faults.shp --spacing 20 \
    --length 60 --width 15

# the trace as the line, two tracks fitted together, all saved
viewinsar model T014A/geo_velocity.h5 T123D/geo_velocity.h5 --comp los \
    --spacing 25 --length 50 --line-file fault_trace.geojson --nodisplay --save fits
```

With `--spacing`: <kbd>u</kbd> / <kbd>Backspace</kbd> removes the last point
of the line, <kbd>n</kbd> starts a new one (its profiles close), <kbd>e</kbd>
exports - every section, `model_sections_ABC….csv` (with each profile's
`point` and `dist_km` along the line), `model_line_…_along.png` and
`model_line_…_info.txt`. The fits take a while: a few seconds to a minute each
with the default sampler (`--steps`, `--walkers`). Profiles closer than twice
`--width` share pixels; the run says so.

**The fault and the data.** The fault is where the profile crosses a `--fault`
trace (the crossing nearest the middle of the profile), with the strike of the
trace there - or `--x0` (km from the start of the profile) and `--strike`.
Distances are taken perpendicular to the fault, so a profile drawn at an angle
to it still gives the right locking depth. `--comp` says what the data are:

| `--comp` | the velocity fitted |
|---|---|
| `auto` (default) | `los` when the line of sight of every file is known (below; nominal angles do not count) and no file is named as a component of a decomposition (`hz`, `east`, `ve` …); `as-is` otherwise, with the GNSS in the LOS (`--vel los`) - the run says which. A file named as the vertical (`up`, `vu`, `vertical`) is a reverse or normal fault's: the run says so |
| `los` | LOS velocity, divided by (LOS unit vector · strike direction), pixel by pixel |
| `par` | fault-parallel velocity, positive along `--strike` |
| `ve`, `vn` | east / north velocity, divided by sin / cos of the strike |
| `as-is` | the data as they are - e.g. the LOS itself; *S* is then the rate seen in the LOS ("slip rate in the LOS" for LOS data) |

One `--comp` for all the files: they all become fault-parallel velocity. The
y axis says what is plotted: `Fault-parallel velocity (mm/yr)` for `par`,
`ve`, `vn` and `los` - all turned into it - and with `as-is` the data as they
are: `LOS velocity (mm/yr)`, or a component by the name of its file (`East
velocity` for `geo_ve.h5` …). `--comp vu` - the vertical - is for a reverse or
a normal fault.

**The line of sight** of every pixel is found next to the data - nothing to
type:

| from | found |
|---|---|
| MintPy | `geometryGeo.h5` / `geometryRadar.h5` (`incidenceAngle`, `azimuthAngle`) beside the file, in `inputs/`, `geo/`, `mintpy/inputs/` or one or two folders up |
| LiCSBAS | `GEOCml*/E.geo`, `N.geo`, `U.geo` with `EQA.dem_par`, for a file anywhere under `TS_GEOCml*/` (e.g. `TS_GEOCml10GACOS/results/vel.filt.mskd.geo.tif` → `GEOCml10GACOS/`) |
| LiCSAR | `{frame}.geo.E.tif`, `.N.tif`, `.U.tif` beside the file, in `metadata/` or in LiCSBAS's `GEOC/` |
| the data | `CENTER_INCIDENCE_ANGLE` and `HEADING`: one value for the whole scene |
| the name | **nominal angles**, as a last resort, when the file's name (or its folder's, or its `ORBIT_DIRECTION`) says the orbit - `asc`, `desc`, `scansardesc`, LiCSAR's `105D_…`, `T125A` - and, if it can, the mission (`alos2`, `scansar`, `s1` …): mid-swath incidence (Sentinel-1 IW 37.5°, ALOS-2 38°) and the heading of the orbit at the scene's latitude, right-looking. One value for the scene, and said so: the incidence changes by 15-25° across a swath. Used for the GNSS in the LOS, never to turn the InSAR into fault-parallel velocity |

`--geom` names it when it is elsewhere, one per input file: any of those
files - also E / N / U GeoTIFFs named `E.geo.tif`, `E.tif`, `*_east.tif` …,
and an incidence + azimuth GeoTIFF pair (`incidenceAngle.tif` +
`azimuthAngle.tif`, `inc.tif` + `az.tif` …) - or **a folder**: it is searched,
with its subfolders (LiCSBAS's `GEOCml*/`, `GEOC/`, `metadata/` first; the
interferogram folders of `GEOC/` skipped), and the first geometry that covers
the data is taken - so `--geom` the LiCSBAS work folder, or a folder of
geometry GeoTIFFs, is enough. Folders fewer than the input files (`--geom
geo_asc/ geo_dsc/` for four tracks) are all searched for every file. Found or
searched, a file whose name says its track (`039A.los.tif`, `105D_09094_…`,
`T14A`) takes that track's geometry only - `039A.E.tif`, `039A_*.geo.E.tif`,
the frames of a track mosaic joined into one - never another track's from the
same folder, and never one of the other pass. `--los-inc` and `--los-az` give one
pair of angles instead (one per file, or one for all; MintPy's convention:
the azimuth of the ground-to-satellite vector, anticlockwise from north -
about 100 for Sentinel-1 ascending, −100 descending).

**LiCSBAS's own files** open as they are, in `swath`, `stat` and `model`:
`TS_GEOCml*/cum_filt.h5` or `cum.h5` (the velocity `vel` by default;
`--dset cum-20200101` an epoch of the cumulative displacement; `viewinsar
list cum_filt.h5` lists them) and the float32 files of `results/` (`vel`,
`vel.filt.mskd`, `vstd` …), placed on the grid of `GEOCml*/EQA.dem_par`.
Their unit (mm/yr, mm for `cum`), LiCSBAS's `results/mask` and the line of
sight are found too.

**The unit.** Every file is read in the unit it says (its `UNIT`), and shown
in view.py's display unit (cm/year for a MintPy velocity). A LiCSBAS
velocity (`vel*`, `vstd*`, also as the GeoTIFF that `LiCSBAS_flt2geotiff.py`
writes) is in mm/yr but does not say so; inside its LiCSBAS folder (or with
its `GEOCml*` as `--geom`) viewinsar reads it in mm/yr. **A velocity that
says no unit at all** (a GeoTIFF; in `model`, with `--gnss`, or with a `vel*`
name) is read in m/yr when 99% of its values are below 1 - a deformation map
in m/yr hardly reaches 1, the noise alone of one in mm/yr exceeds it - and in
mm/yr otherwise; it is then shown in mm/yr, like the GNSS. **`--data-unit`**
gives the unit a file is really in - when it says none, or a wrong one
(e.g. mm/yr values in a file labelled m/year or cm/year, so that the InSAR
looks 10x or 1000x larger than the GNSS): `--data-unit mm/yr` rescales the
map, the profiles, the fit and the comparison with the GNSS. The console
says which unit each file is read in.

`ve`, `vn` and `los` assume the motion is purely fault-parallel (no vertical,
no fault-normal motion), as Hussain et al. did; a LOS nearly perpendicular to
the fault sees little of that motion, and the run warns when it sees less than
a quarter. With them the strike's direction is chosen so that **S > 0 is
right-lateral**, and the sense of slip is written out.

**Several files** (tracks, frames) are fitted together, each with **an offset
of its own** (`a`, `a2`, `a3` ...), so each may be referenced its own way -
an ascending and a descending track each with its own line of sight.
**GNSS** (`--gnss`) within `--gnss-width` km of the line are fitted with them,
with their 1-σ and an offset of their own (`ag`): their fault-parallel
component, from `ve` and `vn` (`--vel` picks another one, to match
`--comp as-is`). `--gnss-nofit` shows them without fitting them.

**The figures**, after Hussain et al.: the pixels as points (red; a colour per
file), the median of each bin with ±1σ (its robust spread), the GNSS in blue
with their error bars (outlined in white), the MAP model as a bold dashed line
in its **95 %** and **68 %** ranges (grey, 25 % and 40 % transparent - 68 %
being ±1σ; `--bands 95 65` shades others), the fault dotted, `A` and `A'` in
the top corners, and the parameters written in a box as Hussain et al. write
them - `slip rate = 20.3 mm/yr (18.4 – 34.0, 95% CI), right-lateral`,
`locking depth = 12.7 km (…)`, `creep rate` in mm/yr, `creep depth` in km,
and the `strain rate` at the fault in µstrain/yr - placed where they hide the
fewest points; the title gives the swath, the pixels and the strike. The
slip rate is given as a rate with its sense; when its 95 % range crosses
zero it keeps its sign (S > 0 right-lateral). Every file is shifted onto the
offset of the first, the GNSS onto the InSAR. `--strain` adds the shear
strain rate of the model, `dv/dx`, below; `--topoxsect` the topography.

The **legend** goes inside the plot where it hides the least, room made for
it; when the plot is too full for it, it goes **outside, on the right**
(beyond the colour bar of a density), and the figure grows to hold it, so
the plot keeps its size. `--legend in` keeps it inside anyway, `--legend
out` puts it outside always, `--legend off` leaves it out; on the profile,
<kbd>l</kbd> moves it - inside, outside, off, in turn - without fitting
again.

A second window holds the **posterior**, drawn as Hussain et al. (2016,
central NAF, Fig. 9): every pair of parameters (*S*, *d₁*, the offset *a*, and *C*, *d₂*,
*x₀* when fitted) as a 2D histogram in `jet`, each scaled to its own peak, the
MAP as a red dot; along the bottom the histogram of each parameter in blue,
the MAP as a red line, the 95 % range dashed and written above it; the
section, `A - A'`, white on black in the corner. The trade-off between *S* and
*d₁* shows there.

**Two views of the pixels.** By default they are points, as above. With
`--cmap-prof CMAP` (alone: `viridis`) they are a **pixel density** instead,
drawn as `viewinsar swath` draws its own - log colour scale, the densest
cells white, a colour bar, cells half a `--bin` wide - with the 95 % and 68 %
ranges in **red** over it, the medians **bigger** (white, a colour per file
with several), the **GNSS red**, outlined in white, and no points:

```bash
viewinsar model TS_GEOCml10GACOS/cum_filt.h5 --fault faults.shp --width 20 \
    --cmap-prof magma --nlook 10
```

`--nlook N` fills in the picture of a multilooked raster, as in swath - the
points or the density - with N × N samples per pixel. **The fit always takes
the real pixels**: interpolated samples are no new data. Before 0.5.5 they went
into the fit, which made the 95 % ranges too narrow (in a test, half as wide
with `--nlook 3`). The title says when the drawing is oversampled, and so does
`_params.txt`.

**Creep: its rate and its depth.** Hussain et al. (2016) constrained the
creep two ways, and `model` does both:

- **In the model, `--creep`**: the creep rate *C* and the depth it reaches,
  *d₂*, are fitted with *S* and *d₁* - creep in the sense of the slip, from
  the surface down to *d₂* ≤ *d₁*, under uniform priors (0 ≤ *C* ≤ 30 mm/yr,
  *d₂* ≤ 40 km: `--c-max`, `--d2-max`; Hussain et al.: 0-15 and 0-40 km at
  Izmit, 0-30 and 0-40 km in the centre). The profile says
  `creep rate = 8.1 mm/yr (7.2 – 8.8, 95% CI)` and `creep depth = 5.5 km
  (…)`, and draws the MAP model without the creep, dotted. The depth is the
  weakly resolved one: Hussain et al. found 6-9 km at the MAP but 95 % ranges
  of ~0-25 km, and creep down to 4, 8 or 12 km fitted Izmit within the noise.
  A prior on *d₁* (`--d1-prior 17 5`, theirs at Izmit) narrows it a little.
- **Straight from the data, `--creep-offset [KM [WIDTH]]`**: a straight line
  through the pixels on each side of the fault, within KM (5) of it and WIDTH
  (2.5) km of the profile line, each with its own slope; their offset at the
  fault is the creep rate, its 95 % range from 2000 bootstrap draws that
  resample each side on its own (Hussain et al. after Bekaert et al. 2015) -
  the method, and the code, of `viewinsar stat`, at one place.
  Drawn as two orange lines either side of the fault, and written in the box.
  With an **ascending and a descending track** the LOS offsets are also
  **split into fault-parallel and vertical**, by least squares and draw by
  draw: Hussain et al. did this at Izmit, where the subsidence of the
  Adapazari basin added to the LOS step and had made an earlier creep rate
  too high. The box then gives `creep offset, fault-parallel = …` and
  `vertical offset = +1.6 mm/yr (…), A' side up`.

```bash
viewinsar model T014A/geo_velocity.h5 T123D/geo_velocity.h5 --comp los \
    --fault faults.shp --labels T014A T123D --creep --creep-offset
```

The offset is model-free, but mind three things. The pixels are taken as
independent, so with spatially correlated noise its range is too narrow. The
lines are straight, while creep that stops a few km down dies away within a
few km of the fault: with *d₂* < KM the offset comes out below *C* (6 mm/yr
down to 4 km gave ~4 mm/yr with ±5 km lines in a test) - a shorter KM
lessens that, with fewer pixels. And the split assumes no fault-normal
motion, with the line of sight of each track its median near the fault.
`--creep-offset` works with or without `--creep`; the two agree when the
creep reaches deeper than KM.

### Reverse and normal faults

`viewinsar model reverse` (or `thrust`) and `viewinsar model normal` fit a
fault that **dips**: locked from the surface down to *d₁*, slipping at *S*
**along its dip** below it - up the dip for a reverse fault, down it for a
normal one. This is the back-slip model of Savage (1983), a 2D edge
dislocation as Segall (2010, eqs. 3.71-3.72) writes it (checked against
Okada 1992 in the tests). With *u* = *t* (*x* − *x₀*) the distance from the
trace towards the dip and *z* = *u*/*d₁* − cot δ:

```
across the fault (+ towards A'):  h(x) = −t·s·S/π · [cos δ · atan z + (sin δ − z cos δ)/(1 + z²)] + a
vertical (+ up):                   v(x) =    s·S/π · [sin δ · atan z + (cos δ + z sin δ)/(1 + z²)] + a
a line of sight:                   ph · h(x) + pu · v(x) + a
```

*s* = +1 for a reverse fault, −1 for a normal one; *t* = +1 when the fault
dips towards `A'` (the end of the profile), −1 towards `A`. Far off, the
hanging wall moves *S* cos δ across the fault (**shortening**, or
**extension**) and *S* sin δ up (**uplift**) or down (**subsidence**). The
free parameters are *S* ≥ 0, *d₁*, the dip δ (and *x₀* with `--free-x0`), an
offset per file and the GNSS's.

```bash
# a reverse fault, from the vertical velocity of an asc + desc decomposition:
# the dip free between 10 and 45 degrees, dipping towards A' (the end of the line)
viewinsar model reverse geo_vu.h5 --comp vu --fault faults.shp --dip 10 45

# a normal fault, from the east and the vertical velocity together - a panel
# each - dipping towards A; the GNSS across the fault
viewinsar model normal geo_ve.h5 geo_vu.h5 --comp ve vu --fault faults.shp \
    --dip-toward start --gnss gnss.csv

# a thrust seen by two tracks, each LOS as it is, the dip fixed at 25 degrees
viewinsar model thrust T014A/geo_velocity.h5 T123D/geo_velocity.h5 --comp los \
    --fault faults.shp --labels T014A T123D --dip 25
```

**The data**, one `--comp` per file (or one for all):

| `--comp` | the velocity fitted | its panel |
|---|---|---|
| `vu` | the vertical velocity, as it is | `Fault vertical velocity (mm/yr)` |
| `ve`, `vn` | the east / north velocity divided by the east / north component of the fault normal: the motion across the fault, + towards `A'`, assuming none along it (the run warns when the component sees less than a quarter of it) | `Fault-normal velocity (mm/yr)` |
| `los` | a line of sight as it is, the model seen through it - *ph* and *pu*, what the LOS sees of the horizontal across the fault and of the vertical, are the medians over the swath; the line of sight must be known (not nominal angles) | `LOS velocity, T014A (mm/yr)`, one per file |
| `auto` (default) | by the file's name: `up` / `vu` / `vertical` → `vu`, `east` / `ve` → `ve`, `north` / `vn` → `vn`; else `los` | |

Each component gets **a panel of its own**, stacked over one distance axis,
with the model in each - the vertical and the horizontal of one fault, fitted
together, pin the dip down; one component alone leaves *S* and the dip trading
off (the posterior shows by how much). The box on the top panel gives
`slip rate = 8.2 mm/yr (7.6 – 9.5, 95% CI), reverse`, `locking depth`, `dip =
29.9° (27.4 – 33.0, 95% CI), towards A'` (`(fixed)` with one `--dip`), then
`shortening rate` and `uplift rate` - `extension rate` and `subsidence rate`
for a normal fault. `A` and `A'` and the title are over the top panel, the
fault dotted in each.

**The dip.** `--dip 30` fixes it; `--dip 10 45` is the range of its uniform
prior (default 5 85). `--dip-toward end` (the default) has the fault dip
towards `A'` - the hanging wall on that side - and `start` towards `A`: draw
the profile, or choose the word, so it matches the fault. The strike is then
given by the right-hand rule (the fault dips to its right).

**GNSS** of a dipping fault are `--vel perp` - across the fault, from `ve`
and `vn` (the default when a file is `ve` / `vn`, or all are LOS) - or
`--vel vu` (the default when the files are all `vu`). When no file has their
component they get a panel of their own, with the model: GNSS horizontals
with InSAR verticals constrain the dip.

What belongs to a strike-slip fault - `--creep`, `--creep-offset`, `--rot`,
`--strain`, `--c-max`, `--d2-max` - is refused here, and `--dip` /
`--dip-toward` with a strike-slip fault. `_params.csv` gives `h_rate` and
`v_rate` - *S* cos δ and *S* sin δ - instead of the strain rate (and the dip
when it is fixed); `_samples.csv` both per sample; `model_sections_ABC.csv`
the fault type, *S*, *d₁*, the dip, *x₀* and both rates.

### Options

| option | default | meaning |
|---|---|---|
| `strike`, `reverse` (`thrust`), `normal` | `strike` | the fault type, the first word |
| `--width` | `20` | half width of the swath, km |
| `--spacing` | off | a line of profiles along the fault (above): a profile every KM |
| `--length` | `50` | with `--spacing`: how far each profile reaches either side of the line, km |
| `--line-lalo` / `--line-yx` / `--line-file` | — | with `--spacing`: the line along the fault (without it, `--line-file` gives the sections) |
| `--bin` | `1` | bin size along the profile, km |
| `--comp` | `auto` | what the data are (above): one for all (strike), or one per file (reverse, normal) |
| `--geom`, `--los-inc`, `--los-az` | found | the line of sight (above), one per file |
| `--data-unit` | the file's own; LiCSBAS: mm/yr | the unit the files are really in, when they say none or a wrong one |
| `--x0`, `--strike` | from the trace | the fault position (km along the profile) and its strike |
| `--dip DEG [DEG]` | `5 85` | reverse / normal: the dip, fixed (one value) or the range of its prior |
| `--dip-toward` | `end` | reverse / normal: the fault dips towards `A'` (`end`) or `A` (`start`) |
| `--creep` | off | fit *C* and *d₂* too |
| `--creep-offset [KM [WIDTH]]` | off; 5, 2.5 | the creep rate from the offset of two lines either side of the fault (above) |
| `--free-x0 [KM]` | off | fit the fault position within ±KM (default 5) |
| `--rot RATE`, `--rot-side` | none, `end` | a fixed linear term on one side, mm/yr per km |
| `--s-max`, `--d1-max`, `--c-max`, `--d2-max` | 60, 60, 30, 40 | the uniform priors: \|S\| ≤ (reverse, normal: 0 ≤ S ≤), 0.1 ≤ d₁ ≤, 0 ≤ C ≤, d₂ ≤ |
| `--d1-prior MEAN SIGMA` | none | a Gaussian prior on *d₁* (Hussain et al. 2016, Izmit: 17 5) |
| `--noise` | `auto` | `auto`: exponential, fitted to the residuals of a least-squares fit, for each file; `white`: each bin's own scatter; `SILL LENGTH`: given, in (mm/yr)² and km |
| `--walkers`, `--steps`, `--burn`, `--seed` | 64, 5000, a third, 0 | the sampler |
| `--gnss`, `--vel`, `--gnss-width`, `--gnss-nofit` | - | GNSS (`--vel perp` / `vu` for reverse and normal faults) |
| `--cmap-prof [CMAP]` | off: points | the pixels as a pixel density in CMAP (above) |
| `--legend [WHERE]` | `auto` | the legend: `auto` inside where it hides nothing, else outside on the right; `in`, `out`, `off` (<kbd>l</kbd> on the profile) |
| `--bands OUTER INNER` | `95 68` | the two ranges of the model shaded, % |
| `--nlook` | `1` | oversample the picture (points or density), not the fit |
| `--strain`, `--ylim`, `--max-points` | off, auto, 200000 | the figure |

The export (<kbd>e</kbd> or `--save [DIR]`, into `DIR` - made if it does not
exist), named after the line, e.g. `model_A_Y470X150_Y130X700_w25km`:

| file | contents |
|---|---|
| `….png` | the fit |
| `…_posterior.png` | the posterior distributions |
| `…_params.csv` | every parameter: `map`, `median`, `lower_95`, `upper_95`; and the strain rate at the fault (reverse, normal: `h_rate`, `v_rate`); with `--creep-offset`, `C_offset` (per file `C_offset_1`, `C_offset_2` …; split, `C_offset_par` and `U_offset`): the offset of the lines as `map`, the bootstrap's median and 95 % range |
| `…_params.txt` | the same, readable, with the sense of slip, the fault position and strike, the noise, the priors, and the sampler and its diagnostics |
| `…_samples.csv` | posterior samples (up to 20 000), for plots of your own |
| `…_bins.csv` | the bins: `x_km`, `median`, `std`, quartiles, `count`, the MAP model at each and the residual |
| `…_points.csv` | every pixel: `dist_km` along the profile, `x_km` from the fault, `perp_km`, `value` as read, `v_model` (the velocity fitted), `v_shifted` (as drawn), row / col, lon / lat, `file_idx` |
| `…_gnss.csv` | the GNSS in the swath, and whether they were fitted |
| `model_map_AB.png` | the map with every section, the fault position of each marked ◇ |
| `model_sections_AB.csv` | two sections or more: one row each - its ends (`start_lat` … `end_lon`), where it crosses the fault (`fault_km`, `fault_lat`, `fault_lon`), `strike_deg`, `sense`, `pixels`, and `S`, `d1` (`C`, `d2`, `x0`), `strain` (and the creep offsets) as `_map`, `_lower_95`, `_upper_95`; reverse, normal: `fault_type`, `S`, `d1`, `dip`, `x0`, `h_rate`, `v_rate` |

**Reading the results.**

- The **strain rate at the fault** is the velocity gradient `dv/dx` at *x* = 0,
  `|S| / (π d₁)` - the engineering shear strain rate; the tensor component
  ε̇ₓᵧ is half of it. For *S* = 26 mm/yr and *d₁* = 16 km it is 0.52 µstrain/yr,
  the ~0.5 µstrain/yr that Hussain et al. (2018) report along the North
  Anatolian Fault.
- *S* and *d₁* trade off: a deeper locking depth with a faster slip rate fits
  a profile almost as well, the more so the shorter the profile. Profiles a few
  times *d₁* long on either side, GNSS in the far field, or a prior on *d₁*
  narrow it - the posterior window shows by how much.
- The bins set the resolution: a creep step, or the fault position with
  `--free-x0`, is resolved to about one bin.
- The run reports the acceptance fraction and the autocorrelation time of the
  chains, and suggests more `--steps` when they are shorter than 50
  autocorrelation times.

---

## Time series: `viewinsar tseries`

LiCSBAS_plot_ts.py, for LiCSBAS **and** MintPy time series, with circles.

```bash
# LiCSBAS, in TS_GEOCml*/: cum_filt.h5, and cum.h5 in red, as LiCSBAS does
viewinsar tseries
viewinsar tseries cum_filt.h5 --radius 2                       # circles of 2 km

# MintPy: velocity.h5, maskTempCoh.h5, temporalCoherence.h5 ... found beside it
viewinsar tseries geo/geo_timeseries_ERA5_demErr.h5 --radius 2 --unit cm

# no window: the time series of one circle, as .csv + figures
viewinsar tseries cum_filt.h5 --radius 2 --p-geo 121.68/-2.51 --nodisplay

# several files - frames of a track - on one map (see below)
viewinsar tseries 061D_*/TS_GEOCml10GACOS/cum.h5 --radius 15
```

Two windows, as in LiCSBAS. The **map** shows the velocity relative to the
reference area, or a noise index from the radio buttons (LiCSBAS `results/`:
`mask`, `coh_avg`, `n_unw`, `vstd` ...; MintPy: mask, temporal coherence,
height), or, with the slider, the cumulative displacement at a date. The
**time series** window has LiCSBAS's fit models (Linear, Annual+L, Quad,
Annual+Q) and its side panel.

| input | action |
|-------|--------|
| left double click | the time series of that pixel - of the circle around it with `--radius` |
| right drag | a new reference area |
| slider | the cumulative displacement at that date |
| `mask` box | the mask on / off |
| <kbd>e</kbd> | export |

**With `--radius KM`** a double click takes every pixel whose centre is within that
distance on the ground. Each one is treated as LiCSBAS treats a single pixel:
minus the reference area's mean at every date, minus its own value at the
reference date. Then, at every date, the plot shows the **median** of the
circle, with the **25th-75th percentiles** as an error bar, and the fits
follow the medians. Masked pixels and pixels without data are left out; the
side panel says how many were used. A circle on a steep gradient gets wide
error bars, which is the point: they show how uniform the area is.

| option | default | meaning |
|--------|---------|---------|
| `--radius` | one pixel | the radius of the circle, km |
| `--i2` / `--no-i2` | `cum.h5` beside `cum_filt.h5` | a second time series, in red |
| `--labels` | the file names | names of the files in the legends |
| `--ref-date` | first date (MintPy: `REF_DATE`; several files: the latest of theirs) | where every series is 0 |
| `--ref X1:X2/Y1:Y2`, `--ref-geo`, `--ref-yx`, `--ref-lalo` | the file's | the reference area |
| `--point X/Y`, `--p-geo LON/LAT`, `--yx`, `--lalo` | the reference area's centre | the first point shown |
| `--results`, `--mask`, `--nomask`, `--vel`, `--lookup` | found | the noise indices, mask, velocity, radar lat/lon |
| `--colormap`, `--unit`, `--vlim`, `--dlim`, `--ylim` | `roma_r`, `mm`, auto | colours, unit, ranges |
| `--save` / `--nodisplay` / `--svg` / `--outfile` / `--outdir` | — | output control |

For LiCSBAS users: `--radius` here is the circle's **radius**. LiCSBAS's `-r` (the
reference area) is `--ref X1:X2/Y1:Y2`; LiCSBAS's `-p X/Y` is `--point X/Y`, and
`--p_geo` and `--ref_geo` are spelled as in LiCSBAS.

A MintPy time series without a `velocity.h5` beside it gets its velocity from
a linear fit of the time series. Radar-coded MintPy data are flipped as MintPy
shows them, and `geometryRadar.h5` gives the lat/lon.

The export (<kbd>e</kbd>, `--save` or `--nodisplay`), named after the point,
e.g. `tseries_Y156X144_r2km`:

| file | contents |
|------|----------|
| `….csv` | `date`, `year`, then per file `median_mm` `q25_mm` `q75_mm` `n_px` (with `--radius`) or `disp_mm` (one pixel) |
| `….png` / `…_map.png` | the time series / the map |
| `…_info.txt` | the file(s), the point, the circle, the reference area and date, the velocity |

### Several files: frames, tracks, versions of one series

```bash
# adjacent LiCSBAS frames of a track, on one map
viewinsar tseries 061D_*/TS_GEOCml10GACOS/cum.h5 --radius 15

# ... with names, and a reference area in their overlap from the start
viewinsar tseries frameA/cum.h5 frameB/cum.h5 --labels A B --radius 15 \
    --ref-geo 121.60/121.70/-2.60/-2.50
```

Several files share **one map**: each is drawn where it lies, with its outline,
its reference area (dashed) and its time series in a colour of its own. A
double click shows the series of **every file that covers the point** - in the
overlap of two frames, both - and the side panel says which files are not
there. The files must be geocoded in one CRS (LiCSBAS GEOC, MintPy `geo_*`),
or radar-coded on one grid (the MintPy time series of one stack, e.g.
`timeseries.h5` and `timeseries_ERA5.h5`, as `tsview.py` takes them).

- **Reference area**: each file starts with its own, so two frames' series at
  one point differ by the motion of their reference areas. Right-drag one in
  their overlap (or give `--ref-geo` / `--ref-lalo`) and every file with data
  there takes it - their series become comparable; a file without data in it
  keeps its own, and says so.
- **Reference date**: `--ref-date`, else the latest of the files' first dates
  (MintPy: `REF_DATE`), so that all start together; a file without that date
  takes its nearest one.
- **Dates**: each file keeps its own. The slider shows each file at its date
  nearest the slider's, and nothing before its first date or after its last.
- `--ref X1:X2/Y1:Y2`, `--ref-yx`, `--yx`, `--point` and `--mask` are in
  pixels: they work for files on one grid only. `--i2`, `--results` and
  `--vel-file` are for one file - several files each bring their own.
- `--labels` names the files in the legends; by default they are told apart by
  their names, or by the folders that differ (`061D_09234_131313`).

The export is named after the point's lat/lon (`tseries_lat-2.5100_lon121.6800_r15km`;
`_N…_E…` on a projected grid). The `.csv` has every date of the files, with a set
of columns per file (`median_mm_1` … `n_px_2`, `nan` where a file has no such
date), and `_info.txt` a block per file: its pixel, circle, reference area (in
pixels and lon/lat), reference date and velocity.

---

## Your own defaults: `viewinsar.ini`

`~/.config/viewinsar/viewinsar.ini` (or `$XDG_CONFIG_HOME/viewinsar/`, or the
file `VIEWINSAR_CONFIG` names) holds rules that apply to every run. Its one
section so far, `[data-scale]`, multiplies the data of a mission - or of files
whose path matches a pattern - before anything else, in `swath`, `stat`,
`model` and `tseries`:

```ini
[data-scale]
# every ALOS-2 product of mine is 10x too large
ALOS-2 = 0.1
# ... or by the path, with * and ?
*/T125D_alos2/*.h5 = 0.1
```

The mission is told by the file's `PLATFORM` attribute (MintPy writes `Alos2`,
`Sen` …); a file without one is told by its name and its folders' (`alos2`,
`scansar`, `s1`, LiCSAR's `105D_…`), and not at all when they name two
missions. Every run says when a rule applies, and `--data-scale` on the command
line comes first (`--data-scale 1` switches the rule off for that run).

When a file's `WAVELENGTH` is not its mission's (ALOS-2 0.238 m, Sentinel-1
0.0555 m) viewinsar says so: a displacement computed with a wrong wavelength is
wrong by the same ratio - a wavelength of 2.4 m for ALOS-2 makes every velocity
10x too large. `--data-scale` corrects what you see; re-running the time
series with the right wavelength corrects the file.

## Python API

The numerical core has no MintPy or matplotlib dependency, so it is usable on
its own:

```python
import numpy as np
from viewinsar import swath_extract, swath_density

swath = swath_extract(data, atr, box=(0, 0, width, length),
                      start_yx=[300, 150], end_yx=[300, 750], half_width=15.)
den = swath_density(swath['dist'], swath['value'], bin_x=0.5, bin_y=1.5,
                    ylim=[-40, 40])

print(swath['dist'], swath['perp'], swath['value'])
print(den['x_centers'], den['median'], den['q25'], den['q75'])
```

The fault helpers work the same way:

```python
from viewinsar import load_faults, profile_crossings

faults = load_faults(['faults.shp'], atr)          # traces in the data's CRS
dist_km, xy = profile_crossings(faults, atr, [300, 150], [300, 750], swath['length'])
```

The creep along a line too, and the two lines of Hussain et al.'s eq. 2:

```python
from viewinsar import creep_along_line, fit_two_lines

res = creep_along_line(data, atr, box=(0, 0, width, length),
                       vertices_yx=[[3670, 744], [3300, 720], [2967, 696]],
                       spacing_km=5, length_km=5, width_km=2.5, factor=10)  # cm/yr -> mm/yr
res['dist_km']                          # the points along the line
prof = res['profiles'][k]               # point k: its pixels, x_km, y_km, values ...
prof['fit']['mean'], prof['fit']['std'] # the creep rate (block 2 - block 1), bootstrap
prof['fit']['offset']                   # ... and least squares (the BLUE)
prof['fit']['boot']                     # every draw: slope_1, slope_2, b1, creep

fit = fit_two_lines(x_1, v_1, x_2, v_2, n_boot=2000)    # any two sides of a fault
```

And the time series, from LiCSBAS or MintPy alike:

```python
from viewinsar.tsdata import TimeSeries, circle_pixels_rc, relative_series, series_stats

ts = TimeSeries('cum_filt.h5')                   # or a MintPy timeseries.h5, in mm
rows, cols = circle_pixels_rc(ts, 156, 144, radius_km=2)
ref = ts.area_mean(ts.ref_area, ts.keep)         # the reference area, every date
rel = relative_series(ts.pixels(rows, cols), ref, ix_ref=0)
stats = series_stats(rel)                        # median, q25, q75, count per date
```

The dislocation model and its sampler, and the GNSS readers:

```python
import numpy as np
from viewinsar import fit_profile, read_gnss
from viewinsar.model import forward, model_band

x = np.random.uniform(-80, 80, 50000)                       # km from the fault
v = forward(x, 20.0, 12.0, 3.0) + np.random.normal(0, 2, x.size)   # mm/yr
res = fit_profile((x, v), bin_km=1.0, walkers=64, steps=5000)      # creep=True ...
res['summary']['S']          # {'map', 'median', 'lo', 'hi'}: the MAP and the 95 % range
res['samples']               # the posterior, columns res['names'] ('S', 'd1', 'a')
best, lo, hi = model_band(res, np.linspace(-80, 80, 401))

gnss = read_gnss('gnss.xlsx')                               # any table or point file
value, sigma = gnss.component('ve')

# a reverse fault dipping 30 degrees towards x > 0: the motion across it and up
from viewinsar.model import dip_slip
h = dip_slip(x, 8.0, 12.0, 30.0, kind='h') + np.random.normal(0, 1, x.size)
z = dip_slip(x, 8.0, 12.0, 30.0, kind='z') + np.random.normal(0, 1, x.size)
res = fit_profile([(x, h), (x, z)], fault='reverse', kinds=['h', 'z'], dip=(5, 85))
res['summary']['dip'], res['derived']['h_rate'], res['derived']['v_rate']
```

`viewinsar.swathProfileViewer`, `viewinsar.StatViewer`,
`viewinsar.TimeSeriesViewer` and `viewinsar.ModelViewer` are the interactive
figures; they are imported lazily, so `import viewinsar` stays cheap.

---

## Notes and caveats

- **Flat caps.** The corridor is a rectangle, not a QGIS round buffer: pixels
  beyond either end point are excluded rather than clamped, so nothing piles up
  at the ends of the profile.
- **Coordinates.** Geographic grids are converted to metres with a local
  equirectangular approximation about the profile's mid-latitude — the same one
  MintPy uses in `transect_yx()`. Fine over tens of km; for continental-scale
  profiles work in UTM.
- **Reading window.** Data is read within the displayed window, so
  `--sub-lat/--sub-lon` bounds the memory. If the corridor pokes outside that
  window, the tool says so instead of silently biasing the density.
- **Zeros.** Exactly-zero pixels are dropped by default, matching
  `plot_transection.py`. Use `--keep-zero` if zero is a real value in your data.
  The DEM behind `--topoxsect` always keeps zeros — 0 m is sea level.
- **`--nlook` interpolates, it does not invent information.** It restores the
  sampling density a multilooked raster lost, which makes the density panel and
  the median readable again, but scatter that multilooking averaged away is
  gone for good. The median it produces is unbiased with respect to the
  full-resolution one.
- **The model is two-dimensional.** An infinitely long, straight, vertical
  strike-slip fault in an elastic half-space: the profile should cross a
  straight stretch of the fault, away from its ends and bends. Vertical
  motion, fault-normal motion and other faults in the swath all end up in the
  residuals - or in the parameters.
- **Fault crossings follow the centre line.** A trace that ends inside the
  corridor without reaching the centre line, or runs along it, is drawn on the
  map but not marked on the profile. Traces are re-projected vertex by vertex,
  so a long straight segment drawn in one CRS bends very slightly in another;
  digitised fault traces are dense enough for this not to matter.

## License

GPL-3.0-or-later, following MintPy.

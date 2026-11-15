# viewinsar

**Interactive swath (buffered) profile viewer for InSAR velocity and time-series
products, built on [MintPy](https://github.com/insarlab/MintPy).**

`plot_transection.py` samples a single line of pixels. `viewinsar` samples a
**corridor**: you draw the same line with two clicks, it gets buffered by
±*W* km, every valid pixel inside the buffer is linearly referenced onto the
line, and the resulting cloud is shown as a 2-D point-density profile with the
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

`viewinsar stat` compares the two sides of a fault: draw a line along it, and
every *SP* km a pair of circles is placed, one on each side, **1 km clear of the
line**; the values in each circle become a violin, block 1 against block 2,
with the block medians joined by dashed lines (see
[Block statistics](#block-statistics-viewinsar-stat)).

`viewinsar tseries` is LiCSBAS_plot_ts.py for LiCSBAS and MintPy time series,
with circles: `--radius 2` and a double click gives, at every date, the median of
all pixels within 2 km, with the 25th-75th percentiles as error bars (see
[Time series](#time-series-viewinsar-tseries)).

Several files at once - the frames of a track, or several tracks - make one
density with a median line for each, and are drawn together on the map; GNSS
velocities (`--gnss`, any table or point format) appear on the map and, inside
a swath, on its profile with their error bars (see
[Several files and GNSS](#several-files-and-gnss)). `viewinsar stat -rate`
turns the two sides of a fault into the rate across it.

`viewinsar model` fits the **slip rate and locking depth** of an elastic
dislocation - and shallow creep - to a profile across a fault, by ensemble
MCMC, as Hussain et al. (2016, 2018) did for the North Anatolian Fault: the
pixels of the swath as points, the maximum a posteriori model with its 95 %
range, and the posterior distributions (see
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

# no DEM file? download one for the map: shaded relief + elevation profile
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

# block statistics on both sides of a fault: circles every 5 km, 2 km radius
viewinsar stat geo_velocity.h5 --spacing 5 --radius 2 --fault faults.shp
# ... or the rate across it (block 1 - block 2), of the east and up velocities
viewinsar stat geo_ve.h5 geo_vu.h5 --rate --spacing 5 --radius 2 --line-file fault_trace.geojson

# slip rate and locking depth of a fault, by MCMC, from a profile across it
viewinsar model geo_velocity.h5 --comp los --fault faults.shp --width 25

# time series, LiCSBAS-style: the median of a 2 km circle at every date
viewinsar tseries cum_filt.h5 --radius 2
viewinsar tseries geo_timeseries_ERA5_demErr.h5 --radius 2

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
| `--dem` | — | DEM for the shaded-relief map background, as in view.py; **without a file, downloaded** for the map extent (see below) |
| `--topoxsect` | — | DEM to profile under the LOS panel, any format / any CRS; **without a file**, the `--dem` DEM or a download |
| `--dem-zoom` | auto | zoom of the downloaded tiles, 0–15 (12 ≈ 38 m, 13 ≈ 19 m) |
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
| `stat`: `-sp`, `-r`, `-rate` | `--spacing`, `--radius`, `--rate` |
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
| `--dem --topoxsect` | downloaded | the same download |
| `--dem dem.tif --topoxsect` | `dem.tif` | `dem.tif` |

The download comes from the [Terrain Tiles](https://registry.opendata.aws/terrain-tiles/)
on AWS Open Data, the "Mapzen Global Terrain" layer of QGIS:
`https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png`.
These are 256-px PNG tiles in Web Mercator, with the height coded in the colour
as R·256 + G + B/256 − 32768 m. The tiles over the map are mosaicked,
re-gridded onto the data's CRS (lat/lon or UTM) and saved as a MintPy geometry
`.h5` (`height`, metres) in `~/.cache/viewinsar` (or `$VIEWINSAR_CACHE`). A
second run over the same area downloads nothing.

The zoom is the finest that keeps the download within about 300 tiles (at
most 12, about 38 m). A big map gets a coarser DEM, which is still plenty for
a median line in 0.5 km bins. `--dem-zoom 13` with `--sub-lat/--sub-lon` gives
about 19 m over a small area. The sources are SRTM, GMTED2010 and ETOPO1 (with
bathymetry) worldwide, plus finer national DEMs where available. Please cite
them: the attribution is printed and written to `_info.txt`.

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

Only a geometry file without `height` makes `--dem` download the Terrain Tiles,
resampled onto the radar grid.

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

Every station is drawn on the map - a **triangle** for a 3D station (one with
`vu`), a **circle** for a 2D one, in black. The stations inside a swath (the
same flat-capped buffer as the pixels) are drawn on its profile, in black with
their 1-σ as error bars, and a **smoothed black line** goes through them (a
smoothing spline weighted by 1/σ²), as a comparison with the InSAR. By default
they are shifted onto the InSAR (`--gnss-ref insar`): by the median, over the
stations, of the InSAR within `--gnss-radius` km of a station minus the
station; `--gnss-ref none` leaves them as read. `--gnss-label` writes the site
names on the map, `--gnss-nofit` leaves the line out.

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

## Block statistics: `viewinsar stat`

How different are the two sides of a fault? Draw a line along it, and

1. every `--spacing` km along the line a **point** is placed (red on the map),
   centred on the line so the first and last sit equally far from its ends;
2. at every point, two **circles** of radius `--radius` km: one on each side, its
   edge `--gap` km (default **1 km**) clear of the line;
3. **block 1** is everything on the **left** of the line, walking from its
   first point to its last, **block 2** everything on the right;
4. the pixels inside each circle become a **violin**, block 1 and block 2 side
   by side at each point, each block in **the colour of its circles on the
   map** (the map and the violins share one legend style, and the map names
   the blocks next to their circles);
5. each violin has **strikes** at its **median** (in the block's colour), its
   **mean** (solid) and its **mode** (dotted: the peak of the kernel density,
   where the violin is widest), and the medians of each block are joined by a
   **dashed line** along the line.

```bash
# draw the line on the map: left-click along the fault, Enter to finish
viewinsar stat geo_velocity.h5 --spacing 5 --radius 2 --fault faults.shp --colormap RdBu_r --vlim -1 1

# or give it: lat/lon or y/x pairs, or the longest line of a vector file
viewinsar stat geo_velocity.h5 --spacing 5 --radius 2 --save --nodisplay \
    --line-lalo -1.20 119.95 -1.55 120.05 -1.80 120.20
viewinsar stat geo_velocity.h5 --spacing 5 --radius 2 --line-file fault_trace.geojson --svg --save

# your own names and colours
viewinsar stat geo_velocity.h5 --spacing 5 --radius 2 --names "West block" "East block"
```

**Every circle stays on its own side.** A circle centre starts `r + gap` out
from its point, along the line's normal. On the inside of a bend that can
bring it within `r + gap` of *another* part of the line, so the circle is
pushed further out along the same normal until it clears **every** part of
the line by the gap; where there is no such spot (a hairpin), it is left out,
and the run says so. A circle is also never placed on the far side of
another strand of the line. So every pixel used is at least the gap away from
the line, and on its own side of it - `_values.csv` lists that distance for
each one. The normal follows the line's direction over ±*r* around the point,
so a hand-digitised trace with small wiggles still gives steady circles.

"Left" and "right" are those of the ground, as the map shows it, also on
radar-coded data, whose grid is a mirror image of the ground (that is why
view.py flips it): there, `geometryRadar.h5` gives the lat/lon of the line
(`--line-lalo`, `--line-file`) and of every pixel in the `.csv` files.

| input | action |
|-------|--------|
| left click | add a point to the line |
| <kbd>Enter</kbd> | finish the line: points, circles on the map, violins |
| <kbd>u</kbd> or <kbd>Backspace</kbd> | remove the last point |
| <kbd>n</kbd> | start a new line (so does a click after <kbd>Enter</kbd>) |
| <kbd>e</kbd> | export |

| option | default | meaning |
|--------|---------|---------|
| `--spacing` | `5` | distance between the points along the line, km |
| `--radius` | `2` | radius of the circles, km |
| `--gap` | `1` | space between the line and every circle, km |
| `--line-yx` / `--line-lalo` | — | the line as row/col or lat/lon pairs, 2 points or more |
| `--line-file` | — | the longest line of a vector file (any format / CRS `--fault` reads) |
| `--names` / `--colors` | `Block 1` `Block 2` / blue, vermilion | the two blocks |
| `--ylim` | auto | value range of the violin plot |
| `--keep-zero` | off | keep exactly-zero pixels |
| `--dem`, `--dem-zoom`, `--lookup`, `--fault` | — | as in `viewinsar swath` |
| `--save` / `--nodisplay` / `--svg` / `--outfile` / `--outdir` | — | output control |

Everything else goes to view.py, as for `swath`. The export (<kbd>e</kbd> or
`--save`) is named after the line's first and last point and the parameters,
e.g. `stat_Y3670X744_Y2967X696_sp5km_r2km`:

| file | contents |
|------|----------|
| `….png` | the violins |
| `…_map.png` | the map with the line, the points and the circles |
| `…_stats.csv` | one row per point: `dist_km`, its row/col and lat/lon, `count` `median` `mean` `mode` `std` `q25` `q75` `min` `max` of each block (`_b1`, `_b2`), the circle offsets `offset_km_b1/_b2`, and `median_b1_minus_b2` |
| `…_values.csv` | every pixel used: `point`, `block`, `value`, `row`, `col`, `dist_to_line_km`, lat/lon |
| `…_info.txt` | the line's vertices, the parameters, which side is which block |

When `--spacing` is smaller than `2r` neighbouring circles overlap and share pixels;
the run notes it.

### Several files, and the rate across the fault

Several files - tracks, frames, or the east and up velocities of a
decomposition - get the same circles, and a row of violins each.

With **`--rate`** the violins give way to the **rate** across the fault: at
every point, the median of block 1 minus the median of block 2, one series per
file, with error bars - the slip (or creep) rate along the fault, as Hussain
et al. (2016) measured the creep along the North Anatolian Fault. The east and
up velocities together give the horizontal and the vertical rate:

```bash
viewinsar stat geo_ve.h5 geo_vu.h5 --labels East Up --rate --spacing 5 --radius 2 \
    --line-file fault_trace.geojson
```

| `--rate-err` | the error bar |
|---|---|
| `iqr` (default) | `sqrt(s₁² + s₂²)`, with `s` = IQR / 1.349 of each circle: the robust spread of the pixels - conservative, since the pixels of a circle are correlated |
| `sem` | the standard error of the two medians, `1.2533 sqrt(s₁²/n₁ + s₂²/n₂)`, as if the pixels were independent |

`…_rate.csv` has one row per point: `dist_km`, its position, and `rate_k`,
`rate_err_k` of every file. A second file's statistics go to `…_2_stats.csv`
and `…_2_values.csv`.

---

## Fault slip models: `viewinsar model`

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

**The fault and the data.** The fault is where the profile crosses a `--fault`
trace (the crossing nearest the middle of the profile), with the strike of the
trace there - or `--x0` (km from the start of the profile) and `--strike`.
Distances are taken perpendicular to the fault, so a profile drawn at an angle
to it still gives the right locking depth. `--comp` says what the data are:

| `--comp` | the velocity fitted |
|---|---|
| `auto` (default) | `los` when the line of sight of every file is known (below; nominal angles do not count) and no file is named as a component of a decomposition (`hz`, `up`, `east`, `ve` …); `as-is` otherwise, with the GNSS in the LOS (`--vel los`) - the run says which |
| `los` | LOS velocity, divided by (LOS unit vector · strike direction), pixel by pixel |
| `par` | fault-parallel velocity, positive along `--strike` |
| `ve`, `vn` | east / north velocity, divided by sin / cos of the strike |
| `as-is` | the data as they are - e.g. the LOS itself; *S* is then the rate seen in the LOS ("slip rate in the LOS" when `auto` chose it for LOS data) |

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
geometry GeoTIFFs, is enough. `--los-inc` and `--los-az` give one
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
with their error bars, the MAP model as a bold dashed line and its 95 % range
in grey, the fault dotted, and the parameters written as Hussain et al. do -
`slip rate = 20.3 mm/yr (18.4 – 34.0, 95% CI), right-lateral`, `locking depth
= 12.7 km (…)`, `creep rate` in mm/yr, `creep depth` in km, and the `strain
rate` at the fault in µstrain/yr - placed where they hide the fewest points.
The slip rate is given as a rate with its sense; when its 95 % range crosses
zero it keeps its sign (S > 0 right-lateral). Every file is shifted
onto the offset of the first, the GNSS onto the InSAR. `--strain` adds the
shear strain rate of the model, `dv/dx`, below; `--topoxsect` the topography.
A second window holds the **posterior**, drawn as Hussain et al. (2016,
Fig. 9): every pair of parameters (*S*, *d₁*, the offset *a*, and *C*, *d₂*,
*x₀* when fitted) as a 2D histogram in `jet`, each scaled to its own peak, the
MAP as a red dot; along the bottom the histogram of each parameter in blue,
the MAP as a red line, the 95 % range dashed and written above it; the
section, `A - A'`, white on black in the corner. The trade-off between *S* and
*d₁* shows there.

| option | default | meaning |
|---|---|---|
| `--width` | `20` | half width of the swath, km |
| `--bin` | `1` | bin size along the profile, km |
| `--comp` | `auto` | what the data are (above) |
| `--geom`, `--los-inc`, `--los-az` | found | the line of sight (above), one per file |
| `--data-unit` | the file's own; LiCSBAS: mm/yr | the unit the files are really in, when they say none or a wrong one |
| `--x0`, `--strike` | from the trace | the fault position (km along the profile) and its strike |
| `--creep` | off | fit *C* and *d₂* too |
| `--free-x0 [KM]` | off | fit the fault position within ±KM (default 5) |
| `--rot RATE`, `--rot-side` | none, `end` | a fixed linear term on one side, mm/yr per km |
| `--s-max`, `--d1-max`, `--c-max`, `--d2-max` | 60, 60, 30, 40 | the uniform priors: \|S\| ≤, 0.1 ≤ d₁ ≤, 0 ≤ C ≤, d₂ ≤ |
| `--d1-prior MEAN SIGMA` | none | a Gaussian prior on *d₁* (Hussain et al. 2016, Izmit: 17 5) |
| `--noise` | `auto` | `auto`: exponential, fitted to the residuals of a least-squares fit, for each file; `white`: each bin's own scatter; `SILL LENGTH`: given, in (mm/yr)² and km |
| `--walkers`, `--steps`, `--burn`, `--seed` | 64, 5000, a third, 0 | the sampler |
| `--gnss`, `--vel`, `--gnss-width`, `--gnss-nofit` | - | GNSS |
| `--strain`, `--ylim`, `--max-points` | off, auto, 200000 | the figure |

The export (<kbd>e</kbd> or `--save [DIR]`, into `DIR` - made if it does not
exist), named after the line, e.g. `model_A_Y470X150_Y130X700_w25km`:

| file | contents |
|---|---|
| `….png` | the fit |
| `…_posterior.png` | the posterior distributions |
| `…_params.csv` | every parameter: `map`, `median`, `lower_95`, `upper_95`; and the strain rate at the fault |
| `…_params.txt` | the same, readable, with the sense of slip, the fault position and strike, the noise, the priors, and the sampler and its diagnostics |
| `…_samples.csv` | posterior samples (up to 20 000), for plots of your own |
| `…_bins.csv` | the bins: `x_km`, `median`, `std`, quartiles, `count`, the MAP model at each and the residual |
| `…_points.csv` | every pixel: `dist_km` along the profile, `x_km` from the fault, `perp_km`, `value` as read, `v_model` (the velocity fitted), `v_shifted` (as drawn), row / col, lon / lat, `file_idx` |
| `…_gnss.csv` | the GNSS in the swath, and whether they were fitted |
| `model_map_AB.png` | the map with every section, the fault position of each marked ◇ |
| `model_sections_AB.csv` | two sections or more: one row each - its ends (`start_lat` … `end_lon`), where it crosses the fault (`fault_km`, `fault_lat`, `fault_lon`), `strike_deg`, `sense`, `pixels`, and `S`, `d1` (`C`, `d2`, `x0`) and `strain` as `_map`, `_lower_95`, `_upper_95` |

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
| `--ref-date` | first date (MintPy: `REF_DATE`) | where every series is 0 |
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

The block statistics too:

```python
from viewinsar import stat_line

res = stat_line(data, atr, box=(0, 0, width, length),
                vertices_yx=[[3670, 744], [3300, 720], [2967, 696]],
                spacing_km=5, radius_km=2, gap_km=1)
res['dist_km']                  # the points along the line
res['stats'][0][k]['median']    # block 1 (left), point k; [1] = block 2
res['stats'][0][k]['mode']      # the peak of its kernel density (also 'mean', 'q25' ...)
res['values'][1][k]             # the pixel values in block 2's circle at point k
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

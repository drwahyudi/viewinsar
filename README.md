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
(`-fault faults.shp`, any vector format, any CRS) and they are drawn on the map
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
with circles: `-r 2` and a double click gives, at every date, the median of
all pixels within 2 km, with the 25th-75th percentiles as error bars (see
[Time series](#time-series-viewinsar-tseries)).

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

MintPy itself is easiest to get from conda-forge; if you already have a MintPy
environment, install `viewinsar` into it.

**Requirements:** Python ≥ 3.9, `mintpy` ≥ 1.5, `numpy` ≥ 1.20, `matplotlib` ≥ 3.5,
`pyproj`, `pyshp`, `shapely`, `h5py`, `pillow` and `scipy` (all already come with MintPy
and matplotlib).
Optional: `pip install "viewinsar[vector]"` adds `pyogrio`, so `--fault` reads
every vector format GDAL knows, not only the common ones it reads by itself.

---

## Quick start

```bash
# interactive - click 2 points on the map to draw the swath centre line
viewinsar swath geo_velocity_msk.h5 -w 15 --bin-x 0.5 --bin-y 1.5 --ylim -40 40

# non-interactive - give the line up front, write the figure + the .csv files
viewinsar swath geo_velocity.h5 \
    --start-lalo 30.125 129.988 --end-lalo 30.250 130.116 \
    -w 15 --save --nodisplay

# a coarse / multilooked raster: oversample the corridor so the profile fills in
viewinsar swath geo_velocity_ml10.h5 --nlook 10

# an elevation profile under the LOS panel, from a DEM in any format or CRS
viewinsar swath geo_velocity.h5 --topoxsect srtm1.dem

# no DEM file? download one for the map: shaded relief + elevation profile
viewinsar swath geo_velocity.h5 -d --topoxsect

# fault traces on the map, and a dashed "Fault line" where a section crosses one
viewinsar swath geo_velocity.h5 -fault faults.shp
viewinsar swath geo_velocity.h5 -fault gem_active_faults.geojson -topoxsect dem.tif

# the point time series, in its own window (right-click a pixel on the map)
viewinsar swath geo_velocity.h5 --timeseries geo_timeseries_ERA5_demErr.h5

# vector output
viewinsar swath geo_velocity.h5 --topoxsect dem.tif --svg --save

# block statistics on both sides of a fault: circles every 5 km, 2 km radius
viewinsar stat geo_velocity.h5 -sp 5 -r 2 -fault faults.shp

# time series, LiCSBAS-style: the median of a 2 km circle at every date
viewinsar tseries cum_filt.h5 -r 2
viewinsar tseries geo_timeseries_ERA5_demErr.h5 -r 2

# list what is inside a file
viewinsar list timeseries.h5
```

`viewinsar FILE ...` is a shortcut for `viewinsar swath FILE ...`.

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
viewinsar swath geo_velocity.h5 -c jet -v -3 3 -u cm -d srtm1.dem \
    --mask geo_maskTempCoh.h5 --sub-lat 30.0 30.5 --sub-lon 129.9 130.3 \
    --coastline --lalo-label --fontsize 12
```

The swath-specific options:

| option | default | meaning |
|--------|---------|---------|
| `-w`, `--width`, `--buffer` | `15` | **half** width of the buffer, in km |
| `--bin-x` | `0.5` | along-profile cell size of the density grid, in km |
| `--bin-y` | auto | value cell size, in the display unit (~55 bins over the Y range) |
| `--ylim` | 0.5–99.5 % | value range of the profile plot |
| `--cmin` / `--cmax` | `1` / auto | limits of the log density colour scale; `--cmax` defaults to the 98th percentile of the occupied cells |
| `--pcmap` | `viridis` | colormap of the density |
| `--x-offset` | `0` | shift the distances so a known feature sits at *x* = 0 |
| `--keep-zero` | off | keep exactly-zero pixels (dropped by default, as in `plot_transection.py`) |
| `--decimate` | `1` | use every *N*th pixel of the corridor |
| `--nlook` | `1` | **oversample** the corridor *N*× per axis (*N*² more samples). The cure for a sparse profile from multilooked data; inverse of `--decimate` |
| `-d`, `--dem` | — | DEM for the shaded-relief map background, as in view.py; **without a file, downloaded** for the map extent (see below) |
| `--topoxsect` | — | DEM to profile under the LOS panel, any format / any CRS; **without a file**, the `-d` DEM or a download |
| `--dem-zoom` | auto | zoom of the downloaded tiles, 0–15 (12 ≈ 38 m, 13 ≈ 19 m) |
| `--lookup` | auto | radar-coded data only: MintPy's `geometryRadar.h5` (found in `inputs/`, `./` or `../inputs/` when not given) |
| `--fault` | — | fault trace file(s), any vector format / any CRS: drawn on the map, and a dashed vertical line + "Fault line" on every profile at each crossing |
| `--fault-label` / `--fault-color` / `--fault-lw` | `Fault line` / `k` / `1.2` | the annotation text (`""` for none), the colour, and the map line width |
| `--svg` | off | write figures as vector `.svg` |
| `--figsize-map` | auto | size of the map window (`--figsize` sizes each section) |
| `--timeseries` | — | time-series file(s) for the 3rd panel |
| `--start-yx` / `--end-yx` | — | line end points in pixel numbers |
| `--start-lalo` / `--end-lalo` | — | line end points in lat/lon |
| `--line-file` | — | GMT-format `lon lat` file with the two end points |
| `--save` / `--nodisplay` / `-o` / `--outdir` | — | output control |

`-nlook`, `-topoxsect`, `-fault` and `-svg` work with one dash or two. Run
`viewinsar swath -h` for the full list.

### Topography: your DEM, or one downloaded for you

`--topoxsect dem.tif` profiles a DEM in any format and any CRS. Its CRS is read
from the file itself, the section is re-projected onto it, and it is sampled
at its own resolution, independently of `--nlook`. A DEM that covers part of
a section gives that part, and a section that gets no topography says why in
its panel.

Without a DEM file, `-d` and `--topoxsect` fetch one:

| command | map background | elevation profile |
|---------|----------------|-------------------|
| `-d` | downloaded | — |
| `--topoxsect` | — | downloaded |
| `-d --topoxsect` | downloaded | the same download |
| `-d dem.tif --topoxsect` | `dem.tif` | `dem.tif` |

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
| `-d`, `--topoxsect` (no file) | its `height`, the DEM of your own processing, already on the radar grid; nothing is downloaded |
| `--topoxsect dem.tif` | a geocoded DEM, reached through its latitude/longitude |
| `--fault faults.shp` | traces mapped from lon/lat onto the radar pixels; `_info.txt` gives each crossing's row/col and lat/lon |

Only a geometry file without `height` makes `-d` download the Terrain Tiles,
resampled onto the radar grid.

### Fault traces

`-fault` takes one or more vector files and reads `.shp` (with its `.prj`),
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

1. every `-sp` km along the line a **point** is placed (red on the map),
   centred on the line so the first and last sit equally far from its ends;
2. at every point, two **circles** of radius `-r` km: one on each side, its
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
viewinsar stat geo_velocity.h5 -sp 5 -r 2 -fault faults.shp -c RdBu_r -v -1 1

# or give it: lat/lon or y/x pairs, or the longest line of a vector file
viewinsar stat geo_velocity.h5 -sp 5 -r 2 --save --nodisplay \
    --line-lalo -1.20 119.95 -1.55 120.05 -1.80 120.20
viewinsar stat geo_velocity.h5 -sp 5 -r 2 --line-file fault_trace.geojson --svg --save

# your own names and colours
viewinsar stat geo_velocity.h5 -sp 5 -r 2 --names "West block" "East block"
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
| `-sp`, `--spacing` | `5` | distance between the points along the line, km |
| `-r`, `--radius` | `2` | radius of the circles, km |
| `--gap` | `1` | space between the line and every circle, km |
| `--line-yx` / `--line-lalo` | — | the line as row/col or lat/lon pairs, 2 points or more |
| `--line-file` | — | the longest line of a vector file (any format / CRS `--fault` reads) |
| `--names` / `--colors` | `Block 1` `Block 2` / blue, vermilion | the two blocks |
| `--ylim` | auto | value range of the violin plot |
| `--keep-zero` | off | keep exactly-zero pixels |
| `-d`, `--dem-zoom`, `--lookup`, `--fault` | — | as in `viewinsar swath` |
| `--save` / `--nodisplay` / `--svg` / `-o` / `--outdir` | — | output control |

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

When `-sp` is smaller than `2r` neighbouring circles overlap and share pixels;
the run notes it.

---

## Time series: `viewinsar tseries`

LiCSBAS_plot_ts.py, for LiCSBAS **and** MintPy time series, with circles.

```bash
# LiCSBAS, in TS_GEOCml*/: cum_filt.h5, and cum.h5 in red, as LiCSBAS does
viewinsar tseries
viewinsar tseries cum_filt.h5 -r 2                       # circles of 2 km

# MintPy: velocity.h5, maskTempCoh.h5, temporalCoherence.h5 ... found beside it
viewinsar tseries geo/geo_timeseries_ERA5_demErr.h5 -r 2 -u cm

# no window: the time series of one circle, as .csv + figures
viewinsar tseries cum_filt.h5 -r 2 --p-geo 121.68/-2.51 --nodisplay
```

Two windows, as in LiCSBAS. The **map** shows the velocity relative to the
reference area, or a noise index from the radio buttons (LiCSBAS `results/`:
`mask`, `coh_avg`, `n_unw`, `vstd` ...; MintPy: mask, temporal coherence,
height), or, with the slider, the cumulative displacement at a date. The
**time series** window has LiCSBAS's fit models (Linear, Annual+L, Quad,
Annual+Q) and its side panel.

| input | action |
|-------|--------|
| left double click | the time series of that pixel - of the circle around it with `-r` |
| right drag | a new reference area |
| slider | the cumulative displacement at that date |
| `mask` box | the mask on / off |
| <kbd>e</kbd> | export |

**With `-r KM`** a double click takes every pixel whose centre is within that
distance on the ground. Each one is treated as LiCSBAS treats a single pixel:
minus the reference area's mean at every date, minus its own value at the
reference date. Then, at every date, the plot shows the **median** of the
circle, with the **25th-75th percentiles** as an error bar, and the fits
follow the medians. Masked pixels and pixels without data are left out; the
side panel says how many were used. A circle on a steep gradient gets wide
error bars, which is the point: they show how uniform the area is.

| option | default | meaning |
|--------|---------|---------|
| `-r`, `--radius` | one pixel | the radius of the circle, km |
| `--i2` / `--no-i2` | `cum.h5` beside `cum_filt.h5` | a second time series, in red |
| `-m`, `--ref-date` | first date (MintPy: `REF_DATE`) | where every series is 0 |
| `--ref X1:X2/Y1:Y2`, `--ref-geo`, `--ref-yx`, `--ref-lalo` | the file's | the reference area |
| `-p X/Y`, `--p-geo LON/LAT`, `--yx`, `--lalo` | the reference area's centre | the first point shown |
| `--results`, `--mask`, `--nomask`, `--vel`, `--lookup` | found | the noise indices, mask, velocity, radar lat/lon |
| `-c`, `-u`, `--vlim`, `--dlim`, `--ylim` | `roma_r`, `mm`, auto | colours, unit, ranges |
| `--save` / `--nodisplay` / `--svg` / `-o` / `--outdir` | — | output control |

For LiCSBAS users: `-r` here is the **radius**. LiCSBAS's `-r` (the reference
area) is `--ref X1:X2/Y1:Y2`. `-p X/Y`, `--p_geo` and `--ref_geo` are spelled
as in LiCSBAS.

A MintPy time series without a `velocity.h5` beside it gets its velocity from
a linear fit of the time series. Radar-coded MintPy data are flipped as MintPy
shows them, and `geometryRadar.h5` gives the lat/lon.

The export (<kbd>e</kbd>, `--save` or `--nodisplay`), named after the point,
e.g. `tseries_Y156X144_r2km`:

| file | contents |
|------|----------|
| `….csv` | `date`, `year`, then per file `median_mm` `q25_mm` `q75_mm` `n_px` (with `-r`) or `disp_mm` (one pixel) |
| `….png` / `…_map.png` | the time series / the map |
| `…_info.txt` | the file(s), the point, the circle, the reference area and date, the velocity |

---

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

`viewinsar.swathProfileViewer`, `viewinsar.StatViewer` and
`viewinsar.TimeSeriesViewer` are the interactive figures; they are imported
lazily, so `import viewinsar` stays cheap.

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
- **Fault crossings follow the centre line.** A trace that ends inside the
  corridor without reaching the centre line, or runs along it, is drawn on the
  map but not marked on the profile. Traces are re-projected vertex by vertex,
  so a long straight segment drawn in one CRS bends very slightly in another;
  digitised fault traces are dense enough for this not to matter.

## License

GPL-3.0-or-later, following MintPy.

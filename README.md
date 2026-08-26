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

The three panels come from three MintPy tools: the map is `view.py`, the
click-to-draw interaction is `plot_transection.py`, and the optional
displacement time series is `tsview.py`.

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

**Requirements:** Python ≥ 3.9, `mintpy` ≥ 1.5, `numpy` ≥ 1.20, `matplotlib` ≥ 3.5.

---

## Quick start

```bash
# interactive - click 2 points on the map to draw the swath centre line
viewinsar swath geo_velocity_msk.h5 -w 15 --bin-x 0.5 --bin-y 1.5 --ylim -40 40

# non-interactive - give the line up front, write the figure + the .csv files
viewinsar swath geo_velocity.h5 \
    --start-lalo 30.125 129.988 --end-lalo 30.250 130.116 \
    -w 15 --save --nodisplay

# add a 3rd panel with the point time series (right-click a pixel on the map)
viewinsar swath geo_velocity.h5 --timeseries geo_timeseries_ERA5_demErr.h5

# list what is inside a file
viewinsar list timeseries.h5
```

`viewinsar FILE ...` is a shortcut for `viewinsar swath FILE ...`.

### Mouse and keyboard

| input | action |
|-------|--------|
| left click ×2 | set the start / end point of the swath centre line |
| right click | plot the time series of that pixel (needs `--timeseries`) |
| <kbd>e</kbd> | export the swath points + the binned stats (`.csv`) and the figure |
| <kbd>n</kbd> | clear the current line and start a new one |

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
| `--cmin` / `--cmax` | `1` / `100` | limits of the log density colour scale |
| `--pcmap` | `viridis` | colormap of the density |
| `--x-offset` | `0` | shift the distances so a known feature sits at *x* = 0 |
| `--keep-zero` | off | keep exactly-zero pixels (dropped by default, as in `plot_transection.py`) |
| `--decimate` | `1` | use every *N*th pixel of the corridor |
| `--timeseries` | — | time-series file(s) for the 3rd panel |
| `--start-yx` / `--end-yx` | — | line end points in pixel numbers |
| `--start-lalo` / `--end-lalo` | — | line end points in lat/lon |
| `--line-file` | — | GMT-format `lon lat` file with the two end points |
| `--save` / `--nodisplay` / `-o` / `--outdir` | — | output control |

Run `viewinsar swath -h` for the full list.

---

## Output

Pressing <kbd>e</kbd> (or passing `--save`) writes four files per profile:

**`swath_..._points.csv`** — one row per pixel, the linear-referencing table:

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
figure draws. The plain statistics are exact and un-clipped; the
density-weighted ones are clipped to the plotted Y range, exactly like the
histogram.

**`swath_..._info.txt`** — everything needed to reproduce the profile.
**`swath_....png`** — the figure.

### Referencing the profile to a fault, afterwards

By default *x* = 0 is the start of the line. To re-zero on a fault you already
know the position of:

```python
import pandas as pd

df = pd.read_csv('swath_Y300X150_Y300X750_w15km_points.csv')
df['dist_from_fault'] = df['dist_km'] - 27.5     # your fault position, in km
df['vel'] = df['value']                          # already in the display unit
```

or do it up front with `--x-offset 27.5`, which also puts a dashed marker at
*x* = 0. The density grid always places a cell **boundary** at *x* = 0, so the
bins never straddle the fault.

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

`viewinsar.swathProfileViewer` is the interactive figure; it is imported lazily,
so `import viewinsar` stays cheap.

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

## License

GPL-3.0-or-later, following MintPy.

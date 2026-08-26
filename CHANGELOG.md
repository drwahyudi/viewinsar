# Changelog

All notable changes to viewinsar are recorded here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project uses [semantic versioning](https://semver.org/).

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

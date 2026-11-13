"""Radar-coded data: the MintPy geometry file as the bridge to lon/lat.

A MintPy product in radar coordinates (no X_FIRST / Y_FIRST) has a companion
geometry file - inputs/geometryRadar.h5 - on the very same grid, holding for
every pixel its `latitude`, `longitude` and `height` (the DEM used in the
processing). That is all the swath viewer needs to put lon/lat things - fault
traces, a DEM file, downloaded terrain - onto the radar grid:

    RadarLookup.to_lonlat(rows, cols)  radar pixel -> lon/lat (bilinear)
    RadarLookup.to_rowcol(lon, lat)    lon/lat -> radar pixel (piecewise
                                       linear over a Delaunay triangulation of
                                       a thinned lookup grid)

find_lookup() searches for the file the way MintPy does, but relative to the
data file as well as to the current folder.
"""

import os

import numpy as np

__all__ = ['find_lookup', 'RadarLookup', 'PIXEL_ATR']

# a "map" whose coordinates are the pixel numbers (col, row), pixel centres on
# integers - lets the map-coordinate helpers work on a radar grid unchanged
PIXEL_ATR = {'X_FIRST': '-0.5', 'Y_FIRST': '-0.5', 'X_STEP': '1.0', 'Y_STEP': '1.0'}

LOOKUP_NAMES = ('geometryRadar.h5',)
LOOKUP_DIRS = ('inputs', '', os.path.join('..', 'inputs'))


def _datasets(fname):
    try:
        import h5py
        with h5py.File(fname, 'r') as fid:
            return {name for name, obj in fid.items() if hasattr(obj, 'shape')}
    except Exception:            # noqa: BLE001 - not an HDF5 file
        return set()


def check_lookup(fname, atr):
    """Why `fname` cannot serve as the lookup for data with metadata `atr`, or None."""
    from .compat import read_attribute

    if not os.path.isfile(fname):
        return f'no such file: {fname}'
    try:
        latr = read_attribute(fname)
    except Exception as err:     # noqa: BLE001
        return f'cannot read {fname}: {err}'
    if 'Y_FIRST' in latr:
        return f'{os.path.basename(fname)} is geocoded; a radar-coded one is needed'
    missing = {'latitude', 'longitude'} - _datasets(fname)
    if missing:
        return f'{os.path.basename(fname)} has no {" / ".join(sorted(missing))} dataset'
    if (latr.get('LENGTH'), latr.get('WIDTH')) != (atr.get('LENGTH'), atr.get('WIDTH')):
        return (f"{os.path.basename(fname)} is {latr.get('LENGTH')} x {latr.get('WIDTH')} "
                f"pixels, the data {atr.get('LENGTH')} x {atr.get('WIDTH')}")
    return None


def find_lookup(data_file, atr, given=None, verbose=True):
    """The radar-coded geometry file for `data_file`, or None.

    `given` (--lookup) is used as it is. Otherwise inputs/geometryRadar.h5,
    ./geometryRadar.h5 and ../inputs/geometryRadar.h5 are tried next to the
    data file, then in the current folder - MintPy's own search order.
    """
    if given:
        problem = check_lookup(given, atr)
        if problem:
            raise ValueError(f'--lookup: {problem}')
        return given

    seen, rejected = set(), []
    for base in (os.path.dirname(os.path.abspath(data_file)), os.getcwd()):
        for sub in LOOKUP_DIRS:
            for name in LOOKUP_NAMES:
                path = os.path.normpath(os.path.join(base, sub, name))
                if path in seen or not os.path.isfile(path):
                    continue
                seen.add(path)
                problem = check_lookup(path, atr)
                if problem is None:
                    if verbose:
                        print(f'lookup: {path}')
                    return path
                rejected.append(problem)
    if verbose and rejected:
        print('lookup: found but not usable - ' + '; '.join(rejected))
    return None


def has_height(fname):
    return 'height' in _datasets(fname)


##############################################################################
class RadarLookup:
    """lon/lat <-> radar pixel, from a MintPy geometryRadar.h5."""

    def __init__(self, fname):
        from .compat import read

        self.file = fname
        self.lat = read(fname, datasetName='latitude', print_msg=False)[0]
        self.lon = read(fname, datasetName='longitude', print_msg=False)[0]
        self.lat = np.asarray(self.lat, dtype=np.float64)
        self.lon = np.asarray(self.lon, dtype=np.float64)
        # no-data: NaN, or 0 / 0 as some processors write it
        bad = ~(np.isfinite(self.lat) & np.isfinite(self.lon)) | ((self.lat == 0) & (self.lon == 0))
        self.lat[bad] = np.nan
        self.lon[bad] = np.nan
        self.length, self.width = self.lat.shape
        self._interp = None

    def extent(self, margin=0.0):
        """(W, S, E, N) of the radar frame, in degrees."""
        west, east = np.nanmin(self.lon), np.nanmax(self.lon)
        south, north = np.nanmin(self.lat), np.nanmax(self.lat)
        pad_x, pad_y = (east - west) * margin, (north - south) * margin
        return west - pad_x, south - pad_y, east + pad_x, north + pad_y

    def to_lonlat(self, rows, cols):
        """lon/lat at fractional radar pixels (bilinear; NaN off the frame)."""
        from .terrain import _bilinear

        rows, cols = np.asarray(rows, dtype=np.float64), np.asarray(cols, dtype=np.float64)
        return _bilinear(self.lon, rows, cols), _bilinear(self.lat, rows, cols)

    def to_rowcol(self, lon, lat):
        """Fractional radar pixels of lon/lat points (NaN outside the frame)."""
        if self._interp is None:
            from scipy.interpolate import LinearNDInterpolator

            # a thinned grid keeps the triangulation quick; the lookup is smooth
            step = max(1, int(np.ceil(np.sqrt(self.length * self.width / 150_000))))
            rr, cc = np.mgrid[0:self.length:step, 0:self.width:step]
            # always include the last row / column, so the frame edge is covered
            if rr[-1, 0] != self.length - 1 or cc[0, -1] != self.width - 1:
                rows = np.r_[np.arange(0, self.length, step), self.length - 1]
                cols = np.r_[np.arange(0, self.width, step), self.width - 1]
                rr, cc = np.meshgrid(np.unique(rows), np.unique(cols), indexing='ij')
            lon_s, lat_s = self.lon[rr, cc], self.lat[rr, cc]
            ok = np.isfinite(lon_s) & np.isfinite(lat_s)
            self._interp = LinearNDInterpolator(np.c_[lon_s[ok], lat_s[ok]],
                                                np.c_[rr[ok], cc[ok]].astype(np.float64))
        out = self._interp(np.c_[np.ravel(lon), np.ravel(lat)])
        shape = np.shape(lon)
        return out[:, 0].reshape(shape), out[:, 1].reshape(shape)

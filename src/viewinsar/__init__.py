"""viewinsar - interactive swath (buffered) profile viewer for InSAR data.

Built on MintPy: the map panel comes from ``mintpy.view``, the click-to-draw
interaction from ``mintpy.plot_transection``, and the optional time-series panel
from ``mintpy.tsview``. What it adds is the swath itself - a flat-capped buffer
around the drawn line, linear referencing of every pixel inside it, and a 2-D
point-density profile with the density-weighted median and quartiles. And
``viewinsar stat``: block statistics in pairs of circles on both sides of a
line drawn along a fault, or their difference (-rate); ``viewinsar tseries``:
LiCSBAS_plot_ts.py for LiCSBAS and MintPy time series, with the median of a
circle at every date; ``viewinsar model``: the slip rate and locking depth of
an elastic dislocation (and shallow creep) fitted to a profile across a fault
by ensemble MCMC, after Hussain et al. (2016, 2018).

    >>> from viewinsar import swath_extract, swath_density, stat_line, fit_profile

The heavy imports (mintpy, matplotlib) live in ``viewinsar.viewer`` and are only
pulled in when the viewer is used, so ``viewinsar.swath`` stays importable with
numpy alone.
"""

__version__ = '0.4.1'
__author__ = 'Donny Wahyudi'
__license__ = 'GPL-3.0-or-later'

from .fault import load_faults, profile_crossings, read_fault_file
from .gnss import read_gnss
from .model import fit_profile
from .radar import RadarLookup, find_lookup
from .stat import stat_line
from .swath import (
    EARTH_RADIUS,
    bilinear_sample,
    get_pixel_size,
    is_projected,
    line_geometry,
    nice_step,
    swath_density,
    swath_extract,
    weighted_quantile,
)
from .terrain import terrain_dem

__all__ = [
    '__version__',
    'EARTH_RADIUS',
    'RadarLookup',
    'bilinear_sample',
    'find_lookup',
    'fit_profile',
    'get_pixel_size',
    'is_projected',
    'line_geometry',
    'load_faults',
    'nice_step',
    'profile_crossings',
    'read_fault_file',
    'read_gnss',
    'stat_line',
    'swath_density',
    'swath_extract',
    'terrain_dem',
    'weighted_quantile',
    'swathProfileViewer',
    'StatViewer',
    'TimeSeriesViewer',
    'ModelViewer',
]


def __getattr__(name):
    """Lazily expose the viewers, so `import viewinsar` stays light."""
    if name == 'swathProfileViewer':
        from .viewer import swathProfileViewer
        return swathProfileViewer
    if name == 'StatViewer':
        from .statview import StatViewer
        return StatViewer
    if name == 'TimeSeriesViewer':
        from .tsviewer import TimeSeriesViewer
        return TimeSeriesViewer
    if name == 'ModelViewer':
        from .modelview import ModelViewer
        return ModelViewer
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

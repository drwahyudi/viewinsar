"""viewinsar - interactive swath (buffered) profile viewer for InSAR data.

Built on MintPy: the map panel comes from ``mintpy.view``, the click-to-draw
interaction from ``mintpy.plot_transection``, and the optional time-series panel
from ``mintpy.tsview``. What it adds is the swath itself - a flat-capped buffer
around the drawn line, linear referencing of every pixel inside it, and a 2-D
point-density profile with the density-weighted median and quartiles.

    >>> from viewinsar import swath_extract, swath_density

The heavy imports (mintpy, matplotlib) live in ``viewinsar.viewer`` and are only
pulled in when the viewer is used, so ``viewinsar.swath`` stays importable with
numpy alone.
"""

__version__ = '0.1.0'
__author__ = 'Donny Wahyudi'
__license__ = 'GPL-3.0-or-later'

from .swath import (
    EARTH_RADIUS,
    get_pixel_size,
    line_geometry,
    nice_step,
    swath_density,
    swath_extract,
    weighted_quantile,
)

__all__ = [
    '__version__',
    'EARTH_RADIUS',
    'get_pixel_size',
    'line_geometry',
    'nice_step',
    'swath_density',
    'swath_extract',
    'weighted_quantile',
    'swathProfileViewer',
]


def __getattr__(name):
    """Lazily expose the viewer, so `import viewinsar` stays light."""
    if name == 'swathProfileViewer':
        from .viewer import swathProfileViewer
        return swathProfileViewer
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

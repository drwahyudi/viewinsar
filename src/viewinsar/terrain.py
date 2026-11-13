"""A DEM for the map extent, downloaded from the Terrain Tiles.

`-d` / `--topoxsect` given WITHOUT a file fetch the elevation from the Terrain
Tiles on AWS Open Data - the "Mapzen Global Terrain" layer of QGIS:

    https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png

256 x 256 PNG tiles in Web Mercator (EPSG:3857), zoom 0-15, the standard XYZ
scheme (y = 0 at the north edge), the height packed into the colour:

    height [m] = R * 256 + G + B / 256 - 32768

The tiles over the map are mosaicked, re-gridded into the DATA's own CRS (so
MintPy's view.py can draw it as the shaded-relief background, and the swath
can sample it like any DEM) and written as a MintPy geometry file - `height`,
in metres - in the cache folder, along with the tiles themselves. A second run
over the same area downloads nothing.

Sources: SRTM, GMTED2010 and ETOPO1 (which includes bathymetry) worldwide,
finer national DEMs where they exist. Using them requires attribution:
ATTRIBUTION below, also written into the DEM file.
"""

import concurrent.futures
import hashlib
import io
import os
import time
import urllib.error
import urllib.request

import numpy as np

__all__ = ['TERRARIUM_URL', 'ATTRIBUTION', 'AUTO_DEM', 'is_auto', 'terrain_dem',
           'choose_zoom', 'lonlat_to_tile', 'decode_terrarium', 'tile_resolution',
           'cache_dir', 'terrain_dem_radar']

AUTO_DEM = 'auto'                            # -d / --topoxsect given without a file
TERRARIUM_URL = 'https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png'
ATTRIBUTION = ('Terrain Tiles (Mapzen/Tilezen, AWS Open Data): SRTM, GMTED2010 and ETOPO1 '
               'courtesy of the U.S. Geological Survey and NOAA, and others - '
               'https://github.com/tilezen/joerd/blob/master/docs/attribution.md')

TILE = 256                                   # pixels per tile side
HALF = 20037508.342789244                    # pi * 6378137: half the Web Mercator world
MAX_LAT = 85.0511287798066                   # the Web Mercator limit
MAX_TILES = 300                              # ~30 MB: the default download budget
MAX_ZOOM_AUTO = 12                           # ~38 m at the equator, about SRTM's 30 m
MAX_SIDE = 8000                              # pixels per side of the re-gridded DEM


def is_auto(value):
    """Was -d / --topoxsect given without a file (or as 'auto')?"""
    return isinstance(value, str) and value.strip().lower() == AUTO_DEM


##############################################################################
#                               tile maths                                   #
##############################################################################
def tile_resolution(zoom, lat=0.0):
    """Ground size of a tile pixel, in metres, at a latitude."""
    return 2 * HALF / (TILE * 2 ** zoom) * np.cos(np.radians(lat))


def lonlat_to_tile(lon, lat, zoom):
    """Fractional XYZ tile coordinates (x east, y south) of lon/lat."""
    num = 2 ** zoom
    lat = np.radians(np.clip(lat, -MAX_LAT, MAX_LAT))
    x = (np.asarray(lon, dtype=np.float64) + 180.0) / 360.0 * num
    y = (1.0 - np.log(np.tan(lat) + 1.0 / np.cos(lat)) / np.pi) / 2.0 * num
    return x, y


def _tile_range(bbox, zoom):
    """First / last tile index along x and y covering a lon/lat box (W, S, E, N)."""
    west, south, east, north = bbox
    x0, y0 = lonlat_to_tile(west, north, zoom)
    x1, y1 = lonlat_to_tile(east, south, zoom)
    last = 2 ** zoom - 1
    return (int(np.clip(np.floor(x0), 0, last)), int(np.clip(np.floor(x1), 0, last)),
            int(np.clip(np.floor(y0), 0, last)), int(np.clip(np.floor(y1), 0, last)))


def choose_zoom(bbox, max_tiles=MAX_TILES, max_zoom=MAX_ZOOM_AUTO):
    """The finest zoom whose tiles over `bbox` stay within the budget."""
    for zoom in range(int(max_zoom), -1, -1):
        tx0, tx1, ty0, ty1 = _tile_range(bbox, zoom)
        if (tx1 - tx0 + 1) * (ty1 - ty0 + 1) <= max_tiles:
            return zoom
    return 0


def decode_terrarium(png):
    """Height in metres from the bytes of a Terrarium PNG tile."""
    from PIL import Image  # Pillow, a matplotlib dependency

    with Image.open(io.BytesIO(png)) as img:
        rgb = np.asarray(img.convert('RGB'), dtype=np.float64)
    return (rgb[..., 0] * 256.0 + rgb[..., 1] + rgb[..., 2] / 256.0 - 32768.0).astype(np.float32)


##############################################################################
#                               download                                     #
##############################################################################
def cache_dir():
    """Where tiles and DEMs are kept: $VIEWINSAR_CACHE, else ~/.cache/viewinsar."""
    return os.environ.get('VIEWINSAR_CACHE') or os.path.join(
        os.path.expanduser('~'), '.cache', 'viewinsar')


def _get(url, timeout=30, retries=2):
    """The body of `url`, or None for a tile that does not exist (404)."""
    from . import __version__

    req = urllib.request.Request(url, headers={'User-Agent': f'viewinsar/{__version__}'})
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code == 404:
                return None
            if attempt == retries or err.code < 500:
                raise
        except (urllib.error.URLError, TimeoutError, OSError):
            if attempt == retries:
                raise
        time.sleep(1.0 + attempt)
    return None


def fetch_tile(zoom, x, y, folder=None, url=TERRARIUM_URL):
    """One tile's PNG bytes, from the cache or the web (None if it does not exist)."""
    folder = folder or cache_dir()
    path = os.path.join(folder, 'terrarium', str(zoom), str(x), f'{y}.png')
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        with open(path, 'rb') as fid:
            return fid.read()
    png = _get(url.format(z=zoom, x=x, y=y))
    if png:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = f'{path}.{os.getpid()}.part'
        with open(tmp, 'wb') as fid:
            fid.write(png)
        os.replace(tmp, path)                # never leave a half-written tile behind
    return png


def mosaic_tiles(zoom, tx0, tx1, ty0, ty1, fetch=None, workers=8, verbose=True):
    """Download and decode a block of tiles into one array (NaN where missing)."""
    fetch = fetch or fetch_tile
    nx, ny = tx1 - tx0 + 1, ty1 - ty0 + 1
    out = np.full((ny * TILE, nx * TILE), np.nan, dtype=np.float32)
    jobs = [(x, y) for y in range(ty0, ty1 + 1) for x in range(tx0, tx1 + 1)]
    failed, done = [], 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch, zoom, x, y): (x, y) for x, y in jobs}
        for fut in concurrent.futures.as_completed(futures):
            x, y = futures[fut]
            try:
                png = fut.result()
            except Exception as err:         # noqa: BLE001 - reported below
                failed.append(err)
                png = None
            if png:
                row, col = (y - ty0) * TILE, (x - tx0) * TILE
                out[row:row + TILE, col:col + TILE] = decode_terrarium(png)
            done += 1
            if verbose and len(jobs) >= 20 and done % max(1, len(jobs) // 5) == 0:
                print(f'   ... {done}/{len(jobs)} tiles')
    if failed and len(failed) == len(jobs):
        raise RuntimeError(f'could not download any Terrain Tile ({failed[0]}). Is there '
                           'internet access (or HTTPS_PROXY)? A DEM file works too: -d dem.tif')
    if failed and verbose:
        print(f'   WARNING: {len(failed)} of {len(jobs)} tiles failed ({failed[0]}); '
              'those areas are left empty.')
    return out


##############################################################################
#                     onto the data's grid                                   #
##############################################################################
def _bilinear(img, rows, cols):
    """Bilinear samples of img at arbitrary fractional (row, col); NaN off the image."""
    nrow, ncol = img.shape
    finite = np.isfinite(rows) & np.isfinite(cols)        # a NaN would index garbage
    inside = (finite & (rows >= -0.5) & (rows <= nrow - 0.5)
              & (cols >= -0.5) & (cols <= ncol - 0.5))
    rows = np.where(finite, np.clip(rows, 0, nrow - 1), 0.0)
    cols = np.where(finite, np.clip(cols, 0, ncol - 1), 0.0)
    r0 = np.minimum(np.floor(rows).astype(np.intp), max(nrow - 2, 0))
    c0 = np.minimum(np.floor(cols).astype(np.intp), max(ncol - 2, 0))
    fr, fc = rows - r0, cols - c0
    r1, c1 = np.minimum(r0 + 1, nrow - 1), np.minimum(c0 + 1, ncol - 1)
    val = (img[r0, c0] * (1 - fr) * (1 - fc) + img[r0, c1] * (1 - fr) * fc
           + img[r1, c0] * fr * (1 - fc) + img[r1, c1] * fr * fc)
    return np.where(inside, val, np.nan)


def _lonlat_box(geo_box, data_crs):
    """lon/lat bounding box (W, S, E, N) of a map box (W, N, E, S) in data_crs."""
    from pyproj import Transformer

    west, north, east, south = geo_box
    xs = np.r_[np.linspace(west, east, 21), np.full(21, east),
               np.linspace(east, west, 21), np.full(21, west)]
    ys = np.r_[np.full(21, north), np.linspace(north, south, 21),
               np.full(21, south), np.linspace(south, north, 21)]
    lon, lat = Transformer.from_crs(data_crs, 4326, always_xy=True).transform(xs, ys)
    ok = np.isfinite(lon) & np.isfinite(lat)
    if not ok.any():
        raise ValueError('cannot place the map extent on the globe (CRS?)')
    return (float(lon[ok].min()), float(max(lat[ok].min(), -MAX_LAT)),
            float(lon[ok].max()), float(min(lat[ok].max(), MAX_LAT)))


def terrain_dem(geo_box, data_crs, projected, zoom=None, max_tiles=MAX_TILES, folder=None,
                fetch=None, margin=0.02, verbose=True):
    """Build (or reuse) a DEM covering a map extent, on the data's own CRS.

    Parameters: geo_box   - (W, N, E, S) of the map, in the data's coordinates
                data_crs  - the data's CRS: EPSG int or WKT
                projected - bool, the data's grid is in metres (else degrees)
                zoom      - int, tile zoom 0-15 (default: the finest within max_tiles)
                max_tiles - int, download budget for the automatic zoom
                folder    - str, cache folder (default: cache_dir())
                fetch     - callable(zoom, x, y) -> PNG bytes, to override the download
                margin    - float, extra border around the map, as a fraction of its size
    Returns:    str, path to a MintPy geometry .h5 with `height` in metres
    """
    import h5py
    from pyproj import Transformer

    folder = folder or cache_dir()
    west, north, east, south = (float(v) for v in geo_box)
    pad_x, pad_y = (east - west) * margin, (north - south) * margin
    box = (west - pad_x, north + pad_y, east + pad_x, south - pad_y)
    bbox = _lonlat_box(box, data_crs)
    if zoom is None:
        zoom = choose_zoom(bbox, max_tiles=max_tiles)
    zoom = int(np.clip(zoom, 0, 15))
    lat_mid = 0.5 * (bbox[1] + bbox[3])
    res_m = float(tile_resolution(zoom, lat_mid))

    # the output grid: square pixels about the tile resolution, in the data's units
    step = res_m if projected else res_m / 111320.0
    width = int(np.ceil((box[2] - box[0]) / step))
    length = int(np.ceil((box[1] - box[3]) / step))
    if max(width, length) > MAX_SIDE:
        step *= max(width, length) / MAX_SIDE
        width = int(np.ceil((box[2] - box[0]) / step))
        length = int(np.ceil((box[1] - box[3]) / step))

    key = hashlib.sha1(repr((zoom, str(data_crs), [round(v, 6) for v in box],
                             round(step, 12))).encode()).hexdigest()[:10]
    out = os.path.join(folder, f'dem_terrarium_z{zoom}_{key}.h5')
    tx0, tx1, ty0, ty1 = _tile_range(bbox, zoom)
    if verbose:
        print(f'terrain: Terrain Tiles at zoom {zoom} (~{res_m:.0f} m), '
              f'{(tx1 - tx0 + 1) * (ty1 - ty0 + 1)} tiles over lon {bbox[0]:.3f}..{bbox[2]:.3f}'
              f', lat {bbox[1]:.3f}..{bbox[3]:.3f}')
    if os.path.isfile(out):
        if verbose:
            print(f'terrain: re-using {out}')
        return out

    tiles = mosaic_tiles(zoom, tx0, tx1, ty0, ty1, fetch=fetch, verbose=verbose)

    # every output pixel centre -> Web Mercator -> its place in the mosaic
    to_merc = Transformer.from_crs(data_crs, 3857, always_xy=True)
    xs = box[0] + (np.arange(width) + 0.5) * step
    dem = np.empty((length, width), dtype=np.float32)
    for r0 in range(0, length, 512):
        r1 = min(r0 + 512, length)
        ys = box[1] - (np.arange(r0, r1) + 0.5) * step
        gx, gy = np.meshgrid(xs, ys)
        dem[r0:r1] = _sample_mosaic(tiles, zoom, tx0, ty0, *to_merc.transform(gx, gy))

    os.makedirs(folder, exist_ok=True)
    atr = {
        'FILE_TYPE': 'geometry', 'UNIT': 'm', 'DATA_TYPE': 'float32',
        'LENGTH': str(length), 'WIDTH': str(width),
        'X_FIRST': repr(box[0]), 'Y_FIRST': repr(box[1]),
        'X_STEP': repr(step), 'Y_STEP': repr(-step),
        'X_UNIT': 'meters' if projected else 'degrees',
        'Y_UNIT': 'meters' if projected else 'degrees',
        'PROCESSOR': 'viewinsar', 'SOURCE': TERRARIUM_URL, 'TILE_ZOOM': str(zoom),
        'ATTRIBUTION': ATTRIBUTION,
    }
    if isinstance(data_crs, int) or str(data_crs).isdigit():
        atr['EPSG'] = str(int(data_crs))
        if 32601 <= int(data_crs) <= 32660 or 32701 <= int(data_crs) <= 32760:
            atr['UTM_ZONE'] = f"{int(data_crs) % 100}{'N' if int(data_crs) < 32700 else 'S'}"
    else:
        atr['CRS_WKT'] = str(data_crs)
    tmp = f'{out}.{os.getpid()}.part'
    with h5py.File(tmp, 'w') as fid:
        fid.create_dataset('height', data=dem, compression='gzip', compression_opts=4)
        for name, value in atr.items():
            fid.attrs[name] = value
    os.replace(tmp, out)
    if verbose:
        finite = np.isfinite(dem)
        rng = (f'{np.nanmin(dem):.0f} to {np.nanmax(dem):.0f} m' if finite.any()
               else 'no data')
        print(f'terrain: wrote {out} ({length} x {width}, {rng})')
    return out


def _sample_mosaic(tiles, zoom, tx0, ty0, mx, my):
    """The mosaic, bilinear, at Web Mercator points (NaN off the mosaic)."""
    scale = TILE * 2 ** zoom / (2 * HALF)                  # mosaic pixels per metre
    cols = (np.asarray(mx, dtype=np.float64) + HALF) * scale - tx0 * TILE - 0.5
    rows = (HALF - np.asarray(my, dtype=np.float64)) * scale - ty0 * TILE - 0.5
    return _bilinear(tiles, rows, cols)


def terrain_dem_radar(lookup, atr, zoom=None, max_tiles=MAX_TILES, folder=None, fetch=None,
                      verbose=True):
    """Terrain Tiles on a RADAR grid, through the lookup's lat/lon of every pixel.

    For radar-coded data whose geometry file has no `height` of its own.

    Parameters: lookup - radar.RadarLookup of the data
                atr    - dict, the data's metadata (copied into the DEM file, so
                         view.py and the swath read it on the same grid)
    Returns:    str, path to a radar-coded MintPy geometry .h5 with `height`
    """
    import h5py
    from pyproj import Transformer

    folder = folder or cache_dir()
    bbox = lookup.extent(margin=0.01)
    if zoom is None:
        zoom = choose_zoom(bbox, max_tiles=max_tiles)
    zoom = int(np.clip(zoom, 0, 15))
    stat = os.stat(lookup.file)
    key = hashlib.sha1(repr((zoom, os.path.abspath(lookup.file), stat.st_size,
                             int(stat.st_mtime))).encode()).hexdigest()[:10]
    out = os.path.join(folder, f'dem_terrarium_radar_z{zoom}_{key}.h5')
    tx0, tx1, ty0, ty1 = _tile_range(bbox, zoom)
    if verbose:
        res_m = float(tile_resolution(zoom, 0.5 * (bbox[1] + bbox[3])))
        print(f'terrain: Terrain Tiles at zoom {zoom} (~{res_m:.0f} m), '
              f'{(tx1 - tx0 + 1) * (ty1 - ty0 + 1)} tiles over the radar frame, '
              'resampled onto it through the lookup')
    if os.path.isfile(out):
        if verbose:
            print(f'terrain: re-using {out}')
        return out

    tiles = mosaic_tiles(zoom, tx0, tx1, ty0, ty1, fetch=fetch, verbose=verbose)
    to_merc = Transformer.from_crs(4326, 3857, always_xy=True)
    dem = np.empty((lookup.length, lookup.width), dtype=np.float32)
    for r0 in range(0, lookup.length, 512):
        r1 = min(r0 + 512, lookup.length)
        dem[r0:r1] = _sample_mosaic(tiles, zoom, tx0, ty0,
                                    *to_merc.transform(lookup.lon[r0:r1], lookup.lat[r0:r1]))

    os.makedirs(folder, exist_ok=True)
    meta = {key: str(val) for key, val in atr.items()}
    meta.update({
        'FILE_TYPE': 'geometry', 'UNIT': 'm', 'DATA_TYPE': 'float32',
        'LENGTH': str(lookup.length), 'WIDTH': str(lookup.width),
        'PROCESSOR': meta.get('PROCESSOR', 'viewinsar'), 'SOURCE': TERRARIUM_URL,
        'TILE_ZOOM': str(zoom), 'ATTRIBUTION': ATTRIBUTION,
        'LOOKUP_FILE': os.path.abspath(lookup.file),
    })
    tmp = f'{out}.{os.getpid()}.part'
    with h5py.File(tmp, 'w') as fid:
        fid.create_dataset('height', data=dem, compression='gzip', compression_opts=4)
        for name, value in meta.items():
            fid.attrs[name] = value
    os.replace(tmp, out)
    if verbose:
        rng = (f'{np.nanmin(dem):.0f} to {np.nanmax(dem):.0f} m' if np.isfinite(dem).any()
               else 'no data')
        print(f'terrain: wrote {out} ({lookup.length} x {lookup.width}, radar grid, {rng})')
    return out

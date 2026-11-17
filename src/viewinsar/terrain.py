"""A DEM for the map extent, downloaded: Copernicus GLO-30, or the Terrain Tiles.

`--dem` / `--topoxsect` given WITHOUT a file fetch the elevation for the map.

**Copernicus DEM GLO-30** (the default, `--dem-source cop30`) - the 1 arc-second
(~30 m) DEM of the TanDEM-X mission, heights above the EGM2008 geoid, as Cloud
Optimized GeoTIFFs on AWS Open Data, one per 1 x 1 degree cell named by its
south-west corner:

    https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif
    name = Copernicus_DSM_COG_10_S01_00_E119_00_DEM   (lat -1..0, lon 119..120)

Each is 3600 rows of 1"; the columns widen with latitude (3600 below 50 deg, 2400
to 60, 1800 to 70 ...), in 1024 x 1024 DEFLATE tiles with overviews at 1/2,
1/4 and 1/8 resolution. Only the parts over the map are read - by HTTP byte
range, at the overview nearest the grid being made - so no GDAL is needed and a
big map downloads its relief at the resolution it is drawn at. The sea has no
tiles: 0 m there, as the product specifies (so do the few countries not
released, Armenia and Azerbaijan). Licence: free, with the notice
COP30_ATTRIBUTION, also written into the DEM file.

**Terrain Tiles** (`--dem-source terrain`, and the fallback when GLO-30 cannot be
had) - the "Mapzen Global Terrain" layer of QGIS on AWS Open Data:

    https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png

256 x 256 PNG tiles in Web Mercator (EPSG:3857), zoom 0-15, the standard XYZ
scheme (y = 0 at the north edge), the height packed into the colour:

    height [m] = R * 256 + G + B / 256 - 32768

Sources: SRTM, GMTED2010 and ETOPO1 (which includes bathymetry) worldwide,
finer national DEMs where they exist; attribution ATTRIBUTION.

Either way the elevation is re-gridded into the DATA's own CRS (so MintPy's
view.py can draw it as the shaded-relief background, and the swath can sample
it like any DEM) and written as a MintPy geometry file - `height`, in metres -
in the cache folder, along with what was downloaded. A second run over the same
area downloads nothing.
"""

import concurrent.futures
import hashlib
import io
import os
import struct
import time
import urllib.error
import urllib.request
import zlib

import numpy as np

__all__ = ['TERRARIUM_URL', 'ATTRIBUTION', 'AUTO_DEM', 'is_auto', 'terrain_dem',
           'choose_zoom', 'lonlat_to_tile', 'decode_terrarium', 'tile_resolution',
           'cache_dir', 'terrain_dem_radar', 'COP30_URL', 'COP30_ATTRIBUTION',
           'DEM_SOURCES', 'cop30_name', 'copernicus_dem', 'copernicus_dem_radar']

AUTO_DEM = 'auto'                            # --dem / --topoxsect given without a file
TERRARIUM_URL = 'https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png'
ATTRIBUTION = ('Terrain Tiles (Mapzen/Tilezen, AWS Open Data): SRTM, GMTED2010 and ETOPO1 '
               'courtesy of the U.S. Geological Survey and NOAA, and others - '
               'https://github.com/tilezen/joerd/blob/master/docs/attribution.md')

COP30_URL = 'https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif'
# the notice of the Copernicus DEM licence (art. 6b: data re-gridded from it)
COP30_ATTRIBUTION = ('Copernicus DEM GLO-30: produced using Copernicus WorldDEM-30 (c) DLR '
                     'e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided '
                     'under COPERNICUS by the European Union and ESA; all rights reserved - '
                     'https://registry.opendata.aws/copernicus-dem')
DEM_SOURCES = ('cop30', 'terrain')           # --dem-source: Copernicus GLO-30, Terrain Tiles

TILE = 256                                   # pixels per tile side
HALF = 20037508.342789244                    # pi * 6378137: half the Web Mercator world
MAX_LAT = 85.0511287798066                   # the Web Mercator limit
MAX_TILES = 300                              # ~30 MB: the default download budget
MAX_ZOOM_AUTO = 12                           # ~38 m at the equator, about SRTM's 30 m
MAX_SIDE = 8000                              # pixels per side of the re-gridded DEM
MAP_SIDE = 3000                              # ... of the one for the map background only

COP30_RES = 1.0 / 3600.0                     # GLO-30 rows: 1 arc-second, in degrees
COP30_MAX_TILES = 400                        # 1 x 1 deg cells: a map of ~20 x 20 deg
HEAD_BYTES = 65536                           # a COG's headers are in its first bytes
M_PER_DEG = 111320.0                         # ~ a degree of latitude, in metres


def is_auto(value):
    """Was --dem / --topoxsect given without a file (or as 'auto')?"""
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


def _get(url, timeout=30, retries=2, byte_range=None):
    """The body of `url` - or of its bytes byte_range = (first, last) - or None for a
    file that does not exist (404, or S3's own 403 for a key it will not list)."""
    from . import __version__

    headers = {'User-Agent': f'viewinsar/{__version__}'}
    if byte_range is not None:
        headers['Range'] = f'bytes={int(byte_range[0])}-{int(byte_range[1])}'
    req = urllib.request.Request(url, headers=headers)
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read()
                if byte_range is not None and resp.status == 200:     # the range ignored
                    body = body[int(byte_range[0]):int(byte_range[1]) + 1]
                return body
        except urllib.error.HTTPError as err:
            if err.code == 404 or (err.code == 403 and 'AmazonS3' in str(
                    err.headers.get('Server', ''))):
                return None
            if err.code == 416:                  # a range past the end of the file
                return b''
            if attempt == retries or err.code < 500:
                raise
        except (urllib.error.URLError, TimeoutError, OSError):
            if attempt == retries:
                raise
        time.sleep(1.0 + attempt)
    return None


def fetch_range(url, start, size):
    """`size` bytes of `url` from byte `start` (fewer at its end); None if it does
    not exist. Every GLO-30 byte comes through here (tests replace it)."""
    return _get(url, byte_range=(start, start + size - 1))


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
                           'internet access (or HTTPS_PROXY)? A DEM file works too: --dem dem.tif')
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


def _outline(west, north, east, south, num=21):
    """x, y of `num` points along each side of a box."""
    xs = np.r_[np.linspace(west, east, num), np.full(num, east),
               np.linspace(east, west, num), np.full(num, west)]
    ys = np.r_[np.full(num, north), np.linspace(north, south, num),
               np.full(num, south), np.linspace(south, north, num)]
    return xs, ys


def _lonlat_box(geo_box, data_crs, limit=MAX_LAT):
    """lon/lat bounding box (W, S, E, N) of a map box (W, N, E, S) in data_crs."""
    from pyproj import Transformer

    lon, lat = Transformer.from_crs(data_crs, 4326, always_xy=True).transform(
        *_outline(*geo_box))
    ok = np.isfinite(lon) & np.isfinite(lat)
    if not ok.any():
        raise ValueError('cannot place the map extent on the globe (CRS?)')
    return (float(lon[ok].min()), float(max(lat[ok].min(), -limit)),
            float(lon[ok].max()), float(min(lat[ok].max(), limit)))


def _grid_size(box, step, max_side):
    """(width, length, step) of a grid of square `step` pixels over a box
    (W, N, E, S), the step widened so neither side passes max_side pixels."""
    width = int(np.ceil((box[2] - box[0]) / step))
    length = int(np.ceil((box[1] - box[3]) / step))
    if max(width, length) > max_side:
        step *= max(width, length) / max_side
        width = int(np.ceil((box[2] - box[0]) / step))
        length = int(np.ceil((box[1] - box[3]) / step))
    return max(width, 1), max(length, 1), step


def _geo_attrs(box, step, length, width, projected, data_crs):
    """The MintPy attributes of a DEM grid on the data's CRS."""
    atr = {
        'FILE_TYPE': 'geometry', 'UNIT': 'm', 'DATA_TYPE': 'float32',
        'LENGTH': str(length), 'WIDTH': str(width),
        'X_FIRST': repr(box[0]), 'Y_FIRST': repr(box[1]),
        'X_STEP': repr(step), 'Y_STEP': repr(-step),
        'X_UNIT': 'meters' if projected else 'degrees',
        'Y_UNIT': 'meters' if projected else 'degrees',
        'PROCESSOR': 'viewinsar',
    }
    if isinstance(data_crs, int) or str(data_crs).isdigit():
        atr['EPSG'] = str(int(data_crs))
        if 32601 <= int(data_crs) <= 32660 or 32701 <= int(data_crs) <= 32760:
            atr['UTM_ZONE'] = f"{int(data_crs) % 100}{'N' if int(data_crs) < 32700 else 'S'}"
    else:
        atr['CRS_WKT'] = str(data_crs)
    return atr


def _write_geometry(out, dem, atr):
    """A MintPy geometry .h5 with `height`, written whole or not at all."""
    import h5py

    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    tmp = f'{out}.{os.getpid()}.part'
    with h5py.File(tmp, 'w') as fid:
        fid.create_dataset('height', data=dem, compression='gzip', compression_opts=4,
                           shuffle=True)
        for name, value in atr.items():
            fid.attrs[name] = value
    os.replace(tmp, out)
    return out


def _padded(geo_box, margin):
    """The map box (W, N, E, S) widened by `margin` of its size on every side."""
    west, north, east, south = (float(v) for v in geo_box)
    pad_x, pad_y = (east - west) * margin, (north - south) * margin
    return (west - pad_x, north + pad_y, east + pad_x, south - pad_y)


def terrain_dem(geo_box, data_crs, projected, zoom=None, max_tiles=MAX_TILES, folder=None,
                fetch=None, margin=0.02, verbose=True, max_side=MAX_SIDE):
    """Build (or reuse) a DEM covering a map extent from the Terrain Tiles, on the
    data's own CRS.

    Parameters: geo_box   - (W, N, E, S) of the map, in the data's coordinates
                data_crs  - the data's CRS: EPSG int or WKT
                projected - bool, the data's grid is in metres (else degrees)
                zoom      - int, tile zoom 0-15 (default: the finest within max_tiles)
                max_tiles - int, download budget for the automatic zoom
                folder    - str, cache folder (default: cache_dir())
                fetch     - callable(zoom, x, y) -> PNG bytes, to override the download
                margin    - float, extra border around the map, as a fraction of its size
                max_side  - int, at most this many pixels a side
    Returns:    str, path to a MintPy geometry .h5 with `height` in metres
    """
    from pyproj import Transformer

    folder = folder or cache_dir()
    box = _padded(geo_box, margin)
    bbox = _lonlat_box(box, data_crs)
    if zoom is None:
        zoom = choose_zoom(bbox, max_tiles=max_tiles)
    zoom = int(np.clip(zoom, 0, 15))
    lat_mid = 0.5 * (bbox[1] + bbox[3])
    res_m = float(tile_resolution(zoom, lat_mid))

    # the output grid: square pixels about the tile resolution, in the data's units
    step = res_m if projected else res_m / M_PER_DEG
    width, length, step = _grid_size(box, step, max_side)

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

    atr = _geo_attrs(box, step, length, width, projected, data_crs)
    atr.update({'SOURCE': TERRARIUM_URL, 'TILE_ZOOM': str(zoom), 'ATTRIBUTION': ATTRIBUTION})
    _write_geometry(out, dem, atr)
    if verbose:
        print(f'terrain: wrote {out} ({length} x {width}, {_range_text(dem)})')
    return out


def _range_text(dem):
    if not np.isfinite(dem).any():
        return 'no data'
    return f'{np.nanmin(dem):.0f} to {np.nanmax(dem):.0f} m'


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

    meta = _radar_attrs(lookup, atr)
    meta.update({'SOURCE': TERRARIUM_URL, 'TILE_ZOOM': str(zoom), 'ATTRIBUTION': ATTRIBUTION})
    _write_geometry(out, dem, meta)
    if verbose:
        print(f'terrain: wrote {out} ({lookup.length} x {lookup.width}, radar grid, '
              f'{_range_text(dem)})')
    return out


def _radar_attrs(lookup, atr):
    """The data's metadata, made a radar-coded geometry file on the same grid."""
    meta = {key: str(val) for key, val in atr.items()}
    meta.update({
        'FILE_TYPE': 'geometry', 'UNIT': 'm', 'DATA_TYPE': 'float32',
        'LENGTH': str(lookup.length), 'WIDTH': str(lookup.width),
        'PROCESSOR': meta.get('PROCESSOR', 'viewinsar'),
        'LOOKUP_FILE': os.path.abspath(lookup.file),
    })
    return meta


##############################################################################
#                       Copernicus DEM GLO-30                                #
##############################################################################
def cop30_name(lat, lon):
    """The GLO-30 tile of the 1 x 1 degree cell whose south-west corner is
    (lat, lon), in whole degrees: (-1, 119) -> Copernicus_DSM_COG_10_S01_00_E119_00_DEM."""
    lat, lon = int(lat), int(lon)
    return (f"Copernicus_DSM_COG_10_{'N' if lat >= 0 else 'S'}{abs(lat):02d}_00_"
            f"{'E' if lon >= 0 else 'W'}{abs(lon):03d}_00_DEM")


# TIFF field types: numpy code, bytes per value; and the sample formats read
_TIFF_TYPES = {1: ('u1', 1), 2: ('S1', 1), 3: ('u2', 2), 4: ('u4', 4), 5: ('u4', 8),
               6: ('i1', 1), 7: ('u1', 1), 8: ('i2', 2), 9: ('i4', 4), 10: ('i4', 8),
               11: ('f4', 4), 12: ('f8', 8), 13: ('u4', 4), 16: ('u8', 8), 17: ('i8', 8),
               18: ('u8', 8)}
_SAMPLE_TYPES = {(1, 8): 'u1', (1, 16): 'u2', (1, 32): 'u4', (2, 8): 'i1', (2, 16): 'i2',
                 (2, 32): 'i4', (3, 32): 'f4', (3, 64): 'f8'}


class _Cog:
    """Just enough of a GeoTIFF reader for the GLO-30 tiles - lon/lat, one band,
    tiled (or stripped), DEFLATE or uncompressed, predictor 1-3, overviews - that
    works from byte ranges: read(offset, size) -> bytes. No GDAL needed.

    levels[k] - dict: width, length, tw, th (internal tile), offsets, counts,
                dtype, comp, pred, sx, sy (pixel size, deg); 0 = full resolution,
                then ever coarser overviews
    x0, y0    - the top-left corner of the image, x1, y1 the bottom-right (deg)
    nodata    - the GDAL_NODATA value, or None
    """

    def __init__(self, read):
        self._read = read
        head = read(0, 16)
        if len(head) < 16 or head[:2] not in (b'II', b'MM'):
            raise ValueError('not a TIFF file')
        self.bo = '<' if head[:2] == b'II' else '>'
        magic = self._unpack('H', head[2:4])
        if magic == 42:
            self.big, ifd = False, self._unpack('I', head[4:8])
        elif magic == 43:                              # BigTIFF
            self.big, ifd = True, self._unpack('Q', head[8:16])
        else:
            raise ValueError('not a TIFF file')
        self.levels, self.nodata, seen = [], None, set()
        while ifd and ifd not in seen and len(seen) < 32:
            seen.add(ifd)
            tags, ifd = self._ifd(ifd)
            if int(self._first(tags, 254, 0)) & 4:     # a transparency mask
                continue
            level = self._level(tags)
            if not self.levels:
                self._georeference(tags, level)
            self.levels.append(level)
        if not self.levels:
            raise ValueError('a TIFF file without an image')
        full = self.levels[0]
        for level in self.levels:                      # overviews: the same extent
            level['sx'] = full['sx'] * full['width'] / level['width']
            level['sy'] = full['sy'] * full['length'] / level['length']
        self.levels.sort(key=lambda lv: lv['sy'])
        self.x1 = self.x0 + full['width'] * full['sx']
        self.y1 = self.y0 - full['length'] * full['sy']

    def _unpack(self, fmt, data):
        return struct.unpack(self.bo + fmt, data)[0]

    @staticmethod
    def _first(tags, tag, default):
        value = tags.get(tag)
        return default if value is None or len(value) == 0 else value[0]

    def _ifd(self, offset):
        """{tag: values} of the directory at `offset`, and the next one's offset."""
        num_fmt, num_size, entry, link = ('Q', 8, 20, 8) if self.big else ('H', 2, 12, 4)
        num = self._unpack(num_fmt, self._read(offset, num_size))
        if num > 4096:
            raise ValueError('a corrupt TIFF directory')
        raw = self._read(offset + num_size, num * entry + link)
        off_fmt = 'Q' if self.big else 'I'
        tags = {}
        for k in range(num):
            item = raw[k * entry:(k + 1) * entry]
            tag, typ = struct.unpack(self.bo + 'HH', item[:4])
            if self.big:
                count, field = self._unpack('Q', item[4:12]), item[12:20]
            else:
                count, field = self._unpack('I', item[4:8]), item[8:12]
            if typ not in _TIFF_TYPES:
                continue
            code, size = _TIFF_TYPES[typ]
            nbytes = count * size
            data = (field[:nbytes] if nbytes <= len(field)
                    else self._read(self._unpack(off_fmt, field), nbytes))
            if typ == 2:
                tags[tag] = data.split(b'\0', 1)[0].decode('latin-1')
                continue
            values = np.frombuffer(data, dtype=np.dtype(code).newbyteorder(self.bo),
                                   count=count * (2 if typ in (5, 10) else 1))
            if typ in (5, 10):                         # rationals
                values = values[0::2] / np.where(values[1::2] == 0, 1, values[1::2])
            tags[tag] = values
        return tags, self._unpack(off_fmt, raw[num * entry:num * entry + link])

    def _level(self, tags):
        width, length = int(self._first(tags, 256, 0)), int(self._first(tags, 257, 0))
        if 322 in tags and 324 in tags:
            tw, th = int(self._first(tags, 322, 0)), int(self._first(tags, 323, 0))
            offsets, counts = tags[324], tags.get(325)
        else:                                          # strips: tiles of the full width
            tw, th = width, min(int(self._first(tags, 278, length)), length)
            offsets, counts = tags.get(273), tags.get(279)
        if not (width and length and tw and th) or offsets is None or counts is None:
            raise ValueError('a TIFF image without its tile table')
        if int(self._first(tags, 277, 1)) != 1:
            raise ValueError('a TIFF image of more than one band')
        fmt, bits = int(self._first(tags, 339, 1)), int(self._first(tags, 258, 1))
        if (fmt, bits) not in _SAMPLE_TYPES:
            raise ValueError(f'TIFF samples of {bits} bits (format {fmt}) are not supported')
        comp, pred = int(self._first(tags, 259, 1)), int(self._first(tags, 317, 1))
        if comp not in (1, 8, 32946):
            raise ValueError(f'TIFF compression {comp} is not supported (DEFLATE is)')
        if pred not in (1, 2, 3):
            raise ValueError(f'TIFF predictor {pred} is not supported')
        return {'width': width, 'length': length, 'tw': tw, 'th': th,
                'offsets': np.asarray(offsets, dtype=np.int64),
                'counts': np.asarray(counts, dtype=np.int64),
                'dtype': np.dtype(_SAMPLE_TYPES[fmt, bits]).newbyteorder(self.bo),
                'comp': comp, 'pred': pred}

    def _georeference(self, tags, level):
        keys, geokeys = tags.get(34735), {}
        if keys is not None and len(keys) >= 4:
            for k in range(int(keys[3])):
                item = [int(v) for v in keys[4 + 4 * k:8 + 4 * k]]
                if len(item) == 4 and item[1] == 0:
                    geokeys[item[0]] = item[3]
        if geokeys.get(1024, 2) != 2:                  # GTModelType: 2 = geographic
            raise ValueError('a GeoTIFF that is not in lon/lat')
        scale, tie, matrix = tags.get(33550), tags.get(33922), tags.get(34264)
        if scale is not None and tie is not None and len(scale) >= 2 and len(tie) >= 6:
            sx, sy = float(scale[0]), float(scale[1])
            x0, y0 = float(tie[3]) - float(tie[0]) * sx, float(tie[4]) + float(tie[1]) * sy
        elif matrix is not None and len(matrix) >= 8 and matrix[1] == 0 and matrix[4] == 0:
            sx, sy = float(matrix[0]), -float(matrix[5])
            x0, y0 = float(matrix[3]), float(matrix[7])
        else:
            raise ValueError('a GeoTIFF without its geo-referencing')
        if not (sx > 0 and sy > 0):
            raise ValueError('a GeoTIFF on a flipped or rotated grid')
        if geokeys.get(1025) == 2:                     # PixelIsPoint: a pixel's centre
            x0, y0 = x0 - 0.5 * sx, y0 + 0.5 * sy
        self.x0, self.y0 = x0, y0
        level['sx'], level['sy'] = sx, sy
        text = tags.get(42113)                         # GDAL_NODATA
        if isinstance(text, str) and text.strip():
            try:
                self.nodata = float(text)
            except ValueError:
                self.nodata = None

    def level_for(self, res):
        """The coarsest level whose rows are at most sqrt(2) x `res` degrees apart:
        the one nearest a grid of that spacing, never much coarser."""
        best = 0
        for k, level in enumerate(self.levels):
            if level['sy'] <= res * np.sqrt(2.0) * (1 + 1e-6):
                best = k
        return best

    def box(self, k, bbox, apron=2):
        """(r0, r1, c0, c1) of level k over the lon/lat box (W, S, E, N), `apron`
        pixels wider on each side; None where they do not overlap."""
        level = self.levels[k]
        west, south = max(bbox[0], self.x0), max(bbox[1], self.y1)
        east, north = min(bbox[2], self.x1), min(bbox[3], self.y0)
        if west >= east or south >= north:
            return None
        c0 = max(int(np.floor((west - self.x0) / level['sx'])) - apron, 0)
        c1 = min(int(np.ceil((east - self.x0) / level['sx'])) + apron, level['width'])
        r0 = max(int(np.floor((self.y0 - north) / level['sy'])) - apron, 0)
        r1 = min(int(np.ceil((self.y0 - south) / level['sy'])) + apron, level['length'])
        return (r0, r1, c0, c1) if r0 < r1 and c0 < c1 else None

    def chunk_ids(self, k, r0, r1, c0, c1):
        """Indices of the internal tiles of level k holding rows r0:r1, cols c0:c1."""
        level = self.levels[k]
        across = -(-level['width'] // level['tw'])
        return [ty * across + tx
                for ty in range(r0 // level['th'], (r1 - 1) // level['th'] + 1)
                for tx in range(c0 // level['tw'], (c1 - 1) // level['tw'] + 1)]

    def decode(self, k, data):
        """One internal tile (or strip) of level k: float32, th x tw, NaN for nodata."""
        level = self.levels[k]
        th, tw, dtype = level['th'], level['tw'], level['dtype']
        raw = zlib.decompress(data) if level['comp'] != 1 else bytes(data)
        rows = min(th, len(raw) // (tw * dtype.itemsize))
        if rows < 1:
            raise ValueError('a truncated TIFF tile')
        raw = raw[:rows * tw * dtype.itemsize]
        if level['pred'] == 3:
            # floating point (TIFF Technical Note 3): each row's bytes differenced,
            # after being shuffled most significant first
            byt = np.frombuffer(raw, np.uint8).reshape(rows, tw * dtype.itemsize)
            byt = np.cumsum(byt, axis=1, dtype=np.uint8)
            byt = np.ascontiguousarray(byt.reshape(rows, dtype.itemsize, tw).transpose(0, 2, 1))
            values = byt.view(dtype.newbyteorder('>')).reshape(rows, tw)
        else:
            values = np.frombuffer(raw, dtype).reshape(rows, tw)
            if level['pred'] == 2:                     # horizontal differencing
                values = np.cumsum(values, axis=1, dtype=dtype)
        out = np.full((th, tw), np.nan, dtype=np.float32)
        out[:rows] = values
        if self.nodata is not None and np.isfinite(self.nodata):
            out[out == np.float32(self.nodata)] = np.nan
        return out

    def window(self, k, r0, r1, c0, c1, chunk):
        """Rows r0:r1, columns c0:c1 of level k (float32); chunk(k, index) gives the
        bytes of an internal tile (empty: a sparse one, left NaN)."""
        level = self.levels[k]
        th, tw = level['th'], level['tw']
        across = -(-level['width'] // tw)
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        for idx in self.chunk_ids(k, r0, r1, c0, c1):
            data = chunk(k, idx)
            if not data:
                continue
            ty, tx = divmod(idx, across)
            tile = self.decode(k, data)
            ra, rb = max(r0, ty * th), min(r1, ty * th + th)
            ca, cb = max(c0, tx * tw), min(c1, tx * tw + tw)
            out[ra - r0:rb - r0, ca - c0:cb - c0] = tile[ra - ty * th:rb - ty * th,
                                                         ca - tx * tw:cb - tx * tw]
        return out


def _save(path, data):
    """Write a cache file whole or not at all."""
    import threading

    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f'{path}.{os.getpid()}.{threading.get_ident()}.part'
    with open(tmp, 'wb') as fid:
        fid.write(data)
    os.replace(tmp, path)


class _Cop30Tile:
    """One GLO-30 tile: its headers and the internal tiles read so far, kept in
    cop30/<name>/ of the cache folder - or `missing` there: no such tile (sea)."""

    def __init__(self, lat, lon, folder, fetch):
        self.lat, self.lon = int(lat), int(lon)
        self.name = cop30_name(lat, lon)
        self.url = COP30_URL.format(name=self.name)
        self.dir = os.path.join(folder, 'cop30', self.name)
        self.fetch = fetch
        self.cog = self.win = self.error = None
        self.level, self.missing, self._head = 0, False, b''

    def open(self):
        """Read the tile's headers (cached); False if there is no such tile."""
        if os.path.isfile(os.path.join(self.dir, 'missing')):
            self.missing = True
            return False
        path = os.path.join(self.dir, 'head.bin')
        cached = os.path.isfile(path) and os.path.getsize(path) > 0
        if cached:
            with open(path, 'rb') as fid:
                self._head = fid.read()
        else:
            self._head = self.fetch(self.url, 0, HEAD_BYTES)
            if self._head is None:
                self.missing = True
                _save(os.path.join(self.dir, 'missing'), b'')
                return False
        size = len(self._head)
        self.cog = _Cog(self._read)
        if not cached or len(self._head) != size:
            _save(path, self._head)
        return True

    def _read(self, offset, size):
        if offset + size > len(self._head):            # headers past the first bytes
            more = self.fetch(self.url, 0, max(2 * len(self._head), offset + size + 4096))
            if not more or len(more) < offset + size:
                raise ValueError(f'{self.name}: a truncated TIFF file')
            self._head = more
        return self._head[offset:offset + size]

    def chunk_path(self, k, idx):
        return os.path.join(self.dir, f'L{k}_{idx}.bin')

    def chunk(self, k, idx):
        """The bytes of internal tile `idx` of level k, from the cache or the web."""
        path = self.chunk_path(k, idx)
        if os.path.isfile(path):
            with open(path, 'rb') as fid:
                return fid.read()
        level = self.cog.levels[k]
        offset, count = int(level['offsets'][idx]), int(level['counts'][idx])
        if count <= 0:
            return b''
        data = self.fetch(self.url, offset, count)
        if data is None or len(data) != count:
            raise ValueError(f'{self.name}: block {idx} came incomplete')
        _save(path, data)
        return data


class _MapGrid:
    """The pixel centres of a DEM grid on the data's CRS - a box (W, N, E, S) of
    square `step` pixels - and the part of it a lon/lat box falls on."""

    def __init__(self, box, step, width, length, crs):
        from pyproj import CRS, Transformer

        self.x0, self.y0, self.step = float(box[0]), float(box[1]), float(step)
        self.shape = (int(length), int(width))
        if CRS.from_user_input(crs).equals(CRS.from_epsg(4326), ignore_axis_order=True):
            self.to_ll = self.from_ll = None
        else:
            self.to_ll = Transformer.from_crs(crs, 4326, always_xy=True)
            self.from_ll = Transformer.from_crs(4326, crs, always_xy=True)

    def valid(self):
        return None                                    # every pixel has a lon/lat

    def blocks(self, bbox, rows=256):
        """(row slice, column slice, lon, lat) of the pixels that may lie in the
        lon/lat box (W, S, E, N), a few rows at a time."""
        xs, ys = _outline(bbox[0], bbox[3], bbox[2], bbox[1])
        if self.from_ll is not None:
            xs, ys = self.from_ll.transform(xs, ys)
            ok = np.isfinite(xs) & np.isfinite(ys)
            if not ok.any():
                return
            xs, ys = xs[ok], ys[ok]
        c0 = max(int(np.floor((np.min(xs) - self.x0) / self.step)) - 1, 0)
        c1 = min(int(np.ceil((np.max(xs) - self.x0) / self.step)) + 1, self.shape[1])
        r0 = max(int(np.floor((self.y0 - np.max(ys)) / self.step)) - 1, 0)
        r1 = min(int(np.ceil((self.y0 - np.min(ys)) / self.step)) + 1, self.shape[0])
        if c0 >= c1 or r0 >= r1:
            return
        gx = self.x0 + (np.arange(c0, c1) + 0.5) * self.step
        for ra in range(r0, r1, rows):
            rb = min(ra + rows, r1)
            lon, lat = np.meshgrid(gx, self.y0 - (np.arange(ra, rb) + 0.5) * self.step)
            if self.to_ll is not None:
                lon, lat = self.to_ll.transform(lon, lat)
            yield slice(ra, rb), slice(c0, c1), np.asarray(lon), np.asarray(lat)


class _LookupGrid:
    """The pixels of a radar grid, by the lon/lat of each (NaN where unknown)."""

    def __init__(self, lon, lat, rows=256):
        self.lon, self.lat, self.rows = lon, lat, rows
        self.shape = lon.shape
        self.boxes = []
        for ra in range(0, self.shape[0], rows):
            lo, la = lon[ra:ra + rows], lat[ra:ra + rows]
            ok = np.isfinite(lo) & np.isfinite(la)
            self.boxes.append((lo[ok].min(), la[ok].min(), lo[ok].max(), la[ok].max())
                              if ok.any() else None)

    def valid(self):
        return np.isfinite(self.lon) & np.isfinite(self.lat)

    def blocks(self, bbox):
        for i, box in enumerate(self.boxes):
            if (box is None or box[0] > bbox[2] or box[2] < bbox[0] or box[1] > bbox[3]
                    or box[3] < bbox[1]):
                continue
            ra, rb = i * self.rows, min((i + 1) * self.rows, self.shape[0])
            yield slice(ra, rb), slice(None), self.lon[ra:rb], self.lat[ra:rb]


def _grid_spacing(lon, lat):
    """The finer of the row and column spacings of a lon/lat grid, in degrees of
    latitude: a median over a sample of its pixels."""
    nrow, ncol = lon.shape
    if nrow < 2 or ncol < 2:
        return COP30_RES
    rr, cc = np.meshgrid(np.arange(0, nrow - 1, max(1, nrow // 200)),
                         np.arange(0, ncol - 1, max(1, ncol // 200)), indexing='ij')
    coslat = np.cos(np.radians(lat[rr, cc]))
    found = []
    for dr, dc in ((1, 0), (0, 1)):
        dist = np.hypot(lat[rr + dr, cc + dc] - lat[rr, cc],
                        (lon[rr + dr, cc + dc] - lon[rr, cc]) * coslat)
        if np.isfinite(dist).any():
            found.append(float(np.nanmedian(dist)))
    return min(found) if found and min(found) > 0 else COP30_RES


def _cop30_cells(bbox, pad=0.01):
    """(lat, lon) of the south-west corners of the GLO-30 cells over a lon/lat box
    (W, S, E, N), and a little around it, north to south."""
    lat0, lat1 = max(int(np.floor(bbox[1] - pad)), -90), min(int(np.floor(bbox[3] + pad)), 89)
    lon0, lon1 = max(int(np.floor(bbox[0] - pad)), -180), min(int(np.floor(bbox[2] + pad)), 179)
    return [(lat, lon) for lat in range(lat1, lat0 - 1, -1) for lon in range(lon0, lon1 + 1)]


def _arcsec(res):
    return f'{res * 3600:.3g}" (~{res * M_PER_DEG:.0f} m)'


def _usable(tile):
    return tile.cog is not None and tile.error is None and tile.win is not None


def _cop30_strip(ref, near, bbox):
    """One row of GLO-30 tiles - the cells of ref's latitude - as a single array
    on ref's pixel lattice, two pixels wider all round, those taken from the
    neighbouring tiles: so the heights interpolate across every seam. 0 m in the
    cells without a tile (the sea), NaN where a tile failed.

    ref  - a usable tile of the row (opened, with a window over the map)
    near - {(lat, lon): tile} of every tile looked for
    Returns (mosaic, ra, ca): mosaic[i, j] is pixel (ra + i, ca + j) of ref's level.
    """
    cog, level = ref.cog, ref.cog.levels[ref.level]
    sx, sy = level['sx'], level['sy']
    per_x, per_y = int(round(1.0 / sx)), int(round(1.0 / sy))       # pixels a degree
    north, south = min(bbox[3], cog.y0), max(bbox[1], cog.y1)
    ra = int(np.floor((cog.y0 - north) / sy)) - 2
    rb = int(np.ceil((cog.y0 - south) / sy)) + 2
    ca = int(np.floor((bbox[0] - cog.x0) / sx)) - 2
    cb = int(np.ceil((bbox[2] - cog.x0) / sx)) + 2
    mosaic = np.full((rb - ra, cb - ca), np.nan, dtype=np.float32)
    for (lat, lon), tile in near.items():
        if abs(lat - ref.lat) > 1:
            continue
        dr, dc = (ref.lat - lat) * per_y, (lon - ref.lon) * per_x     # its pixel 0, 0
        if tile.missing:
            i0, j0 = max(dr - ra, 0), max(dc - ca, 0)
            mosaic[i0:max(dr + per_y - ra, i0), j0:max(dc + per_x - ca, j0)] = 0.0
            continue
        if not _usable(tile):
            continue
        tl = tile.cog.levels[tile.level]
        if (not np.isclose(tl['sx'], sx, rtol=1e-6) or not np.isclose(tl['sy'], sy, rtol=1e-6)
                or abs((tile.cog.x0 - cog.x0) / sx - dc) > 1e-3
                or abs((cog.y0 - tile.cog.y0) / sy - dr) > 1e-3):
            continue                                   # another lattice (a band edge)
        i0, i1 = max(ra - dr, tile.win[0]), min(rb - dr, tile.win[1])
        j0, j1 = max(ca - dc, tile.win[2]), min(cb - dc, tile.win[3])
        if i0 >= i1 or j0 >= j1:
            continue
        try:
            part = tile.cog.window(tile.level, i0, i1, j0, j1, tile.chunk)
        except Exception as err:   # noqa: BLE001 - reported by the caller
            tile.error = err
            continue
        mosaic[i0 + dr - ra:i1 + dr - ra, j0 + dc - ca:j1 + dc - ca] = part
    # rows from beyond the row that could not be had: the row's own edge
    own0, own1 = max(-ra, 0), min(level['length'] - ra, rb - ra)
    for i in [*range(0, own0), *range(own1, rb - ra)]:
        edge = mosaic[own0] if i < own0 else mosaic[own1 - 1]
        gap = np.isnan(mosaic[i])
        mosaic[i, gap] = edge[gap]
    return mosaic, ra, ca


def _cop30_fill(grid, res, bbox, folder, fetch, verbose, workers=8):
    """The GLO-30 heights on `grid` (a _MapGrid or _LookupGrid), float32: 0 m at sea
    (where there is no tile), NaN where a tile could not be had. Also the spacing
    of the levels read, as text."""
    cells = _cop30_cells(bbox)
    if len(cells) > COP30_MAX_TILES:
        raise ValueError(f'the map spans {len(cells)} GLO-30 tiles of 1 x 1 deg, more than '
                         f'{COP30_MAX_TILES}: crop it (--sub-lat / --sub-lon), or use '
                         '--dem-source terrain')
    tiles = [_Cop30Tile(lat, lon, folder, fetch) for lat, lon in cells]

    def opened(tile):
        try:
            tile.open()
        except Exception as err:   # noqa: BLE001 - reported below
            tile.error = err

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(opened, tiles))

    # the blocks each tile needs, from the overview nearest the grid's spacing
    jobs = []
    for tile in tiles:
        if tile.cog is None:
            continue
        tile.level = tile.cog.level_for(res)
        tile.win = tile.cog.box(tile.level, bbox)
        if tile.win is not None:
            jobs += [(tile, tile.level, idx)
                     for idx in tile.cog.chunk_ids(tile.level, *tile.win)
                     if not os.path.isfile(tile.chunk_path(tile.level, idx))]
    if jobs:
        if verbose:
            nbytes = sum(int(t.cog.levels[k]['counts'][i]) for t, k, i in jobs)
            print(f'terrain: downloading {len(jobs)} block(s) of GLO-30, '
                  f'{nbytes / 1e6:.1f} MB ...')

        def download(job):
            tile, k, idx = job
            try:
                tile.chunk(k, idx)
            except Exception as err:   # noqa: BLE001 - reported below
                tile.error = tile.error or err

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            for done, _ in enumerate(pool.map(download, jobs), start=1):
                if verbose and len(jobs) >= 20 and done % max(1, len(jobs) // 5) == 0:
                    print(f'   ... {done}/{len(jobs)} blocks')

    # row by row of tiles: one array for the row, sampled at every pixel in it
    dem = np.full(grid.shape, np.nan, dtype=np.float32)
    covered = np.zeros(grid.shape, dtype=bool)
    near = {(t.lat, t.lon): t for t in tiles}
    spacings = set()
    for lat in sorted({t.lat for t in tiles}, reverse=True):
        row = [t for t in tiles if t.lat == lat and _usable(t)]
        if not row:
            continue
        ref = row[0]
        mosaic, ra, ca = _cop30_strip(ref, near, bbox)
        cog, level = ref.cog, ref.cog.levels[ref.level]
        spacings.add(level['sy'])
        for rs, cs, lon, lat_ in grid.blocks((bbox[0], cog.y1, bbox[2], cog.y0)):
            inside = (lat_ > cog.y1) & (lat_ <= cog.y0) & np.isfinite(lon)
            if not inside.any():
                continue
            rows = (cog.y0 - lat_[inside]) / level['sy'] - 0.5 - ra
            cols = (lon[inside] - cog.x0) / level['sx'] - 0.5 - ca
            dem[rs, cs][inside] = _bilinear(mosaic, rows, cols)
            covered[rs, cs][inside] = True
        del mosaic

    failed = [t for t in tiles if t.error is not None]
    missing = [t for t in tiles if t.missing]
    if not any(_usable(t) for t in tiles) and not missing:
        err = failed[0].error if failed else 'nothing to read'
        raise RuntimeError(f'could not download any Copernicus GLO-30 tile ({err}). Is there '
                           'internet access (or HTTPS_PROXY)? A DEM file works too: --dem dem.tif')

    # rows without a tile at all: the sea, 0 m - except in a tile that failed
    sea = ~covered
    del covered
    valid = grid.valid()
    if valid is not None:
        sea &= valid
    for tile in failed:
        pad = 2 * COP30_RES
        cell = (tile.lon - pad, tile.lat - pad, tile.lon + 1 + pad, tile.lat + 1 + pad)
        for rs, cs, lon, lat_ in grid.blocks(cell):
            sea[rs, cs] &= ~((lon >= cell[0]) & (lon <= cell[2])
                             & (lat_ >= cell[1]) & (lat_ <= cell[3]))
    dem[sea] = 0.0
    del sea
    if verbose:
        if missing:
            print(f'terrain: {len(missing)} of {len(tiles)} GLO-30 cells have no tile - sea '
                  '(or a country not released) - 0 m there')
        if failed:
            print(f'   WARNING: {len(failed)} of {len(tiles)} GLO-30 tiles failed '
                  f'({failed[0].error}); those areas are left empty.')
    return dem, ', '.join(_arcsec(s) for s in sorted(spacings)) or 'none'


def copernicus_dem(geo_box, data_crs, projected, max_side=MAX_SIDE, folder=None, fetch=None,
                   margin=0.02, verbose=True):
    """Build (or reuse) a DEM covering a map extent from the Copernicus GLO-30
    tiles, on the data's own CRS: square pixels of 1" (30.9 m on a projected
    grid), wider where a side would pass max_side pixels - then read from the
    tiles' overview nearest that spacing.

    Parameters: geo_box   - (W, N, E, S) of the map, in the data's coordinates
                data_crs  - the data's CRS: EPSG int or WKT
                projected - bool, the data's grid is in metres (else degrees)
                max_side  - int, at most this many pixels a side
                folder    - str, cache folder (default: cache_dir())
                fetch     - callable(url, start, size) -> bytes or None, the download
                margin    - float, extra border around the map, as a fraction of its size
    Returns:    str, path to a MintPy geometry .h5 with `height` in metres above the
                EGM2008 geoid (0 at sea)
    """
    folder = folder or cache_dir()
    box = _padded(geo_box, margin)
    bbox = _lonlat_box(box, data_crs, limit=90.0)
    step = COP30_RES * (M_PER_DEG if projected else 1.0)
    width, length, step = _grid_size(box, step, max_side)
    res = step / M_PER_DEG if projected else step
    key = hashlib.sha1(repr(('cop30', str(data_crs), [round(v, 6) for v in box],
                             round(step, 12))).encode()).hexdigest()[:10]
    out = os.path.join(folder, f'dem_cop30_{key}.h5')
    if verbose:
        print(f'terrain: Copernicus GLO-30, {len(_cop30_cells(bbox))} tile(s) of 1 x 1 deg over '
              f'lon {bbox[0]:.3f}..{bbox[2]:.3f}, lat {bbox[1]:.3f}..{bbox[3]:.3f}: '
              f'{length} x {width} pixels of {_arcsec(res)}')
    if os.path.isfile(out):
        if verbose:
            print(f'terrain: re-using {out}')
        return out

    grid = _MapGrid(box, step, width, length, data_crs)
    dem, read = _cop30_fill(grid, res, bbox, folder, fetch or fetch_range, verbose)
    atr = _geo_attrs(box, step, length, width, projected, data_crs)
    atr.update({'SOURCE': f'Copernicus DEM GLO-30, {COP30_URL}', 'COP30_READ_AT': read,
                'VERTICAL_DATUM': 'EGM2008 geoid', 'ATTRIBUTION': COP30_ATTRIBUTION})
    _write_geometry(out, dem, atr)
    if verbose:
        print(f'terrain: wrote {out} ({length} x {width}, {_range_text(dem)})')
    return out


def copernicus_dem_radar(lookup, atr, folder=None, fetch=None, verbose=True):
    """Copernicus GLO-30 on a RADAR grid, through the lookup's lat/lon of every
    pixel - read at the overview nearest the grid's own spacing.

    For radar-coded data whose geometry file has no `height` of its own.

    Parameters: lookup - radar.RadarLookup of the data
                atr    - dict, the data's metadata (copied into the DEM file, so
                         view.py and the swath read it on the same grid)
    Returns:    str, path to a radar-coded MintPy geometry .h5 with `height`
    """
    folder = folder or cache_dir()
    bbox = lookup.extent(margin=0.01)
    res = _grid_spacing(lookup.lon, lookup.lat)
    stat = os.stat(lookup.file)
    key = hashlib.sha1(repr(('cop30', os.path.abspath(lookup.file), stat.st_size,
                             int(stat.st_mtime))).encode()).hexdigest()[:10]
    out = os.path.join(folder, f'dem_cop30_radar_{key}.h5')
    if verbose:
        print(f'terrain: Copernicus GLO-30, {len(_cop30_cells(bbox))} tile(s) of 1 x 1 deg over '
              f'the radar frame (pixels ~{_arcsec(res)} apart), resampled onto it through '
              'the lookup')
    if os.path.isfile(out):
        if verbose:
            print(f'terrain: re-using {out}')
        return out

    grid = _LookupGrid(lookup.lon, lookup.lat)
    dem, read = _cop30_fill(grid, res, bbox, folder, fetch or fetch_range, verbose)
    meta = _radar_attrs(lookup, atr)
    meta.update({'SOURCE': f'Copernicus DEM GLO-30, {COP30_URL}', 'COP30_READ_AT': read,
                 'VERTICAL_DATUM': 'EGM2008 geoid', 'ATTRIBUTION': COP30_ATTRIBUTION})
    _write_geometry(out, dem, meta)
    if verbose:
        print(f'terrain: wrote {out} ({lookup.length} x {lookup.width}, radar grid, '
              f'{_range_text(dem)})')
    return out

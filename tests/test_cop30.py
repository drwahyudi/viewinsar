"""--dem / --topoxsect without a file: the Copernicus GLO-30 DEM.

No network here: GLO-30-like Cloud Optimized GeoTIFFs - float32, tiled,
DEFLATE, the floating-point predictor, half and quarter resolution overviews,
PixelIsPoint at the cell's whole-degree corner, as the real ones are - are
written from a known height field (by numpy alone, or by GDAL where rasterio is
installed) and served by byte range in place of the download. The test tiles
are 360 pixels of 10" a side (the real ones 3600 of 1"). A steep plane on top
of the ridge would show any seam between tiles, a north-south gradient any
flipped axis; a cell without a tile is the sea.
"""

import os
import struct
import zlib

import numpy as np
import pytest

from conftest import FAULT_LON, REAL_GET
from viewinsar import cli, terrain

pyproj = pytest.importorskip('pyproj')
h5py = pytest.importorskip('h5py')

RES = 1.0 / 360
LINE = ['--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75']
CELLS = [(29, 129), (29, 130), (30, 129), (30, 130)]
BOX = (129.6, 30.4, 130.4, 29.6)                       # W N E S: four cells meet in it


def height(lon, lat):
    """A ridge on the fixtures' fault, and a north-up gradient (100 m per degree)."""
    km_east = (lon - FAULT_LON) * 111.32 * np.cos(np.radians(30.3))
    return 1500.0 * np.exp(-(km_east / 7.0) ** 2) + 40.0 + 100.0 * (lat - 30.3)


def steep(lon, lat):
    """The ridge on a steep plane: 10 m and 6 m a pixel. Bilinear interpolation
    gives a plane back exactly - unless a seam between tiles breaks it."""
    return height(lon, lat) + 3600.0 * (lon - 130.0) + 2160.0 * (lat - 30.0)


##############################################################################
#                 a GLO-30-like GeoTIFF, with numpy alone                    #
##############################################################################
def _ifd(entries, at, next_ifd):
    """A little-endian TIFF directory placed at byte `at`, long values after it."""
    entries = sorted(entries)
    extra_at = at + 2 + 12 * len(entries) + 4
    head, extra = [struct.pack('<H', len(entries))], b''
    for tag, typ, values in entries:
        if typ == 2:
            data, count = values.encode() + b'\0', len(values) + 1
        else:
            fmt = {3: 'H', 4: 'I', 12: 'd'}[typ]
            data, count = struct.pack('<' + fmt * len(values), *values), len(values)
        if len(data) <= 4:
            field = data.ljust(4, b'\0')
        else:
            extra += b'\0' * (len(extra) % 2)
            field = struct.pack('<I', extra_at + len(extra))
            extra += data
        head.append(struct.pack('<HHI', tag, typ, count) + field)
    head.append(struct.pack('<I', next_ifd))
    return b''.join(head) + extra


def _encode(tile):
    """DEFLATE after the floating-point predictor: each row's bytes, most
    significant first, differenced."""
    rows, cols = tile.shape
    byt = (tile.astype('>f4').view(np.uint8).reshape(rows, cols, 4).transpose(0, 2, 1)
           .reshape(rows, 4 * cols))
    diff = byt.copy()
    diff[:, 1:] = byt[:, 1:] - byt[:, :-1]
    return zlib.compress(diff.tobytes(), 6)


def write_glo30(path, lat, lon, field=height, size=360, block=128, overviews=2, nodata=None):
    """A GLO-30-like tile of the cell (lat, lon): size x size float32 pixels,
    centred on whole multiples of 1/size degree from (lon, lat + 1) - the tie
    point, PixelIsPoint - in block x block DEFLATE tiles, the floating-point
    predictor, and `overviews` levels of 2 x 2 means. Returns the levels."""
    res = 1.0 / size
    rows, cols = np.mgrid[0:size, 0:size]
    pyramid = [field(lon + cols * res, lat + 1 - rows * res).astype(np.float32)]
    for _ in range(overviews):
        z = pyramid[-1]
        pyramid.append(0.25 * (z[0::2, 0::2] + z[1::2, 0::2] + z[0::2, 1::2] + z[1::2, 1::2]))
    blobs = []
    for z in pyramid:
        n_y, n_x = -(-z.shape[0] // block), -(-z.shape[1] // block)
        pad = np.zeros((n_y * block, n_x * block), dtype=np.float32)
        pad[:z.shape[0], :z.shape[1]] = z
        blobs.append([_encode(pad[i * block:(i + 1) * block, j * block:(j + 1) * block])
                      for i in range(n_y) for j in range(n_x)])

    def entries(k, offsets):
        z = pyramid[k]
        out = [(254, 4, [0 if k == 0 else 1]), (256, 4, [z.shape[1]]), (257, 4, [z.shape[0]]),
               (258, 3, [32]), (259, 3, [8]), (262, 3, [1]), (277, 3, [1]), (284, 3, [1]),
               (317, 3, [3]), (322, 3, [block]), (323, 3, [block]), (324, 4, offsets),
               (325, 4, [len(b) for b in blobs[k]]), (339, 3, [3])]
        if k == 0:
            out += [(33550, 12, [res, res, 0.0]), (33922, 12, [0, 0, 0, lon, lat + 1, 0]),
                    (34735, 3, [1, 1, 0, 3, 1024, 0, 1, 2, 1025, 0, 1, 2, 2048, 0, 1, 4326])]
            if nodata is not None:
                out.append((42113, 2, str(nodata)))
        return out

    sizes = [len(_ifd(entries(k, [0] * len(blobs[k])), 0, 0)) for k in range(len(pyramid))]
    ifd_at, pos = [], 8
    for size_k in sizes:
        ifd_at.append(pos)
        pos += size_k + size_k % 2
    offsets = []
    for blob in blobs:
        offsets.append([])
        for data in blob:
            offsets[-1].append(pos)
            pos += len(data)
    with open(path, 'wb') as fid:
        fid.write(b'II' + struct.pack('<HI', 42, 8))
        for k in range(len(pyramid)):
            ifd = _ifd(entries(k, offsets[k]), ifd_at[k],
                       ifd_at[k + 1] if k + 1 < len(pyramid) else 0)
            fid.write(ifd + b'\0' * (len(ifd) % 2))
        for blob in blobs:
            for data in blob:
                fid.write(data)
    return pyramid


def serve_glo30(monkeypatch, folder, cells, field=height):
    """Write GLO-30-like tiles of `cells` into `folder` and serve them by byte
    range instead of the internet; the list of requests is returned."""
    os.makedirs(folder, exist_ok=True)
    for lat, lon in cells:
        write_glo30(os.path.join(folder, f'{terrain.cop30_name(lat, lon)}.tif'), lat, lon,
                    field=field)
    requests = []

    def fetch(url, start, size):
        name = url.rsplit('/', 1)[1]
        requests.append((name, start, size))
        path = os.path.join(folder, name)
        if not os.path.isfile(path):
            return None                                # no tile: the sea
        with open(path, 'rb') as fid:
            fid.seek(start)
            return fid.read(size)

    monkeypatch.setattr(terrain, 'fetch_range', fetch)
    monkeypatch.setattr(terrain, 'COP30_RES', RES)
    return requests


def _read(path):
    with h5py.File(path, 'r') as fid:
        return fid['height'][:], {k: (v.decode() if isinstance(v, bytes) else v)
                                  for k, v in fid.attrs.items()}


def _centres(atr):
    xs = float(atr['X_FIRST']) + (np.arange(int(atr['WIDTH'])) + 0.5) * float(atr['X_STEP'])
    ys = float(atr['Y_FIRST']) + (np.arange(int(atr['LENGTH'])) + 0.5) * float(atr['Y_STEP'])
    return np.meshgrid(xs, ys)


def _whole(cog, blob, k):
    level = cog.levels[k]

    def chunk(kk, idx):
        start = int(cog.levels[kk]['offsets'][idx])
        return blob[start:start + int(cog.levels[kk]['counts'][idx])]

    return cog.window(k, 0, level['length'], 0, level['width'], chunk)


##############################################################################
#                          names, the reader                                 #
##############################################################################
def test_tile_names_and_cells():
    assert terrain.cop30_name(-1, 119) == 'Copernicus_DSM_COG_10_S01_00_E119_00_DEM'
    assert terrain.cop30_name(0, -1) == 'Copernicus_DSM_COG_10_N00_00_W001_00_DEM'
    assert terrain.cop30_name(45, 6) == 'Copernicus_DSM_COG_10_N45_00_E006_00_DEM'
    assert terrain.COP30_URL.format(name='T') == (
        'https://copernicus-dem-30m.s3.amazonaws.com/T/T.tif')
    assert terrain._cop30_cells((119.2, -1.5, 120.5, -0.2)) == [
        (-1, 119), (-1, 120), (-2, 119), (-2, 120)]


def test_reader_reads_back_what_was_written(tmp_path):
    path = str(tmp_path / 't.tif')
    pyramid = write_glo30(path, 30, 130, nodata=-32767)
    blob = open(path, 'rb').read()
    cog = terrain._Cog(lambda offset, size: blob[offset:offset + size])
    assert [level['width'] for level in cog.levels] == [360, 180, 90]
    # PixelIsPoint at (130, 31): the image's corner is half a pixel up and left
    assert (cog.x0, cog.y0) == pytest.approx((130 - RES / 2, 31 + RES / 2), abs=1e-12)
    assert (cog.x1, cog.y1) == pytest.approx((131 - RES / 2, 30 + RES / 2), abs=1e-12)
    assert cog.nodata == -32767
    for k, z in enumerate(pyramid):
        assert np.array_equal(_whole(cog, blob, k), z)
    # the coarsest level no more than sqrt(2) x coarser than asked for
    assert [cog.level_for(r * RES) for r in (1, 1.4, 1.5, 2.1, 3, 9)] == [0, 0, 1, 1, 2, 2]


def test_reader_matches_gdal(tmp_path):
    rasterio = pytest.importorskip('rasterio')
    from rasterio.transform import from_origin

    path, size = str(tmp_path / 'gdal.tif'), 720
    res = 1.0 / size
    rows, cols = np.mgrid[0:size, 0:size]
    data = steep(130 + cols * res, 31 - rows * res).astype(np.float32)
    with rasterio.open(path, 'w', driver='COG', width=size, height=size, count=1,
                       dtype='float32', crs='EPSG:4326',
                       transform=from_origin(130 - res / 2, 31 + res / 2, res, res),
                       BLOCKSIZE=256, COMPRESS='DEFLATE', PREDICTOR='YES',
                       OVERVIEW_RESAMPLING='AVERAGE') as dst:
        dst.update_tags(AREA_OR_POINT='Point')
        dst.write(data, 1)
    blob = open(path, 'rb').read()
    cog = terrain._Cog(lambda offset, size: blob[offset:offset + size])
    with rasterio.open(path) as src:
        assert len(cog.levels) == len(src.overviews(1)) + 1 > 1
        assert (cog.x0, cog.y0) == pytest.approx((src.transform.c, src.transform.f), abs=1e-12)
        assert np.array_equal(_whole(cog, blob, 0), src.read(1))
    for k in range(1, len(cog.levels)):
        with rasterio.open(path, overview_level=k - 1) as ovr:
            assert np.array_equal(_whole(cog, blob, k), ovr.read(1))


def test_range_requests(monkeypatch):
    """Bytes first-last are asked for; a server that sends the whole file is
    coped with; a tile that does not exist (404, or S3's own 403) is told from
    a refusal (anyone else's 403)."""
    import io
    import urllib.error
    import urllib.request

    body, asked = bytes(range(256)) * 4, []

    class Reply(io.BytesIO):
        def __init__(self, data, status):
            super().__init__(data)
            self.status = status

    def urlopen(req, timeout=30):
        asked.append(req.get_header('Range'))
        url = req.full_url
        if url.endswith('whole'):
            return Reply(body, 200)
        if url.endswith('part'):
            return Reply(body[10:20], 206)
        if url.endswith('missing'):
            raise urllib.error.HTTPError(url, 404, 'Not Found', {}, None)
        server = 'AmazonS3' if url.endswith('s3') else 'a proxy'
        raise urllib.error.HTTPError(url, 403, 'Forbidden', {'Server': server}, None)

    monkeypatch.setattr(urllib.request, 'urlopen', urlopen)
    assert REAL_GET('https://x/part', byte_range=(10, 19)) == body[10:20]
    assert asked[-1] == 'bytes=10-19'
    assert REAL_GET('https://x/whole', byte_range=(10, 19)) == body[10:20]
    assert REAL_GET('https://x/missing', byte_range=(0, 9)) is None
    assert REAL_GET('https://x/s3', byte_range=(0, 9)) is None
    with pytest.raises(urllib.error.HTTPError):
        REAL_GET('https://x/refused', byte_range=(0, 9))


##############################################################################
#                           the DEM it makes                                 #
##############################################################################
def test_lonlat_dem_matches_the_field_across_the_seams(monkeypatch, tmp_path):
    serve_glo30(monkeypatch, tmp_path / 'server', CELLS[:3], field=steep)   # (30, 130): sea
    path = terrain.copernicus_dem(BOX, 4326, False, verbose=False)
    dem, atr = _read(path)
    lon, lat = _centres(atr)
    assert atr['ATTRIBUTION'] == terrain.COP30_ATTRIBUTION and atr['EPSG'] == '4326'
    assert atr['VERTICAL_DATUM'] == 'EGM2008 geoid' and float(atr['X_STEP']) == pytest.approx(RES)
    assert np.isfinite(dem).all()
    sea = (lat > 30 + 1.5 * RES) & (lon > 130 + 1.5 * RES)
    land = (lat < 30 - 1.5 * RES) | (lon < 130 - 1.5 * RES)
    assert np.all(dem[sea] == 0)
    err = np.abs(dem - steep(lon, lat))
    # the plane comes back exactly across the seams; what is left is the curve of
    # the ridge between 10" pixels
    assert err[land].max() < 1.0
    seam = land & ((np.abs(lon - 130) < 1.5 * RES) | (np.abs(lat - 30) < 1.5 * RES))
    assert seam.sum() > 100 and err[seam].max() < 1.0
    assert dem[0].mean() > dem[-1].mean()                    # north up


def test_utm_dem(monkeypatch, tmp_path):
    serve_glo30(monkeypatch, tmp_path / 'server', CELLS)
    fwd = pyproj.Transformer.from_crs(4326, 32652, always_xy=True)
    (e0, e1), (n1, n0) = fwd.transform([129.7, 130.3], [30.3, 29.7])
    path = terrain.copernicus_dem((e0, n1, e1, n0), 32652, True, verbose=False)
    dem, atr = _read(path)
    assert atr['UTM_ZONE'] == '52N' and float(atr['X_STEP']) == pytest.approx(
        RES * terrain.M_PER_DEG)
    east, north = _centres(atr)
    lon, lat = pyproj.Transformer.from_crs(32652, 4326, always_xy=True).transform(east, north)
    assert np.isfinite(dem).all() and np.abs(dem - height(lon, lat)).max() < 1.0


def test_a_big_map_is_read_from_an_overview(monkeypatch, tmp_path):
    requests = serve_glo30(monkeypatch, tmp_path / 'server', CELLS)
    path = terrain.copernicus_dem((129.1, 30.9, 130.9, 29.1), 4326, False, max_side=100,
                                  verbose=False)
    dem, atr = _read(path)
    assert max(dem.shape) <= 100 and atr['COP30_READ_AT'].startswith('40"')
    # the headers, then the one block of the quarter-resolution level of each tile
    blocks = [r for r in requests if r[1] > 0]
    assert len(blocks) == 4 and len({r[0] for r in blocks}) == 4
    lon, lat = _centres(atr)
    assert np.abs(dem - height(lon, lat)).max() < 25.0     # 40" pixels on the ridge


def test_a_second_run_downloads_nothing(monkeypatch, tmp_path):
    requests = serve_glo30(monkeypatch, tmp_path / 'server', CELLS)
    first = terrain.copernicus_dem(BOX, 4326, False, verbose=False)
    assert requests
    requests.clear()
    assert terrain.copernicus_dem(BOX, 4326, False, verbose=False) == first and not requests
    os.remove(first)                                         # the blocks are kept, too
    terrain.copernicus_dem(BOX, 4326, False, verbose=False)
    assert not requests
    kept = os.path.join(terrain.cache_dir(), 'cop30', terrain.cop30_name(30, 130))
    assert os.path.isfile(os.path.join(kept, 'head.bin'))


def test_a_failed_tile_is_left_empty_not_made_sea(monkeypatch, tmp_path, capsys):
    serve_glo30(monkeypatch, tmp_path / 'server', CELLS[:3])        # (30, 130): sea
    served = terrain.fetch_range

    def flaky(url, start, size):
        if terrain.cop30_name(29, 130) in url:
            raise OSError('connection reset')
        return served(url, start, size)

    monkeypatch.setattr(terrain, 'fetch_range', flaky)
    dem, atr = _read(terrain.copernicus_dem(BOX, 4326, False))
    out = capsys.readouterr().out
    assert '1 of 4 GLO-30 cells have no tile' in out and '1 of 4 GLO-30 tiles failed' in out
    lon, lat = _centres(atr)
    east = lon > 130 + 1.5 * RES
    assert np.isnan(dem[east & (lat < 30 - 1.5 * RES)]).all()       # unknown: empty
    assert np.all(dem[east & (lat > 30 + 1.5 * RES)] == 0)          # the sea: 0 m
    assert np.isfinite(dem[lon < 130 - 1.5 * RES]).all()


def test_no_network_says_so(monkeypatch, tmp_path):
    serve_glo30(monkeypatch, tmp_path / 'server', [])

    def offline(url, start, size):
        raise OSError('no route to host')

    monkeypatch.setattr(terrain, 'fetch_range', offline)
    with pytest.raises(RuntimeError, match='could not download any Copernicus GLO-30 tile'):
        terrain.copernicus_dem(BOX, 4326, False, verbose=False)


##############################################################################
#                               the CLI                                      #
##############################################################################
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')


def _peak_expected(obj, sec):
    sw = sec.swath_list[0]
    lon0 = float(obj.atr['X_FIRST']) + (sw['start_yx'][1] + 0.5) * float(obj.atr['X_STEP'])
    lon1 = float(obj.atr['X_FIRST']) + (sw['end_yx'][1] + 0.5) * float(obj.atr['X_STEP'])
    return (FAULT_LON - lon0) / (lon1 - lon0) * sw['length']


def test_dem_is_glo30_for_the_map_and_the_profile(velocity_file, monkeypatch, tmp_path, capsys):
    serve_glo30(monkeypatch, tmp_path / 'server', CELLS)
    out = tmp_path / 'out'
    obj = cli.swath_main([velocity_file, *LINE, '--width', '6', '--dem', '--topoxsect',
                          '--save', '--nodisplay', '--outdir', str(out), '--noverbose'])
    # a small map: one download serves both
    assert obj.dem_file == obj.topoxsect == obj.dem_download
    assert os.path.basename(obj.dem_file).startswith('dem_cop30_')
    assert obj.dem_sources[obj.dem_file] == terrain.COP30_ATTRIBUTION
    assert f'DEM: {terrain.COP30_ATTRIBUTION}' in capsys.readouterr().out
    assert obj.transparency == pytest.approx(0.8)          # the relief shows through
    assert len(obj.ax_img.images) >= 2
    sec = obj.sections[0]
    assert sec.topo_error is None
    peak = sec.topo_den['x_centers'][np.nanargmax(sec.topo_den['median'])]
    assert peak == pytest.approx(_peak_expected(obj, sec), abs=1.0)
    info = next(out.glob('*_info.txt')).read_text()
    assert 'DEM source     : Copernicus DEM GLO-30: produced using Copernicus WorldDEM-30' in info


def test_a_big_map_gets_a_lighter_background(velocity_file, monkeypatch, tmp_path):
    """The map's DEM stops at MAP_SIDE pixels a side; the profile's goes on to
    MAX_SIDE - two files, one download."""
    from viewinsar import viewer

    requests = serve_glo30(monkeypatch, tmp_path / 'server', CELLS)
    monkeypatch.setattr(viewer, 'MAP_SIDE', 60)
    obj = cli.swath_main([velocity_file, *LINE, '--width', '6', '--dem', '--topoxsect',
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    assert obj.dem_file != obj.topoxsect and obj.dem_download == obj.topoxsect
    background, _ = _read(obj.dem_file)
    profile, _ = _read(obj.topoxsect)
    assert max(background.shape) <= 60 < max(profile.shape)
    assert all(r[0].startswith('Copernicus_DSM_COG_10_N') for r in requests)


def test_terrain_tiles_when_glo30_cannot_be_had(velocity_file, monkeypatch, tmp_path, capsys):
    from test_terrain import synthetic_tile

    monkeypatch.setattr(terrain, 'fetch_tile',
                        lambda z, x, y, folder=None, url=None: synthetic_tile(z, x, y))

    def offline(url, start, size):
        raise OSError('blocked')

    monkeypatch.setattr(terrain, 'fetch_range', offline)
    obj = cli.swath_main([velocity_file, *LINE, '--width', '6', '--dem', '--topoxsect',
                          '--dem-zoom', '9', '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    out = capsys.readouterr().out
    assert 'WARNING: no Copernicus GLO-30 DEM' in out and 'the Terrain Tiles instead' in out
    assert out.count('no Copernicus GLO-30 DEM') == 1        # not tried again for the map
    assert os.path.basename(obj.topoxsect).startswith('dem_terrarium_')
    assert obj.dem_sources[obj.topoxsect] == terrain.ATTRIBUTION
    assert obj.sections[0].topo_den is not None
    # GLO-30 asked for by name: no other source
    obj = cli.swath_main([velocity_file, *LINE, '--width', '6', '--topoxsect',
                          '--dem-source', 'cop30', '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    out = capsys.readouterr().out
    assert 'no DEM for the map' in out and 'Terrain Tiles instead' not in out
    assert obj.topoxsect is None and obj.sections[0].has_data


def test_dem_source_terrain(velocity_file, monkeypatch, tmp_path, capsys):
    from test_terrain import synthetic_tile

    requests = serve_glo30(monkeypatch, tmp_path / 'server', CELLS)
    monkeypatch.setattr(terrain, 'fetch_tile',
                        lambda z, x, y, folder=None, url=None: synthetic_tile(z, x, y))
    obj = cli.swath_main([velocity_file, *LINE, '--width', '6', '--dem', '--dem-source',
                          'Terrain', '--dem-zoom', '9', '--nodisplay', '--outdir',
                          str(tmp_path), '--noverbose'])
    assert not requests and os.path.basename(obj.dem_file).startswith('dem_terrarium_')
    assert 'NOTE: --dem-zoom' not in capsys.readouterr().out


def test_dem_options(velocity_file, capsys):
    cli.cmd_line_parse([velocity_file, '--dem', '--dem-zoom', '12', '--noverbose'])
    assert 'NOTE: --dem-zoom 12 is for the Terrain Tiles' in capsys.readouterr().out
    inps = cli.cmd_line_parse([velocity_file, '--dem', '--dem-source', 'COP30', '--noverbose'])
    assert inps.dem_source == 'cop30'
    with pytest.raises(SystemExit):
        cli.cmd_line_parse([velocity_file, '--dem', '--dem-source', 'srtm', '--noverbose'])

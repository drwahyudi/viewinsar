"""--dem-source terrain: the Terrain Tiles download, mosaic and re-grid (also
what --dem / --topoxsect fall back to when Copernicus GLO-30 cannot be had).

No network here: the tiles are synthesised with the real Terrarium encoding
and Web Mercator geometry from a known height field, and served in place of
the download. A north-south height gradient catches any flipped axis.
"""

import io
import urllib.error

import numpy as np
import pytest

from conftest import FAULT_LON
from viewinsar import cli, terrain

pyproj = pytest.importorskip('pyproj')
h5py = pytest.importorskip('h5py')
Image = pytest.importorskip('PIL.Image')

LINE = ['--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75']
R_EARTH = 6378137.0


def height(lon, lat):
    """A ridge on the fixtures' fault, plus a north-up gradient (100 m per degree)."""
    km_east = (lon - FAULT_LON) * 111.32 * np.cos(np.radians(30.3))
    return 1500.0 * np.exp(-(km_east / 7.0) ** 2) + 40.0 + 100.0 * (lat - 30.3)


def encode_terrarium(h):
    val = np.asarray(h, dtype=np.float64) + 32768.0
    red = np.floor(val / 256.0)
    green = np.floor(val - red * 256.0)
    blue = np.floor((val - np.floor(val)) * 256.0)
    rgb = np.stack([red, green, blue], axis=-1).clip(0, 255).astype(np.uint8)
    buf = io.BytesIO()
    Image.fromarray(rgb, mode='RGB').save(buf, format='PNG')
    return buf.getvalue()


def synthetic_tile(zoom, x, y, field=height):
    """The PNG the server would send for tile z/x/y of `field`."""
    centre = np.arange(terrain.TILE) + 0.5
    scale = terrain.TILE * 2 ** zoom / (2 * terrain.HALF)
    mx = (x * terrain.TILE + centre[None, :]) / scale - terrain.HALF
    my = terrain.HALF - (y * terrain.TILE + centre[:, None]) / scale
    lon = np.degrees(mx / R_EARTH) + 0 * my
    lat = np.degrees(2 * np.arctan(np.exp(my / R_EARTH)) - np.pi / 2) + 0 * mx
    return encode_terrarium(field(lon, lat))


@pytest.fixture
def served(monkeypatch, tmp_path):
    """Serve synthetic tiles instead of the internet; count the requests."""
    calls = []

    def fake_fetch(zoom, x, y, folder=None, url=terrain.TERRARIUM_URL):
        calls.append((zoom, x, y))
        return synthetic_tile(zoom, x, y)

    monkeypatch.setattr(terrain, 'fetch_tile', fake_fetch)
    monkeypatch.setenv('VIEWINSAR_CACHE', str(tmp_path / 'cache'))
    return calls


##############################################################################
#                               units                                        #
##############################################################################
def test_decode_terrarium_is_exact():
    known = np.array([[-32768.0, 0.0, 1234.5], [8848.0, -10994.75, 0.00390625]])
    assert np.array_equal(terrain.decode_terrarium(encode_terrarium(known)),
                          known.astype(np.float32))


def test_tile_maths_matches_pyproj():
    rng = np.random.default_rng(3)
    lon, lat = rng.uniform(-179, 179, 50), rng.uniform(-80, 80, 50)
    x, y = terrain.lonlat_to_tile(lon, lat, 10)
    mx, my = pyproj.Transformer.from_crs(4326, 3857, always_xy=True).transform(lon, lat)
    assert np.allclose(x, (mx + terrain.HALF) / (2 * terrain.HALF) * 2 ** 10, atol=1e-6)
    assert np.allclose(y, (terrain.HALF - my) / (2 * terrain.HALF) * 2 ** 10, atol=1e-6)
    assert terrain.lonlat_to_tile(0.0, 0.0, 1) == (1.0, 1.0)


def test_zoom_respects_the_budget():
    sulawesi = (118.5, -6.5, 125.5, 2.0)              # about the map in the report
    zoom = terrain.choose_zoom(sulawesi, max_tiles=300)

    def count(z):
        tx0, tx1, ty0, ty1 = terrain._tile_range(sulawesi, z)
        return (tx1 - tx0 + 1) * (ty1 - ty0 + 1)

    assert count(zoom) <= 300 < count(zoom + 1)
    assert terrain.choose_zoom((130.1, 30.1, 130.8, 30.5)) == terrain.MAX_ZOOM_AUTO
    assert terrain.tile_resolution(12) == pytest.approx(38.2, abs=0.1)


def test_tiles_are_cached_on_disk(monkeypatch, tmp_path):
    hits = []

    def fake_get(url, timeout=30, retries=2):
        hits.append(url)
        return synthetic_tile(10, 881, 422)

    monkeypatch.setattr(terrain, '_get', fake_get)
    first = terrain.fetch_tile(10, 881, 422, folder=str(tmp_path))
    second = terrain.fetch_tile(10, 881, 422, folder=str(tmp_path))
    assert first == second and len(hits) == 1
    assert hits[0] == 'https://s3.amazonaws.com/elevation-tiles-prod/terrarium/10/881/422.png'
    assert (tmp_path / 'terrarium' / '10' / '881' / '422.png').is_file()


##############################################################################
#                       mosaic + re-grid                                     #
##############################################################################
def _read(path):
    with h5py.File(path, 'r') as fid:
        return fid['height'][:], {k: (v.decode() if isinstance(v, bytes) else v)
                                  for k, v in fid.attrs.items()}


def _pixel_centres(atr):
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    xs = float(atr['X_FIRST']) + (np.arange(width) + 0.5) * float(atr['X_STEP'])
    ys = float(atr['Y_FIRST']) + (np.arange(length) + 0.5) * float(atr['Y_STEP'])
    return np.meshgrid(xs, ys)


def test_lonlat_grid_matches_the_field(served, tmp_path):
    path = terrain.terrain_dem((130.0, 30.6, 130.9, 30.0), 4326, projected=False, zoom=10,
                               verbose=False)
    dem, atr = _read(path)
    lon, lat = _pixel_centres(atr)
    err = dem - height(lon, lat)
    assert np.isfinite(dem).all()
    assert np.abs(err).max() < 5.0 and abs(err.mean()) < 0.5
    assert atr['EPSG'] == '4326' and atr['FILE_TYPE'] == 'geometry'
    assert 'Terrain Tiles' in atr['ATTRIBUTION']
    # the north-south gradient comes out the right way round
    assert dem[0].mean() > dem[-1].mean()


def test_utm_grid_matches_the_field(served):
    fwd = pyproj.Transformer.from_crs(4326, 32652, always_xy=True)
    (e0, e1), (n1, n0) = fwd.transform([130.0, 130.9], [30.6, 30.0])
    path = terrain.terrain_dem((e0, n1, e1, n0), 32652, projected=True, zoom=10,
                               verbose=False)
    dem, atr = _read(path)
    assert atr['UTM_ZONE'] == '52N' and atr['X_UNIT'] == 'meters'
    east, north = _pixel_centres(atr)
    lon, lat = pyproj.Transformer.from_crs(32652, 4326, always_xy=True).transform(east, north)
    err = dem - height(lon, lat)
    assert np.isfinite(dem).all() and np.abs(err).max() < 5.0


def test_a_second_run_downloads_nothing(served):
    box = (130.0, 30.6, 130.9, 30.0)
    first = terrain.terrain_dem(box, 4326, projected=False, zoom=9, verbose=False)
    n_calls = len(served)
    again = terrain.terrain_dem(box, 4326, projected=False, zoom=9, verbose=False)
    assert again == first and len(served) == n_calls > 0


def test_no_network_says_so(monkeypatch, tmp_path):
    def offline(*args, **kwargs):
        raise urllib.error.URLError('no route to host')

    monkeypatch.setattr(terrain, 'fetch_tile', offline)
    monkeypatch.setenv('VIEWINSAR_CACHE', str(tmp_path))
    with pytest.raises(RuntimeError, match='internet access'):
        terrain.terrain_dem((130.0, 30.6, 130.9, 30.0), 4326, projected=False, zoom=8,
                            verbose=False)


##############################################################################
#                               the CLI                                      #
##############################################################################
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')


@pytest.mark.parametrize('args, dem, topo', [
    (['-d'], 'auto', None),
    (['-d', 'auto'], 'auto', None),
    (['--dem', '--topoxsect'], 'auto', 'auto'),
    (['--topoxsect'], None, 'auto'),
    (['-topoxsect', '-w', '5'], None, 'auto'),
])
def test_dem_without_a_file_parses(velocity_file, args, dem, topo):
    inps = cli.cmd_line_parse([velocity_file, *args, '--noverbose'])
    assert inps.dem_file == dem and inps.topoxsect == topo
    assert '-d' not in inps.view_extra


def test_dem_file_is_still_view_pys_and_feeds_topoxsect(velocity_file, dem_file):
    inps = cli.cmd_line_parse([velocity_file, '-d', dem_file, '--topoxsect', '--noverbose'])
    assert inps.dem_file == dem_file and inps.topoxsect == dem_file
    assert inps.view_extra[-2:] == ['--dem', dem_file]          # view.py's spelling too


def test_bad_dem_options(velocity_file):
    with pytest.raises(ValueError, match='no such file'):
        cli.cmd_line_parse([velocity_file, '-d', 'missing.tif', '--noverbose'])
    with pytest.raises(ValueError, match='0-15'):
        cli.cmd_line_parse([velocity_file, '-d', '--dem-zoom', '16', '--noverbose'])


def test_download_for_the_map_and_the_profile(velocity_file, served, tmp_path):
    out = tmp_path / 'out'
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-d', '--topoxsect',
                          '--dem-source', 'terrain', '--dem-zoom', '10', '--save',
                          '--nodisplay', '--outdir', str(out), '--noverbose'])
    # one download serves both
    assert obj.dem_file == obj.topoxsect == obj.dem_download
    assert obj.dem_file.endswith('.h5') and str(tmp_path / 'cache') in obj.dem_file
    assert obj.transparency == pytest.approx(0.8)          # the relief shows through
    assert len(obj.ax_img.images) >= 2                     # relief + data on the map

    sec = obj.sections[0]
    assert sec.topo_error is None
    peak = sec.topo_den['x_centers'][np.nanargmax(sec.topo_den['median'])]
    sw = sec.swath_list[0]
    lon0 = float(obj.atr['X_FIRST']) + (sw['start_yx'][1] + 0.5) * float(obj.atr['X_STEP'])
    lon1 = float(obj.atr['X_FIRST']) + (sw['end_yx'][1] + 0.5) * float(obj.atr['X_STEP'])
    assert peak == pytest.approx((FAULT_LON - lon0) / (lon1 - lon0) * sw['length'], abs=1.0)
    info = next(out.glob('*_info.txt')).read_text()
    assert 'DEM source     : Terrain Tiles' in info


def test_topoxsect_alone_downloads_without_a_background(velocity_file, served, tmp_path):
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '--topoxsect', '--dem-source',
                          'terrain', '--dem-zoom', '9', '--nodisplay', '--outdir',
                          str(tmp_path), '--noverbose'])
    assert obj.dem_file is None and obj.topoxsect == obj.dem_download
    assert obj.sections[0].topo_den is not None
    assert obj.transparency == pytest.approx(1.0)


def test_alpha_given_by_the_user_wins(velocity_file, served, tmp_path):
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-d', '--dem-source', 'terrain',
                          '--dem-zoom', '9', '--alpha', '0.6', '--nodisplay', '--outdir',
                          str(tmp_path), '--noverbose'])
    assert obj.transparency == pytest.approx(0.6)


def test_offline_keeps_going(velocity_file, monkeypatch, tmp_path, capsys):
    def offline(*args, **kwargs):
        raise urllib.error.URLError('no route to host')

    monkeypatch.setattr(terrain, 'fetch_tile', offline)
    monkeypatch.setenv('VIEWINSAR_CACHE', str(tmp_path))
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-d', '--topoxsect',
                          '--dem-zoom', '8', '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    out = capsys.readouterr().out
    assert 'no Copernicus GLO-30 DEM' in out                # first GLO-30, then
    assert 'no DEM for the map' in out and 'could not download any Terrain Tile' in out
    assert obj.dem_file is None and obj.topoxsect is None
    assert obj.sections[0].has_data                         # the profile still works

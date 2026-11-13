"""Radar-coded data: the DEM, the elevation profile and the faults, all through
MintPy's geometryRadar.h5 (latitude / longitude / height of every pixel).

The synthetic product mimics a DESCENDING right-looking track: rows run south,
columns run WEST, with a slight rotation - so view.py flips the map left-right,
as it does for real data.
"""

import json
import os

import numpy as np
import pytest

from conftest import FAULT_LON
from viewinsar import cli
from viewinsar.radar import RadarLookup, find_lookup
from viewinsar.raster import read_any

h5py = pytest.importorskip('h5py')
pytest.importorskip('mintpy', reason='the viewer needs MintPy')
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

LENGTH, WIDTH = 200, 300
R_EARTH, SAT_HEIGHT = 6371000.0, 628000.0


def lat_of(row, col):
    return 30.62 - 0.0018 * row - 0.0002 * col


def lon_of(row, col):
    return 130.92 - 0.0025 * col + 0.0003 * row


def ridge(lon):
    km_east = (lon - FAULT_LON) * 111.32 * np.cos(np.radians(30.4))
    return 1500.0 * np.exp(-(km_east / 7.0) ** 2) + 40.0


def fault_col(row):
    """The column where the fault (lon = FAULT_LON) crosses a given row."""
    return (130.92 + 0.0003 * row - FAULT_LON) / 0.0025


def _radar_attrs(file_type, unit):
    inc = np.radians(38.0)
    # slant range for a 38 deg incidence from 628 km up
    look = np.arcsin(R_EARTH / (R_EARTH + SAT_HEIGHT) * np.sin(inc))
    srange = (R_EARTH + SAT_HEIGHT) * np.cos(look) - R_EARTH * np.cos(inc)
    return {
        'FILE_TYPE': file_type, 'UNIT': unit, 'DATA_TYPE': 'float32',
        'LENGTH': str(LENGTH), 'WIDTH': str(WIDTH),
        'PLATFORM': 'Alos2', 'PROCESSOR': 'isce', 'ORBIT_DIRECTION': 'DESCENDING',
        'WAVELENGTH': '0.2290', 'HEADING': '-169.0',
        'RANGE_PIXEL_SIZE': str(240.0 * np.sin(inc)),
        'AZIMUTH_PIXEL_SIZE': str(200.0 * (R_EARTH + SAT_HEIGHT) / R_EARTH),
        'STARTING_RANGE': str(srange - 150 * 240.0 * np.sin(inc)),
        'EARTH_RADIUS': str(R_EARTH), 'HEIGHT': str(SAT_HEIGHT),
        'ALOOKS': '1', 'RLOOKS': '1', 'REF_Y': '100', 'REF_X': '150',
        'LAT_REF1': '30.62', 'LON_REF1': '130.92', 'LAT_REF2': '30.56', 'LON_REF2': '130.17',
        'LAT_REF3': '30.26', 'LON_REF3': '130.98', 'LAT_REF4': '30.20', 'LON_REF4': '130.23',
    }


def _write(path, datasets, attrs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with h5py.File(path, 'w') as fid:
        for name, data in datasets.items():
            fid.create_dataset(name, data=np.asarray(data, dtype=np.float32))
        for key, val in attrs.items():
            fid.attrs[key] = val
    return str(path)


def make_radar_product(folder, geometry='inputs', with_height=True):
    """velocity.h5 in `folder`, geometryRadar.h5 at `geometry` relative to it."""
    rows, cols = np.mgrid[0:LENGTH, 0:WIDTH].astype(np.float64)
    lat, lon = lat_of(rows, cols), lon_of(rows, cols)
    km_east = (lon - FAULT_LON) * 111.32 * np.cos(np.radians(30.4))
    rng = np.random.default_rng(5)
    vel = 0.30 * (2 / np.pi) * np.arctan(km_east / 8.0) + rng.normal(0, 0.03, lat.shape)
    vfile = _write(os.path.join(folder, 'velocity.h5'), {'velocity': vel},
                   _radar_attrs('velocity', 'm/year'))
    gfile = None
    if geometry is not None:
        # the real file's datasets, azimuthAngle FIRST: the height must still win
        datasets = {'azimuthAngle': np.full(lat.shape, 100.0), 'incidenceAngle':
                    np.full(lat.shape, 38.0), 'latitude': lat, 'longitude': lon,
                    'slantRangeDistance': np.full(lat.shape, 8e5)}
        if with_height:
            datasets['height'] = ridge(lon)
        gfile = _write(os.path.join(folder, geometry, 'geometryRadar.h5'), datasets,
                       _radar_attrs('geometry', '1'))
    return vfile, gfile


@pytest.fixture
def product(tmp_path):
    return make_radar_product(str(tmp_path / 'proc'))


@pytest.fixture
def fault_geojson(tmp_path):
    obj = {'type': 'FeatureCollection', 'features': [{
        'type': 'Feature', 'properties': {},
        'geometry': {'type': 'LineString',
                     'coordinates': [[FAULT_LON, 30.0], [FAULT_LON, 30.8]]}}]}
    path = tmp_path / 'fault.geojson'
    path.write_text(json.dumps(obj))
    return str(path)


SECTION = ['--start-yx', '100', '120', '--end-yx', '100', '280']   # crosses the fault


def _t_fault():
    """Where along SECTION (0..1) the fault is, from the analytic geometry."""
    return (fault_col(100) - 120) / (280 - 120)


##############################################################################
#                              the lookup                                    #
##############################################################################
def test_lookup_found_like_mintpy_does(tmp_path):
    from mintpy.utils import readfile

    vfile, gfile = make_radar_product(str(tmp_path / 'a'), geometry='inputs')
    atr = readfile.read_attribute(vfile)
    assert find_lookup(vfile, atr, verbose=False) == os.path.normpath(gfile)

    # the user's layout: velocity in a backup folder, inputs/ one level up
    vfile, _ = make_radar_product(str(tmp_path / 'proc' / 'velocity2_backup'), geometry=None)
    _, gfile = make_radar_product(str(tmp_path / 'proc' / 'other'), geometry='../inputs')
    assert find_lookup(vfile, atr, verbose=False) == os.path.normpath(gfile)


def test_a_lookup_of_another_size_is_refused(tmp_path):
    from mintpy.utils import readfile

    vfile, gfile = make_radar_product(str(tmp_path / 'b'))
    atr = dict(readfile.read_attribute(vfile), LENGTH='999')
    assert find_lookup(vfile, atr, verbose=False) is None
    with pytest.raises(ValueError, match='999'):
        find_lookup(vfile, atr, given=gfile, verbose=False)


def test_lookup_round_trip(product):
    lookup = RadarLookup(product[1])
    rows, cols = np.array([3.0, 57.5, 120.25, 198.0]), np.array([1.0, 290.5, 150.0, 77.75])
    lon, lat = lookup.to_lonlat(rows, cols)
    assert np.allclose(lon, lon_of(rows, cols)) and np.allclose(lat, lat_of(rows, cols))
    rr, cc = lookup.to_rowcol(lon, lat)
    assert np.allclose(rr, rows, atol=1e-6) and np.allclose(cc, cols, atol=1e-6)
    assert np.isnan(lookup.to_rowcol(np.array([100.0]), np.array([0.0]))[0]).all()


def test_a_geometry_file_gives_its_height_not_the_first_dataset(product):
    data = read_any(product[1], box=(0, 0, 3, 3))[0]
    assert np.allclose(data, ridge(lon_of(*np.mgrid[0:3, 0:3])), atol=0.01)


##############################################################################
#                              end to end                                    #
##############################################################################
def _peak_km(sec):
    tden = sec.topo_den
    return tden['x_centers'][np.nanargmax(tden['median'])]


def test_faults_on_radar_data(product, fault_geojson, tmp_path):
    out = tmp_path / 'out'
    obj = cli.swath_main([product[0], *SECTION, '-w', '4', '-fault', fault_geojson,
                          '--save', '--nodisplay', '--outdir', str(out), '--noverbose'])
    sec = obj.sections[0]
    length = sec.swath_list[0]['length']
    # float32 lat/lon in the lookup (as MintPy writes them) resolve ~1 m at 130 E
    assert sec.crossings == pytest.approx([_t_fault() * length], abs=0.002)
    # on the map: the trace runs through the radar pixels of lon = FAULT_LON
    xs, ys = obj.fault_obj.get_xdata(), obj.fault_obj.get_ydata()
    ok = np.isfinite(xs)
    assert np.allclose(xs[ok], fault_col(ys[ok]), atol=1e-3)
    info = next(out.glob('*_info.txt')).read_text()
    line = next(ln for ln in info.splitlines() if ln.startswith('fault crossing'))
    assert f'y/x 100.0, {fault_col(100):.1f}' in line
    lat, lon = (float(v) for v in line.split('lat/lon')[1].split(','))
    assert lon == pytest.approx(FAULT_LON, abs=1e-5)     # and its real lat/lon
    assert lat == pytest.approx(lat_of(100, fault_col(100)), abs=1e-5)


def test_dem_from_the_geometry_file(product, tmp_path):
    """-d --topoxsect on radar data: the processing DEM, no download."""
    obj = cli.swath_main([product[0], *SECTION, '-w', '4', '-d', '--topoxsect',
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    assert obj.dem_file == obj.topoxsect == os.path.normpath(product[1])
    assert obj.dem_download is None and obj.topo_mode == 'grid'
    assert len(obj.ax_img.images) >= 2                  # shaded relief behind the data
    sec = obj.sections[0]
    assert sec.topo_error is None
    assert _peak_km(sec) == pytest.approx(_t_fault() * sec.swath_list[0]['length'], abs=1.0)


def test_geocoded_dem_on_radar_data(product, dem_file, tmp_path):
    """A geocoded DEM file (here UTM) reaches the radar section through the lookup."""
    obj = cli.swath_main([product[0], *SECTION, '-w', '4', '--topoxsect', dem_file,
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    sec = obj.sections[0]
    assert obj.topo_mode == 'lookup' and sec.topo_error is None
    assert _peak_km(sec) == pytest.approx(_t_fault() * sec.swath_list[0]['length'], abs=1.5)


def test_download_onto_the_radar_grid(tmp_path, monkeypatch):
    """A geometry file without height: Terrain Tiles resampled onto the radar grid."""
    from test_terrain import synthetic_tile
    from viewinsar import terrain

    vfile, gfile = make_radar_product(str(tmp_path / 'proc'), with_height=False)
    monkeypatch.setattr(terrain, 'fetch_tile',
                        lambda z, x, y, folder=None, url=None: synthetic_tile(z, x, y))
    monkeypatch.setenv('VIEWINSAR_CACHE', str(tmp_path / 'cache'))
    obj = cli.swath_main([vfile, *SECTION, '-w', '4', '-d', '--topoxsect', '--dem-zoom', '11',
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    assert obj.dem_download and obj.topoxsect == obj.dem_download == obj.dem_file
    with h5py.File(obj.dem_download, 'r') as fid:
        assert fid['height'].shape == (LENGTH, WIDTH) and 'Y_FIRST' not in fid.attrs
        dem = fid['height'][:]
    rows, cols = np.mgrid[0:LENGTH, 0:WIDTH]
    # the synthetic tiles also carry a north-south gradient of 100 m per degree
    expect = ridge(lon_of(rows, cols)) + 100.0 * (lat_of(rows, cols) - 30.3)
    assert np.nanmax(np.abs(dem - expect)) < 6.0
    sec = obj.sections[0]
    assert _peak_km(sec) == pytest.approx(_t_fault() * sec.swath_list[0]['length'], abs=1.0)


def test_no_geometry_file_says_what_is_missing(tmp_path, fault_geojson, capsys):
    vfile, _ = make_radar_product(str(tmp_path / 'lonely'), geometry=None)
    obj = cli.swath_main([vfile, *SECTION, '-w', '4', '-d', '--topoxsect',
                          '-fault', fault_geojson, '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    out = capsys.readouterr().out
    assert 'geometryRadar.h5' in out and '--lookup' in out
    assert obj.sections[0].has_data                     # the profile itself still works
    assert obj.topoxsect is None and not obj.faults


def test_lookup_given_by_hand(tmp_path, fault_geojson):
    vfile, _ = make_radar_product(str(tmp_path / 'data'), geometry=None)
    _, gfile = make_radar_product(str(tmp_path / 'elsewhere'), geometry='geom')
    obj = cli.swath_main([vfile, *SECTION, '-w', '4', '--lookup', gfile, '--topoxsect',
                          '-fault', fault_geojson, '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    assert obj.topoxsect == gfile and obj.sections[0].crossings.size == 1
    with pytest.raises(ValueError, match='--lookup'):
        cli.cmd_line_parse([vfile, '--lookup', 'missing.h5', '--noverbose'])

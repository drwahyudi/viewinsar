"""Synthetic MintPy-format products, for the test suite."""

import numpy as np
import pytest

LENGTH, WIDTH = 300, 450
LAT0, LON0, DLAT, DLON = 30.60, 130.00, -0.002, 0.002
FAULT_LON = 130.45
KM_PER_DEG_LON = 111.32 * np.cos(np.deg2rad(30.3))


def _base_attrs():
    return {
        'FILE_TYPE': 'velocity', 'UNIT': 'm/year', 'DATA_TYPE': 'float32',
        'LENGTH': str(LENGTH), 'WIDTH': str(WIDTH),
        'Y_FIRST': str(LAT0), 'X_FIRST': str(LON0),
        'Y_STEP': str(DLAT), 'X_STEP': str(DLON),
        'Y_UNIT': 'degrees', 'X_UNIT': 'degrees',
        'REF_Y': '150', 'REF_X': '225',
        'WAVELENGTH': '0.055465763', 'PLATFORM': 'Sen', 'PROCESSOR': 'isce',
        'ORBIT_DIRECTION': 'ASCENDING', 'ALOOKS': '1', 'RLOOKS': '1',
    }


def make_velocity():
    """An arctan interseismic profile across a N-S fault, plus noise."""
    rng = np.random.default_rng(42)
    lon = LON0 + (np.arange(WIDTH) + 0.5) * DLON
    lat = LAT0 + (np.arange(LENGTH) + 0.5) * DLAT
    lon2d, _ = np.meshgrid(lon, lat)
    xperp = (lon2d - FAULT_LON) * KM_PER_DEG_LON              # km east of the fault
    vel = 0.30 * (2 / np.pi) * np.arctan(xperp / 8.0)         # +/- 0.30 m/year
    vel += rng.normal(0, 0.05, vel.shape)
    vel = vel.astype(np.float32)
    vel[:20, :] = np.nan                                      # a decorrelated strip
    return vel


@pytest.fixture(scope='session')
def velocity_file(tmp_path_factory):
    h5py = pytest.importorskip('h5py')
    path = tmp_path_factory.mktemp('data') / 'velocity.h5'
    with h5py.File(path, 'w') as fid:
        fid.create_dataset('velocity', data=make_velocity())
        for key, val in _base_attrs().items():
            fid.attrs[key] = val
    return str(path)


@pytest.fixture(scope='session')
def timeseries_file(tmp_path_factory, velocity_file):
    h5py = pytest.importorskip('h5py')
    rng = np.random.default_rng(7)
    with h5py.File(velocity_file, 'r') as fid:
        vel = fid['velocity'][:]

    dates = [f'2020{m:02d}12' for m in range(1, 13)] + [f'2021{m:02d}12' for m in range(1, 13)]
    tyr = np.arange(len(dates)) / 12.
    cube = (vel[None, :, :] * tyr[:, None, None]).astype(np.float32)
    cube += rng.normal(0, 0.003, cube.shape).astype(np.float32)
    cube[0] = 0

    path = tmp_path_factory.mktemp('data_ts') / 'timeseries.h5'
    attrs = _base_attrs()
    attrs.update(FILE_TYPE='timeseries', UNIT='m')
    with h5py.File(path, 'w') as fid:
        fid.create_dataset('timeseries', data=cube)
        fid.create_dataset('date', data=np.array(dates, dtype=np.bytes_))
        for key, val in attrs.items():
            fid.attrs[key] = val
    return str(path)


@pytest.fixture(scope='session')
def velocity_attrs():
    return _base_attrs()


@pytest.fixture(scope='session')
def coarse_velocity_file(tmp_path_factory):
    """The same scene multilooked by 6 - the sparse-profile case."""
    h5py = pytest.importorskip('h5py')
    nlk = 6
    vel = make_velocity()
    length, width = (LENGTH // nlk) * nlk, (WIDTH // nlk) * nlk
    blocks = vel[:length, :width].reshape(length // nlk, nlk, width // nlk, nlk)
    import warnings
    with warnings.catch_warnings():          # all-NaN blocks are expected here
        warnings.simplefilter('ignore', RuntimeWarning)
        coarse = np.nanmean(blocks, axis=(1, 3)).astype(np.float32)

    path = tmp_path_factory.mktemp('coarse') / 'velocity_ml6.h5'
    attrs = _base_attrs()
    attrs.update(LENGTH=str(coarse.shape[0]), WIDTH=str(coarse.shape[1]),
                 Y_STEP=str(DLAT * nlk), X_STEP=str(DLON * nlk),
                 ALOOKS=str(nlk), RLOOKS=str(nlk))
    with h5py.File(path, 'w') as fid:
        fid.create_dataset('velocity', data=coarse)
        for key, val in attrs.items():
            fid.attrs[key] = val
    return str(path)


@pytest.fixture(scope='session')
def dem_file(tmp_path_factory):
    """A DEM as a GeoTIFF in a DIFFERENT CRS (UTM) from the data (lat/lon)."""
    rasterio = pytest.importorskip('rasterio')
    pytest.importorskip('pyproj')
    from pyproj import Transformer
    from rasterio.transform import Affine

    # UTM zone 52N covers 126-132 E; the fixtures sit at ~130.4 E, 30.3 N
    fwd = Transformer.from_crs(4326, 32652, always_xy=True)
    corners_e, corners_n = fwd.transform(
        [LON0, LON0 + WIDTH * DLON], [LAT0, LAT0 + LENGTH * DLAT])
    e0, n0 = min(corners_e) - 3000, max(corners_n) + 3000
    step, dl, dw = 200.0, 420, 560
    ee = e0 + (np.arange(dw) + 0.5) * step
    nn = n0 - (np.arange(dl) + 0.5) * step
    EE, NN = np.meshgrid(ee, nn)
    # a ridge over the fault at 130.45 E, so the topography has a clear peak
    fe, fn = fwd.transform(FAULT_LON, 30.3)
    dem = (1500 * np.exp(-((EE - fe) / 7000.) ** 2)
           + 90 * np.sin(EE / 4000.) * np.cos(NN / 5000.) + 40).astype(np.float32)

    path = tmp_path_factory.mktemp('dem') / 'dem_utm.tif'
    with rasterio.open(path, 'w', driver='GTiff', height=dl, width=dw, count=1,
                       dtype='float32', crs='EPSG:32652', nodata=-32768.0,
                       transform=Affine.from_gdal(e0, step, 0, n0, 0, -step)) as dst:
        dst.write(dem, 1)
    return str(path)


@pytest.fixture(autouse=True)
def _close_figures():
    """Viewers built inside a test leave windows behind: close them after each."""
    yield
    try:
        import matplotlib.pyplot as plt
        plt.close('all')
    except ImportError:
        pass

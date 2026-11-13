"""--topoxsect: the DEM is sampled on its own terms, in its own CRS, and says why
when it cannot be used."""

import re

import numpy as np
import pytest

from conftest import FAULT_LON
from viewinsar import cli
from viewinsar.raster import crs_of, utm_epsg
from viewinsar.viewer import TOPO_MAX_SAMPLES, topo_sampling, unit_label

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

LINE = ['--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75']
UTM = 32652


##############################################################################
#                               units                                        #
##############################################################################
@pytest.mark.parametrize('n_pixels, pixel_km, expect', [
    (500_000, 0.03, (1, 1)),               # 30 m DEM, small corridor: every pixel
    (40_000_000, 0.03, (4, 1)),            # 30 m DEM, huge corridor: every 4th pixel
    (10_000, 1.0, (1, 4)),                 # 1 km DEM: oversampled to fill 0.5 km bins
    (3_000_000, 1.0, (1, 1)),              # ... unless that would blow the budget
])
def test_topo_sampling(n_pixels, pixel_km, expect):
    decimate, nlook = topo_sampling(n_pixels, pixel_km, bin_x=0.5)
    assert (decimate, nlook) == expect
    assert n_pixels * nlook ** 2 / decimate ** 2 <= TOPO_MAX_SAMPLES * 1.01


def test_topo_sampling_ignores_the_data_nlook():
    """The DEM never inherits --nlook: that was what emptied the panel."""
    assert topo_sampling(3_200_000, 0.03, bin_x=0.5) == (1, 1)


@pytest.mark.parametrize('zone, code', [
    ('51S', 32751), ('51 S', 32751), ('UTM zone 51S', 32751), ('52N', 32652),
    ('WGS_1984_UTM_Zone_51S', 32751), ('WGS 84 / UTM zone 49N', 32649),
    ('', None), (None, None), ('61S', None),
])
def test_utm_zone_spellings(zone, code):
    assert utm_epsg(zone) == code


def test_epsg_4326_on_a_metric_grid_is_not_trusted():
    """MintPy's GDAL reader can label a projected GeoTIFF EPSG:4326."""
    atr = dict(Y_FIRST='3360000', X_FIRST='640000', Y_STEP='-60', X_STEP='60',
               LENGTH='100', WIDTH='100', EPSG='4326', X_UNIT='meters', Y_UNIT='meters',
               UTM_ZONE='WGS_1984_UTM_Zone_52N')
    assert crs_of(atr) == UTM
    atr.pop('UTM_ZONE')
    assert crs_of(atr) is None                  # unknown, rather than wrong


def test_unit_label():
    assert unit_label('1') == '' and unit_label('') == '' and unit_label(None) == ''
    assert unit_label('cm/year') == ' [cm/year]'


##############################################################################
#                        DEM fixtures, in UTM                                #
##############################################################################
def _utm_dem(path, lon_range, lat_range, step=60.0):
    """A UTM 52N GeoTIFF over a lon/lat box, with a ridge on the fault."""
    rasterio = pytest.importorskip('rasterio')
    from pyproj import Transformer
    from rasterio.transform import Affine

    fwd = Transformer.from_crs(4326, UTM, always_xy=True)
    lons = [lon_range[0], lon_range[1], lon_range[0], lon_range[1]]
    lats = [lat_range[0], lat_range[0], lat_range[1], lat_range[1]]
    east, north = fwd.transform(lons, lats)
    e0, n0 = min(east), max(north)
    width, length = int((max(east) - e0) / step), int((n0 - min(north)) / step)
    ee = e0 + (np.arange(width) + 0.5) * step
    nn = n0 - (np.arange(length) + 0.5) * step
    fault_e = fwd.transform(FAULT_LON, 30.3)[0]
    dem = 1500 * np.exp(-((ee[None, :] - fault_e) / 7000.) ** 2) + 40 + 0 * nn[:, None]
    with rasterio.open(path, 'w', driver='GTiff', height=length, width=width, count=1,
                       dtype='int16', crs=f'EPSG:{UTM}', nodata=-32768,
                       transform=Affine.from_gdal(e0, step, 0, n0, 0, -step)) as dst:
        dst.write(dem.astype('int16'), 1)
    return str(path)


@pytest.fixture(scope='module')
def fine_dem(tmp_path_factory):
    """60 m: with --nlook 10 the old code asked for ~30 million DEM samples."""
    return _utm_dem(tmp_path_factory.mktemp('dem') / 'dem_60m.tif',
                    (130.08, 130.82), (30.15, 30.45))


def _ridge_peak(sec):
    tden = sec.topo_den
    return tden['x_centers'][np.nanargmax(tden['median'])]


def _fault_km(sec):
    """Distance of FAULT_LON along the section, from its own file's grid."""
    sw, atr = sec.swath_list[0], sec.v.atr
    x0, dx = float(atr['X_FIRST']), float(atr['X_STEP'])
    lon0 = x0 + (sw['start_yx'][1] + 0.5) * dx
    lon1 = x0 + (sw['end_yx'][1] + 0.5) * dx
    return (FAULT_LON - lon0) / (lon1 - lon0) * sw['length']


##############################################################################
#                          end to end                                        #
##############################################################################
def test_big_nlook_and_a_fine_dem(coarse_velocity_file, fine_dem, tmp_path):
    """The reported bug: --nlook 10 + a fine DEM gave 'no topography'."""
    out = tmp_path / 'out'
    obj = cli.swath_main([coarse_velocity_file, *LINE, '-w', '10', '-nlook', '10',
                          '-topoxsect', fine_dem, '--save', '--nodisplay',
                          '--outdir', str(out), '--noverbose'])
    sec = obj.sections[0]
    assert sec.topo_error is None and sec.topo_den is not None
    assert sec.topo_swath['nlook'] == 1                      # NOT the data's 10
    assert sec.topo_swath['value'].size <= TOPO_MAX_SAMPLES
    assert _ridge_peak(sec) == pytest.approx(_fault_km(sec), abs=1.5)

    # the line is the EXACT median - no 50 m histogram steps - and so is the csv
    line = sec.ax_topo.get_lines()[0]
    assert np.allclose(line.get_ydata(), sec.topo_den['median'], equal_nan=True)
    stats = np.genfromtxt(next(out.glob('*_stats.csv')), delimiter=',', names=True)
    ok = np.isfinite(stats['elevation_median_m'])
    assert ok.sum() > 50
    assert np.allclose(stats['elevation_median_m'][ok],
                       np.interp(stats['dist_km'][ok], sec.topo_den['x_centers'],
                                 sec.topo_den['median']))
    info = next(out.glob('*_info.txt')).read_text()
    assert 'DEM sampling   : every pixel of a 60 m DEM' in info
    assert f'(EPSG:{UTM})' in info


def test_dem_covering_part_of_the_line(velocity_file, tmp_path):
    """A DEM that stops half way still draws the half it covers."""
    dem = _utm_dem(tmp_path / 'west.tif', (130.08, 130.50), (30.15, 30.45))
    obj = cli.swath_main([velocity_file, *LINE, '-w', '5', '-topoxsect', dem,
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    sec = obj.sections[0]
    assert sec.topo_error is None
    med, xc = sec.topo_den['median'], sec.topo_den['x_centers']
    assert np.isfinite(med[xc < 25]).all()
    assert np.isnan(med[xc > 40]).all() or xc.max() < 40     # nothing past its edge


def test_the_reason_is_shown_in_the_panel(velocity_file, tmp_path):
    dem = _utm_dem(tmp_path / 'elsewhere.tif', (131.5, 131.8), (30.15, 30.45))
    out = tmp_path / 'out'
    obj = cli.swath_main([velocity_file, *LINE, '-w', '5', '-topoxsect', dem, '--save',
                          '--nodisplay', '--outdir', str(out), '--noverbose'])
    sec = obj.sections[0]
    assert sec.topo_den is None and 'does not overlap' in sec.topo_error
    texts = ' '.join(t.get_text() for t in sec.ax_topo.texts)
    assert 'no topography for this profile' in texts and 'does not overlap' in texts
    info = next(out.glob('*_info.txt')).read_text()
    assert 'DEM sampling   : none - the profile does not overlap the DEM' in info


def test_projected_dem_that_mintpy_calls_4326(velocity_file, fine_dem, tmp_path,
                                              monkeypatch):
    """MintPy's GDAL reader gives some UTM GeoTIFFs EPSG=4326 and an ESRI zone
    name; the CRS must come from the file, so the DEM still lines up."""
    from viewinsar import raster, viewer

    real = raster.read_attr_any(fine_dem)
    wrong = dict(real, EPSG='4326', UTM_ZONE='WGS_1984_UTM_Zone_52N',
                 X_UNIT='meters', Y_UNIT='meters')
    monkeypatch.setattr(viewer, 'read_attr_any',
                        lambda f: wrong if f == fine_dem else raster.read_attr_any(f))
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-topoxsect', fine_dem,
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    sec = obj.sections[0]
    assert obj.topo_crs == UTM
    assert sec.topo_error is None
    assert _ridge_peak(sec) == pytest.approx(_fault_km(sec), abs=1.5)


def test_unitless_geotiff_labels(velocity_file, tmp_path):
    """A GeoTIFF has no unit, so view.py displays it in '1'. No '[1]' on the axes."""
    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, *LINE, '-w', '5', '--nodisplay',
                               '--outdir', str(tmp_path), '--noverbose'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.disp_unit = '1'                  # what view.py hands back for a GeoTIFF
    obj.plot()
    sec = obj.sections[0]
    assert sec.ax_swath.get_ylabel() == 'LOS velocity'
    # "x 0.25 cell]", never "x 0.25 1 cell]"
    assert re.search(r'x [0-9.e+-]+ cell\]$', sec.ax_cbar.get_ylabel())
    obj.close()

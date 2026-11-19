"""Tests of the 0.2.0 features: --nlook oversampling and --topoxsect."""

import numpy as np
import pytest

from conftest import FAULT_LON, LENGTH, WIDTH, make_velocity
from viewinsar import cli
from viewinsar.swath import bilinear_sample, swath_density, swath_extract

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')


##############################################################################
#                       bilinear sub-pixel sampling                          #
##############################################################################
def test_bilinear_is_exact_at_the_nodes():
    arr = np.arange(20, dtype=np.float64).reshape(4, 5)
    got = bilinear_sample(arr, np.arange(4.0), np.arange(5.0))
    assert np.allclose(got, arr)


def test_bilinear_is_exact_on_a_linear_ramp():
    """A bilinear scheme must reproduce a bilinear field with zero error."""
    ramp = (2 * np.arange(6)[:, None] + 3 * np.arange(7)[None, :]).astype(np.float64)
    rows = np.arange(0, 5.001, 0.25)
    cols = np.arange(0, 6.001, 0.25)
    want = 2 * rows[:, None] + 3 * cols[None, :]
    assert np.allclose(bilinear_sample(ramp, rows, cols), want, atol=1e-12)


def test_bilinear_renormalises_around_nan():
    """One bad neighbour must not void the sample - the weights rescale."""
    arr = np.array([[1.0, 1.0], [1.0, np.nan]])
    mid = bilinear_sample(arr, np.array([0.5]), np.array([0.5]))[0, 0]
    assert mid == pytest.approx(1.0)
    allnan = bilinear_sample(np.full((2, 2), np.nan), np.array([0.5]), np.array([0.5]))
    assert np.isnan(allnan[0, 0])


def test_bilinear_honours_a_mask():
    arr = np.array([[1.0, 5.0], [1.0, 1.0]])
    mask = np.array([[True, False], [True, True]])
    assert bilinear_sample(arr, np.array([0.5]), np.array([0.5]), mask=mask)[0, 0] \
        == pytest.approx(1.0)


##############################################################################
#                            --nlook                                         #
##############################################################################
def _extract(atr, nlook=1, half_width=6.0):
    data = make_velocity()
    return swath_extract(data, atr, (0, 0, WIDTH, LENGTH), [150, 60], [150, 390],
                         half_width, nlook=nlook)


def test_nlook_multiplies_the_sample_count(velocity_attrs):
    one = _extract(velocity_attrs, nlook=1)
    four = _extract(velocity_attrs, nlook=4)
    ratio = four['value'].size / one['value'].size
    assert 12 < ratio < 18, ratio        # about nlook**2 = 16
    assert four['nlook'] == 4


def test_nlook_keeps_the_corridor_geometry(velocity_attrs):
    """Oversampling must not leak outside the flat-capped buffer."""
    sw = _extract(velocity_attrs, nlook=5, half_width=7.0)
    assert np.abs(sw['perp']).max() <= 7.0 + 1e-6
    assert sw['dist'].min() >= -1e-6
    assert sw['dist'].max() <= sw['length'] + 1e-6


def test_nlook_does_not_change_the_profile(velocity_attrs):
    """The denser sampling must not bias the median - same signal, more points."""
    one, many = _extract(velocity_attrs, nlook=1), _extract(velocity_attrs, nlook=6)
    d1 = swath_density(one['dist'], one['value'], bin_x=1.0, bin_y=0.5, ylim=[-40, 40])
    d6 = swath_density(many['dist'], many['value'], bin_x=1.0, bin_y=0.5, ylim=[-40, 40])
    common = np.interp(d1['x_centers'], d6['x_centers'], d6['median_dw'])
    ok = np.isfinite(d1['median_dw']) & np.isfinite(common)
    assert ok.sum() > 10
    assert np.nanmean(common[ok] - d1['median_dw'][ok]) == pytest.approx(0, abs=0.5)


def test_nlook_rescues_a_multilooked_raster(coarse_velocity_file, velocity_file, tmp_path):
    """The point of the flag: a coarse raster profiles as densely as a fine one."""
    from viewinsar.viewer import swathProfileViewer

    line = ['--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75',
            '-w', '8', '--bin-x', '0.5', '--bin-y', '1.0', '--ylim', '-40', '40',
            '--noverbose', '--nodisplay', '--outdir', str(tmp_path)]

    def profile(fname, extra=()):
        inps = cli.cmd_line_parse([fname, *line, *extra])
        inps.save_fig = False
        obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
        obj.open()
        obj.plot()
        sec = obj.sections[0]
        den = sec.den
        out = (sec.swath_list[0]['value'].size,
               np.count_nonzero(den['counts']) / den['counts'].size, den)
        obj.close()
        return out

    n_fine = profile(velocity_file)[0]
    n_before, fill_before, _ = profile(coarse_velocity_file)
    n_after, fill_after, _ = profile(coarse_velocity_file, ['--nlook', '6'])

    assert n_before < n_fine / 20               # the sparse case, as reported
    assert n_after > 25 * n_before              # ~nlook**2 = 36x more samples
    assert fill_after > 1.3 * fill_before       # a visibly fuller density grid


def test_nlook_and_decimate_are_mutually_exclusive(velocity_file):
    with pytest.raises(ValueError, match='opposites'):
        cli.cmd_line_parse([velocity_file, '--nlook', '4', '--decimate', '2',
                            '--noverbose'])


def test_nlook_guards_against_an_absurd_sample_grid(velocity_attrs):
    with pytest.raises(ValueError, match='limit'):
        swath_extract(make_velocity(), velocity_attrs, (0, 0, WIDTH, LENGTH),
                      [150, 60], [150, 390], 8.0, nlook=200)


def test_auto_cmax_adapts_to_the_sampling(coarse_velocity_file, tmp_path):
    """With --cmax left alone the colour scale must follow the pixel density."""
    from viewinsar.viewer import swathProfileViewer

    base = ['--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75',
            '-w', '8', '--noverbose', '--nodisplay', '--outdir', str(tmp_path)]

    def cmax_for(extra):
        inps = cli.cmd_line_parse([coarse_velocity_file, *base, *extra])
        inps.save_fig = False
        obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
        obj.open()
        obj.plot()
        used = obj.sections[0].cmax_used
        obj.close()
        return used

    assert cmax_for(['--nlook', '6']) > 3 * cmax_for([])
    assert cmax_for(['--cmax', '42']) == 42.0          # an explicit value still wins


def test_the_density_colour_bar_says_pixel_density(coarse_velocity_file, tmp_path):
    """The colour bar of a swath / cross-section profile counts pixels."""
    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([coarse_velocity_file, '--start-lalo', '30.30', '130.15',
                               '--end-lalo', '30.30', '130.75', '-w', '8', '--noverbose',
                               '--nodisplay', '--outdir', str(tmp_path)])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    label = obj.sections[0].ax_cbar.get_ylabel()
    obj.close()
    assert label.startswith('Pixel density\n[count / ') and 'Point' not in label


##############################################################################
#                           --topoxsect                                      #
##############################################################################
def test_topo_panel_profiles_a_dem_in_another_crs(velocity_file, dem_file, tmp_path):
    """The DEM is UTM, the data lat/lon: the line must be re-projected."""
    from matplotlib.backend_bases import MouseButton, MouseEvent

    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, '-w', '8', '--topoxsect', dem_file,
                               '--outdir', str(tmp_path), '--noverbose', '--nodisplay'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    assert obj.topo_atr is not None, 'the DEM was rejected'
    obj.plot()
    obj.fig_map.canvas.draw()

    def click(lon, lat):
        x, y = obj.ax_img.transData.transform((lon, lat))
        obj.select_point(MouseEvent('button_release_event', obj.fig_map.canvas,
                                    x, y, button=MouseButton.LEFT))

    click(130.15, 30.30)
    click(130.75, 30.30)
    sec = obj.sections[0]
    assert sec.ax_topo is not None
    assert sec.topo_den is not None, 'no topography profile'

    # only ONE line, and it is black
    lines = sec.ax_topo.get_lines()
    assert len(lines) == 1
    assert lines[0].get_color() == 'k'

    # it must find the ridge the fixture puts over the fault at 130.45 E
    med = sec.topo_den['median_dw']
    good = np.isfinite(med)
    assert good.mean() > 0.8
    peak_km = sec.topo_den['x_centers'][np.nanargmax(med)]
    fault_km = (FAULT_LON - 130.15) * 111.32 * np.cos(np.deg2rad(30.3))
    assert peak_km == pytest.approx(fault_km, abs=4.0)
    assert np.nanmax(med) > 800 and np.nanmin(med[good]) < 400

    # the topography shares the velocity panel's distance axis
    assert sec.ax_topo.get_xlim() == pytest.approx(sec.ax_swath.get_xlim())
    # and it is written into the stats file
    sec.save()
    stats = np.genfromtxt(next(tmp_path.glob('*_A_*_stats.csv')), delimiter=',', names=True)
    assert 'elevation_median_m' in stats.dtype.names


def test_topo_rescales_onto_the_velocity_profile_length(velocity_file, dem_file, tmp_path):
    """Two grids give slightly different lengths; the panels must still register."""
    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, '--start-lalo', '30.30', '130.15',
                               '--end-lalo', '30.30', '130.75', '-w', '8',
                               '--topoxsect', dem_file, '--outdir', str(tmp_path),
                               '--noverbose', '--nodisplay'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    sec = obj.sections[0]
    assert sec.topo_swath['length'] == pytest.approx(sec.swath_list[0]['length'], rel=1e-9)


def test_bad_dem_is_reported_not_fatal(velocity_file, tmp_path, capsys):
    """A DEM that does not cover the profile must warn, not crash."""
    rasterio = pytest.importorskip('rasterio')
    from rasterio.transform import Affine

    from viewinsar.viewer import swathProfileViewer

    far = tmp_path / 'elsewhere.tif'
    with rasterio.open(far, 'w', driver='GTiff', height=20, width=20, count=1,
                       dtype='float32', crs='EPSG:4326',
                       transform=Affine.from_gdal(0.0, 0.01, 0, 0.0, 0, -0.01)) as dst:
        dst.write(np.ones((20, 20), dtype='float32'), 1)

    inps = cli.cmd_line_parse([velocity_file, '--start-lalo', '30.30', '130.15',
                               '--end-lalo', '30.30', '130.75', '-w', '8',
                               '--topoxsect', str(far), '--outdir', str(tmp_path),
                               '--noverbose', '--nodisplay'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    assert obj.sections[0].has_data               # the velocity profile still works
    assert 'WARNING' in capsys.readouterr().out


def test_missing_dem_is_caught_early(velocity_file):
    with pytest.raises(ValueError, match='no such file'):
        cli.cmd_line_parse([velocity_file, '--topoxsect', '/no/such/dem.tif',
                            '--noverbose'])


##############################################################################
#                      raster / CRS helpers                                  #
##############################################################################
def test_epsg_detection(velocity_attrs):
    from viewinsar.raster import epsg_of
    assert epsg_of(velocity_attrs) == 4326
    assert epsg_of({**velocity_attrs, 'UTM_ZONE': '52N', 'Y_UNIT': 'meters'}) == 32652
    assert epsg_of({**velocity_attrs, 'UTM_ZONE': '52S', 'Y_UNIT': 'meters'}) == 32752
    assert epsg_of({**velocity_attrs, 'EPSG': '3857'}) == 3857
    assert epsg_of({'LENGTH': '3', 'WIDTH': '3'}) is None          # radar coords


def test_line_round_trips_through_another_crs(velocity_attrs, dem_file):
    """Projecting the line onto the DEM grid and back must land where it started."""
    from viewinsar.raster import line_to_other_grid, read_attr_any

    dem_atr = read_attr_any(dem_file)
    start, end = [150, 75], [150, 375]
    there = line_to_other_grid(velocity_attrs, dem_atr, start, end, verbose=False)
    back = line_to_other_grid(dem_atr, velocity_attrs, there[0], there[1], verbose=False)
    for got, want in zip(back, (start, end)):
        assert got[0] == pytest.approx(want[0], abs=2)
        assert got[1] == pytest.approx(want[1], abs=2)


def test_read_any_matches_rasterio(dem_file):
    rasterio = pytest.importorskip('rasterio')
    from viewinsar.raster import read_any

    data = read_any(dem_file)[0]
    with rasterio.open(dem_file) as src:
        assert np.allclose(data, src.read(1), equal_nan=True)
        assert data.shape == src.shape
    sub = read_any(dem_file, box=(10, 20, 60, 70))[0]
    assert sub.shape == (50, 50)
    assert np.allclose(sub, data[20:70, 10:60], equal_nan=True)

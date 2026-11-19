"""`viewinsar stat`: the creep rate along a line, as Hussain et al. (2016).

The geometry is checked on a metric grid (100 m pixels): the points along the
line, the sides, the profile windows across it and the pixels in them. Then
the two-line fit of their eq. 2 against separate fits, and the creep of a
known step recovered along a line. Then the whole command on synthetic
products: a creep step on the N-S fault at 130.45 E of the fixtures, the
radar-coded product, and an ascending and a descending track over an E-W
fault, split into fault-parallel and vertical.
"""

import json
import os

import numpy as np
import pytest

from conftest import DLAT, DLON, FAULT_LON, KM_PER_DEG_LON, LAT0, LON0, _base_attrs
from viewinsar import stat
from viewinsar.model import fault_offset, fit_two_lines, two_line_matrix

METRIC = {'Y_FIRST': '0', 'X_FIRST': '0', 'Y_STEP': '-100', 'X_STEP': '100',
          'Y_UNIT': 'meters', 'X_UNIT': 'meters', 'LENGTH': '800', 'WIDTH': '800'}
BOX = (0, 0, 800, 800)


def frame_km(points_km, atr=METRIC):
    """A frame from vertices given in km (east, north) on the metric grid."""
    yx = yx_km(points_km)
    return stat.polyline_frame(atr, yx)


def yx_km(points_km):
    """[row, col] of points given in km (east, north) on the metric grid."""
    xy = np.asarray(points_km, dtype=float) * 1000.0
    return np.c_[-xy[:, 1] / 100.0, xy[:, 0] / 100.0]


def metric_field(fun):
    """fun(east_km, north_km) on the metric grid."""
    rows, cols = np.mgrid[0:800, 0:800].astype(float)
    return fun(cols * 0.1, -rows * 0.1)


##############################################################################
#                          points along the line                             #
##############################################################################
def test_points_are_spaced_and_centred():
    frame = frame_km([[0, -40], [23.4, -40]])
    dist, points, normals = stat.sample_points(frame, 5.0)
    assert dist == pytest.approx([1.7, 6.7, 11.7, 16.7, 21.7])
    assert np.diff(points[:, 0]) == pytest.approx([5000.0] * 4)
    # walking east, the left is north
    assert normals == pytest.approx(np.tile([0.0, 1.0], (5, 1)))
    # a line shorter than the spacing still gets its middle point
    dist = stat.sample_points(frame_km([[0, -40], [3, -40]]), 5.0)[0]
    assert dist == pytest.approx([1.5])


def test_a_point_on_a_vertex_takes_the_bisector():
    frame = frame_km([[0, -40], [10, -40], [10, -30]])        # a left turn of 90 deg
    dist, points, normals = stat.sample_points(frame, 10.0)
    assert dist == pytest.approx([0.0, 10.0, 20.0])
    assert points[1] == pytest.approx([10000.0, -40000.0])
    assert normals[1] == pytest.approx(np.array([-1.0, 1.0]) / np.sqrt(2))


def test_normals_follow_the_line_not_its_wiggles():
    rng = np.random.default_rng(3)
    east = np.arange(0, 40.01, 0.1)
    north = -40 + rng.normal(0, 0.03, east.size)              # +/- 30 m of digitising noise
    frame = frame_km(np.c_[east, north])
    raw = stat.sample_points(frame, 5.0)[2]
    smooth = stat.sample_points(frame, 5.0, smooth_km=2.0)[2]
    angle = lambda n: np.degrees(np.abs(np.arctan2(n[:, 0], n[:, 1])))  # noqa: E731
    assert angle(smooth).max() < 3.0
    assert angle(raw).max() > angle(smooth).max()


##############################################################################
#                                sides                                       #
##############################################################################
def test_sides_on_a_straight_line():
    frame = frame_km([[0, -40], [40, -40]])                   # walking east
    xy = np.array([[20, -35], [20, -45], [-5, -39], [45, -41]]) * 1000.0
    assert list(stat.side_of_polyline(frame, xy)) == [1, -1, 1, -1]
    assert stat.distance_to_polyline(frame, xy) == pytest.approx([5000, 5000, 5099.0, 5099.0],
                                                                 rel=1e-3)


def test_side_just_outside_a_sharp_bend():
    """Nearest to the vertex of a 150 deg left turn: OUTSIDE, so on the right."""
    turn = np.radians(150)
    a, v = np.array([0.0, -40.0]), np.array([10.0, -40.0])
    b = v + 10 * np.array([np.cos(turn), np.sin(turn)])
    frame = frame_km([a, v, b])
    q = (v + [1.0, 0.5]) * 1000.0
    assert stat.side_of_polyline(frame, q)[0] == -1
    inside = (v + 2 * np.array([np.cos(np.radians(165)), np.sin(np.radians(165))])) * 1000.0
    assert stat.side_of_polyline(frame, inside)[0] == 1


##############################################################################
#                       the profiles across the line                         #
##############################################################################
def test_a_profile_window_and_its_pixels():
    """5 km either side of an E-W line, 2.5 km either side of the profile."""
    frame = frame_km([[0, -40], [40, -40]])                   # walking east: row 400
    dist, points, normals = stat.sample_points(frame, 5.0)
    across, along = stat.profile_axes(normals[3])
    assert across == pytest.approx([0, -1]) and along == pytest.approx([1, 0])
    data = np.ones((800, 800), dtype=np.float32)
    out = stat.profile_pixels(data, BOX, frame, points[3], normals[3], 5.0, 2.5)
    # rows 350 .. 450 but the line's own (x = 0), columns +/- 25 around the point
    assert out['values'].size == 100 * 51
    assert not out['clipped'] and out['dropped'] == 0
    assert np.abs(out['x_km']).max() == pytest.approx(5.0)
    assert np.abs(out['y_km']).max() <= 2.5 + 1e-9
    assert out['x_km'][out['rows'] > 400].min() > 0          # south (block 2) is x > 0
    assert out['x_km'][out['rows'] < 400].max() < 0
    # the outline: 10 km across (north - south), 5 km along (east - west)
    rows, cols = stat.profile_outline(frame, points[3], normals[3], 5.0, 2.5)
    assert np.ptp(rows) == pytest.approx(100.0) and np.ptp(cols) == pytest.approx(50.0)
    rows, _ = stat.profile_outline(frame, points[3], normals[3], 5.0, 2.5, 0.0, 5.0)
    assert rows.min() == pytest.approx(400.0) and rows.max() == pytest.approx(450.0)

    # NaN, masked, zero: left out (zero only by default)
    r, c = 380, int(round(stat.to_yx(frame, points[3])[1]))
    data[r, c], data[r + 1, c], data[r + 2, c] = np.nan, 0.0, 5.0
    mask = np.ones_like(data, dtype=bool)
    mask[r + 3, c] = False
    n_all = out['values'].size
    got = stat.profile_pixels(data, BOX, frame, points[3], normals[3], 5.0, 2.5, mask=mask)
    assert got['values'].size == n_all - 3
    got = stat.profile_pixels(data, BOX, frame, points[3], normals[3], 5.0, 2.5, mask=mask,
                              keep_zero=True)
    assert got['values'].size == n_all - 2
    # a window cut by the data read: only the part inside, and it says so
    sub = data[:420, :]
    got = stat.profile_pixels(sub, (0, 0, 800, 420), frame, points[3], normals[3], 5.0, 2.5)
    assert got['clipped'] and got['rows'].max() == 419
    # no point (outside the other file): nothing
    assert stat.profile_pixels(data, BOX, frame, [np.nan, np.nan], normals[3], 5.0,
                               2.5)['values'].size == 0


def test_pixels_beyond_a_bend_are_left_out():
    """A right angle inside the window: the pixels across the line from where
    their x says are dropped, every pixel kept is on its own side."""
    frame = frame_km([[0, -40], [20, -40], [20, -60]])        # east, then south
    point = np.array([18000.0, -40000.0])                     # 2 km before the corner
    normal = np.array([0.0, 1.0])
    data = np.ones((800, 800), dtype=np.float32)
    out = stat.profile_pixels(data, BOX, frame, point, normal, 5.0, 2.5)
    assert out['dropped'] > 0
    assert stat.trace_offset(frame, point, normal, 5.0, 2.5) == pytest.approx(5.0)
    xy = stat.from_yx(frame, out['rows'], out['cols'])
    sides = stat.side_of_polyline(frame, xy)
    assert np.all(sides == np.where(out['x_km'] < 0, 1, -1))
    # a straight line: nothing to drop, the line at the centre
    straight = frame_km([[0, -40], [40, -40]])
    assert stat.profile_pixels(data, BOX, straight, point, normal, 5.0, 2.5)['dropped'] == 0
    assert stat.trace_offset(straight, point, normal, 5.0, 2.5) == pytest.approx(0.0, abs=1e-9)


##############################################################################
#                the two lines of Hussain et al. (2016, eq. 2)               #
##############################################################################
def test_the_matrix_of_eq2():
    mat = two_line_matrix([-3.0, -1.0], [2.0, 4.0, 5.0])
    assert mat.tolist() == [[-3, 0, 1, 0], [-1, 0, 1, 0], [0, 2, 1, 1], [0, 4, 1, 1],
                            [0, 5, 1, 1]]


def test_two_lines_as_one_system_or_two_fits():
    rng = np.random.default_rng(8)
    xa, xb = rng.uniform(-5, 0, 400), rng.uniform(0, 5, 300)
    va = 1.0 + 0.4 * xa + rng.normal(0, 1, xa.size)
    vb = 1.0 + 3.0 - 0.2 * xb + rng.normal(0, 1, xb.size)
    fit = fit_two_lines(xa, va, xb, vb, n_boot=1500, seed=4)
    slope_a, icpt_a = np.polyfit(xa, va, 1)
    slope_b, icpt_b = np.polyfit(xb, vb, 1)
    assert fit['params'] == pytest.approx([slope_a, slope_b, icpt_a, icpt_b - icpt_a])
    assert fit['offset'] == pytest.approx(icpt_b - icpt_a)
    assert fit['line_a'] == pytest.approx((slope_a, icpt_a))
    assert fit['line_b'] == pytest.approx((slope_b, icpt_b))
    assert fit['offset'] == pytest.approx(3.0, abs=4 * fit['std'])
    # the BLUE's sigma: the scatter about both lines, through (G^T G)^-1
    mat = two_line_matrix(xa, xb)
    resid = np.r_[va, vb] - mat @ fit['params']
    var = resid @ resid / (resid.size - 4)
    assert fit['offset_std'] == pytest.approx(np.sqrt(var * np.linalg.inv(mat.T @ mat)[3, 3]))
    # the bootstrap: every parameter per draw, reproducible, as model --creep-offset's
    assert fit['boot'].shape == (1500, 4)
    again = fit_two_lines(xa, va, xb, vb, n_boot=1500, seed=4)
    assert np.array_equal(fit['boot'], again['boot'])
    assert np.array_equal(fit['boot'][:, 3],
                          fault_offset(xa, va, xb, vb, n_boot=1500, seed=4)['boot'])
    assert fit['lo'] < fit['lo68'] < fit['median'] < fit['hi68'] < fit['hi']
    assert fit['std'] == pytest.approx(fit['offset_std'], rel=0.25)   # independent errors
    # too few points on a side: nothing
    assert fit_two_lines(xa[:2], va[:2], xb, vb) is None


def test_creep_along_a_line_recovers_a_step():
    """A step that grows along the fault, on a tilted background."""
    def field(east, north):
        x = -(north + 40.0)                                   # + towards block 2 (south)
        creep = 2.0 + 0.1 * east
        return 0.3 * x + np.where(x > 0, creep / 2, -creep / 2)

    rng = np.random.default_rng(2)
    data = metric_field(field) + rng.normal(0, 0.5, (800, 800))
    line = yx_km([[2, -40], [57, -40]])
    res = stat.creep_along_line(data, METRIC, BOX, line, 10.0, n_boot=500)
    assert res['length_km'] == pytest.approx(55.0)
    assert res['dist_km'] == pytest.approx([2.5, 12.5, 22.5, 32.5, 42.5, 52.5])
    for dist, prof in zip(res['dist_km'], res['profiles']):
        fit = prof['fit']
        want = 2.0 + 0.1 * (2.0 + dist)
        assert fit['mean'] == pytest.approx(want, abs=max(4 * fit['std'], 0.05))
        assert fit['offset'] == pytest.approx(want, abs=max(4 * fit['offset_std'], 0.05))
        assert fit['line_a'][0] == pytest.approx(0.3, abs=0.03)       # the slopes
        assert prof['n_1'] == prof['n_2'] == 50 * 51
    # the line drawn the other way: the blocks swap, the creep changes sign
    back = stat.creep_along_line(data, METRIC, BOX, line[::-1], 10.0, n_boot=200)
    blue = [p['fit']['offset'] for p in res['profiles']]
    assert [p['fit']['offset'] for p in back['profiles']] == pytest.approx(-np.array(blue[::-1]))
    # every point its own bootstrap: reproducible, and not the same draws
    again = stat.creep_along_line(data, METRIC, BOX, line, 10.0, n_boot=500)
    assert np.array_equal(again['profiles'][2]['fit']['boot'], res['profiles'][2]['fit']['boot'])
    assert res['profiles'][1]['fit']['std'] != res['profiles'][2]['fit']['std']
    # too few pixels: no creep there
    sparse = np.where(rng.random((800, 800)) < 0.0005, data, np.nan)
    few = stat.creep_along_line(sparse, METRIC, BOX, line, 10.0, n_boot=100, min_points=10)
    assert all(p['fit'] is None for p in few['profiles'])
    with pytest.raises(ValueError):
        stat.creep_along_line(data, METRIC, BOX, line, 10.0, width_km=0)


##############################################################################
#                       left / right on a mirrored grid                      #
##############################################################################
def test_mirrored_grids():
    assert not stat.grid_mirrored({'Y_FIRST': '30', 'X_STEP': '0.001', 'Y_STEP': '-0.001'})
    assert stat.grid_mirrored({'Y_FIRST': '30', 'X_STEP': '0.001', 'Y_STEP': '0.001'})
    # radar: right-looking grids are mirror images of the ground; left-looking not
    assert stat.grid_mirrored({'LENGTH': '10'})
    assert not stat.grid_mirrored({'ANTENNA_SIDE': '1'})

    class Lookup:                       # descending: rows run south, columns west
        rows, cols = np.mgrid[0:50, 0:60].astype(float)
        lat, lon = 30.6 - 0.0018 * rows, 130.9 - 0.0025 * cols
    assert stat.grid_mirrored({}, Lookup)
    Lookup.lon = 130.9 + 0.0025 * Lookup.cols          # a grid that is NOT mirrored
    assert not stat.grid_mirrored({}, Lookup)

    # the same line on a mirrored grid: block 1 is still the left of the GROUND
    frame = stat.polyline_frame(dict(METRIC), [[400, 100], [400, 700]], mirrored=True)
    # columns now run west: walking along +col is walking WEST, the left is south
    _, points, normals = stat.sample_points(frame, 10)
    assert normals[0] == pytest.approx([0, -1])
    out = stat.profile_pixels(np.ones((800, 800)), BOX, frame, points[0], normals[0], 5, 2.5)
    assert np.all(out['rows'][out['x_km'] < 0] > 400)        # block 1: south = larger row
    assert stat.to_yx(frame, stat.from_yx(frame, 12.0, 34.0)) == pytest.approx((12.0, 34.0))


##############################################################################
#                           the command line                                 #
##############################################################################
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from viewinsar import cli  # noqa: E402

UP_THE_FAULT = ['--line-lalo', '30.05', str(FAULT_LON), '30.52', str(FAULT_LON)]
CREEP_MM = 4.0                                   # the creep step of creep_file, mm/yr


def _write_h5(path, data, **attrs):
    h5py = pytest.importorskip('h5py')
    with h5py.File(path, 'w') as fid:
        fid.create_dataset('velocity', data=np.asarray(data, dtype=np.float32))
        for key, val in {**_base_attrs(), **attrs}.items():
            fid.attrs[key] = val
    return str(path)


def lonlat_grid():
    lon = LON0 + (np.arange(450) + 0.5) * DLON
    lat = LAT0 + (np.arange(300) + 0.5) * DLAT
    return np.meshgrid(lon, lat)


@pytest.fixture(scope='module')
def creep_file(tmp_path_factory):
    """The N-S fault of the fixtures, creeping: a step of CREEP_MM (east side
    up the LOS) on a linear gradient, in m/year."""
    lon2d, _ = lonlat_grid()
    km_east = (lon2d - FAULT_LON) * KM_PER_DEG_LON
    rng = np.random.default_rng(11)
    vel = (0.4 * km_east + np.where(km_east > 0, CREEP_MM / 2, -CREEP_MM / 2)
           + rng.normal(0, 1.0, km_east.shape)) / 1000.0
    vel[:20, :] = np.nan                                      # a decorrelated strip
    return _write_h5(tmp_path_factory.mktemp('creep') / 'geo_velocity.h5', vel)


def _read_csv(path):
    return np.genfromtxt(path, delimiter=',', names=True, dtype=None, encoding='utf-8')


def test_stat_options(velocity_file):
    inps = cli.stat_cmd_line_parse([velocity_file, '--noverbose', *UP_THE_FAULT])
    assert (inps.spacing, inps.length_km, inps.width_km) == (5.0, 5.0, 2.5)
    assert (inps.n_boot, inps.seed, inps.min_pixels) == (2000, 0, 10)
    assert inps.line_lalo == [[30.05, FAULT_LON], [30.52, FAULT_LON]]
    assert inps.block_names == ['Block 1', 'Block 2']
    inps = cli.stat_cmd_line_parse([velocity_file, '--noverbose', '-sp', '2.5', '--length', '4',
                                    '--width', '1.5', '--boot', '500', '--seed', '7',
                                    '--min-pixels', '20', '--line-yx', '10', '20', '30',
                                    '40', '50', '60', '-c', 'RdBu_r', '--names', 'W', 'E'])
    assert (inps.spacing, inps.length_km, inps.width_km) == (2.5, 4.0, 1.5)
    assert (inps.n_boot, inps.seed, inps.min_pixels) == (500, 7, 20)
    assert inps.line_yx == [[10, 20], [30, 40], [50, 60]]
    assert inps.block_names == ['W', 'E']
    assert inps.view_extra == ['-c', 'RdBu_r']               # the rest goes to view.py

    bad = [['-sp', '0'], ['--length', '-1'], ['--width', '0'], ['--boot', '50'],
           ['--min-pixels', '2'], ['--line-yx', '1', '2', '3'], ['--line-yx', '1', '2'],
           ['--colors', 'red', 'notacolour'], ['--line-file', 'nowhere.shp'], ['--nodisplay'],
           ['--line-yx', '1', '2', '3', '4', '--line-lalo', '30', '130', '30.1', '130'],
           ['--los-inc', '30'], ['--los-inc', '30', '31', '--los-az', '100', '-100'],
           ['--geom', 'nowhere.h5']]
    for extra in bad:
        with pytest.raises(ValueError):
            cli.stat_cmd_line_parse([velocity_file, '--noverbose', *extra])


def test_the_circles_and_the_block_rate_are_gone(velocity_file, capsys):
    """0.5.8: --radius, --gap, --rate and --rate-err (and -r, -rate) are left out,
    with their values and a note, so the commands of before still run."""
    argv = [velocity_file, '--radius', '2', '--gap=0.5', '-rate', '--rate-err', 'sem', '-r',
            '3', '--rate', '--spacing', '5', '-v', '-1', '1']
    assert cli.modern_options(argv, 'stat') == [velocity_file, '--spacing', '5', '-v', '-1',
                                                 '1']
    out = capsys.readouterr().out
    for old in ('--radius', '--gap', '-rate', '--rate-err', '-r', '--rate'):
        assert f'NOTE: {old} is not needed any more' in out
    assert '--length, --width' in out and 'bootstrap' in out
    known = set(cli.create_stat_parser()._option_string_actions)
    assert not known & {'--radius', '--gap', '--rate', '--rate-err', '-r'}
    inps = cli.stat_cmd_line_parse([velocity_file, '--noverbose', '--radius', '2', '--gap',
                                    '1', *UP_THE_FAULT])
    assert inps.view_extra == []


def test_stat_end_to_end(creep_file, tmp_path):
    obj = cli.stat_main([creep_file, '-sp', '10', *UP_THE_FAULT, '--nodisplay', '--noverbose',
                         '--outdir', str(tmp_path)])
    names = sorted(os.listdir(tmp_path))
    base = os.path.basename(obj.outfile_base)
    assert names == sorted(f'{base}{tail}' for tail in (
        '.png', '_map.png', '_creep.csv', '_values.csv', '_info.txt', '_profiles.png'))
    assert base.endswith('_sp10km_L5km_w2.5km')
    assert obj.rate_unit == 'mm/yr' and obj.factor == 10.0      # view.py's cm/year -> mm/yr

    cr = _read_csv(tmp_path / f'{base}_creep.csv')
    assert len(cr) == 6
    assert np.all(cr['label'] == 'geo_velocity') and np.all(cr['orbit'] == 'asc')
    # block 2 (east, the right of a line drawn north) - block 1, in mm/yr
    assert cr['creep'] == pytest.approx(CREEP_MM, abs=0.6)
    assert np.all(np.abs(cr['creep'] - CREEP_MM) < 4 * cr['creep_std'] + 0.05)
    assert cr['creep_blue'] == pytest.approx(cr['creep'], abs=0.1)
    assert np.all((cr['creep_lo95'] < cr['creep']) & (cr['creep'] < cr['creep_hi95']))
    assert np.all((cr['creep_lo95'] < cr['creep_lo68']) & (cr['creep_hi68'] < cr['creep_hi95']))
    assert cr['slope_1'] == pytest.approx(0.4, abs=0.08)     # mm/yr per km, east is +
    assert cr['slope_2'] == pytest.approx(0.4, abs=0.08)
    assert np.all(cr['n_1'] == 26 * 23) and np.all(cr['n_2'] == 26 * 23)   # 5 x 5 km a side
    assert cr['lon'] == pytest.approx(FAULT_LON, abs=0.002)
    assert cr['strike_deg'] == pytest.approx(0.0, abs=0.5)   # the line runs north
    assert cr['trace_km'] == pytest.approx(0.0, abs=1e-6) and np.all(cr['dropped'] == 0)
    assert np.diff(cr['dist_km']) == pytest.approx(10.0)

    val = _read_csv(tmp_path / f'{base}_values.csv')
    assert len(val) == cr['n_1'].sum() + cr['n_2'].sum()
    assert np.all(val['lon'][val['block'] == 1] < FAULT_LON)
    assert np.all(val['lon'][val['block'] == 2] > FAULT_LON)
    assert np.abs(val['x_km']).max() <= 5.0 + 1e-6 and np.abs(val['y_km']).max() <= 2.5 + 1e-6
    # x is the distance east of the line (on the pixel centres of the fault's), across
    km_east = (val['lon'] - cr['lon'][0]) * KM_PER_DEG_LON
    assert val['x_km'] == pytest.approx(km_east, abs=0.01)

    info = (tmp_path / f'{base}_info.txt').read_text()
    assert 'LEFT of the line' in info and 'Hussain et al. (2016), eq. 2' in info
    assert 'creep rate = Block 2 - Block 1' in info

    # the plot: one series, its points the creep rates of the table
    assert len(obj.axes_stat) == 1
    line = [ln for ln in obj.ax_stat.get_lines() if ln.get_gid() == 'creep_1'][0]
    assert np.asarray(line.get_ydata()) == pytest.approx(cr['creep'], abs=1e-6)
    assert line.get_linestyle() == '--'                       # an ascending track
    assert [t.get_text() for t in obj.ax_stat.get_legend().get_texts()] == ['geo_velocity']
    assert obj.ax_stat.get_ylabel() == 'LOS creep rate (mm/yr)'
    # the map: each profile as two halves, the points as black circles, the key
    halves = [p for p in obj.ax_img.patches if str(p.get_gid()).startswith('profile_')]
    assert len(halves) == 12
    assert halves[0].get_edgecolor()[:3] != pytest.approx(halves[1].get_edgecolor()[:3])
    points = [ln for ln in obj.ax_img.get_lines() if ln.get_gid() == 'points'][0]
    assert len(points.get_xdata()) == 6
    keys = [a for a in obj.ax_img.get_children() if isinstance(a, matplotlib.legend.Legend)]
    assert [t.get_text() for t in keys[0].get_texts()] == ['Block 1', 'Block 2',
                                                           'Profile point']


def test_stat_interactive(creep_file, tmp_path):
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    from viewinsar.statview import StatViewer

    inps = cli.stat_cmd_line_parse([creep_file, '-sp', '10', '--boot', '300', '--noverbose',
                                    '--outdir', str(tmp_path)])
    inps.disp_fig = False                     # draw on Agg, without plt.show()
    obj = StatViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    obj.fig_map.canvas.draw()

    def click(lon, lat, button=MouseButton.LEFT):
        x, y = obj.ax_img.transData.transform((lon, lat))
        obj.select_point(MouseEvent('button_release_event', obj.fig_map.canvas, x, y,
                                    button=button))

    def key(name, fig=None):
        event = KeyEvent('key_press_event', (fig or obj.fig_map).canvas, key=name)
        (obj.on_profile_key if fig is not None and fig is obj.fig_prof else obj.on_key)(event)

    click(FAULT_LON, 30.05)
    click(FAULT_LON, 30.30)                    # on the fault: x = 0 is where it creeps
    click(FAULT_LON + 0.05, 30.45)
    assert len(obj.vertices) == 3
    key('u')                                   # undo the last point
    assert len(obj.vertices) == 2
    click(FAULT_LON, 30.52)
    click(FAULT_LON, 30.52, button=MouseButton.RIGHT)          # ignored
    assert len(obj.vertices) == 3
    key('enter')
    assert obj.finished and obj.result is not None
    creep, std, _ = obj.creep_of(0)
    assert np.nanmedian(creep) == pytest.approx(CREEP_MM, abs=1.0)

    # a click on a point of the creep rate: its profile, in its own window
    obj.fig_stat.canvas.draw()
    ax = obj.axes_stat[0]
    x, y = ax.transData.transform((obj.line['dist_km'][2], creep[2]))
    obj.on_stat_click(MouseEvent('button_release_event', obj.fig_stat.canvas, x + 3, y - 3,
                                 button=MouseButton.LEFT))
    assert obj.selected == (0, 2) and obj.fig_prof is not None
    ax1, ax2 = obj.fig_prof.axes[:2]
    pixels = [ln for ln in ax1.get_lines() if ln.get_gid() == 'pixels'][0]
    prof = obj.result['profiles'][2]
    assert len(pixels.get_xdata()) == prof['n_1'] + prof['n_2']
    labels = [t.get_text() for t in ax1.get_legend().get_texts()]
    assert labels == ['BLUE', 'Bootstrap']
    labels = [t.get_text() for t in ax2.get_legend().get_texts()]
    assert labels[0] == 'Bootstrap draws' and labels[1].startswith('Bootstrap: ')
    assert labels[2].startswith('BLUE: ') and '95%' in labels[1]
    assert ax2.get_ylabel() == 'Normalised PDF'
    assert any(ln.get_gid() == 'selected' for ln in obj.ax_img.get_lines())  # on the map
    # far from any point: nothing
    obj.on_stat_click(MouseEvent('button_release_event', obj.fig_stat.canvas, x, y + 200,
                                 button=MouseButton.LEFT))
    assert obj.selected == (0, 2)
    # the arrows: the point after / before (and never past the ends)
    key('right', obj.fig_prof)
    assert obj.selected == (0, 3)
    for _ in range(5):
        key('left', obj.fig_prof)
    assert obj.selected == (0, 0)

    key('e')
    written = {p.name for p in tmp_path.iterdir()}
    assert any(n.endswith('_creep.csv') for n in written)
    assert any(n.endswith('_profile1.png') for n in written)      # the one on screen

    # a click after Enter starts a new line
    click(FAULT_LON, 30.10)
    assert not obj.finished and len(obj.vertices) == 1
    key('n')
    assert obj.vertices == [] and obj.result is None and obj.selected is None
    obj.close()


def test_stat_line_file_takes_the_longest_line(creep_file, tmp_path):
    obj = {'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': {},
         'geometry': {'type': 'LineString', 'coordinates': [[130.2, 30.2], [130.25, 30.22]]}},
        {'type': 'Feature', 'properties': {},
         'geometry': {'type': 'LineString',
                      'coordinates': [[FAULT_LON, 30.05], [FAULT_LON, 30.30],
                                      [FAULT_LON, 30.52]]}}]}
    path = tmp_path / 'trace.geojson'
    path.write_text(json.dumps(obj))
    viewer = cli.stat_main([creep_file, '-sp', '10', '--boot', '300', '--line-file', str(path),
                            '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert viewer.line['length_km'] == pytest.approx(0.47 * 111.19, rel=0.01)
    assert np.nanmedian(viewer.creep_of(0)[0]) == pytest.approx(CREEP_MM, abs=1.0)


def test_stat_on_radar_coded_data(tmp_path):
    """Descending: the grid is a mirror image of the ground; west must stay
    block 1 - and the two lines through an arctan are what its pixels give."""
    from test_radar import make_radar_product

    vfile, _ = make_radar_product(str(tmp_path / 'proc'))
    obj = cli.stat_main([vfile, '-sp', '5', '--boot', '300', '--line-lalo', '30.30',
                         str(FAULT_LON), '30.50', str(FAULT_LON), '--nodisplay',
                         '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert obj.mirrored is True
    cr = _read_csv(next((tmp_path / 'out').glob('*_creep.csv')))
    val = _read_csv(next((tmp_path / 'out').glob('*_values.csv')))
    # (a pixel centre ON the fault - within a metre - may go either way)
    assert np.all(val['lon'][val['block'] == 1] < FAULT_LON + 1e-5)
    assert np.all(val['lon'][val['block'] == 2] > FAULT_LON - 1e-5)
    assert np.all(cr['slope_1'] > 0) and np.all(cr['slope_2'] > 0)   # faster to the east
    # the creep of the noise-free arctan at the same pixels
    for point in cr['point']:
        one = val[val['point'] == point]
        km_east = (one['lon'] - FAULT_LON) * 111.32 * np.cos(np.radians(30.4))
        true = 300 * (2 / np.pi) * np.arctan(km_east / 8.0)          # mm/yr
        a, b = one['block'] == 1, one['block'] == 2
        want = fit_two_lines(one['x_km'][a], true[a], one['x_km'][b], true[b], n_boot=10)
        row = cr[cr['point'] == point]
        assert row['creep'][0] == pytest.approx(want['offset'], abs=4 * row['creep_std'][0] + 0.3)


##############################################################################
#                ascending and descending: the split                         #
##############################################################################
FAULT_LAT = 30.30                     # an E-W fault on the pixel edge of rows 149 | 150
PAR_MM, UP_MM = 6.0, 2.0              # right-lateral creep; the south side up


def los_vector(inc, az):
    i, a = np.radians(inc), np.radians(az)
    return -np.sin(i) * np.sin(a), np.sin(i) * np.cos(a), np.cos(i)


@pytest.fixture(scope='module')
def asc_desc(tmp_path_factory):
    """A descending and an ascending track over an E-W fault that creeps: the
    north side east by PAR_MM (right-lateral), the south side up by UP_MM."""
    folder = tmp_path_factory.mktemp('asc_desc')
    lon2d, lat2d = lonlat_grid()
    north = lat2d > FAULT_LAT
    ve = np.where(north, PAR_MM / 2, -PAR_MM / 2) + 0.2 * (lat2d - FAULT_LAT) * 111.19
    vu = np.where(north, -UP_MM / 2, UP_MM / 2)
    rng = np.random.default_rng(21)
    files = []
    for name, inc, az, orbit in (('T036D', 39.0, -100.0, 'DESCENDING'),
                                 ('T131A', 34.0, 100.0, 'ASCENDING')):
        e, n, u = los_vector(inc, az)
        vel = (ve * e + vu * u + rng.normal(0, 0.8, ve.shape)) / 1000.0
        files.append(_write_h5(folder / f'{name}.h5', vel, ORBIT_DIRECTION=orbit))
    return files


ALONG_THE_FAULT = ['--line-lalo', str(FAULT_LAT), '130.10', str(FAULT_LAT), '130.80']


def test_stat_split_into_fault_parallel_and_vertical(asc_desc, tmp_path):
    obj = cli.stat_main([*asc_desc, '--labels', 'T036D', 'T131A', '--los-inc', '39', '34',
                         '--los-az', '-100', '100', '-sp', '10', '--boot', '500',
                         *ALONG_THE_FAULT, '--names', 'North', 'South', '--nodisplay',
                         '--noverbose', '--outdir', str(tmp_path)])
    assert obj.orbits == ['desc', 'asc'] and obj.split_files == [0, 1]
    base = os.path.basename(obj.outfile_base)
    names = set(os.listdir(tmp_path))
    for tail in ('_creep.csv', '_split.csv', '_values.csv', '_info.txt', '.png', '_map.png',
                 '_profiles.png', '_profiles_2.png'):
        assert f'{base}{tail}' in names, tail

    # the LOS creep of each track: block 2 (south) - block 1 (north)
    cr = _read_csv(tmp_path / f'{base}_creep.csv')
    for label, inc, az in (('T036D', 39.0, -100.0), ('T131A', 34.0, 100.0)):
        e, _, u = los_vector(inc, az)
        rows = cr[cr['label'] == label]
        want = -PAR_MM * e + UP_MM * u
        assert rows['creep'] == pytest.approx(want, abs=0.5)
        assert rows['los_east'] == pytest.approx(e, abs=1e-6)
        assert rows['los_up'] == pytest.approx(u, abs=1e-6)
    assert cr['strike_deg'] == pytest.approx(90.0, abs=0.5)  # walking east

    sp = _read_csv(tmp_path / f'{base}_split.csv')
    assert len(sp) == len(cr) // 2
    assert sp['par_mean'] == pytest.approx(PAR_MM, abs=0.6)      # + right-lateral
    assert sp['up_mean'] == pytest.approx(UP_MM, abs=0.4)        # + South up
    assert sp['par'] == pytest.approx(sp['par_mean'], abs=0.15)
    assert np.all(sp['par_std'] > 0) and np.all(sp['files'] == 'T036D+T131A')
    info = (tmp_path / f'{base}_info.txt').read_text()
    assert 'fault-parallel creep (+ right-lateral)' in info and 'South up' in info

    # two panels: the LOS of each track, and the split in red and blue
    assert len(obj.axes_stat) == 2
    lines = {ln.get_gid(): ln for ax in obj.axes_stat for ln in ax.get_lines() if ln.get_gid()}
    assert lines['creep_1'].get_linestyle() == '-' and lines['creep_2'].get_linestyle() == '--'
    assert np.asarray(lines['split_par'].get_ydata()) == pytest.approx(sp['par_mean'])
    assert matplotlib.colors.to_hex(lines['split_par'].get_color()) == '#d62728'
    assert matplotlib.colors.to_hex(lines['split_up'].get_color()) == '#1f77b4'
    texts = [t.get_text() for t in obj.axes_stat[1].get_legend().get_texts()]
    assert texts == ['Fault-parallel (+ right-lateral)', 'Vertical (+ South up)']
    texts = [t.get_text() for t in obj.axes_stat[0].get_legend().get_texts()]
    assert texts == ['T036D', 'T131A', 'descending', 'ascending']

    # the line drawn the other way: the same right-lateral creep, the other side up
    back = cli.stat_main([*asc_desc, '--los-inc', '39', '34', '--los-az', '-100', '100',
                          '-sp', '10', '--boot', '300', '--line-lalo', str(FAULT_LAT), '130.80',
                          str(FAULT_LAT), '130.10', '--nodisplay', '--noverbose',
                          '--outdir', str(tmp_path / 'back')])
    assert np.nanmedian(back.split['par_mean']) == pytest.approx(PAR_MM, abs=0.6)
    assert np.nanmedian(back.split['up_mean']) == pytest.approx(-UP_MM, abs=0.4)


def test_stat_split_needs_real_geometries(asc_desc, tmp_path, capsys):
    """Nominal angles (of ORBIT_DIRECTION, of a name) are not used for the
    split: one panel, and a note that says how to give the geometry."""
    obj = cli.stat_main([*asc_desc, '-sp', '10', '--boot', '200', *ALONG_THE_FAULT,
                         '--nodisplay', '--outdir', str(tmp_path)])
    out = capsys.readouterr().out
    assert obj.split_files == [] and obj.split is None and len(obj.axes_stat) == 1
    assert 'only nominal angles' in out and '--los-inc / --los-az' in out
    assert not any(p.name.endswith('_split.csv') for p in tmp_path.iterdir())
    assert 'the same profiles across the line on each' in out


def test_help_mentions_stat(capsys):
    cli.main(['--help'])
    out = capsys.readouterr().out
    assert 'stat - fault creep along a line' in out and 'violins' not in out
    with pytest.raises(SystemExit):
        cli.main(['stat', '--help'])
    out = capsys.readouterr().out
    for needle in ('-sp', '--length', '--width', '--boot', '--min-pixels', '--geom',
                   'LEFT side of', 'Hussain et al. (2016)', 'eq. 2', 'BLUE', 'bootstrap',
                   '--line-file', '_creep.csv', '_split.csv', '_values.csv', '_profiles.png',
                   'Backspace', 'fault-parallel', 'right-lateral', 'arrow'):
        assert needle in out, needle
    for gone in ('--radius', '--gap', '--rate', 'violin'):
        assert gone not in out, gone

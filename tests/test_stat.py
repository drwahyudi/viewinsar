"""`viewinsar stat`: block statistics on both sides of a line along a fault.

The geometry is checked on a metric grid (100 m pixels): the points along the
line, the circles - on their own side, every pixel at least the gap from the
line, pushed out or dropped on the inside of a bend - and the pixels counted.
Then the whole command on the synthetic products: the arctan interseismic
field across the N-S fault at 130.45 E gives known medians on either side.
"""

import json
import os

import numpy as np
import pytest

from conftest import DLAT, DLON, FAULT_LON, LAT0, LON0, _base_attrs, make_velocity
from viewinsar import stat

METRIC = {'Y_FIRST': '0', 'X_FIRST': '0', 'Y_STEP': '-100', 'X_STEP': '100',
          'Y_UNIT': 'meters', 'X_UNIT': 'meters', 'LENGTH': '800', 'WIDTH': '800'}


def frame_km(points_km, atr=METRIC):
    """A frame from vertices given in km (east, north) on the metric grid."""
    xy = np.asarray(points_km, dtype=float) * 1000.0
    return stat.polyline_frame(atr, np.c_[-xy[:, 1] / 100.0, xy[:, 0] / 100.0])


def expected_median(km_east):
    """The fixture's noise-free velocity, in cm/year, km east of the fault."""
    return 100 * 0.30 * (2 / np.pi) * np.arctan(km_east / 8.0)


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


def test_ray_hits_the_line():
    frame = frame_km([[0, -40], [40, -40]])
    assert stat.ray_hits_polyline(frame, [20000.0, -30000.0], [0, -1]) == pytest.approx(10000)
    assert stat.ray_hits_polyline(frame, [20000.0, -30000.0], [0, 1]) == np.inf
    assert stat.ray_hits_polyline(frame, [50000.0, -30000.0], [0, -1]) == np.inf


##############################################################################
#                                circles                                     #
##############################################################################
def test_circles_sit_radius_plus_gap_out_on_either_side():
    frame = frame_km([[0, -40], [40, -40]])
    dist, points, normals = stat.sample_points(frame, 5.0)
    centres, offset = stat.place_circles(frame, points, normals, 2.0, 1.0)
    assert offset == pytest.approx(np.full((2, len(dist)), 3.0))
    assert centres[0, :, 1] == pytest.approx(-37000.0)       # block 1: north (left)
    assert centres[1, :, 1] == pytest.approx(-43000.0)       # block 2: south (right)
    assert centres[0, :, 0] == pytest.approx(points[:, 0])


def kinked_case(gap=1.0, radius=2.0, spacing=2.5):
    """The sketch: a kinked line, a pair of circles at every point."""
    vertices = [[2, -10], [14, -22], [22, -24], [31, -38], [44, -44], [58, -47]]
    frame = frame_km(vertices)
    rng = np.random.default_rng(0)
    data = rng.normal(0, 1, (800, 800)).astype(np.float32)
    box = (0, 0, 800, 800)
    yx = np.asarray(vertices, dtype=float) * 1000.0
    yx = np.c_[-yx[:, 1] / 100.0, yx[:, 0] / 100.0]
    res = stat.stat_line(data, METRIC, box, yx, spacing, radius, gap)
    return frame, res


@pytest.mark.parametrize('gap', [1.0, 0.5, 0.0])
def test_every_pixel_is_on_its_side_and_clear_of_the_line(gap):
    frame, res = kinked_case(gap=gap)
    used = 0
    for b, sign in enumerate((1, -1)):
        for k in range(len(res['dist_km'])):
            rows, cols = res['rows'][b][k], res['cols'][b][k]
            if not rows.size:
                continue
            xy = stat.from_yx(frame, rows, cols)
            assert stat.distance_to_polyline(frame, xy).min() / 1000.0 >= gap - 1e-9
            assert np.all(stat.side_of_polyline(frame, xy) == sign)
            # and the circle really is the circle: pixel centres within the radius
            cen = res['centres'][b, k]
            assert np.hypot(*(xy - cen).T).max() <= 2000.0 + 1e-6
            used += rows.size
    assert used > 0


def test_inside_of_a_bend_is_pushed_out_or_dropped():
    frame, res = kinked_case()
    need = res['radius_km'] + res['gap_km']
    off = res['offset_km']
    assert np.isfinite(off).mean() > 0.8                     # most circles stay
    assert np.nanmin(off) == pytest.approx(need)
    assert np.nanmax(off) > need + 0.05                      # some were pushed
    # a pushed centre still clears the whole line by radius + gap
    good = np.isfinite(res['centres'][..., 0])
    clear = stat.distance_to_polyline(frame, res['centres'][good]) / 1000.0
    assert clear.min() >= need - 1e-9


def test_a_hairpin_drops_the_inside_circles():
    # east, then back west 1 km further south: a right-hand hairpin
    frame = frame_km([[0, -40], [30, -40], [0, -41]])
    dist, points, normals = stat.sample_points(frame, 5.0, smooth_km=2.0)
    centres, offset = stat.place_circles(frame, points, normals, 2.0, 1.0)
    # inside the hairpin (the right of both legs, block 2) there is no room at all
    assert np.all(~np.isfinite(offset[1]))
    # outside (block 1), every circle finds its spot
    assert np.all(np.isfinite(offset[0]))
    clear = stat.distance_to_polyline(frame, centres[0]) / 1000.0
    assert clear.min() >= 3.0 - 1e-9


def test_a_circle_never_jumps_over_the_line():
    """A serpentine: south of the first strand, the line comes back twice.

    Beyond the third strand there IS room, on the right of the line - but to
    get there from the first strand a circle would cross the line twice.
    """
    frame = frame_km([[0, -40], [30, -40], [30, -42], [0, -42], [0, -44], [30, -44]])
    dist, points, normals = stat.sample_points(frame, 5.0, smooth_km=2.0)
    centres, offset = stat.place_circles(frame, points, normals, 2.0, 1.0)
    first = (dist > 1) & (dist < 29)
    assert first.sum() >= 4
    assert np.all(~np.isfinite(offset[1, first]))
    # block 1 of the first strand (north, open ground) is untouched
    assert offset[0, first] == pytest.approx(3.0)


def test_pixels_in_a_circle():
    frame = frame_km([[0, -40], [40, -40]])
    data = np.ones((800, 800), dtype=np.float32)
    centre = np.array([20000.0, -30000.0])
    val, rows, cols, clipped = stat.circle_pixels(data, (0, 0, 800, 800), frame, centre, 2.0)
    assert val.size == pytest.approx(np.pi * 20 ** 2, rel=0.02)
    assert not clipped
    # NaN, masked, zero: left out (zero only by default)
    data[300, 200], data[301, 200], data[302, 200] = np.nan, 0.0, 5.0
    mask = np.ones_like(data, dtype=bool)
    mask[303, 200] = False
    n_all = val.size
    val = stat.circle_pixels(data, (0, 0, 800, 800), frame, centre, 2.0, mask=mask)[0]
    assert val.size == n_all - 3
    val = stat.circle_pixels(data, (0, 0, 800, 800), frame, centre, 2.0, mask=mask,
                             keep_zero=True)[0]
    assert val.size == n_all - 2
    # a window that cuts the circle: only the part inside, and it says so
    sub = data[:, :200]
    val, rows, cols, clipped = stat.circle_pixels(sub, (0, 0, 200, 800), frame, centre, 2.0)
    assert clipped and cols.max() == 199 and val.size < n_all * 0.6
    # no centre (a dropped circle): nothing
    assert stat.circle_pixels(data, (0, 0, 800, 800), frame, [np.nan, np.nan], 2.0)[0].size == 0


def test_summarize():
    out = stat.summarize([1, 2, 3, 4, 100])
    assert out['count'] == 5 and out['median'] == 3 and out['max'] == 100
    assert out['q25'] == 2 and out['q75'] == 4 and out['mean'] == 22
    assert np.isfinite(out['mode'])
    empty = stat.summarize([])
    assert empty['count'] == 0 and np.isnan(empty['median']) and np.isnan(empty['mode'])


def test_mode_is_the_peak_of_the_density():
    rng = np.random.default_rng(0)
    skewed = rng.gamma(3.0, 1.0, 20000)                      # the mode of Gamma(3, 1) is 2
    assert stat.kde_mode(skewed) == pytest.approx(2.0, abs=0.15)
    assert stat.kde_mode(skewed) < np.median(skewed) < skewed.mean()
    assert stat.summarize(skewed)['mode'] == stat.kde_mode(skewed)     # what the csv gets
    # two groups: the taller peak wins
    two = np.r_[rng.normal(0, 1, 5000), rng.normal(5, 0.5, 2000)]
    assert stat.kde_mode(two) == pytest.approx(0.0, abs=0.2)
    # too few values, or all equal: the median; nothing: NaN
    assert stat.kde_mode([1.0, 7.0]) == 4.0
    assert stat.kde_mode([3.0, 3.0, 3.0, 3.0]) == 3.0
    assert np.isnan(stat.kde_mode([]))
    assert stat.kde_mode([np.nan, 1.0, 1.0, 2.0, np.inf]) == pytest.approx(1.0, abs=0.2)


def test_mode_sits_where_the_violin_is_widest():
    """The same kernel and bandwidth as matplotlib's violins."""
    mlab = pytest.importorskip('matplotlib.mlab')
    rng = np.random.default_rng(4)
    for values in (rng.gamma(2.0, 1.0, 400), rng.normal(1.0, 0.2, 60),
                   np.r_[rng.normal(-1, 0.3, 300), rng.normal(1, 0.5, 200)]):
        grid = np.linspace(values.min(), values.max(), 200)
        widest = grid[np.argmax(mlab.GaussianKDE(values, 'scott').evaluate(grid))]
        assert stat.kde_mode(values) == pytest.approx(widest, abs=1.5 * (grid[1] - grid[0]))


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
    atr = dict(METRIC)
    frame = stat.polyline_frame(atr, [[400, 100], [400, 700]], mirrored=True)
    # columns now run west: walking along +col is walking WEST, the left is south
    _, points, normals = stat.sample_points(frame, 10)
    assert normals[0] == pytest.approx([0, -1])
    centres = stat.place_circles(frame, points, normals, 2, 1)[0]
    rows, _ = stat.to_yx(frame, centres[0])
    assert np.all(rows > 400)                                 # south = larger row
    assert stat.to_yx(frame, stat.from_yx(frame, 12.0, 34.0)) == pytest.approx((12.0, 34.0))


##############################################################################
#                       the whole line, on the fixture field                 #
##############################################################################
def test_stat_line_recovers_the_block_velocities():
    """A line up the fault: block 1 (west) moves -, block 2 (east) moves +."""
    atr = _base_attrs()
    vel = make_velocity()
    col = (FAULT_LON - LON0) / DLON - 0.5
    row = lambda lat: (lat - LAT0) / DLAT - 0.5               # noqa: E731
    line = [[row(30.05), col], [row(30.52), col]]             # south -> north
    res = stat.stat_line(vel, atr, (0, 0, vel.shape[1], vel.shape[0]), line, 5.0, 2.0, 1.0)
    assert res['length_km'] == pytest.approx(0.47 * 111.19, rel=0.01)
    assert len(res['dist_km']) == 11
    med = np.array([[s['median'] for s in res['stats'][b]] for b in range(2)]) * 100
    want = expected_median(3.0)                               # the circle centres: 3 km out
    assert np.median(med[0]) == pytest.approx(-want, abs=0.6)
    assert np.median(med[1]) == pytest.approx(want, abs=0.6)
    assert np.all(med[0] < 0) and np.all(med[1] > 0)
    # every pixel of block 1 is west of the fault, of block 2 east, by more than 1 km
    for b, sign in enumerate((-1, 1)):
        cols = np.concatenate(res['cols'][b])
        km = ((cols + 0.5) * DLON + LON0 - FAULT_LON) * 111.32 * np.cos(np.radians(30.3))
        assert np.all(sign * km > 1.0 - 0.01)


##############################################################################
#                           the command line                                 #
##############################################################################
mintpy = pytest.importorskip('mintpy', reason='the viewer needs MintPy')
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from viewinsar import cli  # noqa: E402

UP_THE_FAULT = ['--line-lalo', '30.05', str(FAULT_LON), '30.52', str(FAULT_LON)]


def test_stat_options(velocity_file):
    inps = cli.stat_cmd_line_parse([velocity_file, '--noverbose', *UP_THE_FAULT])
    assert (inps.spacing, inps.radius, inps.gap) == (5.0, 2.0, 1.0)
    assert inps.line_lalo == [[30.05, FAULT_LON], [30.52, FAULT_LON]]
    assert inps.block_names == ['Block 1', 'Block 2']
    inps = cli.stat_cmd_line_parse([velocity_file, '--noverbose', '-sp', '2.5', '-r', '1',
                                    '--gap', '0.5', '--line-yx', '10', '20', '30', '40',
                                    '50', '60', '-c', 'RdBu_r', '--names', 'W', 'E'])
    assert (inps.spacing, inps.radius, inps.gap) == (2.5, 1.0, 0.5)
    assert inps.line_yx == [[10, 20], [30, 40], [50, 60]]
    assert inps.block_names == ['W', 'E']
    assert inps.view_extra == ['-c', 'RdBu_r']               # the rest goes to view.py

    bad = [['-sp', '0'], ['-r', '-1'], ['--gap', '-0.1'], ['--line-yx', '1', '2', '3'],
           ['--line-yx', '1', '2'], ['--colors', 'red', 'notacolour'],
           ['--line-file', 'nowhere.shp'], ['--nodisplay'],
           ['--line-yx', '1', '2', '3', '4', '--line-lalo', '30', '130', '30.1', '130']]
    for extra in bad:
        with pytest.raises(ValueError):
            cli.stat_cmd_line_parse([velocity_file, '--noverbose', *extra])


def _read_csv(path):
    return np.genfromtxt(path, delimiter=',', names=True)


def test_stat_end_to_end(velocity_file, tmp_path):
    obj = cli.stat_main([velocity_file, '-sp', '10', '-r', '2', *UP_THE_FAULT,
                         '--nodisplay', '--noverbose', '--outdir', str(tmp_path)])
    names = sorted(os.listdir(tmp_path))
    base = os.path.basename(obj.outfile_base)
    assert names == sorted(f'{base}{tail}' for tail in
                           ('.png', '_map.png', '_stats.csv', '_values.csv', '_info.txt'))
    assert base.endswith('_sp10km_r2km')

    st = _read_csv(tmp_path / f'{base}_stats.csv')
    assert len(st) == 6
    want = expected_median(3.0)                               # cm/year, view.py's unit
    assert np.median(st['median_b1']) == pytest.approx(-want, abs=0.8)
    assert np.median(st['median_b2']) == pytest.approx(want, abs=0.8)
    assert st['median_b1_minus_b2'] == pytest.approx(st['median_b1'] - st['median_b2'])
    assert st['offset_km_b1'] == pytest.approx(3.0) and st['offset_km_b2'] == pytest.approx(3.0)
    assert np.all(st['count_b1'] > 200) and np.all(st['count_b2'] > 200)
    assert st['lon'] == pytest.approx(FAULT_LON, abs=0.002)
    assert np.diff(st['dist_km']) == pytest.approx(10.0)

    val = _read_csv(tmp_path / f'{base}_values.csv')
    assert len(val) == st['count_b1'].sum() + st['count_b2'].sum()
    assert val['dist_to_line_km'].min() >= 1.0 - 1e-6
    assert np.all(val['lon'][val['block'] == 1] < FAULT_LON)
    assert np.all(val['lon'][val['block'] == 2] > FAULT_LON)

    info = (tmp_path / f'{base}_info.txt').read_text()
    assert 'LEFT of the line' in info and 'gap            : 1.0 km' in info

    # the mode of each circle, next to its median and mean
    for b in (1, 2):
        assert np.all(np.isfinite(st[f'mode_b{b}']))
        assert np.all(np.abs(st[f'mode_b{b}'] - st[f'median_b{b}']) < 3.0)

    # the plot: two violins per point, two dashed median lines
    ax = obj.ax_stat
    dashed = [ln for ln in ax.get_lines() if ln.get_linestyle() == '--']
    assert len(dashed) == 2
    assert len(dashed[0].get_xdata()) == 6
    assert dashed[0].get_color() != dashed[1].get_color()
    # just the block names in the legend - which side is which is on the map
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels == ['Block 1', 'Block 2', 'Median', 'Mean', 'Mode']
    # a strike per violin at its median, mean and mode
    strikes = {c.get_gid(): c for c in ax.collections if c.get_gid()}
    for b in (1, 2):
        for key in ('median', 'mean', 'mode'):
            ys = [seg[0][1] for seg in strikes[f'{key}_b{b}'].get_segments()]
            assert ys == pytest.approx(st[f'{key}_b{b}'], abs=1e-5)
    assert strikes['mode_b1'].get_linestyle() != strikes['mean_b1'].get_linestyle()

    # the colours: a block's circles on the map are filled and edged as its violins
    bodies = [c for c in ax.collections if not c.get_gid()]
    assert len(bodies) == 12
    circles = obj.ax_img.patches
    assert len(circles) == 12
    for body, circle in zip(bodies[::6], circles[::6]):
        assert body.get_facecolor()[0][:3] == pytest.approx(circle.get_facecolor()[:3])
        assert body.get_edgecolor()[0][:3] == pytest.approx(circle.get_edgecolor()[:3])
        assert body.get_facecolor()[0][3] == 1.0              # not see-through
    assert bodies[0].get_facecolor()[0][:3] != pytest.approx(bodies[6].get_facecolor()[0][:3])
    # and the map has the key
    keys = [a for a in obj.ax_img.get_children() if isinstance(a, matplotlib.legend.Legend)]
    assert [t.get_text() for t in keys[0].get_texts()] == ['Block 1', 'Block 2']


def test_stat_interactive(velocity_file, tmp_path):
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    from viewinsar.statview import StatViewer

    inps = cli.stat_cmd_line_parse([velocity_file, '-sp', '10', '--noverbose',
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

    def key(name):
        obj.on_key(KeyEvent('key_press_event', obj.fig_map.canvas, key=name))

    click(FAULT_LON, 30.05)
    click(FAULT_LON + 0.01, 30.30)
    click(FAULT_LON, 30.45)
    assert len(obj.vertices) == 3
    key('u')                                   # undo the last point
    assert len(obj.vertices) == 2
    click(FAULT_LON, 30.52)
    click(FAULT_LON, 30.52, button=MouseButton.RIGHT)          # ignored
    assert len(obj.vertices) == 3
    key('enter')
    assert obj.finished and obj.result is not None
    med = [np.nanmedian([s['median'] for s in obj.result['stats'][b]]) for b in range(2)]
    assert med[0] < -3 and med[1] > 3
    # the map shows the circles: one patch per circle
    n_circles = int(np.isfinite(obj.result['offset_km']).sum())
    assert len(obj.ax_img.patches) >= n_circles

    key('e')
    assert any(p.name.endswith('_stats.csv') for p in tmp_path.iterdir())

    # a click after Enter starts a new line
    click(FAULT_LON, 30.10)
    assert not obj.finished and len(obj.vertices) == 1
    key('n')
    assert obj.vertices == [] and obj.result is None
    obj.close()


def test_stat_line_file_takes_the_longest_line(velocity_file, tmp_path):
    obj = {'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': {},
         'geometry': {'type': 'LineString', 'coordinates': [[130.2, 30.2], [130.25, 30.22]]}},
        {'type': 'Feature', 'properties': {},
         'geometry': {'type': 'LineString',
                      'coordinates': [[FAULT_LON, 30.05], [FAULT_LON, 30.30],
                                      [FAULT_LON, 30.52]]}}]}
    path = tmp_path / 'trace.geojson'
    path.write_text(json.dumps(obj))
    viewer = cli.stat_main([velocity_file, '-sp', '10', '--line-file', str(path), '--nodisplay',
                            '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert viewer.result['length_km'] == pytest.approx(0.47 * 111.19, rel=0.01)
    med = [np.nanmedian([s['median'] for s in viewer.result['stats'][b]]) for b in range(2)]
    assert med[0] < -3 and med[1] > 3


def test_stat_on_radar_coded_data(tmp_path):
    """Descending: the grid is a mirror image of the ground; west must stay block 1."""
    from test_radar import make_radar_product

    vfile, _ = make_radar_product(str(tmp_path / 'proc'))
    obj = cli.stat_main([vfile, '-sp', '5', '-r', '2', '--line-lalo', '30.30',
                         str(FAULT_LON), '30.50', str(FAULT_LON), '--nodisplay',
                         '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert obj.mirrored is True
    st = _read_csv(next((tmp_path / 'out').glob('*_stats.csv')))
    assert np.all(st['median_b1'] < 0) and np.all(st['median_b2'] > 0)
    assert np.median(st['median_b2']) == pytest.approx(expected_median(3.0), abs=1.5)
    val = _read_csv(next((tmp_path / 'out').glob('*_values.csv')))
    assert np.all(val['lon'][val['block'] == 1] < FAULT_LON)
    assert np.all(val['lon'][val['block'] == 2] > FAULT_LON)


def test_help_mentions_stat(capsys):
    cli.main(['--help'])
    out = capsys.readouterr().out
    assert 'stat' in out and 'violins' in out
    with pytest.raises(SystemExit):
        cli.main(['stat', '--help'])
    out = capsys.readouterr().out
    for needle in ('-sp', '--radius', '--gap', 'LEFT side of the line', '--line-file',
                   '_values.csv', 'median_b1_minus_b2', 'Backspace'):
        assert needle in out, needle

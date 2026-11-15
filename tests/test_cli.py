"""End-to-end tests of the CLI and the viewer. These need MintPy installed."""

import os
import sys

import numpy as np
import pytest

import viewinsar
from viewinsar import cli

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')


##############################################################################
#                              plumbing                                      #
##############################################################################
def test_version_and_help(capsys):
    assert cli.main(['--version']) == 0
    assert viewinsar.__version__ in capsys.readouterr().out
    assert cli.main([]) == 0
    out = capsys.readouterr().out
    assert 'swath' in out and 'list' in out


def test_list_subcommand(velocity_file, capsys):
    cli.main(['list', velocity_file])
    assert capsys.readouterr().out.strip().splitlines() == ['velocity  <- default']


def test_bare_file_implies_swath(velocity_file, monkeypatch):
    """`viewinsar FILE` must behave like `viewinsar swath FILE`."""
    seen = {}
    monkeypatch.setitem(cli.SUBCOMMANDS, 'swath',
                        (lambda argv: seen.setdefault('argv', argv), 'x'))
    cli.main([velocity_file, '-w', '3'])
    assert seen['argv'] == [velocity_file, '-w', '3']


def test_view_options_are_forwarded(velocity_file):
    inps = cli.cmd_line_parse([velocity_file, '-w', '5', '-c', 'RdBu', '-v', '-3', '3',
                               '--noverbose'])
    assert inps.view_extra == ['-c', 'RdBu', '-v', '-3', '3']
    cmd = cli.get_view_cmd(inps)
    assert cmd.startswith('view.py ') and '-c RdBu' in cmd and '--noverbose' in cmd
    assert inps.width == 5.0            # our --width must survive view.py's namesake


def test_lalo_is_converted_to_yx(velocity_file):
    inps = cli.cmd_line_parse([velocity_file, '--start-lalo', '30.30', '130.15',
                               '--end-lalo', '30.30', '130.75', '--noverbose'])
    assert inps.start_yx is not None and inps.end_yx is not None
    assert inps.start_yx[0] == inps.end_yx[0]          # same latitude -> same row
    assert inps.start_yx[1] < inps.end_yx[1]


def test_line_file(velocity_file, tmp_path):
    path = tmp_path / 'line.xy'
    path.write_text('>\n130.1500 30.2000\n130.7500 30.4000\n')
    inps = cli.cmd_line_parse([velocity_file, '--line-file', str(path), '--noverbose'])
    assert inps.start_lalo == [30.2, 130.15]
    assert inps.end_lalo == [30.4, 130.75]


def test_several_profiles_up_front(velocity_file, tmp_path):
    """A-A' from --start-lalo, then one per --line-file segment, then --profile
    and --profile-yx, in that order."""
    path = tmp_path / 'lines.xy'
    path.write_text("# two sections\n>  B-B'\n130.20 30.50\n130.30 30.10\n"
                    ">  C-C'\n130.40 30.50\n130.45 30.30\n130.50 30.10\n")
    inps = cli.cmd_line_parse([velocity_file, '--start-lalo', '30.30', '130.15',
                               '--end-lalo', '30.30', '130.75', '--line-file', str(path),
                               '--profile', '30.45', '130.60', '30.15', '130.70',
                               '-profile', '30.40', '130.20', '30.20', '130.20',
                               '--profile-yx', '10', '20', '60', '80', '--noverbose'])
    profs = inps.profiles
    assert len(profs) == 6
    assert profs[0]['start_yx'] == inps.start_yx and inps.start_lalo == [30.30, 130.15]
    # a segment of three points runs from its first to its last
    assert profs[1]['start_lalo'] == [30.5, 130.2] and profs[1]['end_lalo'] == [30.1, 130.3]
    assert profs[2]['start_lalo'] == [30.5, 130.4] and profs[2]['end_lalo'] == [30.1, 130.5]
    assert profs[3]['start_lalo'] == [30.45, 130.60]
    assert profs[5] == dict(start_yx=[10, 20], end_yx=[60, 80], start_lalo=None,
                            end_lalo=None)
    # the conversion: same longitude -> same column, north -> a smaller row
    assert profs[4]['start_yx'][1] == profs[4]['end_yx'][1]
    assert profs[4]['start_yx'][0] < profs[4]['end_yx'][0]


def test_line_file_vector_and_errors(velocity_file, tmp_path):
    """A GeoJSON of two lines (one in UTM) is two sections; bad input is caught."""
    import json

    from pyproj import Transformer

    tfm = Transformer.from_crs(4326, 32652, always_xy=True)
    e0, n0 = tfm.transform(130.20, 30.50)
    e1, n1 = tfm.transform(130.30, 30.10)
    utm = tmp_path / 'utm.geojson'
    utm.write_text(json.dumps({
        'type': 'FeatureCollection',
        'crs': {'type': 'name', 'properties': {'name': 'urn:ogc:def:crs:EPSG::32652'}},
        'features': [{'type': 'Feature', 'properties': {}, 'geometry': {
            'type': 'MultiLineString', 'coordinates': [[[e0, n0], [e1, n1]],
                                                       [[e1, n1], [e0, n0]]]}}]}))
    lines = cli.read_profile_file(str(utm))
    assert len(lines) == 2
    assert lines[0][0] == pytest.approx([30.50, 130.20], abs=1e-6)
    assert lines[1][1] == pytest.approx([30.50, 130.20], abs=1e-6)
    inps = cli.cmd_line_parse([velocity_file, '--line-file', str(utm), '--noverbose'])
    assert len(inps.profiles) == 2 and inps.start_lalo == pytest.approx([30.50, 130.20])

    with pytest.raises(ValueError, match='needs --end-lalo'):
        cli.cmd_line_parse([velocity_file, '--start-lalo', '30.3', '130.2', '--noverbose'])
    with pytest.raises(ValueError, match='same pixel'):
        cli.cmd_line_parse([velocity_file, '--profile', '30.3', '130.2', '30.3', '130.2',
                            '--noverbose'])
    with pytest.raises(ValueError, match='no such file'):
        cli.cmd_line_parse([velocity_file, '--line-file', str(tmp_path / 'none.xy'),
                            '--noverbose'])
    empty = tmp_path / 'empty.xy'
    empty.write_text('>\n130.2 30.5\n')
    with pytest.raises(ValueError, match='no line of two or more points'):
        cli.cmd_line_parse([velocity_file, '--line-file', str(empty), '--noverbose'])


@pytest.mark.parametrize('bad', [['-w', '0'], ['--bin-x', '-1']])
def test_bad_arguments_raise(velocity_file, bad):
    with pytest.raises(ValueError):
        cli.cmd_line_parse([velocity_file, *bad, '--noverbose'])


##############################################################################
#                             end to end                                     #
##############################################################################
def _run(args):
    obj = cli.swath_main(args)
    return obj


def test_headless_run_writes_everything(velocity_file, tmp_path):
    out = tmp_path / 'out'
    obj = _run([velocity_file, '--start-lalo', '30.30', '130.15',
                '--end-lalo', '30.30', '130.75', '-w', '10',
                '--bin-x', '0.5', '--bin-y', '1.5', '--ylim', '-40', '40',
                '--save', '--nodisplay', '--outdir', str(out), '--noverbose'])

    files = sorted(p.name for p in out.iterdir())
    assert any(f.endswith('_points.csv') for f in files)
    assert any(f.endswith('_stats.csv') for f in files)
    assert any(f.endswith('_info.txt') for f in files)
    assert any(f.endswith('.png') for f in files)
    assert any(f.startswith('swath_map_') for f in files)      # the map is saved too

    pts = np.genfromtxt(next(out.glob('*_points.csv')), delimiter=',', names=True)
    assert set(pts.dtype.names) >= {'dist_km', 'perp_km', 'value', 'row', 'col', 'lat', 'lon'}
    assert pts.size == obj.sections[0].swath_list[0]['value'].size
    assert np.abs(pts['perp_km']).max() <= 10.0 + 1e-4
    assert pts['dist_km'].min() >= -1e-4

    stats = np.genfromtxt(next(out.glob('*_stats.csv')), delimiter=',', names=True)
    assert stats['count'].sum() == pts.size
    # the arctan step must cross zero near the fault, ~ 28.8 km along this line
    crossing = stats['dist_km'][np.nanargmin(np.abs(stats['median']))]
    assert crossing == pytest.approx(28.8, abs=2.0)


def test_headless_several_profiles(velocity_file, tmp_path):
    """Three sections up front: three windows, three sets of files, one map."""
    out = tmp_path / 'multi'
    obj = _run([velocity_file, '--start-lalo', '30.30', '130.15', '--end-lalo', '30.30',
                '130.75', '--profile', '30.45', '130.20', '30.45', '130.70',
                '--profile', '30.15', '130.20', '30.15', '130.70', '-w', '5',
                '--save', '--nodisplay', '--outdir', str(out), '--noverbose'])
    assert [s.name for s in obj.sections] == ["A-A'", "B-B'", "C-C'"]
    assert all(s.has_data for s in obj.sections)
    assert len({id(s.fig) for s in obj.sections}) == 3          # a window each
    files = sorted(p.name for p in out.iterdir())
    for label in 'ABC':
        assert any(f.startswith(f'swath_{label}_') and f.endswith('_points.csv')
                   for f in files), label
        assert any(f.startswith(f'swath_{label}_') and f.endswith('.png') for f in files)
    assert 'swath_map_ABC.png' in files
    # each section is its own line: 0.6 deg, 0.5 deg and 0.5 deg of longitude
    for sec, deg, lat in zip(obj.sections, (0.6, 0.5, 0.5), (30.30, 30.45, 30.15)):
        want = deg * 111.32 * np.cos(np.deg2rad(lat))
        assert sec.swath_list[0]['length'] == pytest.approx(want, rel=0.01), sec.name
    texts = [t.get_text() for t in obj.ax_img.texts]
    for label in ('A', "A'", 'B', "B'", 'C', "C'"):
        assert label in texts, label
    obj.close()


def test_profile_length_is_physical(velocity_file, tmp_path):
    obj = _run([velocity_file, '--start-lalo', '30.30', '130.15',
                '--end-lalo', '30.30', '130.75', '-w', '5',
                '--save', '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    expected = 0.6 * 111.32 * np.cos(np.deg2rad(30.3))       # 0.6 deg of longitude
    assert obj.sections[0].swath_list[0]['length'] == pytest.approx(expected, rel=0.01)


def test_x_offset_is_applied_to_the_csv(velocity_file, tmp_path):
    out = tmp_path / 'off'
    _run([velocity_file, '--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75',
          '-w', '5', '--x-offset', '28.8', '--save', '--nodisplay',
          '--outdir', str(out), '--noverbose'])
    pts = np.genfromtxt(next(out.glob('*_points.csv')), delimiter=',', names=True)
    assert pts['dist_km'].min() == pytest.approx(-28.8, abs=0.2)
    assert pts['dist_km'].max() > 0


def test_yx_input_and_subset(velocity_file, tmp_path):
    obj = _run([velocity_file, '--start-yx', '150', '75', '--end-yx', '150', '375',
                '-w', '5', '--sub-lat', '30.15', '30.45', '--sub-lon', '130.05', '130.85',
                '--save', '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    assert obj.sections[0].swath_list[0]['value'].size > 0
    assert obj.pix_box[3] - obj.pix_box[1] < 300      # the subset really was applied


def test_interactive_clicks(velocity_file, tmp_path):
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, '-w', '8', '--bin-x', '0.5',
                              '--outdir', str(tmp_path), '--noverbose', '--nodisplay'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    obj.fig_map.canvas.draw()

    # the map and the section live in SEPARATE windows
    assert obj.fig_map is not obj.sections[0].fig
    assert obj.fig_map.get_label() == 'Map'
    assert obj.sections[0].label == 'A'

    def click(lon, lat, button=MouseButton.LEFT):
        x, y = obj.ax_img.transData.transform((lon, lat))
        evt = MouseEvent('button_release_event', obj.fig_map.canvas, x, y, button=button)
        assert evt.inaxes is obj.ax_img
        obj.select_point(evt)

    def key(k):
        x, y = obj.ax_img.transData.transform((130.5, 30.3))
        obj.on_key(KeyEvent('key_press_event', obj.fig_map.canvas, k, x, y))

    sec = obj.sections[0]
    assert not sec.has_data
    click(130.18, 30.25)
    assert not sec.has_data                  # one click is not enough
    click(130.72, 30.38)

    swath = sec.swath_list[0]
    assert np.abs(swath['perp']).max() <= 8.0 + 1e-6
    assert swath['dist'].min() >= -1e-6
    assert swath['dist'].max() <= swath['length'] + 1e-6
    assert sec.poly_obj is not None and sec.line_obj is not None
    assert sec.txt0.get_text() == 'A' and sec.txt1.get_text() == "A'"

    # "e" exports, "n" clears only the current section
    key('e')
    assert list(tmp_path.glob('*_A_*_points.csv'))
    assert list(tmp_path.glob('swath_map_*'))
    key('n')
    assert not sec.has_data and sec.start_yx is None and sec.poly_obj is None


def test_multiple_sections(velocity_file, tmp_path):
    """"=" starts B-B' in a new window and leaves A-A' on the map."""
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, '-w', '8', '--svg',
                               '--outdir', str(tmp_path), '--noverbose', '--nodisplay'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    obj.fig_map.canvas.draw()

    def click(lon, lat):
        x, y = obj.ax_img.transData.transform((lon, lat))
        obj.select_point(MouseEvent('button_release_event', obj.fig_map.canvas,
                                    x, y, button=MouseButton.LEFT))

    def key(k, fig=None):
        x, y = obj.ax_img.transData.transform((130.5, 30.3))
        obj.on_key(KeyEvent('key_press_event', (fig or obj.fig_map).canvas, k, x, y))

    click(130.15, 30.22)
    click(130.78, 30.40)
    key('=')
    assert len(obj.sections) == 2 and obj.sections[1].label == 'B'
    click(130.20, 30.42)
    click(130.80, 30.18)
    # "+" works too, and from a section window
    key('+', fig=obj.sections[1].fig)
    click(130.10, 30.30)
    click(130.85, 30.33)

    assert [s.label for s in obj.sections] == ['A', 'B', 'C']
    assert len({id(s.fig) for s in obj.sections}) == 3          # three windows
    assert all(s.has_data for s in obj.sections)
    # every earlier line is still drawn on the map, with its own labels
    assert all(s.line_obj is not None for s in obj.sections)
    assert [s.txt1.get_text() for s in obj.sections] == ["A'", "B'", "C'"]

    key('e')
    svg = sorted(p.name for p in tmp_path.glob('*.svg'))
    assert len(svg) == 4, svg                       # map + 3 sections
    assert any(f.startswith('swath_map_ABC') for f in svg)
    assert sum(1 for p in tmp_path.glob('*_points.csv')) == 3
    assert next(tmp_path.glob('*_A_*.svg')).read_text()[:200].lstrip().startswith('<?xml')


def test_new_section_needs_a_drawn_one(velocity_file, tmp_path, capsys):
    from matplotlib.backend_bases import KeyEvent

    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, '-w', '8', '--outdir', str(tmp_path),
                               '--noverbose', '--nodisplay'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    x, y = obj.ax_img.transData.transform((130.5, 30.3))
    obj.on_key(KeyEvent('key_press_event', obj.fig_map.canvas, '=', x, y))
    assert len(obj.sections) == 1                   # refused: A has no profile yet
    assert 'no profile yet' in capsys.readouterr().out


def test_timeseries_panel(velocity_file, timeseries_file, tmp_path):
    from matplotlib.backend_bases import MouseButton, MouseEvent

    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, '-w', '8', '--timeseries', timeseries_file,
                               '--outdir', str(tmp_path), '--noverbose', '--nodisplay'])
    inps.save_fig = False
    obj = swathProfileViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    if obj.ts_obj is None:
        pytest.skip('tsview could not be initiated in this environment')
    obj.plot()
    obj.fig_map.canvas.draw()
    assert obj.ax_ts is not None
    assert obj.fig_ts is not obj.fig_map                # its own window

    x, y = obj.ax_img.transData.transform((130.55, 30.33))
    obj.select_point(MouseEvent('button_release_event', obj.fig_map.canvas, x, y,
                                button=MouseButton.RIGHT))
    assert obj.ax_ts.lines or obj.ax_ts.collections


def test_console_entry_point_is_installed():
    """`viewinsar` must be on PATH once the package is installed."""
    import shutil
    import subprocess
    exe = shutil.which('viewinsar')
    if exe is None:
        pytest.skip('the package is not pip-installed in this environment')
    out = subprocess.run([sys.executable, '-m', 'viewinsar.cli', '--version'],
                         capture_output=True, text=True,
                         env={**os.environ, 'MPLBACKEND': 'Agg'})
    assert viewinsar.__version__ in out.stdout


##############################################################################
#        0.5.0: one --long spelling per option; the old ones still work      #
##############################################################################
@pytest.mark.parametrize('command, old, new', [
    (cmd, old, new) for cmd, table in cli.OLD_SPELLINGS.items() for old, new in table.items()])
def test_every_old_spelling_is_an_option_of_today(command, old, new):
    """Each spelling of before maps onto an option the command has now."""
    parser = {'swath': cli.create_swath_parser, 'stat': cli.create_stat_parser,
              'model': cli.create_model_parser, 'tseries': cli.create_tseries_parser,
              'list': cli.create_list_parser}[command]()
    known = set(parser._option_string_actions)
    assert new in known and old not in known


def test_old_spellings_are_rewritten_with_a_note(capsys):
    argv = ['v.h5', '-fault', 'f.shp', '-gnss', 'g.csv', '-vel', 'los', '-d', '-w15',
            '--buffer=8', '-v', '-3', '3', '-u', 'mm/yr', '--start-lalo', '-0.8', '121.1']
    got = cli.modern_options(argv, 'swath')
    assert got == ['v.h5', '--fault', 'f.shp', '--gnss', 'g.csv', '--vel', 'los', '--dem',
                   '--width', '15', '--width=8', '-v', '-3', '3', '-u', 'mm/yr',
                   '--start-lalo', '-0.8', '121.1']        # view.py's -v, -u untouched
    out = capsys.readouterr().out
    assert out.count('NOTE:') == 1 and '--fault for -fault' in out and '--width for -w' in out
    # today's spellings: nothing to say
    assert cli.modern_options(got, 'swath') == got and not capsys.readouterr().out
    # tseries' own letters are tseries' own
    assert cli.modern_options(['cum.h5', '-r', '2', '-m', '20170105', '-p', '10/20', '-c',
                               'jet', '-u', 'cm'], 'tseries') == [
        'cum.h5', '--radius', '2', '--ref-date', '20170105', '--point', '10/20',
        '--colormap', 'jet', '--unit', 'cm']


def test_old_and_new_spellings_parse_the_same(velocity_file):
    old = cli.cmd_line_parse([velocity_file, '-w', '12', '-nlook', '2', '-svg', '--noverbose'])
    new = cli.cmd_line_parse([velocity_file, '--width', '12', '--nlook', '2', '--svg',
                              '--noverbose'])
    assert (old.width, old.nlook, old.svg) == (new.width, new.nlook, new.svg) == (12.0, 2, True)
    assert old.argv == new.argv[:1] + ['--width', '12', '--nlook', '2', '--svg', '--noverbose']


def test_swath_gnss_nofit_is_gone_but_harmless(capsys):
    """0.5.1: swath draws no line through the GNSS, so --gnss-nofit is left out
    there (with a note); model keeps its own - the GNSS shown, not fitted."""
    assert '--gnss-nofit' not in cli.create_swath_parser()._option_string_actions
    assert '--gnss-nofit' in cli.create_model_parser()._option_string_actions
    assert cli.modern_options(['v.h5', '--gnss-nofit', '--width', '5'], 'swath') == [
        'v.h5', '--width', '5']
    assert 'not needed any more' in capsys.readouterr().out
    assert cli.modern_options(['v.h5', '--gnss-nofit'], 'model') == ['v.h5', '--gnss-nofit']

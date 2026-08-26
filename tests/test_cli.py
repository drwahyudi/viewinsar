"""End-to-end tests of the CLI and the viewer. These need MintPy installed."""

import os
import sys

import numpy as np
import pytest

import viewinsar
from viewinsar import cli

pytest.importorskip('mintpy', reason='the viewer needs MintPy')
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
    assert capsys.readouterr().out.strip().splitlines() == ['velocity']


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

    pts = np.genfromtxt(next(out.glob('*_points.csv')), delimiter=',', names=True)
    assert set(pts.dtype.names) >= {'dist_km', 'perp_km', 'value', 'row', 'col', 'lat', 'lon'}
    assert pts.size == obj.swath_list[0]['value'].size
    assert np.abs(pts['perp_km']).max() <= 10.0 + 1e-4
    assert pts['dist_km'].min() >= -1e-4

    stats = np.genfromtxt(next(out.glob('*_stats.csv')), delimiter=',', names=True)
    assert stats['count'].sum() == pts.size
    # the arctan step must cross zero near the fault, ~ 28.8 km along this line
    crossing = stats['dist_km'][np.nanargmin(np.abs(stats['median']))]
    assert crossing == pytest.approx(28.8, abs=2.0)


def test_profile_length_is_physical(velocity_file, tmp_path):
    obj = _run([velocity_file, '--start-lalo', '30.30', '130.15',
                '--end-lalo', '30.30', '130.75', '-w', '5',
                '--save', '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    expected = 0.6 * 111.32 * np.cos(np.deg2rad(30.3))       # 0.6 deg of longitude
    assert obj.swath_list[0]['length'] == pytest.approx(expected, rel=0.01)


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
    assert obj.swath_list[0]['value'].size > 0
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
    obj.fig.canvas.draw()

    def click(lon, lat, button=MouseButton.LEFT):
        x, y = obj.ax_img.transData.transform((lon, lat))
        evt = MouseEvent('button_release_event', obj.fig.canvas, x, y, button=button)
        assert evt.inaxes is obj.ax_img
        obj.select_point(evt)

    assert obj.swath_list is None
    click(130.18, 30.25)
    assert obj.swath_list is None            # one click is not enough
    click(130.72, 30.38)

    swath = obj.swath_list[0]
    assert np.abs(swath['perp']).max() <= 8.0 + 1e-6
    assert swath['dist'].min() >= -1e-6
    assert swath['dist'].max() <= swath['length'] + 1e-6
    assert obj.poly_obj is not None and obj.line_obj is not None

    # "e" exports, "n" resets
    x, y = obj.ax_swath.transData.transform((5, 0))
    obj.on_key(KeyEvent('key_press_event', obj.fig.canvas, 'e', x, y))
    assert list(tmp_path.glob('*_points.csv'))
    obj.on_key(KeyEvent('key_press_event', obj.fig.canvas, 'n', x, y))
    assert obj.swath_list is None and obj.start_yx is None and obj.poly_obj is None


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
    obj.fig.canvas.draw()
    assert obj.ax_ts is not None

    x, y = obj.ax_img.transData.transform((130.55, 30.33))
    obj.select_point(MouseEvent('button_release_event', obj.fig.canvas, x, y,
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

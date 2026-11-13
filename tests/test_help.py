"""The help text must document every feature, and every example must parse."""

import shlex

import pytest

from viewinsar import cli

pytest.importorskip('mintpy', reason='the examples forward options to view.py')


def _swath_help():
    return cli.create_swath_parser().format_help()


def _flat(text):
    """Collapse whitespace, so a phrase wrapped over two lines still matches."""
    return ' '.join(text.split())


@pytest.mark.parametrize('needle', [
    # 0.2.0 options
    '--nlook', '--topoxsect', '--svg', '--figext', '--dpi', '--figsize-map',
    # the multi-section workflow
    '"=" or "+"', '"n"', '"e"', "A-A'", "B-B'", 'step by step',
    # what gets written
    'output', '_points.csv', '_stats.csv', '_info.txt', 'swath_map_ABC',
    'elevation_median_m', 'dist_km', 'perp_km',
    # behaviour worth knowing
    'northing', 'never "intercept"', 'LAST epoch', 'sea level',
    'not information', 'Cannot be combined with --decimate', 'OWN window',
    # --fault
    '--fault', '-fault', '--fault-label', '--fault-color', '--fault-lw', 'Fault line',
    '.shp', '.geojson', '.gpkg', '.kml', 'any CRS', 'dashed vertical line',
    'crosses each fault trace', '--faultline',
    # the single-dash spellings the options were asked for with
    '-nlook', '-topoxsect', '-svg', 'one dash or two',
    # 0.3.0: the DEM, downloaded or given once
    '-d, --dem [DEM]', '--topoxsect, -topoxsect, --topo [DEM]', '--dem-zoom',
    'Terrain Tiles', 'Mapzen Global Terrain', 'terrarium/{z}/{x}/{y}.png',
    '~/.cache/viewinsar', 'VIEWINSAR_CACHE', 'not --nlook', 'says why in its panel',
    'the exact median', 'DEM covering part of a section', 'please cite them',
    # 0.4.0: radar-coded data, through MintPy's geometry file
    '--lookup', 'geometryRadar.h5', 'radar-coded data', 'RADAR-coded data only',
    '../inputs/', 'nothing downloaded', 'mapped onto the radar grid',
    "each crossing's row/col and lat/lon",
])
def test_swath_help_documents(needle):
    assert needle in _flat(_swath_help()), f'`viewinsar swath -h` does not mention {needle!r}'


def test_top_level_help_points_at_the_new_features(capsys):
    assert cli.main(['--help']) == 0
    out = capsys.readouterr().out
    for needle in ('--nlook', '--topoxsect', '--svg', '--timeseries', '-fault', '-d --topoxsect',
                   '"="', '"n"', '"e"', "A-A'", 'swath', 'list'):
        assert needle in out, f'`viewinsar -h` does not mention {needle!r}'


def test_list_help_explains_the_default_marker():
    with pytest.raises(SystemExit):
        cli.list_main(['-h'])


def test_list_help_text(capsys):
    with pytest.raises(SystemExit):
        cli.list_main(['--help'])
    assert '<- default' in capsys.readouterr().out


def _example_commands():
    """Every `viewinsar swath ...` command in the help epilog, joined across \\."""
    cmds, buf = [], ''
    for raw in cli.EXAMPLE.splitlines():
        line = raw.strip()
        if buf:
            buf += ' ' + line.rstrip('\\').strip()
            if not line.endswith('\\'):
                cmds.append(buf)
                buf = ''
        elif line.startswith('viewinsar swath '):
            if line.endswith('\\'):
                buf = line.rstrip('\\').strip()
            else:
                cmds.append(line)
    return cmds


def test_there_are_examples():
    assert len(_example_commands()) >= 10


@pytest.mark.parametrize('cmd', _example_commands())
def test_every_example_parses(cmd):
    """No example may use an option that neither viewinsar nor view.py knows."""
    from mintpy.cli.view import create_parser as view_parser

    args = shlex.split(cmd)[2:]                       # drop "viewinsar swath"
    _, unknown = cli.create_swath_parser().parse_known_args(args)
    if unknown:
        # whatever viewinsar does not recognise is forwarded to view.py
        _, still = view_parser().parse_known_args(['placeholder.h5', *unknown])
        assert not still, f'example uses unknown option(s) {still}:\n  {cmd}'


def test_help_lines_fit_a_terminal():
    """RawTextHelpFormatter does not wrap, so keep every line readable."""
    for name, text in (('swath -h', _swath_help()), ('-h', cli.top_level_help())):
        long = [ln for ln in text.splitlines() if len(ln) > 92]
        assert not long, f'`viewinsar {name}` has lines over 92 columns:\n' + '\n'.join(long)

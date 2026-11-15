"""The help text must document every feature, and every example must parse."""

import re
import shlex

import pytest

from viewinsar import cli


def _swath_help():
    return cli.create_swath_parser().format_help()


def _flat(text):
    """Collapse whitespace, so a phrase wrapped over two lines still matches -
    and write `-d [DEM], --dem [DEM]` (Python < 3.13) as `-d, --dem [DEM]`."""
    return re.sub(r' \[DEM\],', ',', ' '.join(text.split()))


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
    # 0.5.0: one spelling for every option, --long; the old ones still work
    'written with two dashes', '-fault, -gnss, -vel, -d, -w', 'still work, with a note',
    # 0.3.0: the DEM, downloaded or given once
    '--dem [DEM]', '--topoxsect [DEM]', '--dem-zoom',
    'Terrain Tiles', 'Mapzen Global Terrain', 'terrarium/{z}/{x}/{y}.png',
    # 0.5.3: Copernicus GLO-30 by default
    'Copernicus GLO-30', '--dem-source', 'cop30', 'copernicus-dem-30m', 'EGM2008',
    'The sea has no tiles', 'no GDAL needed', 'overviews', 'the fallback when GLO-30',
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
    for needle in ('--nlook', '--topoxsect', '--svg', '--timeseries', '--fault', '--dem --topoxsect',
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
    view_cli = pytest.importorskip('mintpy.cli.view', reason='checked against view.py')
    view_parser = view_cli.create_parser

    args = shlex.split(cmd, comments=True)[2:]   # drop "viewinsar swath"; a shell drops # comments
    _, unknown = cli.create_swath_parser().parse_known_args(args)
    if unknown:
        # whatever viewinsar does not recognise is forwarded to view.py
        _, still = view_parser().parse_known_args(['placeholder.h5', *unknown])
        assert not still, f'example uses unknown option(s) {still}:\n  {cmd}'


def test_help_lines_fit_a_terminal():
    """RawTextHelpFormatter does not wrap, so keep every line readable."""
    helps = [('swath -h', _swath_help()), ('-h', cli.top_level_help())]
    helps += [(f'{name} -h', cli.command_help(name)) for name in ('stat', 'model', 'tseries')]
    for name, text in helps:
        long = [ln for ln in text.splitlines() if len(ln) > 92]
        assert not long, f'`viewinsar {name}` has lines over 92 columns:\n' + '\n'.join(long)


def _sections_of_top_help():
    """{command: the text of its section} of `viewinsar -h`."""
    text = cli.top_level_help()
    out, name = {}, None
    for line in text.splitlines():
        head = re.match(r'^(swath|stat|model|tseries) - ', line)
        if head:
            name = head.group(1)
            out[name] = []
        elif line.startswith('common options') or line.startswith('install'):
            name = None
        elif name:
            out[name].append(line)
    return {k: '\n'.join(v) for k, v in out.items()}


def test_top_level_help_covers_every_command():
    """`viewinsar -h` has a section on each command, with what it does, its keys,
    its output and examples - and every option it names really exists."""
    parsers = {'swath': cli.create_swath_parser, 'stat': cli.create_stat_parser,
               'model': cli.create_model_parser, 'tseries': cli.create_tseries_parser}
    sections = _sections_of_top_help()
    assert set(sections) == set(parsers)
    for name, text in sections.items():
        for needle in ('keys', 'output', 'examples'):
            assert needle in text, (name, needle)
        known = set(parsers[name]()._option_string_actions)
        named = set(re.findall(r'(?<![\w/+])(--?[a-zA-Z][\w-]*)', text))
        assert named, name
        assert not named - known, f'{name}: unknown option(s) {sorted(named - known)}'
    top = cli.top_level_help()
    for needle in ('--profile', '--line-file', '--rate', '--creep', 'posteriors',
                   'model_sections_ABC.csv', '--save [DIR]', '--radius KM', 'LiCSBAS',
                   'pip install "viewinsar[mintpy]"', '--help-all'):
        assert needle in top, needle


def test_help_all_and_help_command(capsys):
    assert cli.main(['--help-all']) == 0
    out = capsys.readouterr().out
    for name in ('swath', 'stat', 'model', 'tseries', 'list'):
        assert f'usage: viewinsar {name}' in out, name
    assert out.startswith(cli.top_level_help()[:40])
    assert cli.main(['help', 'model']) == 0
    out = capsys.readouterr().out
    assert out.startswith('usage: viewinsar model') and '--profile' in out
    assert cli.main(['help', 'all']) == 0
    assert 'usage: viewinsar tseries' in capsys.readouterr().out

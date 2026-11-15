"""Your own defaults: ~/.config/viewinsar/viewinsar.ini.

One section so far, [data-scale]: multiply the data of a mission - or of the
files whose path matches a pattern - by a factor before anything else (the
map, the profiles, the fit, the statistics, the time series):

    [data-scale]
    # every ALOS-2 product of mine is 10x too large
    ALOS-2 = 0.1
    # ... or by the path (with * and ?)
    */ALOS2_*/geo_velocity.h5 = 0.1

A mission is told by the file's PLATFORM attribute (MintPy: Alos2, Sen ...),
or, without one, by its name and its folders' (alos2, scansar, s1 ...).
--data-scale on the command line comes first; --data-scale 1 switches a rule
off for that run. The file is VIEWINSAR_CONFIG if that is set, else
$XDG_CONFIG_HOME/viewinsar/viewinsar.ini, else ~/.config/viewinsar/viewinsar.ini.
"""

import configparser
import fnmatch
import functools
import os

__all__ = ['config_path', 'scale_rules', 'config_scale', 'mission_of_file',
           'wavelength_note']

# the radar wavelength of each mission, m (L-band ALOS-2 PALSAR-2, C-band S1)
WAVELENGTH = {'ALOS-2': 0.2384, 'Sentinel-1': 0.05547}


def config_path():
    """Where your viewinsar.ini is (or would be)."""
    if os.environ.get('VIEWINSAR_CONFIG'):
        return os.path.expanduser(os.environ['VIEWINSAR_CONFIG'])
    base = os.environ.get('XDG_CONFIG_HOME') or os.path.join(os.path.expanduser('~'), '.config')
    return os.path.join(base, 'viewinsar', 'viewinsar.ini')


@functools.lru_cache(maxsize=4)
def _rules(path, mtime):
    parser = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    parser.optionxform = str                         # paths keep their case
    try:
        parser.read(path, encoding='utf-8')
    except configparser.Error as err:
        print(f'WARNING: {path}: {err} - its rules are not used')
        return ()
    if not parser.has_section('data-scale'):
        return ()
    out = []
    for key, value in parser.items('data-scale'):
        try:
            factor = float(value)
        except ValueError:
            print(f'WARNING: {path}: [data-scale] {key} = {value!r} is not a number - not used')
            continue
        if not factor or factor != factor:
            print(f'WARNING: {path}: [data-scale] {key} = {value} - not used')
            continue
        out.append((key.strip(), factor))
    return tuple(out)


def scale_rules():
    """[(key, factor)] of the [data-scale] section of your viewinsar.ini."""
    path = config_path()
    try:
        mtime = os.path.getmtime(path)
    except OSError:
        return ()
    return _rules(path, mtime)


def _is_pattern(key):
    return any(c in key for c in '*?[/\\')


def mission_of_file(fname, atr=None):
    """(mission, how): 'ALOS-2' / 'Sentinel-1' from the PLATFORM attribute, or
    - for a file without one - from its name and its folders'; (None, ...)
    when it cannot be told."""
    from .gnss import _mission_of, _name_tokens

    platform = str((atr or {}).get('PLATFORM', '') or '').strip()
    if platform:
        flat = platform.lower().replace('-', '').replace('_', '').replace(' ', '')
        mission = _mission_of([flat])
        if mission is None and flat.startswith(('sen', 's1')):
            mission = 'Sentinel-1'
        return mission, f'PLATFORM {platform}'
    found = set()
    for tokens in _name_tokens(fname):
        for i in range(len(tokens)):              # every word, and a LiCSAR frame id
            found.update(m for m in (_mission_of([tokens[i]]), _mission_of(tokens[i:i + 3]))
                         if m)
    if len(found) == 1:
        return found.pop(), 'its name'
    return None, None


def config_scale(fname, atr=None):
    """(factor, why) of the first rule of viewinsar.ini that applies to the
    file - by mission or by path - else (1.0, None)."""
    rules = scale_rules()
    if not rules:
        return 1.0, None
    from .gnss import _mission_of

    path = os.path.abspath(str(fname))
    mission, how = None, None
    for key, factor in rules:
        if _is_pattern(key):
            if fnmatch.fnmatch(path, os.path.expanduser(key)) or fnmatch.fnmatch(
                    os.path.basename(path), key):
                return factor, f'{config_path()}: {key} = {factor:g}'
            continue
        want = _mission_of([key.lower().replace('-', '').replace('_', '')])
        if want is None and key.lower().startswith(('sen', 's1')):
            want = 'Sentinel-1'
        if want is None:
            continue
        if how is None:
            mission, how = mission_of_file(fname, atr)
            how = how or 'unknown'
        if mission == want:
            return factor, f'{config_path()}: {key} = {factor:g}, as {how} says {mission}'
    return 1.0, None


def wavelength_note(fname, atr):
    """A note when the file's WAVELENGTH is not its mission's - velocities
    computed with a wrong one are wrong by the same ratio - else None."""
    try:
        lam = float((atr or {}).get('WAVELENGTH'))
    except (TypeError, ValueError):
        return None
    mission, _ = mission_of_file(fname, atr)
    want = WAVELENGTH.get(mission)
    if want is None or not lam > 0:
        if lam > 1.0:                               # no SAR flies a metre-long wave
            return (f'NOTE: {os.path.basename(fname)} says WAVELENGTH = {lam:g} m, more than '
                    'any SAR: phase turned into distance with it is wrong by as much')
        return None
    ratio = lam / want
    if 0.67 < ratio < 1.5:
        return None
    return (f'NOTE: {os.path.basename(fname)} says WAVELENGTH = {lam:g} m; {mission} has '
            f'{want:g} m. A displacement or velocity computed with it is {ratio:.3g}x what it '
            f'should be: --data-scale {1 / ratio:.3g} shows it right, re-running the '
            'time series with the right WAVELENGTH corrects the file')

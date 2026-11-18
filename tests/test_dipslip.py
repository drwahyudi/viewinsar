"""`viewinsar model reverse|normal`: dip-slip faults.

The forward model against Okada (1992) - a reverse fault dipping 30 deg and a
normal fault dipping 60 deg, each slipping along its dip below its locking
depth - its limits and symmetries, the fit, then the whole command on
synthetic products: a fault striking N10E, the vertical and east velocity of
a reverse fault (8 mm/yr, 12 km, 30 deg towards the east) and of a normal
fault (4 mm/yr, 10 km, 60 deg towards the west), two LOS tracks of the
normal fault, GNSS, the trace.
"""

import json

import numpy as np
import pytest

from test_model import FAST, _attrs, _csv, _write
from viewinsar import cli
from viewinsar import model as mdl

FLON, FLAT, STRIKE = 130.45, 30.30, 10.0
KM_LAT = 111.32
KM_LON = KM_LAT * np.cos(np.radians(FLAT))
NORMAL = np.array([np.cos(np.radians(STRIKE)), -np.sin(np.radians(STRIKE))])  # to the ESE
LINE = ['--start-lalo', '30.30', '130.20', '--end-lalo', '30.30', '130.70', '-w', '10']
REVERSE = dict(S=8.0, d1=12.0, dip=30.0, toward=1)       # dipping east: towards A'
NORMAL_FAULT = dict(S=4.0, d1=10.0, dip=60.0, toward=-1)  # dipping west: towards A


##############################################################################
#                              the forward model                             #
##############################################################################
# Okada (1992), with pyrocko: the surface displacement of a fault 10000 km long
# whose top edge is at d1 and which slips 1 along its dip down to 2000 km
# (rake 90: reverse, -90: normal), at x km from the trace, relative to the
# point above its top edge - the far edges of the fault give an offset only
OKADA = {
    'reverse': dict(d1=12.0, dip=30.0, sense=1,
                    x=[-19.2154, 5.7846, 15.7846, 20.7846, 25.7846, 35.7846, 60.7846],
                    h=[0.4172, 0.2075, 0.0338, 0.0, 0.0133, -0.0133, -0.125],
                    z=[-0.4972, -0.3872, -0.1597, 0.0, 0.0782, 0.051, -0.0088]),
    'normal': dict(d1=10.0, dip=60.0, sense=-1,
                   x=[-34.2265, -9.2265, 0.7735, 5.7735, 10.7735, 20.7735, 45.7735],
                   h=[-0.4273, -0.2717, -0.0646, 0.0, -0.0457, -0.11, -0.0915],
                   z=[0.5707, 0.5048, 0.2687, 0.0, -0.2051, -0.2844, -0.2709]),
}


@pytest.mark.parametrize('name', sorted(OKADA))
def test_dip_slip_matches_okada(name):
    ref = OKADA[name]
    x = np.asarray(ref['x'])
    for kind in ('h', 'z'):
        v = mdl.dip_slip(x, 1.0, ref['d1'], ref['dip'], kind=kind, sense=ref['sense'])
        assert v - v[3] == pytest.approx(np.asarray(ref[kind]), abs=0.012), kind


def test_dip_slip_limits_and_symmetries():
    S, d1, dip = 8.0, 12.0, 30.0
    far = np.array([-1e7, 1e7])
    cos, sin = np.cos(np.radians(dip)), np.sin(np.radians(dip))
    # a reverse fault dipping towards x > 0: the hanging wall there moves
    # S cos(dip) towards the footwall and S sin(dip) up, far off
    h = mdl.dip_slip(far, S, d1, dip, kind='h')
    z = mdl.dip_slip(far, S, d1, dip, kind='z')
    assert h[1] - h[0] == pytest.approx(-S * cos, rel=1e-5)
    assert z[1] - z[0] == pytest.approx(S * sin, rel=1e-5)
    # a normal fault: the other way; dipping towards x < 0: the mirror image
    x = np.linspace(-60, 60, 121)
    for kind in ('h', 'z'):
        assert mdl.dip_slip(x, S, d1, dip, kind=kind, sense=-1) == pytest.approx(
            -mdl.dip_slip(x, S, d1, dip, kind=kind))
    assert mdl.dip_slip(x, S, d1, dip, kind='h', toward=-1) == pytest.approx(
        -mdl.dip_slip(-x, S, d1, dip, kind='h'))
    assert mdl.dip_slip(x, S, d1, dip, kind='z', toward=-1) == pytest.approx(
        mdl.dip_slip(-x, S, d1, dip, kind='z'))
    # the trace elsewhere, an offset, a line of sight: ph h + pu z + a
    assert mdl.dip_slip(x + 3, S, d1, dip, a=2.0, x0=3.0, kind='z') == pytest.approx(
        mdl.dip_slip(x, S, d1, dip, kind='z') + 2.0)
    los = mdl.dip_slip(x, S, d1, dip, kind=(0.6, 0.75))
    assert los == pytest.approx(0.6 * mdl.dip_slip(x, S, d1, dip, kind='h')
                                + 0.75 * mdl.dip_slip(x, S, d1, dip, kind='z'))
    # K parameter sets at once: (K, x.size)
    many = mdl.dip_slip(x, [S, 2 * S], [d1, d1], [dip, 45.0], kind='z')
    assert many.shape == (2, x.size)
    assert many[0] == pytest.approx(mdl.dip_slip(x, S, d1, dip, kind='z'))
    with pytest.raises(ValueError, match='unknown kind'):
        mdl.dip_slip(x, S, d1, dip, kind='q')


def test_dip_slip_model_parameters():
    free = mdl.DipSlipModel(sense=1, kinds=['z', 'h'], dip=(10.0, 45.0), gnss=True,
                            gnss_kind='h', free_x0=(-3.0, 3.0))
    assert free.fault == 'reverse' and free.names == ['S', 'd1', 'dip', 'x0', 'a', 'a2', 'ag']
    assert free.bounds['dip'] == (10.0, 45.0) and free.bounds['S'][0] == 0.0
    fixed = mdl.DipSlipModel(sense=-1, kinds=['z'], dip=30.0)
    assert fixed.fault == 'normal' and fixed.names == ['S', 'd1', 'a'] and fixed.dip_fixed == 30
    theta = np.array([[4.0, 10.0, 1.0], [6.0, 12.0, 0.0]])
    assert fixed.unpack(theta)['dip'] == pytest.approx([30.0, 30.0])
    h, z = fixed.rates(theta)
    assert h == pytest.approx([4 * np.cos(np.radians(30)), 6 * np.cos(np.radians(30))])
    assert z == pytest.approx([4 * np.sin(np.radians(30)), 6 * np.sin(np.radians(30))])
    # the prior: inside the bounds only
    assert np.isfinite(free.log_prior(np.array([5.0, 10.0, 30.0, 0.0, 0.0, 0.0, 0.0])))[0]
    assert not np.isfinite(free.log_prior(np.array([5.0, 10.0, 60.0, 0.0, 0.0, 0.0, 0.0])))[0]
    assert not np.isfinite(free.log_prior(np.array([-1.0, 10.0, 30.0, 0.0, 0.0, 0.0, 0.0])))[0]


def test_fit_dip_slip_recovers_slip_depth_and_dip():
    rng = np.random.default_rng(4)
    x = rng.uniform(-60, 60, 30000)
    truth = dict(S=4.0, d1=10.0, dip=60.0)
    h = mdl.dip_slip(x, **truth, a=1.5, kind='h', sense=-1, toward=-1)
    z = mdl.dip_slip(x, **truth, a=-0.5, kind='z', sense=-1, toward=-1)
    res = mdl.fit_profile([(x, h + rng.normal(0, 1.0, x.size)),
                           (x, z + rng.normal(0, 1.0, x.size))],
                          fault='normal', kinds=['h', 'z'], toward=-1, dip=(5.0, 85.0),
                          walkers=32, steps=1500, seed=2, verbose=False)
    s = res['summary']
    assert res['model'].fault == 'normal' and res['names'] == ['S', 'd1', 'dip', 'a', 'a2']
    assert s['S']['map'] == pytest.approx(4.0, abs=0.6)
    assert s['d1']['map'] == pytest.approx(10.0, abs=2.5)
    assert s['dip']['map'] == pytest.approx(60.0, abs=6.0)
    assert s['a']['map'] == pytest.approx(1.5, abs=0.3)
    assert s['a2']['map'] == pytest.approx(-0.5, abs=0.3)
    for name in ('S', 'd1', 'dip'):
        assert s[name]['lo'] <= truth[name] <= s[name]['hi'], name
    d = res['derived']
    dip = np.radians(s['dip']['map'])
    assert d['h_rate']['map'] == pytest.approx(s['S']['map'] * np.cos(dip), rel=0.02)
    assert d['v_rate']['map'] == pytest.approx(s['S']['map'] * np.sin(dip), rel=0.02)
    assert 'strain' not in d
    # the curve of each component, and its ranges
    xx = np.linspace(-50, 50, 11)
    best, bands = mdl.model_bands(res, xx, offset='a2', kind='z')
    assert best == pytest.approx(mdl.dip_slip(xx, **truth, a=-0.5, kind='z', sense=-1,
                                              toward=-1), abs=0.3)
    assert np.all(bands[95.0][0] <= bands[68.0][0] + 1e-9)


##############################################################################
#                               the products                                 #
##############################################################################
def _frame(lon, lat):
    """km from the fault, perpendicular to it (+ to the ESE)."""
    east, north = (lon - FLON) * KM_LON, (lat - FLAT) * KM_LAT
    return east * NORMAL[0] + north * NORMAL[1]


def _motion(fault, x):
    """(across the fault, vertical) of a reverse or normal fault (mm/yr)."""
    par = REVERSE if fault == 'reverse' else NORMAL_FAULT
    sense = 1 if fault == 'reverse' else -1
    kw = dict(toward=par['toward'], sense=sense)
    return (mdl.dip_slip(x, par['S'], par['d1'], par['dip'], kind='h', **kw),
            mdl.dip_slip(x, par['S'], par['d1'], par['dip'], kind='z', **kw))


@pytest.fixture(scope='module')
def dip_products(tmp_path_factory):
    pytest.importorskip('h5py')
    root = tmp_path_factory.mktemp('dipslip')
    length, width, step = 150, 250, 0.004
    lon = 130.0 + (np.arange(width) + 0.5) * step
    lat = 30.60 - (np.arange(length) + 0.5) * step
    lon2, lat2 = np.meshgrid(lon, lat)
    x = _frame(lon2, lat2)
    rng = np.random.default_rng(7)
    out = {}
    for fault in ('reverse', 'normal'):
        h, z = _motion(fault, x)
        ve = h * NORMAL[0] + 3.0 + rng.normal(0, 0.8, x.shape)
        vu = z - 1.0 + rng.normal(0, 0.8, x.shape)
        out[f'{fault}_ve'] = _write(str(root / fault / 'geo_ve.h5'), {'velocity': ve / 1000},
                                    _attrs(30.60, 130.0, length, width, step))
        out[f'{fault}_vu'] = _write(str(root / fault / 'geo_vu.h5'), {'velocity': vu / 1000},
                                    _attrs(30.60, 130.0, length, width, step))
    # the normal fault seen by two tracks, each LOS with its geometry
    h, z = _motion('normal', x)
    for name, inc0, az, head in (('asc', 39.0, 102.0, -12.0), ('desc', 37.0, -102.0, -168.0)):
        inc = inc0 + 4.0 * (np.arange(width) / width - 0.5)[None, :] * np.ones((length, 1))
        ue = -np.sin(np.radians(inc)) * np.sin(np.radians(az))
        un = np.sin(np.radians(inc)) * np.cos(np.radians(az))
        los = (h * (ue * NORMAL[0] + un * NORMAL[1]) + z * np.cos(np.radians(inc))
               + rng.normal(0, 0.8, x.shape))
        atr = _attrs(30.60, 130.0, length, width, step, HEADING=head,
                     CENTER_INCIDENCE_ANGLE=inc0,
                     ORBIT_DIRECTION='ASCENDING' if abs(head) < 90 else 'DESCENDING')
        out[name] = _write(str(root / name / 'geo_velocity.h5'), {'velocity': los / 1000}, atr)
        _write(str(root / name / 'geo_geometryGeo.h5'),
               {'incidenceAngle': inc, 'azimuthAngle': np.full(inc.shape, az),
                'height': np.full(inc.shape, 300.0)}, dict(atr, FILE_TYPE='geometry', UNIT='1'))
    # GNSS along the line, of the reverse fault: east, north and up
    rows = ['site,lon,lat,ve,vn,vu,se,sn,su']
    for i, lon_g in enumerate(np.linspace(130.22, 130.68, 9)):
        lat_g = 30.30 + 0.01 * np.sin(i)
        hh, zz = _motion('reverse', _frame(lon_g, lat_g))
        rows.append(f'G{i:02d},{lon_g:.4f},{lat_g:.4f},{hh * NORMAL[0] + 2:.3f},'
                    f'{hh * NORMAL[1] - 1:.3f},{zz:.3f},0.4,0.4,1.2')
    out['gnss'] = str(root / 'gnss.csv')
    with open(out['gnss'], 'w') as fid:
        fid.write('\n'.join(rows) + '\n')
    t = np.linspace(-90, 90, 31)
    coords = [[FLON + v * np.sin(np.radians(STRIKE)) / KM_LON,
               FLAT + v * np.cos(np.radians(STRIKE)) / KM_LAT] for v in t]
    out['fault'] = str(root / 'fault.geojson')
    with open(out['fault'], 'w') as fid:
        json.dump({'type': 'FeatureCollection', 'features': [{
            'type': 'Feature', 'properties': {},
            'geometry': {'type': 'LineString', 'coordinates': coords}}]}, fid)
    return out


##############################################################################
#                               the command                                  #
##############################################################################
def test_model_fault_types_options(dip_products):
    vu, ve = dip_products['reverse_vu'], dip_products['reverse_ve']
    q = ['--noverbose']

    def parse(*args):
        return cli.model_cmd_line_parse([*args, *q])

    # the first word: the fault type; strike when it is none
    assert parse(ve, '--comp', 've').fault_type == 'strike'
    assert parse('strike', ve, '--comp', 've').fault_type == 'strike'
    assert parse('thrust', vu).fault_type == 'reverse'
    rev = parse('reverse', vu, '--comp', 'vu')
    assert (rev.fault_type, rev.comps, rev.comp, rev.file) == ('reverse', ['vu'], 'vu', [vu])
    assert rev.dip_range == (5.0, 85.0) and rev.dip_toward == 'end' and rev.gnss_comp == 'vu'
    assert rev.argv[0] == 'reverse'
    # one --comp per file: the files after them stay inputs
    nor = parse('normal', '--comp', 've', 'vu', ve, vu, '--dip-toward', 'start')
    assert nor.file == [ve, vu] and nor.comps == ['ve', 'vu'] and nor.comp == 'mixed'
    assert nor.gnss_comp == 'perp' and nor.dip_toward == 'start'
    assert parse('normal', ve, vu, '--comp', 've,vu').comps == ['ve', 'vu']
    assert parse('reverse', ve, vu, '--comp', 'vu').comps == ['vu', 'vu']
    # --dip: fixed, or a range (either order); the files after it stay inputs
    assert parse('reverse', '--dip', '30', vu).dip_range == (30.0, 30.0)
    assert parse('reverse', vu, '--dip', '45', '10').dip_range == (10.0, 45.0)
    for bad in ('0', '95'):
        with pytest.raises(ValueError, match='--dip'):
            parse('reverse', vu, '--dip', bad)
    # what belongs to the other kind of fault
    with pytest.raises(ValueError, match='reverse or a normal'):
        parse(vu, '--comp', 'vu')
    with pytest.raises(ValueError, match='ONE --comp'):
        parse(ve, vu, '--comp', 've', 'los')
    with pytest.raises(ValueError, match='--dip'):
        parse(ve, '--dip', '30')
    for opt in (['--creep'], ['--creep-offset'], ['--rot', '0.02'], ['--strain']):
        with pytest.raises(ValueError, match='strike-slip'):
            parse('reverse', vu, *opt)
    with pytest.raises(ValueError, match='--comp par'):
        parse('reverse', vu, '--comp', 'par')
    with pytest.raises(ValueError, match='components for'):
        parse('normal', ve, vu, vu, '--comp', 've', 'vu')
    with pytest.raises(ValueError, match='not a component'):
        parse('normal', ve, '--comp', 'up')
    with pytest.raises(ValueError, match='comes first'):
        parse(vu, 'reverse')
    # the GNSS: across the fault (perp) or vertical (vu)
    with pytest.raises(ValueError, match='perp'):
        parse('reverse', vu, '--gnss', dip_products['gnss'], '--vel', 'par')
    with pytest.raises(ValueError, match='--vel perp'):
        parse(ve, '--comp', 've', '--gnss', dip_products['gnss'], '--vel', 'perp')
    assert parse('reverse', vu, '--comp', 'vu', '--gnss', dip_products['gnss'], '--vel',
                 'perp').gnss_comp == 'perp'


def test_model_dip_slip_names_the_data(dip_products, tmp_path):
    """--comp auto for a reverse or normal fault: by the name of each file
    (up / vu / vertical, east / ve, north / vn), else its line of sight."""
    import shutil

    from viewinsar.modelview import ModelViewer

    hz = tmp_path / 'geo_hz.h5'
    shutil.copy(dip_products['reverse_ve'], hz)
    plain = tmp_path / 'geo_rate.h5'
    shutil.copy(dip_products['reverse_ve'], plain)

    def resolved(*args):
        inps = cli.model_cmd_line_parse([*args, '--noverbose'])
        inps.disp_fig = False
        obj = ModelViewer(inps, cli.get_view_cmd(inps))
        try:
            obj.open()
            return obj.comps, obj.comp, obj.gnss_comp
        finally:
            obj.close()

    assert resolved('reverse', dip_products['reverse_vu'], dip_products['reverse_ve']) == \
        (['vu', 've'], 'mixed', 'perp')
    assert resolved('normal', dip_products['asc'], dip_products['desc'])[:2] == \
        (['los', 'los'], 'los')
    assert resolved('reverse', dip_products['reverse_vu'])[2] == 'vu'
    with pytest.raises(ValueError, match='ve or vn'):
        resolved('reverse', str(hz))
    with pytest.raises(ValueError, match='--comp vu, ve or vn'):
        resolved('reverse', str(plain))
    # a strike-slip fault and a file named as the vertical: said
    with pytest.raises(ValueError, match='reverse or a normal'):
        resolved(dip_products['reverse_vu'])


def _panels(sec):
    return [ax.get_ylabel() for ax in sec.ax_panels]


def test_model_reverse_fault_from_the_vertical(dip_products, tmp_path):
    """The vertical velocity of a reverse fault and the GNSS across it: two
    panels - the InSAR vertical, the GNSS across the fault - the model in
    each, the dip and the far-field rates written out."""
    from matplotlib.collections import PathCollection
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.text import Annotation

    out = tmp_path / 'out'
    obj = cli.model_main(['reverse', dip_products['reverse_vu'], '--comp', 'vu', '--gnss',
                          dip_products['gnss'], '--vel', 'perp', '-fault', dip_products['fault'],
                          *LINE, *FAST, '--nodisplay', '--noverbose', '--outdir', str(out)])
    sec = obj.sections[0]
    res = sec.result
    s = res['summary']
    assert res['model'].fault == 'reverse' and res['model'].kinds == ['z']
    assert res['model'].gnss_kind == 'h' and res['names'] == ['S', 'd1', 'dip', 'a', 'ag']
    assert s['S']['map'] == pytest.approx(REVERSE['S'], abs=1.0)
    assert s['d1']['map'] == pytest.approx(REVERSE['d1'], abs=2.5)
    assert s['dip']['map'] == pytest.approx(REVERSE['dip'], abs=4.0)
    # the strike by the right-hand rule: the fault dips to its right (east)
    assert sec.strike == pytest.approx(STRIKE, abs=1.0)
    # a panel for the InSAR vertical, one for the GNSS across the fault
    assert _panels(sec) == ['Fault vertical velocity\n(mm/yr)', 'Fault-normal velocity\n(mm/yr)']
    assert [spec['kind'] for spec in sec.specs] == ['z', 'h']
    gnss_ax = sec.ax_panels[1]                     # the GNSS only: no pixel, no median
    assert len(gnss_ax.containers) == 1
    assert not [c for c in gnss_ax.collections if type(c) is PathCollection]
    assert [c for c in sec.ax_swath.collections if type(c) is PathCollection]
    texts = [a for a in sec.ax_swath.artists if isinstance(a, AnchoredText)]
    lines = texts[0].txt._text.get_text().split('\n')
    assert lines[0].startswith('slip rate = ') and lines[0].endswith(', reverse')
    assert lines[2].startswith('dip = ') and "°" in lines[2] and lines[2].endswith("towards A'")
    assert lines[-2].startswith('shortening rate = ') and lines[-1].startswith('uplift rate = ')
    # A and A' over the top panel, the title; the fault dotted in each panel
    top = {a.get_text() for a in sec.ax_swath.texts if isinstance(a, Annotation)}
    assert {'A', "A'"} <= top
    assert sec.ax_swath.get_title().startswith("A-A': swath ±10 km")
    assert 'reverse fault dipping towards A\'' in sec.ax_swath.get_title()
    for ax in sec.ax_panels:
        assert [ln for ln in ax.get_lines() if ln.get_linestyle() == ':']
        assert [ln for ln in ax.get_lines() if ln.get_linestyle() == '--']   # the MAP
    # the files: the rates, the samples, the GNSS
    params = _csv(next(out.glob('*_params.csv')))
    rows = {str(n): r for n, r in zip(params['name'], params)}
    assert {'S', 'd1', 'dip', 'h_rate', 'v_rate'} <= set(rows) and 'strain' not in rows
    dip = np.radians(rows['dip']['map'])
    assert rows['h_rate']['map'] == pytest.approx(rows['S']['map'] * np.cos(dip), rel=1e-3)
    with open(next(out.glob('*_samples.csv'))) as fid:
        assert fid.readline().strip() == 'S,d1,dip,a,ag,h_rate,v_rate'
    bins = _csv(next(out.glob('*_bins.csv')))
    assert np.median(np.abs(bins['residual'])) < 0.6
    report = next(out.glob('*_params.txt')).read_text()
    for needle in ('fault type     : reverse - dipping towards A\'', 'Savage 1983',
                   'priors         : uniform, S 0 to', 'dip 5 to 85', '--comp vu',
                   'right-hand rule'):
        assert needle in report, needle
    assert next(out.glob('*_gnss.csv')).exists()


def test_model_normal_fault_east_and_vertical(dip_products, tmp_path):
    """The east and the vertical velocity of a normal fault, its dip fixed,
    two sections: a panel each, the files first in the legend, the table of
    the sections with the dip and the rates."""
    out = tmp_path / 'out'
    obj = cli.model_main(['normal', '--comp', 've', 'vu', dip_products['normal_ve'],
                          dip_products['normal_vu'], '--labels', 'East', 'Up', '--dip', '60',
                          '--dip-toward', 'start', '-fault', dip_products['fault'],
                          '--profile', '30.30', '130.20', '30.30', '130.70',
                          '--profile', '30.40', '130.22', '30.40', '130.72', '-w', '10',
                          *FAST, '--nodisplay', '--noverbose', '--outdir', str(out)])
    assert len(obj.sections) == 2
    sec = obj.sections[0]
    res = sec.result
    s = res['summary']
    assert res['names'] == ['S', 'd1', 'a', 'a2'] and res['model'].kinds == ['h', 'z']
    assert s['S']['map'] == pytest.approx(NORMAL_FAULT['S'], abs=0.6)
    assert s['d1']['map'] == pytest.approx(NORMAL_FAULT['d1'], abs=2.5)
    assert sec.strike == pytest.approx(STRIKE + 180.0, abs=1.0)     # dipping west
    assert _panels(sec) == ['Fault-normal velocity\n(mm/yr)', 'Fault vertical velocity\n(mm/yr)']
    legend = sec.ax_swath.get_legend()
    assert [t.get_text() for t in legend.get_texts()][:2] == ['East', 'Up']
    lines = sec.parameter_lines()
    assert lines[0].endswith(', normal') and lines[2] == 'dip = 60° (fixed), towards A'
    assert lines[-2].startswith('extension rate = ') and lines[-1].startswith('subsidence rate')
    # the vertical bins sit on the model, each brought onto its own panel's offset
    bins = _csv(next(out.glob('model_A_*_bins.csv')))
    assert set(bins['file_idx']) == {0, 1}
    assert np.median(np.abs(bins['residual'])) < 0.6
    assert np.allclose(bins['median_shifted'], bins['median'])         # one file a panel
    params = _csv(next(out.glob('model_A_*_params.csv')))
    rows = {str(n): r for n, r in zip(params['name'], params)}
    assert rows['dip']['map'] == 60.0 and rows['dip']['lower_95'] == 60.0
    table = _csv(next(out.glob('model_sections_AB.csv')))
    assert list(table['fault_type']) == ['normal', 'normal']
    assert list(table['dip_map']) == [60.0, 60.0]
    assert 'h_rate_map' in table.dtype.names and 'strain_map' not in table.dtype.names
    assert table['S_map'] == pytest.approx([NORMAL_FAULT['S']] * 2, abs=0.6)


def test_model_dip_slip_lines_of_sight(dip_products, tmp_path):
    """Two tracks of a normal fault, each LOS as it is - a panel each, the
    model seen through each line of sight; drawn as a density."""
    obj = cli.model_main(['normal', dip_products['asc'], dip_products['desc'], '--labels', 'asc',
                          'desc', '--dip-toward', 'start', '-fault', dip_products['fault'],
                          *LINE, '--cmap-prof', 'magma', *FAST, '--nodisplay', '--noverbose',
                          '--outdir', str(tmp_path)])
    sec = obj.sections[0]
    res = sec.result
    s = res['summary']
    kinds = res['model'].kinds
    assert all(isinstance(k, tuple) for k in kinds)
    # ascending and descending both look east-ish for an N-S fault: ph opposite
    assert kinds[0][0] * kinds[1][0] < 0 and all(0.7 < k[1] < 0.85 for k in kinds)
    assert s['S']['map'] == pytest.approx(NORMAL_FAULT['S'], abs=0.6)
    assert s['dip']['map'] == pytest.approx(NORMAL_FAULT['dip'], abs=6.0)
    assert _panels(sec) == ['LOS velocity, asc\n(mm/yr)', 'LOS velocity, desc\n(mm/yr)']
    assert all(cax.axison for cax in sec.ax_cbars)                    # a density each
    report = next(tmp_path.glob('*_params.txt')).read_text()
    assert report.count('--comp los: the LOS velocity, as it is - it sees') == 2


def test_model_help_mentions_fault_types(capsys):
    with pytest.raises(SystemExit):
        cli.main(['model', 'reverse', '--help'])
    out = capsys.readouterr().out
    assert out.startswith('usage: viewinsar model [strike|reverse|normal]')
    for needle in ('--dip DEG', '--dip-toward', 'Savage (1983)', 'perp', 'reverse geo_vu.h5',
                   '--comp ve vu', 'h_rate'):
        assert needle in out, needle
    cli.main(['--help'])
    assert 'model [strike|reverse|normal]' in capsys.readouterr().out


def test_model_dip_slip_draws_every_panel_once_again(dip_products, tmp_path):
    """The "l" key on a figure of two panels: the legend moves, the fit stays."""
    from matplotlib.backend_bases import KeyEvent

    from viewinsar.modelview import ModelViewer

    inps = cli.model_cmd_line_parse(['normal', '--comp', 've', 'vu', dip_products['normal_ve'],
                                     dip_products['normal_vu'], '--dip-toward', 'start',
                                     '-fault', dip_products['fault'], *LINE, *FAST,
                                     '--noverbose', '--outdir', str(tmp_path)])
    inps.disp_fig = False
    obj = ModelViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    sec = obj.sections[0]
    result = sec.result
    sec.fig.canvas.draw()
    seen = [sec.legend_at]
    for _ in range(3):
        x, y = sec.ax_swath.transAxes.transform((0.5, 0.5))
        sec.fig.canvas.callbacks.process('key_press_event', KeyEvent(
            'key_press_event', sec.fig.canvas, key='l', x=x, y=y))
        seen.append(sec.legend_at)
    assert sorted(seen[1:]) == ['in', 'off', 'out'] and sec.result is result
    assert len(sec.ax_panels) == 2 and all(ax.lines for ax in sec.ax_panels)
    obj.close()

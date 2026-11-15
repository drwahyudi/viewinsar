"""--data-scale and your viewinsar.ini: a product whose values are a factor
off (an ALOS-2 product 10x too large) - multiplied before anything, on the
map, in the profiles, the fit, the statistics and the time series."""

import os

import numpy as np
import pytest

from viewinsar import cli
from viewinsar import settings as st


def write_ini(tmp_path, monkeypatch, text):
    path = tmp_path / 'viewinsar.ini'
    path.write_text(text)
    monkeypatch.setenv('VIEWINSAR_CONFIG', str(path))
    return str(path)


def write_velocity(path, platform=None, wavelength=None, scale=10.0):
    """A MintPy geo_velocity.h5: an arctan step of +/-10 mm/yr, in m/year,
    times `scale` (an ALOS-2 product 10x too large)."""
    h5py = pytest.importorskip('h5py')
    length, width, step = 120, 150, 0.004
    lon = 130.15 + (np.arange(width) + 0.5) * step
    vel = np.tile(20.0 / np.pi * np.arctan((lon - 130.45) * 111.32 / 5.0), (length, 1))
    os.makedirs(os.path.dirname(str(path)) or '.', exist_ok=True)
    with h5py.File(str(path), 'w') as fid:
        fid.create_dataset('velocity', data=(vel * scale / 1000.0).astype(np.float32))
        atr = {'FILE_TYPE': 'velocity', 'UNIT': 'm/year', 'LENGTH': str(length),
               'WIDTH': str(width), 'Y_FIRST': '30.50', 'X_FIRST': '130.15',
               'Y_STEP': str(-step), 'X_STEP': str(step), 'Y_UNIT': 'degrees',
               'X_UNIT': 'degrees', 'REF_Y': '60', 'REF_X': '75'}
        if platform:
            atr['PLATFORM'] = platform
        if wavelength:
            atr['WAVELENGTH'] = str(wavelength)
        fid.attrs.update(atr)
    return str(path)


SWATH = ['--start-lalo', '30.30', '130.20', '--end-lalo', '30.30', '130.70', '--width', '5',
         '--nodisplay', '--noverbose']


def _swath(path, tmp_path, *extra):
    obj = cli.swath_main([path, *SWATH, '--outdir', str(tmp_path / 'out'), *extra])
    obj.close()
    return obj


def _range_mm(obj):
    """The 2-98% range of the profile values, in mm/yr."""
    from viewinsar.viewer import velocity_factor

    vals = obj.data_list[0][np.isfinite(obj.data_list[0])]
    lo, hi = np.percentile(vals, [2, 98])
    return (hi - lo) * velocity_factor(obj.disp_unit, 'mm/yr')


def test_data_scale_divides_everything(tmp_path):
    vel = write_velocity(tmp_path / 'alos2' / 'geo_velocity.h5')      # 10x too large
    big = _swath(vel, tmp_path)
    fixed = _swath(vel, tmp_path, '--data-scale', '0.1')
    assert _range_mm(big) == pytest.approx(10 * _range_mm(fixed), rel=1e-3)
    assert 15 < _range_mm(fixed) < 21                                 # +/-10 mm/yr: right
    # the map and the profiles alike
    assert np.nanmax(np.abs(fixed.data_img)) == pytest.approx(
        np.nanmax(np.abs(big.data_img)) / 10, rel=1e-4)


def test_data_scale_checks(tmp_path):
    vel = write_velocity(tmp_path / 'geo_velocity.h5')
    for bad in ('0', 'nan'):
        with pytest.raises(ValueError, match='no factor'):
            cli.cmd_line_parse([vel, '--data-scale', bad, '--noverbose'])
    with pytest.raises(ValueError, match='one per file'):
        cli.cmd_line_parse([vel, vel, vel, '--data-scale', '0.1', '1', '--noverbose'])


def test_mission_of_file(tmp_path):
    assert st.mission_of_file('x.h5', {'PLATFORM': 'Alos2'}) == ('ALOS-2', 'PLATFORM Alos2')
    assert st.mission_of_file('x.h5', {'PLATFORM': 'Sen'})[0] == 'Sentinel-1'
    assert st.mission_of_file(str(tmp_path / 'ALOS2_desc' / 'velocity.h5'))[0] == 'ALOS-2'
    # a platform it does not know is not second-guessed from the name
    assert st.mission_of_file(str(tmp_path / 'alos2' / 'v.h5'), {'PLATFORM': 'TSX'})[0] is None
    # a name that says two missions says none
    assert st.mission_of_file(str(tmp_path / 'alos2_vs_s1' / 'velocity.h5'))[0] is None


def test_ini_rule_by_mission(tmp_path, monkeypatch, capsys):
    ini = write_ini(tmp_path, monkeypatch, '[data-scale]\n# mine are 10x too large\n'
                                           'ALOS-2 = 0.1\n')
    alos = write_velocity(tmp_path / 'a' / 'geo_velocity.h5', platform='Alos2')
    sen = write_velocity(tmp_path / 's' / 'geo_velocity.h5', platform='Sen', scale=1.0)
    assert st.config_scale(alos, {'PLATFORM': 'Alos2'})[0] == 0.1
    obj = _swath(alos, tmp_path)
    assert 15 < _range_mm(obj) < 21
    out = capsys.readouterr().out
    assert 'its values x 0.1' in out and ini in out and 'PLATFORM Alos2' in out
    assert 15 < _range_mm(_swath(sen, tmp_path)) < 21                   # untouched
    # --data-scale comes first: 1 switches the rule off for a run
    assert _range_mm(_swath(alos, tmp_path, '--data-scale', '1')) > 150


def test_ini_rule_by_path(tmp_path, monkeypatch):
    write_ini(tmp_path, monkeypatch, '[data-scale]\n*/T125D_alos/*.h5 = 0.1\nbad = x\n')
    vel = write_velocity(tmp_path / 'T125D_alos' / 'geo_velocity.h5')
    assert st.config_scale(vel)[0] == 0.1
    other = write_velocity(tmp_path / 'T087A' / 'geo_velocity.h5', scale=1.0)
    assert st.config_scale(other) == (1.0, None)


def test_wavelength_note():
    note = st.wavelength_note('geo_velocity.h5', {'PLATFORM': 'Alos2', 'WAVELENGTH': '2.384'})
    assert '10x' in note and '--data-scale 0.1' in note
    assert st.wavelength_note('v.h5', {'PLATFORM': 'Alos2', 'WAVELENGTH': '0.2424'}) is None
    assert st.wavelength_note('v.h5', {'PLATFORM': 'Sen', 'WAVELENGTH': '0.05546576'}) is None
    assert 'more than any SAR' in st.wavelength_note('v.h5', {'WAVELENGTH': '2.4'})


def test_model_with_data_scale(tmp_path):
    """The fit is on the scaled data: S of a product 10x too large, x 0.1."""
    import test_licsbas as tl

    vel, fault = tl.licsbas_tree(str(tmp_path / 'work'))
    import rasterio

    with rasterio.open(vel) as src:
        prof, data = src.profile, src.read(1)
    big = tmp_path / 'alos2_vel_mskd.tif'
    with rasterio.open(str(big), 'w', **prof) as dst:
        dst.write((data * 10).astype(np.float32), 1)
    args = [str(big), '--fault', fault, *tl.PROFILE, '--geom',
            str(tmp_path / 'work' / 'GEOCml10GACOS'), *tl.FAST, '--nodisplay', '--noverbose',
            '--outdir', str(tmp_path / 'out')]
    obj = cli.model_main(args + ['--data-scale', '0.1'])
    assert obj.comp == 'los'
    assert obj.sections[0].result['summary']['S']['map'] == pytest.approx(-tl.S_TRUE, abs=1.5)


def test_tseries_data_scale(tmp_path, monkeypatch):
    from test_tseries import write_mintpy
    from viewinsar.tsdata import TimeSeries

    path = write_mintpy(str(tmp_path / 'mintpy'))
    plain = TimeSeries(path, verbose=False)
    scaled = TimeSeries(path, verbose=False, scale=0.1)
    assert scaled.to_mm == pytest.approx(plain.to_mm * 0.1)
    assert scaled.velocity() == pytest.approx(plain.velocity() * 0.1, abs=1e-6)
    assert scaled.slice(-1) == pytest.approx(plain.slice(-1) * 0.1, abs=1e-6)
    plain.close(), scaled.close()
    # ... and by a rule of the ini, by the path
    write_ini(tmp_path, monkeypatch, f'[data-scale]\n{tmp_path}/mintpy/*.h5 = 0.5\n')
    ruled = TimeSeries(path, verbose=False)
    assert ruled.to_mm == pytest.approx(plain.to_mm * 0.5)
    ruled.close()

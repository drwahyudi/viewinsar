"""`viewinsar tseries`: LiCSBAS_plot_ts.py with circles.

The synthetic time series are exact: every pixel moves at a rate that grows
linearly from west to east, plus a signal common to the whole scene at each
date (the "atmosphere"). Relative to the reference area and the reference
date, a pixel's series is then (rate - rate_ref) * (t - t_ref), with no
atmosphere left, and the median and quartiles of a circle follow from the
rates of its pixels.
"""

import datetime as dt
import json
import os

import numpy as np
import pytest

from viewinsar import tsdata
from viewinsar.tsdata import (
    TimeSeries,
    circle_pixels_rc,
    linear_velocity,
    model_curve,
    parse_geo_range,
    parse_range,
    relative_series,
    series_stats,
)

h5py = pytest.importorskip('h5py')

LENGTH, WIDTH, N_DATE = 60, 80, 30
LAT1, LON1, DLAT, DLON = -2.0, 121.0, -0.002, 0.002        # LiCSBAS: centre of pixel 0
DATES = [dt.date(2018, 1, 1) + dt.timedelta(days=24 * i) for i in range(N_DATE)]
YEARS = np.array([(d - DATES[0]).days for d in DATES]) / 365.25
REF = (60, 65, 10, 15)                                      # x1, x2, y1, y2


def rate_field():
    """mm/yr: 0.5 * (x - 40), so the reference area (x 60-64) moves at 11."""
    return np.broadcast_to(0.5 * (np.arange(WIDTH) - 40.0), (LENGTH, WIDTH)).copy()


def common_signal():
    rng = np.random.default_rng(1)
    sig = rng.normal(0, 8, N_DATE)
    return sig - sig[0]


def make_cube():
    cube = rate_field()[None] * YEARS[:, None, None] + common_signal()[:, None, None]
    cube[:, :, :5] = np.nan                                 # no data west of x = 5
    return cube.astype(np.float32)


def write_licsbas(folder, name='cum.h5', filtered=False, gap=True):
    os.makedirs(os.path.join(folder, 'results'), exist_ok=True)
    cube = make_cube()
    with h5py.File(os.path.join(folder, name), 'w') as f:
        f.create_dataset('cum', data=cube)
        f.create_dataset('vel', data=rate_field().astype(np.float32))
        f.create_dataset('imdates', data=np.array([int(d.strftime('%Y%m%d')) for d in DATES]))
        f.create_dataset('refarea', data='{}:{}/{}:{}'.format(*REF))
        for key, val in (('corner_lat', LAT1), ('corner_lon', LON1), ('post_lat', DLAT),
                         ('post_lon', DLON)):
            f.create_dataset(key, data=val)
        if gap:
            g = np.zeros((N_DATE - 1, LENGTH, WIDTH), dtype=np.int8)
            g[10, :, 40:] = 1                               # a network gap, east half
            f.create_dataset('gap', data=g)
        if filtered:
            f.create_dataset('filtwidth_km', data=2.0)
            f.create_dataset('filtwidth_yr', data=0.066)
            f.create_dataset('deramp_flag', data='')
    mask = np.ones((LENGTH, WIDTH), np.float32)
    mask[40:, :] = 0                                        # the south masked
    layers = dict(mask=mask, coh_avg=np.full((LENGTH, WIDTH), 0.3),
                  n_unw=np.full((LENGTH, WIDTH), 450.0), hgt=np.full((LENGTH, WIDTH), 800.0))
    for key, arr in layers.items():
        np.asarray(arr, dtype='<f4').tofile(os.path.join(folder, 'results', key))
    return os.path.join(folder, name)


@pytest.fixture
def licsbas(tmp_path):
    folder = str(tmp_path / 'TS_GEOCml10')
    write_licsbas(folder, 'cum.h5')
    return write_licsbas(folder, 'cum_filt.h5', filtered=True)


def expected_series(rows, cols, ix_ref=0, ref=REF):
    """The exact relative series of these pixels: (n_date, n_px)."""
    rates = rate_field()[rows, cols]
    x1, x2, y1, y2 = ref
    ref_rate = rate_field()[y1:y2, x1:x2].mean()
    t = YEARS - YEARS[ix_ref]
    return (rates - ref_rate)[None, :] * t[:, None]


##############################################################################
#                                 the core                                   #
##############################################################################
def test_ranges():
    assert parse_range('140:145/125:130', 200, 150) == (140, 145, 125, 130)
    assert parse_range('0:0/0:0', 200, 150) == (0, 200, 0, 150)          # LiCSBAS: all
    for bad in ('140:145', '145:140/1:2', '0:300/0:10', 'a:b/c:d'):
        with pytest.raises(ValueError):
            parse_range(bad, 200, 150)
    assert parse_geo_range('121.1/121.2/-2.1/-2.0') == (121.1, 121.2, -2.1, -2.0)
    with pytest.raises(ValueError):
        parse_geo_range('121.1/121.2')


def test_series_relative_to_the_reference():
    values = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 9.0]])          # 3 dates, 2 pixels
    ref = np.array([1.0, 2.0, 3.0])
    rel = relative_series(values, ref, 0)
    assert rel[0] == pytest.approx([0, 0])
    assert rel[:, 1] == pytest.approx([0, 3, 5])
    assert relative_series(values, ref, 1)[1] == pytest.approx([0, 0])


def test_series_stats():
    values = np.array([[0, 0, 0, 0, np.nan], [1, 2, 3, 4, 100]], dtype=float)
    out = series_stats(values)
    assert out['median'] == pytest.approx([0, 3])
    assert out['q25'] == pytest.approx([0, 2]) and out['q75'] == pytest.approx([0, 4])
    assert list(out['count']) == [4, 5]
    one = series_stats(values[:, :1])                                    # one pixel
    assert one['median'] == pytest.approx(one['q25'])
    assert one['median'] == pytest.approx(one['q75'])
    empty = series_stats(np.full((2, 3), np.nan))
    assert np.isnan(empty['median']).all() and list(empty['count']) == [0, 0]


def test_fit_models_recover_a_known_signal():
    t = np.linspace(0, 5, 80)
    y = 3.0 - 12.0 * t + 4.0 * np.sin(2 * np.pi * t + 0.7)
    out = np.linspace(0, 5, 11)
    want = 3.0 - 12.0 * out + 4.0 * np.sin(2 * np.pi * out + 0.7)
    assert model_curve(t, y, out, 1) == pytest.approx(want, abs=1e-8)        # Annual+L
    assert model_curve(t, y, out, 3) == pytest.approx(want, abs=1e-8)        # Annual+Q
    lin = model_curve(t, 2.0 + 5.0 * t, out, 0)
    assert lin == pytest.approx(2.0 + 5.0 * out)
    y[::3] = np.nan                                                           # gaps: dropped
    assert model_curve(t, y, out, 1) == pytest.approx(want, abs=1e-8)
    assert np.isnan(model_curve(t[:2], y[:2], out, 3)).all()                  # too few


def test_linear_velocity_is_nan_aware():
    cube = rate_field()[None] * YEARS[:, None, None] + 7.0
    cube[3:9, 0, 0] = np.nan
    cube[:, 1, 1] = np.nan
    vel = linear_velocity(cube, YEARS, max_values=5000)                       # many bands
    assert vel[2:] == pytest.approx(rate_field()[2:])
    assert vel[0, 0] == pytest.approx(rate_field()[0, 0])
    assert np.isnan(vel[1, 1])


##############################################################################
#                                 reading                                    #
##############################################################################
def test_licsbas_file(licsbas):
    ts = TimeSeries(licsbas, verbose=False)
    assert ts.kind == 'LiCSBAS' and ts.to_mm == 1.0
    assert (ts.n_date, ts.length, ts.width) == (N_DATE, LENGTH, WIDTH)
    assert ts.dates[0] == '20180101' and ts.years[-1] == pytest.approx(YEARS[-1])
    assert ts.ref_area == REF
    # LiCSBAS's corner_lat/lon is the centre of the first pixel
    assert ts.yx2lalo(0, 0) == pytest.approx((LAT1, LON1))
    assert ts.yx2lalo(10, 20) == pytest.approx((LAT1 + 10 * DLAT, LON1 + 20 * DLON))
    assert ts.lalo2yx(LAT1 + 10 * DLAT, LON1 + 20 * DLON) == (10, 20)
    x_km, y_km = ts.pixel_size_km(30)
    assert y_km == pytest.approx(0.2224, rel=0.01) and x_km == pytest.approx(0.2222, rel=0.01)
    assert ts.label.startswith('1: s=2.0km, t=0.07yr (24d)')
    assert set(ts.layers) == {'mask', 'coh_avg', 'n_unw', 'hgt'}
    assert ts.keep[39].all() and not ts.keep[40].any()
    assert ts.velocity() == pytest.approx(rate_field())
    assert ts.gap_fraction(np.array([0, 0]), np.array([10, 50]))[10] == pytest.approx(0.5)
    ts.close()


def write_mintpy(folder, geocoded=True, with_velocity=True):
    os.makedirs(os.path.join(folder, 'inputs'), exist_ok=True)
    ref_y, ref_x = 12, 62
    rates = rate_field()
    cube = (rates - rates[ref_y, ref_x])[None] * YEARS[:, None, None] / 1000.0   # metres
    atr = {'FILE_TYPE': 'timeseries', 'UNIT': 'm', 'LENGTH': str(LENGTH), 'WIDTH': str(WIDTH),
           'REF_Y': str(ref_y), 'REF_X': str(ref_x), 'REF_DATE': DATES[2].strftime('%Y%m%d')}
    if geocoded:
        atr.update(Y_FIRST=str(LAT1 - DLAT / 2), X_FIRST=str(LON1 - DLON / 2),
                   Y_STEP=str(DLAT), X_STEP=str(DLON), Y_UNIT='degrees', X_UNIT='degrees')

    def write(name, dsets, extra=None):
        with h5py.File(os.path.join(folder, name), 'w') as f:
            for key, val in dsets.items():
                f.create_dataset(key, data=val)
            for key, val in {**atr, **(extra or {})}.items():
                f.attrs[key] = val

    write('timeseries.h5', {'timeseries': cube.astype(np.float32),
                            'date': np.array([d.strftime('%Y%m%d') for d in DATES], dtype='S8')})
    if with_velocity:
        write('velocity.h5', {'velocity': ((rates - rates[ref_y, ref_x]) / 1000.0)
                               .astype(np.float32)}, {'FILE_TYPE': 'velocity', 'UNIT': 'm/year'})
    mask = np.ones((LENGTH, WIDTH), bool)
    mask[:, 70:] = False
    write('maskTempCoh.h5', {'mask': mask}, {'FILE_TYPE': 'mask', 'UNIT': '1'})
    write('temporalCoherence.h5', {'temporalCoherence': np.full((LENGTH, WIDTH), 0.9,
                                                                np.float32)})
    write(os.path.join('inputs', 'geometryGeo.h5'), {'height': np.full((LENGTH, WIDTH), 120.0,
                                                                       np.float32)})
    return os.path.join(folder, 'timeseries.h5')


def test_mintpy_file(tmp_path):
    path = write_mintpy(str(tmp_path / 'mintpy'))
    ts = TimeSeries(path, verbose=False)
    assert ts.kind == 'MintPy' and ts.to_mm == 1000.0
    assert ts.ref_area == (62, 63, 12, 13) and ts.ref_date == DATES[2].strftime('%Y%m%d')
    assert list(ts.layers) == ['mask', 'temporalCoherence', 'hgt']
    assert ts.keep[:, :70].all() and not ts.keep[:, 70:].any()
    assert ts.vel_source.endswith('velocity.h5')
    assert ts.velocity() == pytest.approx(rate_field() - rate_field()[12, 62], abs=1e-4)
    assert ts.slice(-1)[0, 0] == pytest.approx((rate_field()[0, 0] - 11) * YEARS[-1], rel=1e-5)
    ts.close()
    # no velocity.h5: a linear fit of the time series
    path = write_mintpy(str(tmp_path / 'novel'), with_velocity=False)
    ts = TimeSeries(path, verbose=False)
    assert ts.velocity() == pytest.approx(rate_field() - rate_field()[12, 62], abs=1e-3)
    assert ts.vel_source == 'linear fit of the time series'
    ts.close()


def test_hdfeos5_file(tmp_path):
    path = str(tmp_path / 'S1_test.he5')
    with h5py.File(path, 'w') as f:
        obs = f.create_group('HDFEOS/GRIDS/timeseries/observation')
        obs.create_dataset('displacement', data=make_cube() / 1000.0)
        obs.create_dataset('date', data=np.array([d.strftime('%Y%m%d') for d in DATES],
                                                 dtype='S8'))
        f.create_dataset('HDFEOS/GRIDS/timeseries/quality/mask',
                         data=np.ones((LENGTH, WIDTH), bool))
        f.create_dataset('HDFEOS/GRIDS/timeseries/geometry/height',
                         data=np.zeros((LENGTH, WIDTH), np.float32))
        for key, val in {'UNIT': 'm', 'REF_Y': '12', 'REF_X': '62', 'Y_FIRST': '-1.999',
                         'X_FIRST': '120.999', 'Y_STEP': '-0.002', 'X_STEP': '0.002'}.items():
            f.attrs[key] = val
    ts = TimeSeries(path, verbose=False)
    assert ts.kind == 'MintPy' and ts.n_date == N_DATE and ts.to_mm == 1000.0
    assert list(ts.layers) == ['mask', 'hgt'] and ts.keep.all()
    ts.close()


def test_not_a_time_series(tmp_path):
    path = str(tmp_path / 'other.h5')
    with h5py.File(path, 'w') as f:
        f.create_dataset('velocity', data=np.zeros((3, 3)))
    with pytest.raises(ValueError, match='neither'):
        TimeSeries(path)


def test_circle_pixels(licsbas):
    ts = TimeSeries(licsbas, verbose=False)
    x_km, y_km = ts.pixel_size_km(30)
    rows, cols = circle_pixels_rc(ts, 30, 40, 2.0)
    assert rows.size == pytest.approx(np.pi * 4 / (x_km * y_km), rel=0.05)
    dist = np.hypot((rows - 30) * y_km, (cols - 40) * x_km)
    assert dist.max() <= 2.0 and (30, 40) in set(zip(rows, cols))
    # every pixel within the radius is there: none just outside left behind
    rr, cc = np.mgrid[0:LENGTH, 0:WIDTH]
    inside = np.hypot((rr - 30) * y_km, (cc - 40) * x_km) <= 2.0
    assert inside.sum() == rows.size
    # clipped at the edge of the grid
    rows, cols = circle_pixels_rc(ts, 0, 0, 2.0)
    assert rows.min() == 0 and cols.min() == 0 and rows.size < np.pi * 4 / (x_km * y_km) / 3
    ts.close()


##############################################################################
#                              the viewer                                    #
##############################################################################
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from viewinsar import cli  # noqa: E402


def make_viewer(args):
    from viewinsar.tsviewer import TimeSeriesViewer

    inps = cli.tseries_cmd_line_parse(args)
    inps.disp_fig = False                      # Agg, no plt.show()
    obj = TimeSeriesViewer(inps)
    obj.open()
    obj.plot()
    return obj


def test_tseries_options(licsbas, tmp_path, monkeypatch):
    inps = cli.tseries_cmd_line_parse([licsbas, '--noverbose'])
    assert inps.radius is None and inps.unit == 'mm' and inps.colormap == 'roma_r'
    assert inps.file2.endswith('cum.h5')                       # LiCSBAS: cum.h5 beside
    inps = cli.tseries_cmd_line_parse([licsbas, '--no-i2', '-r', '1.5', '-p', '40/30',
                                       '--ref', '0:10/0:10', '-m', '20180125'])
    assert inps.file2 is None and inps.radius == 1.5 and inps.point_xy == [40, 30]
    inps = cli.tseries_cmd_line_parse([licsbas, '--p_geo', '121.08/-2.06',
                                       '--ref_geo', '121.0/121.01/-2.0/-2.01'])
    assert inps.point_geo == [121.08, -2.06]
    # LiCSBAS's -r is the reference area: say so
    with pytest.raises(SystemExit):
        cli.tseries_cmd_line_parse([licsbas, '-r', '140:145/125:130'])
    for bad in (['-r', '0'], ['-p', '40'], ['--vlim', '5', '1'], ['-m', '2018-01-01'],
                ['--yx', '1', '2', '-p', '1/2'], ['--ref', '0:1/0:1', '--ref-yx', '1', '2'],
                ['--i2', 'nowhere.h5'], ['--auto-crange', '30']):
        with pytest.raises((ValueError, SystemExit)):
            cli.tseries_cmd_line_parse([licsbas, *bad])
    # no file named: cum_filt.h5 in the current folder
    monkeypatch.chdir(os.path.dirname(licsbas))
    assert cli.tseries_cmd_line_parse([]).file == 'cum_filt.h5'
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match='no time series given'):
        cli.tseries_cmd_line_parse([])


def test_radius_hint(licsbas, capsys):
    with pytest.raises(SystemExit):
        cli.tseries_cmd_line_parse([licsbas, '-r', '140:145/125:130'])
    assert '--ref X1:X2/Y1:Y2' in capsys.readouterr().err


def test_circle_end_to_end(licsbas, tmp_path):
    """--nodisplay: the median and quartiles of a circle, exactly."""
    out = tmp_path / 'out'
    obj = cli.tseries_main([licsbas, '--no-i2', '-r', '1', '-p', '40/30', '--nodisplay',
                            '--noverbose', '--outdir', str(out)])
    base = 'tseries_Y30X40_r1km'
    assert sorted(os.listdir(out)) == sorted(f'{base}{t}' for t in
                                             ('.csv', '.png', '_map.png', '_info.txt'))
    rows, cols = circle_pixels_rc(obj.ts, 30, 40, 1.0)
    want = expected_series(rows, cols)
    table = np.genfromtxt(out / f'{base}.csv', delimiter=',', names=True, dtype=None,
                          encoding='utf-8')
    assert list(table.dtype.names) == ['date', 'year', 'median_mm', 'q25_mm', 'q75_mm', 'n_px']
    assert table['median_mm'] == pytest.approx(np.median(want, axis=1), abs=1e-3)
    assert table['q25_mm'] == pytest.approx(np.percentile(want, 25, axis=1), abs=1e-3)
    assert table['q75_mm'] == pytest.approx(np.percentile(want, 75, axis=1), abs=1e-3)
    assert (table['n_px'] == rows.size).all()
    assert table['year'][0] == pytest.approx(2018.0)
    # the atmosphere common to all pixels is gone; the error bars are a real spread
    assert (table['q75_mm'] - table['q25_mm'])[-1] > 1.0
    info = (out / f'{base}_info.txt').read_text()
    assert f'{rows.size} pixels, {rows.size} used' in info
    assert 'reference area : 60:65/10:15' in info
    # the velocity of the circle: the median rate, relative to the reference area
    assert obj.result['vel'][0] == pytest.approx(np.median(rate_field()[rows, cols]) - 11.0,
                                                 abs=1e-4)


def test_one_pixel_and_the_reference_date(licsbas):
    obj = make_viewer([licsbas, '--no-i2', '--yx', '20', '50', '-m',
                       DATES[5].strftime('%Y%m%d'), '--noverbose'])
    assert obj.ix_ref == 5
    stats = obj.result['series'][0]
    want = expected_series(np.array([20]), np.array([50]), ix_ref=5)[:, 0]
    assert stats['median'] == pytest.approx(want, abs=1e-3)
    assert stats['q25'] == pytest.approx(stats['q75'])                 # one pixel: no bar
    assert stats['median'][5] == 0.0
    # no error bars drawn for one pixel, a scatter as LiCSBAS
    assert not obj.ax_ts.containers
    obj.close()


def test_interactive(licsbas, tmp_path):
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    obj = make_viewer([licsbas, '-r', '1', '--noverbose', '--outdir', str(tmp_path)])
    assert obj.ts2 is not None and len(obj.result['series']) == 2       # cum.h5 in red
    obj.fig_map.canvas.draw()
    canvas = obj.fig_map.canvas

    def click(col, row, dbl=True, button=MouseButton.LEFT):
        x, y = obj.ax_map.transData.transform((col, row))
        obj.on_click(MouseEvent('button_press_event', canvas, x, y, button=button,
                                dblclick=dbl))

    click(40, 30)
    assert obj.point == (30, 40)
    bars = obj.ax_ts.containers
    assert len(bars) == 2                                    # an error bar per file
    rows, cols = circle_pixels_rc(obj.ts, 30, 40, 1.0)
    want = np.median(expected_series(rows, cols), axis=1)
    assert bars[0].lines[0].get_ydata() == pytest.approx(want, abs=1e-3)
    assert 'median of' in obj.ax_ts.get_title()
    click(20, 20, dbl=False)                                 # a single click: nothing
    click(20, 20, button=MouseButton.RIGHT)
    assert obj.point == (30, 40)

    # the gap of the network: in the east half, so a circle at x 40 has half
    click(60, 30)
    assert any(c.get_label() == 'Gap of ifg network' for c in obj.ax_ts.collections)

    # right drag: a new reference area, and the series follow it
    class Ev:
        def __init__(self, x, y):
            self.xdata, self.ydata = x, y

    obj.on_ref_select(Ev(10, 5), Ev(20, 15))
    assert obj.ref == (10, 20, 5, 15)
    new_ref_rate = rate_field()[5:15, 10:20].mean()
    rows, cols = circle_pixels_rc(obj.ts, 30, 60, 1.0)
    want = np.median((rate_field()[rows, cols] - new_ref_rate)[None] * YEARS[:, None], axis=1)
    assert obj.result['series'][0]['median'] == pytest.approx(want, abs=1e-3)
    obj.on_ref_select(Ev(0, 0), Ev(3, 3))                    # no data there: refused
    assert obj.ref == (10, 20, 5, 15)

    # the slider: the cumulative displacement at a date, relative to the reference
    obj.slider.set_val(obj.date_nums[-1])
    assert obj.cum_mode and obj.date_index == N_DATE - 1
    shown = obj.image.get_array()
    assert shown[30, 60] == pytest.approx((rate_field()[30, 60] - new_ref_rate) * YEARS[-1],
                                          abs=1e-3)
    assert '(Ref: 2018/01/01)' in obj.ax_map.get_title()
    # the layers
    obj.radio.set_active(list(obj.layers).index('hgt'))
    assert not obj.cum_mode and obj.ax_map.get_title() == 'hgt'
    obj.radio.set_active(0)
    assert obj.ax_map.get_title() == 'vel(1)'

    # the mask (south masked): a circle across it uses fewer pixels with it on
    obj.select(39, 40)
    kept_on = obj.result['n_kept']
    obj.toggle_mask()
    assert not obj.use_mask and obj.result['n_kept'] > kept_on
    obj.toggle_mask()
    assert obj.result['n_kept'] == kept_on

    # the fit models
    visible = [ln.get_visible() for ln in obj.fit_lines[0]]
    assert visible == [True, True, False, False]
    obj.toggle_fit('Quad')
    assert obj.fit_lines[0][2].get_visible() and obj.fit_lines[1][2].get_visible()

    # export
    obj.on_key(KeyEvent('key_press_event', canvas, key='e'))
    names = os.listdir(tmp_path)
    assert 'tseries_Y39X40_r1km.csv' in names and 'tseries_Y39X40_r1km_map.png' in names
    header = (tmp_path / 'tseries_Y39X40_r1km.csv').read_text().splitlines()[0]
    assert header == ('date,year,median_mm_1,q25_mm_1,q75_mm_1,n_px_1,'
                      'median_mm_2,q25_mm_2,q75_mm_2,n_px_2')
    obj.close()


def test_a_circle_with_no_data(licsbas):
    obj = make_viewer([licsbas, '--no-i2', '-r', '0.3', '--yx', '30', '1', '--noverbose'])
    assert obj.result['n_kept'] == 0 and obj.result['series'] == []
    assert obj.ax_ts.get_title().startswith('NaN @(1, 30)')
    assert obj.save_outputs() == []
    obj.close()


def test_mintpy_viewer_in_cm_and_by_lalo(tmp_path):
    path = write_mintpy(str(tmp_path / 'mintpy'))
    lat, lon = LAT1 + 30 * DLAT, LON1 + 40 * DLON
    obj = make_viewer([path, '-r', '1', '--lalo', str(lat), str(lon), '-u', 'cm',
                       '--noverbose'])
    assert obj.point == (30, 40) and obj.ix_ref == 2                   # REF_DATE
    rows, cols = circle_pixels_rc(obj.ts, 30, 40, 1.0)
    rates = rate_field()[rows, cols] - rate_field()[12, 62]
    want = np.median(rates[None] * (YEARS - YEARS[2])[:, None], axis=1) / 10.0
    assert obj.result['series'][0]['median'] == pytest.approx(want, abs=1e-4)
    assert obj.ax_ts.get_ylabel() == 'Displacement (cm)'
    assert 'temporalCoherence' in obj.layers and 'hgt' in obj.layers
    obj.close()


def test_mintpy_utm_by_lalo(tmp_path):
    """A UTM grid (as your geo_timeseries): --lalo in degrees still finds the pixel."""
    pyproj = pytest.importorskip('pyproj')
    north0, east0, step = 9780000.0, 840000.0, 100.0          # UTM 51S, Sulawesi
    atr = {'FILE_TYPE': 'timeseries', 'UNIT': 'm', 'Y_FIRST': str(north0),
           'X_FIRST': str(east0), 'Y_STEP': str(-step), 'X_STEP': str(step),
           'Y_UNIT': 'meters', 'X_UNIT': 'meters', 'EPSG': '32751', 'UTM_ZONE': '51S',
           'REF_Y': '5', 'REF_X': '5'}
    path = str(tmp_path / 'geo_timeseries.h5')
    cube = rate_field()[None] * YEARS[:, None, None] / 1000.0
    with h5py.File(path, 'w') as f:
        f.create_dataset('timeseries', data=cube.astype(np.float32))
        f.create_dataset('date', data=np.array([d.strftime('%Y%m%d') for d in DATES],
                                               dtype='S8'))
        for key, val in atr.items():
            f.attrs[key] = val
    east, north = east0 + 25.5 * step, north0 - 20.5 * step    # the centre of pixel (20, 25)
    lon, lat = pyproj.Transformer.from_crs(32751, 4326, always_xy=True).transform(east, north)
    obj = make_viewer([path, '-r', '0.5', '--lalo', f'{lat:.7f}', f'{lon:.7f}', '--noverbose'])
    assert obj.point == (20, 25)
    assert obj.result['lalo'] == pytest.approx((north, east))
    assert obj.side_text.get_text().startswith(f'N:{north:.1f}\nE:{east:.1f}')
    # 100 m pixels: a 0.5 km circle holds about pi * 25 of them
    assert obj.result['n_total'] == pytest.approx(np.pi * 25, rel=0.1)
    obj.close()


def test_radar_coded_mintpy(tmp_path):
    """Descending: the map is flipped left-right as MintPy does; lat/lon via
    geometryRadar.h5."""
    from test_radar import _radar_attrs, lat_of, lon_of, make_radar_product

    vfile, _ = make_radar_product(str(tmp_path / 'proc'))
    with h5py.File(vfile, 'r') as f:
        vel = f['velocity'][()]
    years = YEARS[:12]
    path = os.path.join(os.path.dirname(vfile), 'timeseries.h5')
    with h5py.File(path, 'w') as f:
        f.create_dataset('timeseries', data=(vel[None] * years[:, None, None]).astype(np.float32))
        f.create_dataset('date', data=np.array([d.strftime('%Y%m%d') for d in DATES[:12]],
                                               dtype='S8'))
        for key, val in _radar_attrs('timeseries', 'm').items():
            f.attrs[key] = val
    obj = make_viewer([path, '-r', '1', '--lalo', f'{lat_of(100, 150):.5f}',
                       f'{lon_of(100, 150):.5f}', '--noverbose'])
    assert obj.point == pytest.approx((100, 150), abs=1)
    # columns run west on a descending track: east goes to the right; rows stay
    # top-down, as for any image (matplotlib calls that 'inverted' too)
    assert obj.ax_map.xaxis_inverted() and obj.ax_map.yaxis_inverted()
    assert obj.result['lalo'] == pytest.approx((lat_of(*obj.point), lon_of(*obj.point)),
                                               abs=1e-3)
    assert obj.ts.vel_source.endswith('velocity.h5')                 # found beside it
    obj.close()


def test_help_and_top_level(capsys):
    cli.main(['--help'])
    out = capsys.readouterr().out
    assert 'tseries' in out and 'LiCSBAS_plot_ts.py' in out
    with pytest.raises(SystemExit):
        cli.main(['tseries', '--help'])
    out = capsys.readouterr().out
    for needle in ('--radius KM', 'median', '25th', '--ref X1:X2/Y1:Y2', '--p-geo', '--i2',
                   'cum_filt.h5', 'timeseries', 'HDF-EOS5', 'right drag'):
        assert needle in out, needle


def test_module_is_light():
    """`import viewinsar` still pulls in neither matplotlib nor MintPy."""
    import subprocess
    import sys

    code = ('import sys, viewinsar.tsdata; '
            'print(sorted(m for m in ("matplotlib", "mintpy", "h5py") if m in sys.modules))')
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True,
                         check=True).stdout
    assert json.loads(out.replace("'", '"')) == []
    assert tsdata.MODELS == ('Linear', 'Annual+L', 'Quad', 'Annual+Q')

"""`viewinsar tseries` with several files: frames of a track on one map.

Two synthetic LiCSBAS frames of one track, with test_tseries' exact series:
the west one, and the east one 50 columns further east - they overlap on 30
columns, 25 of them with data - which starts two dates later. The ground moves
at the same rate in both, so once their reference areas are on the same
ground, their series at a point are the same.
"""

import os

import numpy as np
import pytest

from test_tseries import (
    DATES,
    DLAT,
    DLON,
    LAT1,
    LENGTH,
    LON1,
    N_DATE,
    WIDTH,
    YEARS,
    rate_field,
    write_licsbas,
    write_mintpy,
)

h5py = pytest.importorskip('h5py')
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from viewinsar import cli  # noqa: E402
from viewinsar.tsdata import circle_pixels_rc  # noqa: E402
from viewinsar.viewer import file_labels  # noqa: E402

SHIFT = 50                    # column 0 of the east frame is column 50 of the west one
CSV_HEADER = ('date,year,median_mm_1,q25_mm_1,q75_mm_1,n_px_1,'
              'median_mm_2,q25_mm_2,q75_mm_2,n_px_2')


@pytest.fixture
def frames(tmp_path):
    """The west frame (all dates) and the east one (from the third date)."""
    west = write_licsbas(str(tmp_path / '061D_09234_131313' / 'TS_GEOCml10GACOS'))
    east = write_licsbas(str(tmp_path / '061D_09012_161516' / 'TS_GEOCml10GACOS'),
                         shift=SHIFT, index=range(2, N_DATE))
    return west, east


def lalo(row, col):
    """lat, lon of a pixel centre of the west frame's grid (col may lie east of it)."""
    return LAT1 + row * DLAT, LON1 + col * DLON


def ground_rate(cols):
    """mm/yr of the ground at columns of the west grid (east of it too)."""
    return 0.5 * (np.asarray(cols, dtype=float) - 40.0)


def make_viewer(args):
    from viewinsar.tsmulti import MultiTimeSeriesViewer

    inps = cli.tseries_cmd_line_parse(args)
    inps.disp_fig = False                      # Agg, no plt.show()
    obj = MultiTimeSeriesViewer(inps)
    obj.open()
    obj.plot()
    return obj


def shown(frame):
    """What the map shows of a file, NaN where nothing."""
    return np.ma.filled(np.ma.asarray(frame.image.get_array(), dtype=float), np.nan)


class Ev:
    def __init__(self, lat, lon):
        self.xdata, self.ydata = lon, lat


##############################################################################
def test_file_labels():
    assert file_labels(['061D_09234_131313/TS_GEOCml10GACOS/cum.h5',
                        '061D_09012_161516/TS_GEOCml10GACOS/cum.h5']) == [
        '061D_09234_131313', '061D_09012_161516']
    assert file_labels(['x/a/b/vel.h5', 'y/a/b/vel.h5']) == ['x', 'y']
    assert file_labels(['a/vel.tif', 'b/vel.tif']) == ['a/vel', 'b/vel']        # as before
    assert file_labels(['TS/cum.h5', 'TS/cum_filt.h5']) == ['cum', 'cum_filt']


def test_options(frames, tmp_path):
    west, east = frames
    inps = cli.tseries_cmd_line_parse([west, east, '-r', '15'])
    assert inps.files == [west, east] and inps.file == west and inps.file2 is None
    assert inps.radius == 15
    # a quoted wildcard (no shell to expand it on Windows): sorted, each file once
    pattern = str(tmp_path / '061D_*' / 'TS_GEOCml10GACOS' / 'cum.h5')
    assert cli.tseries_cmd_line_parse([pattern, east]).files == [east, west]
    inps = cli.tseries_cmd_line_parse([west, east, '--labels', 'W', 'E',
                                       '--data-scale', '1', '0.5'])
    assert inps.labels == ['W', 'E'] and inps.data_scale == [1, 0.5]
    for bad, match in ((['--labels', 'W'], '--labels'), (['--i2', west], '--i2'),
                       (['--results', str(tmp_path)], '--results'),
                       (['--vel-file', west], '--vel-file'),
                       (['--data-scale', '1', '2', '3'], '--data-scale')):
        with pytest.raises(ValueError, match=match):
            cli.tseries_cmd_line_parse([west, east, *bad])
    with pytest.raises(ValueError, match='no such file'):
        cli.tseries_cmd_line_parse([west, str(tmp_path / 'nowhere.h5')])
    with pytest.raises(ValueError, match='no such file'):
        cli.tseries_cmd_line_parse([str(tmp_path / 'nowhere_*' / 'cum.h5')])


def test_two_frames_end_to_end(frames, tmp_path):
    """--nodisplay, a reference area on the same ground in both: the same
    series, each on its own dates."""
    west, east = frames
    lat, lon = lalo(30, 65)                      # in the overlap: east column 15
    (lat1, lon1), (lat2, lon2) = lalo(5, 60), lalo(15, 70)
    out = tmp_path / 'out'
    obj = cli.tseries_main([west, east, '--radius', '1', '--lalo', str(lat), str(lon),
                            '--ref-geo', f'{lon1}/{lon2}/{lat2}/{lat1}', '--nodisplay',
                            '--noverbose', '--outdir', str(out)])
    fw, fe = obj.frames
    assert fw.ref == (60, 71, 5, 16) and fe.ref == (10, 21, 5, 16)
    # the common reference date: the latest first date, the east frame's
    assert (fw.ix_ref, fe.ix_ref) == (2, 0)
    base = f'tseries_lat{lat:.4f}_lon{lon:.4f}_r1km'
    assert sorted(os.listdir(out)) == sorted(f'{base}{t}' for t in
                                             ('.csv', '.png', '_map.png', '_info.txt'))
    rows, cols = circle_pixels_rc(fw.ts, 30, 65, 1.0)
    want = np.median((ground_rate(cols) - ground_rate(65))[None]
                     * (YEARS - YEARS[2])[:, None], axis=1)
    assert (out / f'{base}.csv').read_text().splitlines()[0] == CSV_HEADER
    table = np.genfromtxt(out / f'{base}.csv', delimiter=',', names=True, dtype=None,
                          encoding='utf-8')
    assert len(table) == N_DATE                                  # every date of either
    assert table['median_mm_1'] == pytest.approx(want, abs=1e-3)
    assert table['median_mm_2'][2:] == pytest.approx(want[2:], abs=1e-3)
    assert np.isnan(table['median_mm_2'][:2]).all() and list(table['n_px_2'][:2]) == [0, 0]
    assert (table['n_px_1'] == rows.size).all() and (table['n_px_2'][2:] == rows.size).all()
    info = (out / f'{base}_info.txt').read_text()
    for needle in ('[1] 061D_09234_131313', '[2] 061D_09012_161516',
                   'reference area : 60:71/5:16', 'reference area : 10:21/5:16',
                   'reference date : 20180218', f'{rows.size} pixels, {rows.size} used',
                   'velocity       : 0.000 mm/yr (median of the circle)'):
        assert needle in info, needle


def test_two_frames_interactive(frames, tmp_path, capsys):
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    west, east = frames
    obj = make_viewer([west, east, '--radius', '1', '--outdir', str(tmp_path)])
    assert 'its own reference area' in capsys.readouterr().out
    fw, fe = obj.frames
    assert [fr.label for fr in obj.frames] == ['061D_09234_131313', '061D_09012_161516']
    assert obj.ref_day.strftime('%Y%m%d') == '20180218' and (fw.ix_ref, fe.ix_ref) == (2, 0)
    # each file drawn where it lies, with its outline: the pixels with data, x 5 - 79
    assert fw.image.get_extent() == pytest.approx(
        [LON1 - DLON / 2, LON1 + (WIDTH - 0.5) * DLON, LAT1 + (LENGTH - 0.5) * DLAT,
         LAT1 - DLAT / 2])
    assert fe.image.get_extent()[0] == pytest.approx(LON1 + (SHIFT - 0.5) * DLON)
    xs, ys = fe.outline.get_data()
    assert min(xs) == pytest.approx(LON1 + (SHIFT + 4.5) * DLON)
    assert max(xs) == pytest.approx(LON1 + (SHIFT + 79.5) * DLON)
    assert max(ys) == pytest.approx(LAT1 - DLAT / 2)
    assert np.isnan(shown(fw)[:, :5]).all() and np.isfinite(shown(fw)[:40, 5:]).all()

    obj.fig_map.canvas.draw()
    canvas = obj.fig_map.canvas

    def click(lat, lon, dbl=True, button=MouseButton.LEFT):
        x, y = obj.ax_map.transData.transform((lon, lat))
        obj.on_click(MouseEvent('button_press_event', canvas, x, y, button=button,
                                dblclick=dbl))

    # the overlap: both files, each relative to its own reference area
    click(*lalo(30, 65))
    here = obj.result['frames']
    assert [r['frame'] for r in here] == [fw, fe]
    assert [(r['row'], r['col']) for r in here] == [(30, 65), (30, 15)]
    rows, cols = circle_pixels_rc(fw.ts, 30, 65, 1.0)
    t = YEARS - YEARS[2]

    def want(ref_rate):
        return np.median((ground_rate(cols) - ref_rate)[None] * t[:, None], axis=1)

    assert here[0]['stats']['median'] == pytest.approx(want(11.0), abs=1e-3)   # x 60-64
    assert here[1]['stats']['median'] == pytest.approx(want(36.0)[2:], abs=1e-3)
    assert obj.ax_ts.get_title().startswith('vel(1) = 1.5, vel(2) = -23.5 mm/yr @(-2.06')
    assert len(obj.ax_ts.containers) == 2                        # an error bar each
    click(*lalo(30, 20), dbl=False)                              # a single click: nothing
    click(*lalo(30, 20), button=MouseButton.RIGHT)
    assert obj.result['frames'][0]['col'] == 65

    # a right drag in the overlap: both files on that ground - the same series
    obj.on_ref_select(Ev(*lalo(5, 60)), Ev(*lalo(15, 70)))
    assert fw.ref == (60, 71, 5, 16) and fe.ref == (10, 21, 5, 16)
    a, b = (r['stats']['median'] for r in obj.result['frames'])
    assert a == pytest.approx(want(12.5), abs=1e-3)
    assert b == pytest.approx(a[2:], abs=1e-4)
    # ... one where the east file has no data (its x 0-3): it keeps its own
    capsys.readouterr()
    obj.on_ref_select(Ev(*lalo(5, 50)), Ev(*lalo(15, 53)))
    assert fw.ref == (50, 54, 5, 16) and fe.ref == (10, 21, 5, 16)
    assert '2: 061D_09012_161516 keeps its reference area' in capsys.readouterr().out
    # ... and one with data in neither (the south is masked): refused
    obj.on_ref_select(Ev(*lalo(45, 50)), Ev(*lalo(55, 70)))
    assert fw.ref == (50, 54, 5, 16) and fe.ref == (10, 21, 5, 16)
    assert 'has all nan' in capsys.readouterr().out

    # a point in one file only, and one in neither
    click(*lalo(30, 20))
    assert [r['frame'] for r in obj.result['frames']] == [fw]
    assert 'not here: 2' in obj.side_text.get_text()
    click(*lalo(LENGTH + 0.2, 30))                               # the map's margin
    assert obj.result['frames'][0]['col'] == 20

    # the slider: the east file has nothing before its first date
    obj.slider.set_val(obj.date_nums[0])
    assert obj.cum_mode and np.isfinite(shown(fw)).any() and not np.isfinite(shown(fe)).any()
    obj.slider.set_val(obj.date_nums[-1])
    ref_rate = ground_rate(np.arange(50, 54)).mean()
    assert shown(fw)[30, 65] == pytest.approx((ground_rate(65) - ref_rate) * t[-1], abs=1e-3)
    assert '(Ref: 2018/02/18)' in obj.ax_map.get_title()
    # the layers: the velocity, then those of any file
    assert list(obj.layers) == ['vel', 'mask', 'coh_avg', 'n_unw', 'hgt']
    obj.radio.set_active(list(obj.layers).index('hgt'))
    assert not obj.cum_mode and obj.ax_map.get_title() == 'hgt'

    # the mask (the south of both): a circle across it uses more pixels without
    lat, lon = lalo(39, 65)
    obj.select(lon, lat)
    kept = [r['n_kept'] for r in obj.result['frames']]
    obj.toggle_mask()
    assert [r['n_kept'] > k for r, k in zip(obj.result['frames'], kept)] == [True, True]
    obj.toggle_mask()
    assert [r['n_kept'] for r in obj.result['frames']] == kept
    # the fit models, for every file
    obj.toggle_fit('Quad')
    assert len(obj.fit_lines) == 2 and all(lines[2].get_visible() for lines in obj.fit_lines)

    # export
    obj.on_key(KeyEvent('key_press_event', canvas, key='e'))
    base = f'tseries_lat{lat:.4f}_lon{lon:.4f}_r1km'
    assert f'{base}_map.png' in os.listdir(tmp_path)
    assert (tmp_path / f'{base}.csv').read_text().splitlines()[0] == CSV_HEADER
    obj.close()


def test_versions_of_one_time_series(tmp_path):
    """MintPy's several files: versions of one series, on one grid - the pixel
    options and --mask apply to every file."""
    one, two = write_mintpy(str(tmp_path / 'a')), write_mintpy(str(tmp_path / 'b'))
    mask = os.path.join(str(tmp_path / 'a'), 'maskTempCoh.h5')
    obj = make_viewer([one, two, '--radius', '1', '--yx', '30', '40', '--ref', '60:65/10:15',
                       '--mask', mask, '--noverbose'])
    assert obj.one_grid and [fr.ref for fr in obj.frames] == [(60, 65, 10, 15)] * 2
    assert [fr.label for fr in obj.frames] == ['a/timeseries', 'b/timeseries']
    assert all(fr.ts.mask_source == mask for fr in obj.frames)
    rows, cols = circle_pixels_rc(obj.frames[0].ts, 30, 40, 1.0)
    want = np.median((rate_field()[rows, cols] - 11.0)[None] * (YEARS - YEARS[2])[:, None],
                     axis=1)                                     # MintPy: REF_DATE
    a, b = (r['stats']['median'] for r in obj.result['frames'])
    assert a == pytest.approx(want, abs=1e-3) and b == pytest.approx(a)
    obj.close()


def test_radar_coded_on_one_grid(tmp_path):
    """Radar-coded files of one stack: the map in its pixels, flipped as MintPy
    does; lat/lon through geometryRadar.h5."""
    from test_radar import _radar_attrs, lat_of, lon_of, make_radar_product

    vfile, _ = make_radar_product(str(tmp_path / 'proc'))
    with h5py.File(vfile, 'r') as f:
        vel = f['velocity'][()]
    paths = []
    for name, factor in (('timeseries.h5', 1.0), ('timeseries_ERA5.h5', 0.5)):
        path = os.path.join(os.path.dirname(vfile), name)
        with h5py.File(path, 'w') as f:
            f.create_dataset('timeseries', data=(factor * vel[None]
                                                 * YEARS[:12, None, None]).astype(np.float32))
            f.create_dataset('date', data=np.array([d.strftime('%Y%m%d') for d in DATES[:12]],
                                                   dtype='S8'))
            for key, val in _radar_attrs('timeseries', 'm').items():
                f.attrs[key] = val
        paths.append(path)
    obj = make_viewer([*paths, '--radius', '1', '--lalo', f'{lat_of(60, 220):.5f}',
                       f'{lon_of(60, 220):.5f}', '--noverbose'])
    assert not obj.geocoded and obj.one_grid
    assert obj.ax_map.xaxis_inverted() and obj.ax_map.yaxis_inverted()   # descending
    here = obj.result['frames']
    assert len(here) == 2 and (here[0]['row'], here[0]['col']) == pytest.approx((60, 220),
                                                                                abs=1)
    assert obj.result['lalo'] == pytest.approx((lat_of(here[0]['row'], here[0]['col']),
                                                lon_of(here[0]['row'], here[0]['col'])),
                                               abs=1e-3)
    a, b = (r['stats']['median'] for r in here)
    assert np.nanmax(np.abs(a)) > 1 and b == pytest.approx(0.5 * a, abs=1e-6)
    assert obj.outfile_base.endswith(f'tseries_Y{here[0]["row"]}X{here[0]["col"]}_r1km')
    obj.close()


def test_files_that_cannot_share_a_map(frames, tmp_path):
    from viewinsar.tsmulti import MultiTimeSeriesViewer

    west, east = frames

    def opened(args):
        inps = cli.tseries_cmd_line_parse([*args, '--noverbose'])
        inps.disp_fig = False
        obj = MultiTimeSeriesViewer(inps)
        try:
            obj.open()
        finally:
            obj.close()

    # what is in the pixels of one grid
    mask = os.path.join(os.path.dirname(west), 'results', 'mask')
    for opt in (['--yx', '30', '40'], ['--point', '40/30'], ['--ref', '0:10/0:10'],
                ['--ref-yx', '3', '3'], ['--mask', mask]):
        with pytest.raises(ValueError, match='different grids'):
            opened([west, east, *opt])
    # a point outside every file; a reference area without data in any
    with pytest.raises(ValueError, match='outside every file'):
        opened([west, east, '--lalo', '-3.0', '121.0'])
    (lat1, lon1), (lat2, lon2) = lalo(50, 5), lalo(55, 10)          # the masked south
    with pytest.raises(ValueError, match='no file has a valid pixel'):
        opened([west, east, '--ref-geo', f'{lon1}/{lon2}/{lat2}/{lat1}'])

    # radar-coded with geocoded; lat/lon with UTM
    radar = str(tmp_path / 'timeseries.h5')
    with h5py.File(radar, 'w') as f:
        f.create_dataset('timeseries', data=np.zeros((3, 10, 12), np.float32))
        f.create_dataset('date', data=np.array([d.strftime('%Y%m%d') for d in DATES[:3]],
                                               dtype='S8'))
        f.attrs['UNIT'] = 'm'
    with pytest.raises(ValueError, match='all be geocoded'):
        opened([west, radar])
    utm = str(tmp_path / 'geo_timeseries.h5')
    with h5py.File(utm, 'w') as f:
        f.create_dataset('timeseries', data=np.zeros((3, 10, 12), np.float32))
        f.create_dataset('date', data=np.array([d.strftime('%Y%m%d') for d in DATES[:3]],
                                               dtype='S8'))
        for key, val in {'UNIT': 'm', 'Y_FIRST': '9780000', 'X_FIRST': '840000',
                         'Y_STEP': '-100', 'X_STEP': '100', 'Y_UNIT': 'meters',
                         'X_UNIT': 'meters', 'EPSG': '32751'}.items():
            f.attrs[key] = val
    with pytest.raises(ValueError, match='one CRS'):
        opened([west, utm])


def test_help_tells_of_several_files(capsys):
    with pytest.raises(SystemExit):
        cli.main(['tseries', '--help'])
    out = capsys.readouterr().out
    for needle in ('--labels', 'several files', 'one map', 'TS_GEOCml10GACOS/cum.h5'):
        assert needle in out, needle


def test_one_file_with_a_label(frames):
    """--labels with one file: its name in the legend; the viewer is the one
    of a single file, as before."""
    from test_tseries import make_viewer as make_single
    from viewinsar.tsviewer import TimeSeriesViewer

    obj = make_single([frames[0], '--labels', 'west', '--noverbose'])
    assert isinstance(obj, TimeSeriesViewer) and obj.ts.label == '1: west'
    assert obj.ax_ts.get_legend().get_texts()[0].get_text() == '1: west'
    obj.close()

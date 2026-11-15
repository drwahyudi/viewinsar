"""GNSS on the swath profiles, and several files in one swath.

The readers (text with any header, GMT psvelo, .xlsx, GeoJSON, KML, shapefile),
the components - the LOS against MintPy's own; then
`viewinsar swath -gnss` on the synthetic velocity (an arctan across the N-S
fault at 130.45 E) and `viewinsar swath FILE1 FILE2` on two frames of it.
"""

import json
import zipfile

import numpy as np
import pytest

from conftest import FAULT_LON, KM_PER_DEG_LON, _base_attrs, make_velocity
from viewinsar import gnss as gn

h5py = pytest.importorskip('h5py')


##############################################################################
#                                 reading                                    #
##############################################################################
def test_read_csv_with_units_in_the_header(tmp_path):
    path = tmp_path / 'a.csv'
    path.write_text('Station,Longitude,Latitude,Ve(cm/yr),Vn(cm/yr),Vu(cm/yr),sig_e,sig_n,sig_u\n'
                    'AAAA,120.1,-2.0,1.05,-0.32,0.11,0.05,0.06,0.15\n'
                    'BBBB,120.3,-2.1,1.20,-0.20,,0.04,0.05,\n')
    g = gn.read_gnss(str(path), verbose=False)
    assert g.unit == 'cm/yr' and len(g) == 2
    assert list(g.site) == ['AAAA', 'BBBB']
    assert g.ve == pytest.approx([1.05, 1.20]) and g.sn == pytest.approx([0.06, 0.05])
    assert list(g.has_vu) == [True, False]
    mm = g.scaled(10.0, 'mm/yr')
    assert mm.ve == pytest.approx([10.5, 12.0]) and mm.se == pytest.approx([0.5, 0.4])
    assert mm.unit == 'mm/yr'


def test_read_headerless_psvelo_and_commented_header(tmp_path):
    vel = tmp_path / 'b.vel'
    vel.write_text('120.1 -2.0 10.5 -3.2 0.5 0.6 0.01 AAAA\n120.3 -2.1 12.0 -2.0 0.4 0.5 0.0 BBBB\n')
    g = gn.read_gnss(str(vel), verbose=False)
    assert list(g.site) == ['AAAA', 'BBBB'] and g.vn == pytest.approx([-3.2, -2.0])
    assert not g.has_vu.any() and g.unit == 'mm/yr'

    txt = tmp_path / 'c.txt'
    txt.write_text('# lon lat ve vn vu se sn su site\n120.1 -2.0 10.5 -3.2 1 0.5 0.6 1 X1\n')
    g = gn.read_gnss(str(txt), verbose=False)
    assert g.vu == pytest.approx([1.0]) and list(g.site) == ['X1']


def _write_xlsx(path, rows):
    """A minimal workbook - shared strings for the text, numbers inline."""
    strings = sorted({c for r in rows for c in r if isinstance(c, str)})
    main = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    rel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    cells = []
    for i, row in enumerate(rows, start=1):
        out = []
        for j, val in enumerate(row):
            ref = f'{chr(65 + j)}{i}'
            if isinstance(val, str):
                out.append(f'<c r="{ref}" t="s"><v>{strings.index(val)}</v></c>')
            else:
                out.append(f'<c r="{ref}"><v>{val}</v></c>')
        cells.append(f'<row r="{i}">{"".join(out)}</row>')
    with zipfile.ZipFile(path, 'w') as z:
        z.writestr('xl/workbook.xml', f'<workbook xmlns="{main}" xmlns:r="{rel}"><sheets>'
                   '<sheet name="gnss" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels',
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
                   'relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml" '
                   'Type="x"/></Relationships>')
        z.writestr('xl/sharedStrings.xml', f'<sst xmlns="{main}">'
                   + ''.join(f'<si><t>{s}</t></si>' for s in strings) + '</sst>')
        z.writestr('xl/worksheets/sheet1.xml',
                   f'<worksheet xmlns="{main}"><sheetData>{"".join(cells)}</sheetData>'
                   '</worksheet>')


def test_read_xlsx(tmp_path):
    path = tmp_path / 'gnss.xlsx'
    _write_xlsx(path, [['site', 'lon', 'lat', 'Ve [mm/yr]', 'Vn [mm/yr]', 'Se', 'Sn'],
                       ['P1', 130.1, 30.2, 3.5, -1.0, 0.4, 0.5],
                       ['P2', 130.3, 30.25, 4.5, -1.5, 0.3, 0.4]])
    g = gn.read_gnss(str(path), verbose=False)
    assert list(g.site) == ['P1', 'P2'] and g.lon == pytest.approx([130.1, 130.3])
    assert g.ve == pytest.approx([3.5, 4.5]) and g.se == pytest.approx([0.4, 0.3])


def test_read_geojson_and_kml(tmp_path):
    feats = [{'type': 'Feature', 'geometry': {'type': 'Point', 'coordinates': [130.1, 30.2]},
              'properties': {'name': 'P1', 've': 3.5, 'vn': -1.0, 'vu': 0.2}},
             {'type': 'Feature', 'geometry': {'type': 'Point', 'coordinates': [130.3, 30.25]},
              'properties': {'name': 'P2', 've': 4.5, 'vn': -1.5}}]
    path = tmp_path / 'gnss.geojson'
    path.write_text(json.dumps({'type': 'FeatureCollection', 'features': feats}))
    g = gn.read_gnss(str(path), verbose=False)
    assert list(g.site) == ['P1', 'P2'] and g.vn == pytest.approx([-1.0, -1.5])
    assert list(g.has_vu) == [True, False]

    kml = tmp_path / 'gnss.kml'
    marks = ''.join(
        f'<Placemark><name>{f["properties"]["name"]}</name><ExtendedData>'
        f'<Data name="ve"><value>{f["properties"]["ve"]}</value></Data>'
        f'<Data name="vn"><value>{f["properties"]["vn"]}</value></Data></ExtendedData>'
        f'<Point><coordinates>{f["geometry"]["coordinates"][0]},'
        f'{f["geometry"]["coordinates"][1]},0</coordinates></Point></Placemark>'
        for f in feats)
    kml.write_text(f'<kml xmlns="http://www.opengis.net/kml/2.2"><Document>{marks}</Document>'
                   '</kml>')
    g = gn.read_gnss(str(kml), verbose=False)
    assert list(g.site) == ['P1', 'P2'] and g.ve == pytest.approx([3.5, 4.5])


def test_read_shapefile(tmp_path):
    shapefile = pytest.importorskip('shapefile')
    path = str(tmp_path / 'gnss.shp')
    with shapefile.Writer(path, shapeType=shapefile.POINT) as w:
        for name in ('site', 've', 'vn', 'se', 'sn'):
            w.field(name, 'C' if name == 'site' else 'N', decimal=0 if name == 'site' else 3)
        w.point(130.1, 30.2)
        w.record('P1', 3.5, -1.0, 0.4, 0.5)
        w.point(130.3, 30.25)
        w.record('P2', 4.5, -1.5, 0.3, 0.4)
    g = gn.read_gnss(path, verbose=False)
    assert list(g.site) == ['P1', 'P2'] and g.ve == pytest.approx([3.5, 4.5])
    assert g.lat == pytest.approx([30.2, 30.25])


def test_unreadable_tables_say_why(tmp_path):
    path = tmp_path / 'bad.csv'
    path.write_text('name,x_coord,y_coord,speed\nA,1,2,3\n')
    with pytest.raises(ValueError, match='lon'):
        gn.read_gnss(str(path), verbose=False)


##############################################################################
#                            the components                                  #
##############################################################################
def test_los_matches_mintpy():
    ut = pytest.importorskip('mintpy.utils.utils')
    inc, az = 38.5, gn.heading_to_azimuth(-12.0)
    assert az == pytest.approx(102.0)
    assert gn.heading_to_azimuth(-168.0) == pytest.approx(-102.0)
    g = gn.GNSS(['A', 'B'], [130.0, 130.1], [30.0, 30.1], [10.0, -4.0], [3.0, 2.0],
                vu=[1.0, np.nan], se=[1, 1], sn=[1, 1], su=[2, np.nan])
    value, sigma = g.component('los', gn.los_unit_vector(inc, az))
    want = ut.enu2los(np.array([10.0, -4.0]), np.array([3.0, 2.0]), np.array([1.0, 0.0]),
                      inc_angle=inc, az_angle=az)
    assert value == pytest.approx(want)
    ue, un, uu = gn.los_unit_vector(inc, az)
    assert sigma[1] == pytest.approx(np.hypot(ue, un))           # 2D: no vertical sigma
    assert sigma[0] == pytest.approx(np.sqrt(ue ** 2 + un ** 2 + (2 * uu) ** 2))
    assert g.component('vn')[0] == pytest.approx([3.0, 2.0])
    with pytest.raises(ValueError, match='unit vector'):
        g.component('los')


def test_los_geometry_from_a_file_or_the_attributes(tmp_path):
    atr = _base_attrs()
    atr.update(HEADING='-12.0', CENTER_INCIDENCE_ANGLE='39.0')
    geom = gn.LOSGeometry.for_file(str(tmp_path / 'velocity.h5'), atr, verbose=False)
    assert geom.constant and geom.az == pytest.approx(102.0)
    assert geom.at_lonlat([130.2], [30.3])[0] == pytest.approx([39.0])

    # a geocoded geometry file beside the data is found and used, pixel by pixel
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    inc = np.tile(np.linspace(30, 45, width), (length, 1))
    with h5py.File(tmp_path / 'geometryGeo.h5', 'w') as fid:
        fid.create_dataset('incidenceAngle', data=inc.astype(np.float32))
        fid.create_dataset('azimuthAngle', data=np.full(inc.shape, 102.0, np.float32))
        for key, val in dict(atr, FILE_TYPE='geometry').items():
            fid.attrs[key] = val
    geom = gn.LOSGeometry.for_file(str(tmp_path / 'velocity.h5'), atr, verbose=False)
    assert not geom.constant
    lon0, dlon = float(atr['X_FIRST']), float(atr['X_STEP'])
    col = 200
    got, _ = geom.at_lonlat([lon0 + (col + 0.5) * dlon], [30.3])
    assert got[0] == pytest.approx(inc[0, col], abs=1e-4)
    off, _ = geom.at_lonlat([lon0 - 1.0], [30.3])
    assert np.isnan(off[0])


##############################################################################
#                         swath: GNSS, several files                         #
##############################################################################
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from viewinsar import cli  # noqa: E402

ALONG_30_3 = ['--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75', '--width', '8']


def _east_mm(lon):
    """The fixture's noise-free velocity, in mm/yr, read as an east velocity."""
    return 300.0 * (2 / np.pi) * np.arctan((lon - FAULT_LON) * KM_PER_DEG_LON / 8.0)


@pytest.fixture(scope='module')
def gnss_csv(tmp_path_factory):
    path = tmp_path_factory.mktemp('gnss') / 'gnss.csv'
    rows = ['site,lon,lat,ve,vn,vu,se,sn,su']
    for i, lon in enumerate(np.linspace(130.20, 130.70, 6)):
        vu, su = ('0.5', '2.0') if i % 2 else ('', '')
        rows.append(f'G{i},{lon:.3f},{30.30 + 0.01 * (i % 2):.3f},{_east_mm(lon):.2f},1.0,'
                    f'{vu},3.0,3.0,{su}')
    rows.append('FAR,130.45,30.50,0.0,0.0,,3.0,3.0,')            # 22 km off the line
    path.write_text('\n'.join(rows) + '\n')
    return str(path)


def test_swath_with_gnss(velocity_file, gnss_csv, tmp_path):
    out = tmp_path / 'out'
    obj = cli.swath_main([velocity_file, *ALONG_30_3, '-gnss', gnss_csv, '-vel', 've',
                          '--gnss-radius', '2', '--gnss-label', '--nodisplay', '--noverbose',
                          '--outdir', str(out)])
    sec = obj.sections[0]
    gin = sec.gnss_in
    assert gin['index'].size == 6 and 'FAR' not in gin['site']
    # the data are in cm/year here: the GNSS were converted, and agree with the InSAR
    assert obj.disp_unit.startswith('cm')
    assert gin['value'] == pytest.approx(
        [_east_mm(lon) / 10 for lon in obj.gnss.lon[gin['index']]], abs=1e-3)
    # (the pixels carry 5 cm/yr of noise each: the shift is that of a few medians)
    assert gin['n_ref'] == 6 and abs(gin['offset']) < 0.7
    # triangles for the 3D stations, circles for the 2D ones - and no line through them
    labels = [t.get_text() for t in sec.ax_swath.get_legend().get_texts()]
    assert {'GNSS (3D)', 'GNSS (2D)'} <= set(labels)
    assert not any('GNSS' in lab and 'smooth' in lab for lab in labels)
    # all seven on the map (labelled), as triangles and circles
    markers = [c for c in obj.ax_img.collections if c.get_offsets().shape[0] in (3, 4)]
    assert sum(c.get_offsets().shape[0] for c in markers) == 7
    table = np.genfromtxt(next(out.glob('*_gnss.csv')), delimiter=',', names=True,
                          dtype=None, encoding='utf-8')
    assert len(table) == 6
    assert table['value_shifted'] == pytest.approx(table['value'] + gin['offset'], abs=1e-4)
    assert 'GNSS in swath  : 6 station(s)' in next(out.glob('*_info.txt')).read_text()


def test_swath_gnss_in_the_los(velocity_file, gnss_csv, tmp_path):
    """-vel los: through the geometry of the data (here its attributes)."""
    import shutil

    vfile = str(tmp_path / 'velocity.h5')
    shutil.copy(velocity_file, vfile)
    with h5py.File(vfile, 'a') as fid:
        fid.attrs['HEADING'] = '-12.0'
        fid.attrs['CENTER_INCIDENCE_ANGLE'] = '39.0'
    obj = cli.swath_main([vfile, *ALONG_30_3, '-gnss', gnss_csv, '--gnss-ref', 'none',
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'o')])
    gin = obj.sections[0].gnss_in
    ue, un, uu = gn.los_unit_vector(39.0, 102.0)
    idx = gin['index']
    want = (obj.gnss.ve[idx] * ue + obj.gnss.vn[idx] * un
            + np.nan_to_num(obj.gnss.vu[idx]) * uu)
    assert gin['value'] == pytest.approx(want, abs=1e-6)
    assert gin['offset'] == 0.0 and 'attributes' in obj.gnss_note


@pytest.fixture(scope='module')
def two_frames(tmp_path_factory):
    """The fixture velocity cut into two overlapping frames; the east one 1 cm/yr up."""
    folder = tmp_path_factory.mktemp('frames')
    vel = make_velocity()
    files = []
    for name, (c0, c1), shift in (('west', (0, 260), 0.0), ('east', (200, 450), 0.01)):
        atr = _base_attrs()
        atr.update(WIDTH=str(c1 - c0),
                   X_FIRST=str(float(atr['X_FIRST']) + c0 * float(atr['X_STEP'])))
        path = folder / name / 'velocity.h5'
        path.parent.mkdir()
        with h5py.File(path, 'w') as fid:
            fid.create_dataset('velocity', data=vel[:, c0:c1] + shift)
            for key, val in atr.items():
                fid.attrs[key] = val
        files.append(str(path))
    return files


def test_swath_of_two_frames(two_frames, velocity_file, tmp_path):
    from matplotlib.image import AxesImage

    one = cli.swath_main([velocity_file, *ALONG_30_3, '--nodisplay', '--noverbose',
                          '--outdir', str(tmp_path / 'one')])
    obj = cli.swath_main([*two_frames, '--labels', 'West', 'East', *ALONG_30_3,
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'two')])
    assert obj.labels == ['West', 'East']
    sec = obj.sections[0]
    assert len(sec.swath_list) == 2
    # the second frame is sampled along the same line, on its own grid
    east = sec.swath_list[1]
    assert east['lon'].min() >= float(_base_attrs()['X_FIRST']) + 200 * 0.002 - 1e-9
    assert east['dist'].max() == pytest.approx(sec.swath_list[0]['length'], abs=0.6)
    # one density for both: every sample of both frames
    n_all = sum(s['value'].size for s in sec.swath_list)
    assert sec.den['count'].sum() == n_all
    assert n_all > sum(s['value'].size for s in one.sections[0].swath_list)
    # a median line per frame, named, and both frames on the map
    labels = [t.get_text() for t in sec.ax_swath.get_legend().get_texts()]
    assert 'Median - West' in labels and 'Median - East' in labels
    assert len([a for a in obj.ax_img.get_images() if isinstance(a, AxesImage)]) == 2
    pts = np.genfromtxt(next((tmp_path / 'two').glob('*_points.csv')), delimiter=',',
                        names=True)
    assert set(pts['file_idx']) == {0, 1}
    info = next((tmp_path / 'two').glob('*_info.txt')).read_text()
    assert "['West', 'East']" in info


def test_swath_labels_must_match_the_files(two_frames):
    with pytest.raises(ValueError, match='--labels'):
        cli.cmd_line_parse([*two_frames, '--labels', 'only-one', '--noverbose'])


def test_swath_with_several_gnss_components(velocity_file, gnss_csv, tmp_path):
    """--vel ve vn vu: every component on the profile, each in its colour and
    shifted onto the InSAR on its own; a row per station and component."""
    out = tmp_path / 'out'
    obj = cli.swath_main([velocity_file, *ALONG_30_3, '--gnss', gnss_csv, '--vel', 've', 'vn',
                          'vu', 've', '--gnss-radius', '2', '--nodisplay', '--noverbose',
                          '--outdir', str(out)])
    assert obj.gnss_comps == ['ve', 'vn', 'vu'] and obj.gnss_comp == 've'   # each once
    sec = obj.sections[0]
    gin = sec.gnss_in
    idx = gin['index']
    assert list(gin['series']) == ['ve', 'vn', 'vu'] and idx.size == 6
    assert gin['series']['ve']['value'] == pytest.approx(
        [_east_mm(lon) / 10 for lon in obj.gnss.lon[idx]], abs=1e-3)
    assert gin['series']['vn']['value'] == pytest.approx(np.full(6, 0.1), abs=1e-6)  # cm/yr
    vu = gin['series']['vu']['value']
    assert np.isfinite(vu).sum() == 3                                 # the 3D stations
    # the first component is still where it always was
    assert gin['value'] == pytest.approx(gin['series']['ve']['value'])
    assert gin['offset'] == gin['series']['ve']['offset']
    offsets = {c: one['offset'] for c, one in gin['series'].items()}
    assert offsets['ve'] != offsets['vn'] and all(one['n_ref'] for one in gin['series'].values())
    labels = [t.get_text() for t in sec.ax_swath.get_legend().get_texts()]
    assert {'GNSS ve', 'GNSS vn', 'GNSS vu', '3D station (2D: circle)'} <= set(labels)
    table = np.genfromtxt(next(out.glob('*_gnss.csv')), delimiter=',', names=True,
                          dtype=None, encoding='utf-8')
    assert len(table) == 6 + 6 + 3 and set(table['component']) == {'ve', 'vn', 'vu'}
    info = next(out.glob('*_info.txt')).read_text()
    assert '--vel ve vn vu' in info and 'each component shifted' in info


def test_gnss_vu_without_3d_stations(velocity_file, tmp_path, capsys):
    path = tmp_path / 'g2d.csv'
    path.write_text('site,lon,lat,ve,vn\nA,130.30,30.30,1,2\nB,130.60,30.30,3,4\n')
    args = [velocity_file, *ALONG_30_3, '--gnss', str(path), '--nodisplay', '--noverbose',
            '--outdir', str(tmp_path / 'o')]
    obj = cli.swath_main(args + ['--vel', 'vu', 've'])                # left out, said
    assert obj.gnss_comps == ['ve'] and 'no station has a vertical' in capsys.readouterr().out
    obj = cli.swath_main(args + ['--vel', 'vu'])                       # alone: no GNSS
    assert obj.gnss is None and 'no station has a vertical' in capsys.readouterr().out


def test_several_components_survive_a_missing_los(velocity_file, gnss_csv, tmp_path, capsys):
    """No look angles anywhere (not even ORBIT_DIRECTION for nominal ones):
    --vel los ve keeps ve, and says why los is gone."""
    import shutil

    vfile = str(tmp_path / 'velocity.h5')
    shutil.copy(velocity_file, vfile)
    with h5py.File(vfile, 'a') as fid:
        del fid.attrs['ORBIT_DIRECTION']
    obj = cli.swath_main([vfile, *ALONG_30_3, '--gnss', gnss_csv, '--vel', 'los', 've',
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'o')])
    assert obj.gnss_comps == ['ve'] and 'los left out' in capsys.readouterr().out

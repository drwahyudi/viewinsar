"""--fault: vector readers, re-projection, crossings, and what gets drawn."""

import json
import sqlite3
import struct
import zipfile

import numpy as np
import pytest

from conftest import DLAT, DLON, FAULT_LON, LAT0, LENGTH, LON0, WIDTH, make_velocity
from viewinsar import cli, fault
from viewinsar.swath import line_geometry, swath_extract

pyproj = pytest.importorskip('pyproj')

LINE = ['--start-lalo', '30.30', '130.15', '--end-lalo', '30.30', '130.75']
UTM = 32652                                  # UTM 52N covers the fixtures (~130.4 E)


def _trace_lonlat():
    """The fixtures' fault: N-S at FAULT_LON, running past the data (30.0-30.6 N)."""
    lat = np.linspace(29.9, 30.7, 17)
    return np.c_[np.full_like(lat, FAULT_LON), lat]


def _to_utm(xy):
    tfm = pyproj.Transformer.from_crs(4326, UTM, always_xy=True)
    return np.c_[tfm.transform(xy[:, 0], xy[:, 1])]


def _esri_wkt(epsg):
    return pyproj.CRS.from_epsg(epsg).to_wkt(version='WKT1_ESRI')


##############################################################################
#                       one fault, every format                              #
##############################################################################
@pytest.fixture
def geojson_file(tmp_path):
    obj = {'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': {'name': 'main'},
         'geometry': {'type': 'LineString', 'coordinates': _trace_lonlat().tolist()}},
        {'type': 'Feature', 'properties': {'name': 'a point, not a trace'},
         'geometry': {'type': 'Point', 'coordinates': [130.3, 30.3]}}]}
    path = tmp_path / 'faults.geojson'
    path.write_text(json.dumps(obj))
    return str(path)


@pytest.fixture
def geojson_legacy_utm_file(tmp_path):
    """Pre-RFC-7946 GeoJSON: a named CRS, here UTM."""
    obj = {'type': 'FeatureCollection',
           'crs': {'type': 'name', 'properties': {'name': f'urn:ogc:def:crs:EPSG::{UTM}'}},
           'features': [{'type': 'Feature', 'properties': {},
                         'geometry': {'type': 'MultiLineString',
                                      'coordinates': [_to_utm(_trace_lonlat()).tolist()]}}]}
    path = tmp_path / 'faults_utm.geojson'
    path.write_text(json.dumps(obj))
    return str(path)


def _write_shp(path, lines, prj=None):
    shapefile = pytest.importorskip('shapefile')
    with shapefile.Writer(str(path), shapeType=shapefile.POLYLINE) as wrt:
        wrt.field('name', 'C')
        for i, line in enumerate(lines):
            wrt.line([line.tolist()])
            wrt.record(f'f{i}')
    if prj:
        path.with_suffix('.prj').write_text(prj)
    return str(path.with_suffix('.shp'))


@pytest.fixture
def shp_utm_file(tmp_path):
    """A UTM shapefile, with the ESRI .prj a GIS writes."""
    return _write_shp(tmp_path / 'faults_utm', [_to_utm(_trace_lonlat())], _esri_wkt(UTM))


@pytest.fixture
def shp_noprj_file(tmp_path):
    """lon/lat numbers, but no .prj."""
    return _write_shp(tmp_path / 'faults_noprj', [_trace_lonlat()])


@pytest.fixture
def gmt_file(tmp_path):
    xy = _trace_lonlat()
    text = '# fault traces, GMT multi-segment\n> -Z0\n'
    text += '\n'.join(f'{x} {y}' for x, y in xy[:9]) + '\n>\n'
    text += '\n'.join(f'{x} {y}' for x, y in xy[8:]) + '\n'
    path = tmp_path / 'faults.gmt'
    path.write_text(text)
    return str(path)


@pytest.fixture
def gmt_utm_file(tmp_path):
    """OGR-GMT output: the CRS lives in a '# @Je' header."""
    xy = _to_utm(_trace_lonlat())
    text = f'# @VGMT1.0 @GLINESTRING\n# @Je{UTM}\n# FEATURE_DATA\n>\n'
    text += '\n'.join(f'{x:.3f} {y:.3f}' for x, y in xy) + '\n'
    path = tmp_path / 'faults_utm.gmt'
    path.write_text(text)
    return str(path)


def _kml_text(lines):
    marks = ''.join(
        '<Placemark><LineString><coordinates>'
        + ' '.join(f'{x},{y},0' for x, y in line)
        + '</coordinates></LineString></Placemark>' for line in lines)
    return ('<?xml version="1.0" encoding="UTF-8"?>'
            '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>' + marks
            + '<Placemark><Point><coordinates>130.3,30.3,0</coordinates></Point></Placemark>'
            '</Document></kml>')


@pytest.fixture
def kml_file(tmp_path):
    path = tmp_path / 'faults.kml'
    path.write_text(_kml_text([_trace_lonlat()]))
    return str(path)


@pytest.fixture
def kmz_file(tmp_path):
    path = tmp_path / 'faults.kmz'
    with zipfile.ZipFile(path, 'w') as zfile:
        zfile.writestr('doc.kml', _kml_text([_trace_lonlat()]))
    return str(path)


def _wkb_linestring(line, big_endian=False, z=False):
    """WKB of a LineString - as shapely would write it, without shapely."""
    order = '>' if big_endian else '<'
    pts = np.c_[line, np.zeros(len(line))] if z else np.asarray(line, dtype=float)
    return (struct.pack(order + 'BII', 0 if big_endian else 1, 1002 if z else 2, len(pts))
            + np.asarray(pts, dtype=np.dtype(order + 'f8')).tobytes())


def _write_gpkg(path, lines, epsg):
    """A minimal GeoPackage, the way GDAL lays one out."""
    con = sqlite3.connect(path)
    con.executescript("""
        CREATE TABLE gpkg_spatial_ref_sys (srs_name TEXT, srs_id INTEGER PRIMARY KEY,
            organization TEXT, organization_coordsys_id INTEGER, definition TEXT,
            description TEXT);
        CREATE TABLE gpkg_contents (table_name TEXT PRIMARY KEY, data_type TEXT,
            identifier TEXT, srs_id INTEGER);
        CREATE TABLE gpkg_geometry_columns (table_name TEXT, column_name TEXT,
            geometry_type_name TEXT, srs_id INTEGER, z INTEGER, m INTEGER);
        CREATE TABLE "fault lines" (fid INTEGER PRIMARY KEY, geom BLOB, name TEXT);
    """)
    con.execute('INSERT INTO gpkg_spatial_ref_sys VALUES (?, ?, ?, ?, ?, ?)',
                ('utm', epsg, 'EPSG', epsg, pyproj.CRS.from_epsg(epsg).to_wkt(), ''))
    con.execute("INSERT INTO gpkg_contents VALUES ('fault lines', 'features', 'f', ?)", (epsg,))
    con.execute("INSERT INTO gpkg_geometry_columns VALUES "
                "('fault lines', 'geom', 'LINESTRING', ?, 0, 0)", (epsg,))
    for i, line in enumerate(lines):
        # header: magic, version, flags (little endian + an xy envelope), srs_id
        header = b'GP' + bytes([0, 0b0000_0011]) + struct.pack('<i', epsg)
        envelope = struct.pack('<4d', line[:, 0].min(), line[:, 0].max(),
                               line[:, 1].min(), line[:, 1].max())
        blob = header + envelope + _wkb_linestring(line)
        con.execute('INSERT INTO "fault lines" VALUES (?, ?, ?)', (i, blob, f'f{i}'))
    con.commit()
    con.close()
    return str(path)


@pytest.fixture
def gpkg_utm_file(tmp_path):
    return _write_gpkg(tmp_path / 'faults_utm.gpkg', [_to_utm(_trace_lonlat())], UTM)


ALL_FORMATS = ['geojson_file', 'geojson_legacy_utm_file', 'shp_utm_file', 'shp_noprj_file',
               'gmt_file', 'gmt_utm_file', 'kml_file', 'kmz_file', 'gpkg_utm_file']


@pytest.mark.parametrize('name', ALL_FORMATS)
def test_every_format_lands_on_the_fault(name, request, velocity_attrs):
    """Whatever the format and CRS, the trace ends up at FAULT_LON on the data grid."""
    path = request.getfixturevalue(name)
    traces = fault.load_faults([path], velocity_attrs, verbose=False)
    pts = np.vstack(traces)
    assert np.allclose(pts[:, 0], FAULT_LON, atol=1e-6), name
    assert pts[:, 1].min() < 30.0 and pts[:, 1].max() > 30.6     # untouched, not clipped


def test_points_are_not_traces(geojson_file, kml_file):
    for path in (geojson_file, kml_file):
        layers = fault.read_fault_file(path)
        assert sum(len(lines) for lines, _ in layers) == 1


def test_crs_is_read_from_the_file(shp_utm_file, gmt_utm_file, gpkg_utm_file,
                                   geojson_legacy_utm_file, shp_noprj_file):
    assert pyproj.CRS(fault.read_fault_file(shp_utm_file)[0][1]).to_epsg() == UTM
    assert fault.read_fault_file(gmt_utm_file)[0][1] == UTM
    assert fault.read_fault_file(gpkg_utm_file)[0][1] == UTM
    assert fault.read_fault_file(geojson_legacy_utm_file)[0][1] == UTM
    assert fault.read_fault_file(shp_noprj_file)[0][1] is None


def _utm_attrs():
    """A UTM 52N grid over the same area as the lat/lon fixtures."""
    e0, n0 = _to_utm(np.array([[LON0, LAT0]]))[0]
    return dict(LENGTH='330', WIDTH='430', X_FIRST=str(round(e0, -3)),
                Y_FIRST=str(round(n0, -3)), X_STEP='200.0', Y_STEP='-200.0',
                X_UNIT='meters', Y_UNIT='meters', EPSG=str(UTM))


def test_lonlat_file_on_a_utm_grid(geojson_file, kml_file):
    """WGS84 faults on a UTM map: re-projected into easting / northing."""
    atr = _utm_attrs()
    expect = _to_utm(_trace_lonlat())
    for path in (geojson_file, kml_file):
        pts = np.vstack(fault.load_faults([path], atr, verbose=False))
        assert np.allclose(pts, expect, atol=1e-3)


def test_file_without_crs_in_metres_is_taken_as_the_data_crs(tmp_path):
    """A UTM shapefile that lost its .prj, on a UTM grid: used as it is."""
    xy = _to_utm(_trace_lonlat())
    path = _write_shp(tmp_path / 'lost_prj', [xy])
    pts = np.vstack(fault.load_faults([path], _utm_attrs(), verbose=False))
    assert np.allclose(pts, xy)


def test_far_away_traces_are_dropped(tmp_path, velocity_attrs, capsys):
    path = tmp_path / 'far.gmt'
    path.write_text('>\n10.0 45.0\n10.5 45.5\n>\n' + '\n'.join(
        f'{x} {y}' for x, y in _trace_lonlat()) + '\n')
    traces = fault.load_faults([str(path)], velocity_attrs, verbose=True)
    assert len(traces) == 1
    assert '1 of 2 trace(s)' in capsys.readouterr().out


def test_antimeridian(tmp_path):
    """A lat/lon grid running past 180 E still meets faults written as -175."""
    atr = dict(LENGTH='100', WIDTH='200', X_FIRST='175.0', Y_FIRST='-10.0',
               X_STEP='0.05', Y_STEP='-0.05', Y_UNIT='degrees', X_UNIT='degrees')
    path = tmp_path / 'dateline.gmt'
    path.write_text('>\n-176.0 -11.0\n-176.0 -14.0\n')
    pts = np.vstack(fault.load_faults([str(path)], atr, verbose=False))
    assert np.allclose(pts[:, 0], 184.0)


def test_radar_coordinates_are_refused():
    atr = {'LENGTH': '10', 'WIDTH': '10'}
    with pytest.raises(ValueError, match='geocoded'):
        fault.load_faults(['whatever.shp'], atr, verbose=False)


def test_unknown_format_says_what_to_install(tmp_path):
    path = tmp_path / 'faults.fgb'
    path.write_bytes(b'not a flatgeobuf')
    with pytest.raises(Exception):         # noqa: B017 - OGR, or our own advice
        fault.read_fault_file(str(path))


##############################################################################
#                             crossings                                      #
##############################################################################
ROW, COL0, COL1 = 150, 75, 375     # lat 30.299, lon 130.151 -> 130.751


def test_crossing_matches_the_swath_distance_axis(velocity_attrs):
    """The crossing sits exactly where the swath puts pixels at FAULT_LON."""
    traces = [_trace_lonlat()]
    geom = line_geometry(velocity_attrs, [ROW, COL0], [ROW, COL1], 1.0)
    length_km = geom['length'] * geom['dist_scale']
    dists, pts = fault.profile_crossings(traces, velocity_attrs, [ROW, COL0], [ROW, COL1],
                                         length_km)
    assert dists.size == 1
    t = (FAULT_LON - (LON0 + (COL0 + 0.5) * DLON)) / ((COL1 - COL0) * DLON)
    assert dists[0] == pytest.approx(t * length_km, rel=1e-12)
    assert pts[0] == pytest.approx([FAULT_LON, LAT0 + (ROW + 0.5) * DLAT])

    # the two pixel columns either side of the fault straddle that distance
    swath = swath_extract(make_velocity(), velocity_attrs, (0, 0, WIDTH, LENGTH),
                          [ROW, COL0], [ROW, COL1], 0.05)
    col_fault = (FAULT_LON - LON0) / DLON - 0.5                 # 224.5
    near = (swath['row'] == ROW) & np.isclose(np.abs(swath['col'] - col_fault), 0.5)
    assert near.sum() == 2
    assert swath['dist'][near].mean() == pytest.approx(dists[0], abs=1e-4)


def test_direction_matters(velocity_attrs):
    traces = [_trace_lonlat()]
    fwd, _ = fault.profile_crossings(traces, velocity_attrs, [ROW, COL0], [ROW, COL1], 60.0)
    bwd, _ = fault.profile_crossings(traces, velocity_attrs, [ROW, COL1], [ROW, COL0], 60.0)
    assert fwd[0] + bwd[0] == pytest.approx(60.0)


def test_no_crossing(velocity_attrs):
    traces = [_trace_lonlat() + [0.5, 0.0]]                       # east of the line's end
    dists, pts = fault.profile_crossings(traces, velocity_attrs, [ROW, COL0], [ROW, COL1], 60.)
    assert dists.size == 0 and pts.shape == (0, 2)


def test_zigzag_and_shared_vertex(velocity_attrs):
    """Three crossings; a vertex ON the line is counted once, not twice."""
    lat_line = LAT0 + (ROW + 0.5) * DLAT
    zigzag = np.array([[130.30, 30.45], [130.34, 30.15], [130.38, 30.45], [130.42, 30.15]])
    vertex_on_line = np.array([[130.60, 30.45], [130.62, lat_line], [130.64, 30.15]])
    dists, pts = fault.profile_crossings([zigzag, vertex_on_line], velocity_attrs,
                                         [ROW, COL0], [ROW, COL1], 60.0)
    assert dists.size == 4
    assert np.all(np.diff(dists) > 0)
    assert pts[-1] == pytest.approx([130.62, lat_line])


def test_to_display_in_pixels(velocity_attrs):
    xs, ys = fault.to_display([_trace_lonlat()], velocity_attrs, 'yx')
    assert np.allclose(xs, (FAULT_LON - LON0) / DLON - 0.5)
    xs, ys = fault.to_display([_trace_lonlat(), _trace_lonlat()], velocity_attrs, 'geo')
    assert np.isnan(xs).sum() == 1                      # one NaN between the two traces


##############################################################################
#                          end to end                                        #
##############################################################################
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')


def _vlines(ax):
    """The x of every dashed vertical line drawn by axvline."""
    return sorted(float(ln.get_xdata()[0]) for ln in ax.get_lines()
                  if ln.get_linestyle() == '--' and len(ln.get_xdata()) == 2
                  and ln.get_xdata()[0] == ln.get_xdata()[1])


def _fault_texts(ax, text='Fault line'):
    return [t for t in ax.texts if t.get_text() == text]


def test_fault_on_the_map_and_the_profile(velocity_file, geojson_file, tmp_path):
    out = tmp_path / 'out'
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-fault', geojson_file,
                          '--save', '--nodisplay', '--outdir', str(out), '--noverbose'])
    sec = obj.sections[0]

    # the map: the trace is drawn, and did not stretch the map past the data
    assert obj.fault_obj is not None
    assert np.nanmean(obj.fault_obj.get_xdata()) == pytest.approx(FAULT_LON)
    ylo, yhi = sorted(obj.ax_img.get_ylim())
    assert ylo >= 30.0 - 0.01 and yhi <= 30.6 + 0.01

    # the profile: one dashed line at the crossing, annotated
    length = sec.swath_list[0]['length']
    t = (FAULT_LON - (LON0 + (COL0 + 0.5) * DLON)) / ((COL1 - COL0) * DLON)
    assert sec.crossings == pytest.approx([t * length])
    assert _vlines(sec.ax_swath) == pytest.approx([t * length])
    assert len(_fault_texts(sec.ax_swath)) == 1

    # the arctan step is centred on the fault: the median crosses zero there
    den = sec.den
    ok = np.isfinite(den['median_dw'])
    x0 = den['x_centers'][ok][np.argmin(np.abs(den['median_dw'][ok]))]
    assert x0 == pytest.approx(sec.crossings[0], abs=1.5)

    info = next(out.glob('*_info.txt')).read_text()
    assert f'fault crossing : {t * length:.4f} km' in info
    assert 'faults.geojson' in info


def test_fault_with_topography(velocity_file, shp_utm_file, dem_file, tmp_path):
    """A UTM shapefile + a UTM DEM on lat/lon data: the line is on both panels."""
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-fault', shp_utm_file,
                          '-topoxsect', dem_file, '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    sec = obj.sections[0]
    assert sec.crossings.size == 1
    assert _vlines(sec.ax_swath) == pytest.approx(list(sec.crossings))
    assert _vlines(sec.ax_topo) == pytest.approx(list(sec.crossings))
    assert len(_fault_texts(sec.ax_swath)) == 1 and not _fault_texts(sec.ax_topo)
    # the DEM ridge was built over the fault: the topography peaks at the line
    tden = sec.topo_den
    peak = tden['x_centers'][np.nanargmax(tden['median_dw'])]
    assert peak == pytest.approx(sec.crossings[0], abs=2.0)


def test_label_colour_and_x_offset(velocity_file, geojson_file, tmp_path):
    out = tmp_path / 'off'
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '--fault', geojson_file,
                          '--fault-label', 'Main fault', '--fault-color', 'darkred',
                          '--x-offset', '10', '--save', '--nodisplay',
                          '--outdir', str(out), '--noverbose'])
    sec = obj.sections[0]
    texts = _fault_texts(sec.ax_swath, 'Main fault')
    assert len(texts) == 1 and not _fault_texts(sec.ax_swath)
    assert matplotlib.colors.same_color(texts[0].get_color(), 'darkred')
    assert texts[0].xy[0] == pytest.approx(sec.crossings[0] - 10)
    assert pytest.approx(sec.crossings[0] - 10) in _vlines(sec.ax_swath)
    info = next(out.glob('*_info.txt')).read_text()
    assert f'fault crossing : {sec.crossings[0] - 10:.4f} km' in info


def test_empty_label_keeps_the_line(velocity_file, geojson_file, tmp_path):
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-fault', geojson_file,
                          '--fault-label', '', '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    sec = obj.sections[0]
    assert len(_vlines(sec.ax_swath)) == 1
    assert [t.get_text() for t in sec.ax_swath.texts] == ['A', "A'"]     # no fault label


def test_crossing_over_no_data_widens_the_axis(velocity_file, tmp_path):
    """A fault in the NaN strip at the top of the scene, past the last data: still shown."""
    path = tmp_path / 'ew.gmt'
    path.write_text('>\n130.0 30.585\n131.0 30.585\n')     # rows < 20 are all NaN
    obj = cli.swath_main([velocity_file, '--start-lalo', '30.10', '130.30',
                          '--end-lalo', '30.595', '130.30', '-w', '4', '-fault', str(path),
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    sec = obj.sections[0]
    assert sec.crossings.size == 1
    assert sec.crossings[0] > sec.den['x_edges'][-1]         # beyond the last data
    lo, hi = sec.ax_swath.get_xlim()
    assert lo < sec.crossings[0] < hi


def test_pixel_coordinate_map(velocity_file, geojson_file, tmp_path):
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-fault', geojson_file,
                          '--coord', 'yx', '--nodisplay', '--outdir', str(tmp_path),
                          '--noverbose'])
    assert np.nanmean(obj.fault_obj.get_xdata()) == pytest.approx((FAULT_LON - LON0) / DLON - 0.5)
    assert obj.sections[0].crossings.size == 1


def test_unreadable_file_does_not_stop_the_viewer(velocity_file, tmp_path, capsys):
    path = tmp_path / 'broken.geojson'
    path.write_text('{ not json')
    obj = cli.swath_main([velocity_file, *LINE, '-w', '6', '-fault', str(path),
                          '--nodisplay', '--outdir', str(tmp_path), '--noverbose'])
    assert 'cannot use --fault' in capsys.readouterr().out
    assert obj.sections[0].has_data and obj.sections[0].crossings.size == 0


def test_new_section_finds_its_own_crossings(velocity_file, geojson_file, tmp_path):
    """Sections B, C ... each get their own crossings; "n" clears them."""
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    from viewinsar.viewer import swathProfileViewer

    inps = cli.cmd_line_parse([velocity_file, '-w', '5', '-fault', geojson_file,
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

    def key(k):
        x, y = obj.ax_img.transData.transform((130.5, 30.3))
        obj.on_key(KeyEvent('key_press_event', obj.fig_map.canvas, k, x, y))

    click(130.15, 30.30)
    click(130.75, 30.30)                 # A: crosses the fault
    key('=')
    click(130.55, 30.20)
    click(130.85, 30.40)                 # B: east of it, no crossing
    a, b = obj.sections
    assert a.crossings.size == 1 and b.crossings.size == 0
    assert _vlines(a.ax_swath) and not _vlines(b.ax_swath)
    key('n')                             # clears the CURRENT section, B
    assert b.crossings.size == 0 and a.crossings.size == 1


##############################################################################
#                               the CLI                                      #
##############################################################################
def test_single_dash_spellings(velocity_file, geojson_file):
    inps = cli.cmd_line_parse([velocity_file, '-nlook', '3', '-svg', '-fault', geojson_file,
                               '-topoxsect', velocity_file, '--noverbose'])
    assert inps.nlook == 3 and inps.fig_ext == '.svg'
    assert inps.fault_file == [geojson_file] and inps.topoxsect == velocity_file
    assert inps.view_extra == []


def test_view_options_are_not_taken_for_ours(velocity_file):
    """No abbreviations: view.py's -n / --faultline must reach view.py untouched."""
    inps = cli.cmd_line_parse([velocity_file, '-n', '0', '--faultline', 'f.shp',
                               '--noverbose'])
    assert inps.nlook == 1 and inps.fault_file is None
    assert inps.view_extra == ['-n', '0', '--faultline', 'f.shp']


def test_several_fault_files_then_more_options(velocity_file, geojson_file, kml_file):
    inps = cli.cmd_line_parse([velocity_file, '-fault', geojson_file, kml_file, '-w', '4',
                               '-c', 'RdBu', '--noverbose'])
    assert inps.fault_file == [geojson_file, kml_file]
    assert inps.width == 4.0 and inps.view_extra == ['-c', 'RdBu']


@pytest.mark.parametrize('bad', [['-fault', 'no_such_file.shp'], ['--fault-lw', '0']])
def test_bad_fault_options(velocity_file, geojson_file, bad):
    args = [velocity_file, *bad, '--noverbose']
    if bad[0] == '--fault-lw':
        args += ['-fault', geojson_file]
    with pytest.raises(ValueError):
        cli.cmd_line_parse(args)


def test_wkb_reader():
    """ISO and extended WKB, both byte orders, Z, multi-part and polygons."""
    import struct as st

    from viewinsar.wkb import wkb_to_geojson

    line = [[130.0, 30.0], [130.5, 30.5], [131.0, 30.2]]
    for big in (False, True):
        for z in (False, True):
            geom = wkb_to_geojson(_wkb_linestring(line, big_endian=big, z=z))
            assert geom['type'] == 'LineString' and geom['coordinates'] == line
    two = st.pack('<BII', 1, 5, 2) + _wkb_linestring(line) + _wkb_linestring(line[:2])
    geom = wkb_to_geojson(two)
    assert geom['type'] == 'MultiLineString' and len(geom['coordinates'][1]) == 2
    ring = line + [line[0]]
    poly = st.pack('<BIII', 1, 3, 1, len(ring)) + np.asarray(ring, '<f8').tobytes()
    assert wkb_to_geojson(poly)['coordinates'][0][-1] == line[0]
    ewkb = st.pack('<BII', 1, 0x20000002, 4326) + st.pack('<I', 2) + np.asarray(
        line[:2], '<f8').tobytes()                     # PostGIS: with an SRID
    assert wkb_to_geojson(ewkb)['coordinates'] == line[:2]
    point = st.pack('<BI', 1, 1) + np.asarray([1.5, 2.5], '<f8').tobytes()
    assert wkb_to_geojson(point) == {'type': 'Point', 'coordinates': [1.5, 2.5]}

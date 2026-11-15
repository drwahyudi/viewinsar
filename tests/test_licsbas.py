"""LiCSBAS and LiCSAR products: the line-of-sight geometry from their east,
north and up (GEOCml*/E.geo ... and {frame}.geo.E.tif), found by itself; the
mm/yr of a LiCSBAS velocity; and `viewinsar model` on a LiCSBAS folder with
nothing but the file and the fault."""

import json
import os

import numpy as np
import pytest

from viewinsar import cli
from viewinsar import licsbas as lb
from viewinsar.gnss import LOSGeometry, los_unit_vector

W, L, D = 300, 260, 0.004
LAT_C, LON_C = -0.75, 120.70              # centre of the first pixel (corner_lat / lon)
FLON, FLAT, STRIKE = 121.30, -1.27, 110.0
S_TRUE, D1_TRUE = 20.0, 12.0              # left-lateral
KM_LAT = 111.32
KM_LON = KM_LAT * np.cos(np.radians(FLAT))


def grid():
    lon = LON_C + np.arange(W) * D
    lat = LAT_C - np.arange(L) * D
    return np.meshgrid(lon, lat)


LON_MAX = LON_C + (W - 1) * D
LAT_MID = LAT_C - (L - 1) * D / 2


def angles(lon, lat):
    """Ascending: the incidence grows west to east (31 to 45), the azimuth near 100.5."""
    inc = 31.0 + 14.0 * (np.asarray(lon) - LON_C) / (LON_MAX - LON_C)
    az = 100.5 + 0.3 * (np.asarray(lat) - LAT_MID)
    return inc, az


def fault_x(lon, lat):
    """km from the fault, positive to the right of the strike (SSW)."""
    phi = np.radians(STRIKE)
    east, north = (np.asarray(lon) - FLON) * KM_LON, (np.asarray(lat) - FLAT) * KM_LAT
    return east * np.cos(phi) - north * np.sin(phi)


def strike_vec():
    phi = np.radians(STRIKE)
    return np.sin(phi), np.cos(phi)


def write_dem_par(folder, width, length, lat_c, lon_c, dlat, dlon):
    with open(os.path.join(folder, 'EQA.dem_par'), 'w') as fid:
        fid.write('Gamma DIFF&GEO DEM/MAP parameter file\ntitle:\tsynthetic\n'
                  'DEM_projection:     EQA\ndata_format:        REAL*4\n'
                  f'width:                 {width}\nnlines:                {length}\n'
                  f'corner_lat:    {lat_c:.8f}  decimal degrees\n'
                  f'corner_lon:   {lon_c:.8f}  decimal degrees\n'
                  f'post_lat:   {dlat:.8e} decimal degrees\n'
                  f'post_lon:    {dlon:.8e} decimal degrees\n\nellipsoid_name: WGS 84\n')


def write_geoc(folder, dtype='<f4'):
    """GEOCml*/E.geo, N.geo, U.geo and EQA.dem_par, as LiCSBAS02 writes them."""
    os.makedirs(folder, exist_ok=True)
    lon2, lat2 = grid()
    ue, un, uu = los_unit_vector(*angles(lon2, lat2))
    for comp, arr in zip('ENU', (ue, un, uu)):
        arr.astype(dtype).tofile(os.path.join(folder, f'{comp}.geo'))
    write_dem_par(folder, W, L, LAT_C, LON_C, -D, D)
    return folder


def licsbas_tree(root):
    """A LiCSBAS work folder: GEOCml10GACOS/ (the geometry), a decoy GEOCml20/
    (other angles, another grid) and TS_GEOCml10GACOS/results/vel.filt.mskd.geo.tif
    - the LOS velocity in mm/yr with no unit, across a left-lateral fault."""
    rasterio = pytest.importorskip('rasterio')
    from rasterio.transform import from_origin

    write_geoc(os.path.join(root, 'GEOCml10GACOS'))
    decoy = os.path.join(root, 'GEOCml20')
    os.makedirs(decoy)
    for comp, val in zip('ENU', (0.6, -0.1, 0.7937)):
        np.full((L // 2, W // 2), val, '<f4').tofile(os.path.join(decoy, f'{comp}.geo'))
    write_dem_par(decoy, W // 2, L // 2, LAT_C, LON_C, -2 * D, 2 * D)

    lon2, lat2 = grid()
    ue, un, _ = los_unit_vector(*angles(lon2, lat2))
    sx, sy = strike_vec()
    v_s = (S_TRUE / np.pi) * np.arctan(fault_x(lon2, lat2) / D1_TRUE)
    los = v_s * (ue * sx + un * sy) + 4.0 + np.random.default_rng(5).normal(0, 1.0, v_s.shape)
    los[:10, :25] = np.nan
    results = os.path.join(root, 'TS_GEOCml10GACOS', 'results')
    os.makedirs(results)
    vel = os.path.join(results, 'vel.filt.mskd.geo.tif')
    with rasterio.open(vel, 'w', driver='GTiff', height=L, width=W, count=1, dtype='float32',
                       crs='EPSG:4326', nodata=np.nan,
                       transform=from_origin(LON_C - D / 2, LAT_C + D / 2, D, D)) as dst:
        dst.write(los.astype(np.float32), 1)

    return vel, write_fault(root)


def write_fault(root):
    """The trace: a straight line along the strike, well past the data."""
    sx, sy = strike_vec()
    t = np.linspace(-120, 120, 41)
    coords = [[FLON + v * sx / KM_LON, FLAT + v * sy / KM_LAT] for v in t]
    fault = os.path.join(root, 'fault.geojson')
    with open(fault, 'w') as fid:
        json.dump({'type': 'FeatureCollection', 'features': [{
            'type': 'Feature', 'properties': {},
            'geometry': {'type': 'LineString', 'coordinates': coords}}]}, fid)
    return fault


PROFILE = ['--start-lalo', '-0.80', '121.10', '--end-lalo', '-1.70', '121.45', '-w', '20']
FAST = ['--walkers', '32', '--steps', '1500', '--seed', '1']


##############################################################################
#                       E, N, U and the EQA.dem_par grid                     #
##############################################################################
@pytest.mark.parametrize('inc, az', [(39.0, 100.5), (36.0, -100.0), (44.0, 81.0)])
def test_enu_to_angles_inverts_the_unit_vector(inc, az):
    e, n, u = los_unit_vector(inc, az)
    got_inc, got_az = lb.enu_to_angles(e, n, u)
    assert got_inc == pytest.approx(inc) and got_az == pytest.approx(az)
    # ascending (az ~100): the satellite to the west - the LOS points west and up
    if az == 100.5:
        assert e < 0 and u > 0


def test_enu_to_angles_no_data():
    inc, az = lb.enu_to_angles([0.0, np.nan], [0.0, 0.1], [0.0, 0.8])
    assert np.isnan(inc).all() and np.isnan(az).all()


def test_dem_par_grid_is_the_pixel_centre(tmp_path):
    write_dem_par(tmp_path, 10, 8, -0.75, 120.70, -0.004, 0.004)
    atr = lb.grid_attrs(lb.read_dem_par(tmp_path / 'EQA.dem_par'))
    assert atr['WIDTH'] == '10' and atr['LENGTH'] == '8'
    # corner_lat / corner_lon: the centre of the first pixel; Y_FIRST: its edge
    assert float(atr['Y_FIRST']) == pytest.approx(-0.748)
    assert float(atr['X_FIRST']) == pytest.approx(120.698)
    assert float(atr['Y_STEP']) == pytest.approx(-0.004)


@pytest.mark.parametrize('dtype', ['<f4', '>f4'])
def test_read_enu_geo_either_byte_order(tmp_path, dtype):
    folder = write_geoc(str(tmp_path / 'GEOCml10'), dtype=dtype)
    east, north, up, atr = lb.read_enu_geo(folder)
    lon2, lat2 = grid()
    inc, az = lb.enu_to_angles(east, north, up)
    want_inc, want_az = angles(lon2, lat2)
    assert np.allclose(inc, want_inc, atol=1e-3) and np.allclose(az, want_az, atol=1e-3)
    assert atr['EPSG'] == '4326'


def test_read_enu_geo_checks_the_size(tmp_path):
    folder = write_geoc(str(tmp_path / 'GEOCml10'))
    write_dem_par(folder, W + 1, L, LAT_C, LON_C, -D, D)
    with pytest.raises(ValueError, match='bytes'):
        lb.read_enu_geo(folder)


def test_tif_triples(tmp_path):
    for comp in 'ENU':
        (tmp_path / f'105A_09094_131313.geo.{comp}.tif').write_bytes(b'')
    (tmp_path / 'SLOPE.tif').write_bytes(b'')                   # ends in E.tif: no set
    triple = lb.enu_tif_triple(str(tmp_path / '105A_09094_131313.geo.N.tif'))
    assert [os.path.basename(p) for p in triple] == [
        '105A_09094_131313.geo.E.tif', '105A_09094_131313.geo.N.tif',
        '105A_09094_131313.geo.U.tif']
    assert lb.enu_tif_triple(str(tmp_path / 'SLOPE.tif')) is None
    assert lb.find_enu_tifs(str(tmp_path)) == triple


##############################################################################
#                             found by themselves                            #
##############################################################################
def test_candidates_prefer_the_geoc_of_the_time_series(tmp_path):
    work = tmp_path / 'work'
    write_geoc(str(work / 'GEOCml20'))
    write_geoc(str(work / 'GEOCml10GACOS'))
    (work / 'GEOC').mkdir()
    for comp in 'ENU':
        (work / 'GEOC' / f'frame.geo.{comp}.tif').write_bytes(b'')
    data = work / 'TS_GEOCml10GACOS' / 'results' / 'vel.filt.mskd.geo.tif'
    data.parent.mkdir(parents=True)
    data.write_bytes(b'')
    cands = [os.path.basename(p) for p in lb.geometry_candidates(str(data))]
    assert cands == ['GEOCml10GACOS', 'GEOCml20', 'GEOC']
    assert lb.in_licsbas_tree(str(data))
    # data copied out of the tree: nothing is guessed
    outside = tmp_path / 'elsewhere' / 'vel.tif'
    outside.parent.mkdir()
    outside.write_bytes(b'')
    assert lb.geometry_candidates(str(outside)) == [] and not lb.in_licsbas_tree(str(outside))


def test_licsbas_unit(tmp_path):
    work = tmp_path / 'work'
    geoc = write_geoc(str(work / 'GEOCml10'))
    res = work / 'TS_GEOCml10' / 'results'
    res.mkdir(parents=True)
    assert lb.licsbas_unit(str(res / 'vel.filt.mskd.geo.tif')) == 'mm/year'
    assert lb.licsbas_unit(str(res / 'vstd.geo.tif')) == 'mm/year'
    assert lb.licsbas_unit(str(res / 'coh_avg.geo.tif')) is None
    elsewhere = str(tmp_path / 'vel.tif')
    assert lb.licsbas_unit(elsewhere) is None
    assert lb.licsbas_unit(elsewhere, geom_hint=geoc) == 'mm/year'       # --geom GEOCml10
    assert lb.licsbas_unit(elsewhere, geom_hint=[geoc]) == 'mm/year'


def test_los_geometry_from_any_path(tmp_path):
    folder = write_geoc(str(tmp_path / 'GEOCml10'))
    lon2, lat2 = grid()
    want_inc, want_az = angles(lon2, lat2)
    for path in (folder, os.path.join(folder, 'E.geo'), os.path.join(folder, 'EQA.dem_par')):
        geom = LOSGeometry.from_path(path)
        assert geom.kind == 'enu' and geom.inc.shape == (L, W)
        inc, az = geom.at_lonlat(lon2[100, [10, 200]], lat2[100, [10, 200]])
        assert inc == pytest.approx(want_inc[100, [10, 200]], abs=1e-3)
        assert az == pytest.approx(want_az[100, [10, 200]], abs=1e-3)
    with pytest.raises(ValueError, match='not a geometry'):
        LOSGeometry.from_path(str(tmp_path / 'nothing.dat'))


def test_los_geometry_from_licsar_tifs(tmp_path):
    rasterio = pytest.importorskip('rasterio')
    from rasterio.transform import from_origin

    lon2, lat2 = grid()
    for comp, arr in zip('ENU', los_unit_vector(*angles(lon2, lat2))):
        with rasterio.open(tmp_path / f'105A.geo.{comp}.tif', 'w', driver='GTiff', height=L,
                           width=W, count=1, dtype='float32', crs='EPSG:4326',
                           transform=from_origin(LON_C - D / 2, LAT_C + D / 2, D, D)) as dst:
            dst.write(arr.astype(np.float32), 1)
    geom = LOSGeometry.from_path(str(tmp_path / '105A.geo.U.tif'))
    inc, az = geom.at_lonlat(lon2[5, [3, 250]], lat2[5, [3, 250]])
    want_inc, want_az = angles(lon2, lat2)
    assert inc == pytest.approx(want_inc[5, [3, 250]], abs=1e-3)
    assert az == pytest.approx(want_az[5, [3, 250]], abs=1e-3)


def test_geometry_that_misses_the_data_is_skipped(tmp_path):
    vel, _ = licsbas_tree(str(tmp_path / 'work'))
    from viewinsar.compat import read_attribute
    atr = read_attribute(vel)
    geom = LOSGeometry.for_file(vel, atr, verbose=False)
    assert geom.source.startswith(str(tmp_path / 'work' / 'GEOCml10GACOS'))
    far = dict(atr, X_FIRST=str(float(atr['X_FIRST']) + 5.0))      # 5 degrees east
    assert not geom.covers(far)


##############################################################################
#                                end to end                                  #
##############################################################################
def test_model_on_a_licsbas_folder(tmp_path):
    """The file and the fault, nothing else: the unit (mm/yr), the LOS of every
    pixel (GEOCml10GACOS, not the decoy) and --comp los are found."""
    vel, fault = licsbas_tree(str(tmp_path / 'work'))
    out = tmp_path / 'out'
    obj = cli.model_main([vel, '-fault', fault, *PROFILE, '--strain', *FAST, '--nodisplay',
                          '--noverbose', '--outdir', str(out)])
    assert obj.disp_unit == 'mm/year' and obj.model_unit == 'mm/yr'
    assert obj.comp == 'los'
    assert obj._los[0].source.startswith(str(tmp_path / 'work' / 'GEOCml10GACOS'))
    sec = obj.sections[0]
    s = sec.result['summary']
    assert s['S']['map'] == pytest.approx(-S_TRUE, abs=1.5)          # S < 0: left-lateral
    assert s['d1']['map'] == pytest.approx(D1_TRUE, abs=2.0)

    # every pixel turned into fault-parallel velocity with its own LOS
    pts = np.genfromtxt(next(out.glob('*_points.csv')), delimiter=',', names=True)
    sx, sy = strike_vec()
    want = (S_TRUE / np.pi) * np.arctan(fault_x(pts['lon'], pts['lat']) / D1_TRUE)
    ue, un, _ = los_unit_vector(*angles(pts['lon'], pts['lat']))
    offset = 4.0 / np.median(ue * sx + un * sy)
    # along the oriented strike (290 deg) the left-lateral motion is -v_s
    assert np.median(pts['v_model'] + want) == pytest.approx(-offset, abs=0.5)

    # the parameters as Hussain et al. write them
    lines = sec.parameter_lines()
    assert lines[0].startswith('slip rate = 2') and ' mm/yr (' in lines[0]
    assert lines[0].endswith('95% CI), left-lateral')
    assert lines[1].startswith('locking depth = ') and ' km (' in lines[1]
    assert lines[-1].startswith('strain rate = ') and 'µstrain/yr' in lines[-1]
    assert sec.ax_swath.get_ylabel() == 'Fault-parallel velocity [mm/yr]'
    assert 'µstrain/yr' in sec.ax_strain.get_ylabel()
    obj.close()


def test_data_unit_option(tmp_path):
    """A velocity without a unit, out of its LiCSBAS folder: --data-unit."""
    vel, fault = licsbas_tree(str(tmp_path / 'work'))
    lone = tmp_path / 'velocity_asc.tif'
    lone.write_bytes(open(vel, 'rb').read())
    args = [str(lone), '-fault', fault, *PROFILE, '--comp', 'as-is', *FAST, '--nodisplay',
            '--noverbose', '--outdir', str(tmp_path / 'out')]
    obj = cli.model_main(args)
    # no unit, out of its LiCSBAS folder: its numbers (tens) say mm/yr
    assert obj.model_unit == 'mm/yr' and obj.disp_unit == 'mm/year'
    assert obj._unit_guess[str(lone)][0] == 'mm/year'
    obj = cli.model_main(args + ['--data-unit', 'mm/yr'])
    assert not obj._unit_guess
    assert obj.disp_unit == 'mm/year' and obj.model_unit == 'mm/yr'
    assert obj.sections[0].parameter_lines()[0].split(' (')[0].endswith('mm/yr')
    with pytest.raises(ValueError, match='unknown unit'):
        cli.model_cmd_line_parse([str(lone), '--data-unit', 'furlongs'])


def test_comp_auto_without_geometry_fits_as_is(tmp_path):
    """No geometry anywhere: --comp auto fits the data as they are."""
    vel, fault = licsbas_tree(str(tmp_path / 'work'))
    lone = tmp_path / 'velocity_asc.tif'
    lone.write_bytes(open(vel, 'rb').read())
    obj = cli.model_main([str(lone), '-fault', fault, *PROFILE, '--data-unit', 'mm/yr',
                          *FAST, '--nodisplay', '--noverbose', '--outdir', str(tmp_path)])
    assert obj.comp == 'as-is'
    # ... and --geom points at the geometry: los again, one --geom per file
    obj = cli.model_main([str(lone), '-fault', fault, *PROFILE, '--data-unit', 'mm/yr',
                          '--geom', str(tmp_path / 'work' / 'GEOCml10GACOS'), *FAST,
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path)])
    assert obj.comp == 'los'
    assert obj.sections[0].result['summary']['S']['map'] == pytest.approx(-S_TRUE, abs=1.5)
    # ... or one angle pair per file
    obj = cli.model_main([str(lone), '-fault', fault, *PROFILE, '--data-unit', 'mm/yr',
                          '--los-inc', '38', '--los-az', '100.5', *FAST, '--nodisplay',
                          '--noverbose', '--outdir', str(tmp_path)])
    assert obj.comp == 'los' and obj._los[0].constant


def test_swath_gnss_in_the_los_of_a_licsbas_file(tmp_path):
    """-vel los (the default) takes the LOS of each station from E, N, U."""
    vel, _ = licsbas_tree(str(tmp_path / 'work'))
    gnss = tmp_path / 'gnss.csv'
    gnss.write_text('site,lon,lat,ve,vn,se,sn\nA,121.0,-0.9,10,0,1,1\nB,121.6,-1.5,0,10,1,1\n')
    obj = cli.swath_main([vel, '-gnss', str(gnss), '--start-lalo', '-0.8', '121.0',
                          '--end-lalo', '-1.6', '121.6', '-w', '10', '--nodisplay',
                          '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert obj.disp_unit == 'mm/year'
    assert 'GEOCml10GACOS' in obj.gnss_note
    ue, un, _ = los_unit_vector(*angles(np.array([121.0, 121.6]), np.array([-0.9, -1.5])))
    assert obj.gnss_value == pytest.approx([10 * ue[0], 10 * un[1]], abs=0.05)
    obj.close()


def test_mintpy_without_gdal_reads_geotiff_with_rasterio(tmp_path, monkeypatch):
    from viewinsar import compat

    if not compat.have_mintpy():
        assert not compat.mintpy_reads('x.tif') and not compat.mintpy_reads('x.h5')
        return
    monkeypatch.setattr(compat, '_gdal_importable', lambda: False)
    assert not compat.mintpy_reads('vel.geo.tif') and compat.mintpy_reads('velocity.h5')
    monkeypatch.setattr(compat, '_gdal_importable', lambda: True)
    assert compat.mintpy_reads('vel.geo.tif')


##############################################################################
#                       a file that states a wrong unit                      #
##############################################################################
def write_mislabelled(path, unit):
    """A MintPy velocity file whose values are mm/yr, labelled `unit`: an
    arctan step of +/-10 mm/yr across a north-south line at 130.45 E."""
    h5py = pytest.importorskip('h5py')
    length, width, step = 120, 150, 0.004
    lon = 130.15 + (np.arange(width) + 0.5) * step
    vel = np.tile(20.0 / np.pi * np.arctan((lon - 130.45) * 111.32 / 5.0), (length, 1))
    with h5py.File(path, 'w') as fid:
        fid.create_dataset('velocity', data=vel.astype(np.float32))
        fid.attrs.update({'FILE_TYPE': 'velocity', 'UNIT': unit, 'LENGTH': str(length),
                          'WIDTH': str(width), 'Y_FIRST': '30.50', 'X_FIRST': '130.15',
                          'Y_STEP': str(-step), 'X_STEP': str(step), 'Y_UNIT': 'degrees',
                          'X_UNIT': 'degrees', 'REF_Y': '60', 'REF_X': '75'})
    return str(path)


@pytest.mark.parametrize('stated, factor', [('m/year', 1000.0), ('cm/year', 10.0)])
def test_data_unit_corrects_a_wrong_unit(tmp_path, stated, factor):
    """mm/yr labelled m/year (or cm/year): the map, the profiles and the GNSS
    were `factor` x off; --data-unit mm/yr puts them right."""
    vel = write_mislabelled(tmp_path / 'geo_velocity.h5', stated)
    gnss = tmp_path / 'gnss.csv'
    gnss.write_text('site,lon,lat,ve,vn,vu,se,sn,su\nG1,130.30,30.30,0,0,8,1,1,1\n')
    common = [vel, '--start-lalo', '30.30', '130.20', '--end-lalo', '30.30', '130.70',
              '-w', '5', '-gnss', str(gnss), '-vel', 'vu', '--nodisplay', '--noverbose',
              '--outdir', str(tmp_path / 'out')]
    wrong = cli.swath_main(common)
    right = cli.swath_main(common + ['--data-unit', 'mm/yr'])
    assert right.disp_unit == wrong.disp_unit                   # the same display unit
    to_disp = 1.0 / {'cm/year': 10.0, 'mm/year': 1.0}[str(right.disp_unit)]
    # the profile: the file's numbers in mm/yr, in the display unit
    import h5py
    with h5py.File(vel, 'r') as fid:
        raw = float(np.abs(fid['velocity'][()]).max())
    amp = np.nanmax(np.abs(right.data_list[0]))
    assert amp == pytest.approx(raw * to_disp, rel=1e-4)
    assert np.nanmax(np.abs(wrong.data_list[0])) == pytest.approx(amp * factor, rel=1e-4)
    # the map and its colour range, scaled the same
    assert np.nanmax(np.abs(right.data_img)) == pytest.approx(
        np.nanmax(np.abs(wrong.data_img)) / factor, rel=1e-4)
    assert right.vlim[1] == pytest.approx(wrong.vlim[1] / factor, rel=1e-4)
    # the GNSS stay what they are: 8 mm/yr up
    assert right.gnss_value[0] == pytest.approx(8.0 * to_disp)
    assert wrong.gnss_value[0] == pytest.approx(8.0 * to_disp)
    right.close()
    wrong.close()


def test_data_unit_must_convert(tmp_path):
    vel = write_mislabelled(tmp_path / 'geo_velocity.h5', 'radian')
    with pytest.raises(ValueError, match='does not convert'):
        cli.swath_main([vel, '--data-unit', 'mm/yr', '--nodisplay', '--noverbose',
                        '--start-lalo', '30.30', '130.20', '--end-lalo', '30.30', '130.70'])


def test_map_follows_the_unit_the_file_says(tmp_path):
    """A velocity labelled mm/year: view.py takes any "velocity" as m/year; the
    map is put back in mm/year, as the profiles read it."""
    vel = write_mislabelled(tmp_path / 'geo_velocity.h5', 'mm/year')
    obj = cli.swath_main([vel, '--start-lalo', '30.30', '130.20', '--end-lalo', '30.30',
                          '130.70', '-w', '5', '--nodisplay', '--noverbose',
                          '--outdir', str(tmp_path / 'out')])
    import h5py
    with h5py.File(vel, 'r') as fid:
        raw = float(np.abs(fid['velocity'][()]).max())
    to_disp = 1.0 / {'cm/year': 10.0, 'mm/year': 1.0}[str(obj.disp_unit)]
    assert np.nanmax(np.abs(obj.data_list[0])) == pytest.approx(raw * to_disp, rel=1e-4)
    assert np.nanmax(np.abs(obj.data_img)) == pytest.approx(raw * to_disp, rel=1e-4)
    obj.close()


##############################################################################
#                     LiCSBAS's own files, opened directly                   #
##############################################################################
def licsbas_native(root):
    """TS_GEOCml10GACOS/cum_filt.h5 (vel, cum, the grid as LiCSBAS writes it),
    results/vel.filt.mskd and results/mask (float32, no header), next to
    GEOCml10GACOS/ (E, N, U, EQA.dem_par)."""
    h5py = pytest.importorskip('h5py')
    write_geoc(os.path.join(root, 'GEOCml10GACOS'))
    lon2, lat2 = grid()
    ue, un, _ = los_unit_vector(*angles(lon2, lat2))
    sx, sy = strike_vec()
    v_s = (S_TRUE / np.pi) * np.arctan(fault_x(lon2, lat2) / D1_TRUE)
    vel = v_s * (ue * sx + un * sy) + 4.0 + np.random.default_rng(6).normal(0, 1.0, v_s.shape)
    mask = np.ones(vel.shape, dtype=np.float32)
    mask[:40, :60] = 0.0                                   # LiCSBAS: 0 = masked
    ts = os.path.join(root, 'TS_GEOCml10GACOS')
    os.makedirs(os.path.join(ts, 'results'))
    dates = [20180101, 20190101, 20200101]
    years = np.array([0.0, 1.0, 2.0])
    with h5py.File(os.path.join(ts, 'cum_filt.h5'), 'w') as h5:
        h5['vel'] = vel.astype(np.float32)
        h5['vintercept'] = np.zeros(vel.shape, np.float32)
        h5['cum'] = (years[:, None, None] * vel[None]).astype(np.float32)
        h5['imdates'] = np.array(dates, dtype=np.int32)
        h5['corner_lat'], h5['corner_lon'] = LAT_C, LON_C
        h5['post_lat'], h5['post_lon'] = -D, D
        h5['refarea'] = b'10:11/10:11'
    vel_mskd = np.where(mask > 0, vel, np.nan).astype('<f4')
    vel_mskd.tofile(os.path.join(ts, 'results', 'vel.filt.mskd'))
    mask.astype('<f4').tofile(os.path.join(ts, 'results', 'mask'))
    return (os.path.join(ts, 'cum_filt.h5'), os.path.join(ts, 'results', 'vel.filt.mskd'),
            write_fault(root))


def test_licsbas_files_are_recognised(tmp_path):
    h5, raw, _ = licsbas_native(str(tmp_path / 'work'))
    assert lb.kind_of(h5)[0] == 'h5'
    kind, par = lb.kind_of(raw)
    assert kind == 'raw' and par.endswith(os.path.join('GEOCml10GACOS', 'EQA.dem_par'))
    names = lb.slices(h5)
    assert names[0] == 'vel' and 'cum-20200101' in names and 'vintercept' in names
    atr = lb.attributes(h5)
    assert atr['UNIT'] == 'mm/year' and float(atr['Y_FIRST']) == pytest.approx(LAT_C + D / 2)
    assert lb.attributes(h5, 'cum-20190101')['UNIT'] == 'mm'
    data, _ = lb.read(h5, 'cum-20200101', box=(5, 6, 15, 16))
    vel, _ = lb.read(h5, 'vel', box=(5, 6, 15, 16))
    assert np.allclose(data, 2.0 * vel)
    raw_vel = lb.read(raw, box=(100, 100, 110, 110))[0]
    assert np.allclose(raw_vel, vel_at(h5, (100, 100, 110, 110)), equal_nan=True)
    assert lb.mask_for(h5).endswith(os.path.join('results', 'mask'))
    assert lb.mask_for(raw).endswith(os.path.join('results', 'mask'))
    # MintPy's own reader is never asked: it knows neither their grid nor their unit
    from viewinsar import compat
    assert not compat.mintpy_reads(h5) and not compat.mintpy_reads(raw)
    assert compat.get_slice_list(h5)[0] == 'vel'
    assert cli.default_dset(h5, print_msg=False) == 'vel'


def vel_at(h5, box):
    return lb.read(h5, 'vel', box=box)[0]


def test_list_licsbas(tmp_path, capsys):
    h5, _, _ = licsbas_native(str(tmp_path / 'work'))
    cli.main(['list', h5])
    out = capsys.readouterr().out.splitlines()
    assert out[0] == 'vel  <- default' and 'cum-20180101' in out


@pytest.mark.parametrize('which', ['h5', 'raw'])
def test_model_on_licsbas_files(tmp_path, which):
    """cum_filt.h5 or results/vel.filt.mskd as they are: the velocity in mm/yr,
    LiCSBAS's mask, the LOS of every pixel and --comp los, found by themselves."""
    h5, raw, fault = licsbas_native(str(tmp_path / 'work'))
    obj = cli.model_main([h5 if which == 'h5' else raw, '-fault', fault, *PROFILE, *FAST,
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert obj.model_unit == 'mm/yr' and obj.comp == 'los'
    assert obj._los[0].source.startswith(str(tmp_path / 'work' / 'GEOCml10GACOS'))
    # LiCSBAS's mask: the corner it masks is not used (NaN in vel.filt.mskd, or
    # results/mask found for cum_filt.h5)
    msk = obj.msk_list[0]
    off = np.isnan(obj.data_list[0]) | (~msk if msk is not None else False)
    assert off[:40, :60].all() and not off[100:, 100:].any()
    assert (msk is not None) == (which == 'h5')
    s = obj.sections[0].result['summary']
    assert s['S']['map'] == pytest.approx(-S_TRUE, abs=1.5)
    assert s['d1']['map'] == pytest.approx(D1_TRUE, abs=2.0)
    obj.close()


def test_swath_on_a_licsbas_cube_epoch(tmp_path):
    """--dset cum-YYYYMMDD: an epoch of the cumulative displacement, in mm."""
    h5, _, _ = licsbas_native(str(tmp_path / 'work'))
    obj = cli.swath_main([h5, '--dset', 'cum-20200101', '--start-lalo', '-0.80', '121.10',
                          '--end-lalo', '-1.70', '121.45', '-w', '5', '--nodisplay',
                          '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert obj.disp_unit == 'mm'
    vel = vel_at(h5, (0, 0, W, L))
    assert np.nanmax(np.abs(obj.data_list[0])) == pytest.approx(
        2.0 * np.nanmax(np.abs(np.where(np.isnan(obj.data_list[0]), np.nan, vel))), rel=1e-4)
    obj.close()


##############################################################################
#        a file that says nothing about itself: no unit, no geometry         #
##############################################################################
def write_bare_tif(path, scale=1e-3):
    """A descending LOS velocity GeoTIFF with no unit and no geometry anywhere
    near it - as ALOS-2 ScanSAR velocities are passed around - in `scale` x
    mm/yr (1e-3: m/yr), masked pixels 0."""
    rasterio = pytest.importorskip('rasterio')
    from rasterio.transform import from_origin

    lon2, lat2 = grid()
    ue, un, _ = los_unit_vector(38.0, -101.7)
    sx, sy = strike_vec()
    v_s = (S_TRUE / np.pi) * np.arctan(fault_x(lon2, lat2) / D1_TRUE)
    los = v_s * (ue * sx + un * sy) + 4.0 + np.random.default_rng(7).normal(0, 1.0, v_s.shape)
    los[:10, :25] = 0.0
    with rasterio.open(str(path), 'w', driver='GTiff', height=L, width=W, count=1,
                       dtype='float32', crs='EPSG:4326',
                       transform=from_origin(LON_C - D / 2, LAT_C + D / 2, D, D)) as dst:
        dst.write((los * scale).astype(np.float32), 1)
    return str(path)


BARE_GNSS = ('site,lon,lat,ve,vn,vu\nA,121.15,-0.95,20,5,\nB,121.30,-1.30,25,8,2\n'
             'C,121.40,-1.60,30,10,\n')


@pytest.mark.parametrize('scale, unit', [(1e-3, 'm/year'), (1.0, 'mm/year')])
def test_guess_velocity_unit(scale, unit):
    from viewinsar.compat import guess_velocity_unit

    v = np.random.default_rng(1).normal(0, 8, 10_000) * scale       # mm/yr, in `scale`
    v[:100], v[100:200] = 0.0, np.nan                               # masked: not counted
    got, p99 = guess_velocity_unit(v)
    assert got == unit and p99 == pytest.approx(np.percentile(np.abs(v[200:]), 99))
    assert guess_velocity_unit([0.0, np.nan]) == (None, None)


@pytest.mark.parametrize('name, dset, yes', [
    ('vel.filt.mskd.geo.tif', None, True), ('geo_velocity.h5', None, True),
    ('slip_rate.tif', None, True), ('scansardesc_masked.tif', None, False),
    ('scansardesc_masked.tif', 'velocity', True), ('elevation.tif', None, False),
    ('coherence.tif', None, False)])
def test_velocity_name(name, dset, yes):
    from viewinsar.compat import velocity_name

    assert velocity_name(name, dset) is yes


def test_nominal_heading_of_sentinel1():
    from viewinsar.gnss import heading_to_azimuth, nominal_heading

    asc, desc = nominal_heading(0.0, True, 98.18, 6.77), nominal_heading(0.0, False, 98.18, 6.77)
    assert asc == pytest.approx(-12.0, abs=0.5) and desc == pytest.approx(-168.0, abs=0.5)
    assert heading_to_azimuth(asc) == pytest.approx(102.0, abs=0.5)
    assert heading_to_azimuth(desc) == pytest.approx(-102.0, abs=0.5)
    assert nominal_heading(50.0, True, 98.18, 6.77) < asc - 2      # away from north, poleward


@pytest.mark.parametrize('rel, orbit, mission', [
    ('scansardesc_masked.tif', 'desc', 'ALOS-2'),
    ('alos2_asc_vel.tif', 'asc', 'ALOS-2'),
    ('105D_09094_131313.geo.vel.tif', 'desc', 'Sentinel-1'),
    ('S1_descending_vel.tif', 'desc', 'Sentinel-1'),
    (os.path.join('T125A', 'velocity.tif'), 'asc', None),          # the folder says it
    ('cascadia_vel.tif', None, None),
    ('velocity_2d.tif', None, None),
    ('asc_desc_up.tif', None, None)])
def test_nominal_geometry_from_the_name(tmp_path, rel, orbit, mission):
    from viewinsar.gnss import nominal_geometry

    atr = {'Y_FIRST': '-0.75', 'Y_STEP': '-0.004', 'LENGTH': '260'}
    got = nominal_geometry(str(tmp_path / rel), atr)
    if orbit is None:
        assert got is None
        return
    inc, az, what = got
    assert (az > 0) is (orbit == 'asc') and 95 < abs(az) < 105
    assert inc == pytest.approx(37.5 if mission == 'Sentinel-1' else 38.0)
    assert (mission or 'polar-orbiting') in what
    # the attributes say it too: ORBIT_DIRECTION beats the name
    other = 'ASCENDING' if orbit == 'desc' else 'DESCENDING'
    assert (nominal_geometry(str(tmp_path / rel), {**atr, 'ORBIT_DIRECTION': other})[1] > 0) \
        is (other == 'ASCENDING')


def test_los_geometry_nominal_last(tmp_path):
    atr = {'Y_FIRST': '-0.75', 'Y_STEP': '-0.004', 'LENGTH': '260'}
    geom = LOSGeometry.for_file(str(tmp_path / 'scansardesc_masked.tif'), atr, verbose=False)
    assert geom.kind == 'nominal' and geom.constant and geom.az < 0
    with pytest.raises(ValueError, match='ascending or descending'):
        LOSGeometry.for_file(str(tmp_path / 'velocity.tif'), atr, verbose=False)
    # the attributes' one value comes first
    geom = LOSGeometry.for_file(str(tmp_path / 'scansardesc_masked.tif'),
                                {**atr, 'CENTER_INCIDENCE_ANGLE': '33', 'HEADING': '-170'},
                                verbose=False)
    assert geom.kind == 'attributes' and geom.inc == 33


def _nominal_uvec():
    from viewinsar.gnss import heading_to_azimuth, nominal_heading

    lat = (LAT_C + D / 2) - D * L / 2
    az = float(heading_to_azimuth(nominal_heading(lat, False, 97.92, 6.85)))
    return los_unit_vector(38.0, az)


def test_model_on_a_bare_geotiff_with_gnss(tmp_path, capsys):
    """scansardesc_masked.tif, in m/yr, says nothing: no option needed - the
    unit from its numbers, the GNSS in a nominal LOS, the LOS fitted as is."""
    tif = write_bare_tif(tmp_path / 'scansardesc_masked.tif')
    fault = write_fault(str(tmp_path))
    gnss = tmp_path / 'GNSS_obs.csv'
    gnss.write_text(BARE_GNSS)
    obj = cli.model_main([tif, '-fault', fault, '-gnss', str(gnss), '--gnss-nofit', *PROFILE,
                          *FAST, '--nodisplay', '--outdir', str(tmp_path / 'out')])
    out = capsys.readouterr().out
    assert 'read in m/year' in out and 'not a velocity' not in out
    assert obj.disp_unit == 'mm/year' and obj.model_unit == 'mm/yr'
    assert 'nominal angles instead, for ALOS-2, descending, right-looking' in out
    assert obj.comp == 'as-is' and obj.gnss_comp == 'los'
    gin = obj.sections[0].gnss_in
    ue, un, uu = _nominal_uvec()
    want = np.array([20 * ue + 5 * un, 25 * ue + 8 * un + 2 * uu, 30 * ue + 10 * un])
    assert gin['index'].size == 3 and gin['value'] == pytest.approx(want, abs=0.05)
    # the LOS step across the fault, in mm/yr: S x (LOS . strike)
    sx, sy = strike_vec()
    lu = los_unit_vector(38.0, -101.7)
    s_los = abs(obj.sections[0].result['summary']['S']['map'])
    assert s_los == pytest.approx(S_TRUE * abs(lu[0] * sx + lu[1] * sy), abs=2.0)
    lines = obj.sections[0].parameter_lines()
    assert lines[0].startswith('slip rate in the LOS = ') and lines[0].split(' (')[0].endswith(
        'mm/yr') and lines[-1].startswith('strain rate in the LOS = ')
    # never turned into fault-parallel velocity with nominal angles
    with pytest.raises(ValueError, match='only nominal angles'):
        cli.model_main([tif, '-fault', fault, *PROFILE, '--comp', 'los', *FAST, '--nodisplay',
                        '--noverbose', '--outdir', str(tmp_path / 'out')])


def test_swath_on_a_bare_geotiff_with_gnss(tmp_path):
    tif = write_bare_tif(tmp_path / 'scansardesc_masked.tif')
    gnss = tmp_path / 'GNSS_obs.csv'
    gnss.write_text(BARE_GNSS)
    obj = cli.swath_main([tif, '-gnss', str(gnss), '--start-lalo', '-0.8', '121.1',
                          '--end-lalo', '-1.7', '121.45', '-w', '20', '--nodisplay',
                          '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert obj.disp_unit == 'mm/year' and 'nominal angles' in obj.gnss_note
    # the map is in mm/yr as well - the same numbers as the profiles
    assert np.nanmedian(obj.data_img) == pytest.approx(np.nanmedian(obj.data_list[0]), rel=1e-4)
    assert 2 < np.nanmedian(obj.data_list[0]) < 8                    # mm/yr, not m/yr
    ue, un, uu = _nominal_uvec()
    assert obj.gnss_value == pytest.approx(
        [20 * ue + 5 * un, 25 * ue + 8 * un + 2 * uu, 30 * ue + 10 * un], abs=0.05)
    obj.close()


def test_no_unit_guessed_for_what_is_not_a_velocity(tmp_path):
    from viewinsar.compat import is_unitless

    tif = write_bare_tif(tmp_path / 'coherence.tif')
    obj = cli.swath_main([tif, '--start-lalo', '-0.8', '121.1', '--end-lalo', '-1.7', '121.45',
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert not obj._unit_guess and is_unitless(obj.disp_unit)
    obj.close()


##############################################################################
#        --geom FOLDER: the geometry GeoTIFFs of LiCSBAS / LiCSAR in it      #
##############################################################################
def write_tif(path, arr, lon_c=LON_C, lat_c=LAT_C, step=D):
    rasterio = pytest.importorskip('rasterio')
    from rasterio.transform import from_origin

    with rasterio.open(str(path), 'w', driver='GTiff', height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype='float32', crs='EPSG:4326',
                       transform=from_origin(lon_c - step / 2, lat_c + step / 2, step, step)) as dst:
        dst.write(arr.astype(np.float32), 1)
    return str(path)


def write_enu_tifs(folder, names=('{c}.geo.tif',), lon_c=LON_C, lat_c=LAT_C):
    """E, N, U GeoTIFFs of the test geometry, named after `names` ({c}: E N U)."""
    os.makedirs(folder, exist_ok=True)
    lon2, lat2 = grid()
    lon2, lat2 = lon2 - LON_C + lon_c, lat2 - LAT_C + lat_c
    enu = los_unit_vector(*angles(lon2 - lon_c + LON_C, lat2 - lat_c + LAT_C))
    for pattern in names:
        for comp, arr in zip('ENU', enu):
            write_tif(os.path.join(folder, pattern.format(c=comp)), arr, lon_c, lat_c)
    return folder


@pytest.mark.parametrize('pattern', ['{c}.geo.tif', '105A_09094_131313.geo.{c}.tif', '{c}.tif',
                                     'los_{lc}.tif', 'los_{word}.tif', 'LOS_{WORD}.tiff'])
def test_enu_tif_names(tmp_path, pattern):
    words = {'E': 'east', 'N': 'north', 'U': 'up'}
    for c in 'ENU':
        name = pattern.format(c=c, lc=c.lower(), word=words[c], WORD=words[c].upper())
        (tmp_path / name).write_bytes(b'')
    (tmp_path / 'hgt.geo.tif').write_bytes(b'')
    (tmp_path / 'vel.filt.mskd.geo.tif').write_bytes(b'')
    found = lb.all_enu_tifs(str(tmp_path))
    assert len(found) == 1 and len(found[0]) == 3
    names = [os.path.basename(p).lower() for p in found[0]]
    assert ('e' in names[0] or 'east' in names[0]) and ('u' in names[2] or 'up' in names[2])


def test_angle_tif_pairs(tmp_path):
    for name in ('geo_incidenceAngle.tif', 'geo_azimuthAngle.tif', 'inc.tif', 'az.tif',
                 'princess.tif', 'cum_inc_only.tif'):
        (tmp_path / name).write_bytes(b'')
    pairs = [tuple(os.path.basename(p) for p in pair) for pair in lb.all_angle_tifs(str(tmp_path))]
    assert sorted(pairs) == [('geo_incidenceAngle.tif', 'geo_azimuthAngle.tif'),
                             ('inc.tif', 'az.tif')]
    assert lb.is_angle_tif(str(tmp_path / 'az.tif')) == (str(tmp_path / 'inc.tif'),
                                                          str(tmp_path / 'az.tif'))


def test_geometry_in_folder_order(tmp_path):
    root = tmp_path / 'work'
    write_geoc(str(root / 'GEOCml10'))
    write_enu_tifs(str(root / 'GEOC'), names=('105A_09094_131313.geo.{c}.tif',))
    os.makedirs(root / 'GEOC' / '20180101_20180113')              # interferograms: skipped
    write_enu_tifs(str(root / 'GEOC' / '20180101_20180113'))
    found = lb.geometry_in_folder(str(root))
    rel = [os.path.relpath(p, root) for p in found]
    assert rel == ['GEOCml10', os.path.join('GEOC', '105A_09094_131313.geo.E.tif')]


def test_geom_folder_takes_the_geometry_that_covers_the_data(tmp_path):
    vel, _ = licsbas_tree(str(tmp_path / 'work'))
    from viewinsar.compat import read_attribute

    geoms = tmp_path / 'geometry'
    write_enu_tifs(str(geoms / 'other_frame'), lon_c=LON_C + 5.0)    # 5 degrees east
    write_enu_tifs(str(geoms / 'this_frame'), names=('frame.geo.{c}.tif',))
    geom = LOSGeometry.for_file(vel, read_attribute(vel), geom_file=str(geoms), verbose=False)
    assert 'this_frame' in geom.source
    with pytest.raises(ValueError, match='none of its geometries covers'):
        LOSGeometry.for_file(vel, read_attribute(vel), geom_file=str(geoms / 'other_frame'),
                             verbose=False)
    with pytest.raises(ValueError, match='no LOS geometry in it'):
        os.makedirs(tmp_path / 'empty')
        LOSGeometry.for_file(vel, read_attribute(vel), geom_file=str(tmp_path / 'empty'),
                             verbose=False)


@pytest.mark.parametrize('kind', ['enu', 'angles'])
def test_model_with_a_geom_folder_of_tifs(tmp_path, kind):
    """A LiCSBAS velocity out of its folder, and --geom a folder of geometry
    GeoTIFFs: --comp auto turns the LOS into fault-parallel velocity."""
    vel, fault = licsbas_tree(str(tmp_path / 'work'))
    lone = tmp_path / 'vel_mskd.tif'
    lone.write_bytes(open(vel, 'rb').read())
    geoms = tmp_path / 'my_geometry'
    if kind == 'enu':
        write_enu_tifs(str(geoms / 'GEOC'), names=('105A_09094_131313.geo.{c}.tif',))
    else:
        os.makedirs(geoms)
        lon2, lat2 = grid()
        inc, az = angles(lon2, lat2)
        write_tif(geoms / 'incidenceAngle.tif', inc)
        write_tif(geoms / 'azimuthAngle.tif', az)
    obj = cli.model_main([str(lone), '--fault', fault, *PROFILE, '--geom', str(geoms), *FAST,
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'out')])
    assert obj.comp == 'los' and obj.model_unit == 'mm/yr'
    assert obj._los[0].kind == ('enu' if kind == 'enu' else 'angles')
    assert obj.sections[0].result['summary']['S']['map'] == pytest.approx(-S_TRUE, abs=1.5)
    assert 'left-lateral' in obj.sections[0].parameter_lines()[0]

"""`viewinsar model`: elastic dislocation fits of profiles across a fault.

The core - the forward model, the bins, the noise covariance, the ensemble
sampler and the fit - on profiles of known slip rate and locking depth; then
the whole command on synthetic products: LOS velocities of an ascending and a
descending track across a right-lateral fault striking N80E (20 mm/yr below
12 km), an east-velocity map, a fault trace and GNSS.
"""

import json
import os

import numpy as np
import pytest
from matplotlib import colors as mcolors

from viewinsar import cli
from viewinsar import model as mdl
from viewinsar.modelview import (
    azimuth_of,
    ci_text,
    fault_normal,
    local_en,
    orient_strike,
    slip_sense,
    strike_of_points,
)

S_TRUE, D1_TRUE = 20.0, 12.0                 # mm/yr, km
FLON, FLAT, STRIKE = 130.45, 30.30, 80.0
KM_LAT = 111.32
KM_LON = KM_LAT * np.cos(np.radians(FLAT))
PROFILE = ['--start-lalo', '30.58', '130.40', '--end-lalo', '30.02', '130.52', '-w', '15']
FAST = ['--walkers', '32', '--steps', '1500', '--seed', '1']


##############################################################################
#                              the products                                  #
##############################################################################
def _fault_frame(lon, lat):
    """km from the fault, perpendicular to it (+ to the SSE), for lon / lat."""
    east, north = (lon - FLON) * KM_LON, (lat - FLAT) * KM_LAT
    normal = np.array([np.cos(np.radians(STRIKE)), -np.sin(np.radians(STRIKE))])
    return east * normal[0] + north * normal[1]


def _vpar_n80e(x):
    """The velocity along N80E: right-lateral, so the SSE block moves WSW."""
    return -(S_TRUE / np.pi) * np.arctan(x / D1_TRUE)


def _write(path, datasets, attrs):
    import h5py

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with h5py.File(path, 'w') as fid:
        for name, data in datasets.items():
            fid.create_dataset(name, data=np.asarray(data, dtype=np.float32))
        for key, val in attrs.items():
            fid.attrs[key] = val
    return path


def _attrs(lat0, lon0, length, width, step, **extra):
    atr = {'FILE_TYPE': 'velocity', 'UNIT': 'm/year', 'DATA_TYPE': 'float32',
           'LENGTH': str(length), 'WIDTH': str(width), 'Y_FIRST': str(lat0),
           'X_FIRST': str(lon0), 'Y_STEP': str(-step), 'X_STEP': str(step),
           'Y_UNIT': 'degrees', 'X_UNIT': 'degrees', 'REF_Y': '10', 'REF_X': '10',
           'WAVELENGTH': '0.055465763', 'PLATFORM': 'Sen', 'PROCESSOR': 'isce',
           'EPSG': '4326', 'ALOOKS': '1', 'RLOOKS': '1'}
    atr.update({k: str(v) for k, v in extra.items()})
    return atr


def _creep_n80e(x, rate, depth):
    """Shallow creep, in the sense of the slip, along N80E (mm/yr)."""
    return -rate * (np.heaviside(x, 0.5) - 0.5 - np.arctan(x / depth) / np.pi)


def make_track(folder, inc0, az, offset, lat0=30.60, lon0=130.0, seed=0, heading=-12.0,
               creep=None, up_step=0.0):
    """A LOS velocity (m/year) and its geometryGeo.h5, in `folder`.

    creep   - (rate, depth): shallow creep on the fault, mm/yr and km
    up_step - mm/yr: the SSE side moves up by that much relative to the NNW
              side - a vertical step at the fault, which the LOS sees
    """
    length, width, step = 150, 250, 0.004
    lon = lon0 + (np.arange(width) + 0.5) * step
    lat = lat0 - (np.arange(length) + 0.5) * step
    lon2, lat2 = np.meshgrid(lon, lat)
    inc = inc0 + 4.0 * (np.arange(width) / width - 0.5)[None, :] * np.ones((length, 1))
    ue = -np.sin(np.radians(inc)) * np.sin(np.radians(az))
    un = np.sin(np.radians(inc)) * np.cos(np.radians(az))
    p_e, p_n = np.sin(np.radians(STRIKE)), np.cos(np.radians(STRIKE))
    x = _fault_frame(lon2, lat2)
    vpar = _vpar_n80e(x) + (_creep_n80e(x, *creep) if creep else 0.0)
    los = (vpar * (ue * p_e + un * p_n) + up_step * (np.heaviside(x, 0.5) - 0.5)
           * np.cos(np.radians(inc)) + offset)
    los += np.random.default_rng(seed).normal(0, 1.0, los.shape)
    atr = _attrs(lat0, lon0, length, width, step, HEADING=heading,
                 CENTER_INCIDENCE_ANGLE=inc0,
                 ORBIT_DIRECTION='ASCENDING' if abs(heading) < 90 else 'DESCENDING')
    vfile = _write(os.path.join(folder, 'geo_velocity.h5'), {'velocity': los / 1000.0}, atr)
    _write(os.path.join(folder, 'geo_geometryGeo.h5'),
           {'incidenceAngle': inc, 'azimuthAngle': np.full(inc.shape, az),
            'height': 500 + 300 * np.exp(-(_fault_frame(lon2, lat2) / 20) ** 2)},
           dict(atr, FILE_TYPE='geometry', UNIT='1'))
    return vfile


@pytest.fixture(scope='module')
def products(tmp_path_factory):
    pytest.importorskip('h5py')
    root = tmp_path_factory.mktemp('model')
    out = {'asc': make_track(str(root / 'asc'), 39.0, 102.0, 2.0, seed=1),
           'desc': make_track(str(root / 'desc'), 37.0, -102.0, -3.0, lat0=30.62,
                              lon0=130.10, seed=2, heading=-168.0)}

    # an east-velocity map: the fault-parallel motion seen east, + 4 mm/yr
    length, width, step = 150, 250, 0.004
    lon = 130.0 + (np.arange(width) + 0.5) * step
    lat = 30.60 - (np.arange(length) + 0.5) * step
    lon2, lat2 = np.meshgrid(lon, lat)
    ve = _vpar_n80e(_fault_frame(lon2, lat2)) * np.sin(np.radians(STRIKE)) + 4.0
    ve += np.random.default_rng(3).normal(0, 1.0, ve.shape)
    out['east'] = _write(str(root / 'east' / 'geo_ve.h5'), {'velocity': ve / 1000.0},
                         _attrs(30.60, 130.0, length, width, step))

    # the trace: a straight line along the strike, well past the map
    t = np.linspace(-90, 90, 31)
    coords = [[FLON + v * np.sin(np.radians(STRIKE)) / KM_LON,
               FLAT + v * np.cos(np.radians(STRIKE)) / KM_LAT] for v in t]
    out['fault'] = str(root / 'fault.geojson')
    with open(out['fault'], 'w') as fid:
        json.dump({'type': 'FeatureCollection', 'features': [{
            'type': 'Feature', 'properties': {},
            'geometry': {'type': 'LineString', 'coordinates': coords}}]}, fid)

    # GNSS along the profile (Eurasia-fixed-like: + 5 mm/yr east), one far off it
    rows = ['site,lon,lat,ve,vn,se,sn']
    sites = [(130.40, 30.56), (130.42, 30.48), (130.44, 30.40), (130.455, 30.33),
             (130.47, 30.25), (130.48, 30.18), (130.50, 30.10), (130.52, 30.04),
             (130.85, 30.50)]
    for i, (lon_g, lat_g) in enumerate(sites):
        vp = _vpar_n80e(_fault_frame(lon_g, lat_g))
        rows.append(f'S{i:02d},{lon_g},{lat_g},{vp * np.sin(np.radians(STRIKE)) + 5:.3f},'
                    f'{vp * np.cos(np.radians(STRIKE)):.3f},0.5,0.5')
    out['gnss'] = str(root / 'gnss.csv')
    with open(out['gnss'], 'w') as fid:
        fid.write('\n'.join(rows) + '\n')
    return out


def _csv(path):
    return np.genfromtxt(path, delimiter=',', names=True, dtype=None, encoding='utf-8')


##############################################################################
#                              the forward model                             #
##############################################################################
def test_forward_far_field_and_shapes():
    x = np.array([-1e6, -12.0, 0.0, 12.0, 1e6])
    v = mdl.forward(x, 20.0, 12.0, 3.0)
    assert v.shape == x.shape
    assert v[0] == pytest.approx(3.0 - 10.0, abs=1e-3)
    assert v[-1] == pytest.approx(3.0 + 10.0, abs=1e-3)
    assert v[2] == pytest.approx(3.0)
    assert v[3] - 3.0 == pytest.approx(5.0)                   # atan(1) = pi / 4
    assert v[1] - 3.0 == pytest.approx(-(v[3] - 3.0))
    # K parameter sets at once -> (K, x.size)
    many = mdl.forward(x, np.array([20.0, 10.0]), np.array([12.0, 5.0]), np.zeros(2))
    assert many.shape == (2, x.size)
    assert many[0] == pytest.approx(v - 3.0)


def test_creep_is_a_step_that_dies_away():
    x = np.array([-1e5, -1e-6, 1e-6, 1e5])
    creep = mdl.forward(x, 20.0, 12.0, 0.0, creep=6.0, d2=3.0) - mdl.forward(x, 20.0, 12.0)
    assert creep[0] == pytest.approx(0.0, abs=1e-3) and creep[-1] == pytest.approx(0.0, abs=1e-3)
    assert creep[2] - creep[1] == pytest.approx(6.0, abs=1e-3)
    # in the sense of S: a negative slip rate steps down
    neg = mdl.forward(x, -20.0, 12.0, 0.0, creep=6.0, d2=3.0) - mdl.forward(x, -20.0, 12.0)
    assert neg[2] - neg[1] == pytest.approx(-6.0, abs=1e-3)


def test_rotation_on_one_side():
    x = np.array([-50.0, 50.0])
    base = mdl.forward(x, 20.0, 12.0)
    end = mdl.forward(x, 20.0, 12.0, rot=0.02) - base
    start = mdl.forward(x, 20.0, 12.0, rot=0.02, rot_side=-1) - base
    assert end == pytest.approx([0.0, 1.0])
    assert start == pytest.approx([-1.0, 0.0])


def test_model_parameters_and_priors():
    model = mdl.ProfileModel(n_sets=3, creep=True, free_x0=(-2, 2), gnss=True)
    assert model.names == ['S', 'd1', 'C', 'd2', 'x0', 'a', 'a2', 'a3', 'ag']
    assert mdl.describe('a3') == 'offset of file 3' and mdl.describe('d1') == 'locking depth'
    theta = np.array([[20, 12, 5, 3, 0, 0, 0, 0, 0],
                      [20, 12, 5, 13, 0, 0, 0, 0, 0],       # d2 > d1
                      [20, 12, 5, 3, 3, 0, 0, 0, 0],        # x0 out of bounds
                      [99, 12, 5, 3, 0, 0, 0, 0, 0]], dtype=float)
    lp = model.log_prior(theta)
    assert lp[0] == 0.0 and np.all(np.isinf(lp[1:]))
    pri = mdl.ProfileModel(d1_prior=(17.0, 5.0))
    assert pri.log_prior(np.array([20.0, 22.0, 0.0]))[0] == pytest.approx(-0.5)


def test_bins_take_the_median():
    rng = np.random.default_rng(0)
    x = rng.uniform(-10, 10, 20000)
    v = 2.0 * np.floor(x) + rng.normal(0, 0.5, x.size)
    bins = mdl.bin_profile(x, v, 1.0)
    assert bins['x'].size == 20
    assert bins['median'] == pytest.approx(2.0 * np.floor(bins['x']), abs=0.05)
    assert bins['std'] == pytest.approx(0.5, abs=0.08)
    assert np.all(bins['sem'] < bins['std'] / 10)
    with pytest.raises(ValueError, match='bin'):
        mdl.bin_profile(x[:5], v[:5], 1.0)


def test_noise_covariance_is_recovered():
    rng = np.random.default_rng(4)
    x = np.arange(0, 300, 1.0)
    phi = np.exp(-1.0 / 8.0)                         # exponential, length 8 km
    r = np.zeros(x.size)
    for i in range(1, x.size):
        r[i] = phi * r[i - 1] + rng.normal(0, np.sqrt(1 - phi ** 2) * 2.0)
    sill, length = mdl.fit_noise_covariance(x, r, 1.0)
    assert 2.0 < sill < 6.0 and 4.0 < length < 16.0
    white = mdl.fit_noise_covariance(x, rng.normal(0, 1, x.size), 1.0)
    assert white[1] < 1.0
    cov = mdl.exponential_covariance([0, 8, 16], 4.0, 8.0, nugget=1.0)
    assert cov[0, 0] == pytest.approx(5.0) and cov[0, 1] == pytest.approx(4.0 / np.e)


##############################################################################
#                                the sampler                                 #
##############################################################################
def test_ensemble_sampler_draws_a_gaussian():
    cov = np.array([[4.0, 1.8], [1.8, 1.0]])
    prec = np.linalg.inv(cov)
    mean = np.array([3.0, -1.0])

    def log_prob(theta):
        d = np.atleast_2d(theta) - mean
        return -0.5 * np.einsum('ij,jk,ik->i', d, prec, d)

    rng = np.random.default_rng(0)
    p0 = mean + 0.01 * rng.standard_normal((32, 2))
    chain, lp, accept = mdl.ensemble_sample(log_prob, p0, 3000, rng=rng)
    assert chain.shape == (3000, 32, 2) and lp.shape == (3000, 32)
    assert 0.3 < accept < 0.9
    flat = chain[500:].reshape(-1, 2)
    assert flat.mean(axis=0) == pytest.approx(mean, abs=0.15)
    assert np.cov(flat.T) == pytest.approx(cov, rel=0.15, abs=0.15)
    assert np.all(mdl.autocorr_time(chain[500:]) < 100)
    with pytest.raises(ValueError, match='walkers'):
        mdl.ensemble_sample(log_prob, p0[:3], 10)


def test_autocorrelation_time():
    rng = np.random.default_rng(1)
    white = rng.normal(size=(4000, 8, 1))
    assert mdl.autocorr_time(white)[0] == pytest.approx(1.0, abs=0.3)
    phi, ar = 0.9, np.zeros((4000, 8, 1))
    for i in range(1, 4000):
        ar[i] = phi * ar[i - 1] + rng.normal(size=(8, 1))
    assert mdl.autocorr_time(ar)[0] == pytest.approx((1 + phi) / (1 - phi), rel=0.3)


##############################################################################
#                                  the fit                                   #
##############################################################################
def _profile(rng, n=40000, half=80.0, noise=1.5, **kw):
    x = rng.uniform(-half, half, n)
    return x, mdl.forward(x, kw.pop('S', S_TRUE), kw.pop('d1', D1_TRUE), **kw) \
        + rng.normal(0, noise, n)


def test_fit_recovers_slip_rate_and_locking_depth():
    rng = np.random.default_rng(10)
    x, v = _profile(rng, a=3.0)
    res = mdl.fit_profile((x, v), walkers=32, steps=2000, seed=2, verbose=False)
    s = res['summary']
    assert s['S']['map'] == pytest.approx(S_TRUE, abs=1.0)
    assert s['d1']['map'] == pytest.approx(D1_TRUE, abs=2.0)
    assert s['a']['map'] == pytest.approx(3.0, abs=0.3)
    for name, true in (('S', S_TRUE), ('d1', D1_TRUE)):
        assert s[name]['lo'] < true < s[name]['hi']
    assert res['samples'].shape[1] == 3 and res['names'] == ['S', 'd1', 'a']
    gamma = res['derived']['strain']
    assert gamma['map'] == pytest.approx(s['S']['map'] / (np.pi * s['d1']['map']))
    # the bounds are the priors: nothing outside
    lo, hi = res['model'].lo_hi()
    assert np.all(res['samples'] >= lo) and np.all(res['samples'] <= hi)
    # the band holds the MAP, the strain rate peaks at the fault
    xx = np.linspace(-50, 50, 101)
    best, band_lo, band_hi = mdl.model_band(res, xx)
    assert np.all(band_lo <= best + 1e-6) and np.all(best <= band_hi + 1e-6)
    rate, r_lo, r_hi = mdl.strain_rate(res, xx)
    assert np.argmax(rate) == 50 and rate[50] == pytest.approx(gamma['map'])


def test_fit_creep_offsets_gnss_and_fault_position():
    rng = np.random.default_rng(11)
    kw = dict(x0=1.3, creep=6.0, d2=3.0)
    set1 = _profile(rng, n=30000, half=60.0, noise=1.0, a=2.0, **kw)
    set2 = _profile(rng, n=30000, half=60.0, noise=1.0, a=-4.0, **kw)
    gx = np.array([-55, -30, -15, -6, 8, 20, 40, 58.0])
    gnss = dict(x=gx, v=mdl.forward(gx, S_TRUE, D1_TRUE, 9.0, x0=1.3) + rng.normal(0, 0.3, 8),
                sigma=np.full(8, 0.3))
    res = mdl.fit_profile([set1, set2], gnss=gnss, creep=True, free_x0=(-3.0, 3.0),
                          walkers=40, steps=3000, seed=3, verbose=False)
    s = res['summary']
    assert res['names'] == ['S', 'd1', 'C', 'd2', 'x0', 'a', 'a2', 'ag']
    assert s['S']['map'] == pytest.approx(S_TRUE, abs=1.5)
    assert s['d1']['map'] == pytest.approx(D1_TRUE, abs=3.0)
    assert s['C']['map'] == pytest.approx(6.0, abs=1.5)
    assert s['x0']['map'] == pytest.approx(1.3, abs=0.6)
    assert s['a']['map'] - s['a2']['map'] == pytest.approx(6.0, abs=0.3)
    assert s['ag']['map'] == pytest.approx(9.0, abs=0.6)
    assert np.all(res['samples'][:, 3] <= res['samples'][:, 1])        # d2 <= d1
    assert len(res['noise']) == 2 and res['sets'] == [0, 1]


def test_fit_prior_noise_modes_and_errors():
    rng = np.random.default_rng(12)
    x, v = _profile(rng, n=20000, half=30.0, noise=2.0)
    # a tight prior on d1 pulls it there
    free = mdl.fit_profile((x, v), walkers=24, steps=800, verbose=False)
    res = mdl.fit_profile((x, v), d1_prior=(25.0, 0.5), walkers=24, steps=800, verbose=False)
    assert res['summary']['d1']['median'] > free['summary']['d1']['median'] + 3.0
    white = mdl.fit_profile((x, v), noise='white', walkers=24, steps=600, verbose=False)
    assert white['noise'][0][0] == 0.0 and white['noise_mode'] == 'white'
    given = mdl.fit_profile((x, v), noise=(4.0, 10.0), walkers=24, steps=600, verbose=False)
    assert given['noise'][0] == (4.0, 10.0) and given['noise_mode'] == 'given'
    with pytest.raises(ValueError, match='burn'):
        mdl.fit_profile((x, v), steps=100, burn=100, verbose=False)
    with pytest.raises(ValueError, match='too few bins'):
        mdl.fit_profile((x[:4], v[:4]), verbose=False)


##############################################################################
#                               the geometry                                 #
##############################################################################
def test_fault_offset_of_two_lines():
    """Hussain et al. (2016): a line each side, each with its own slope, and
    their offset at the fault - its spread the standard error of the two
    intercepts (2 sigma / sqrt(n) each, for x uniform over one side)."""
    rng = np.random.default_rng(4)
    x_a, x_b = rng.uniform(-5, 0, 3000), rng.uniform(0, 5, 3000)
    v_a = -1.0 + 0.4 * x_a + rng.normal(0, 1.0, x_a.size)
    v_b = 6.0 + 0.1 * x_b + rng.normal(0, 1.0, x_b.size)
    fit = mdl.fault_offset(x_a, v_a, x_b, v_b, n_boot=1000, seed=1)
    assert fit['offset'] == pytest.approx(7.0, abs=0.2)
    assert fit['lo'] < 7.0 < fit['hi'] and fit['lo'] < fit['offset'] < fit['hi']
    assert fit['line_a'][0] == pytest.approx(0.4, abs=0.06)
    assert fit['line_b'][0] == pytest.approx(0.1, abs=0.06)
    assert fit['boot'].shape == (1000,) and fit['n_a'] == fit['n_b'] == 3000
    assert fit['std'] == pytest.approx(np.hypot(2 / np.sqrt(3000), 2 / np.sqrt(3000)), rel=0.15)
    # the same seed, the same draws - in chunks or not
    again = mdl.fault_offset(x_a, v_a, x_b, v_b, n_boot=1000, seed=1)
    assert np.array_equal(again['boot'], fit['boot'])
    # too few points on a side
    assert mdl.fault_offset(x_a[:2], v_a[:2], x_b, v_b) is None


def test_split_offsets_into_horizontal_and_vertical():
    """The LOS offsets of an ascending and a descending track as one fault-
    parallel and one vertical offset - draw by draw too; none from two alike."""
    h, u = np.array([-0.63, 0.56]), np.array([0.78, 0.80])
    los = h * 4.0 + u * 1.5
    par, up = mdl.split_offsets(los, h, u)
    assert par[0] == pytest.approx(4.0) and up[0] == pytest.approx(1.5)
    draws = np.vstack([los[0] + np.arange(5.0) * h[0], los[1] + np.arange(5.0) * h[1]])
    par, up = mdl.split_offsets(draws, h, u)
    assert par == pytest.approx(4.0 + np.arange(5.0)) and up == pytest.approx(np.full(5, 1.5))
    assert mdl.split_offsets(los, [0.60, 0.62], [0.78, 0.79]) is None
    assert mdl.split_offsets(los[:1], h[:1], u[:1]) is None


def test_parameters_written_as_before():
    """The parameters as the box writes them: '20.1 mm/yr (19.4 – 22.2, 95% CI)'
    - a true minus sign, an angle in degrees."""
    assert ci_text(20.14, 19.43, 22.21, 'mm/yr') == '20.1 mm/yr (19.4 – 22.2, 95% CI)'
    assert ci_text(-12.3, -16.2, -8.9, 'mm/yr') == \
        '\u221212.3 mm/yr (\u221216.2 – \u22128.9, 95% CI)'
    assert ci_text(0.514, 0.371, 0.549, 'µstrain/yr') == \
        '0.51 µstrain/yr (0.37 – 0.55, 95% CI)'
    assert ci_text(29.87, 27.39, 32.95, 'deg') == '29.9° (27.4 – 33.0, 95% CI)'
    assert ci_text(3.0, 3.0, 3.0, '') == '3.0 (3.0 – 3.0, 95% CI)'


def test_geometry_helpers():
    assert azimuth_of(1, 0) == pytest.approx(90.0) and azimuth_of(0, -1) == pytest.approx(180.0)
    assert azimuth_of(-1, 1) == pytest.approx(315.0)
    east, north = local_en([130.1], [30.1], 130.0, 30.0)
    assert east[0] == pytest.approx(0.1 * np.cos(np.radians(30)) * 111.32, rel=0.01)
    assert north[0] == pytest.approx(11.13, rel=0.01)
    t = np.linspace(-5, 5, 11)
    assert strike_of_points(t * np.sin(np.radians(80)), t * np.cos(np.radians(80))) \
        == pytest.approx(80.0)
    # a profile running south: the strike points west, so that S > 0 is dextral
    assert orient_strike(80.0, 170.0) == pytest.approx(260.0)
    assert orient_strike(260.0, 170.0) == pytest.approx(260.0)
    assert slip_sense(5.0, 260.0, 170.0) == 'right-lateral'
    assert slip_sense(-5.0, 260.0, 170.0) == 'left-lateral'
    # a N-S fault, a W -> E profile: north is to the left of east
    assert slip_sense(5.0, orient_strike(0.0, 90.0), 90.0) == 'right-lateral'
    normal, angle = fault_normal(80.0, 170.0)
    assert normal == pytest.approx([np.cos(np.radians(80)), -np.sin(np.radians(80))])
    assert angle == pytest.approx(0.0, abs=1e-6)
    assert fault_normal(80.0, 125.0)[1] == pytest.approx(45.0)


##############################################################################
#                               the command                                  #
##############################################################################
def test_model_options(products):
    asc = products['asc']
    inps = cli.model_cmd_line_parse([asc, '--comp', 'los', '--noise', '3', '12', '--creep',
                                     '--free-x0', '--d1-prior', '17', '5', '--noverbose'])
    assert inps.noise == (3.0, 12.0) and inps.free_x0 == 5.0 and inps.creep
    assert inps.gnss_comp == 'par' and inps.bin_x == inps.bin_km == 1.0
    # --comp auto (the default) settles the GNSS component once the files are read
    auto = cli.model_cmd_line_parse([asc, '-gnss', products['gnss'], '--noverbose'])
    assert auto.comp == 'auto' and auto.gnss_comp is None
    with pytest.raises(ValueError, match='-vel'):
        cli.model_cmd_line_parse([asc, '--comp', 'as-is', '-gnss', products['gnss'],
                                  '--noverbose'])
    with pytest.raises(ValueError, match='--strike'):
        cli.model_cmd_line_parse([asc, '--comp', 'par', '-gnss', products['gnss'],
                                  '--noverbose'])
    with pytest.raises(ValueError, match='--noise'):
        cli.model_cmd_line_parse([asc, '--noise', 'pink', '--noverbose'])
    with pytest.raises(ValueError, match='--burn'):
        cli.model_cmd_line_parse([asc, '--steps', '100', '--burn', '200', '--noverbose'])
    with pytest.raises(ValueError, match='--nodisplay'):
        cli.model_cmd_line_parse([asc, '--nodisplay', '--noverbose'])
    with pytest.raises(ValueError, match='--geom'):
        cli.model_cmd_line_parse([asc, '--geom', 'nowhere.h5', '--noverbose'])
    # --cmap-prof: a colormap, viridis when left out - a file after it is an input
    assert cli.model_cmd_line_parse([asc, '--noverbose']).cmap_prof is None
    assert cli.model_cmd_line_parse([asc, '--cmap-prof', '--noverbose']).cmap_prof == 'viridis'
    assert cli.model_cmd_line_parse([asc, '--cmap-prof', 'magma_r']).cmap_prof == 'magma_r'
    bare = cli.model_cmd_line_parse(['--cmap-prof', asc, '--noverbose'])
    assert bare.cmap_prof == 'viridis' and bare.file == [asc]
    with pytest.raises(ValueError, match='--cmap-prof: no colormap'):
        cli.model_cmd_line_parse([asc, '--cmap-prof', 'no_such_map', '--noverbose'])
    # --creep-offset [KM [WIDTH]]: 5 and 2.5 km when left out - a file after it is an input
    assert cli.model_cmd_line_parse([asc, '--noverbose']).creep_offset is None
    assert cli.model_cmd_line_parse([asc, '--creep-offset']).creep_offset == (5.0, 2.5)
    assert cli.model_cmd_line_parse([asc, '--creep-offset', '4']).creep_offset == (4.0, 2.5)
    assert cli.model_cmd_line_parse([asc, '--creep-offset', '4', '1.5',
                                     '--noverbose']).creep_offset == (4.0, 1.5)
    bare = cli.model_cmd_line_parse(['--creep-offset', asc, '--comp', 'los', '--noverbose'])
    assert bare.creep_offset == (5.0, 2.5) and bare.file == [asc] and bare.comp == 'los'
    with pytest.raises(ValueError, match='KM and WIDTH at most'):
        cli.model_cmd_line_parse([asc, '--creep-offset', '4', '2', '1', '--noverbose'])
    with pytest.raises(ValueError, match='must be positive'):
        cli.model_cmd_line_parse([asc, '--creep-offset', '0', '--noverbose'])
    # the legend: inside, room made for it, else beside (auto); in, out, off
    plain = cli.model_cmd_line_parse([asc, '--noverbose'])
    assert plain.model_legend == 'auto'
    asked = cli.model_cmd_line_parse([asc, '--legend', 'out', '--noverbose'])
    assert asked.model_legend == 'out'
    assert cli.model_cmd_line_parse([asc, '--legend', '--noverbose']).model_legend == 'auto'
    bare = cli.model_cmd_line_parse(['--legend', asc, '--noverbose'])
    assert bare.model_legend == 'auto' and bare.file == [asc]
    with pytest.raises(SystemExit):
        cli.model_cmd_line_parse([asc, '--legend', 'above', '--noverbose'])
    # --comp vu: the vertical is a reverse or a normal fault's (viewinsar model
    # reverse|normal), not a strike-slip one's; one --comp for all the files
    with pytest.raises(ValueError, match='reverse or a normal'):
        cli.model_cmd_line_parse([asc, '--comp', 'vu', '-gnss', products['gnss'],
                                  '--noverbose'])
    with pytest.raises(ValueError, match='ONE --comp'):
        cli.model_cmd_line_parse([asc, products['desc'], '--comp', 'los', 'as-is',
                                  '--noverbose'])
    two = cli.model_cmd_line_parse([asc, products['desc'], '--comp', 'los', '--noverbose'])
    assert two.comps == ['los', 'los'] and two.fault_type == 'strike'
    # --bands OUTER INNER: 95 68 (1 sigma), or as asked
    assert cli.model_cmd_line_parse([asc, '--noverbose']).bands == [95.0, 68.0]
    assert cli.model_cmd_line_parse([asc, '--bands', '95', '65', '--noverbose']).bands == \
        [95.0, 65.0]
    with pytest.raises(ValueError, match='--bands'):
        cli.model_cmd_line_parse([asc, '--bands', '60', '70', '--noverbose'])


def test_model_end_to_end(products, tmp_path):
    out = tmp_path / 'out'
    obj = cli.model_main([products['asc'], '--comp', 'los', '-fault', products['fault'],
                          '-gnss', products['gnss'], *PROFILE, '--strain', *FAST,
                          '--nodisplay', '--noverbose', '--outdir', str(out)])
    sec = obj.sections[0]
    res = sec.result
    s = res['summary']
    assert s['S']['map'] == pytest.approx(S_TRUE, abs=2.0)
    assert s['d1']['map'] == pytest.approx(D1_TRUE, abs=4.0)
    assert s['ag']['map'] == pytest.approx(-5.0 * np.sin(np.radians(80)), abs=1.0)
    # the fault where the trace crosses, its strike, oriented: S > 0 = right-lateral
    assert sec.fault_how.startswith('where the profile crosses')
    assert sec.strike == pytest.approx(260.0, abs=1.0) and sec.strike_how.startswith('the fault')
    assert sec.sense_known and slip_sense(s['S']['map'], sec.strike, sec.prof_az) == \
        'right-lateral'
    assert sec.prof_az == pytest.approx(170.0, abs=2.0)

    base = os.path.basename(sec.outfile_base)
    assert base.startswith('model_A_')
    names = set(os.listdir(out))
    for tail in ('.png', '_posterior.png', '_params.csv', '_params.txt', '_samples.csv',
                 '_bins.csv', '_points.csv', '_gnss.csv'):
        assert f'{base}{tail}' in names, tail
    assert 'model_map_A.png' in names

    par = _csv(out / f'{base}_params.csv')
    assert list(par['name']) == ['S', 'd1', 'a', 'ag', 'strain']
    assert par['map'][0] == pytest.approx(s['S']['map'], rel=1e-5)
    assert np.all(par['lower_95'] <= par['upper_95'])
    text = (out / f'{base}_params.txt').read_text()
    for needle in ('right-lateral', 'Goodman & Weare', 'exponential covariance',
                   'the fault trace at the fault', 'where the profile crosses the fault trace'):
        assert needle in text, needle

    pts = _csv(out / f'{base}_points.csv')
    assert pts.dtype.names == ('dist_km', 'x_km', 'perp_km', 'value', 'v_model',
                               'v_shifted', 'row', 'col', 'lon', 'lat', 'file_idx')
    assert len(pts) == res['n_pixels']
    # the distance from the fault is the distance perpendicular to it
    assert pts['x_km'] == pytest.approx(_fault_frame(pts['lon'], pts['lat']), abs=0.15)
    # and the fault-parallel velocity is the LOS divided by its projection
    want = -_vpar_n80e(_fault_frame(pts['lon'], pts['lat']))      # along N260E
    assert np.median(pts['v_model'] - want) == pytest.approx(2.0 / 0.63, abs=0.6)

    gn = _csv(out / f'{base}_gnss.csv')
    assert len(gn) == 8 and 'S08' not in list(gn['site'])
    assert np.all(gn['fitted'] == 1)
    samples = _csv(out / f'{base}_samples.csv')
    assert samples.dtype.names == ('S', 'd1', 'a', 'ag', 'strain')

    # the window: red points, the medians, the MAP model in bold dashes and its
    # 95 % and 68 % ranges, the GNSS blue, the fault dotted, A and A', the
    # parameters in a box, the legend, the title; the strain panel
    _check_look(sec)
    ax = sec.ax_swath
    assert ax.get_title().startswith("A-A': swath ±15 km") and 'GNSS' in ax.get_title()
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels[:3] == ['InSAR', 'Median per 1 km (±1σ)', 'MAP model']
    for label in ('95% range', '68% range', 'GNSS', 'Fault'):
        assert label in labels, label
    _check_legend_clear(sec)
    gnss = [ln for ln in ax.get_lines() if ln.get_marker() == 'o' and ln.get_mec() == 'w']
    assert gnss and all(mcolors.to_hex(ln.get_mfc()) == '#1f77b4' for ln in gnss)
    assert sec.info.startswith("A-A': swath ±15 km") and 'GNSS' in sec.info
    assert ax.get_ylabel() == 'Fault-parallel velocity (mm/yr)'
    assert sec.ax_strain.get_xlabel() == 'Distance from the fault (km)'
    strain = [ln for ln in sec.ax_strain.get_lines() if ln.get_linestyle() == '--']
    assert len(strain) == 1
    assert sec.ax_strain is not None and sec.fig_post is not None
    _check_posterior(sec, ['S', 'd1', 'a'])


def _check_look(sec, density=False):
    """The profile as viewinsar draws it: A and A' in the top corners, the
    parameters in a box, the MAP model in bold black dashes in its 95 % and
    68 % ranges (25 % and 40 % transparent; red on a density), the medians, the
    fault dotted."""
    from matplotlib.collections import PolyCollection
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.text import Annotation

    ax = sec.ax_swath
    notes = {a.get_text(): a for a in ax.texts if isinstance(a, Annotation)}
    assert {'A', "A'"} <= set(notes) and notes["A'"].get_bbox_patch() is not None
    texts = [a for a in ax.artists if isinstance(a, AnchoredText)]
    assert len(texts) == 1
    text = texts[0].txt._text
    assert text.get_style() == 'normal' and '95% CI)' in text.get_text()
    assert text.get_text().startswith('slip rate')
    model = [ln for ln in ax.get_lines() if ln.get_linestyle() == '--' and ln.get_color() == 'k'
             and len(ln.get_xdata()) > 100]
    assert len(model) == 1                                      # the MAP: bold dashes
    assert [ln for ln in ax.get_lines() if ln.get_linestyle() == ':'
            and len(ln.get_xdata()) == 2]                       # the fault, dotted
    medians = [ln for ln in ax.get_lines() if ln.get_marker() == 'o' and ln.get_mec() != 'w']
    assert medians
    polys = [c for c in ax.collections if isinstance(c, PolyCollection)]
    shades = [mcolors.to_hex(c.get_facecolor()[0]) for c in polys]
    assert shades == (['#d62728', '#d62728'] if density else ['#9e9e9e', '#dbdbdb'])
    # 25 % and 40 % transparent (lighter on a density, so that it shows)
    assert [c.get_alpha() for c in polys] == ([0.4, 0.3] if density else [0.75, 0.6])


def _check_legend_clear(sec):
    """The legend hides neither the parameters nor the marks: inside the plot,
    or beside it on the right."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.offsetbox import AnchoredText

    ax = sec.ax_swath
    renderer = FigureCanvasAgg(ax.figure).get_renderer()     # the run closed its figures
    leg = ax.get_legend().get_window_extent(renderer)
    box = ax.get_window_extent(renderer)
    others = [a.get_window_extent(renderer) for a in ax.artists if isinstance(a, AnchoredText)]
    others += [a.get_window_extent(renderer) for a in ax.texts]
    for other in others:
        assert not leg.overlaps(other)
    inside = leg.x0 >= box.x0 - 1 and leg.x1 <= box.x1 + 1 and leg.y0 >= box.y0 - 1 \
        and leg.y1 <= box.y1 + 1
    assert inside or leg.x0 >= box.x1, (leg, box)
    return inside


def _check_posterior(sec, names):
    """The posterior window in the style of Hussain et al. (2016, Fig. 9)."""
    from matplotlib.collections import QuadMesh

    fig = sec.fig_post
    num = len(names)
    cbar = [ax for ax in fig.axes if 'frequency' in ax.get_xlabel()]
    joint = [ax for ax in fig.axes if any(isinstance(c, QuadMesh) for c in ax.collections)
             and ax not in cbar]
    hist = [ax for ax in fig.axes if ax.patches and not ax.collections]
    assert len(joint) == num * (num - 1) // 2             # every pair, jet
    assert all(ax.collections[0].get_cmap().name == 'jet' for ax in joint)
    assert all(ax.collections[0].get_clim() == (0, 1) for ax in joint)   # to its peak
    assert len(hist) == num                               # the marginals, blue
    assert [ax.get_xlabel().split(' [')[0] for ax in hist][:2] == ['Slip rate',
                                                                  'Locking depth']
    assert hist[0].get_ylabel() == 'Frequency'
    # the MAP: a red line on every marginal, a red dot on every pair
    for ax in hist:
        assert any(ln.get_color() == '#d62728' for ln in ax.get_lines())
    for ax in joint:
        assert any(ln.get_marker() == 'o' for ln in ax.get_lines())
    texts = [t.get_text() for t in fig.texts]
    assert f" {sec.label} - {sec.label}' " in texts      # white on black, top right
    assert len(cbar) == 1


def test_model_two_tracks(products, tmp_path):
    obj = cli.model_main([products['asc'], products['desc'], '--labels', 'T1A', 'T2D',
                          '--comp', 'los', '-fault', products['fault'], *PROFILE, *FAST,
                          '--nodisplay', '--noverbose', '--outdir', str(tmp_path)])
    res = obj.sections[0].result
    s = res['summary']
    assert res['names'] == ['S', 'd1', 'a', 'a2'] and res['sets'] == [0, 1]
    assert s['S']['map'] == pytest.approx(S_TRUE, abs=2.0)
    # each track keeps its own reference: +2 and -3 mm/yr in the LOS, divided by
    # the LOS . strike of each track (along N260E: asc +0.63, desc -0.56)
    asc = 2.0 / (-np.cos(np.radians(182)) * np.sin(np.radians(39)))
    desc = -3.0 / (-np.cos(np.radians(22)) * np.sin(np.radians(37)))
    assert s['a']['map'] - s['a2']['map'] == pytest.approx(asc - desc, abs=0.7)
    labels = [t.get_text() for t in obj.sections[0].ax_swath.get_legend().get_texts()]
    assert labels[:2] == ['T1A', 'T2D']
    bins = _csv(next(tmp_path.glob('*_bins.csv')))
    assert set(bins['file_idx']) == {0, 1}
    # the bins of both tracks, shifted, sit on the model
    assert np.median(np.abs(bins['residual'])) < 1.0


def test_model_creep_offset(tmp_path):
    """--creep --creep-offset on an ascending and a descending track with 6 mm/yr
    of creep down to 4 km and the A' (SSE) side rising 2 mm/yr: the offset of
    each track, biased each its own way by the vertical step; the two split into
    the fault-parallel creep and the vertical - as Hussain et al. (2016) found
    at Izmit; on the profile, in _params.csv / .txt and the box."""
    pytest.importorskip('h5py')
    tracks = dict(creep=(6.0, 4.0), up_step=2.0)
    asc = make_track(str(tmp_path / 'asc'), 39.0, 102.0, 2.0, seed=1, **tracks)
    desc = make_track(str(tmp_path / 'desc'), 37.0, -102.0, -3.0, lat0=30.62, lon0=130.10,
                      seed=2, heading=-168.0, **tracks)
    out = tmp_path / 'out'
    obj = cli.model_main([asc, desc, '--labels', 'asc', 'desc', '--comp', 'los', *PROFILE,
                          '--creep', '--creep-offset', '--legend', 'out', '--bands', '95',
                          '65', *FAST, '--nodisplay', '--noverbose', '--outdir', str(out)])
    sec = obj.sections[0]
    off = sec.offset
    assert off['km'] == 5.0 and off['width'] == 2.5
    f_asc, f_desc = off['files']
    assert f_asc['label'] == 'asc' and f_asc['n_a'] > 100 and f_asc['n_b'] > 100
    # the step up adds to the creep in the ascending LOS, takes from it in the descending
    assert f_asc['rate'] > f_desc['rate'] + 2.0
    for one in (f_asc, f_desc):
        assert one['lo'] < one['rate'] < one['hi']
    # split: the creep (less than 6 - the lines' 5 km span a creep that dies away
    # within 4 km) and the A' side up
    sp = off['split']
    assert sp is not None and sp['n'] == 2
    assert 2.5 < sp['par'] < 6.0 and sp['par_lo'] < sp['par'] < sp['par_hi']
    assert 0.5 < sp['up'] < 3.5 and sp['up_lo'] > 0.0

    base = os.path.basename(sec.outfile_base)
    par = _csv(out / f'{base}_params.csv')
    names = list(par['name'])
    for name in ('C', 'd2', 'C_offset_1', 'C_offset_2', 'C_offset_par', 'U_offset'):
        assert name in names, name
    row = names.index('C_offset_par')
    assert par['map'][row] == pytest.approx(sp['par'], rel=1e-5)
    assert par['lower_95'][row] == pytest.approx(sp['par_lo'], rel=1e-5)
    text = (out / f'{base}_params.txt').read_text()
    for needle in ('creep offset   : as Hussain et al. (2016)', 'split of the 2 LOS offsets',
                   "+ = the A' side up"):
        assert needle in text, needle

    # the profile: the box, the lines either side, the model without the creep
    ax = sec.ax_swath
    box = '\n'.join(sec.parameter_lines())
    for needle in ('creep rate = ', 'creep depth = ', 'creep offset, fault-parallel = ',
                   "A' side up"):
        assert needle in box, needle
    # --legend out: everything, beside the plot - the figure grown for it, the
    # plot as wide (in inches) as ever
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    for label in ('asc', 'desc', 'MAP model', '95% range', '65% range', 'MAP without creep',
                  'Offset lines, ±5 km', 'Fault'):
        assert label in labels, label
    assert not _check_legend_clear(sec)                       # beside, not inside
    (_, home), = [(a, box) for a, box in sec.home if a is ax]
    narrow, wide = sec._narrow[0], ax.figure.get_size_inches()[0]
    assert wide > narrow
    assert ax.get_position().width * wide == pytest.approx(home.width * narrow, rel=1e-6)
    assert sum(ln.get_color() == '#ff7f0e' for ln in ax.get_lines()) == 4
    obj.close()


def test_model_density_view(products, tmp_path):
    """--cmap-prof: the pixels as a point density, as swath draws it - no points,
    the 95 % range red, bigger medians, a colour bar; --nlook fills the density
    in but the fit takes the real pixels, exactly as without it."""
    from matplotlib.collections import PathCollection, PolyCollection, QuadMesh

    args = [products['asc'], products['desc'], '--labels', 'T1A', 'T2D', '--comp', 'los',
            '-fault', products['fault'], '-gnss', products['gnss'], *PROFILE, *FAST,
            '--nodisplay', '--noverbose']
    plain = cli.model_main([*args, '--outdir', str(tmp_path / 'plain')])
    dense = cli.model_main([*args, '--cmap-prof', 'magma', '--nlook', '3',
                            '--outdir', str(tmp_path / 'dense')])
    sp, sd = plain.sections[0], dense.sections[0]
    assert sd.result['n_pixels'] == sp.result['n_pixels']
    for name in ('S', 'd1', 'a', 'a2', 'ag'):
        for key in ('map', 'lo', 'hi'):
            assert sd.result['summary'][name][key] == pytest.approx(
                sp.result['summary'][name][key], rel=1e-9), (name, key)
    drawn = sum(f['x'].size for f in sd.fp_draw.values())
    assert drawn > 5 * sd.result['n_pixels'] and sp.fp_draw is None
    assert '(drawn 3x oversampled)' in sd.ax_swath.get_title()
    assert 'oversampled what is drawn only' in next(
        (tmp_path / 'dense').glob('*_params.txt')).read_text()
    pts = _csv(next((tmp_path / 'dense').glob('*_points.csv')))
    assert len(pts) == sd.result['n_pixels'] and np.all(pts['row'] == np.round(pts['row']))

    def bands(ax):
        return [mcolors.to_hex(c.get_facecolor()[0]) for c in ax.collections
                if isinstance(c, PolyCollection)]

    def median_size(ax):
        return [ln for ln in ax.get_lines() if ln.get_marker() == 'o'
                and ln.get_mfc() != 'k' and ln.get_mec() != 'w'][0].get_markersize()

    def gnss_face(ax):
        return {mcolors.to_hex(ln.get_mfc()) for ln in ax.get_lines()
                if ln.get_marker() == 'o' and ln.get_mec() == 'w'}

    ax = sd.ax_swath
    mesh = [c for c in ax.collections if isinstance(c, QuadMesh)]
    assert len(mesh) == 1 and mesh[0].get_cmap().name.startswith('magma')
    # every sample drawn, but for the few beyond the plot's limits
    assert 0.99 * drawn < mesh[0].get_array().sum() <= drawn
    assert not [c for c in ax.collections if isinstance(c, PathCollection)]
    assert bands(ax) == ['#d62728', '#d62728'] and bands(sp.ax_swath) == ['#9e9e9e', '#dbdbdb']
    _check_look(sd, density=True)
    _check_look(sp)
    # the medians: on a density twice the size
    assert median_size(ax) > 1.8 * median_size(sp.ax_swath)
    # the GNSS: red on a density, blue on the points - outlined in white either way
    assert gnss_face(ax) == {'#d62728'} and gnss_face(sp.ax_swath) == {'#1f77b4'}
    assert sd.ax_cbar.get_ylabel().startswith('Point density')
    # the legend: the files first, the GNSS red in it too
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels[:3] == ['T1A', 'T2D', 'Median per 1 km (±1σ)'] and 'GNSS' in labels
    _check_legend_clear(sd)
    # the default view: the points, no density
    assert [c for c in sp.ax_swath.collections if isinstance(c, PathCollection)]
    assert not [c for c in sp.ax_swath.collections if isinstance(c, QuadMesh)]
    plain.close()
    dense.close()


def test_model_several_profiles(products, tmp_path):
    """--profile and a two-line --line-file: A-A', B-B', C-C' fitted one by one,
    each with its own figures, and one table of all of them."""
    lines = tmp_path / 'lines.xy'
    lines.write_text(">  B-B'\n130.45 30.58\n130.57 30.02\n")
    out = tmp_path / 'out'
    obj = cli.model_main([products['asc'], '--comp', 'los', '-fault', products['fault'],
                          *PROFILE, '--line-file', str(lines),
                          '--profile', '30.58', '130.35', '30.02', '130.47', *FAST,
                          '--nodisplay', '--noverbose', '--outdir', str(out)])
    assert [s.label for s in obj.sections] == ['A', 'B', 'C']
    names = set(os.listdir(out))
    for sec in obj.sections:
        assert sec.has_data, sec.name
        assert sec.result['summary']['S']['map'] == pytest.approx(S_TRUE, abs=3.0), sec.name
        assert sec.fig_post is not None and sec.fig_post is not obj.sections[0].fig_post \
            or sec is obj.sections[0]
        _check_posterior(sec, ['S', 'd1', 'a'])
        base = os.path.basename(sec.outfile_base)
        assert base.startswith(f'model_{sec.label}_')
        for tail in ('.png', '_posterior.png', '_params.csv', '_params.txt'):
            assert f'{base}{tail}' in names, (sec.name, tail)
    assert 'model_map_ABC.png' in names
    table = _csv(out / 'model_sections_ABC.csv')
    assert list(table['section']) == ['A', 'B', 'C']
    assert list(table['sense']) == ['right-lateral'] * 3
    for k, sec in enumerate(obj.sections):
        assert table['S_map'][k] == pytest.approx(sec.result['summary']['S']['map'],
                                                  rel=1e-5)
    # where each profile crosses the fault trace (N80E through 130.45 E, 30.30 N)
    assert table['fault_lat'] == pytest.approx(
        FLAT + (table['fault_lon'] - FLON) * KM_LON / KM_LAT / np.tan(np.radians(STRIKE)),
        abs=0.01)
    assert np.all(table['S_lower_95'] <= table['S_upper_95'])
    obj.close()


def test_model_east_velocity_given_fault(products, tmp_path):
    """No trace: --x0 and --strike place the fault; the east velocity is divided
    by sin(strike)."""
    obj = cli.model_main([products['east'], '--comp', 've', '--x0', '31.45', '--strike', '80',
                          *PROFILE, *FAST, '--legend', 'off', '--nodisplay', '--noverbose',
                          '--outdir', str(tmp_path)])
    sec = obj.sections[0]
    assert sec.ax_swath.get_legend() is None                  # --legend off: none
    assert sec.fault_how == '--x0' and sec.strike_how == '--strike'
    assert sec.strike == pytest.approx(260.0)            # turned: S > 0 is right-lateral
    s = sec.result['summary']
    assert s['S']['map'] == pytest.approx(S_TRUE, abs=2.0)
    assert s['a']['map'] == pytest.approx(-4.0 / np.sin(np.radians(80)), abs=0.8)


def test_model_as_is_needs_no_strike(products, tmp_path):
    obj = cli.model_main([products['asc'], '--comp', 'as-is', *PROFILE, *FAST, '--nodisplay',
                          '--noverbose', '--outdir', str(tmp_path), '--noise', 'white'])
    sec = obj.sections[0]
    assert not sec.sense_known and sec.strike_how == 'perpendicular to the profile'
    assert sec.fault_how.startswith('the middle')
    # the LOS as it is: the slip rate seen in the LOS (x ~0.63)
    assert abs(sec.result['summary']['S']['map']) == pytest.approx(S_TRUE * 0.63, abs=2.0)
    # what is plotted: the LOS, as it is - and the rates are those it sees
    assert sec.ax_swath.get_ylabel() == 'LOS velocity (mm/yr)'
    assert sec.parameter_lines()[0].startswith('slip rate in the LOS = ')


def test_model_names_the_data(products, tmp_path):
    """The y axis says what is plotted: fault-parallel velocity once the data
    are turned into it; a component as it is, by the name of its file; the
    LOS as it is. A file named as the vertical is a reverse or a normal
    fault's (test_dipslip)."""
    import shutil

    from viewinsar.modelview import VEL_LABEL, ModelViewer

    up = tmp_path / 'geo_up.h5'
    shutil.copy(products['east'], up)
    cases = ((products['east'], 'auto', 'as-is', 'east'), (products['east'], 've', 've', 'par'),
             (products['asc'], 'as-is', 'as-is', 'los'), (products['asc'], 'auto', 'los', 'par'))
    for fname, comp, resolved, plotted in cases:
        inps = cli.model_cmd_line_parse([fname, '--comp', comp, '--noverbose'])
        inps.disp_fig = False
        obj = ModelViewer(inps, cli.get_view_cmd(inps))
        obj.open()
        assert (obj.comp, obj.plotted) == (resolved, plotted), (fname, comp)
        assert obj.comps == [resolved]
        obj.close()
    inps = cli.model_cmd_line_parse([str(up), '--noverbose'])
    inps.disp_fig = False
    obj = ModelViewer(inps, cli.get_view_cmd(inps))
    with pytest.raises(ValueError, match='reverse or a normal'):
        obj.open()
    obj.close()
    assert VEL_LABEL['z'] == 'Fault vertical velocity' and VEL_LABEL['par'] == \
        'Fault-parallel velocity' and VEL_LABEL['los'] == 'LOS velocity'
    assert VEL_LABEL['h'] == 'Fault-normal velocity'


def test_model_interactive(products, tmp_path):
    from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent

    from viewinsar.modelview import ModelViewer

    inps = cli.model_cmd_line_parse([products['asc'], '--comp', 'los', '-fault',
                                     products['fault'], '-w', '15', *FAST, '--noverbose',
                                     '--outdir', str(tmp_path)])
    inps.disp_fig = False
    obj = ModelViewer(inps, cli.get_view_cmd(inps))
    obj.open()
    obj.plot()
    obj.fig_map.canvas.draw()

    def click(lon, lat):
        x, y = obj.ax_img.transData.transform((lon, lat))
        obj.select_point(MouseEvent('button_release_event', obj.fig_map.canvas, x, y,
                                    button=MouseButton.LEFT))

    def key(name):
        obj.on_key(KeyEvent('key_press_event', obj.fig_map.canvas, key=name))

    click(130.40, 30.58)
    click(130.52, 30.02)
    sec = obj.sections[0]
    assert sec.has_data and sec.result['summary']['S']['map'] == pytest.approx(S_TRUE, abs=2.5)
    assert len(sec.map_marks) == 1                       # the fault on the map
    # "l" on the profile: the legend in, out, off, in turn - drawn again, not
    # fitted again, and not the log y axis of matplotlib's own "l"
    result = sec.result
    sec.fig.canvas.draw()
    seen = [sec.legend_at]
    for _ in range(3):
        # the mouse over the plot: where matplotlib's "l" would act
        x, y = sec.ax_swath.transAxes.transform((0.5, 0.5))
        sec.fig.canvas.callbacks.process('key_press_event', KeyEvent(
            'key_press_event', sec.fig.canvas, key='l', x=x, y=y))
        seen.append(sec.legend_at)
    order = ['in', 'out', 'off']
    assert seen[1:] == [order[(order.index(seen[0]) + k) % 3] for k in (1, 2, 3)]
    assert sec.result is result and sec.ax_swath.get_yscale() == 'linear'
    assert (sec.ax_swath.get_legend() is None) == (sec.legend_at == 'off')
    key('e')
    assert any(p.name.endswith('_params.csv') for p in tmp_path.iterdir())
    key('=')
    assert len(obj.sections) == 2 and not obj.sections[1].has_data
    key('n')
    assert not obj.sections[1].has_data
    obj.close()


def test_model_on_radar_coded_data(tmp_path):
    """A radar grid: the lon / lat of every pixel come from geometryRadar.h5."""
    from conftest import FAULT_LON
    from test_radar import SECTION, make_radar_product

    vfile, _ = make_radar_product(str(tmp_path / 'proc'))
    fault = tmp_path / 'fault.geojson'
    fault.write_text(json.dumps({'type': 'FeatureCollection', 'features': [{
        'type': 'Feature', 'properties': {}, 'geometry': {
            'type': 'LineString', 'coordinates': [[FAULT_LON, 30.0], [FAULT_LON, 30.8]]}}]}))
    obj = cli.model_main([vfile, *SECTION, '--comp', 'as-is', '-w', '5', '-fault', str(fault),
                          '--s-max', '1500',
                          *FAST, '--nodisplay', '--noverbose', '--outdir', str(tmp_path / 'out')])
    sec = obj.sections[0]
    assert sec.fault_how.startswith('where the profile crosses')
    assert sec.strike % 180 == pytest.approx(0.0, abs=2.0) or \
        sec.strike % 180 == pytest.approx(180.0, abs=2.0)
    # the fixture: 0.30 m/yr x 2/pi atan(x / 8 km) -> S = 600 mm/yr, d1 = 8 km
    s = sec.result['summary']
    assert abs(s['S']['map']) == pytest.approx(600.0, rel=0.1)
    assert s['d1']['map'] == pytest.approx(8.0, rel=0.25)


def test_help_mentions_model(capsys):
    cli.main(['--help'])
    out = capsys.readouterr().out
    assert 'model' in out and 'locking depth' in out
    with pytest.raises(SystemExit):
        cli.main(['model', '--help'])
    out = capsys.readouterr().out
    for needle in ('--comp', '--creep', '--free-x0', '--d1-prior', '--noise', '--walkers',
                   '-gnss', '_posterior.png', '_samples.csv', 'Goodman & Weare',
                   'Hussain et al.', '--cmap-prof', 'point DENSITY', '--creep-offset',
                   'bootstrap', '--legend', '--bands'):
        assert needle in out, needle

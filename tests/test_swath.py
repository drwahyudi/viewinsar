"""Tests of the numerical core - no MintPy, no matplotlib needed."""

import numpy as np
import pytest

from conftest import DLAT, DLON, LAT0, LENGTH, LON0, WIDTH, make_velocity
from viewinsar.swath import (
    EARTH_RADIUS,
    get_pixel_size,
    line_geometry,
    nice_step,
    swath_density,
    swath_extract,
    weighted_quantile,
)

KM = 1e-3


##############################################################################
#                            pixel size                                      #
##############################################################################
def test_pixel_size_geographic(velocity_attrs):
    lat_c = 30.3
    x_step, y_step, unit = get_pixel_size(velocity_attrs, lat_c=lat_c)
    assert unit == 'm'
    assert x_step == pytest.approx(
        abs(DLON) * np.pi / 180 * EARTH_RADIUS * np.cos(np.deg2rad(lat_c)), rel=1e-12)
    assert y_step == pytest.approx(abs(DLAT) * np.pi / 180 * EARTH_RADIUS, rel=1e-12)


def test_pixel_size_utm():
    atr = {'Y_FIRST': '4e6', 'X_FIRST': '5e5', 'Y_STEP': '-30', 'X_STEP': '30',
           'UTM_ZONE': '11N', 'LENGTH': '10', 'WIDTH': '10'}
    assert get_pixel_size(atr) == (30.0, 30.0, 'm')


def test_pixel_size_radar():
    atr = {'LENGTH': '10', 'WIDTH': '10', 'PROCESSOR': 'isce',
           'RANGE_PIXEL_SIZE': '2.33', 'AZIMUTH_PIXEL_SIZE': '14.1',
           'HEIGHT': '693000.0', 'EARTH_RADIUS': '6371000.0',
           'INCIDENCE_ANGLE': '39.0'}
    x_step, y_step, unit = get_pixel_size(atr)
    assert unit == 'm'
    assert x_step == pytest.approx(2.33 / np.sin(np.deg2rad(39.0)), rel=1e-9)
    assert y_step == pytest.approx(14.1 * 6371000.0 / (6371000.0 + 693000.0), rel=1e-9)


def test_pixel_size_unknown_falls_back_to_pixel():
    assert get_pixel_size({'LENGTH': '10', 'WIDTH': '10'}) == (1.0, 1.0, 'pixel')


##############################################################################
#                          line geometry                                     #
##############################################################################
def test_line_geometry_length_and_frame(velocity_attrs):
    start, end = [150, 75], [150, 375]           # due east, 300 columns
    geom = line_geometry(velocity_attrs, start, end, half_width=10.)
    x_step = geom['x_step']
    assert geom['length'] == pytest.approx(300 * x_step, rel=1e-12)
    assert geom['uvec'] == pytest.approx([1.0, 0.0], abs=1e-12)
    # +Y is north, so the left normal of an eastward line points north
    assert geom['nvec'] == pytest.approx([0.0, 1.0], abs=1e-12)
    assert geom['corners_yx'].shape == (4, 2)


def test_line_geometry_rejects_zero_length(velocity_attrs):
    with pytest.raises(ValueError, match='identical'):
        line_geometry(velocity_attrs, [10, 10], [10, 10], half_width=5.)


##############################################################################
#                    buffer + linear referencing                             #
##############################################################################
def _brute_force(atr, start, end, row, col):
    """Along/off-profile distance of one pixel, computed from scratch, in km."""
    la0, lo0 = LAT0 + (start[0] + .5) * DLAT, LON0 + (start[1] + .5) * DLON
    la1, lo1 = LAT0 + (end[0] + .5) * DLAT, LON0 + (end[1] + .5) * DLON
    lap, lop = LAT0 + (row + .5) * DLAT, LON0 + (col + .5) * DLON
    ky = EARTH_RADIUS * np.pi / 180 * KM
    kx = ky * np.cos(np.deg2rad((la0 + la1) / 2))
    pnt = np.array([(lop - lo0) * kx, (lap - la0) * ky])
    uvec = np.array([(lo1 - lo0) * kx, (la1 - la0) * ky])
    uvec = uvec / np.hypot(*uvec)
    nvec = np.array([-uvec[1], uvec[0]])
    return float(pnt @ uvec), float(pnt @ nvec)


@pytest.mark.parametrize('start,end', [
    ([150, 60], [150, 400]),      # due east
    ([40, 225], [260, 225]),      # due south
    ([60, 60], [250, 400]),       # oblique
    ([250, 400], [60, 60]),       # the same line, reversed
])
def test_linear_referencing_matches_brute_force(velocity_attrs, start, end):
    data = make_velocity()
    box = (0, 0, WIDTH, LENGTH)
    sw = swath_extract(data, velocity_attrs, box, start, end, half_width=8.)

    assert sw['value'].size > 1000
    rng = np.random.default_rng(0)
    for i in rng.integers(0, sw['value'].size, 25):
        dist, perp = _brute_force(velocity_attrs, start, end, sw['row'][i], sw['col'][i])
        assert dist == pytest.approx(sw['dist'][i], abs=2e-3)
        assert perp == pytest.approx(sw['perp'][i], abs=2e-3)
        # and the value really is the pixel it claims to be
        assert sw['value'][i] == pytest.approx(data[sw['row'][i], sw['col'][i]], rel=1e-6)


def test_flat_cap_bounds(velocity_attrs):
    data = make_velocity()
    sw = swath_extract(data, velocity_attrs, (0, 0, WIDTH, LENGTH),
                       [60, 60], [250, 400], half_width=7.5)
    assert np.abs(sw['perp']).max() <= 7.5 + 1e-6      # inside the buffer
    assert sw['dist'].min() >= -1e-6                   # flat cap at the start
    assert sw['dist'].max() <= sw['length'] + 1e-6     # flat cap at the end
    # a rectangle: the corridor area is recovered by counting the pixels
    x_step, y_step = get_pixel_size(velocity_attrs, lat_c=30.3)[:2]
    area_pix = sw['value'].size * x_step * y_step * KM * KM
    assert area_pix == pytest.approx(2 * 7.5 * sw['length'], rel=0.10)


def test_lat_lon_columns_match_the_row_col(velocity_attrs):
    sw = swath_extract(make_velocity(), velocity_attrs, (0, 0, WIDTH, LENGTH),
                       [150, 75], [150, 375], half_width=5.)
    assert np.allclose(sw['lat'], LAT0 + (sw['row'] + 0.5) * DLAT)
    assert np.allclose(sw['lon'], LON0 + (sw['col'] + 0.5) * DLON)


def test_nan_and_zero_handling(velocity_attrs):
    data = make_velocity()
    data[150, 100:110] = 0.0
    kwargs = dict(atr=velocity_attrs, box=(0, 0, WIDTH, LENGTH),
                  start_yx=[150, 75], end_yx=[150, 375], half_width=5.)
    drop = swath_extract(data, keep_zero=False, **kwargs)
    keep = swath_extract(data, keep_zero=True, **kwargs)
    assert keep['value'].size - drop['value'].size == 10
    assert np.isfinite(drop['value']).all()             # the NaN strip never appears
    assert not np.any(drop['value'] == 0)


def test_mask_is_applied(velocity_attrs):
    data = make_velocity()
    mask = np.ones(data.shape, dtype=bool)
    mask[:, 200:220] = False
    kwargs = dict(atr=velocity_attrs, box=(0, 0, WIDTH, LENGTH),
                  start_yx=[150, 75], end_yx=[150, 375], half_width=5.)
    full = swath_extract(data, **kwargs)
    part = swath_extract(data, mask=mask, **kwargs)
    assert part['value'].size < full['value'].size
    assert not np.any((part['col'] >= 200) & (part['col'] < 220))


def test_decimate(velocity_attrs):
    kwargs = dict(atr=velocity_attrs, box=(0, 0, WIDTH, LENGTH),
                  start_yx=[150, 75], end_yx=[150, 375], half_width=5.)
    full = swath_extract(make_velocity(), **kwargs)
    dec = swath_extract(make_velocity(), decimate=3, **kwargs)
    assert dec['value'].size == pytest.approx(full['value'].size / 9, rel=0.15)
    assert np.abs(dec['perp']).max() <= 5.0 + 1e-6


def test_clipping_flag_and_offset_box(velocity_attrs):
    """A read window that does not cover the corridor must be flagged."""
    data = make_velocity()
    full = swath_extract(data, velocity_attrs, (0, 0, WIDTH, LENGTH),
                         [150, 75], [150, 375], half_width=5.)
    assert full['clipped'] is False

    # a window that cuts into the +/-5 km corridor (~23 rows) on all sides
    box = (100, 135, 350, 165)
    sub = data[box[1]:box[3], box[0]:box[2]]
    part = swath_extract(sub, velocity_attrs, box, [150, 75], [150, 375], half_width=5.)
    assert part['clipped'] is True
    assert part['value'].size < full['value'].size
    # the row/col are still full-file indices, and still point at the right values
    assert part['row'].min() >= box[1] and part['col'].min() >= box[0]
    assert np.allclose(part['value'], data[part['row'], part['col']], equal_nan=True)


def test_corridor_entirely_outside_raises(velocity_attrs):
    data = make_velocity()[:20, :20]
    with pytest.raises(ValueError, match='outside'):
        swath_extract(data, velocity_attrs, (400, 250, 420, 270),
                      [10, 10], [20, 200], half_width=1.)


##############################################################################
#                       density grid + statistics                            #
##############################################################################
def test_zero_is_always_a_cell_boundary():
    rng = np.random.default_rng(1)
    dist = rng.uniform(0, 57.3, 20000)
    val = rng.normal(0, 5, 20000)
    for offset in (0.0, 27.5, -3.3, 12.25):
        den = swath_density(dist, val, bin_x=0.5, bin_y=1.5, ylim=[-40, 40], x_offset=offset)
        assert np.min(np.abs(den['x_edges'])) == pytest.approx(0.0, abs=1e-12)


def test_density_counts_and_coverage():
    rng = np.random.default_rng(2)
    dist = rng.uniform(0, 30, 50000)
    val = rng.normal(0, 4, 50000)
    den = swath_density(dist, val, bin_x=0.5, bin_y=1.0, ylim=[-40, 40])
    assert den['counts'].sum() + den['num_outside'] == 50000
    assert den['count'].sum() == 50000            # the exact stats keep every point
    assert den['x_edges'][0] <= dist.min() and den['x_edges'][-1] >= dist.max()


def test_density_recovers_a_known_step():
    """A clean step at 20 km must show up in the median within one bin."""
    rng = np.random.default_rng(3)
    dist = rng.uniform(0, 40, 200000)
    val = np.where(dist < 20, -10.0, 10.0) + rng.normal(0, 1.0, dist.size)
    den = swath_density(dist, val, bin_x=0.5, bin_y=0.5, ylim=[-30, 30])
    left = den['median_dw'][den['x_centers'] < 19]
    right = den['median_dw'][den['x_centers'] > 21]
    assert np.nanmedian(left) == pytest.approx(-10.0, abs=0.3)
    assert np.nanmedian(right) == pytest.approx(10.0, abs=0.3)
    # the exact and the density-weighted medians agree away from the step
    assert np.nanmax(np.abs(den['median'] - den['median_dw'])[den['x_centers'] < 19]) < 0.6


def test_x_offset_shifts_the_axis():
    rng = np.random.default_rng(4)
    dist = rng.uniform(0, 40, 20000)
    val = rng.normal(0, 3, 20000)
    plain = swath_density(dist, val, bin_x=0.5, bin_y=1.0, ylim=[-20, 20])
    moved = swath_density(dist, val, bin_x=0.5, bin_y=1.0, ylim=[-20, 20], x_offset=20.)
    assert moved['x_centers'].min() == pytest.approx(plain['x_centers'].min() - 20, abs=0.5)
    assert moved['counts'].sum() == plain['counts'].sum()


def test_ylim_clipping_is_reported():
    dist = np.linspace(0, 10, 1000)
    val = np.linspace(-100, 100, 1000)
    den = swath_density(dist, val, bin_x=1.0, bin_y=1.0, ylim=[-10, 10])
    assert den['num_outside'] == np.sum((val < den['y_edges'][0]) | (val > den['y_edges'][-1]))
    assert den['num_outside'] > 0
    assert den['count'].sum() == 1000            # exact stats are NOT clipped


def test_auto_bin_y_and_ylim():
    rng = np.random.default_rng(5)
    val = rng.normal(0, 5, 30000)
    den = swath_density(rng.uniform(0, 10, 30000), val, bin_x=0.5)
    assert 30 <= den['y_edges'].size <= 90
    assert den['ylim'][0] < np.median(val) < den['ylim'][1]


def test_weighted_quantile():
    values = np.array([1., 2., 3., 4.])
    assert weighted_quantile(values, np.array([1., 1., 1., 1.]), 0.5) == 2.
    assert weighted_quantile(values, np.array([0., 0., 10., 0.]), 0.5) == 3.
    assert np.isnan(weighted_quantile(values, np.zeros(4), 0.5))
    # against numpy on a flat weighting
    rng = np.random.default_rng(6)
    vals = rng.normal(size=999)
    for q in (0.25, 0.5, 0.75):
        assert weighted_quantile(vals, np.ones(999), q) == pytest.approx(
            np.quantile(vals, q, method='inverted_cdf'), abs=1e-9)


@pytest.mark.parametrize('vrange,expect', [(80, 1.5), (8, 0.15), (800, 15.0), (0, 1.0)])
def test_nice_step(vrange, expect):
    assert nice_step(vrange) == pytest.approx(expect, rel=1e-9)


def test_density_needs_points():
    with pytest.raises(ValueError, match='no valid pixel'):
        swath_density(np.array([]), np.array([]), bin_x=0.5)

"""Elastic dislocation models of interseismic velocity profiles - the core of
`viewinsar model`.

Strike-slip faults. The fault-parallel velocity at distance x along a
fault-normal profile (Savage & Burford 1973), with optional shallow creep
from the surface down to d2 and a rotation term (Hussain et al. 2016, 2018):

    v(x) = S/pi arctan(u / d1) + a
         + C [H(u) - 1/2 - 1/pi arctan(u / d2)]       creep, d2 <= d1
         + rot u   on one side of the fault           rotation (fixed)
    u = x - x0

S is the slip rate below the locking depth d1, C the creep rate on the fault
above d2 - in the same sense as S, so 0 <= C - a the offset and x0 the fault
position. Every input file (frame, track) has an offset of its own, a, a2,
a3 ..., and the GNSS one more, ag: each is referenced on its own.

Reverse and normal faults. A fault dipping at dip, locked from the surface to
d1 and slipping steadily at S (along its dip) below: the 2D edge dislocation
of Freund & Barnett (1976), the back-slip model of Savage (1983) - Segall
(2010, eqs. 3.71-3.72). With u the distance from the trace towards the side
the fault dips to (the hanging wall) and z = u / d1 - cot(dip), at the
surface:

    horizontal (+ towards the dip)  -S/pi [cos(dip) atan z + (sin(dip) - z cos(dip)) / (1 + z^2)]
    vertical   (+ up)               +S/pi [sin(dip) atan z + (cos(dip) + z sin(dip)) / (1 + z^2)]

for a reverse fault (the hanging wall moving up the dip); a normal fault is
the opposite. Far from the fault the hanging wall moves S cos(dip) towards
the footwall and S sin(dip) up, relative to it: the shortening and uplift
rates (extension and subsidence for a normal fault). Checked against Okada
(1992) for a long, deep fault.

Each file's swath is binned along the profile (the median of every bin, as
Ou et al. 2022 do before fitting). The bins are fitted under uniform priors
with the affine-invariant ensemble Markov chain Monte Carlo sampler of
Goodman & Weare (2010) - the sampler of emcee, used by Hussain et al. - and a
noise covariance that decays exponentially with distance, sill exp(-r / L),
estimated from the residuals of a least-squares fit (or given). The result:
the maximum a posteriori (MAP) model and the 95 % range of every parameter
and of the model curve.

Numpy and scipy only.
"""

import warnings

import numpy as np

__all__ = ['ProfileModel', 'DipSlipModel', 'forward', 'dip_slip', 'bin_profile',
           'exponential_covariance',
           'fit_noise_covariance', 'ensemble_sample', 'autocorr_time', 'fit_profile',
           'summarize_samples', 'model_band', 'model_bands', 'strain_rate', 'fault_offset',
           'fit_two_lines', 'two_line_matrix', 'split_offsets',
           'describe', 'NAMES', 'UNITS']

NAMES = {'S': 'slip rate', 'd1': 'locking depth', 'C': 'creep rate', 'd2': 'creep depth',
         'x0': 'fault position', 'a': 'offset', 'ag': 'GNSS offset', 'dip': 'dip'}
UNITS = {'S': 'mm/yr', 'd1': 'km', 'C': 'mm/yr', 'd2': 'km', 'x0': 'km', 'a': 'mm/yr',
         'ag': 'mm/yr', 'dip': 'deg'}
FAULTS = ('strike', 'reverse', 'normal')


def describe(name):
    """'slip rate', 'offset of file 2' ..."""
    if name in NAMES:
        return NAMES[name]
    if name.startswith('a') and name[1:].isdigit():
        return f'offset of file {name[1:]}'
    return name


##############################################################################
#                               the model                                    #
##############################################################################
class ProfileModel:
    """Which parameters are free, their bounds and priors.

    n_sets   - number of InSAR data sets (files), each with its own offset
    creep    - add C and d2 (shallow creep above d2 <= d1)
    free_x0  - (lo, hi) in km: fit the fault position within that range;
               None: the fault is at x = 0
    rot      - mm/yr per km, a linear term on one side of the fault (fixed)
    rot_side - +1: on the x > 0 side, -1: on the x < 0 side
    gnss     - add the offset of the GNSS, ag
    bounds   - {name: (lo, hi)}; offsets default to the data's range (set by
               fit_profile)
    d1_prior - (mean, sigma): a Gaussian prior on d1 on top of its bounds
    """

    def __init__(self, n_sets=1, creep=False, free_x0=None, rot=0.0, rot_side=1, gnss=False,
                 bounds=None, d1_prior=None):
        self.n_sets, self.creep, self.free_x0 = int(n_sets), bool(creep), free_x0
        self.rot, self.rot_side, self.gnss = float(rot or 0.0), 1 if rot_side >= 0 else -1, gnss
        self.offsets = ['a'] + [f'a{k + 1}' for k in range(1, self.n_sets)]
        self.names = ['S', 'd1']
        if self.creep:
            self.names += ['C', 'd2']
        if free_x0 is not None:
            self.names.append('x0')
        self.names += self.offsets
        if gnss:
            self.names.append('ag')
        self.bounds = {'S': (-60.0, 60.0), 'd1': (0.1, 60.0), 'C': (0.0, 30.0),
                       'd2': (0.01, 40.0)}
        if free_x0 is not None:
            self.bounds['x0'] = (float(free_x0[0]), float(free_x0[1]))
        for name in self.offsets + (['ag'] if gnss else []):
            self.bounds[name] = (-1e3, 1e3)
        self.bounds.update({k: tuple(map(float, v)) for k, v in (bounds or {}).items()
                            if k in self.names})
        self.d1_prior = d1_prior
        self.fault = 'strike'
        self.kinds, self.gnss_kind = [None] * self.n_sets, None   # one component: no kinds

    @property
    def ndim(self):
        return len(self.names)

    def index(self, name):
        return self.names.index(name)

    def lo_hi(self):
        lo = np.array([self.bounds[n][0] for n in self.names], dtype=np.float64)
        hi = np.array([self.bounds[n][1] for n in self.names], dtype=np.float64)
        return lo, hi

    def unpack(self, theta):
        """{name: array} from (..., ndim) parameters; x0 = 0 when it is fixed."""
        theta = np.asarray(theta, dtype=np.float64)
        out = {n: theta[..., i] for i, n in enumerate(self.names)}
        out.setdefault('x0', np.zeros(theta.shape[:-1]))
        return out

    def predict(self, theta, x, offset='a', kind=None):
        """The model at x (km): (x.size,) for one parameter vector, (K, x.size)
        for K of them (kind: for the models of dip-slip faults)."""
        par = self.unpack(theta)
        return forward(x, par['S'], par['d1'], par[offset], x0=par['x0'],
                       creep=par['C'] if self.creep else None,
                       d2=par['d2'] if self.creep else None,
                       rot=self.rot, rot_side=self.rot_side)

    def log_prior(self, theta):
        theta = np.atleast_2d(theta)
        lo, hi = self.lo_hi()
        ok = np.all((theta >= lo) & (theta <= hi), axis=1)
        if self.creep:
            ok &= theta[:, self.index('d2')] <= theta[:, self.index('d1')]
        out = np.where(ok, 0.0, -np.inf)
        if self.d1_prior is not None:
            mean, sig = self.d1_prior
            out = out - 0.5 * ((theta[:, self.index('d1')] - mean) / sig) ** 2
        return out


class DipSlipModel:
    """The model of a reverse or a normal fault: its parameters, bounds, priors.

    sense    - +1: reverse (the hanging wall moves up the dip), -1: normal
    kinds    - per data set (file): 'h' (the horizontal velocity across the
               fault, + towards the end of the profile), 'z' (vertical, + up)
               or (ph, pu): a line of sight, ph h + pu z
    dip      - (lo, hi): the bounds of the dip (deg), or one value: fixed
    toward   - +1: the fault dips towards x > 0 (the end of the profile), -1
    free_x0  - (lo, hi) in km: fit the position of the trace; None: at 0
    gnss     - add the offset of the GNSS, ag; gnss_kind: their kind
    bounds, d1_prior - as ProfileModel; S >= 0, its sense given by `sense`
    """

    def __init__(self, sense=1, kinds=('z',), dip=(5.0, 85.0), toward=1, free_x0=None,
                 gnss=False, gnss_kind='z', bounds=None, d1_prior=None):
        self.sense = 1 if sense >= 0 else -1
        self.fault = 'reverse' if self.sense > 0 else 'normal'
        self.toward = 1 if toward >= 0 else -1
        self.kinds = list(kinds)
        self.n_sets = len(self.kinds)
        self.free_x0, self.gnss, self.gnss_kind = free_x0, gnss, gnss_kind
        self.creep, self.rot, self.rot_side = False, 0.0, 1    # strike-slip only
        dip = np.atleast_1d(np.asarray(dip, dtype=np.float64))
        self.dip_fixed = float(dip[0]) if dip.size == 1 or dip[0] == dip[-1] else None
        self.offsets = ['a'] + [f'a{k + 1}' for k in range(1, self.n_sets)]
        self.names = ['S', 'd1'] + (['dip'] if self.dip_fixed is None else [])
        if free_x0 is not None:
            self.names.append('x0')
        self.names += self.offsets
        if gnss:
            self.names.append('ag')
        self.bounds = {'S': (0.0, 60.0), 'd1': (0.1, 60.0)}
        if self.dip_fixed is None:
            self.bounds['dip'] = (float(dip[0]), float(dip[-1]))
        if free_x0 is not None:
            self.bounds['x0'] = (float(free_x0[0]), float(free_x0[1]))
        for name in self.offsets + (['ag'] if gnss else []):
            self.bounds[name] = (-1e3, 1e3)
        self.bounds.update({k: tuple(map(float, v)) for k, v in (bounds or {}).items()
                            if k in self.names})
        self.d1_prior = d1_prior

    @property
    def ndim(self):
        return len(self.names)

    def index(self, name):
        return self.names.index(name)

    def lo_hi(self):
        lo = np.array([self.bounds[n][0] for n in self.names], dtype=np.float64)
        hi = np.array([self.bounds[n][1] for n in self.names], dtype=np.float64)
        return lo, hi

    def unpack(self, theta):
        """{name: array}; x0 = 0 when it is fixed, dip the fixed one."""
        theta = np.asarray(theta, dtype=np.float64)
        out = {n: theta[..., i] for i, n in enumerate(self.names)}
        out.setdefault('x0', np.zeros(theta.shape[:-1]))
        if self.dip_fixed is not None:
            out['dip'] = np.full(theta.shape[:-1], self.dip_fixed)
        return out

    def predict(self, theta, x, offset='a', kind=None):
        """The model at x (km), of a kind ('h', 'z', (ph, pu)) - by default
        that of the first data set."""
        par = self.unpack(theta)
        return dip_slip(x, par['S'], par['d1'], par['dip'], par[offset], x0=par['x0'],
                        kind=self.kinds[0] if kind is None else kind, toward=self.toward,
                        sense=self.sense)

    def log_prior(self, theta):
        theta = np.atleast_2d(theta)
        lo, hi = self.lo_hi()
        ok = np.all((theta >= lo) & (theta <= hi), axis=1)
        out = np.where(ok, 0.0, -np.inf)
        if self.d1_prior is not None:
            mean, sig = self.d1_prior
            out = out - 0.5 * ((theta[:, self.index('d1')] - mean) / sig) ** 2
        return out

    def rates(self, theta):
        """The far-field rates of the hanging wall relative to the footwall:
        horizontal S cos(dip) (shortening, or extension) and vertical S
        sin(dip) (uplift, or subsidence)."""
        par = self.unpack(theta)
        dip = np.radians(par['dip'])
        return par['S'] * np.cos(dip), par['S'] * np.sin(dip)


def dip_slip(x, S, d1, dip, a=0.0, x0=0.0, kind='z', toward=1, sense=1):
    """The surface velocity across a dip-slip fault locked from the surface to
    d1, slipping at S along its dip below (see the module's docstring).

    kind   - 'h': horizontal, across the fault (+ towards x > 0); 'z': up;
             (ph, pu): the line of sight that sees ph of h and pu of z
    toward - +1: the fault dips towards x > 0, -1: towards x < 0
    sense  - +1: reverse, -1: normal
    Scalar parameters give an array shaped like x; arrays of K parameters give
    (K, x.size).
    """
    x = np.asarray(x, dtype=np.float64)
    pars = [np.asarray(p, dtype=np.float64) for p in (S, d1, dip, a, x0)]
    if any(p.ndim for p in pars):
        pars = [p[..., None] for p in pars]
    S, d1, dip, a, x0 = pars
    rad = np.radians(dip)
    cos, sin = np.cos(rad), np.sin(rad)
    u = (x - x0) * toward                     # from the trace, + towards the dip
    z = u / d1 - cos / sin                    # from above the deep end of the locked part
    atz, den = np.arctan(z), 1.0 + z * z
    amp = sense * S / np.pi
    if isinstance(kind, str):
        if kind == 'h':
            return -amp * (cos * atz + (sin - z * cos) / den) * toward + a
        if kind == 'z':
            return amp * (sin * atz + (cos + z * sin) / den) + a
        raise ValueError(f'unknown kind {kind!r}: h, z or (ph, pu)')
    ph, pu = (float(c) for c in kind)
    horizontal = -amp * (cos * atz + (sin - z * cos) / den) * toward
    vertical = amp * (sin * atz + (cos + z * sin) / den)
    return ph * horizontal + pu * vertical + a


def forward(x, S, d1, a=0.0, x0=0.0, creep=None, d2=None, rot=0.0, rot_side=1):
    """Fault-parallel velocity of a screw dislocation (+ shallow creep, rotation).

    Scalar parameters give an array shaped like x; arrays of K parameters give
    (K, x.size) - K models at once.
    """
    x = np.asarray(x, dtype=np.float64)
    pars = [np.asarray(p, dtype=np.float64) for p in (S, d1, a, x0)]
    if creep is not None:
        pars += [np.asarray(creep, dtype=np.float64), np.asarray(d2, dtype=np.float64)]
    if any(p.ndim for p in pars):
        pars = [p[..., None] for p in pars]
    S, d1, a, x0 = pars[:4]
    u = x - x0
    v = S / np.pi * np.arctan(u / d1) + a
    if creep is not None:
        C, d2 = pars[4], pars[5]
        v = v + C * np.sign(S) * (np.heaviside(u, 0.5) - 0.5 - np.arctan(u / d2) / np.pi)
    if rot:
        v = v + rot * np.where(u * rot_side > 0, u, 0.0)
    return v


##############################################################################
#                             the data and noise                             #
##############################################################################
def bin_profile(x, v, bin_km=1.0, min_count=3, min_bins=4):
    """Medians of a profile in bins along it.

    Returns dict: x (the median distance of the points in each bin), median,
    std (1.4826 x the median absolute deviation), q25, q75, count, sem
    (1.2533 std / sqrt(count), the standard error of a median) - only the bins
    with at least min_count points.
    """
    x, v = np.asarray(x, dtype=np.float64), np.asarray(v, dtype=np.float64)
    good = np.isfinite(x) & np.isfinite(v)
    x, v = x[good], v[good]
    if not x.size:
        raise ValueError('no data in the swath')
    edges = np.arange(np.floor(x.min() / bin_km) * bin_km, x.max() + bin_km, bin_km)
    if edges.size < 2:
        edges = np.array([x.min(), x.min() + bin_km])
    which = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    order = np.argsort(which, kind='stable')
    which, xs, vs = which[order], x[order], v[order]
    splits = np.flatnonzero(np.diff(which)) + 1
    out = {k: [] for k in ('x', 'median', 'std', 'count', 'q25', 'q75')}
    for part_x, part_v in zip(np.split(xs, splits), np.split(vs, splits)):
        if part_v.size < min_count:
            continue
        med = np.median(part_v)
        q25, q75 = np.percentile(part_v, [25, 75])
        out['x'].append(np.median(part_x))
        out['median'].append(med)
        out['std'].append(1.4826 * np.median(np.abs(part_v - med)))
        out['count'].append(part_v.size)
        out['q25'].append(q25)
        out['q75'].append(q75)
    out = {k: np.asarray(val, dtype=np.float64) for k, val in out.items()}
    if out['x'].size < min_bins:
        raise ValueError(f'only {out["x"].size} bin(s) with {min_count}+ pixels - a longer '
                         'profile, a wider swath or a smaller bin')
    # a bin of identical values has no spread: give it the typical one
    pos = out['std'] > 0
    floor = np.median(out['std'][pos]) if pos.any() else 1e-3
    out['std'] = np.where(pos, out['std'], floor)
    out['sem'] = 1.2533 * out['std'] / np.sqrt(out['count'])
    return out


def exponential_covariance(x, sill, length, nugget=None):
    """sill exp(-|xi - xj| / length), plus a diagonal nugget."""
    x = np.asarray(x, dtype=np.float64)
    lag = np.abs(x[:, None] - x[None, :])
    cov = sill * np.exp(-lag / max(float(length), 1e-6))
    if nugget is not None:
        cov = cov + np.diag(np.broadcast_to(np.asarray(nugget, dtype=np.float64), x.shape))
    return cov


def fit_noise_covariance(x, resid, bin_km=1.0):
    """(sill, length) of an exponential covariance, from residuals along a profile.

    The empirical covariance of the residuals in classes of separation (one
    bin wide, up to a third of the profile) is fitted by sill exp(-r / L).
    Residuals that do not decay give white noise (L ~ 0).
    """
    from scipy.optimize import curve_fit

    x, r = np.asarray(x, dtype=np.float64), np.asarray(resid, dtype=np.float64)
    r = r - r.mean()
    var = float(np.mean(r ** 2))
    white = (max(var, 1e-12), bin_km * 0.1)
    if r.size < 8 or var <= 0:
        return white
    lag = np.abs(x[:, None] - x[None, :])
    prod = r[:, None] * r[None, :]
    max_lag = (x.max() - x.min()) / 3.0
    centres, values = [0.0], [var]
    for lo in np.arange(bin_km * 0.5, max_lag, bin_km):
        sel = (lag > lo) & (lag <= lo + bin_km)
        if sel.sum() >= 6:
            centres.append(lo + bin_km / 2)
            values.append(float(prod[sel].mean()))
    centres, values = np.asarray(centres), np.asarray(values)
    if centres.size < 4:
        return white
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            (sill, length), _ = curve_fit(
                lambda h, s, ell: s * np.exp(-h / ell), centres, values,
                p0=(var, max(3 * bin_km, max_lag / 5)),
                bounds=([var * 0.05, bin_km * 0.1], [var * 3.0, max(max_lag * 3, bin_km)]))
        return float(sill), float(length)
    except (RuntimeError, ValueError):
        return white


##############################################################################
#                               the sampler                                  #
##############################################################################
def ensemble_sample(log_prob, p0, n_steps, stretch=2.0, rng=None, progress=None):
    """The affine-invariant ensemble sampler of Goodman & Weare (2010).

    Each walker of one half of the ensemble makes a "stretch move" along the
    line to a random walker of the other half, the two halves in turn - as
    emcee does.

    log_prob - function (K, n) -> (K,), vectorized over walkers
    p0       - (K, n) start positions, K even and >= 2 n
    Returns: chain (n_steps, K, n), log_prob (n_steps, K), acceptance fraction
    """
    rng = np.random.default_rng(rng)
    pos = np.array(p0, dtype=np.float64)
    n_walk, ndim = pos.shape
    if n_walk % 2 or n_walk < 2 * ndim:
        raise ValueError(f'need an even number of walkers, at least {2 * ndim}; got {n_walk}')
    lp = np.asarray(log_prob(pos), dtype=np.float64)
    if not np.all(np.isfinite(lp)):
        raise ValueError('every walker must start where the posterior is finite')
    chain = np.empty((n_steps, n_walk, ndim))
    lp_chain = np.empty((n_steps, n_walk))
    halves = (np.arange(n_walk // 2), np.arange(n_walk // 2, n_walk))
    accepted = 0
    for step in range(n_steps):
        for move, other in (halves, halves[::-1]):
            m = move.size
            z = ((stretch - 1.0) * rng.random(m) + 1.0) ** 2 / stretch
            partner = pos[other[rng.integers(0, other.size, m)]]
            prop = partner + z[:, None] * (pos[move] - partner)
            lp_new = log_prob(prop)
            take = np.log(rng.random(m)) < (ndim - 1) * np.log(z) + lp_new - lp[move]
            pos[move[take]] = prop[take]
            lp[move[take]] = lp_new[take]
            accepted += int(take.sum())
        chain[step], lp_chain[step] = pos, lp
        if progress is not None and (step + 1) % max(1, n_steps // 10) == 0:
            progress(step + 1, n_steps)
    return chain, lp_chain, accepted / float(n_steps * n_walk)


def autocorr_time(chain, c=5.0):
    """Integrated autocorrelation time of every parameter, in steps.

    As emcee: the autocorrelation function of each walker, averaged over the
    walkers, summed up to Sokal's automatic window (c times the estimate).
    chain - (n_steps, n_walkers, ndim)
    """
    n_steps, _, ndim = chain.shape
    nfft = 1 << int(2 * n_steps - 1).bit_length()
    tau = np.ones(ndim)
    for i in range(ndim):
        y = chain[:, :, i] - chain[:, :, i].mean(axis=0)
        spec = np.fft.rfft(y, n=nfft, axis=0)
        acf = np.fft.irfft(spec * np.conj(spec), n=nfft, axis=0)[:n_steps].mean(axis=1)
        if acf[0] <= 0:
            continue
        taus = 2.0 * np.cumsum(acf / acf[0]) - 1.0
        beyond = np.arange(n_steps) >= c * taus
        tau[i] = taus[np.argmax(beyond)] if beyond.any() else taus[-1]
    return tau


def summarize_samples(samples, best=None, names=None):
    """{name: dict(map, median, lo, hi)} - lo / hi the 2.5 / 97.5 percentiles."""
    names = names or [f'p{i}' for i in range(samples.shape[1])]
    out = {}
    for i, name in enumerate(names):
        lo, med, hi = np.percentile(samples[:, i], [2.5, 50, 97.5])
        out[name] = dict(map=float(best[i]) if best is not None else float(med),
                         median=float(med), lo=float(lo), hi=float(hi))
    return out


##############################################################################
#                                the fit                                     #
##############################################################################
def _as_sets(data):
    """[(x, v), ...] from one (x, v) pair or a list of them."""
    if (isinstance(data, tuple) and len(data) == 2 and np.ndim(data[0]) == 1
            and np.ndim(data[1]) == 1 and np.size(data[0]) == np.size(data[1])):
        return [tuple(np.asarray(d, dtype=np.float64) for d in data)]
    return [tuple(np.asarray(d, dtype=np.float64) for d in pair) for pair in data]


def _rot_term(model, u):
    if not model.rot:
        return 0.0
    return model.rot * np.where(u * model.rot_side > 0, u, 0.0)


def _grid_start(model, sets, gnss):
    """A robust start: for a grid of d1 (and x0), S and the offsets are linear -
    solve them by least squares and keep the best."""
    lo, hi = model.lo_hi()
    d1_lo, d1_hi = model.bounds['d1']
    d1_lo = max(d1_lo, 0.3)
    d1_grid = np.geomspace(d1_lo, max(min(d1_hi, 60.0), d1_lo * 1.01), 14)
    x0_grid = (np.linspace(*model.bounds['x0'], 9) if model.free_x0 is not None else [0.0])
    n_off = len(model.offsets) + (1 if gnss is not None else 0)
    # the scatter of the InSAR bins, to weigh the GNSS against them
    scatter = float(np.median(np.concatenate([b['std'] for b in sets]))) or 1.0
    best = None
    for d1 in d1_grid:
        for x0 in x0_grid:
            rows, rhs = [], []
            for k, b in enumerate(sets):
                u = b['x'] - x0
                col = np.zeros((u.size, 1 + n_off))
                col[:, 0] = np.arctan(u / d1) / np.pi
                col[:, 1 + k] = 1.0
                rows.append(col)
                rhs.append(b['median'] - _rot_term(model, u))
            if gnss is not None:
                u = gnss['x'] - x0
                w = scatter / gnss['sigma']
                col = np.zeros((u.size, 1 + n_off))
                col[:, 0] = np.arctan(u / d1) / np.pi
                col[:, -1] = 1.0
                rows.append(col * w[:, None])
                rhs.append((gnss['v'] - _rot_term(model, u)) * w)
            amat, bvec = np.vstack(rows), np.concatenate(rhs)
            sol = np.linalg.lstsq(amat, bvec, rcond=None)[0]
            rss = float(np.sum((amat @ sol - bvec) ** 2))
            if best is None or rss < best[0]:
                best = (rss, d1, x0, sol)
    _, d1, x0, sol = best
    start = {'S': sol[0], 'd1': d1, 'x0': x0}
    for k, name in enumerate(model.offsets):
        start[name] = sol[1 + k]
    if gnss is not None:
        start['ag'] = sol[-1]
    if model.creep:
        start['C'] = 0.05 * abs(sol[0])
        start['d2'] = min(2.0, 0.5 * d1)
    theta = np.array([start[n] for n in model.names], dtype=np.float64)
    margin = 1e-4 * (hi - lo)
    return np.clip(theta, lo + margin, hi - margin)


def _grid_start_dip(model, sets, gnss):
    """The start of a dip-slip fit: for a grid of d1, dip (and x0), S and the
    offsets are linear - solved by least squares, S >= 0, the best kept."""
    lo, hi = model.lo_hi()
    d1_lo, d1_hi = model.bounds['d1']
    d1_lo = max(d1_lo, 0.3)
    d1_grid = np.geomspace(d1_lo, max(min(d1_hi, 60.0), d1_lo * 1.01), 12)
    dip_grid = ([model.dip_fixed] if model.dip_fixed is not None
                else np.linspace(*model.bounds['dip'], 9)[1:-1])
    x0_grid = (np.linspace(*model.bounds['x0'], 7) if model.free_x0 is not None else [0.0])
    n_off = len(model.offsets) + (1 if gnss is not None else 0)
    scatter = float(np.median(np.concatenate([b['std'] for b in sets]))) or 1.0
    best = None
    for d1 in d1_grid:
        for dip in dip_grid:
            for x0 in x0_grid:
                rows, rhs = [], []
                for k, b in enumerate(sets):
                    col = np.zeros((b['x'].size, 1 + n_off))
                    col[:, 0] = dip_slip(b['x'], 1.0, d1, dip, 0.0, x0=x0, kind=model.kinds[k],
                                         toward=model.toward, sense=model.sense)
                    col[:, 1 + k] = 1.0
                    rows.append(col)
                    rhs.append(b['median'])
                if gnss is not None:
                    w = scatter / gnss['sigma']
                    col = np.zeros((gnss['x'].size, 1 + n_off))
                    col[:, 0] = dip_slip(gnss['x'], 1.0, d1, dip, 0.0, x0=x0,
                                         kind=model.gnss_kind, toward=model.toward,
                                         sense=model.sense)
                    col[:, -1] = 1.0
                    rows.append(col * w[:, None])
                    rhs.append(gnss['v'] * w)
                amat, bvec = np.vstack(rows), np.concatenate(rhs)
                sol = np.linalg.lstsq(amat, bvec, rcond=None)[0]
                if sol[0] < 0:                 # the wrong sense: S = 0, the offsets alone
                    sol = np.r_[0.0, np.linalg.lstsq(amat[:, 1:], bvec, rcond=None)[0]]
                rss = float(np.sum((amat @ sol - bvec) ** 2))
                if best is None or rss < best[0]:
                    best = (rss, d1, dip, x0, sol)
    _, d1, dip, x0, sol = best
    start = {'S': max(sol[0], 0.05 * hi[model.index('S')]), 'd1': d1, 'dip': dip, 'x0': x0}
    for k, name in enumerate(model.offsets):
        start[name] = sol[1 + k]
    if gnss is not None:
        start['ag'] = sol[-1]
    theta = np.array([start[n] for n in model.names], dtype=np.float64)
    margin = 1e-4 * (hi - lo)
    return np.clip(theta, lo + margin, hi - margin)


def fit_profile(data, bin_km=1.0, gnss=None, noise='auto', walkers=64, steps=5000, burn=None,
                seed=0, max_samples=48000, creep=False, free_x0=None, rot=0.0, rot_side=1,
                bounds=None, d1_prior=None, min_bins=3, verbose=True, fault='strike', kinds=None,
                gnss_kind='z', dip=(5.0, 85.0), toward=1):
    """Fit the dislocation model to fault-parallel velocity profiles.

    data   - (x, v), or a list of them - one per input file, each with its
             own offset. x: km from the fault, v: mm/yr, one value per pixel
    gnss   - None, or dict(x, v, sigma): GNSS, fitted with an offset of their own
    noise  - 'auto': exponential, from the residuals of a least-squares fit,
             one per file; 'white': the scatter of each bin, uncorrelated;
             or (sill, length) in (mm/yr)^2 and km, for every file
    walkers, steps, burn, seed - the ensemble sampler (burn: default steps / 3)
    creep, free_x0, rot, rot_side, bounds, d1_prior - see ProfileModel
    fault  - 'strike' (the default: ProfileModel), or 'reverse' / 'normal':
             DipSlipModel, with kinds (one per file: 'h', 'z' or (ph, pu)),
             gnss_kind, dip and toward - see DipSlipModel
    Returns dict: model, names, sets (indices of the files used), bins (one
    dict per file used), noise [(sill, length) per file], samples, log_prob,
    best (the MAP), summary, derived (strike-slip: the strain rate at the
    fault; dip-slip: the horizontal and vertical rates), acceptance, tau,
    lsq, n_pixels, walkers, steps, burn, gnss
    """
    from scipy.linalg import solve_triangular
    from scipy.optimize import least_squares, minimize

    say = print if verbose else (lambda *a, **k: None)
    data = _as_sets(data)
    bins, used = [], []
    for k, (x, v) in enumerate(data):
        try:
            bins.append(bin_profile(x, v, bin_km, min_bins=min_bins))
            used.append(k)
        except ValueError as err:
            say(f'model: file {k + 1} left out of the fit: {err}')
    if not bins or sum(b['x'].size for b in bins) < 5:
        raise ValueError('too few bins with data to fit a model - a longer profile, a wider '
                         'swath (-w) or a smaller --bin')
    if gnss is not None:
        gnss = {k: np.asarray(gnss[k], dtype=np.float64).ravel() for k in ('x', 'v', 'sigma')}
        sig = np.where(np.isfinite(gnss['sigma']) & (gnss['sigma'] > 0), gnss['sigma'], np.nan)
        fill = float(np.nanmedian(sig)) if np.isfinite(sig).any() else 1.0
        gnss['sigma'] = np.nan_to_num(sig, nan=fill)
        keep = np.isfinite(gnss['x']) & np.isfinite(gnss['v'])
        gnss = {k: val[keep] for k, val in gnss.items()} if keep.any() else None

    if fault not in FAULTS:
        raise ValueError(f'unknown fault type {fault!r}: one of {", ".join(FAULTS)}')
    if fault == 'strike':
        model = ProfileModel(n_sets=len(bins), creep=creep, free_x0=free_x0, rot=rot,
                             rot_side=rot_side, gnss=gnss is not None, bounds=bounds,
                             d1_prior=d1_prior)
    else:
        kinds = list(kinds) if kinds is not None else ['z'] * len(data)
        model = DipSlipModel(sense=1 if fault == 'reverse' else -1,
                             kinds=[kinds[k] for k in used], dip=dip, toward=toward,
                             free_x0=free_x0, gnss=gnss is not None, gnss_kind=gnss_kind,
                             bounds=bounds, d1_prior=d1_prior)
    # the offsets: anywhere the data could put them (unless given)
    s_lim = max(abs(model.bounds['S'][0]), abs(model.bounds['S'][1]))
    for name, b in zip(model.offsets, bins):
        if name not in (bounds or {}):
            span = max(float(np.ptp(b['median'])), 1.0)
            model.bounds[name] = (float(b['median'].min()) - 2 * span - s_lim,
                                  float(b['median'].max()) + 2 * span + s_lim)
    if gnss is not None and 'ag' not in (bounds or {}):
        span = max(float(np.ptp(gnss['v'])), 1.0)
        model.bounds['ag'] = (float(gnss['v'].min()) - 2 * span - s_lim,
                              float(gnss['v'].max()) + 2 * span + s_lim)
    lo, hi = model.lo_hi()

    # 1. least squares, from the best of a grid
    theta0 = (_grid_start if fault == 'strike' else _grid_start_dip)(model, bins, gnss)
    scatter = float(np.median(np.concatenate([b['std'] for b in bins]))) or 1.0

    def resid_lsq(theta):
        out = [model.predict(theta, b['x'], offset=name, kind=kind) - b['median']
               for name, b, kind in zip(model.offsets, bins, model.kinds)]
        if gnss is not None:
            out.append((model.predict(theta, gnss['x'], offset='ag', kind=model.gnss_kind)
                        - gnss['v']) / gnss['sigma'] * scatter)
        return np.concatenate(out)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        lsq = least_squares(resid_lsq, theta0, bounds=(lo, hi), x_scale='jac')
    theta_lsq = lsq.x.copy()
    if model.creep:
        i1, i2 = model.index('d1'), model.index('d2')
        theta_lsq[i2] = min(theta_lsq[i2], theta_lsq[i1] * 0.99)

    # 2. the noise of each file, applied by whitening through its Cholesky factor
    noises, whiten = [], []
    for name, b, kind in zip(model.offsets, bins, model.kinds):
        if isinstance(noise, (tuple, list)):
            sill, length = float(noise[0]), float(noise[1])
            nugget = b['sem'] ** 2
        elif noise == 'white':
            sill, length, nugget = 0.0, bin_km * 0.1, b['std'] ** 2
        else:
            resid = b['median'] - model.predict(theta_lsq, b['x'], offset=name, kind=kind)
            sill, length = fit_noise_covariance(b['x'], resid, bin_km)
            nugget = b['sem'] ** 2
        cov = exponential_covariance(b['x'], sill, length, nugget + 1e-8)
        try:
            chol = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            chol = np.linalg.cholesky(cov + np.eye(len(cov)) * 1e-6 * np.trace(cov) / len(cov))
        whiten.append(solve_triangular(chol, np.eye(len(cov)), lower=True))
        noises.append((sill, length))
    n_pix = int(sum(np.size(data[k][0]) for k in used))
    say(f'model: {sum(b["x"].size for b in bins)} bins of {bin_km:g} km from {n_pix:,} '
        f'pixels in {len(bins)} file(s)' + (f' + {gnss["x"].size} GNSS' if gnss else ''))
    for k, (sill, length) in zip(used, noises):
        say(f'   noise of file {k + 1}: ' + (
            'white (the scatter of each bin)' if noise == 'white' else
            f'exponential, sill {sill:.3g} (mm/yr)^2, length {length:.3g} km'))

    def log_prob(theta):
        theta = np.atleast_2d(theta)
        lp = model.log_prior(theta)
        ok = np.isfinite(lp)
        out = np.full(theta.shape[0], -np.inf)
        if not ok.any():
            return out
        th = theta[ok]
        ll = np.zeros(th.shape[0])
        for name, b, winv, kind in zip(model.offsets, bins, whiten, model.kinds):
            pred = np.atleast_2d(model.predict(th, b['x'], offset=name, kind=kind))
            z = (b['median'][None, :] - pred) @ winv.T
            ll -= 0.5 * np.sum(z ** 2, axis=1)
        if gnss is not None:
            gp = np.atleast_2d(model.predict(th, gnss['x'], offset='ag', kind=model.gnss_kind))
            ll -= 0.5 * np.sum(((gnss['v'][None, :] - gp) / gnss['sigma'][None, :]) ** 2,
                               axis=1)
        out[ok] = lp[ok] + ll
        return out

    # 3. the ensemble, started in a small ball around the mode of the posterior
    #    (the least-squares fit, moved by the priors and the noise covariance)
    start = theta_lsq
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        opt0 = minimize(lambda t: -log_prob(t)[0], theta_lsq, method='Nelder-Mead',
                        options=dict(xatol=1e-5, fatol=1e-6, maxiter=300 * model.ndim))
    if np.isfinite(opt0.fun) and -opt0.fun > log_prob(theta_lsq)[0]:
        start = opt0.x
    rng = np.random.default_rng(seed)
    walkers = max(int(walkers) + int(walkers) % 2, 2 * model.ndim + 2)
    scale = np.maximum(np.abs(start) * 0.01, (hi - lo) * 1e-3)
    p0 = start + scale * rng.standard_normal((walkers, model.ndim))
    margin = 1e-6 * (hi - lo)
    p0 = np.clip(p0, lo + margin, hi - margin)
    if model.creep:
        p0[:, model.index('d2')] = np.minimum(p0[:, model.index('d2')],
                                              p0[:, model.index('d1')] * 0.999)
    bad = ~np.isfinite(log_prob(p0))
    if bad.any():
        p0[bad] = start
    burn = steps // 3 if burn is None else int(burn)
    if not 0 <= burn < steps:
        raise ValueError(f'--burn ({burn}) must be at least 0 and below --steps ({steps})')

    def progress(done, total):
        say(f'   MCMC {done}/{total} steps ({walkers} walkers)')

    chain, lp_chain, accept = ensemble_sample(log_prob, p0, steps, rng=rng,
                                              progress=progress if verbose else None)
    tau = autocorr_time(chain[burn:])
    samples = chain[burn:].reshape(-1, model.ndim)
    lps = lp_chain[burn:].reshape(-1)
    if samples.shape[0] > max_samples:
        keep = rng.choice(samples.shape[0], max_samples, replace=False)
        samples, lps = samples[keep], lps[keep]

    # 4. the MAP: the best position the walkers visited, polished
    best = chain.reshape(-1, model.ndim)[np.argmax(lp_chain)].copy()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        opt = minimize(lambda t: -log_prob(t)[0], best, method='Nelder-Mead',
                       options=dict(xatol=1e-7, fatol=1e-9, maxiter=200 * model.ndim))
    if np.isfinite(opt.fun) and -opt.fun > log_prob(best)[0]:
        best = opt.x
    summary = summarize_samples(samples, best, model.names)

    def spread(values, at_map):
        lo_v, med_v, hi_v = np.percentile(values, [2.5, 50, 97.5])
        return dict(map=float(at_map), median=float(med_v), lo=float(lo_v), hi=float(hi_v))

    if fault == 'strike':
        # the shear strain rate at the fault, dv/dx at x = 0 = |S| / (pi d1)
        i_s, i_d = model.index('S'), model.index('d1')
        derived = {'strain': spread(np.abs(samples[:, i_s]) / (np.pi * samples[:, i_d]),
                                    abs(best[i_s]) / (np.pi * best[i_d]))}
    else:
        # the hanging wall against the footwall, far off: horizontal and vertical
        h_all, v_all = model.rates(samples)
        h_map, v_map = model.rates(best)
        derived = {'h_rate': spread(h_all, h_map), 'v_rate': spread(v_all, v_map)}

    say(f'   acceptance {accept:.2f}, autocorrelation time {np.max(tau):.0f} steps')
    for name in model.names:
        s = summary[name]
        say(f'   {name:3s} {s["map"]:9.3f}   95%: {s["lo"]:9.3f} to {s["hi"]:9.3f}'
            f'   {describe(name)}')
    if steps - burn < 50 * np.max(tau):
        say(f'   NOTE: {steps - burn} steps after the burn-in, under 50 autocorrelation times'
            f' ({50 * np.max(tau):.0f}): more --steps would make the ranges more reliable')
    return dict(model=model, names=list(model.names), sets=used, bins=bins, noise=noises,
                noise_mode=noise if isinstance(noise, str) else 'given',
                samples=samples, log_prob=lps, best=np.asarray(best, dtype=np.float64),
                summary=summary, derived=derived, acceptance=accept, tau=tau,
                lsq=theta_lsq, n_pixels=n_pix, walkers=walkers, steps=steps, burn=burn,
                gnss=gnss, bin_km=bin_km)


##############################################################################
#                the creep rate straight from the data                       #
##############################################################################
def _lines(x, v):
    """Least-squares slope and intercept of v(x) along the last axis."""
    mx, mv = x.mean(axis=-1), v.mean(axis=-1)
    dx = x - mx[..., None]
    sxx = (dx * dx).sum(axis=-1)
    sxv = (dx * (v - mv[..., None])).sum(axis=-1)
    slope = np.where(sxx > 0, sxv / np.where(sxx > 0, sxx, 1.0), 0.0)
    return slope, mv - slope * mx


TWO_LINES = ('slope_a', 'slope_b', 'intercept_a', 'offset')   # fit_two_lines' parameters


def two_line_matrix(x_a, x_b):
    """The design matrix of Hussain et al. (2016, eq. 2): a row per point -
    [x, 0, 1, 0] for side a, [0, x, 1, 1] for side b - for the parameters
    (slope_a, slope_b, intercept_a, offset): side a is the line
    slope_a x + intercept_a, side b slope_b x + intercept_a + offset."""
    x_a = np.asarray(x_a, dtype=np.float64).ravel()
    x_b = np.asarray(x_b, dtype=np.float64).ravel()
    mat = np.zeros((x_a.size + x_b.size, 4))
    mat[:x_a.size, 0] = x_a
    mat[x_a.size:, 1] = x_b
    mat[:, 2] = 1.0
    mat[x_a.size:, 3] = 1.0
    return mat


def fit_two_lines(x_a, v_a, x_b, v_b, n_boot=2000, seed=0, min_points=3):
    """The offset at x = 0 between two straight lines, one either side of a
    fault - side a (x < 0) and side b (x > 0) - as Hussain et al. (2016, eq.
    2) measure fault creep:

        v_a = slope_a x + intercept_a                 (side a)
        v_b = slope_b x + intercept_a + offset        (side b)

    solved two ways. Least squares (two_line_matrix): the best linear unbiased
    estimator (BLUE) for independent errors of one variance, its covariance
    sigma^2 (G^T G)^-1 from the residuals - too narrow when the errors are
    correlated or large (their Fig. S3). And a percentile bootstrap that
    resamples each side on its own, n_a of n_a and n_b of n_b (Bekaert et al.
    2015): every draw gives all four parameters, their spread the one
    Hussain et al. use.

    Returns dict - params (4,) and cov (4, 4): the BLUE, as TWO_LINES; sigma:
    the scatter about the lines; boot (n_boot, 4): the bootstrap draws;
    offset (the BLUE) and offset_std (its own 1 sigma); mean, median, std,
    lo / hi (2.5 / 97.5 %), lo68 / hi68 (16 / 84 %) of the offset draws;
    n_a, n_b; line_a, line_b as (slope, intercept). None with fewer than
    min_points a side.
    """
    x_a, v_a = np.asarray(x_a, dtype=np.float64), np.asarray(v_a, dtype=np.float64)
    x_b, v_b = np.asarray(x_b, dtype=np.float64), np.asarray(v_b, dtype=np.float64)
    n_a, n_b = x_a.size, x_b.size
    if min(n_a, n_b) < min_points:
        return None
    mat, data = two_line_matrix(x_a, x_b), np.r_[v_a, v_b]
    params = np.linalg.lstsq(mat, data, rcond=None)[0]
    resid = data - mat @ params
    var = float(resid @ resid) / max(data.size - 4, 1)
    try:
        cov = var * np.linalg.inv(mat.T @ mat)
    except np.linalg.LinAlgError:
        cov = np.full((4, 4), np.nan)

    rng = np.random.default_rng(seed)
    boot = np.empty((n_boot, 4))
    chunk = max(1, min(n_boot, 1_000_000 // max(n_a, n_b)))
    for start in range(0, n_boot, chunk):
        m = min(chunk, n_boot - start)
        ia, ib = rng.integers(0, n_a, (m, n_a)), rng.integers(0, n_b, (m, n_b))
        slope_a, icpt_a = _lines(x_a[ia], v_a[ia])
        slope_b, icpt_b = _lines(x_b[ib], v_b[ib])
        boot[start:start + m] = np.c_[slope_a, slope_b, icpt_a, icpt_b - icpt_a]
    off = boot[:, 3]
    lo, lo68, med, hi68, hi = np.percentile(off, [2.5, 16, 50, 84, 97.5])
    return dict(params=params, cov=cov, sigma=float(np.sqrt(var)), boot=boot, names=TWO_LINES,
                offset=float(params[3]), offset_std=float(np.sqrt(cov[3, 3])),
                mean=float(off.mean()), median=float(med),
                std=float(off.std(ddof=1)) if n_boot > 1 else 0.0, lo=float(lo),
                hi=float(hi), lo68=float(lo68), hi68=float(hi68), n_a=n_a, n_b=n_b,
                line_a=(float(params[0]), float(params[2])),
                line_b=(float(params[1]), float(params[2] + params[3])))


def fault_offset(x_a, v_a, x_b, v_b, n_boot=2000, seed=0, min_points=3):
    """The step at x = 0 between two straight lines, as Hussain et al. (2016)
    measure fault creep (fit_two_lines - their eq. 2): one line through the
    points on one side (x < 0), one through those on the other (x > 0), each
    with its own slope; the offset is the second line minus the first at
    x = 0. Its spread: a percentile bootstrap that resamples each side on its
    own, n_a of n_a and n_b of n_b (Hussain et al. 2016, after Bekaert et al.
    2015) - the points taken as independent.

    Returns dict(offset, lo, hi, std, boot, line_a, line_b, n_a, n_b), the
    lines as (slope, intercept); None with fewer than min_points a side.
    """
    fit = fit_two_lines(x_a, v_a, x_b, v_b, n_boot=n_boot, seed=seed, min_points=min_points)
    if fit is None:
        return None
    return dict(offset=fit['offset'], lo=fit['lo'], hi=fit['hi'], std=fit['std'],
                boot=fit['boot'][:, 3].copy(), n_a=fit['n_a'], n_b=fit['n_b'],
                line_a=fit['line_a'], line_b=fit['line_b'])


def split_offsets(los_offsets, horizontal, up, min_angle=10.0):
    """LOS offsets of several geometries (ascending, descending) as one
    horizontal and one vertical offset, as Hussain et al. (2016) split creep:
    los_k = P * horizontal_k + U * up_k, by least squares - draw by draw when
    the offsets are bootstrap draws, (n_geometries, n_draw).

    horizontal - per geometry, the LOS unit vector . the horizontal direction
    up         - per geometry, the LOS unit vector's up component
    Returns (P, U), arrays like a row of los_offsets - or None when the
    geometries are too alike (their (horizontal, up) less than min_angle
    degrees apart) to tell the two apart.
    """
    g = np.column_stack([np.asarray(horizontal, dtype=np.float64),
                         np.asarray(up, dtype=np.float64)])
    if g.shape[0] < 2 or not np.all(np.isfinite(g)):
        return None
    unit = g / np.linalg.norm(g, axis=1)[:, None]
    cos = np.clip(np.abs(unit @ unit.T), 0.0, 1.0)
    if np.degrees(np.arccos(cos.min())) < min_angle:
        return None
    sol = np.linalg.lstsq(g, np.asarray(los_offsets, dtype=np.float64).reshape(len(g), -1),
                          rcond=None)[0]
    return sol[0], sol[1]


def _draws(result, n_draw, seed):
    samples = result['samples']
    rng = np.random.default_rng(seed)
    return samples[rng.choice(samples.shape[0], min(n_draw, samples.shape[0]), replace=False)]


def model_bands(result, x, levels=(95.0, 68.0), n_draw=1000, seed=1, offset='a', kind=None):
    """The MAP curve and, for each level, the central level % envelope of the
    posterior curves at x - {95.0: (lo, hi), 68.0: (lo, hi)}: the 95 % and
    68 % confidence bounds that Hussain et al. (2016, Fig. 9) shade. kind:
    the component of a dip-slip model ('h', 'z', (ph, pu))."""
    model = result['model']
    curves = np.atleast_2d(model.predict(_draws(result, n_draw, seed), x, offset=offset,
                                         kind=kind))
    bands = {}
    for level in levels:
        tail = (100.0 - float(level)) / 2.0
        lo, hi = np.percentile(curves, [tail, 100.0 - tail], axis=0)
        bands[float(level)] = (lo, hi)
    return np.asarray(model.predict(result['best'], x, offset=offset, kind=kind)), bands


def model_band(result, x, n_draw=1000, seed=1, offset='a'):
    """The MAP curve and the 2.5 / 97.5 % envelope of the posterior curves at x."""
    best, bands = model_bands(result, x, levels=(95.0,), n_draw=n_draw, seed=seed,
                              offset=offset)
    return (best, *bands[95.0])


def strain_rate(result, x, n_draw=1000, seed=1):
    """The shear strain rate dv/dx (mm/yr per km = microstrain/yr) of the
    model at x - the MAP and the 2.5 / 97.5 % envelope. The step of a creeping
    fault at x0 itself is left out."""
    model = result['model']

    def rate(theta):
        par = model.unpack(np.atleast_2d(theta))
        u = np.asarray(x, dtype=np.float64)[None, :] - par['x0'][:, None]
        S, d1 = par['S'][:, None], par['d1'][:, None]
        out = S / np.pi * d1 / (d1 ** 2 + u ** 2)
        if model.creep:
            C, d2 = (par['C'] * np.sign(par['S']))[:, None], par['d2'][:, None]
            out = out - C / np.pi * d2 / (d2 ** 2 + u ** 2)
        if model.rot:
            out = out + model.rot * (u * model.rot_side > 0)
        return out

    curves = rate(_draws(result, n_draw, seed))
    lo, hi = np.percentile(curves, [2.5, 97.5], axis=0)
    return rate(result['best'])[0], lo, hi

"""Combine-style CLs limits for a single-bin cut-and-count.

The counting version of pyfitter's sensitivity_scan/simple_limit_combine_card.py, with the
same conventions, so the two limits can be compared directly:

- POI is the signal strength mu: N_CE = mu * (CE yield on the card).
- lnN nuisances as in Combine: each yield scales by exp(sum_k theta_k * ln kappa_k), with a
  standard-normal constraint on every theta_k.
- Test statistic: one-sided profile likelihood q_mu = 2 (NLL(mu) - NLL(mu_hat)) for mu_hat < mu,
  else 0, with mu_hat >= 0.
- asymptotic: Combine AsymptoticLimits expected limits and bands from the background-only Asimov
  dataset, plus a per-toy asymptotic CLs cross-check (as the pyfitter default).
- empirical: Combine HybridNew-style toy CLs scan over a mu grid,
  CL_sb = P(q >= q_obs | S+B), CL_b = P(q >= q_obs | B).

As in pyfitter, toys are thrown with the nuisances at their nominal values and the constraints
centred on zero. The likelihood is a single Poisson term, so q_mu depends only on the observed
count n and is cached per n: the toys then cost almost nothing.

Usage:
    python counting_stats.py --card <tag>/cutcount/count_card.yaml [--method empirical]
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import scipy.optimize as opt
import scipy.stats as stats
import yaml

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

THETA_BOUND = 7.0  # same range as the theta_ parameters in simple_limit_combine_12d.py


class CountingModel:
    """Single-bin Poisson likelihood with Combine lnN nuisances.

    lambda(mu, theta) = mu * s * exp(theta . ln kappa_s) + sum_b y_b * exp(theta . ln kappa_b)
    NLL = lambda - n ln(lambda) + 0.5 * sum(theta^2)
    """

    def __init__(self, yields, signal, systematics=None, freeze_nuisances=False):
        if signal not in yields:
            raise ValueError(f'Signal process "{signal}" not in the card processes {list(yields)}')
        self.processes = list(yields)
        self.signal = signal
        self.yields = np.array([float(yields[p]) for p in self.processes])
        self.is_signal = np.array([p == signal for p in self.processes])
        self.ref_signal_rate = float(yields[signal])
        if self.ref_signal_rate <= 0:
            raise ValueError(f'Signal yield must be > 0, got {self.ref_signal_rate}')

        self.nuisances = []
        rows = []
        for name, syst in (systematics or {}).items():
            if syst.get('type', 'lnN') != 'lnN':
                logger.warning('Systematic %s has type %s, only lnN is used in a counting card, skipped', name, syst.get('type'))
                continue
            effects = syst.get('effects', {})
            row = np.array([np.log(float(effects[p])) if p in effects else 0.0 for p in self.processes])
            if np.any(row != 0):
                self.nuisances.append(name)
                rows.append(row)
        self.log_kappa = np.array(rows) if rows else np.zeros((0, len(self.processes)))
        if freeze_nuisances:
            self.nuisances = []
            self.log_kappa = np.zeros((0, len(self.processes)))
        self._q_cache = {}

    @property
    def n_theta(self):
        return len(self.nuisances)

    def background(self):
        return float(np.sum(self.yields[~self.is_signal]))

    def expected(self, mu, theta=None):
        """Expected count lambda(mu, theta)."""
        scale = np.where(self.is_signal, mu, 1.0)
        if theta is None or self.n_theta == 0:
            return float(np.sum(scale * self.yields))
        return float(np.sum(scale * self.yields * np.exp(theta @ self.log_kappa)))

    def _nll_and_grad(self, mu, theta, n):
        """NLL and its gradient wrt (mu, theta)."""
        per_process = self.yields * (np.exp(theta @ self.log_kappa) if self.n_theta else 1.0)
        scale = np.where(self.is_signal, mu, 1.0)
        lam = max(float(np.sum(scale * per_process)), 1e-300)
        nll = lam - n * np.log(lam) + 0.5 * float(theta @ theta)
        common = 1.0 - n / lam
        d_mu = common * float(np.sum(per_process[self.is_signal]))
        d_theta = common * (self.log_kappa @ (scale * per_process)) + theta if self.n_theta else np.zeros(0)
        return nll, d_mu, d_theta

    def _minimize(self, n, mu_fixed=None, mu_start=0.0):
        """Minimum NLL over theta (and mu >= 0 unless mu_fixed). Returns (nll, mu_hat)."""
        if mu_fixed is not None and self.n_theta == 0:
            return self._nll_and_grad(mu_fixed, np.zeros(0), n)[0], mu_fixed
        if mu_fixed is None and self.n_theta == 0:
            # Analytic: lambda = n at the minimum, clipped at mu = 0
            mu_hat = max(0.0, (n - self.background()) / self.ref_signal_rate)
            return self._nll_and_grad(mu_hat, np.zeros(0), n)[0], mu_hat

        theta_bounds = [(-THETA_BOUND, THETA_BOUND)] * self.n_theta
        if mu_fixed is not None:
            def fun(x):
                nll, _, d_theta = self._nll_and_grad(mu_fixed, x, n)
                return nll, d_theta
            result = opt.minimize(fun, np.zeros(self.n_theta), jac=True, method='L-BFGS-B', bounds=theta_bounds)
            return float(result.fun), mu_fixed

        def fun(x):
            nll, d_mu, d_theta = self._nll_and_grad(x[0], x[1:], n)
            return nll, np.concatenate([[d_mu], d_theta])
        x0 = np.concatenate([[max(0.0, mu_start)], np.zeros(self.n_theta)])
        result = opt.minimize(fun, x0, jac=True, method='L-BFGS-B', bounds=[(0.0, None)] + theta_bounds)
        return float(result.fun), float(result.x[0])

    def q_mu(self, mu_test, n):
        """One-sided profile likelihood test statistic for an observed count n (may be non-integer: Asimov)."""
        key = (round(float(mu_test), 10), float(n))
        if key in self._q_cache:
            return self._q_cache[key]
        mu_start = (n - self.background()) / self.ref_signal_rate
        best_nll, mu_hat = self._minimize(n, mu_start=mu_start)
        if mu_hat >= mu_test:
            q = 0.0
        else:
            conditional_nll, _ = self._minimize(n, mu_fixed=mu_test)
            q = max(0.0, 2.0 * (conditional_nll - best_nll))
        self._q_cache[key] = q
        return q

    def asimov_count(self):
        """Background-only Asimov count: the nominal expected background."""
        return self.expected(0.0)

    def sqrt_q_asimov(self, mu_test):
        return np.sqrt(self.q_mu(mu_test, self.asimov_count()))


def _bracket(func, target, start=1.0, limit=1e5):
    """Grow an upper bound until func(upper) passes target; func is increasing."""
    upper = start
    while func(upper) < target and upper < limit:
        upper *= 1.5
    if upper >= limit:
        raise RuntimeError(f'Could not bracket the limit below mu = {limit:g}')
    return upper


def asimov_expected_limits(model, cl):
    """Combine AsymptoticLimits expected limits: sqrt(q_mu,A) = Phi^-1(1 - alpha * Phi(N)) + N."""
    alpha = 1.0 - cl
    limits = {}
    for band in (-2, -1, 0, 1, 2):
        target = stats.norm.ppf(1.0 - alpha * stats.norm.cdf(band)) + band
        upper = _bracket(model.sqrt_q_asimov, target)
        limits[band] = opt.brentq(lambda mu: model.sqrt_q_asimov(mu) - target, 1e-6, upper, xtol=1e-4)
    return limits


def asymptotic_cls(model, mu_test, n):
    """Cowan et al. asymptotic CLs for count n, with sigma from the Asimov dataset at mu_test (as in Combine)."""
    q = model.q_mu(mu_test, n)
    sqrt_q = np.sqrt(q)
    sqrt_q_asimov = model.sqrt_q_asimov(mu_test)
    if sqrt_q_asimov <= 0:
        return 1.0
    if sqrt_q <= sqrt_q_asimov:
        p_sb = stats.norm.sf(sqrt_q)
        one_minus_pb = stats.norm.cdf(sqrt_q_asimov - sqrt_q)
    else:
        p_sb = stats.norm.sf((q + sqrt_q_asimov**2) / (2.0 * sqrt_q_asimov))
        one_minus_pb = stats.norm.sf((q - sqrt_q_asimov**2) / (2.0 * sqrt_q_asimov))
    return p_sb / one_minus_pb if one_minus_pb > 1e-12 else 1.0


def asymptotic_limit(model, n, cl):
    """Asymptotic CLs upper limit on mu for an observed count n."""
    target = 1.0 - cl
    mu_hat = max(0.0, (n - model.background()) / model.ref_signal_rate)
    lower = max(mu_hat, 1e-6)
    upper = _bracket(lambda mu: -asymptotic_cls(model, mu, n), -target, start=max(1.0, 2.0 * lower))
    return opt.brentq(lambda mu: asymptotic_cls(model, mu, n) - target, lower, upper, xtol=1e-4)


def toy_asymptotic_limits(model, n_toys, cl, seed):
    """Per-toy asymptotic limits on background-only Poisson toys (the pyfitter toy cross-check)."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(model.asimov_count(), size=n_toys)
    per_n = {n: asymptotic_limit(model, n, cl) for n in np.unique(counts)}
    return np.array([per_n[n] for n in counts]), counts


def empirical_cls_scan(model, mu_grid, n_toys, cl, seed, n_obs=None):
    """HybridNew-style toy CLs on a mu grid.

    Returns the CLs curves for the expected quantiles (from the background-only q distribution)
    and, if n_obs is given, the observed curve.
    """
    quantiles = {-2: 0.025, -1: 0.16, 0: 0.5, 1: 0.84, 2: 0.975}
    curves = {band: np.full(len(mu_grid), np.nan) for band in quantiles}
    observed = np.full(len(mu_grid), np.nan) if n_obs is not None else None
    for index, mu in enumerate(mu_grid):
        if mu <= 0:
            for band in curves:
                curves[band][index] = 1.0
            if observed is not None:
                observed[index] = 1.0
            continue
        rng_b = np.random.default_rng(seed + index * 20000)
        rng_sb = np.random.default_rng(seed + index * 20000 + n_toys)
        n_b = rng_b.poisson(model.expected(0.0), size=n_toys)
        n_sb = rng_sb.poisson(model.expected(mu), size=n_toys)
        q_by_n = {n: model.q_mu(mu, n) for n in np.unique(np.concatenate([n_b, n_sb]))}
        b_q = np.array([q_by_n[n] for n in n_b])
        sb_q = np.array([q_by_n[n] for n in n_sb])

        def cls_at(q_obs):
            cl_sb = float(np.mean(sb_q >= q_obs))
            cl_b = max(float(np.mean(b_q >= q_obs)), 0.5 / len(b_q))
            return cl_sb / cl_b

        # A limit at the lower quantile comes from a downward fluctuation, i.e. a large q
        for band, p in quantiles.items():
            curves[band][index] = cls_at(float(np.quantile(b_q, 1.0 - p)))
        if observed is not None:
            observed[index] = cls_at(model.q_mu(mu, n_obs))
        logger.info('Empirical scan point %d/%d: mu=%.5g, CLs(expected median)=%.5g',
                    index + 1, len(mu_grid), mu, curves[0][index])
    return curves, observed


def crossing(mu_grid, cls_values, target):
    """First mu where CLs falls through target, linearly interpolated (as in pyfitter)."""
    below = np.where((cls_values[:-1] > target) & (cls_values[1:] <= target))[0]
    if not len(below):
        return np.nan
    i = below[0]
    return float(np.interp(target, cls_values[i:i + 2][::-1], mu_grid[i:i + 2][::-1]))


def load_card(path):
    with open(path) as f:
        card = yaml.safe_load(f)
    yields = {name: float(info['yield']) for name, info in card['processes'].items()}
    signal = card.get('poi', {}).get('name') or card.get('metadata', {}).get('signal_process', 'CE')
    return card, yields, signal


def rmue_conversion(card, args):
    """Effective muons and signal efficiency for R_mue = N_CE / (N_mu * capture * efficiency)."""
    norm = card.get('normalization', {}) or {}
    n_muons = args.n_muons if args.n_muons is not None else norm.get('stopped_muons', 5.588e15)
    capture = args.capture_rate if args.capture_rate is not None else norm.get('captures_per_stopped_muon', 0.609)
    efficiency = args.signal_efficiency if args.signal_efficiency is not None else norm.get('signal_efficiency', 0.1535)
    return float(n_muons), float(capture), float(efficiency)


def run(card_path, args):
    card, yields, signal = load_card(card_path)
    model = CountingModel(yields, signal, card.get('systematics'), freeze_nuisances=args.freeze_nuisances)
    n_muons, capture, efficiency = rmue_conversion(card, args)
    effective_muons = n_muons * capture
    to_rmue = 1.0 / (effective_muons * efficiency)
    s_ref = model.ref_signal_rate
    n_obs = args.n_obs if args.n_obs is not None else card.get('observed')

    logger.info('Card: %s', Path(card_path).resolve())
    logger.info('Yields: %s', ', '.join(f'{p}={y:.4g}' for p, y in yields.items()))
    logger.info('Signal %s = %.4g events, total background = %.4g events', signal, s_ref, model.background())
    logger.info('Nuisances (lnN): %s', ', '.join(model.nuisances) if model.nuisances else 'none (frozen or absent)')
    logger.info('R_mu_e normalization: N_muons=%.6g, capture_rate=%.4g, signal_efficiency=%.6g, SES=%.4e',
                n_muons, capture, efficiency, to_rmue)

    out = {
        'card': str(Path(card_path).resolve()),
        'method': args.method,
        'confidence_level': args.cl,
        'signal_yield': s_ref,
        'background_yield': model.background(),
        'nuisances': model.nuisances,
        'n_muons': n_muons,
        'capture_rate': capture,
        'signal_efficiency': efficiency,
        'single_event_sensitivity': to_rmue,
    }

    # Asimov expected limits are cheap and always reported
    asimov_mu = asimov_expected_limits(model, args.cl)
    asimov_events = {band: mu * s_ref for band, mu in asimov_mu.items()}
    logger.info('Asimov expected limit (Combine AsymptoticLimits): %.4f events, R_mu_e = %.6e',
                asimov_events[0], asimov_events[0] * to_rmue)
    logger.info('Asimov 68%% band: [%.4f, %.4f] events, [%.6e, %.6e] R_mu_e',
                asimov_events[-1], asimov_events[1], asimov_events[-1] * to_rmue, asimov_events[1] * to_rmue)
    logger.info('Asimov 95%% band: [%.4f, %.4f] events, [%.6e, %.6e] R_mu_e',
                asimov_events[-2], asimov_events[2], asimov_events[-2] * to_rmue, asimov_events[2] * to_rmue)
    out['asimov'] = {f'band_{band:+d}': {'mu': float(asimov_mu[band]), 'events': float(asimov_events[band]),
                                         'rmue': float(asimov_events[band] * to_rmue)} for band in asimov_mu}

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.method == 'asymptotic':
        limits, counts = toy_asymptotic_limits(model, args.toys, args.cl, args.seed)
        median, low, high = np.percentile(limits, [50.0, 15.87, 84.13])
        logger.info('Toy cross-check (%d toys): median %.4f events, 68%% band [%.4f, %.4f] events',
                    len(limits), median * s_ref, low * s_ref, high * s_ref)
        logger.info('Toy median expected limit: R_mu_e = %.6e', median * s_ref * to_rmue)
        out['toys'] = {'n': int(len(limits)), 'median_events': float(median * s_ref),
                       'low_events': float(low * s_ref), 'high_events': float(high * s_ref),
                       'median_rmue': float(median * s_ref * to_rmue)}
        if n_obs is not None:
            obs_mu = asymptotic_limit(model, float(n_obs), args.cl)
            logger.info('Observed (n=%g) asymptotic limit: %.4f events, R_mu_e = %.6e', n_obs, obs_mu * s_ref, obs_mu * s_ref * to_rmue)
            out['observed'] = {'n': float(n_obs), 'mu': float(obs_mu), 'events': float(obs_mu * s_ref),
                               'rmue': float(obs_mu * s_ref * to_rmue)}
        np.savez(output_dir / 'count_cls_limits.npz', limits=limits, event_limits=limits * s_ref,
                 rmue_limits=limits * s_ref * to_rmue, toy_counts=counts,
                 asimov_bands=np.array([-2, -1, 0, 1, 2]),
                 asimov_events=np.array([asimov_events[b] for b in (-2, -1, 0, 1, 2)]))
        plot_toy_limits(limits * s_ref, asimov_events, output_dir / 'count_cls_limits.png')
    else:
        mu_max = args.mu_max if args.mu_max is not None else 2.0 * asimov_mu[2]
        mu_grid = np.linspace(args.mu_min, mu_max, args.mu_points)
        logger.info('Empirical toy CLs scan: %d points from mu=%.3f to %.3f, %d toys per ensemble',
                    len(mu_grid), args.mu_min, mu_max, args.toys)
        curves, observed = empirical_cls_scan(model, mu_grid, args.toys, args.cl, args.seed,
                                              n_obs=float(n_obs) if n_obs is not None else None)
        target = 1.0 - args.cl
        emp_mu = {band: crossing(mu_grid, curve, target) for band, curve in curves.items()}
        if not np.isfinite(emp_mu[0]):
            logger.warning('CLs scan did not cross %.5g; try a higher --mu-max', target)
        logger.info('Empirical toy CLs expected limit: mu=%.6g, N_CE=%.6g events, R_mu_e = %.6e',
                    emp_mu[0], emp_mu[0] * s_ref, emp_mu[0] * s_ref * to_rmue)
        logger.info('Empirical 68%% band: [%.4f, %.4f] events', emp_mu[-1] * s_ref, emp_mu[1] * s_ref)
        out['empirical'] = {f'band_{band:+d}': {'mu': float(mu), 'events': float(mu * s_ref),
                                                'rmue': float(mu * s_ref * to_rmue)} for band, mu in emp_mu.items()}
        if observed is not None:
            obs_mu = crossing(mu_grid, observed, target)
            logger.info('Observed (n=%g) empirical limit: %.4f events, R_mu_e = %.6e', n_obs, obs_mu * s_ref, obs_mu * s_ref * to_rmue)
            out['observed'] = {'n': float(n_obs), 'mu': float(obs_mu), 'events': float(obs_mu * s_ref),
                               'rmue': float(obs_mu * s_ref * to_rmue)}
        np.savez(output_dir / 'count_empirical_cls_scan.npz', mu_grid=mu_grid,
                 **{f'cls_expected_{band:+d}': curve for band, curve in curves.items()},
                 cls_observed=observed if observed is not None else np.array([]))
        plot_cls_scan(mu_grid * s_ref, curves, observed, target, args.cl, output_dir / 'count_empirical_cls_scan.png')

    results_path = output_dir / f'count_limits_{args.method}.yaml'
    with open(results_path, 'w') as f:
        yaml.safe_dump(out, f, sort_keys=False)
    logger.info('Wrote %s', results_path)
    return out


def plot_toy_limits(event_limits, asimov_events, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    median, low, high = np.percentile(event_limits, [50.0, 15.87, 84.13])
    plt.figure(figsize=(9, 6))
    plt.hist(event_limits, bins=35, color='mediumpurple', alpha=0.6, edgecolor='indigo')
    plt.axvline(median, color='black', linestyle='--', linewidth=2, label=f'Toy median: {median:.2f} events')
    plt.axvspan(low, high, color='green', alpha=0.2, label='68% band')
    plt.axvline(asimov_events[0], color='red', linewidth=2, label=f'Asimov expected: {asimov_events[0]:.2f} events')
    plt.xlabel('Expected upper limit on CE signal yield')
    plt.ylabel('Toys / bin')
    plt.title('Cut-and-count asymptotic CLs limits')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    logger.info('Saved limit distribution plot: %s', path)


def plot_cls_scan(events_grid, curves, observed, target, cl, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(9, 6))
    plt.fill_between(events_grid, curves[-2], curves[2], color='yellow', label='Expected 95%')
    plt.fill_between(events_grid, curves[-1], curves[1], color='limegreen', label='Expected 68%')
    plt.plot(events_grid, curves[0], 'k--', linewidth=2, label='Expected median')
    if observed is not None:
        plt.plot(events_grid, observed, 'o-', color='black', label='Observed')
    plt.axhline(target, color='red', linestyle='--', linewidth=2, label=f'CL={cl}')
    plt.xlabel('CE signal yield [events]')
    plt.ylabel('CLs')
    plt.title('Cut-and-count empirical toy CLs scan (Combine style)')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.ylim([0, 1.1])
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    logger.info('Saved CLs scan plot: %s', path)


def add_stats_arguments(parser):
    """Limit options, shared with run_cutcount.py."""
    parser.add_argument('--method', choices=['asymptotic', 'empirical'], default='asymptotic',
                        help='asymptotic: Asimov expected limits + per-toy asymptotic cross-check; empirical: toy CLs scan')
    parser.add_argument('--toys', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--cl', type=float, default=0.90)
    parser.add_argument('--mu-min', type=float, default=0.0)
    parser.add_argument('--mu-max', type=float, default=None, help='Upper end of the empirical mu scan (default: 2x the Asimov +2 sigma limit)')
    parser.add_argument('--mu-points', type=int, default=33)
    parser.add_argument('--n-obs', type=float, default=None, help='Observed count in the window, for an observed limit (default: card "observed", else expected only)')
    parser.add_argument('--n-muons', type=float, default=None, help='Stopped muons (default: card normalization, else 5.588e15)')
    parser.add_argument('--capture-rate', type=float, default=None, help='Captures per stopped muon (default: card normalization, else 0.609)')
    parser.add_argument('--signal-efficiency', type=float, default=None, help='CE efficiency for R_mue (default: card normalization, else 0.1535)')
    parser.add_argument('--freeze-nuisances', action='store_true')


def main():
    parser = argparse.ArgumentParser(description='Combine-style CLs limits for a single-bin cut-and-count card')
    parser.add_argument('--card', required=True, help='Counting card (written by run_cutcount.py)')
    parser.add_argument('--output-dir', default=None, help='Default: the directory of the card')
    add_stats_arguments(parser)
    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = str(Path(args.card).resolve().parent)
    run(args.card, args)


if __name__ == '__main__':
    main()

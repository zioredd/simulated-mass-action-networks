"""Generator for hidden mass-action chemical reaction networks (prototype, CC0 creator data).

A system has n species S0..S{n-1}; h of them are never measured. m reactions, each with
1-2 reactant molecules and 0-2 product molecules (elementary mass action, order <= 2).
E experiments with different initial concentrations are simulated with a stiff solver;
observed species are sampled at T random times with multiplicative + additive noise.

Identifiability safeguards (applied before release):
- every reaction must carry >= FLUX_MIN of the largest reaction's integrated flux in
  at least one experiment, otherwise it is resampled;
- no two reactions have identical net stoichiometry AND identical reactants;
- every species takes part in at least one reaction;
- each reactant complex (multiset) occurs in at most one reaction. Without this, mass-action
  data only identify the summed net change per complex and the reaction set is not unique;
- each reaction's net-change direction is unique among valid product sets for its reactants
  (e.g. S2+S2->S2 and S2+S2->0 are proportional and therefore excluded);
- systems whose simulation exceeds 5 s wall time (pathological stiffness) are resampled;
- every hidden species is a reactant in at least one reaction (added 2026-09-26 after review
  RN-17 found 14.1% of hidden systems where hidden production was unobservable).

Difficulty knobs live in CONFIG and are tuned on development systems only.
"""
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp

CONFIG = dict(n_species=(4, 7), n_hidden=(0, 2), n_reactions=(3, 7), n_exp=4, n_times=12,
              t_max=10.0, log_k=(-1.5, 1.0), noise_mult=0.05, noise_add=0.01, flux_min=0.05,
              p_zero_product=0.15)


def canon(r):
    reac, prod = r
    s = lambda ms: "+".join(f"S{i}" for i in sorted(ms)) or "0"
    return f"{s(reac)}->{s(prod)}"


def sample_reaction(rng, n):
    k = rng.integers(1, 3)
    reac = tuple(sorted(rng.integers(0, n, size=k)))
    if rng.random() < CONFIG["p_zero_product"]:
        prod = ()
    else:
        prod = tuple(sorted(rng.integers(0, n, size=rng.integers(1, 3))))
    if sorted(reac) == sorted(prod) or not unique_direction(reac, prod, n):
        return None
    return reac, prod


def valid_products(reac, n):
    reac = tuple(sorted(reac))
    out = [()] + [(i,) for i in range(n)] + [tuple(sorted((i, j))) for i in range(n) for j in range(i, n)]
    return [p for p in out if sorted(p) != sorted(reac)]


def net_vec(reac, prod, n):
    v = np.zeros(n)
    for i in reac:
        v[i] -= 1
    for i in prod:
        v[i] += 1
    return v


def unique_direction(reac, prod, n):
    """True if no other valid product set gives a positively proportional net vector."""
    v = net_vec(reac, prod, n)
    for p in valid_products(reac, n):
        if tuple(sorted(p)) == tuple(sorted(prod)):
            continue
        w = net_vec(reac, p, n)
        if np.allclose(v / np.linalg.norm(v), w / np.linalg.norm(w)):
            return False
    return True


def stoich(reactions, n):
    A = np.zeros((len(reactions), n))  # reactant orders
    N = np.zeros((len(reactions), n))  # net change
    for j, (reac, prod) in enumerate(reactions):
        for i in reac:
            A[j, i] += 1
            N[j, i] -= 1
        for i in prod:
            N[j, i] += 1
    return A, N


def rhs_factory(A, N, k):
    def f(t, x):
        x = np.maximum(x, 0.0)
        rate = k * np.prod(np.power(x[None, :], A), axis=1)
        return rate @ N
    return f


def simulate(reactions, k, n, x0s, times):
    A, N = stoich(reactions, n)
    f = rhs_factory(A, N, k)
    out, flux = [], []
    for x0, ts in zip(x0s, times):
        sol = solve_ivp(f, (0, CONFIG["t_max"]), x0, method="LSODA", t_eval=np.sort(ts), rtol=1e-6, atol=1e-9)
        if not sol.success or sol.y.shape[1] != len(ts):
            return None, None
        out.append(sol.y.T)
        dense = solve_ivp(f, (0, CONFIG["t_max"]), x0, method="LSODA", t_eval=np.linspace(0, CONFIG["t_max"], 200),
                          rtol=1e-6, atol=1e-9)
        X = np.maximum(dense.y.T, 0)
        rates = k[None, :] * np.prod(np.power(X[:, None, :], A[None]), axis=2)
        flux.append(np.trapezoid(rates, dense.t, axis=0))
    return out, np.array(flux)


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def guarded_simulate(reactions, k, n, x0s, times, seconds=5.0):
    """Discard pathologically stiff systems: simulation must finish within `seconds` wall time."""
    import signal
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return simulate(reactions, k, n, x0s, times)
    except _Timeout:
        return None, None
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def make_system(rng):
    c = CONFIG
    n = int(rng.integers(c["n_species"][0], c["n_species"][1] + 1))
    h = int(min(rng.integers(c["n_hidden"][0], c["n_hidden"][1] + 1), n - 2))
    m = int(rng.integers(c["n_reactions"][0], c["n_reactions"][1] + 1))
    for _ in range(200):
        reactions, keys = [], set()
        tries = 0
        while len(reactions) < m and tries < 500:
            tries += 1
            r = sample_reaction(rng, n)
            if r is None or r[0] in keys:
                continue
            reactions.append(r)
            keys.add(r[0])
        if len(reactions) < m:
            continue
        used = {i for r in reactions for i in r[0] + r[1]}
        if len(used) < n:
            continue
        k = 10 ** rng.uniform(*c["log_k"], size=m)
        x0s = [rng.uniform(0.2, 2.0, size=n) * (rng.random(n) < 0.8) for _ in range(c["n_exp"])]
        times = [rng.uniform(0.05, c["t_max"], size=c["n_times"]) for _ in range(c["n_exp"])]
        traj, flux = guarded_simulate(reactions, k, n, x0s, times)
        if traj is None or not np.all(np.isfinite(flux)) or flux.max() <= 0:
            continue
        rel = flux / flux.max(axis=1, keepdims=True).clip(1e-12)
        if (rel.max(axis=0) < c["flux_min"]).any():
            continue
        if max(float(np.abs(t).max()) for t in traj) > 1e3:
            continue
        # fairness: a hidden species must be a reactant in some (flux-passing) reaction, otherwise
        # reactions that merely produce it are observationally identical to ones that do not
        reactant_species = sorted({i for r in reactions for i in r[0]})
        h_eff = min(h, len(reactant_species))
        hidden = sorted(rng.choice(reactant_species, size=h_eff, replace=False).tolist()) if h_eff else []
        obs = [i for i in range(n) if i not in hidden]
        exps = []
        for x0, ts, tr in zip(x0s, times, traj):
            order = np.argsort(ts)
            y = tr[:, obs]
            y = y * (1 + c["noise_mult"] * rng.standard_normal(y.shape)) + c["noise_add"] * rng.standard_normal(y.shape)
            exps.append(dict(x0=[float(v) for v in x0], t=[float(v) for v in np.sort(ts)],
                             y=[[float(v) for v in row] for row in y]))
        return dict(n_species=n, hidden=hidden, observed=obs, experiments=exps,
                    reactions=sorted(canon(r) for r in reactions), k={canon(r): float(kk) for r, kk in zip(reactions, k)})
    return None


def main():
    seed, count, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    rng = np.random.default_rng(seed)
    systems = []
    while len(systems) < count:
        s = make_system(rng)
        if s is not None:
            s["system_id"] = f"sys{seed}_{len(systems):05d}"
            systems.append(s)
    json.dump(systems, open(out, "w"))
    print(out, len(systems), "systems; mean reactions", np.mean([len(s["reactions"]) for s in systems]))


if __name__ == "__main__":
    main()

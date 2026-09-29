"""Mixed-difficulty system set: each system draws its own knobs (frozen candidate ranges)."""
import json, sys
import numpy as np
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gen

RANGES = dict(n_exp=(6, 10), n_times=(15, 25), noise_mult=(0.02, 0.03))


def main():
    seed, count, out, prefix = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    rng = np.random.default_rng(seed)
    S = []
    while len(S) < count:
        gen.CONFIG.update(n_exp=int(rng.integers(RANGES["n_exp"][0], RANGES["n_exp"][1] + 1)),
                          n_times=int(rng.integers(RANGES["n_times"][0], RANGES["n_times"][1] + 1)),
                          noise_mult=float(rng.uniform(*RANGES["noise_mult"])), noise_add=0.004, n_hidden=(0, 1))
        s = gen.make_system(rng)
        if s:
            s["system_id"] = f"{prefix}_{len(S):04d}"
            s["knobs"] = {k: gen.CONFIG[k] for k in ("n_exp", "n_times", "noise_mult")}
            S.append(s)
    json.dump(S, open(out, "w"))
    print(out, len(S))


if __name__ == "__main__":
    main()

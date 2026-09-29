# Simulated Mass-Action Reaction Network Experiments

## Overview

A collection of 3,150 randomly generated chemical reaction networks with simulated experiments. Each network has 4 to 7 species and 3 to 7 elementary mass-action reactions (one or two reactant molecules, zero to two product molecules). For each network, 6 to 10 experiments were simulated from different initial concentrations, and the concentrations of the observed species were sampled at 15 to 25 random times in (0, 10] with noise. In about half of the networks one species is never measured.

The data were created by the uploader with a purpose-built generator, which is included in this upload. They are not derived from, or a re-upload of, any existing dataset, and they are simulated rather than laboratory measurements.

## Intended Use

Benchmarking methods that infer the structure of a chemical reaction network (the exact set of reactions) from sparse, noisy, partially observed time series across several experiments: sparse regression, ODE fitting and model selection, amortized inference trained on the labelled networks, and combinatorial search. The independent unit is the whole network; splits should keep every network's experiments together.

## Creation Procedure

1. Draw the number of species (4-7), reactions (3-7) and hidden species (0-1), and per-network difficulty settings: 6-10 experiments, 15-25 samples per experiment, 2-3% multiplicative noise plus 0.004 additive Gaussian noise.
2. Sample reactions with 1-2 reactant molecules and 0-2 product molecules, keeping at most one reaction per reactant combination and only reactions whose net change is not proportional to another valid reaction with the same reactants. Every species must take part in a reaction.
3. Draw rate constants log-uniformly in [10^-1.5, 10^1] and initial concentrations uniformly in [0.2, 2.0], with each species independently absent (0) with probability 0.2.
4. Simulate every experiment over t in [0, 10] with the LSODA stiff solver (rtol 1e-6, atol 1e-9); discard networks whose simulation fails, exceeds 5 seconds, or reaches concentrations above 1000.
5. Keep a network only if every reaction carries at least 5% of the largest reaction's integrated flux in some experiment.
6. Choose hidden species only among species that are consumed by some reaction; remove their measurements.
7. Sample observation times uniformly in (0, 10], add the noise, and store the true reaction list. Identifiers are random and carry no information.

## Known Limitations

- Simulated data: deterministic mass-action kinetics in a well-mixed volume, no temperature effects, no catalysis beyond elementary reactions, no measurement bias or drift.
- The networks are random, not drawn from real chemistry, so chemical plausibility carries no information.
- Reaction rates are not included, only the reaction list.
- Some measured values are slightly negative because of additive noise.
- Even with perfect methods a small fraction of reactions may be practically indistinguishable at the noisiest settings.

## Files

This repository contains these files at its root (the README is this description).

- **`systems.jsonl` — JSON Lines, 3,150 lines:** one generated network per line.
- **`gen.py` — Python:** network sampling, identifiability filters, stiff ODE simulation and noise model.
- **`make_mixed.py` — Python:** draws each network's difficulty settings and calls `gen.py`.
- **`make_raw.py` — Python:** assembled `systems.jsonl` from the generator output and assigned random identifiers.
- **`requirements.txt` — text:** generator environment (Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1).
- **`RELEASE.json` — JSON:** release version, record count and SHA-256 hashes of every file.
- **`RAW_UPLOAD.md` — Markdown:** generator usage and what can be reproduced.
- **`LICENSE.txt` — text:** CC0 1.0 dedication for the data and the code.

The generator can produce a fresh demonstration with `python make_mixed.py 12345 5 demo.json demo` after installing `requirements.txt`; that demonstration is not the frozen corpus. The generation seeds of the frozen corpus are retained by the creator and are not included.

## Fields

Each line of `systems.jsonl` is a JSON object with these fields.

- **`raw_id` — string:** random identifier of the generated network.
- **`n_species` — integer:** number of species, named `S0` to `S(n_species-1)`.
- **`hidden` — list of integers:** indices of species that are never measured.
- **`experiments` — list of objects:** one object per experiment, with `x0` (initial concentration of every species), `t` (sorted sample times) and `y` (one row of measured values per time, observed species in increasing index order).
- **`reactions` — list of strings:** the true reactions, for example `S0+S3->S1` or `S4->0`.

## Source and Licence

- Source: creator-generated for this dataset; this repository is the source page.
- Licence: CC0 1.0 Universal (public domain dedication) for the data and the generator code, https://creativecommons.org/publicdomain/zero/1.0/

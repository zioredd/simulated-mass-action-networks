# Raw upload: Simulated Mass-Action Reaction Network Experiments

This upload holds the dataset and the exact code that generated it.

- `systems.jsonl`: 3,150 generated networks, one JSON object per line (see the dataset description for the fields).
- `gen.py`: network sampling, identifiability filters, stiff ODE simulation and noise model.
- `make_mixed.py`: draws each network's difficulty settings (6-10 experiments, 15-25 samples per experiment, 2-3% multiplicative noise, 0.004 additive noise, 0-1 hidden species) and calls `gen.py`. Usage: `python make_mixed.py <seed> <count> <out.json> <id_prefix>`.
- `make_raw.py`: wrote `systems.jsonl` from `make_mixed.py` output. It keeps `n_species`, `hidden`, `experiments` and `reactions`, drops rate constants and difficulty settings, and assigns each line a `raw_id` of fresh random entropy. In the original workspace it lived in `src/`, with `gen.py` and `make_mixed.py` in `prototype/`.
- `requirements.txt`: the generator environment (Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1).
- `RELEASE.json`: counts and SHA-256 hashes of every file in this upload.
- `LICENSE.txt`: CC0 1.0 for the data and the code.

A fresh demonstration set can be generated with `pip install -r requirements.txt` and then `python make_mixed.py 12345 5 demo.json demo`. That demonstration is not the frozen corpus. The generation seeds of the frozen corpus are retained by the creator for audits and are not included; with them, `make_mixed.py` reproduces every network, experiment and measurement in `systems.jsonl` exactly (verified for all 3,150 lines). The `raw_id` values are random and cannot be regenerated.

The challenge preparation program reads only `systems.jsonl`. Participant tables are created separately by that program and are not part of this upload.

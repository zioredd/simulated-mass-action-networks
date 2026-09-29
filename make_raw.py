"""Build raw/systems.jsonl (the uploaded raw dataset) from the frozen generator.
Truth fields kept in raw: reactions (rates and difficulty knobs are dropped). raw_id is secret
random entropy (not derived from seeds), so knowing the generator and seeds does not map to ids.
Usage: python src/make_raw.py <out_dir> <seed1,seed2,...> <per_seed>
"""
import json, os, secrets, sys, subprocess
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "prototype"))


def main():
    out, seeds, per = sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3])
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "systems.jsonl"), "w", encoding="utf-8") as f:
        for seed in seeds:
            tmp = os.path.join(out, f".tmp_{seed}.json")
            subprocess.run([sys.executable, os.path.join(ROOT, "prototype", "make_mixed.py"), str(seed), str(per), tmp, f"r{seed}"],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for s in json.load(open(tmp)):
                f.write(json.dumps({"raw_id": secrets.token_hex(12), "n_species": s["n_species"], "hidden": s["hidden"],
                                    "experiments": s["experiments"], "reactions": s["reactions"]}) + "\n")
            os.remove(tmp)
    print("wrote", os.path.join(out, "systems.jsonl"))


if __name__ == "__main__":
    main()

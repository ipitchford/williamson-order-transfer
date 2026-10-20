"""Reject corrupted certificate obligations in normal and optimised Python.

All mutations are made in a temporary directory. Published inputs are read-only.
"""
import ast
import copy
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def change_term(polynomial, forest, delta):
    for term in polynomial:
        if term["forest"] == forest:
            value = Fraction(term["coefficient"]) + delta
            if value:
                term["coefficient"] = str(value)
            else:
                polynomial.remove(term)
            return
    polynomial.append({"forest": forest, "coefficient": str(delta)})


def cases(certificate, ranks):
    bad = copy.deepcopy(certificate)
    bad["coefficients"][0]["polynomial"][0]["coefficient"] = "2"
    yield "coefficient", bad, ranks, "Root recursion"

    bad = copy.deepcopy(certificate)
    bad["coefficients"][0]["exact_flow_value"] = "999"
    yield "exact-flow-value", bad, ranks, "check failed"

    bad = copy.deepcopy(certificate)
    # Preserve every fixed-first polarisation sum while violating a shuffle.
    for family in bad["tensor_formulas"]:
        words = [ast.literal_eval(entry["word"]) for entry in family["formulas"]]
        if len(words) == 3 and len(words[0]) == 3:
            pair = next(((i, j) for i in range(3) for j in range(i + 1, 3)
                         if words[i][0] == words[j][0]), None)
            if pair:
                left, right = (family["formulas"][index]["polynomial"] for index in pair)
                forest = left[0]["forest"]
                change_term(left, forest, Fraction(1))
                change_term(right, forest, Fraction(-1))
                break
    else:
        raise RuntimeError("No aab family found for shuffle corruption")
    yield "shuffle", bad, ranks, "Shuffle mismatch"

    bad_ranks = copy.deepcopy(ranks)
    bad_ranks[0]["rank"] += 1
    yield "rank", certificate, bad_ranks, "check failed"


def main():
    certificate = json.loads((HERE / "order6_symbolic_certificate.json").read_text())
    ranks = json.loads((HERE / "tensor_ranks.json").read_text())
    outcomes = []
    with tempfile.TemporaryDirectory(prefix="williamson-corruption-") as temporary:
        directory = Path(temporary)
        for label, data, rank_data, expected in cases(certificate, ranks):
            certificate_path = directory / (label + ".json")
            rank_path = directory / (label + "-ranks.json")
            certificate_path.write_text(json.dumps(data))
            rank_path.write_text(json.dumps(rank_data))
            for optimisation in (0, 1):
                command = [sys.executable, "-B", "-E"]
                if optimisation:
                    command.append("-O")
                command += [str(HERE / "verify_certificate.py"),
                            "--certificate", str(certificate_path),
                            "--ranks", str(rank_path)]
                result = subprocess.run(command, text=True, capture_output=True, timeout=60)
                if result.returncode == 0 or '"status": "PASS"' in result.stdout:
                    raise RuntimeError(f"Corruption accepted: {label}, optimise={optimisation}")
                if expected not in result.stderr:
                    raise RuntimeError(f"Unexpected rejection for {label}: {result.stderr}")
                outcomes.append({"corruption": label, "optimise": optimisation,
                                 "returncode": result.returncode,
                                 "expected_obligation_rejected": expected})
    print(json.dumps({"status": "PASS", "negative_controls": outcomes}, indent=2))


if __name__ == "__main__":
    main()

"""Replay all exact checks and corruption controls without modifying inputs.

Run `python3 -B code/replay.py` and `python3 -B -O code/replay.py`.
The interpreter optimisation flag is explicitly passed to each child and read
back from that child's runtime. NumPy is required for the direct matrix checks.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

HERE = Path(__file__).resolve().parent
CHILD = (
    "import json,runpy,sys;from pathlib import Path;"
    "sys.path.insert(0,str(Path(sys.argv[1]).resolve().parent));"
    "print('EP_RUNTIME '+json.dumps({'optimise':sys.flags.optimize,"
    "'python':sys.version,'executable':sys.executable}),flush=True);"
    "target=sys.argv.pop(1);runpy.run_path(target,run_name='__main__')"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    paths = sorted(path for path in HERE.rglob("*")
                   if path.suffix in (".py", ".json"))
    before = {str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in paths}
    for path in paths:
        if path.suffix == ".py" and any(isinstance(node, ast.Assert)
                                        for node in ast.walk(ast.parse(path.read_text()))):
            raise RuntimeError(f"Removable proof assertion remains: {path.name}")
    children = {}
    for relative in ("verify_certificate.py", "verify.py",
                     "independent/audit_order6.py", "independent/audit_tensor_ranks.py",
                     "test_corruptions.py", "order7_criterion.py"):
        command = [sys.executable, "-B", "-E"]
        if sys.flags.optimize:
            command.append("-" + "O" * sys.flags.optimize)
        command += ["-c", CHILD, str(HERE / relative)]
        environment = dict(os.environ)
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        result = subprocess.run(command, cwd=HERE, env=environment,
                                text=True, capture_output=True, timeout=180)
        if result.returncode:
            raise RuntimeError(f"{relative} failed:\n{result.stdout}\n{result.stderr}")
        runtime_line, body = result.stdout.split("\n", 1)
        if not runtime_line.startswith("EP_RUNTIME "):
            raise RuntimeError(f"Missing child runtime evidence: {relative}")
        runtime = json.loads(runtime_line.removeprefix("EP_RUNTIME "))
        if runtime["optimise"] != sys.flags.optimize:
            raise RuntimeError(f"Child optimisation mode differs: {relative}")
        report = json.loads(body)
        if report.get("status") != "PASS":
            raise RuntimeError(f"No verified PASS: {relative}")
        children[relative] = {"runtime": runtime, "report": report}
    after = {str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest()
             for path in paths}
    if before != after:
        raise RuntimeError("Replay changed source/certificate files")
    report = {"status": "PASS", "python": sys.version, "platform": platform.platform(),
              "optimise": sys.flags.optimize, "source_sha256": before,
              "input_files_unchanged": True, "children": children}
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        # Exclusive creation prevents overwriting earlier evidence.
        with args.output.open("x") as output:
            output.write(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()

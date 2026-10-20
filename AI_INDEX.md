# Agent navigation

## Identity and starting point

- Release slug: `williamson-order-transfer`
- Version: `1.0.0-candidate`
- Scholarly creator: Anonymous
- Main statement and proof: [paper.md](paper.md), [PDF](paper.pdf)
- Machine-readable claims: [CLAIMS.json](CLAIMS.json)
- Status: [STATUS.md](STATUS.md)
- Provenance/input hash: [PROVENANCE.md](PROVENANCE.md)
- Evidence limits: [ASSURANCE.md](ASSURANCE.md)
- Prior-art boundary: [novelty_report.md](novelty_report.md), [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md)
- Citation metadata: [CITATION.cff](CITATION.cff), [Zenodo metadata](.zenodo.json)
- Licensing: [LICENSES.md](LICENSES.md)

## Scientific map

Theorem 3 proves the all-degree repeated-tail contraction identity.
Theorem 6 proves order transfer through six for all finite real Williamson
lists satisfying classical order conditions. Theorem 7 is a sufficient
order-seven criterion with three **additional** tensor assumptions.

Do not infer automatic seventh-order transfer, historical priority,
prescribed-stage existence, stability or performance from these statements.
Do not reinterpret the rank-deficiency boundary as a counterexample.

## Replay

Run from repository root:

```sh
python3 -B code/replay.py
python3 -B -O code/replay.py
```

See [code/README.md](code/README.md) for exact scope and
[code/replay.py](code/replay.py) for the executable entry point.
The launcher checks each child's actual
optimisation mode and input hashes. `--output` creates a new receipt exclusively.
No existing receipt should be overwritten.

The universal theorem is a written mathematical argument. Finite exact checks,
AI review, formal verification and independent reproduction are distinct
assurance categories; only the first two are described by this package.

## Useful next research question

Determine whether the three Theorem 7 tensor equalities follow automatically
from the classical order-seven conditions on genuine Williamson lists.
The family-rank computation does not answer this question or prove independence.

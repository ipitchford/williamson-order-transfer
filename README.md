# Order transfer through six for Williamson low-storage Runge–Kutta methods

**Research candidate · Anonymous · version 1.0.0-candidate**

This paper explains when a low-storage Runge–Kutta method keeps its order after
its additive updates are replaced by matrix exponentials. The answer proved
here is uniform over every finite real Williamson coefficient list through
classical order six. No positivity or distinct-stage-node assumption is needed.

Start with [the manuscript](paper.md), especially its proof roadmap and the
fourth-order worked argument. [AI_INDEX.md](AI_INDEX.md) provides machine-agent
navigation; [CLAIMS.json](CLAIMS.json) records individual claim boundaries.

## Main results

- An all-degree repeated-tail tensor identity follows directly from the
  Williamson recurrence, without division or nonsingularity assumptions.
- Polarisation and shuffle relations imply unconditional order transfer
  through six, assuming the corresponding classical order conditions.
- At order seven, three additional explicit tensor equalities suffice.
  Their automatic validity is **not proved**. The work does not settle
  unrestricted-order transfer or provide a seventh-order counterexample.

The theorem does not construct a coefficient list of a prescribed order and
stage count, nor establish superior stability, error constants or speed.

## Reproduce the finite checks

Python 3.10+ and NumPy are needed for the complete suite:

```sh
python3 -B code/replay.py
python3 -B -O code/replay.py
```

Both commands should return JSON with `status: PASS`. They run exact rational
certificate checks, independent implementations, matrix-orientation checks and
deliberately corrupted-data rejection tests. They do not modify source inputs.
See [code/README.md](code/README.md) for dependencies, evidence counts and
exclusive receipt output options.

Passing finite checks supports the written proof; it is not proof-assistant
formalisation or external human peer review. See [ASSURANCE.md](ASSURANCE.md).

## Publication and provenance

Repository: <https://github.com/ipitchford/williamson-order-transfer>.
Version DOI: <https://doi.org/10.5281/zenodo.23501588>.
Citation metadata: [CITATION.cff](CITATION.cff). The canonical provider records
identify the public availability and exact assets of this candidate version.

See [STATUS.md](STATUS.md), [PROVENANCE.md](PROVENANCE.md),
[REVISION_RESPONSE.md](REVISION_RESPONSE.md) and [SOURCES.md](SOURCES.md).
Original prose and structured records are CC0-1.0; original code is MIT.
Third-party source works retain their own terms: [LICENSES.md](LICENSES.md).

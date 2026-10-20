# Bounded computational replay

On 20 October 2026, the publication-preparation replay completed successfully
under Python 3.14.7 with NumPy 2.4.2, in both normal and optimized modes.
Each of the six child processes reported the expected interpreter mode;
source and certificate hashes were unchanged after execution.

The six checks are documented in [code/README.md](code/README.md).
They cover the emitted order-six certificate, rational identities and ordered
series, separately implemented prefix coefficients and shuffle ranks, deliberate
certificate corruptions, and the conditional weight-seven criterion.
Four corruption classes were each rejected in normal and optimized mode.

The weight-seven check covers 115 families, identifies three unresolved
directions, verifies that the selected scalar coordinate removes each rank
deficiency, and checks the targets 1/280, 1/84 and 1/168. It does not show
automatic seventh-order transfer or construct a seventh-order coefficient list.

The executable entry point is [code/replay.py](code/replay.py). It can emit a
new JSON receipt with `--output PATH`; it refuses to overwrite an existing
receipt. The package manifest identifies shipped inputs. CI repeats the two
interpreter modes; its actual result is visible in the repository Actions tab.

These finite checks support the written argument and detect specified
transcription errors. They are not proof-assistant verification, external
human peer review, an exhaustive novelty search or empirical performance evidence.

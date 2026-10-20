# Exact checks

From the repository root:

```sh
python3 -B code/replay.py
python3 -B -O code/replay.py
```

Python 3.10 or later is required. The direct matrix and tableau checks use
NumPy; the remaining checks use only the Python standard library. The recorded
publication replay used Python 3.14.7 and NumPy 2.4.2.

To preserve a receipt, add `--output /path/to/new-receipt.json`. Receipt creation
is exclusive: an existing file will not be overwritten. The launcher checks
that every child actually uses the requested optimisation mode and that source
and certificate files remain unchanged. Correctness checks are explicit
exceptions, not removable Python assertions.

## What runs

- `verify_certificate.py` checks the emitted formulas without importing their
  generator: 197 ordered-tree polynomials, 151 tensor formulas, 87 polarisation
  equations, 212 shuffle equations, 47 nontrivial family ranks, and root
  recursions/exact-flow values.
- `verify.py` checks 1,360 rational identities and 788 direct ordered-series
  comparisons. Its noncommuting-matrix control distinguishes the specified
  multiplication convention from its reversal.
- `independent/audit_order6.py` uses a separate implementation to compare
  9,653 prefix coefficients, with additional row and tensor checks.
- `independent/audit_tensor_ranks.py` independently enumerates positional
  shuffles and uses rational elimination and Bareiss determinants.
- `test_corruptions.py` mutates temporary copies of coefficient, exact-flow,
  shuffle, and rank data. Each mutation must be rejected under normal and
  optimised Python. No published input is modified.
- `order7_criterion.py` checks all 115 weight-seven elementary-vector
  families, the three unresolved family directions, the selected-coordinate
  rank increases, and the three exact continuous targets. This supports the
  **conditional** order-seven criterion; it does not prove that classical
  order-seven conditions automatically imply those additional tests.

## Regeneration

The submitted certificate can be regenerated to a new path:

```sh
python3 -B code/universal_tree_reduction.py --order 6 --output /path/to/new-certificate.json
python3 -B code/verify_certificate.py --certificate /path/to/new-certificate.json
```

The order-seven generator stops at its first undetermined family. The stored
`order7_boundary.json` is a boundary receipt, not a certificate of order-seven
transfer or a counterexample.

These are exact finite checks supporting the written argument. They are not a
proof-assistant formalisation, independent experimental reproduction, or a claim
of human peer review.

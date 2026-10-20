# Weil window certificate (Arb / python-flint 0.9)

Setup: `python3 -m venv venv && venv/bin/pip install python-flint mpmath`
- `check_ef.py`: checks explicit-formula conventions against zeta zeros.
- `check_zeros.py`: checks the Galerkin matrix against the sum over zeros of |f̂(γ)|².
- `weil.py`: x-space Galerkin matrix of Q, giving λ_min upper bounds (`run_table.py` → `table.txt`).
- `envelope.py`: floating-point test of the envelope reduction.
- `check_consistency.py`: checks M_R ≤ Q against the x-space matrix.
- `cert.py`: rigorous certificate; the proof chain is in its docstring. Run: `venv/bin/python cert.py 192 1.75 800 1050 64 1e-21`
  → `CERTIFIED Q(f) >= 9.99e-22 ||f||^2` (log in `cert_1.75.log`; about 8 min on 1 core, 330 MB).
L must be a dyadic rational, so that a·r is exact (interval Bessel recurrences blow up otherwise).

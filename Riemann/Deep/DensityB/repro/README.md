# repro: numerics behind PROOF.md section 4 (our scripts; ANTEDB code not included)

- best_beta.py: runs ANTEDB's compute_best_beta_bounds (exponent pairs to depth 5) and writes best_beta.json (pieces x0, x1, c, m: beta(a) <= c + m a on [x0, x1]). Needs EXPDB_PY = path of github.com/teorth/expdb blueprint/src/python and pycddlib<3 (gmp headers).
- explore.py: own numeric pipeline (all ANTEDB literature large-value estimates, raising to powers, subdivision, HM with best beta, LV_zeta via beta, Heath-Brown 12th moment) for Corollary zero-large-cor-0, with and without the new Montgomery range. Usage: python3 explore.py 0.86 0.87
- antedb_check.py: ANTEDB's own Corollary 11.7 LV computation (compute_large_value_region + compute_sup_rho_on_tau) with and without the new hypothesis. Usage: EXPDB_PY=... uv run --no-project --python 3.11 --with "pycddlib<3" --with numpy --with sympy --with scipy --with matplotlib python3 antedb_check.py 0.865 0.87 2.34

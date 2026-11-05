# ooc/ -- memory-lean (out-of-core) Weil-positivity certificate (weilooc, 2026-11-05)

Same theorem and proof chain as ../PROOF.md; the changes are proved in DERIVATION.md.  Every process stays far below
1 GB (blocks of <= 480 x 480 fixed-point entries, primes processed in groups, files read by explicit seek/read), runs
single-threaded, and stops cleanly at a time budget; chain.sh restarts it (exit 3 = resume) under a hard alarm.

| file | role |
|---|---|
| fxmm.py | exact fixed-point matrix products: 45 primes < 2^21, float64 BLAS (exact), explicit CRT, one rounding |
| oocla.py | disk matrices (FileMat), views, exact block linear combinations, resumable block loops, norm bounds |
| gvals.c | Legendre values at the Gram nodes by the norm-one rotation recurrence (DERIVATION.md section 1) |
| bandplan.py, bvals.c | band nodes / error budget (as cert4.band_lean) and spherical Bessel values (as cert2.sbessel) |
| fxkern.c, fxkern.py | dense midpoint kernels on blocks: Cholesky, triangular solves (data only, never trusted) |
| oocert.py | driver: H1-H5, Grams, G1/G2, X (iterative refinement), W', Z, H', Cholesky, exact residual, margins |
| wgen.py | Schur weight (numerics; the rigorous bound is bc.schur, re-run in every certificate) |
| chain.sh | runner: one process at a time, perl alarm 2400 s, resume loop |
| make_config.py | run parameters for a tag (s from the rigorous Schur bound, gamma, Omega, u, N) |
| runs/ | config.json, state.json (all margins and error terms), log.txt, chain.out of every certificate run; runs/weights
  holds the Schur weights (the configs' wfile paths pointed to the scratch copies) |

Build (in a scratch directory; binaries are not kept in the repo):
  cc -O2 -o gvals gvals.c -I/opt/homebrew/include -L/opt/homebrew/lib -lflint -lgmp -lmpfr
  cc -O2 -o bvals bvals.c (same flags);  cc -O2 -shared -fPIC -o libfxkern.dylib fxkern.c (same flags)
  python3 -m venv venv && venv/bin/pip install python-flint==0.9.0 numpy   (numpy must use Accelerate/BLAS)
Run:  OOC_BIN=<scratch dir> chain.sh <rundir>   with <rundir>/config.json as in runs/.
Validation: mu = 5 and mu = 16 (log 4) reproduce cert4/bc.py (same H: lambda_min(H) 6.705e-20 / 4.825e-16 at mu = 5);
Gram and band matrices agree with grams_c / cert4.band_lean within their radii (L = 31/20, N = 800).

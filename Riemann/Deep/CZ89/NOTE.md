CZ89 audit note (2026-10-21)
Paper: Cobeli-Zaharescu, More than 8/9 zeros of zeta are simple or on the critical line, arXiv 2610.23467 (TeX read in full; Lamzouri 2609.02882 TeX checked for Lemma 3.1 and Prop 2.1).
Verdict: proof correct. The seven-point lemma is true but tight (margin 6.45e-5). The paper's computer evidence is not reproducible. I supplied an independent rigorous certificate.
Rerun (native python 3.14, numpy, python-flint 0.9.0; under 130 MB; about 1 min in total):
  cd cert; python table.py 13   (10 s, writes tables.npz; python table.py 15 writes tables_15.npz)
  python bnb.py full            (certificate; CZ_TARGET=0.04407 is the negative control; CZ_GLOBALM2=40 runs the paper's Taylor test)
  python bnb.py local 300 1e8 0.002   (cube around g*)
  python checks.py              (Arb scalars: C, U, f_* > 0.79, K_*(1/2), K_*(0.743), S(g*))
  python validate.py; python basins.py 11 400 4000 240 0.046   (sanity checks; basin census)
Scratch (.scratch/cz89: venv, tables, sources) deleted after the run.


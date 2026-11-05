THEOREM: For every odd m ≥ 1 and every n ≥ 0 the n-th Taylor coefficient c_n(m,m) of 1/q_{m,m}(z) is positive, where q_{m,m} is the denominator of the diagonal [m/m] Padé approximant to e^z. This is the diagonal case k = m of Fasi–Higham, SIMAX 40 (2019), Conj. 3.2. The proof is computer-assisted, using Arb ball arithmetic.

**Scope of the proof**
- **m ≥ 203:** analytic proof with explicit constants (PROOF.md §1–8). Every finite inequality is checked in ball arithmetic and holds uniformly for ν = m+½ ≥ 203.5. It covers all n, from n ≤ 2m through 2m ≤ n ≤ Cm² up to n ≥ Cm².
- **3 ≤ m ≤ 201:** certified separately for each m with the earlier ../FasiHigham/arb_cert.py. It uses certified root balls, then partial fractions for the tail, then an Arb power series; every m returns POS.
- **m = 1:** c_n = 2^{−n}.

**Method (m ≥ 203)**
1. **Exact expansion.** Let p(z) = q(−z). The Padé remainder gives e^z q = p − z^{2m+1}E with E_j > 0. Hence 1/q = (e^z/p)·Σ_l X^l, where X = z^{2m+1}E/p = πI_ν(z/2)/K_ν(z/2). Write ζ = z/(2ν). Then c_n(2ν)^n/Λ = I_L + Σ_l κ^l J_l. Each J_l is a contour integral with its own saddle point, at ζ + (2l+1)(√(1+ζ²)−1) = (n − l(2m+1))/ν. No localisation of the zeros of q is needed.
2. **Error bounds.** Liouville–Green (Olver) bounds for K_ν and I_ν: |h| ≤ e^{2V(U1)/ν} − 1 on explicit progressive paths. They are proved by a Volterra argument and spot-checked against exact Bessel values (ratio ≤ 0.50).
3. **Good term.** A rigorous saddle-point lower bound for the dominant J_l, with tabulated relative error e ≤ 0.49. Every other term gets an explicit exponential upper bound.
4. **Regimes.**
   - n ≤ 12: exact rational functions of m.
   - n ≥ 8 and n/ν < 60: 3000 τ-boxes. The worst ratio of bad terms to the good term is 5.5e−18, and the rate gap is ≥ 0.18.
   - n/ν ≥ 60: linear-in-τ bounds, with total ratio ≤ 8.9e−435.

**Files** (Riemann/Deep/FasiHigham2): PROOF.md (full proof); verify.py (re-runs every step; verify.log holds its last run); proof code fh2lib.py, good.py, goodsmall.py, mktab.py (builds tab203.json), assembly.py, large_tau.py, smalln.py, identities.py; cross-checks at m = 203 sanity.py, check_good2.py, check_neg_IL.py, check_bounds.py. Dependencies: python-flint ≥ 0.9, mpmath, sympy.

**Verification run** (2026-11-05, one core, python3 verify.py --sanity): every step passed in 1471 s. Proof steps: A identities 20 s, B n ≤ 12 3 s, C tables 429 s, D Liouville–Green grid 58 s, E τ-sweeps 82 s, F τ ≥ 60 <1 s, G m ≤ 201 424 s.
Cross-checks at m = 203: S1 lower bound ≤ exact c_n for n ≤ 4000 (7 s); S2 good-term bound ≤ exact J_l (219 s); S3 I_L and N<0 bounds ≥ numeric values (224 s); S4 bad-term bounds ≥ numeric values (6 s).

**Novelty**
- Scardino (arXiv:2610.39051) proves k = 0 and states that the diagonal remains open.
- Goodchild (Zenodo 23919117, github JacobGoodchild/findformula3) certifies odd m ≤ 361 by computer, one m at a time.
- A search on 2026-11-05 of arXiv for Nov 2026, Zenodo and the GitHub history found nothing newer.
- To our knowledge this is the first proof for all odd m.

**Caveats**
- The proof is computer-assisted. It relies on Arb/python-flint, on the classical Liouville–Green theorem (Olver; reproved in outline) and, for m ≤ 201, on the earlier certificate.
- Only one quick review of §3–7 has been done. It found two presentation slips, both now fixed. A full independent review is still pending.
- One apparent bad-term violation came from a 50-digit check with too little precision; recomputed at 1800 bits it holds.

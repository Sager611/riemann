AUDIT PASSED: WeilSimple's certificate is a faithful, rigorous proof that QW_λ (CCM eq. (3.10)) has a simple even ground state for every λ ∈ [√2, √20], with QW_λ^N even-simple for N ≥ N₀ at the 52 listed windows and the eight OEIS A401172 terms correct; no false step found.

# REPORT: line-by-line audit of "QW_λ has a simple even ground state for λ ∈ [√2, √20]" (wsgsaudit, 2026-11-20)

Audited: discoveries/weil-simple-ground-state.md, Riemann/Deep/WeilSimple/ (DERIVATION.md, weil_matrix.py, cert_simple.py, rscan.c, run_windows.py,
chain.py, table.py, bracket.py, abracket.py, a401172.py, uimprove.py, logs/), against the CCM TeX (arXiv 2511.22755) and the OEIS entry.
Proof chain with source lines: FULLPROOF.md. Own code and logs: this folder. Nothing in WeilSimple/ was edited; all runs used a /tmp copy.

## 1. Fidelity to CCM's form: correct
- Pole 2Re(f^(i/2) conj f^(-i/2)) = +2<f,cosh(x/2)>² (even) / -2<f,sinh(x/2)>² (odd); archimedean term ∫|f^|² 2θ'/2π with 2θ' = -log π + Re ψ(1/4+it/2);
  primes Λ(n)n^{-1/2}, n < c (the atom n = c = λ² sits at C_f(L) = 0). Basis φ_j, ψ_j = orthonormal basis of CCM's E_N = span{V_n : |n| ≤ N}.
- New analytic check: WS's x-space archimedean term equals -W_R of CCM's own matrix formula (weinfty, TeX l. 711) identically (FULLPROOF 1b).
- Numerics (own code): 25 entries (j ≤ 9, both blocks, c = 5) by quadrature with numerically computed autocorrelations: 2.5e-14 (double);
  CCM's t-integral itself for the (0,0) archimedean entry: 2.9e-11; the lattice form with my own μ(k) vs the matrix: 3e-13 relative (smooth f),
  1.9e-10 for the discontinuous f = φ_0 (F ~ 1/k, the regime the L² claim needs).

## 2. Lower bounds on all of L²: correct
- Lattice identity (all f ∈ L²(I), values in (-∞, ∞]), sector identity, minorant (dropping nodes with μ ≥ 1 / μ ≥ 0, Arb lower bounds of the
  weights), P0 = deff·Id + finite rank, Woodbury with Chekmachev's Lemma 5.3 bound, Cholesky residual criterion, and the new g column:
  all checked against the code line by line (FULLPROOF 2-3). Reviewer verdict on this chain, pasted with exact formulas: CORRECT.
- Infinitely many nodes: the tail bound uses h increasing and b/(b²+r²) decreasing (re-proved); rscan.c's a-priori rounding bound E holds
  empirically with ratio 0.385 at c = 19 up to k = 2·10^9 (rscan_dbg.c); flagged nodes are kept exactly, so the scan only decides membership.
- The "201-node bug" is a completeness failure, not a soundness one: K2 = 201 at c = 7 (odd) gives margin -2.66e-24 (no claim), K2 = 401 passes.
- Since (3.10) is lsc, CCM's closure (E a core) is ≥ it pointwise, so the bounds hold for CCM's A_λ whatever closure is meant.

## 3. Logic: correct
- U < min(L_odd, L_2) ⇒ the spectral subspace below min(L_odd, L_2) is 1-dimensional and even (needs neither discreteness nor CCM Thm 3.6).
- Parity commutes with Q; every per-sector min-max value is non-increasing in c (extension by zero keeps Q, including endpoint atoms);
  chain step U(c_i) < ε(c_{i+1}) covers [c_i, c_{i+1}]; 37 contiguous steps 2 → 20 (own exact checker chaincheck.py and WS's table.py).
- All-N claim: compressions have per-sector min-max values ≥ the full ones and g ∈ E_N for N ≥ NG, so QW_λ^N is even-simple; that is exactly
  CCM's hypothesis of Thm 5.10 (finmain, l. 1087), whose δ_N(ξ) = 1 normalization CCM derive from even-simplicity; WS's basis spans CCM's E_N.
  52 windows checked, N₀ ∈ [30, 160]; CCM's λ = 3, √12 have N₀ = 80 ≤ 120 (rerun from scratch).

## 4. OEIS A401172 terms: correct
- Entry fetched read-only: a(n) = ⌊-log10 κ(n)⌋, κ(n) = inf Q(f)/||f||² over smooth f supported in (-log(n)/2, log(n)/2); data n ≤ 9 = 2..37.
- Lower bound of κ(n): plain even certificate at ε plus L_odd(n) > U(n) ≥ ε; upper: rigorous Rayleigh quotients. All eight decisions
  (a(10..13,15,16,18,19) = 42, 48, 53, 58, 69, 74, 85, 90) and n = 2..9 recomputed from the logs; n = 16 rerun (identical margin).
  Wiktorowicz's bound Q(g)/||g||² ≤ 8.62132e-49 at n = 11 (OEIS comment) lies inside our bracket [2.7411e-49, 9.45177e-49].
- Gap filled here (not stated by WS): U bounds the OEIS infimum over C_c^∞ of the open interval because C_c^∞(I°) is a form core
  (dilation + mollification, FULLPROOF 5b).

## 5. Falsification: every attack failed in the way a sound certificate must
- ε above a truncated eigenvalue (an upper bound of the true one) at c = 5 and 31/4, odd and deflated: 4/4 certificates fail.
- Just outside the first step (ε = 2.00001 U(2), odd, c ∈ (9/4, 3], logs/negctl2.log): passes for c ≤ 29/12, where λ_1^o(E_60) ≥ 3.6e-3 > ε, and fails
  from c = 5/2, where λ_1^o(E_60) = 1.78e-3 < ε already forces failure: the certificate tracks the true boundary. Stretched to 2 → 3: fails. K2 = 201 at c = 7: fails (as logged). ε = 0.5·λ^N: passes (minorant within a factor 2).
- The new deflation code: with g ⊥ ξ_N (|<g, ξ_N>| = 1e-57) or M = 1e-40 at ε > U the deflated certificate fails; with the true g it passes.
- Smallest failing case found: none in the chain. The tightest certified margins are tiny on the O(1) scale of Id - W_up (e.g. 3.5e-28 at
  c = 31/4, odd) but are proved with residuals below 1e-70 at 256-576 bits; precision is not the bottleneck.

## 6. Reruns (FULLPROOF 6)
c = 2 (rerun.sh), steps 29/4 → 31/4 and 16 → 133/8, c = 9 and 12 from scratch (NG search included), c = 19 (L = 1.472, inside Botha's unaudited
(1.19, 1.498]), OEIS n = 16: all reproduce the logged U, ε, flags and margins. Peak RSS 268 MB (c = 19). Not rerun (above the 300 MB cap):
windows 18, 20 and the three steps 37/2 → 20 (312-676 MB); same code path, logs internally consistent, reproduced by the orchestrator on 6 Nov.

## 7. Nits (none affects a conclusion)
- U_up strings are 12-digit roundings of Arb upper endpoints and may sit up to 5e-12 (relative) below them; every comparison has margin ≥ 2.
- cert_simple.py's own default is still K2 = max(201, ⌈5c²/3⌉) (l. 238); every logged run passes K2 explicitly (401 or more), so it is unused.
- memcap.sh misses short runs (peak 0 MB); /usr/bin/time -l was used for footprints.
- Integer rows c = 2-6, 8, 10 ran with K2 = 201 (valid, weaker minorant); the discovery file's "401-node fix" applies to the other rows.

## 8. What survives
Everything claimed: the continuum statement on λ ∈ [√2, √20] (an independent confirmation of Botha's range; on log half-widths (1.19, 1.498]
now audited line by line), the explicit thresholds N₀ at 52 windows, the eight OEIS terms. Bounded range only; not formalized.

Files: FULLPROOF.md; indep.py, tcheck.py (fidelity); rscan_dbg.c, rscan_test.py (rounding); negctl.py, negctl2.py, defl_falsify.py (falsification);
rerun_step.sh, chaincheck.py (reruns, chain); latrough.py (lattice identity, rough f); logs/.

RH relevance: N13, N1

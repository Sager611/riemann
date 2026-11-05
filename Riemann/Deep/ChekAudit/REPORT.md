SOUND: no false or unjustified step in the proof (secs. 2-5 re-derived; the identity checked to 1e-28). Our own code re-certifies c = 5, 10, 16 to every printed digit, and all c = 45/47 checks agree. The full c = 45, 47 and 71 runs rest on the author's unpublished code.

# ChekAudit: V. Chekmachev, "Rigorous lower bounds for the truncated Weil form up to support length ln 71 via a lattice minorant" (Zenodo 23848597, 1 Nov 2026)

RH relevance: N1

## Verdict in short
- Claim: Q_c(f) >= eps_c ||f||^2 for all complex f in L^2(-L/2, L/2), L = ln c, for c = 2..31, 45, 47; even f also at c = 71. With monotonicity in c this gives Weil positivity on half-width (ln 47)/2 = 1.925 for all f, and (ln 71)/2 = 2.131 for even f.
- Q_c is exactly our functional (PROOF.md sec. 0) with mu = c and half-width a = (ln c)/2. Our own independent code computes their Galerkin values to every printed digit (c = 5, 6, 7, 8, 14, 31, 45).
- The proof chain (sec. 2-5) is correct. I re-derived every identity and checked the identity numerically to 1e-28.
- My own implementation of their certificate (Theorem 5.10 with the rule of Definition 6.1, Arb balls, rigorous scan) reproduces their margins and every set size of Table A.2 at c = 5, 10 and 16. It proves Q_16 >= 1.5128e-75 ||f||^2 in both parity sectors. We had certified only 1e-76.
- At c = 45 our independent Galerkin upper bounds (N up to 1200) are 1.79e-231. That is still 1.55 x eps_45, and they are still falling slowly (see sec. 3). No contradiction.
- At c = 47 our own scan of all 2e8 nodes with k <= 2e8 finds exactly their selected nodes: 5719 odd and 2134 even, none missing and none extra. This range holds 95% and 98% of all their selected nodes.
- Not re-run by us: the full c = 45, 47 and 71 certificates (matrix orders 4424-38377 at 1200-1450 bits, and scans of 3.5e10 and 1.7e13 nodes done on a GPU). The paper lists its trust assumptions in sec. 6.3. The code is unpublished ("after journal submission").
- Consequence for N1: the record passes to Chekmachev, with half-width 1.925 for all f and 2.131 for even f. Our 1.55 and 1.6 are superseded. His c = 25 and 31 certificates alone (both re-checked in Arb by the author) already beat 1.6: (ln 25)/2 = 1.609 and (ln 31)/2 = 1.717.


## 1. Definitions: Q_c is exactly our functional
| item | paper (sec. 2.2, Lemma A) | ours (WeilRecord/PROOF.md sec. 0) | same? |
|---|---|---|---|
| support | f in L^2(I), I = [-L/2, L/2], L = ln c | f supported in [-a, a] | yes: a = (ln c)/2, mu = e^{2a} = c |
| Fourier transform | f^(r) = int f e^{-irx} | F(t) = int f e^{itx} | yes: F(t) = f^(-t); every term is even in t |
| archimedean | (1/2pi) int h+(r) abs(f^(r))^2 dr, h+(r) = Re psi(1/4 + ir/2) - ln pi | (1/2pi) int a(t) abs(F)^2, a(t) = Re psi(1/4 + it/2) - log pi | yes. The x-space form (2.1) is checked against psi in Lemma A(a); I re-derived it from Gauss's formula |
| primes | - sum_{n<c} Lambda(n) n^{-1/2} (g(ln n) + g(-ln n)), g = f * conj f(-x) | - <f, C_p f>, with C_p = sum_{n<mu} Lambda(n) n^{-1/2}(S_{log n} + S_{-log n}) | yes: <f, (S_s + S_-s) f> = 2 Re g(s). The term n = c sits at the edge of supp g, where g = 0 |
| poles | 2 Re[f^(i/2) conj f^(-i/2)], equal to 2<f, cosh(x/2)>^2 for even f and -2<f, sinh(x/2)>^2 for odd f (Lemma 2.3) | 2 abs(int f cosh(x/2))^2 - 2 abs(int f sinh(x/2))^2 | yes. pdftotext drops the bar; Lemma 2.3 fixes the meaning |
| norm, order | Q_c >= eps Id means Q_c(f) >= eps int abs(f)^2 | Q(f) >= tau norm(f)^2 | yes |
| map | c = 16 is our mu = 16 (a = log 4); c = e^{3.1} = 22.2 is a = 1.55; c = e^{3.2} = 24.5 is a = 1.6; c = 47 is a = 1.925; c = 71 is a = 2.131 | | |

I also checked the definitions numerically. Independent x-space code (galerkin_x.py) builds the Galerkin matrix of Q_c in their cosine basis from eq. (2.1) by quadrature, with closed-form autocorrelations and no lattice, digamma or Cauchy formulas. It reproduces every printed digit of their Table A.3(c):

| c, N | paper lambda_min(Q_{c,N}) | ours |
|---|---|---|
| 5, 50 | 1.172002576e-17 | 1.17200257596e-17 |
| 5, 100 | 1.005022530e-17 | 1.00502253020e-17 |
| 6, 50 | 8.651689900e-23 | 8.65168990015e-23 |
| 7, 100 | 7.730541166e-28 | 7.73054116631e-28 |
| 8, 100 | 4.581271669e-33 | 4.58127166921e-33 |
| 14, 50 | 5.137368428e-61 | 5.13736842847e-61 |
| 14, 100 | 1.642099861e-64 | 1.64209986084e-64 |
| 31, 512 | 1.0007e-155 | 1.000672907e-155 (galerkin_fast.py: the same integrals in closed form) |
| 45, 300 | 2.40e-229 | 2.403779694e-229 |
| 45, 400 | 3.07e-231 | 3.066554641e-231 |

## 2. Step-by-step verdict
| step | paper | verdict | how checked |
|---|---|---|---|
| form domain, closedness, core | Lemma 2.2 | correct | standard (log-weighted Fourier space) |
| complex f and parity sectors | Lemma 2.3 | correct | real symmetric kernel, even measure, pole split |
| monotonicity in c | Prop. 2.6 | correct | for f supported in I_{c'} the atoms with c' <= n < c see g(+-ln n) = 0. This is the same argument as our REPORT5 |
| even sector suffices (all c) | Prop. 2.5 | not used for any finite-c claim | the paper says so (sec. 2.3) |
| lattice identity | Thm 3.1 | correct, exact for every f in L^2(I) | g is continuous on the circle of length 2L. Its Fourier coefficients abs(f^(r_k))^2 >= 0 sum to 2L norm(f)^2, so the series converges absolutely and uniformly (Fejer). The archimedean part follows by Tonelli. I re-derived m(r) = h+(r) - T(r) and the tail series (3.1)-(3.2). Numerically (identity.py, c = 13): relative difference 3e-28 (even) and 1e-32 (odd) for random f with three boundary moments; for generic f (no moments, K = 6000) 1.3e-5 and 4e-9, equal to the omitted tail; with the full symbol h+ instead of h+ - T the identity fails by 3e-4 to 8e-2 |
| two orthogonal node systems | Thm 3.2, Cor. 3.3 | correct | even nodes: Neumann cosines of period L; odd nodes: Dirichlet cosines; each system is complete in the even sector, sum (1/L) e_k e_k^T = Id/2 |
| odd-node replacement (the minorant) | Thm 4.1 | correct | termwise for mu(r_k) >= t1 outside the finite set K, then the Parseval identity of Thm 3.2 |
| monotone bound beyond the scan | Lemma 4.2 | correct | abs(T(r)) <= 2 sum c^{-a_n}(1/r + a_n/r^2), which decreases in r; h+' > 0 from Im psi' < 0. At nodes sin(r_k L) = 0 |
| dual check | Thm 5.1 | correct | standard: N = B B^T <= P iff lambda_max(B^T P^{-1} B) <= 1 |
| Cauchy forms of the inner products | Prop. 5.2, 5.9(c) | correct | I re-derived <phi_j, e_k>, <phi_j, v>, <e_k, v>, norm(v)^2 and the odd-sector versions |
| Woodbury upper bound W_up | Lemma 5.3 | correct | completion of the square, S > 0 |
| positivity tests | Lemma 5.4 (margins only), Lemma 5.5 (proof) | correct | Lemma 5.5 is elementary: row-sum bound on the residual. Lemma 5.4 cites Higham Thm 10.3 for completed factorisations; it is only used for printed margins |
| Stirling remainder | Lemma 5.6 | correct | re-derived from Binet's formula |
| rigorous scan | Def. 5.7, Lemma 5.8, sec. 6.3 | logically correct | CPU interval scan described in detail (outward rounding, FMA phase split). The GPU scan for c = 45, 47 and 71 is only stated ("organized differently, but all rounding errors are bounded above"). It is not repeated on a CPU, but is cross-checked by LLL/Fincke-Pohst enumeration (c = 47, 71) |
| certificate theorem | Thm 5.10 | correct | R - eps >= P - N is compared termwise; d_eff > 0 |
| computation c <= 31 | sec. 6.4 (i)-(v) | reproduced by us at c = 5, 10, 16 (both sectors at 5 and 16) | see sec. 3. The author re-assembled and re-checked all c <= 31 in Arb. A second implementation by another person matched c = 5, 7, 23, 29 |
| computation c = 45, 47 | sec. 6.4 (vi), 6.5 | not reproduced in full; consistent with every check we ran | the final assembly and Lemma 5.5 were run only in the main (Rust/astro-float) implementation; symbol, constants, weights and Cauchy sums were re-computed in Arb; no ECC memory; GPU integer GEMM assumed exact |
| c = 71, even only | Prop. 1.3 | same method; ad hoc precision (1450 bits) and eps = 1.5 Delta | the final Arb check (v) was not run. Margin 2.375e-373 at eps_71; matrix radii <= 1e-424 |
| eps_c >= 0.97 Theta_c | Thm 1.1(b) | correct as arithmetic | Theta_c recomputed with mpmath (theta.py) matches Table A.1 to 1e-10 at c = 2, 3, 5, 10, 16, 22, 23, 25, 31, 45, 47, 71 |

First false or unjustified step: none found. The weakest links are not mathematical but practical:
- (a) the GPU scan for c = 45, 47 and 71 is not described in detail (sec. 6.3, p. 21);
- (b) the final positivity check at c = 45 and 47 ran in a single unpublished implementation (sec. 9, "Limitations");
- (c) the hardware has no ECC memory.
The paper states all three openly as assumptions.

## 3. Cross-checks
Each line says what was computed, the result, and what it proves.
- Scale Theta_c (theta.py, mpmath, 60 digits): matches Table A.1 to <= 3.8e-10 relative at all 12 tested c. Hence the ratio claim eps_c >= 0.97 Theta_c is correct arithmetic.
- Their Galerkin values versus our independent x-space code: equal to every printed digit (table in sec. 1). This confirms the definitions and normalisation.
- Lattice identity, Thm 3.1 (identity.py, c = 13, 200 bits):
  - random even and odd f with three boundary moments: |lattice - direct|/direct = 3e-28, 8e-28 (even) and 1e-32, 4e-33 (odd);
  - generic f without moments (K = 6000 nodes): 1.3e-5 and 4e-9, which is exactly the size of the omitted tail;
  - the full symbol h+ in place of h+ - T: off by 3e-4 to 8e-2.
- Independent certificate, even sector (cert.py: Arb balls, rigorous 80-bit scan, Lemma 5.5). Every set size equals Table A.2.

  | c | Delta (ours, eps = 0) | Delta (paper) | margin at eps_c (ours, rigorous) | paper's rigorous margin |
  |---|---|---|---|---|
  | 5 | 4.002816064e-18 | 4.002816e-18 | 4.0028e-19 > 0 at eps_5 = 7.447489666e-18 | 4.002e-19 |
  | 10 | 5.026355931e-44 | 5.026356e-44 | 5.0264e-45 > 0 at eps_10 = 9.133835606e-44 | 5.026e-45 |
  | 16 | 7.293623148e-76 | 7.293623e-76 | 7.2936e-77 > 0 at eps_16 = 1.512800554e-75 | 7.293e-77 |

- Independent certificate, odd sector:
  - c = 5: margin 8.28e-16 > 0 at eps5- = 5.693e-15;
  - c = 16: margin 2.99e-72 > 0 at eps16- = 2.106e-71, with |J-| = 27 and the pole column, as in Table A.4.
  - So our code proves Q_5 >= 7.447e-18 Id and Q_16 >= 1.5128e-75 Id on all of L^2. cert2.py, the fast Cauchy-sum version, reproduces c = 16 exactly.
- Our own certificates against theirs (lower bounds; no conflict):
  - mu = 5: ours 1e-20 (smoke test; H-margin 6.705e-20), theirs 7.447e-18;
  - mu = 16: ours 1e-76 (T-margin 2.541e-78 after tau), theirs 1.513e-75;
  - both stay below the Galerkin upper bounds: 9.16e-18 at mu = 5 and 2.37e-75 at mu = 16.
- The brief's test "eps_c below our lambda_min at e^3.1 and e^3.2":
  - eps_23 = 6.02e-113 <= 4.6e-110 and eps_25 = 1.17e-123 <= 8.4e-123 (lambda_min falls with c, so eps_c must sit below lambda_min at every c' <= c): passed.
  - Caveat: our 4.6e-110 and 8.4e-123 are eigenvalues of our own certificate matrix H, built from a minorant of Q. They are lower-bound-type numbers, not upper bounds on lambda_min(Q), so this test cannot refute anything. The real upper-bound test is the next item.
- Galerkin upper bounds (galerkin_fast.py: closed forms of the same x-space integrals, validated against galerkin_x.py) against eps_c. A valid certificate needs eps_c <= lambda_min <= every Galerkin value.

  | c | N: lambda_min(Q_{c,N}) | eps_c | verdict |
  |---|---|---|---|
  | 5 | 200: 9.34496e-18; 400: 9.19890e-18; 800: 9.17004e-18; 1200: 9.16523e-18; 1600: 9.16354e-18 (about N^-2, limit about 9.161e-18) | 7.447e-18 (0.81 x limit) | consistent |
  | 16 | 400: 2.6195e-75; 800: 2.4402e-75; 1200: 2.4129e-75; 1600: 2.3746e-75 (irregular steps) | 1.5128e-75 (certified by us) | consistent |
  | 45 | 300: 2.4038e-229; 400: 3.0666e-231; 600: 2.2408e-231; 800: 1.9748e-231; 1000: 1.8676e-231; 1200: 1.7885e-231 | 1.1506e-231 (0.64 x the N = 1200 value) | consistent so far |

  At c = 45 the last step fell 4.4%. A naive N^-0.5 fit through the last three points would extrapolate below eps_45. That fit is not evidence: at c = 16 the steps are just as irregular (the drop from 1200 to 1600 is larger than from 800 to 1200), and there the lower bound eps_16 is proved by our own code. Settling c = 45 by Galerkin needs N >> 1200. Each N costs O(N^3) at 1100 bits (1200 took 7 min), which did not fit this run.
- c = 47 data audit (Zenodo data record 23847003, certificate zip):
  - manifest: k_scan = 50,977,471,206; kscan_e = 18,753,563,601; |K+| = 6063, |K-| = 6301, bad even 869, |G| = 371000, order 7170, rigorous margin 6.5486e-244. These equal Table A.2.
  - The counts of retained nodes are consistent: 12364 odd and 371869 even, giving 6047 odd and 2172 even selected nodes.
  - The scan-node counts quoted for c = 31, 45, 47 and 71 equal the counts predicted by the monotone bound of Lemma 4.2: k_scan ~ 2L e^{1+2 sigma_c}, with sigma_31 = 8.32, sigma_47 = 10.81, sigma_71 = 13.85.
  - scan47.py: a float64 scan with exact uint64 phases (error < 1e-10) over all 1e8 odd and 1e8 even nodes with k <= 2e8 flags exactly their selected nodes (5719 odd with mu < 1, 2134 even with mu < 0), with 0 missing and 0 extra. Our values agree with their 40-digit symbol values to <= 2.2e-10 on 356717 retained nodes.
- Not done: an independent certificate at c = 25 (cert2.py; matrix order 558, S of order 1058 at 800 bits). It did not finish within its 1700 s limit. Nothing was learned from it.

## 4. Companion paper (Zenodo 23848535)
It is a theory paper about how lambda(c) decays: the boundary expansion of the logarithmic Laplacian, the Hadamard formula d lambda/da = -2 beta^2, and a reduction of positivity for all c to a bounded "boundary balance" r_c. It claims no new positivity: the boundedness of r_c is stated as open. Not audited in depth. It changes nothing in the N1 record.

## 5. What this means for N1
- If the c = 47 computation is right, the record is Chekmachev's: Weil positivity on half-width (ln 47)/2 = 1.925 for all f (eps_47 = 1.874e-242), and (ln 71)/2 = 2.131 for even f. Everything we could check agrees.
- Already c = 25 and 31 of the same paper (a = 1.609 and 1.717) beat our 1.6. The author re-checked both in Arb, and a second implementation by another person matched c = 23 and 29. Our certificates at 1.55 and 1.6 are superseded either way.
- Open trust items, all stated in the paper:
  - the code is unpublished;
  - the c = 45, 47 and 71 final assembly and positivity ran in one Rust/astro-float implementation;
  - the scan for c = 45, 47 and 71 ran on a GPU only;
  - the hardware has no ECC memory.
- Next audit step: a full independent certificate at c = 47. It needs a scan of 3.5e10 nodes (about 10 h single-core with scan47.py) and a dense matrix of order 7170 at 1200 bits (about 1.2e11 multiplications).

## 6. Files
- REPORT.md (this file).
- evidence/: galerkin_x.log, gfast.log (all Galerkin runs), identity13.log, cert_a.log, cert_16.log, cert_16b.log, scan47_a.log. Scripts: theta.py, galerkin_x.py, galerkin_fast.py, identity.py, cert.py, cert2.py, scan47.py.
- Run any script with: uv run --no-project --with mpmath --with numpy --with python-flint python3 <script>, from the directory that holds the scripts.
- Downloads (PDFs, data zip) were deleted after use. State: .scratch/chekaudit/STATE.md.

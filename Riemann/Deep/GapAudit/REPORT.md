AUDIT PASSED: Λ > 3.324 holds as written — Bui–Milinovich's Lemma main admits the midpoint mollifier as printed, the exact identity holds, and the certificate reruns identically (every printed digit) and matches an independent integral derived from BM's TeX alone.

# gapaudit: line-by-line referee of our record Λ > 3.324 (node N8, falsification first)

Sources: BM = Bui–Milinovich arXiv 1410.3635v3 (latest; TeX `.scratch/gapaudit/bm/main.tex`, "l." = line, deleted after use);
BBLR = Bettin–Bui–Li–Radziwiłł arXiv 1609.02539v1 (TeX l.103-136, Thm 1.1; both TeX files deleted after use). Λ = limsup (t_{n+1}−t_n)·log t_n/(2π) over
ordinates of zeros ON the critical line (BM l.53, l.63-66): unconditional (lesson 127), not the RH quantity λ.

## Verdicts by item

| Item | Verdict | Evidence |
|---|---|---|
| Criterion: h1 > 1 ⇒ Λ > κ for F_mid = e^{ivtL}ζ(s)ζ(s+iκπ/L)M(s+iκπ/(2L)) | holds | BM l.164-210 uses only: F vanishes at t_n and t_n − κπ/L, Wirtinger (i) for complex f (l.139-155), end pieces; all unchanged |
| Lemma main admits α3 = −β3 = iκπ/(2L), same error, same length y = T^ϑ | holds, as printed | l.368 "shifts α_j, β_j ≪ L⁻¹", l.384 "uniformly"; α3, β3 tracked through l.531-612; BBLR Thm 1.1 is bilinear in two independent sequences ≪ a^ε |
| Midpoint factor = e^{−iκπϑ(x1+x2−x3−x4)/2} = e^{ic(r2−r1)/2}; c_1 keeps (v−S)² | holds | l.377 factor y^{−(α3+β3)x−α3(x1+x2)−β3(x3+x4)}; l.422-433 Σ∂ gives −S independent of the evaluation point |
| GapBar identity e^{ic(r2−r1)/2}U_L conj U_R = 4 sin(c(X1−X2+σ)/2) sin(c(X3−X4+σ)/2) | holds (exact) | algebra checked; 2000 random points at 50 digits: max error 1.9·10⁻⁵⁰ (`logs/idmaj.txt`) |
| Code computes BM's functional for real P, r = 2, κ, v, ϑ | holds | independent 7-D integral (`indep7.py`, BM TeX only) = certmid's Arb h1 to 10 digits at degrees 8 and 16 |
| Certified run certmid_r2_d80_k3.324 | holds, rerun identical | h1 = 1.00056020750 ± 4.5·10⁻¹⁵, H = 6.5058025·10⁻³⁷ ± 5.6·10⁻⁴⁷ (47 s, 207 MB) |
| Strict inequality decided by Arb | holds | `ok = H > 0 and Q1b > 0 and theta < 1/4` on Arb balls with tail radii added (certr.py l.89-104) |
| Float ceiling 3.32658 (r = 2) not claimed as proved | holds | RECORDS.md l.24 "Method ceiling (float, not proved)"; GapBar REPORT "3.32658 (r = 2, float)" |

## (1) Referee of DERIVATION.md against the TeX
- **§0, reduction.** BM l.171: F = e^{ivtL}ζ(½+it)ζ(½+it+iκπ/L)M(½+it). Under Λ ≤ κ consecutive critical-line zeros are ≤ (1+o(1))2κπ/L
  apart, so {t_n} ∪ {t_n − κπ/L} has gaps ≤ (1+o(1))κπ/L (l.182). The mollifier only adds zeros. Wirtinger (i) holds for complex f
  (l.155). End pieces (l.194): |M(½+it+iκπ/(2L))| ≤ Σ|d_r(h)P[h]|h^{−½} ≪ T^{ϑ/2+ε}, the same bound as unshifted, so
  |F^{(k)}|² ≪ T^{2/3+ϑ+ε} ≪ T^{1−ε}. Conclusion l.198-202: h1 > 1 ⇒ Λ > κ. Inequality directions: all ≤ as printed; correct.
- **§0, Lemma main (l.370-384).** Statement: "Suppose ϑ < 1/4", I_{α,β}(M) (l.364-367, mollifier shifts α3, β3 included,
  "where the shifts α_j, β_j ≪ L⁻¹") = c(α,β)·A_{r+2}(log y)^{r²+4r}L⁴ŵ(0)T/(2(r²−1)!((r−1)!)⁴) + O(TL^{(r+2)²−1}), "uniformly for
  α_j, β_j ≪ L⁻¹". With real P, conj M(½+α3+it) = M(½+ᾱ3−it), so |M(½+it+iκπ/(2L))|² is the case α3 = iκπ/(2L), β3 = −iκπ/(2L);
  |α3| = κπ/(2L) < |α2| = κπ/L, which BM themselves use (l.391). Same error terms, same mollifier length.
- **Proof of Lemma main.** α3, β3 enter only through h^{−α3}, k^{−β3} (l.531, l.538); Lemmas 4.1/4.2 (l.440-492, uniform in shifts
  ≪ 1/log y) turn them into y^{−(α3+β3)x−α3(x1+x2)−β3(x3+x4)} (l.598-605); I_2..I_6 arise from I_1 by swaps of α1, α2, β1, β2 only
  (l.610), so the factor is common to all six terms (l.612, l.650-666).
- **BBLR as quoted (l.500-517)** shows Σ a_h conj(a_k) with no coefficient hypotheses. The original (BBLR Thm 1.1, TeX l.103-136):
  A, B independent Dirichlet polynomials of length T^ϑ with α_a ≪ a^ε, β_b ≪ b^ε; shifts complex ≪ 1/log T; Φ smooth on [1,2],
  Φ^{(j)} ≪ T^ε; error O(T^{1/2+2ϑ+ε} + T^{3/4+ϑ+ε}) = o(T) for ϑ < 1/4. For the midpoint, a_h = d_r(h)P[h]h^{−iκπ/(2L)} has
  |a_h| = |d_r(h)P[h]| ≪ h^ε, and it is even the quadratic case a_h conj(a_k). Hypotheses met; ϑ = 2499/10000 < 1/4.
- **Derivative step (l.393-434).** F_mid'/(ie^{ivtL}) = L·Q(L⁻¹(∂_{α1}+∂_{α2}+∂_{α3}))[…] with Q(x) = v + x; Cauchy circles of radius
  ≍ 1/L around the shifted point stay in the uniform range. On the exponentials of c(α,β), L⁻¹Σ_j∂_{α_j} and L⁻¹Σ_j∂_{β_j} both give −S,
  S = ϑ(x+Σx_i) + (1−ϑ(x1+x3))t1 + (1−ϑ(x2+x4))t2, independent of the evaluation point, so c_1 = ∫(…)(v−S)². Correct as in DERIVATION §0.
- **The midpoint factor.** y^{−α3(x1+x2)−β3(x3+x4)} = exp(−i(κπ/2)(log y/L)(x1+x2−x3−x4)); log y/L = ϑ + O(1/L), and the
  replacement costs O(TL^{(r+2)²−1}) — the same simplification BM make for their T- and y-powers (l.434). = e^{ic(r2−r1)/2}, c = κπ.
- **Passage to the real kernel (§3, GapCeil §1-2).** Own derivation: integrating t3, t4 exactly turns BM's weights into
  (κπ)⁻²(1−e^{−ia})(1−e^{ib}), a = c(X3−X4+σ), b = c(X1−X2+σ); times e^{ic(X2−X4)} this is GapCeil's
  e^{ic(X1−X3)} + e^{ic(X2−X4)} − e^{ic(X1−X4)}e^{icσ} − e^{ic(X2−X3)}e^{−icσ}, in BM's own labels (X_i = ϑx_i); τ_i = A_i t_i absorbs
  BM's (1−ϑ(x1+x3))(1−ϑ(x2+x4)). So the twist e^{ic(X3+X4−X1−X2)/2} in optmid.py/gtw.py is applied to the right variables, and the
  identity and product form hold exactly (2cos(a−b) − 2cos(a+b+2θ) = 4 sin(a+θ) sin(b+θ)).

## (2) Code audit (certmid.py → certr.main with gtw.gamma_arrays_tw; gcert.py helpers; gcore.py float path)
- **Functional.** certr.main computes Q0 = ⟨φ, M0φ⟩, Q1 = v²Q0 − 2v·MS + MSS (MS = a1 + b0, MSS = a2 + 2b1 + c0), φ = I^{2r}f,
  f(U) = P(1−U/ϑ) = Σ A_n P_n(2U/ϑ−1) (checked: I[0][0] = ϑ/2, I[1][0] = −ϑ/2 is ∫_U^ϑ 1 = ϑ − U; Gram ϑ/(2k+1)), so
  P(w) = Σ(−1)^n A_n P_n(2w−1), real, degree 80 (81 coefficients). h1 = Q0/(κ²Q1) is scale-free (c_j = const·Q_j, same const).
- **Twist.** gtw adds (−c/2, +c/2) to the (X1, X3) and (X2, X4) phases of the four terms of certr.gamma_arrays_r = multiplication by
  e^{ic(X3+X4−X1−X2)/2}; optmid.py's float rho_tw multiplies by the same factor. Matches §1 above.
- **Exactness.** κ = 831/250, v = 683629/500000, ϑ = 2499/10000 exact rationals; A_n exact decimals (printed in the log); Legendre
  window n2 = N + D + 8 holds every I^iφ (i < D) exactly; X^q blocks taken from T^q in a larger window (exact).
- **Rounding/tail.** Arb balls throughout (prec 366). Tail i or j ≥ D: |Γ_ij| ≤ 4(i+j+8)²((i+4)(j+4))⁴γ^{i+j}, γ = 24 ≥ 2c+2 = 22.89,
  with ‖I^i‖ ≤ ϑ^i/i!, X^{q0+p} ≤ ϑ^{q0}; the sum is cut at i, j < 340 plus 10⁻⁶⁰ (the omitted terms are < 10⁻⁴⁰⁰). Base bound
  |û_ab| ≤ (2c+2)^{a+b} re-derived (Leibniz: phases |α| ≤ c, |∂^q m_k(1−s, ±ic)| ≤ (c+2)^q). On the computed arrays (D = 80) the
  ratio |Γ_ij|/majorant is ≤ 1.3·10⁻⁷ and falls to 4·10⁻¹⁰³ at i+j ≥ 140 for all three arrays a, b, c (`logs/idmaj.txt`).
  The certificate has slack: rad(H) = 5.6·10⁻⁴⁷ (which contains the tail terms e0 + κ²e1) against H = 6.5·10⁻³⁷, so the tail bound
  could be 10¹⁰ times too small and the certificate would still hold.
- **Reality.** certr symmetrises Re Γ. This is exact for real φ because c_0, c_1 are real (leading constants of real mean values);
  indep7.py's Im parts are ≤ 10⁻¹² relative. (gtw's IMCHECK tests the antisymmetric part of Im Γ, which matters only for complex φ.)
- **Decision.** `ok = bool(H > 0 and Q1b > 0 and theta < fmpq(1, 4))` (certr.py l.104): Arb comparisons are True only if the
  whole ball is positive; H > 0 and Q1 > 0 give h1 > 1 and c_0 > 0.

## (3) Falsification runs (all under memcap 290 MB, one process at a time)
| Run | Result | Peak |
|---|---|---|
| indep7.py, BM's own data (r = 1, ϑ = 1/4, v = 1.26, κ = 3.18, P = 1−5.8x+6.4x², no shift), n = 12 | h1 = 1.0021778143 (BM: > 1.0002; GapCeil/GapLambda 1.00217781433) | < 45 MB |
| indep7.py, midpoint μ = c/2, r = 2, deg 8 (optmid_r2_d8_k3.25.json), n = 12/16/20 | 1.0177219 / 1.01560266 / **1.0156026615** | 85 MB |
| GapBar float (h1_float of that JSON) | 1.015602661541 | — |
| certmid.py, same JSON, κ = 3.25, D = 60 (Arb) | **1.01560266154 ± 4.5·10⁻¹³**, certified | 86 MB |
| certmid.py, d80 P reconstructed from the log (0 mismatches), κ = 3.324, D = 80 | Q0, Q1, h1 = 1.00056020750 ± 4.5·10⁻¹⁵, H, tail: identical to `GapBar/logs/certmid_r2_d80_k3.324.txt`; CERTIFIED True | 207 MB |
| certmid.py, same P, κ = 3.3266 (above the float ceiling 3.32658) | h1 = 0.998453611563, H < 0, CERTIFIED False (consistent) | 185 MB |
| identity at 2000 random (X1..X4, σ, c), 50 digits | max error 1.9·10⁻⁵⁰ | — |
| test P = the d80 Legendre series cut at degree 16 (not optimal), κ = 3.324: certmid D = 80 / indep7 n = 24, 28 | 0.450910418099 ± 1.9·10⁻¹³ / 0.4509104178, **0.4509104181** | 209 / 272 MB top (RSS 298 MB) |
| symmetry c_j(μ) = c_j(c − μ) (DERIVATION §4), indep7, deg 8, μ/c = 0.2 and 0.8, n = 16 | c0, c1 equal to 13 digits; h1 = 1.0065489444 both | 45 MB |
- The degree-8 agreement tests the whole Γ-array/Legendre chain (every array entry i, j < D enters through I^iφ), the twist and
  the r = 2 weights against a computation that uses only BM's TeX; the degree-16 test repeats this with certmid's exact D = 80,
  κ = 3.324, v settings. Degree 80 itself is out of reach of 7-D quadrature (the n = 28 run already peaked at 272 MB); it uses the same code.

## (4) What holds, minor findings, dependents
- **Holds as written:** Λ > 3.324 unconditionally; the identity; the certified run; the float ceiling 3.32658 is correctly labelled
  float. Nothing false found; nothing needs the unshifted fallback (Λ > 3.30, GapCeil, survives anyway).
- **Minor (no effect):** (a) BM l.395 divides by ie^{2ivθ(t)}, a leftover of the θ(t) version; should be ie^{ivtL} (|F'|² unaffected).
  (b) BM Lemma 4.2's proof (l.488) quotes Σ_{n≤y}d_r(n) = y(log y)^{r−1}/(r−1)! + O((log y)^{r−2}); the error should carry a factor y.
  (c) BM's numerics use ϑ = 1/4, outside "ϑ < 1/4"; our runs use 2499/10000. (d) BM's quote of BBLR omits the coefficient
  hypotheses; the original has them and they hold. (e) The d80 input JSON is not in GapBar/logs (only the log's A_n list); the
  exact reconstruction is `logs/optmid_r2_d80_k3.324_reconstructed.json` (orchestrator may copy it there).
  (f) DERIVATION §1 (complex P) is an extension of BM's proof, not the printed lemma: Lemma main has M, not M̄, on the right
  (l.366), so for complex P it is not the mean square; BBLR's bilinear form makes the extension routine. Not load-bearing for 3.324
  (real P), only for "complex P is equivalent" in RECORDS.md l.24 and the barrier certificates of DERIVATION §7.
- **Dependents** (rg gapbar|3.324): RECORDS.md l.24 (the record), NOGO.md l.630 (gapbar: NO BARRIER), TREE.md l.483 and
  Riemann/Tree/R447.lean (round 17 input check). None needs a change.

## Review
- reviewer subagent, inline (BM reduction l.164-210, Lemma main l.364-384 and its proof, BBLR Thm 1.1, derivative step):
  **CORRECT** — Lemma main permits complex α3, β3 = O(1/L); κπ/2 is a fixed constant inside the uniform range together with the
  Cauchy circles; real P gives the conjugation; the twisted coefficients have modulus |d_r(h)P[h]|; differentiate the exact uniform
  formula before replacing log y/L by ϑ; c_0, c_1 are real; needs c_1 > 0 (certified).

## Files
`indep7.py` (independent BM integral), `idmaj.py` (identity, majorant), `FULLPROOF.md` (the chain with line numbers),
`logs/`: indep7_r2_d8_k3.25.txt, certmid_d8_k3.25_D60.txt, rerun_certmid_r2_d80_k3.324.txt, above_ceiling_d80_k3.3266.txt,
idmaj.txt, optmid_r2_d80_k3.324_reconstructed.json, test_r2_d16_trunc.json, certmid_test_d16_k3.324.txt,
indep7_test_d16_k3.324.txt, indep7_symmetry_d8.txt.

RH relevance: N8

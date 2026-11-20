# FULLPROOF: Λ_f > 218/125 for primitive GL(2) L-functions (audit chain, gl2gapaudit, 20 Nov 2026)
Sources: TeX of arXiv 1410.7765 (BMMRTW), line numbers of .scratch/gl2gapaudit/tex/1410.7765.tex (single file, same as GapFamily's);
per-gap tool: ExploreShape CERT.md Lemma 4 (l.110–121), cert5.py (copy in GapFamily, identical to ExploreShape's: checked below).

## 0. Objects
- f primitive holomorphic or Maass cusp form, level q fixed (l.130–132). L(s,f) = Σ a_f(n)n^{−s}, a_f(1) = 1, critical line ½ (l.133–138).
- c_f := Res_{s=1} Σ|a_f(n)|²n^{−s} (l.230–234); c_f > 0 (a_f(1) = 1, positive coefficients, simple pole). Only c_f > 0 is used:
  c_f cancels in the criterion.
- Λ_f := limsup (γ̃(n+1) − γ̃(n))/(π/log(√q γ̃(n))), γ̃ = ordinates of zeros ON the line, repeated by multiplicity (l.156–162).
  Normaliser: N(T,f) = (T/π)log(√q T/2πe) + O(log 𝔮) (l.144–147, IK Thm 5.8 halved: IK count |γ| ≤ T, degree 2), so the mean
  spacing at height t is π/(log t + O_q(1)). On [T, 2T]: log(√q γ) = log T + O_q(1); and log(√q γ) ≥ log T since q ≥ 1.

## 1. Mean values (Thm mixed, l.239–247, proved l.400–478 from Thm shiftedmoment l.255–264)
- Thm shiftedmoment (hypotheses: primitive, GL(2)/Q, s = ½+it, T large, α, β ∈ C, |α|, |β| ≪ 1/log T):
  ∫_T^{2T} L(s+α,f)L(1−s+β,f̄)dt = ∫_T^{2T}{D(1+α+β) + (t/2π)^{−2(α+β)}D(1−α−β)}dt + O(T), D(s) := Σ|a_f(n)|²n^{−s}.
  [Audit] With conductor q the dual factor is (√q t/2π)^{−2(α+β)}; the paper drops √q: the difference is
  D(1−z)(t/2π)^{−2z}(q^{−z} − 1) = O(log q) for |z| ≪ 1/log T, i.e. O(T) after integration: absorbed (q fixed). Harmless.
- l.405–417: (t/2π)^{−2z} = T^{−2z}(1 + O(1/log T)), D(1 ± z) = ±c_f/z + O(1), z = α+β.
  [Audit] As printed, (c_f/z)·O(1/log T) is not O(1) for z → 0. True bound: (t/2πT)^{−2z} − 1 = O(|z|) on [T, 2T],
  so the error is O(1) uniformly, including z = 0 (removable singularity). Sloppy wording, statement right.
- F(z;T) = c_f T(1 − T^{−2z})/z = c_f T Σ_n (−1)^n 2^{n+1} z^n (log T)^{n+1}/(n+1)! (l.418–425): checked (series of (1 − e^{−2zℓ})/z).
- R(α,β;T) := I − F (l.426–431). [Audit] Printed with ∫_0^T; must be ∫_T^{2T} (typo; the final "dyadic intervals" sentence l.477
  is a leftover: the theorem is stated for [T, 2T], which is what the Cauchy step gives directly). Harmless.
- Cauchy (l.448–462): R = O(T) on |α|, |β| ≤ 2/log T ⇒ ∂_α^μ∂_β^ν R(0,0) = O(μ!ν! T(log T)^{μ+ν}) (two circles of radius
  1/log T; needs the O(T) of Thm shiftedmoment uniform in complex shifts of size ≤ 2/log T: "|α|, |β| ≪ 1/log T" as printed).
- ∂^{μ+ν}F at 0 = c_f T (−1)^{μ+ν}2^{μ+ν+1}(log T)^{μ+ν+1}/(μ+ν+1) (l.470–475): checked (n!·coefficient of z^n, n = μ+ν).
- Needed strength: R = o(T log T) uniformly on |α|, |β| ≤ 2/log T already gives ∂^μ∂^νR(0,0) = o(T L^{μ+ν+1}) (same Cauchy step),
  enough for §2–3; the printed O(T) is stronger than needed.
- Result (all fixed μ, ν ≥ 0; we use μ, ν ≤ 6): ∫_T^{2T}L^{(μ)}(½+it,f)L^{(ν)}(½−it,f̄)dt = (−1)^{μ+ν}2^{μ+ν+1}/(μ+ν+1)·c_f T L^{μ+ν+1}
  + O_{μ,ν}(T L^{μ+ν}), L := log T. L^{(ν)}(½−it,f̄) = conj L^{(ν)}(½+it,f) (coefficients conjugate). Uniformity in μ, ν: not needed.

## 2. Mean values of g^{(j)}, g(t) := e^{itL}L(½+it,f), L = log T (twist ρ = 1 of l.285–287)
- d/dt L(½+it,f) = iL'(½+it,f), so g^{(j)} = i^j e^{itL}Σ_k C(j,k)L^{j−k}L^{(k)}(½+it,f) and
  |g^{(j)}|² = Σ_{k,l} C(j,k)C(j,l)L^{2j−k−l}L^{(k)}(½+it,f)·conj L^{(l)}(½+it,f).
- §1 with (μ,ν) = (k,l): ∫_T^{2T}|g^{(j)}|² = c_f T L^{2j+1}Σ_{k,l}C(j,k)C(j,l)(−1)^{k+l}2^{k+l+1}/(k+l+1) + O_j(T L^{2j})
  = c_f T L^{2j+1}∫_0^2(1−s)^{2j}ds + O_j(TL^{2j}) = 2c_f T L^{2j+1}/(2j+1)·(1 + O_j(1/L)). Identity checked in exact
  rationals for j ≤ 6 (indep.py (c)). j = 0, 1 reproduce the paper's three integrals (l.354–362) and its 3ρ²−6ρ+4 at ρ = 1.
- Spectral reading: ∫|g^{(j)}|²/∫|g|² = L^{2j}E[U^{2j}], U uniform on [−1, 1]: L(½+it) has frequencies −log n spread
  uniformly (weight c_f d log n, Rankin–Selberg) over [−2ℓ(t), 0], ℓ(t) = log(√q t/2π) (the two halves of the approximate
  functional equation, lengths √q t/2π each); the twist e^{itL} centres them. Scale: ℓ(t) = L + O_q(1) on [T, 2T], so L = log T
  and log(√q T/2π) give the same leading moments; the zero density (1/π)log(√q t/2π) (§0) has the same scale L. In the zero
  normaliser and in the moments only L·(1 + O(1/L)) enters; both corrections vanish in the limit. (Numerical check: §5.)

## 3. Reduction (if Λ_f ≤ Λ₀ := 218/125 then contradiction)
- H(A) := −Σ_{j≤6}a_jA^{2j}/(2j+1); H(Λ₀) = 7.7047·10⁻⁴ > 0 (exact), H(1.7445) = 2.58·10⁻⁴ > 0, H(1.745) < 0 (indep.py (b)).
  Fix A ∈ (Λ₀, 1.7445] with H(A) > 0.
- Λ_f ≤ Λ₀ < A ⇒ ∃T₀: every gap between consecutive distinct critical zeros γ̃(n) ≥ T₀ is ≤ A·π/log(√q γ̃(n)) ≤ Aπ/L for
  γ̃(n) ∈ [T, 2T] (log(√q γ̃) ≥ log T, q ≥ 1). Multiplicity: repeated ordinates give gap 0; only distinct zeros are used
  (as l.292–297). Straddling gaps (reviewer's clarification, adopted): take Λ₀ < B < A; eventually every gap is ≤ Bπ/log(√q γ̃(n));
  the zero before T is T − o(1), so its gap is < Aπ/L too. Hence γ(1) − T, 2T − γ(N) ≤ Aπ/L (stronger than (zeroassumption) l.300–303).
- Rescale x = λt, λ := L/A: gaps in x are ≤ π. F_R(x) := Re g(x/λ), F_I := Im g(x/λ) are real-analytic, vanish at the
  zeros x_n = λγ(n). On each [x_n, x_{n+1}] apply (V₆) (§4) to f = F_R(x_n + ·), F_I(x_n + ·) ∈ H⁶, f(0) = f(g_n) = 0,
  g_n ≤ π. Sum over n: the six boundary terms φ_pq[f^{(p)}f^{(q)}] telescope (one global rescaling, F smooth: jets agree
  at shared endpoints), leaving φ_pq(F^{(p)}F^{(q)}(x_N) − F^{(p)}F^{(q)}(x_1)).
- ∫_{x_1}^{x_N}(F_R^{(j)})² + (F_I^{(j)})² dx = λ^{1−2j}∫_{γ(1)}^{γ(N)}|g^{(j)}|²dt. So, with §2,
  0 ≤ Σ_j a_jλ^{1−2j}[∫_T^{2T}|g^{(j)}|² − (end pieces)] + (end jets) = λ·2c_f T L·(−H(A) + O(1/L)) + E,
  since λ^{−2j}L^{2j+1} = A^{2j}L. [Audit] The REPORT's "−2c_f H(A)TL" omits the common factor λ = L/A > 0: sign unaffected.
- E (end pieces, end jets): |L^{(k)}(½+it,f)| ≪ (log T)^k·max_{|s−½−it| = 1/log T}|L(s,f)| ≪ T^{1/2−δ+ε} (Cauchy +
  t-aspect subconvexity near the line, l.330–337: Good, Meurman (level 1), Blomer–Harcos (any level, nebentypus), Michel–Venkatesh;
  any δ > 0 suffices). So |g^{(j)}| ≪ T^{1/2−δ+2ε}, end pieces ≪ (Aπ/L)T^{1−2δ+4ε}, end jets λ^{−p−q}|g^{(p)}g^{(q)}| ≪ T^{1−2δ+4ε}:
  E = O(T^{1−2δ+5ε}) = o(λ T L). Hence H(A) ≤ O(1/L) + o(1): contradiction for large T. So Λ_f > 218/125. ∎
- Without subconvexity (convexity bound only) E would be O(T^{1+ε}), not o(TL²): subconvexity is load-bearing but is a theorem
  for every f in scope (Blomer–Harcos covers all levels and Maass forms).

## 4. Per-gap inequality (V₆) and its certificate
- (V₆): ∀g ∈ (0, π], ∀ real f ∈ H⁶(0, g), f(0) = f(g) = 0: Σ_{j≤6}a_j‖f^{(j)}‖² + φ12[f'f''] + φ14[f'f''''] + φ32[f'''f''] + φ34[f'''f'''']
  + φ52[f⁽⁵⁾f''] + φ54[f⁽⁵⁾f''''] ≥ 0, [·] = [·]_0^g. Sufficient (Lemma 4, CERT.md l.110–121; re-derived, GapCertAudit FULLPROOF F1 +
  lemma46.py): p := Σa_jx^j > 0 on [1, ∞), S_odd(g) ⪰ 0, S_even(g) + C/2 ⪰ 0; reviewer condition deg p = 6.
- Code = audited code: run/cert5.py, exact6.py, exact5.py sha256 prefixes 0271dc508c4b5aa8, a2da4920ed882b7e, 940bc3a0840718f2 =
  GapCertAudit REPORT2 l.6–7; GapFamily copies = ExploreShape = git HEAD (cmp, git diff).
- Rerun on a copy (run/): cert JSON and criterion JSON byte-identical to GapFamily's (cmp). (1) p(1) = 0.0146345, no real root ≥ 1;
  a6 = 321577562/367354223721 = 8.75·10⁻⁴ > 0 (deg p = 6). (2) C11 = 5.71·10⁻⁶, det C = 5.82·10⁻¹¹ > 0. (3) leading coefficients of
  N11, N22, DET > 0 (4.9·10⁻⁶, 7.5·10⁻⁷, 1.3·10⁻¹⁵); X* = 2.2195379 (DET's largest possibly-real root).
- Cover: g < (π/√X*).lower() = 2.108718: x_k = (kπ/g)² > X* for all k ≥ 1 ⇒ every R(x_k) ≻ 0, C ≻ 0: termwise valid, down to g → 0
  (arbitrarily small gaps, which occur: zeros may be close). [2.108718, π.upper()]: 110 Arb cells, k ≤ 61, dropped modes
  x_k ≥ 62² > X* (PSD); min det S_odd ≥ 2.64·10⁻⁷ (g ≈ 3.1414), min det(S_even + C/2) ≥ 6.01·10⁻⁷ (g ≈ 2.1139). So the cover is
  all of (0, π] = every gap the reduction needs (§3). CERTIFIED.
- Data flow: cert5.py reads SRC = logs/gl2_J6_1.744m.json, rationalises a_j, φ_pq (limit_denominator 10¹²; D0 = E12 = E32 = 0) and
  writes exactly these to logs/cert5_gl2J6_1.744m.json; gl2crit.py reads 'a' from that file (Fractions), asserts a1 = 1.
  indep.py (a): own rationalisation of SRC = the certified a and φ (all 13 equal); max |float − rational| = 3.1·10⁻¹⁴.
  The φ_pq do not enter H (they telescope); the a_j of H are the certified ones.
- H exact (Fraction; 94/97-digit numerator/denominator) = 7.70468593084270331952·10⁻⁴; second formula H = −(1/Λ)∫_0^Λ p(v²)dv
  (uniform measure, mpmath Gauss–Legendre, 60 digits): equal to 4·10⁻⁶³.
- Falsification of (V₆), independent of Lemma 4 (indep.py (d)): the form itself on f(x) = h(x/g), h = u^{m+1}(1 − u), m < M,
  exact rational matrices; exact LDLᵀ: all pivots > 0 at g = 0.01, 0.1, 0.5, 1, 1.5, 2, 2.108718, 2.5, 3, 3.1, 3.14, 3.141592,
  3.1415926; min form/‖f‖² (M = 22) = 9.8·10⁴, 922, 34.6, 6.96, 2.26, 1.07, 0.916, 0.478, 0.0860, 0.0246, 1.39·10⁻³,
  4.77·10⁻⁴, 4.77·10⁻⁴ (M = 8: 4.92·10⁻⁴ at g ≈ π; converging from above, no negative direction). For f = sin(kπx/g) all six
  boundary terms vanish (f'', f'''' vanish at the ends) and the form is (g/2)p((kπ/g)²): check (1).

## 5. Numerical test of the mean-value input (Δ, level 1, weight 12; meanval.gp, PARI 2.19 lfun, logs/meanval_T{200,400}.log)
- L(½+it, Δ) = PARI lfun(lfunetaquo([1,24]), 6+it) on t ∈ [T, 2T], step 1/50 (trapezoid; spectrally accurate); L', L'' by 8th-order
  central differences (6th-order difference: relative change ≤ 3·10⁻⁷); Im ∫L'L̄ = −½[|L|²] (boundary) reproduced (−4.2528 vs −4.2443).
- c_Δ = Res D = 0.384084013, Laurent D(1+z) = c/z + d + ez + fz², d = 0.7960286, e = 0.051777, f = −0.00764, from smoothed
  sums Σ|a(n)|²n^{−1}(log n)^k e^{−n/X}, n ≤ 1.6·10⁵, X = 1000, 2000, 4000, with the Γ(w)-pole 1/X corrections (fits agree to 10⁻⁷).
- Recipe = Thm shiftedmoment's main term differentiated exactly at the local scale ℓ(t) = log(t/2π):
  ∂^n_z[D(1+z) + e^{−2zℓ}D(1−z)]_{z=0} = c(−1)^n(2ℓ)^{n+1}/(n+1) + d(1+δ_{n0})(−2ℓ)^n − en(−2ℓ)^{n−1} + f(2δ_{n2} + n(n−1)(−2ℓ)^{n−2}).
- I_μν := ∫_T^{2T}L^{(μ)}conj L^{(ν)}, μ ≤ ν ≤ 2. I/recipe: T = 200: 0.977, 0.975, 0.981, 0.981, 0.986, 0.988; T = 400: 1.0087, 1.0080,
  1.0057, 1.0049, 1.0031, 1.0017 (fluctuation of the second moment, size ~ Good's (T log T)^{2/3}; sign changes with T).
  I/(Thm mixed's leading term with log T): T = 400: 1.113, 0.844, 0.744, 0.744, 0.644, 0.552 (T = 200: 1.091 … 0.486): the expected
  O((μ+ν+1)/log T) lower-order terms; all of them are the recipe's (ℓ = log T − log 2π − O(1) and the d-, e-terms).
- g = e^{it log T}L: ∫|g^{(j)}|² / (2c T L^{2j+1}/(2j+1)), j = 0, 1, 2: T = 200: 1.091, 1.269, 1.578; T = 400: 1.113, 1.251, 1.495
  (→ 1 like 1/log T: the twist log T sits log 2π above the spectral centre ℓ(t)); against the recipe: 1.0087, 1.0026, 1.0019.
- Verdict of the test: the normalisation (c_f = Res Σ|a(n)|²n^{−s}, factor (−1)^{μ+ν}2^{μ+ν+1}/(μ+ν+1), scale 2ℓ(t) with
  ℓ(t) = log T + O(1)) is right; the leading asymptotics converge slowly (1/log T), which the limit argument does not need.
- Zero count (zeros.gp, logs/zeros_T200.log): Z(t) = L(6+it)e^{i(arg Γ(6+it) − t log 2π)} is real (|Im Z| ≤ 4·10⁻³⁵); 245 sign
  changes on [200, 400] vs N(2T) − N(T) = (T/π)log(2T/πe) = 244.89: density (1/π)log(t/2π), as §0. Largest gap there: 2.886 in
  Λ_f's units π/log γ (1.981 in local units π/log(γ/2π)); the mean gap in Λ_f's units is 1.48 at this height (→ 1 like 1/log γ).
- The spectral half-width ℓ(t) = log(√q t/2π) is the same function as π × zero density (§0): in local units the moments are
  exactly the uniform ones; using L = log T in both places is a uniform (1 + O(1/L)) approximation.

## 6. Used as printed (not re-refereed line by line)
- Thm shiftedmoment's proof (§7–8 of 1410.7765, l.1177–1544: Ramachandra-type approximate functional equation, J_i error terms,
  Lemmas weightsum, truncatedsum, stirling, GeneralMV). Checked: its main-term computation (S_1, S_2, l.1283–1360) holds for complex
  shifts ≪ 1/log T (weights n^{−2Re α} = O(1); the X-dependent terms c_fTΓ(−z)(X/2)^{−z} + (c_fX^z/z)∫(t/2π)^{−2z} are O(T)
  uniformly, z = α + β, incl. z → 0). Its Stirling factor (t/2π)^{−2(α+β)} lacks √q (O(T) for fixed q, §1). Its error section
  (l.1372–1430) treats one "representative" term and says the other J_i follow similarly (a sketch; published: JMAA 2015,
  DOI 10.1016/j.jmaa.2015.04.007, = arXiv v4, the TeX used here).
- Lemma 4 and cert5.py: audited in GapCertAudit REPORT2 (same files, sha256); here rerun and falsification-tested (§4).
- Subconvexity in the t-aspect for all levels, holomorphic and Maass (Blomer–Harcos; Michel–Venkatesh), as cited l.330–337.

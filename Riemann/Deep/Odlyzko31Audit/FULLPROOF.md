# odlaudit FULLPROOF: the chain behind discriminant-bounds-strip, with source line numbers

Sources: D = Riemann/Deep/Odlyzko31/DERIVATION.md (line numbers), C = certify.py, L = Challenge/DedekindZetaChallenge.lean (AxiomMath/DedekindZetaZeroFree@09afbb4, fetched 19 Nov 2026), M = Mathlib@0493595e Mathlib/NumberTheory/NumberField/DedekindZeta.lean.

## Step 0. Input H(15/32) (D l.8-9)
- L: `theorem main (K) [Field K] [NumberField K] : ∃ f : ℂ → ℂ, DifferentiableOn ℂ f {1}ᶜ ∧ (∀ s, 1 < s.re → f s = NumberField.dedekindZeta K s) ∧ ∀ s, 31 / 32 < s.re → s ≠ 1 → f s ≠ 0`.
- M l.49-50: dedekindZeta K s = LSeries (n ↦ #{I : absNorm I = n}) s: the Dirichlet series, junk 0 where it diverges (Re s ≤ 1). The existential f avoids the junk: f is holomorphic on the open connected set ℂ∖{1} and equals ζ_K on Re s > 1, so by the identity theorem f is THE continuation of ζ_K on ℂ∖{1}; hence ζ_K(s) ≠ 0 for Re s > 31/32 (incl. 31/32 < Re s ≤ 1). Not vacuous (f exists: the classical continuation). 31/32 elaborates in ℝ (s.re : ℝ).
- Reflection: Λ_K(s) = Λ_K(1 − s) (classical, not in Lean) maps a nontrivial zero ρ to 1 − ρ; so 1/32 ≤ β ≤ 31/32, i.e. |β − 1/2| ≤ 15/32. Zeros of Λ_K = nontrivial zeros of ζ_K (Γ-factors have no zeros; ζ_K ≠ 0 on Re s ≥ 1 and by FE none with Re s < 0 except trivial ones cancelled by Γ poles). CHECKED.
- Proof status of L (HPAudit/REPORT.md l.79-110, not re-run here): Comparator match, standard 3 axioms, three kernels in third-party CI.

## Step 1. Explicit formula (P) (D l.13-33; (P) at l.31, E at l.32)
Re-derived by hand (this audit):
- (1.1) residues of Φ·Λ'/Λ on [1−c, c]×[−T, T]: Σ_ρ Φ(ρ) − Φ(0) − Φ(1) = (1/2πi)[∫_(c) − ∫_(1−c)]; z ↦ 1 − z with Φ(1 − z) = Φ(z) (F even) and Λ'/Λ(1 − z) = −Λ'/Λ(z) turns ∫_(1−c) into −∫_(c). CHECKED.
- constant (1/2)log|d_K| → F(0)·log|d_K| = log|d_K| (Fourier inversion on Re z = 1/2). CHECKED.
- ζ_K'/ζ_K → −2 Σ_a Λ_K(a) Na^{−1/2} F(log Na) ((1/2πi)∫_(c) Φ(z)e^{−uz}dz = e^{−u/2}F(u)). CHECKED.
- Γ pieces: Γ_C'/Γ_C = −log 2π + ψ, Γ_R'/Γ_R(z) = −(1/2)log π + (1/2)ψ(z/2); with ψ(z) = −γ + ∫_0^∞ (e^{−x} − e^{−zx})/(1 − e^{−x}) dx and (1/π)∫F̂(t)cos(tx)dt = 2F(x): (1/π)∫F̂ Re ψ(1/2+it) = 2[ψ(1/2) + I_c], (1/π)∫F̂ Re ψ(1/4+it/2) = 2[ψ(1/4) + I_c + I_r] (2e^{−x/2}/(1 − e^{−2x}) = e^{x/2}/sinh x = 1/(2 sinh(x/2)) + 1/(2 cosh(x/2))). CHECKED.
- Collect: r2-part 2r2(γ + log 8π − I_c) (ψ(1/2) = −γ − 2 log 2), r1-part r1(γ + log 8π − I_c) + r1(π/2 − I_r) (ψ(1/4) = −γ − π/2 − 3 log 2). Total n(γ + log 8π − I_c) + r1(π/2 − I_r): (P) as in D l.31. CHECKED. Independent numerical identity: see Step 6.
- Numerical identity, independent PARI/GP code (identity.gp, identity2.gp; logs/identity*.log): (P) with F = H(x/1.5)·sech(15x/32), zeros from lfunzeros, primes from idealprimedec, I_c, I_r, E by intnum with numerically integrated tails: residual 9.2e-13 (ℚ, 138 zeros), 5.4e-12 (ℚ(√5), r1 = 2), 2.5e-11 (cubic d = −23, r1 = r2 = 1), 5.2e-13 (ℚ(√−3), r2 = 1, zeros of ζ ∪ L(χ₋₃) up to 400). Second kernel G(x/2)·sech(w_ab x) on ℚ(√−3): 1.0e-8 (tail-limited, T/2→T change 6e-8). Control: dropping the r1 term moves ℚ(√5) by 0.888. So every Γ-term (Γ_R and Γ_C), the pole term E and the prime normalisation are right. Artifact noted: lfunzeros on the product ζ_K of ℚ(√−3) misses close pairs from the two factors (residual 7.5e-9 / 1.8e-6 at T = 250), cured by the factorwise zeros.

## Step 2. Positivity under H(w) (D l.35-42)
- Lemma A (D l.36-37): ∫ e^{zx}/cosh x dx = π/cos(πz/2) (|Re z| < 1); sum over z = (±a + it)/w gives (π/2w)·[1/cos(α+iτ) + 1/cos(−α+iτ)] = (π/w)cos α cosh τ/(cosh²τ − sin²α) = (2π/w)cos α cosh τ/(cosh 2τ + cos 2α), α = πa/2w, τ = πt/2w. Re-derived by hand; > 0 for |α| < π/2. Numerically (consts_mp.py, 256-bit quad/quadosc) at 5 (a, w, t) incl. a = 0.46 vs w = 15/32: diff ≤ 1.4e-76. CHECKED.
- Lemma B (D l.38): G, H are the autocorrelations of cos(πt/2) and 1 + cos πt on [−1, 1] with norms 1 and 3: re-derived by hand (both closed forms reproduced), and FT(S)(t) = |ĥ(t)|²/‖h‖² checked at t = 0, 1.7, 5.3, 12 (diff ≤ 7e-80). h ≥ 0 ⇒ S ≥ 0, S ≤ S(0) = 1, Ŝ ≥ 0. CHECKED.
- Prop C (D l.39-41): Re Φ(1/2 + a + iγ) = ∫ F(x) e^{ax} cos γx dx = ∫ F(x) cosh(ax) cos γx dx (F even); F·cosh(a·) = S(·/b) · cosh(a·)/cosh(w·); its FT is (1/2π)(b Ŝ(b·)) * (Lemma A kernel) ≥ 0 for |a| < w, and = b Ŝ(bγ) ≥ 0 for |a| = w. Real zeros (γ = 0) included. Conjugation symmetry of the zero multiset makes Σ Φ(ρ) = Σ Re Φ(ρ); absolute convergence from Φ(c+it) ≪ (1+|t|)^{-2} (F ∈ C², compact support) and N(T+1) − N(T) ≪ log(|d|T^n). F ≥ 0 ⇒ prime terms ≥ 0. Hence (T): log|d| ≥ r1 log A + 2r2 log B − E. CHECKED.
- Signature monotonicity (D l.42, Remark 2): I_r(F) ≤ ∫_0^∞ dx/(2cosh(x/2)) = π/2 (needs only F ≥ 0) ⇒ log A ≥ log B ⇒ r1 = 0 value bounds every signature; E > 0 ⇒ (log B − E/n), (log A − E/n) increase in n ⇒ a value certified at n holds for all larger n. CHECKED.

## Step 3. Extremality and limit (D l.44-49)
- Prop D: φ = F·cosh(w·) continuous, integrable (F cosh(x/2) integrable, w ≤ ½), φ̂ ≥ 0 (admissibility at a = w) ⇒ φ̂ ∈ L¹ and φ(x) = (1/2π)∫φ̂ e^{−itx} ⇒ |φ(x)| ≤ φ(0) = 1 ⇒ F ≤ sech(wx) ⇒ I_c(F) ≥ I_c(w), I_r(F) ≥ I_r(w); E > 0 ⇒ (T)/n < log(A_r^{r1/n} A_c^{2r2/n}). CHECKED.
- Converse: F_b = H(x/b)sech(wx) admissible (Step 2); 1 − F_b ≤ (1 − sech wx) + min(Cx²/b², 1) with C = sup(1 − H(t))/t² (H = 1 − (π²/6)t² + O(t⁴)); dominated convergence ⇒ I(F_b) → I(w); E(F_b)/n → 0 for fixed b. Order: n → ∞ first, then b → ∞ (or slow diagonal b(n)). CHECKED (reviewer agrees).
- Scope of optimality: within the admissible class only (zeros may sit anywhere in the strip at any height). The page says so.

## Step 4. Constants (D l.50-55; certify.py limit_consts)
- consts_mp.py (mpmath, 256-bit, [0, ∞) split at powers of 2, own code): w = ½: I_c = log 2, I_r = π/2 − 1 to 15 digits, A_c = 22.3816160955 = 4πe^γ, A_r = 60.8395403239 = 4πe^{1+γ}; w = 15/32: I_c = 0.649168544494663, I_r = 0.538769195117049, A_c = 23.3878942309, A_r = 65.6439656548; w_ab: 26.9662525665, 84.2999165937; Liu (arXiv): 26.9032357855, 83.949939377; w = 0: 44.763232191, 215.332515953 = 8πe^γ, 8πe^{γ+π/2}. All agree with certify.log. CHECKED.

## Step 5. Thresholds and table (D l.66-72; certify.py)
- 256-bit recomputation (logs/consts_mp_special.log): n = 2732, r1 = 0, H, b = 20.412: 22.3816960154 > 22.38169 > 4πe^γ = 22.3816160955; n = 1618, r1 = n, b = 20.742: 60.8402198645 > 60.84021 > 60.8395403239; n = 1581, r1 = 0, b = 17.398: 22.0002274248 > 22; n = 1223, r1 = n, b = 19.129: 60.0003991246 > 60; abelian lines 482/285/35448 also pass. Extension to all larger n and all signatures: Step 2 monotonicity. CHECKED.
- Scans (make_kernels.py crossover): start at n = 2, double, then bisect; the certified n equals the scan's n0. CHECKED.
- 40 random table entries (seed 20261119, logs/consts_mp.log): all exceed their 3-decimal certified values; logA, logB, E agree with certified.json to ≤ 3.9e-15 relative. CHECKED.

## Step 6. Sharpness control (posit.py, logs/posit.log)
- For the certified crossover kernel F = H(x/20.412)·sech(15x/32), min over γ ∈ [0, 3] of Re Φ(½ + a + iγ): ≥ 0 for a = 0, 0.3, 0.46 and −2.8e-15 (rounding) at a = 15/32; −0.25 at a = 0.475, −1.1 at 0.49, −1.9 at 0.5. So the bound genuinely needs β ≤ 31/32 for every zero: the Cureton–Thorner input is load-bearing at exactly its stated width.


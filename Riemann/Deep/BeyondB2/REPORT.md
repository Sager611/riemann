IDEA: modified system gives B = B(σ) < B_2 (exact closed form; B(0.0075) = 0.8725612664, dB/dσ(0) = −0.1078) if on every zero-detected row class 𝓑 of count U^R the signed row sum of the correction, Σ_{u∈𝓑} c_u L(w,χ_u), is U^{−σR} below the paper's absolute bound max_𝓑|c_u|·(#𝓑)^{3/4}(Σ_𝓑|L(w,χ_u)|⁴)^{1/4}; and B = B_2 − λ exactly if the direct Cauchy pairing (NF3 Prop direct-estimate) saves Z^{−λ}.

Status: neither inequality is proved, and none of the scheme's inputs can prove them. Every input is blind to the phase of an individual row. The one provable signed tool, a large sieve in the Möbius variable, loses Z^{0.0333} at the contact. Finite models satisfy the dual inequality generically (σ_eff ≈ 0.10–0.19, more than 10× the 0.008 needed). It fails on phase-aligned subclasses that carry the same magnitude data.

# BeyondB2: below the contact B_2 = (777 − √1081)/852 by cross-row cancellation (node N7, 2026-11-05)

Sources: Liu, Hexagon 2611.00082v1 (`number_field_new_v2-1.tex`, "l."); its input NF3 = Liu, arXiv 2610.26835v2 (`number_fields_v4.tex`, "NF3 l."). The scalar model is the audited one in [../Hexagon2](../Hexagon2/REPORT.md) (99 checks). All new computations are in `verify.py` (17 checks, 0 failed, 6 s, 115 MB), `infimum_scan.py` (250 s, 90 MB) and `finite_model.py` (16 s, 270 MB). Each runs as one single-threaded process.

## Short answer
1. **Where sign information is lost (§1).** Two places matter at the contact. Everything else is positive by nature (moments, counts) or inactive.
   - **Dual side:** the triangle inequality over the frequency rows u of a zero-detected class. NF3 (v:cmp:11), l.9389; 2611 Lemma sc:reference, l.3497–3501 and l.3534–3538. It drops the phases of 1/L(s,ηχ̄_u), of the selected prime product, of the Euler correction H_u and of the root number.
   - **Direct side:** the Cauchy–Schwarz pairing over the theta rows k. NF3 Prop direct-estimate, l.8520. It drops the correlation between the theta sum B_k and the Gauss-sum Gram factor Ã_k. The constant 11/3 in ℓ = 11/3 − 4B comes from this pairing.
2. **What a signed estimate would buy (§2), computed exactly.**

| lever | hypothesis | new infimum | to reach 0.873 | to reach 0.87 | to reach 13/15 |
|---|---|---|---|---|---|
| dual: signed class sum | saving U^{−σR} on each class | B(σ), closed form (§2.1); dB/dσ = −0.10784 | σ = 0.003518 | σ ≈ 0.0296 | σ ≈ 0.06 (then 13/15 binds) |
| direct: signed pairing | saving Z^{−λ} | B_2 − λ exactly; the 13/15 floor also moves to 13/15 − λ | λ = 0.000382 | λ = 0.003382 | λ = 0.006715 |
| signed large sieve (provable type) | none (conjectural diagonal sextic large sieve) | B_2: no gain; it loses Z^{0.0333} at the contact | — | — | — |
| pairing rows (unit classes, Galois conjugates) | — | B_2: paired rows have different denominators; conjugate pairs give only 2 Re | — | — | — |

   Per unit of exponent saved at the contact, the direct lever is 4.7 times stronger than the dual one: dB = −1 against −1/(3+4δ_2) = −0.2126.
3. **Tests (§3).** The finite model uses exact L-values of every Dirichlet character modulo a prime q ≤ 786433.
   - Classes selected by size cancel at a power rate: σ_eff = 0.10–0.19 for the largest corrections and for the contact-like class.
   - Two kinds of class show no power saving: classes selected by prime-sum amplitude alone (σ_eff ≈ −0.03), and large random classes, which sit on a family-mean bias of about 0.1.
   - Phase-aligned quarter-arc subsets keep r ≈ 0.90 with about 28% of the rows.
4. **What a proof needs (§4).** One new arithmetic input.
   - Dual side: phase information about 1/L(s,ηχ̄_u) on the rows where zeros are detected, which amounts to excluding phase conspiracies among hypothetical zeros. No known method gives this.
   - Direct side: a bilinear or spectral bound for a concrete pairing of theta coefficients with the Gauss-sum Gram factor, beating Cauchy–Schwarz by any power. This is the cleaner target, because it concerns an actual arithmetic sum rather than hypothetical zeros. The cost: such a bound combined with the existing correction bounds is equivalent to the new zero-free region. So it has to come from cancellation that does not presuppose zero-free information, for example the spectral theory of the sextic theta function.

## 1. Where the scheme takes absolute values across rows
"Row" here means a frequency u of the dual side, written Σ_u in NF3 (v:cmp:6) at l.8775: Σ_u^{(6)} q_u^{−z} ζ_F^S(6z) L(w,χ_u) H_{η,u,Z}(s,w,z) / L(s,ηχ̄_u). It can also mean a theta row k of the direct side, an amplitude class or a parameter dyad ρ.

| # | step (label, line) | across | sign information discarded | active at the contact? |
|---|---|---|---|---|
| S1 | NF3 (v:cmp:11), l.9386–9391: Σ_u \|L(w,χ_u)H_u^{piece}\| | rows u of one class | phase of c_u = q_u^{−z}ζ(6z)H_u/L(s,ηχ̄_u), of the selected product Q_sel(u) and of L(w,χ_u) | **yes**: it produces the class first moment |
| S2 | 2611 Lemma sc:reference, l.3497–3501, and γ_2, l.3534–3538: Σ_{u∈𝓑}\|S_m\| ≤ (#𝓑)^{3/4}(Σ\|S_m\|⁴)^{1/4} | rows u ∈ 𝓑 | none beyond S1 (Hölder on moduli); this is where (3R+1)/4 enters E(R) | **yes**: E(R_d), G_N |
| S3 | NF3 l.10406–10409 "root numbers are removed by absolute values"; 2611 l.3352–3354 | rows u | the root number ε_u (a sextic Gauss-sum phase) | inside S1 |
| S4 | NF3 l.10369–10372: triangle inequality over the S-part 𝔞 | finite S-divisors | phases ψ̄_u*(𝔞) | no (bounded factor) |
| S5 | NF3 l.9306–9310: local six-row table, "every phase coefficient has modulus one … bound all its monomials individually" | local factors at p \| u | local Gauss phases a, b_0 | no (local errors only) |
| S6 | 2611 l.3855–3856 and l.4396–4398: the whole holomorphic numerator is moved *before* absolute values | — | nothing: this is where signs are still kept (both AFE halves, deletions) | it marks the boundary of S1 |
| S7 | NF3 Prop direct-estimate l.8504–8538: "Cauchy for the exact pairing Q'^{−1/2}ΣÃB", with (v:cmp:2) l.8247 and (v:cmp:3) l.8252 | theta rows [k], q_k ≪ X'Y' | correlation of B_k^{(J)} (theta sum, target η) with Ã_k (Gauss sums g_s(−k), s ~ Y') | **yes**: B_low = 11/12 − ℓ/4 = B, the source of the 11/3 |
| S8 | NF3 l.8372–8375: Hilbert-space triangle inequality inside ‖Ã‖² | columns s | profile phases (the paper says these "cancel in each positive square") | no (a norm bound) |
| S9 | 2611 l.3303: detector "Mellin shift and dyadic separation", one large dyad chosen | detector dyads (r,m) | the signed identity Σ_{(r,m)} M_rS_m(ρ) = −1 | no (the count only needs one large piece) |
| S10 | 2611 l.3438 (Markov), (sc:sobolev) l.3325 | rows u | none: positive moments Σ\|M_rS_mQ\|² | no |
| S11 | 2611 (sc:lower-rows) l.3727, (sc:row-slope) l.3736: dyads ρ summed by the triangle inequality | parameter dyads ρ | cross-scale cancellation | no: the top dyad ρ = h dominates (slope > 0.11), so no cancellation is possible |
| S12 | 2611 l.4181–4184 (unequal factors): Cauchy–Schwarz on the counted parameters | rows u | phases of the mixed product | it reproduces the contact exactly (neutrality), so moduli alone cannot help |

Liu marks the same boundary: 2611 l.4256–4260 ("neither supplies nor rules out a phase-sensitive mixed inverse–ordinary estimate") and NF3 (sec:v-detail-37), l.11558–11563 ("Additional mean cancellation in the coupled numerator/correction, a sharper Gram-matrix estimate for the original sum … changes an input"). S1/S2 is the first of these; S7 is the second.

## 2. Modified systems

### 2.1 Dual lever: signed class sum, exact infimum
**Hypothesis H(σ).** Let 𝓑 be a row class of count U^R, fix (s,w,z) on the contours and write N_u = L(w,χ_u). Assume the signed row sum satisfies |Σ_{u∈𝓑} c_u N_u| ≤ U^{−σR+ε} times the absolute bound the paper proves for Σ_{u∈𝓑}|c_u N_u|. That bound is max_𝓑|c_u|·(#𝓑)^{3/4}(Σ_𝓑|N_u|⁴)^{1/4}, the Hölder step with M2 in 2611 Lemma sc:reference and NF3 (eq:v-refined-block) l.10350. Here σ = 1/2 is square-root cancellation over the class, and σ = 0 is the paper's bound.

**Effect on the system.** The class first moment (3R+1)/4 becomes (3R+1)/4 − σR. In 2611 (sc:reference-functions) and (sc:numerator-tests) this gives:
- E(R) → E(R) − hσR on every class and every numerator band. The gains G_N, G_L, G_R and the short-sum term D_sh are relative to the reference, so they are unchanged.
- The row slope becomes S(R) → S(R) − σR.
- The floor bins (trivial count R = 1) gain hσ, and the small-amplitude bins gain at least (7/8)hσ.

**KKT computation (exact, sympy; `verify.py` A1–A6).** The active margin becomes g_σ = G_N − E_σ(R_d) = g + hσR_d. At x = 1/2 it is affine in B and in b, with no B·b term.
- Its b-coefficient is [(33−36σ)δ² − (61−98σ)δ + (20−60σ)] / (24(4δ−5)).
- The contact amplitude is the smaller root: δ*(σ) = (61 − 98σ − r)/(2(33 − 36σ)), with r = √(964σ² − 1156σ + 1081).
- At δ* the margin is free of b. Setting it to 0 gives

  **B(σ) = (103528σ² + 7316σr + 51710σ − 3823r − 51257) / (12(8720σ² + 640σr + 4912σ − 356r − 5179)).**

- B(0) = 259/284 − √1081/852 = B_2 exactly, and δ*(0) = δ_2.
- dB/dσ(0) = −h_2R_2/(3 + 4δ_2) = −0.1078420638.

| σ | 0 | 0.001 | 0.0035 | 0.005 | 0.0075 | σ_f = 0.007959 |
|---|---|---|---|---|---|---|
| B(σ) | 0.8733819666 | 0.8732739158 | 0.8730019505 | 0.8728374989 | 0.8725612664 | 0.8725102987 |

**Validity.**
- For σ ≤ 0.0075 the contact (δ*(σ), 1/2) is the only active margin. All other entries of 2611 Prop sc:rational-certificate and Lemma sc:outside stay positive on 301² grids (A7). In that range B(σ) is the exact infimum.
- At σ_f = 0.007959 the floor requirement (sc:floor-small-rows) starts to bind (A8).
- `infimum_scan.py` re-optimizes b with all the structural constraints. The paper's numerator-band constants are re-choosable and are not imposed; the principal box is relaxed to Re s > 0.67 (§2.2). The scan stays within 1.2·10⁻⁵ of B(σ) up to σ = 0.03 (0.8722824, 0.8711429, 0.8699600 at σ = 0.01, 0.02, 0.03). Beyond that it reaches 0.8676503 at σ = 0.05 and saturates at **13/15** for σ ≳ 0.06, because the branch ℓ ≤ 1/5 of the direct estimate binds.
- So the dual lever alone cannot pass the 13/15 floor.

### 2.2 Direct lever: signed Cauchy pairing, B_2 − λ exactly
- The direct exponent is l_x/2 + b/12 = C(B_low), with B_low = 11/12 − ℓ/4 for every b (B2).
- Every margin of the correction system is a function of (ℓ, b, δ, x) only (B1). B itself appears only through the direct relation ℓ = 11/3 − 4B and through the principal box B > 0.87.
- Hence a bound |I_mod| ≪ Z^{l_x/2+b/12−λ}, uniform near (ℓ_2, b_2), gives B_low = 11/12 − ℓ/4 − λ under the same margins. The infimum is **B_2 − λ exactly**. Uniform means: for the empty subset and for every selected subset of length d < 8λ/3; longer subsets already have the slack 3d/8 of NF3 l.8522–8527.
- This is R310's dB/dC = 1/4 for the constant C = 11/3, now located: C is the Cauchy–Schwarz pairing S7.
- The principal box is not a floor. With the same corner (Re w, Re z) = (0.95, 0.165), all 12 decay exponents of (vtwo:principal-polynomial) exceed 1, and all five selected-factor exponents are positive, for Re s ≥ 0.67. The binding monomial is R, with exponent 6 Re s − 3.01 (B3).
- The 13/15 floor of R109 is the ℓ ≤ 1/5 branch of the same direct estimate (the (d+5ℓ−1)_+/8 term), so a uniform saving λ moves it to 13/15 − λ.
- What has to be beaten: the d = 0 (empty-subset) term, i.e. the pairing of the selected-prime-weighted theta sum with Ã. Nonempty subsets already have slack of at least 3d/8.

### 2.3 Signed large sieve (the provable kind of signed estimate): no gain
Write row u's s-integral as the Möbius sum F_u = Σ_{m~Z} μ(m)ηχ̄_u(m)w(m/Z). The scheme bounds |F_u| ≤ Z^{a+ε} pointwise. A diagonal sextic large sieve would give Σ_u|F_u|² ≪ Z^{2+ε}, which is conjectural in this range: (UZ)^{2/3} > Z, and the selected large sieve of 2611 Thm op:selected needs prime-amplitude slope > 1/2, while the contact has δ_2 = 0.426.
- **LS1** gives |Σ_{u∈𝓑} b_uF_u| ≤ Z‖b‖_2. Compared with the absolute route Z^a‖b‖_1 it loses Z^{1−a−hR/2}. At the contact this is **1 − a_2 − h_2R_2/2 = −251/2343 + 21635√1081/5065566 = 0.0332965** (C1). It is positive on the whole of [1/8,5/6]×[0,1/2] (min 0.0072, C3).
- No Hölder interpolation between the pointwise bound and the large sieve helps: the exponent is linear in 1/p with slope 2Δ_1 > 0.
- **LS2** also replaces the pointwise prime bound by M2 for S_mQ_full (allowed, since m_0 + 2ℓ/h = 0.855 ≤ 1). It loses Z^{0.0782397} (C2).
- In H(σ) units, LS1 corresponds to σ = −0.0656. The system that takes min(absolute, LS) per class keeps the contact, so its infimum is B_2.
- Mean-square averaging over the Mellin scale Z is the same estimate (Plancherel ↔ large sieve in m). It cannot resolve zeros whose ordinates lie within 1/log Z of each other, and #𝓑 = U^{0.63} ≫ H forces such clusters.

### 2.4 Pairing rows, positivity
- **Unit classes εu.** χ_{εu} = χ_εχ_u, so the paired rows have the different denominators L(s,ηχ̄_εχ̄_u), and the orthogonality over ε is already used in (v:cmp:6).
- **Galois conjugates c(u)**, when F is CM. For a c-invariant target, corr_{c(u)} = conj(corr_u). The pair gives 2 Re corr_u: a factor 2, no power saving.
- **Positivity, η = 1.** After reflection the row ratio is L(1−w,χ̄_u)/L(s,χ̄_u), the same L-function. Zeros would cancel only on the diagonal s = 1 − w, but the s- and w-integrals are independent. Coupling them is a different operator (NF3 l.11566–11569, "combining theta averages with cancellation"). The bootstrap also needs every target, because β_* is a supremum over the whole family.

## 3. Finite-model test of H(σ) (`finite_model.py`; summary also in `verify.py` D)
**Model.**
- Rows: all nonprincipal Dirichlet characters modulo a prime q, as proxies for χ_u at U = q.
- Exact L-values: L(s,χ_j) = q^{−s}·FFT_k[ζ(s, g^k/q)], with the Hurwitz ζ computed by Euler–Maclaurin. This matches `mpmath.dirichlet` to 10⁻¹³ (D0).
- Per row:
  - T = L(1/2−θ_0, χ)·Q(χ)·c(χ);
  - Q = P^{−1/2}Σ_{P<p≤2P} χ̄(p)w(p/P), with P = q^{0.21};
  - c = 1/L(0.75, ηχ̄), with η = χ_7 (generic) or η = 1.
- Classes of size K = q^{0.6266} = q^{R_2}.
- Statistic: r = |Σ_𝓑T|/Σ_𝓑|T|, with σ_eff = −ln r/ln K fitted over K = q^{0.4…0.8}. H(σ) asks for r ≤ K^{−σ}; random phases would give 1/2.

| class (q = 786433, η = χ_7, K = 4946) | r | fitted σ_eff |
|---|---|---|
| top \|T\| (largest corrections) | 0.302 | 0.163 |
| contact-like: top \|c\|·\|Q\| | 0.179 | 0.096 |
| zero-like: top \|c\| (smallest \|L(0.75)\|) | 0.150 | 0.112 |
| amplitude only: top \|Q\| | 0.145 | −0.028 (no power saving) |
| random | 0.105 | 0.030 (saturates at family-mean ratio 0.104) |
| phase-aligned quarter-arc subset of the contact-like class (1396 rows) | 0.905 | → 0 |

- q = 65537 and 163841 give the same picture: top |T| σ_eff = 0.186 and 0.180; contact-like 0.171 and 0.113.
- With η = 1 the zero-like and contact-like classes cancel strongly (σ_eff up to 0.8): numerator and denominator are the same L-function.

**Reading.**
1. Wherever the corrections are large, they cancel at a power rate well above the σ ≥ 0.0035–0.008 that is needed.
2. Selecting by prime-sum amplitude alone induces a phase bias (an Euler-product coupling of Q and L). Large random classes carry the family mean, which comes from the diagonal n = pm. In the scheme the diagonal is a computable secondary term, smaller by about Z^{−(g−θ_0)ℓ} on the contact class, but it would have to be split off before H(σ) can hold. That last point is heuristic.
3. The quarter-arc subsets have the same magnitude profile as their class and show no cancellation. Every input of the scheme is a magnitude statement that stays valid on subsets, so H(σ) cannot follow from those inputs.

These are proxies. Characters modulo q are not sextic Kummer characters, the prime window holds only 4 primes, and no row has a zero off the line.

## 4. Adversarial check: what a proof would need
1. **The inputs cannot give it.** M1, M2, the inverse moments, the pointwise prime bound, the detector lower bounds and the Sobolev step are all invariant under u ↦ e^{iφ_u}·(row data), and monotone under passing to subclasses. Liu's saturation model (2611 l.4196–4262) with phases e^{iφ_u} aligned to the row correction satisfies them all and keeps r = 1. So any proof needs an input that ties phases to arithmetic.
2. **Dual side.** One would need a phase-sensitive estimate for Σ_{u∈𝓑} Q_sel(u)L(w,χ_u)/L(s,ηχ̄_u) over rows where zeros of L(·,ηχ̄_u) are detected. Candidates: a "ratios over thin subfamilies" theorem, equidistribution of arg L'(ρ_u)·ε_u, or repulsion between zeros of the twists L(s,ηχ̄_u). None is known. A further problem: the class is formal and may be empty, so the estimate has to hold for every hypothetical configuration consistent with the counts, and the quarter-arc example shows that magnitude counts alone allow conspiracies.
3. **Direct side.** One would need |Σ_k Ã_kB_k^{(∅)}| ≤ Z^{−λ}‖Ã‖‖B‖ for the empty-subset pairing, uniformly near (ℓ_2, b_2).
   - Under GRH it should hold with λ up to about B_2 − 1/2 (heuristic: the theta average is then about Z^{C(1/2)}).
   - Unconditionally it is equivalent, given the correction bounds, to the new half-plane. So the saving has to come from structure that is independent of zero-free information: Ã_k is essentially a sextic theta coefficient (Kubota–Patterson), and the pairing is a Rankin–Selberg-type integral of theta against an η-twisted sum.
   - A nontrivial bound on its cuspidal part, beyond Cauchy–Schwarz, would suffice; λ = 0.0034 already gives 0.87.
   - Nothing of this kind is known for the sextic cover. The only signed tool here, the large sieve, loses (§2.3).
4. **Proof pieces already in place.** B(σ) and B_2 − λ are exact consequences of the audited scalar system. The principal box re-derivation down to Re s ≥ 0.67 is a finite monomial check. The rest of the chain is unchanged: contours, Mellin argument and descent.
5. **Risks.**
   - Numerics A7 and `infimum_scan.py` use float grids, not Bernstein certificates. A proof at a fixed σ would need the Appendix A certificates redone with E − hσR.
   - The validity range σ ≤ 0.0075 depends on the paper's floor buffer of −1/1000.
   - The finite model is a proxy (§3).

## 5. Files
- `verify.py`: exact sympy checks A1–A8 (dual lever), B1–B3 (direct lever, principal box), C1–C3 (large-sieve deficits), D0–D2 (finite model at q = 65537, 163841). 17 checks, 0 failed.
- `infimum_scan.py`: numerical infimum of the modified system for σ ≤ 0.5, with b re-optimized.
- `finite_model.py`: the §3 table at q = 65537, 163841, 786433, for η = χ_7 and η = 1.

Run: `uv run --no-project --with numpy --with sympy --with mpmath python verify.py`

Sources: [Hexagon 2611.00082v1](https://hexagonmath.org/abs/2611.00082) · [arXiv 2610.26835v2](https://arxiv.org/abs/2610.26835v2) · [Hexagon2 audit](../Hexagon2/REPORT.md) · [R310](../../Tree/R310.lean) · [R109](../../Tree/R109.lean)

OPEN: grhsub's Alladi–Johnson survivor (2(2 − θ)/3 − ε)φ(k) = 0.7511φ(k) (N11) is Montgomery's density exponent, not the method's limit: Chen–Gupta–Li's (qT)^{7(1−σ)/3} (arXiv 2507.08296) in the same mean-value step gives (6/7 − ε)φ(k) = 0.857φ(k) and the tail from X ≥ k^{7/6+ε} instead of k^{1.331} (exponent arithmetic, reviewer CORRECT; routine, off-tree). ROBUST: grhharv4's θ < 3/4 (fails already without a resonator; large sieve X^{1+ε}); morirate's rate 4π (o(1/log T) merges move it); explore36's 2πX² (floats support r ≤ 0.95 for close pairs); R393's Pólya–Wiman (an n₀ needs zeros within 2.3n₀/log γ of the line).

# CeilAudit4 — falsification of our own N1, N4, N6, N10, N11, N12 barriers (ceilaudit4, Orca, 20 Nov 2026)

Selection (RH relevance: criterion, rate, threshold or exponent): (A) N11 grhsub's Alladi–Johnson threshold; (B) N11 grhharv4's θ < 3/4 for 2606.15635; (C) N1 morirate's rate 4π; (D) N1 explore36's height 2πX²; (E) N6 R393's Pólya–Wiman limit (added: N6 had no other candidate). Not chosen: mertensch's exponent floor (Mertens78 REPORT l.14, 53: signed LP certificates already push θ toward 7/8 at rising cost, so the floor 0.8992 is a cost record, not a barrier); explore30 (its block of ≈ 100–750 zeros at δ = 0.37 is already checked far inside explicit zero density); odlyzko31, robinch, weildefl (constants, ranges, cost).

## A. N11 grhsub — Alladi–Johnson signs: "only (2(2 − θ)/3 − ε)φ(k) survives"; tail needs X ≥ k^{1.331}

**Step fixed.** GRHSub NOTES.md l.84–100, Claim (MV): Σ_{χ≠χ₀}|P_y(χ)|² ≪ y^{−η} for y ≥ k^A, P_y(χ) = Σ_{p>y} χ(p)/p, via the explicit formula and a zero-density bound N_k(σ,T) = Σ_{χ mod k} N(σ,T,χ) ≪ (kT)^{c(σ)(1−σ)}; the dyadic zero sum needs (i) c_k(σ) < 2A for the k-exponent and (ii) T-exponent c_T(σ)(1−σ) ≤ 1 (absorbed by 1/|ρ|), on σ ∈ [1/2, θ] (zeros above θ are absent). Markov then gives #{l: C(k,l) < 0} ≥ (1/A − ε)φ(k). So A₀ = sup_{[1/2,θ]} c_eff(σ)/2, fraction 1/A₀.
**Class.** One density input: Montgomery's large-sieve exponent c(σ) = 3/(2 − σ) (Topics Thm 12.1), maximal at σ = θ: A₀ = 3/(2(2 − θ)) = 1.33138, fraction 0.75110.
**Freedoms the lemma admits.** Any single-modulus density bound of the same shape (the proof uses only its k- and T-exponents): Montgomery's large-σ 2/σ (σ ≥ 4/5), Huxley's 12/5, and Chen–Gupta–Li, arXiv 2507.08296 (v2, July 2025), Thm "Zero-density estimate", TeX l.158–185, "in all cases": Σ_{χ mod q} N(σ,T,χ) ≤ (qT)^{o(1)}(q^{7(1−σ)/3}T^{2(1−σ)} + (qT)^{30(1−σ)/13}), 1/2 < σ < 1. Its T-exponents are ≤ 1 for σ ≥ 17/30; below, Montgomery's 3/(2 − σ) ≤ 90/43 < 7/3.
**Test** (a_density.py, logs/a_density.log; 15.6 MB; c_eff = best admissible k-exponent, sup on a 2·10⁵-point grid):

| inputs | sup c_eff (at σ) | A₀ | fraction 1/A₀ |
|---|---|---|---|
| Montgomery only (the barrier) | 2.66275 (θ) | 1.33138 | 0.75110 |
| + large-σ 2/σ | 2.5 (4/5) | 1.25 | 0.8 |
| + Huxley 12/5 | 2.4 (3/4) | 1.2 | 0.83333 |
| + Chen–Gupta–Li 7/3 | 7/3 on [5/7, 6/7] ⊂ [1/2, θ] | 7/6 | 6/7 = 0.857143 |

The check with the CGL terms by hand: y^{−2(1−σ)}k^{7(1−σ)/3}U^{2(1−σ)−1} ≤ [y^{7/(3A)−2}]^{1−σ} (needs A > 7/6); (kU)^{30(1−σ)/13}y^{−2(1−σ)}/U ≤ [y^{30/(13A)−2}]^{1−σ} for σ ≥ 17/30 (needs A > 15/13); with 1 − σ ≥ 1 − θ both give a power saving. The (qT)^{o(1)} with T ≤ T₀ = y^{20} is absorbed.
**Verdict.** OPEN (input upgrade, arithmetic, no new idea): the stated survivor (2(2 − θ)/3)φ(k) is Montgomery's exponent, not the method's limit; the same proof gives C(k,l) < 0 for ≥ (6/7 − ε)φ(k) classes and the tail from X ≥ k^{7/6+ε}, and 6/7 no longer depends on θ (the same for every θ ≥ 5/7). ROBUST part: "φ(k)(1 − o(1)) does not transfer" in the class Parseval + Markov with known density exponents: 1 − o(1) needs A → 1, i.e. c_k ≤ 2 on [1/2, θ], the density hypothesis in the k-aspect on [5/7, 4/5] (open; best 7/3).
Certificate steps: restate (ZD) piecewise (Montgomery on [1/2, 17/30], CGL on [17/30, θ]); re-check NOTES l.94's dyadic bound with separate k- and T-exponents; Markov unchanged; reviewer. Value: routine, off-tree (COROLLARIES), the NOGO line needs its class.

## B. N11 grhharv4 — 2606.15635 (long resonator, L(1/2, χ_{8q}), q prime): "needs θ < 3/4"

**Step fixed.** TeX l.285–319: S₁ = 2Σ_{m,n∈R}ψ(m)ψ(n)Σ_l l^{−1/2}Σ*_q χ_{8q}(lmn)V(l/√q)Φ(q/X); Lemma logd (l.202–208, GRH: Σ_q (log q)χ_{8q}(c)Φ(q/X) = δ_{c=□}Φ̂(1)X + O(X^{1/2+ε}log(c+2))) is applied to each c = lmn and the error summed with absolute values over l ≤ X^{1/2+ε}: X^{3/4+2ε}(Σψ)² ≤ X^{3/4+2ε}|R|Σψ², |R| ≤ N = X^{1/4−5ε}. With Re ρ ≤ θ the lemma gives X^{θ+ε}log(c+2), so the error is X^{θ+1/4+2ε}N Σψ² against the main term X·e^{(...)}Σψ².
**Class.** Twist-uniform prime-sum inputs: the error may depend on the conductor c = lmn only through log c, because the long resonator's moduli are huge: log m ≤ Σ_k a log N(k + 1 + log(log N log₂N))/(k² log₃N), i.e. m up to N^{7.6…19.5} for log N = 10²…10¹⁶ (b_resonator.py), growing like N^{log₂N/log₃N}; any power of c is fatal.
**Freedoms and tests** (b_resonator.py, logs/b_resonator.log, 19.7 MB):
1. Unbalanced AFE (root number +1, lengths L₁L₂ = X): X^θ(√L₁ + √L₂) ≥ 2X^{θ+1/4}; exponents 1.1233, 1.1733, 1.2483 for L₁ = X^{1/2}, X^{0.6}, X^{0.75}: no gain.
2. Shorter resonator: the AFE-summed error X^{θ+1/4} = X^{1.1233} exceeds X already for |R| = 1 (no resonator), so no choice of N helps.
3. Cancellation over l (phase information the lemma discards): the twist-uniform tool is Heath-Brown's quadratic large sieve (coefficients absorb ψ_{mn}): Σ_{l∼L₁}|E_l|² ≪ X^ε(L₁ + X)X log X, so each dyadic l-block contributes X^{1+ε} (exponent 1.0000 at L₁ = X^{0.1}, X^{0.25}, X^{0.5}), against the needed X^{1−η}/N: fails by X^ε·N. The family is shorter (L₁ ≤ X^{1/2}) than the q-sum (X), so the large sieve is Parseval-tight here.
4. Family zero density (Heath-Brown: Σ_{|d|≤Q}N(σ,T,χ_d) ≪ (QT)^ε Q^{3(1−σ)/(2−σ)}T^{(3−2σ)(1−σ)/(2−σ)}): transfers the untwisted first moment (m = n = 1): max_σ [σ − 1/4 + 3(1−σ)/(2(2−σ))] = 0.80051 at σ = 0.7753 < 1, error X^{0.80}; with the twist ψ_{mn} the conductors are ≈ l·mn and the bound is a power of mn: useless. An adversary with one real zero at θ for every twisted character ψ_{lmn} (Siegel-zero exclusion does not reach θ) makes every E_c ≈ −X^θΦ̂(θ) of one sign, so Σ_l l^{−1/2}E_c ≍ X^{θ+1/4}: consistent with every family density bound, since those are polynomial in the conductor.
**Verdict.** ROBUST. Exact class: AFE-expanded first moments twisted by a resonator of moduli X^{Θ(log₂X/log₃X)}, with twist-uniform (log-conductor) inputs: the half-plane's pointwise X^θ log c and the quadratic large sieve; inside it θ < 3/4 is necessary already without a resonator. Outside: a log-conductor-uniform family zero-density estimate for twisted quadratic families (unknown), or a short resonator, which only gives the unconditional shape e^{c√(log X/log₂X)} (FHX26, l.115–119).

## C. N1 morirate — RH + zero statistics cannot force lim (−log λ_min)/μ = 4π

**Step fixed.** MoriRate (ii): λ_min is a monotone functional of the zero multiset; under RH (W(u) = Σ_ρ m_ρ|û(γ_ρ)|², Mori L2.1) deleting one sample at height γ₀ < 4πμ costs e^{−4g₂(γ₀/2πμ)}, g₂(y) = log((2 + √(4 − y²))/y); a pair closer than δ* ≈ e^{−3.24μ+3.6} (l.79) acts as one sample. The step: from RH plus zero statistics to the rate.
**Class.** Zero-side information stable under merging a proportion ε/log T of the zeros at height T into δ*-close pairs: RH, N(T) and S, S₁ bounds (each merged pair moves its zeros by ≤ ½ spacing, so S changes by ≤ ½), every asymptotic statistic with o(1) error (Montgomery's F(α), |α| ≤ 1, changes by ≈ ε(1 − cos 2πα)/log T; simple/distinct proportions by ε/log T), and verification to any fixed height.
**Freedoms examined.**
1. Prime side under RH: primes ≤ X see the zeros only through Σ_ρ e^{iγy}, |y| < log X (explore36); a merge to distance δ* changes it by ≤ δ* log X per pair, invisible unless log X ≳ e^{3.2μ}.
2. Explicit error terms in the statistics: the merge changes F(α) by O(ε/log T), inside the o(1) of Montgomery's theorem, and any ε > 0 already moves the limit (to ≈ 4π(1 + ε)).
3. Moments of ζ′ at zeros: Gonek's Σ|ζ′(ρ)|² ~ T log⁴T/24π² is insensitive (mass from the other zeros); his lower bound Σ|ζ′(ρ)|^{−2} ≫ T only grows under merging. An upper bound Σ_{γ≤T}|ζ′(ρ)|^{−2} ≪ T (Gonek–Hejhal conjecture) would cap the zeros with |ζ′(ρ)| ≤ e^{−cμ} at ≪ Te^{−2cμ} = o(μ) and pass the barrier, but it is a conjecture.
**Test** (c_lostsample.py, logs/c_lostsample.log, 73 MB): mean of g₂ on (0, 2) = 1.570796326795 = π/2 (quad); merging ε/log T of the zeros below 4πμ loses (2εμ)·log(2μ/e)/log 4πμ samples and, density-weighted, adds 0.0306·4πμ nats at μ = 10, ε = 0.1, rising to 0.0883·4πμ at μ = 10¹² (→ ε·4πμ as log μ → ∞): the rate changes by a fixed factor while every statistic in the class moves by o(1).
**Verdict.** ROBUST (model-based, as stated: it obstructs proofs, not implications; the lost-sample law's transfer to λ_min as μ → ∞ remains open, MoriRate (ii′)). Exact class: zero-side inputs stable under merging o(1/log T)-proportions into e^{−3.2μ}-close pairs, i.e. RH + counting + S-bounds + o(1)-asymptotic statistics + finite verification. Outside: a negative-moment upper bound for ζ′(ρ) (conjectural) or any statement giving all but o(μ) of the zeros below 4πμ simple and e^{−3.2μ}-separated.

## D. N1 explore36 — primes ≤ X cannot exclude near-line pairs above ≈ 2πX²

**Step fixed.** ExploreMinimax §1–2: Le Cam two-point lemma (D2) with the Weil-window distance η_c(Z₀, Z₁) = ‖T_{ΔG}‖ on L²(−a, a), 2a = log c ≤ log X; Z₁ merges two consecutive ordinates into ½ ± δ + iγ̄ (δ = 1/(2L*), L* = ½log(t/2π)) and moves the M nearest ordinates by ε_k (Levenberg–Marquardt on the Hilbert–Schmidt norm, |ε_k| ≤ CAP_SP spacings); ρ = η_c(Z₀, Z₁^ε)/η_c(Z₀, Z₁⁰). Threshold variable r = log X/L* (r = 1 ⟺ t = 2πX²).
**Class.** Arguments that use primes ≤ X only through windows c ≤ X and keep a margin above ρ·η_bare; local, order-preserving moves; the endpoint "≈ 2πX²" is the Landau-density heuristic (their D4, "a conjecture, not a theorem"); exact invisibility needs infinitely many moves (D3, conjectural).
**Freedoms.** (1) Margin below ρ·η_bare, i.e. exact use of the prime data: outside the class by construction; N1's actual window certificates have margins ≈ e^{15}μ⁵e^{−4πμ} but by explore19 positivity is blind to single pairs above ≈ 13.9c anyway, so for certificates the binding limit is 13.9c, not 2πc². (2) The adversary's own freedoms near r → 1: more neighbours (M), larger moves (cap 0.45 → 0.95 spacing), choice of pair. Tested with explore36's own compensate/opnorm (lecam1.py) on ζ's zeros at t ≈ 10⁴ (L* = 3.686), two pairs (gaps 0.69 and 0.87 mean spacings):

| pair gap | cap | r = 0.8 | r = 0.9 | r = 0.95 | r = 1.0 |
|---|---|---|---|---|---|
| 0.69 sp | 0.45 | – | 1.7e-2 / 4.6e-3 / 2.3e-3 | 1.8e-2 / 4.8e-3 / 2.4e-3 | 2.1e-2 / 1.0e-2 / 8.1e-3 |
| 0.69 sp | 0.95 | 5.8e-3 / 2.2e-5 / – | 1.7e-2 / 4.6e-3 / – | – | 2.1e-2 / 1.0e-2 / – |
| 0.87 sp | 0.45 | – | 0.26 / 0.24 / 0.23 | 0.43 / 0.42 / 0.42 | 0.55 / 0.52 / 0.51 |
| 0.87 sp | 0.95 | 6.5e-3 / 1.8e-3 / – | 0.19 / 0.13 / – | – | 0.54 / 0.52 / – |

(ρ_op at M = 100 / 400 / 800; logs/d_threshold_cap045.log, d_threshold_cap095.log; 16 + 3.6 min, peak 203 MB by time -l.)
**Reading.** A close pair (gap ≲ 0.69 spacing: 22.7 % of the 11,661 normalized gaps on 5·10³–1.5·10⁴) is hidden to ρ ≈ 2·10⁻³ up to r = 0.95 and to ≈ 10⁻² at r = 1 (still falling in M); a typical pair (0.87 spacing) is hidden only for r ≤ 0.8 (ρ = 1.8·10⁻³) and stays visible at r = 0.9 (ρ ≈ 0.13 even with moves up to 0.73 spacing). Since the adversary chooses the pair, the barrier stands, but its numerical support is r ≤ 0.95 (t ≥ 2πX^{2.11}) at t = 10⁴, not r = 1.
**Verdict.** ROBUST (information barrier, floats). Exact class: zero-side window functionals for c ≤ X with margin ≥ ρ·η_bare, ρ ≈ 2·10⁻³ at log X ≤ 0.95·½log(t/2π) (close pairs, M ≤ 800, moves < 0.5 spacing); the endpoint 2πX² is heuristic. No freedom of the lemma passes it; it does not bind N1's certificates, whose positivity horizon is ≈ 13.9c.

## E. N6 R393 — "Pólya–Wiman needs finitely many nonreal zeros, so no unconditional n₀ for all d"

**Step fixed.** R393 (Wild): J^{d,n} hyperbolic for all d ⟺ F^{(n)} ∈ LP, F(x) = Σγ(k)x^k/k! (genus 0, Ξ(z) ∝ F(−z²)); Pólya–Wiman (Craven–Csordas–Smith 1987, Kim 1990) gives F^{(n)} ∈ LP for n ≥ n₀ only if F has finitely many nonreal zeros. Holland's exponent 5/3 against the two-thirds wedge's 3/2 (same round) was not re-examined (superseded by savvak, R393).
**Class.** Zero-location inputs: the Lean θ₁₄₁ strip (|δ| ≤ 0.3733, δ = Re ρ − ½, i.e. the parabola |Im x| ≤ 0.75|x|^{1/2} for F), RH to height T, genus 0.
**Freedom.** The derivative operator itself: x-derivatives might heal infinitely many off-line zeros (as derivatives of Ξ push zeros toward the line, Conrey 1983), which would give an n₀ without finiteness. Mechanism: on an oscillating stretch, d/dx acts on a local factor r as (1 + (Δx/iπ)d/dx), i.e. n derivatives translate r by ∓inΔx/π, so a pair with offset δ_z (in local z-spacings) is healed after n ≈ πδ_z derivatives.
**Test** (e_healing.py, logs/e_healing.log; mpmath 40 digits, Cauchy-integral derivatives, Newton continuation; 15 MB): G(x) = sinh(π√x)/(π√x) (zeros z = k), the zeros z = K, K + 1 replaced by K + ½ ± iδ_z. |Im z| of the tracked zero of F^{(n)}: δ_z = 0.3: healed at n = 1; δ_z = 1: n = 2 (K = 40) or 1 (K = 200); δ_z = 3: 3, 2.61, 2.23, 1.85, 1.48, 1.10, 0.75, 0.21–0.35, then real at n = 8, for K = 40 and 200 alike (≈ 0.38 per derivative): n_heal ≈ 2.7δ_z, independent of K.
**Consequence.** For ζ, δ_z = δ log(γ/2π)/2π, so F^{(n)} keeps a nonreal zero for each off-line zero with δ ≳ 2.3n/log(γ/2π). An unconditional n₀ would force every zero into σ ≤ ½ + 2.3n₀/log(γ/2π), a zero-free region of width O(1/log γ) beside the line: RH-strength, and every strip or verification input admits models (off-line pairs with fixed δ at all heights) in which no n₀ exists.
**Verdict.** ROBUST. Exact class: zero-location inputs (θ₁₄₁ strip/parabola, finite verification, genus 0); a large-n result for all d needs off-line zeros within O(n₀/log γ) of the line, i.e. near-RH input. (Heuristic plus model test, not a proof.)

## Verification and novelty
- Strongest claim (A): one reviewer subagent, formulas pasted inline (MV dyadic bound, CGL's two terms, T-exponents, o(1) loss, arithmetic 7/6, 6/7, 17/30, 5/7): CORRECT ("within the stated density-based MV reduction"; Montgomery needs A > 45/43 on [1/2, 17/30], CGL A > max(7/6, 15/13); the 6/7 is a limit of these inputs and this reduction, not absolute).
- CGL statement read in its TeX (l.158–185, "In all cases"); arXiv abs: Chen, Gupta, Li, 11 Jul 2025, v2.
- Novelty (A): python3 .scratch/lst.py 'Alladi|smallest prime factor|least prime factor|zero[- ]density|Dirichlet L-functions.*zeros' math.NT --zenodo '"Alladi" "smallest prime factor"': arXiv 2026-10/11 titles only 2610.15989, 2610.37371 (zero density, other families); Zenodo nothing relevant (98 unrelated). Our files: 2507.08296 is already used for Linnik's 7/3 (COROLLARIES l.3) and densityq, never for Alladi–Johnson (Alladi appears only in GRHSub, COROLLARIES l.70, TREE N11, NOGO). Verdict: routine input upgrade (COROLLARIES), not a discovery.
- Not certified here (brief): the 6/7 needs the NOTES l.84–100 write-up with the piecewise (ZD) and the orchestrator's COROLLARIES edit.

## Suggested NOGO lines
- grhsub (replace "only (2(2 − θ)/3 − ε)φ(k) survives"): with the mean-value route (Parseval + Markov) the survivor is (1/A₀ − ε)φ(k), A₀ = sup_{[1/2,θ]} c_k(σ)/2 over k-aspect density exponents with T-exponent ≤ 1: 0.7511 (Montgomery), 6/7 (Chen–Gupta–Li 7/3); φ(k)(1 − o(1)) needs the k-aspect density hypothesis on [5/7, 4/5] [CeilAudit4](Riemann/Deep/CeilAudit4/REPORT.md)
- grhharv4: 2606.15635 under Re ρ ≤ θ needs θ < 3/4 in the class of AFE-expanded, resonator-twisted first moments with log-conductor-uniform inputs (pointwise X^θ log c; quadratic large sieve gives X^{1+ε}); it fails already with |R| = 1 (X^{θ+1/4}); Heath-Brown's family density transfers only the untwisted moment (X^{0.80}) [CeilAudit4]
- morirate: class = zero-side inputs stable under merging o(1/log T)-proportions into e^{−3.2μ}-close pairs (RH, counting, S-bounds, o(1) statistics, finite verification); outside: a negative-moment bound Σ|ζ′(ρ)|^{−2} ≪ T (conjectural) [CeilAudit4]
- explore36: class = window functionals c ≤ X with margin ≥ ρ·η_bare; floats: ρ ≈ 2·10⁻³ for log X ≤ 0.95·½log(t/2π) (pairs with gap ≲ 0.69 spacing, M ≤ 800, moves < 0.5 spacing, t ≈ 10⁴); typical pairs only for ≤ 0.8; 2πX² itself is the Landau-density heuristic; N1 certificates are blind above ≈ 13.9c anyway [CeilAudit4]
- R393 Pólya–Wiman: class = zero-location inputs (θ₁₄₁ strip, verification, genus 0); x-derivatives heal an off-line pair after ≈ 2.7δ_z derivatives (δ_z in local spacings; model, K-independent), so an unconditional n₀ for all d needs every zero in σ ≤ ½ + 2.3n₀/log(γ/2π) [CeilAudit4]

RH relevance: N11 (A: GRH-to-half-plane transfer threshold A₀; B: the θ < 3/4 transfer limit), N1 (C: the λ_min rate; D: the prime-length/height threshold), N6 (E: the large-n half of the Jensen criterion).

## Files
- Riemann/Deep/CeilAudit4/REPORT.md (this file)
- a_density.py (A: A₀ and fractions per density input), b_resonator.py (B: exponents, resonator moduli), c_lostsample.py (C: g₂ mean, merge bookkeeping), d_threshold.py (D: explore36's compensate/opnorm near r → 1; imports Riemann/Deep/ExploreMinimax/lecam1.py), e_healing.py (E: derivative healing model)
- logs/: a_density.log, b_resonator.log, c_lostsample.log, d_threshold_cap045.log, d_threshold_cap095.log, e_healing.log (+ .err with memcap peaks and /usr/bin/time -l footprints: 15.6, 19.7, 73, 203/93, 35 MB)
- Compute: ≈ 22 min total, one process at a time, single-threaded, every process < 300 MB.

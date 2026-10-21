# Density3: can the ring class zero-density exponent be improved? Proofs and barriers

Date 2026-10-21. Brief: `.scratch/orca/density3-brief.md`. Scratch: `.scratch/density3/` (`qa1.py` → `qa1.out`
numerical checks; `axq.py`, `axq2.py` → `axq.out`, `axq2.out` arXiv searches; `src/yang/` TeX of arXiv 2607.04476).
References: [C] = [Combos NOTE](../Combos/NOTE.md) (Lemmas 2–5, Thm 1), [D2] = [Density2 NOTE](../Density2/NOTE.md)
(Lemma 4, Prop 6), [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md) (Prop 1, Thm 2′, §4). Notation as there: D < 0, G = Cl(O),
h = |G| ≍ |D|^{1/2±ε}, N_D(σ;t) = Σ_{χ∈Ĝ} #{ρ of L(s,χ*) : β ≥ σ, |γ − t| ≤ 1}, T = |t| + 3,
e(σ) := 2(1−σ)/(3−2σ) (the current exponent, [R2] Thm 1′(a)). All exponents are of |D| unless "h-units" is said.

## Verdicts

| # | Item | Result | Status |
|---|---|---|---|
| 1 | Is e(σ) improvable with the present tools? | No: over {weighted large sieve on k-th powers of dyadic blocks, first moment, twisted second moment [C] Lemma 4, any mollifier length} the optimum is exactly e(σ). The pure Carlson route (no class II) also ends at e(σ): its worst block is N = h^{1/(3−2σ)} (Prop 1) | barrier, proved |
| 2 | Fourth moment input | If Σ_{χ≠1}\|L(½+iu,χ)\|⁴ ≪ \|D\|^{a+ε}(\|u\|+3)^B, then N_D(σ;t) ≪ \|D\|^{min(e(σ), (1+a)(1−σ)/(2−σ))+ε}T^{O(1)}; this beats e(σ) iff a < 1/(3−2σ) (Thm 2). Only a = 1 is known | conditional, proved |
| 3 | Subconvexity input | If L(½+iu,χ) ≪ (\|D\|(\|u\|+3)²)^{¼−δ+ε} for all χ ≠ 1, then N_D(σ;t) ≪ (\|D\|^{1−2δ}T)^{2(1−σ)+ε}, better than e(σ) iff σ > (1−3δ)/(1−2δ) (Thm 3). Unconditional instances: δ > 0 (Michel–Venkatesh 2010), δ = 1/1889 (Blomer–Harcos–Michel 2007, as quoted in arXiv 2607.04476), δ = 1/120 claimed (Yang, arXiv 2607.04476, July 2026): σ > 0.99947, resp. σ > 0.99153 | new statement, routine, microscopic |
| 4 | Barrier for the fourth-moment route | Improving e at σ₀ this way needs a < 1/(3−2σ₀), hence max_χ\|L(½,χ)\| ≪ \|D\|^{1/(4(3−2σ₀))}: at σ₀ = 5/6, Burgess strength \|D\|^{3/16} (best known ¼ − 1/120, claimed). Halász–Montgomery at the critical block needs \|D\|^{(4σ₀−3)/(4(3−2σ₀))} = \|D\|^{1/16} at σ₀ = 5/6 (beyond Weyl) (Prop 4) | barrier, proved |
| 5 | Almost-all threshold A₀ = 3/2 of [R2] | unchanged: the gains of items 2–3 live at σ > 0.99 and the binding range is σ ∈ [0.81, 0.84]. A per-D fourth moment lowers A₀ iff a < 3/4 (A₀ = 1.4919 at a = 0.74, 1.2944 at a = ½) | no change |
| 6 | Halász–Montgomery with the lattice structure | below \|D\|^{1/2} the off-diagonal is a sum over primitive ideals, ≈ N^{1/2}\|L(½,ψ)\|; no Pólya–Vinogradov-type saving exists; the group structure gives no subdivision (all ψ ≠ 1 have conductor ≍ \|D\| for fundamental D); an exceptional-set version with the known fourth moment is exactly neutral (Lemma 5, §4) | NOGO |
| 7 | Lead | λ_π(p²) = 1 at every inert p for every θ_χ, so one amplifier is maximal on the whole ring class family inside F₀(\|D\|, χ_D); an amplified fourth moment of length L for that GL₂ family gives Σ_χ\|L(½,χ)\|⁴ ≪ \|D\|^{1+ε}/L (§6). Not verified | lead |

## 0. Inputs used

- (LS) [D2] Lemma 4 (M = 0): Σ_{χ∈Ĝ}|Σ_𝔞 c_𝔞χ(𝔞)|² ≤ ε^{−1}L_D(2N + (w/2π)|D|^{1/2})‖c‖²_ε for c supported on (𝔞,f) = 1,
  N𝔞 ≤ N, ‖c‖²_ε = Σ m(𝔞)^{1+2ε}|c_𝔞|². Below N = √|D|/2 it is an identity on primitive ideals (Lemma 5(a)).
- (M2) [D2] (2.1) with X = 1, and with X = V: for 1 ≤ V ≤ |D| and u ∈ R,
  Σ_{χ≠1}|L(½+iu,χ)|² ≤ C(ε)|D|^{½+ε}(|u|+3)^{1+ε};  Σ_{χ∈Ĝ}|M_V(½+iu,χ)|² ≤ C(ε)(V + |D|^{½})V^{ε}
  (the second by (LS) and [C] Lemma 3; M_V(s,χ) = Σ_{N𝔡≤V,(𝔡,f)=1} μ_K(𝔡)χ(𝔡)N𝔡^{−s}).
- (Det) [C] Lemma 5 (zero detector) with mollifier length V ≤ Y, and the class I count of [R2] Prop 1:
  if V = |D|^{1/2}, Y ≥ V, Z = 2Y log²Y, the class I characters number ≤ C(ε)(|D|T)^{ε}Z^{2−2σ+ε}.
- Class II condition of (Det): ∫|L M_V(½+i(γ_χ+u),χ)|e^{−|u|}du ≥ c₁Y^{σ−½}, with |γ_χ − t| ≤ 1.

## 1. The present toolkit stops exactly at e(σ)

Write N = h^α for a dyadic block (h-units). The weighted large sieve applied to the k-th power of a block with coefficients
|b_𝔞| ≤ d(𝔞)N^{−σ} gives at most h^{f_k(α)+o(1)} large values, f_k(α) := max(kα, 1) + kα(1−2σ) (constant N^k + h, then
[C] Lemma 3). Put g(α) := min_{k≥1} f_k(α) and α* := 1/(3−2σ) ∈ (½, 1).

**Proposition 1.** Let ½ < σ < 1.
- (a) g(α) ≤ g(α*) = 4(1−σ)/(3−2σ) = 2e(σ) for all α ∈ (0, 1], with equality only at α*.
- (b) [C]'s Carlson argument with k chosen per block as the minimiser of f_k, mollifier length |D|^{ε}, class II by [C]
  Lemma 4, and Y := |D|^{1/(2(3−2σ))} (instead of [C]'s |D|^{σ−½}) gives N_D(σ;t) ≪ |D|^{e(σ)+ε}T^{O(1)}: a second proof of
  [R2] Thm 1′(a) in the |D|-aspect.
- (c) No combination of: the detector (Det) with any mollifier length X = h^x, class I blocks bounded by g, class II bounded by
  the first moment (≤ h·U^{−1}) or by [C] Lemma 4 (≤ hX·U^{−2}), U = Y^{σ−½}, and any Y, yields an exponent below e(σ).

*Proof.* (a) For kα ≤ 1, f_k = 1 + kα(1−2σ) decreases in k; for kα ≥ 1, f_k = kα(2−2σ) increases in k. So
g(α) = min(f_{⌊1/α⌋}, f_{⌈1/α⌉}), and the local maxima of g sit at the crossings 1 + kα(1−2σ) = (k+1)α(2−2σ), i.e.
α_k = 1/(k+2−2σ), where g(α_k) = (2−2σ)(k+1)/(k+2−2σ). This decreases in k because 2−2σ < 1; k = 1 gives α₁ = α* and 2e(σ).
(b) Class II: [D2] (2.1) with X = |D|^{ε} gives ≪ |D|^{½+ε}Y^{1−2σ}T^{1+ε} = |D|^{e(σ)+ε}T^{1+ε}, since
½ + (1−2σ)/(2(3−2σ)) = (4−4σ)/(2(3−2σ)). Class I: every block has N ≤ Z ≪ Y log²Y = h^{α*+o(1)}, so by (a) and the
Sobolev step of [C] Thm 1 each of the O(log|D|) blocks contributes ≪ |D|^{e(σ)+ε}; k ≤ 2/ε, so the constants are C(ε).
(c) In h-units let Y = h^y, X = h^x, E = the resulting exponent. If y < α*: the first-moment bound is 1 + y(½−σ) > 3/2 − σ
≥ 2e(σ) (the difference (3/2−σ) − 2e(σ), times (3−2σ), is 2(σ−½)² ≥ 0), and the Lemma 4 bound is > 1 + α*(1−2σ) = 2e(σ).
So y ≥ α*. If x ≤ α*, the block N = h^{α*} is a class I block and costs 2e(σ) by (a). If x > α*: with the first moment
and y < 1, the class II term alone exceeds 3/2 − σ ≥ 2e(σ); with the first moment and y ≥ 1, E ≥ max(y(2−2σ), 1 + y(½−σ)) ≥
its value at the crossing y = 2/(3−2σ), which is 2e(σ); with Lemma 4 and y ≥ 1,
E ≥ max(y(2−2σ), 1 + x + y(1−2σ)) ≥ (1+x)(2−2σ) > (1+α*)(2−2σ) > 2e(σ); with Lemma 4 and α* ≤ y < 1, the class II term is
≥ 1 + x + (1−2σ) > 2 − 2σ + α*, and (2−2σ+α*)(3−2σ) − (4−4σ) = 4σ² − 6σ + 3 > 0. □

qa1.py (1)–(2) confirms numerically: the worst block sits at α = 1/(3−2σ) with value e(σ) for σ ∈ {0.55, …, 0.99}, and a
grid search over (x, y, class II option) returns e(σ) to 5 digits for σ ∈ {0.6, 0.75, 5/6, 0.9, 0.97}.

So the two routes ([C] re-optimised, and [R2]) end at the same exponent for the same reason: at length N = |D|^{1/(2(3−2σ))}
the large sieve bound hN^{1−2σ} and the squared bound N^{4−4σ} coincide, and the long-mollifier route meets class II at
Y = |D|^{1/(3−2σ)} = N². Any improvement needs an input that is not a second moment.

## 2. Fourth moments

**Hypothesis (M_a).** For every u ∈ R and ε > 0: Σ_{χ∈Ĝ, χ≠1}|L(½+iu,χ)|⁴ ≤ C(ε)|D|^{a+ε}(|u|+3)^{B}, with B absolute.

(M₁) holds: approximate functional equation of L(s,χ*)² (conductor ≤ |D|²(|u|+3)⁴, length |D|(|u|+3)²), the factors E_χ as in
[C] Lemma 4, (LS) and [C] Lemma 3. The full GL₂ family gives the same a = 1 by positivity (Yang, arXiv 2607.04476, Thm E with
F = Q, q = |D|, ω = χ_D; it covers the cuspidal θ_χ, χ² ≠ 1). Lindelöf on average would be a = ½.

**Theorem 2.** Assume (M_a), ½ ≤ a ≤ 1. For 0.6 ≤ σ ≤ 1, t ∈ R, ε > 0:

  N_D(σ;t) ≤ C(ε)|D|^{E_a(σ)+ε}T^{B+2},  E_a(σ) := min(e(σ), (1+a)(1−σ)/(2−σ)).

The second entry is the smaller one iff a < 1/(3−2σ), i.e. iff σ > (3 − 1/a)/2.

*Proof.* The first entry is [R2] Thm 1′(a). For the second, run (Det) with V = |D|^{½} and Y ≥ V free. Class I: (Det) gives
≤ C(ε)(|D|T)^{ε}Y^{2−2σ+ε}. Class II: with dμ(v) := e^{−|v−t|}dv, |γ_χ − t| ≤ 1, and Hölder in (χ, v) with exponents 4, 4/3,

  c₁Y^{σ−½}|𝒳_II| ≤ e∫Σ_{χ∈𝒳_II}|L M_V(½+iv,χ)|dμ ≤ e(∫Σ_{χ≠1}|L(½+iv,χ)|⁴dμ)^{1/4}(∫Σ_{χ∈𝒳_II}|M_V(½+iv,χ)|^{4/3}dμ)^{3/4}.

By Hölder in χ, Σ_{𝒳_II}|M_V|^{4/3} ≤ |𝒳_II|^{1/3}(Σ_χ|M_V|²)^{2/3}; by Hölder in v the last factor is
≤ |𝒳_II|^{1/4}2^{1/4}(∫Σ_χ|M_V|²dμ)^{1/2} ≤ C(ε)|𝒳_II|^{1/4}|D|^{¼+ε} by (M2). With (M_a) and ∫(|v|+3)^{B}dμ ≤ C T^{B}:

  |𝒳_II|^{3/4} ≤ C(ε)Y^{½−σ}|D|^{(1+a)/4+ε}T^{B/4},  so  |𝒳_II| ≤ C(ε)Y^{(2/3)(1−2σ)}|D|^{(1+a)/3+ε}T^{B/3}.

Choose Y := |D|^{(1+a)/(4−2σ)}T^{B/(4−2σ)} (≥ |D|^{½} because a ≥ ½ ≥ 1 − σ). Both classes are then
≤ C(ε)|D|^{(1+a)(1−σ)/(2−σ)+2ε}T^{B+1}; each χ has ≤ c₀log(|D|T²) zeros in the window. Finally
(1+a)/(2−σ) < 2/(3−2σ) ⟺ a < 1/(3−2σ). □

With a = 1 the second entry is 2(1−σ)/(2−σ) > e(σ): the known fourth moment never helps.

## 3. Subconvexity: the only unconditional gain, near σ = 1

**Hypothesis (H_δ).** For every χ ∈ Ĝ, χ ≠ 1, u ∈ R, ε > 0: |L(½+iu,χ*)| ≤ C(ε)(|D|(|u|+3)²)^{¼−δ+ε}.

Known cases. (i) Some absolute δ > 0: for χ² ≠ 1, L(s,χ*) = L(s, AI(χ*)) with AI(χ*) cuspidal on GL₂/Q of conductor
A(χ*) ≤ |D|; apply the uniform GL₂ subconvexity of Michel–Venkatesh (Publ. Math. IHÉS 111 (2010)) to AI(χ*) ⊗ |·|^{iu}. For χ² = 1, χ ≠ 1, L(s,χ*)
is a product of two Dirichlet L-functions whose conductors multiply to A(χ*), and hybrid Burgess gives δ = 1/16.
(ii) δ = 1/1889: Blomer–Harcos–Michel, Ann. Sci. ÉNS 40 (2007), Cor 1, as quoted in arXiv 2607.04476 (its uniformity in u and
in non-maximal orders was not checked here). (iii) δ = 1/120 for all orders of quadratic extensions, uniformly in u: Yang,
arXiv 2607.04476 (July 2026), Cor 1.3 — unrefereed.

**Theorem 3.** Assume (H_δ), 0 < δ < ¼. For 0.6 ≤ σ ≤ 1, t ∈ R, ε > 0:

  N_D(σ;t) ≤ C(ε)(|D|^{1−2δ}T)^{2(1−σ)}(|D|T)^{ε}.

It beats e(σ) iff σ > σ_δ := (1−3δ)/(1−2δ). So N_D(σ;t) ≪ |D|^{min(e(σ), 2(1−2δ)(1−σ))+ε}T^{O(1)}, and Theorem 2 with
a = 1 − 2δ (implied by (H_δ) and (M2)) adds nothing.

*Proof.* (Det) with V = |D|^{½}, Y ≥ V. Class I as before. Class II: |L(½+iv,χ)| ≤ d(f)²|L(½+iv,χ*)| ≤ C(ε)|D|^{¼−δ+2ε}
(|v|+3)^{½}, so for χ ∈ 𝒳_II, by Cauchy–Schwarz in u,
c₁Y^{σ−½} ≤ C(ε)|D|^{¼−δ+2ε}T^{½}(∫|M_V(½+i(γ_χ+u),χ)|²e^{−|u|}du)^{½}. Square, sum over 𝒳_II, use e^{−|u|} ≤ e^{1−|v−t|}
after v = γ_χ + u, and (M2):

  c₁²Y^{2σ−1}|𝒳_II| ≤ C(ε)|D|^{½−2δ+4ε}T·e∫e^{−|v−t|}Σ_χ|M_V(½+iv,χ)|²dv ≤ C(ε)|D|^{1−2δ+5ε}T.

Take Y := |D|^{1−2δ}T. Then |𝒳_I| + |𝒳_II| ≤ C(ε)(|D|^{1−2δ}T)^{2−2σ+6ε}, and the window count finishes.
Comparison: 2(1−2δ)(1−σ) < 2(1−σ)/(3−2σ) ⟺ (1−2δ)(3−2σ) < 1 ⟺ σ > σ_δ. Endpoints: Hölder with exponent p ∈ [2, ∞] on L
interpolates the class II bounds hU^{−1} and |D|^{1−2δ}U^{−2} geometrically; after balancing with Y^{2−2σ}, the exponent is
monotone in the interpolation parameter, with the constant sign of (1−σ) − δ(3−2σ). The same computation shows that Theorem 2
with a = 1 − 2δ crosses e(σ) at the same σ_δ and is weaker beyond it. □

Numbers (qa1.py (3)): σ_δ = 0.999957 (δ = 1/23041, Duke–Friedlander–Iwaniec 1995), 0.999470 (1/1889), 0.991525 (1/120),
0.928571 (1/16, Burgess strength, unknown for this family). At σ = 1 − δ/4 the gain is about 1% of the exponent. This is the
classical "sup bound × mean value" route: new as a statement for this family, routine, and of no use below σ = 0.99.

## 4. Barriers

**Lemma 5 (the Halász–Montgomery off-diagonal).** Let N ≥ 1, w ∈ C_c^∞((½,2)), ψ ∈ Ĝ, ψ ≠ 1, and
O_N(ψ) := Σ_{𝔟 primitive, (𝔟,f)=1} ψ(𝔟)w(N𝔟/N).
- (a) If N < √|D|/2, distinct primitive ideals prime to f of norm ≤ N lie in distinct classes. So for c supported on them,
  Σ_χ|Σ_𝔟 c_𝔟χ(𝔟)|² = h Σ_𝔟|c_𝔟|² (Parseval), and the weighted large sieve is an identity below the threshold.
- (b) L(s,ψ) = ζ^{(f)}(2s)·Σ_{𝔟 prim} ψ(𝔟)N𝔟^{−s}, ζ^{(f)}(2s) := Σ_{(m,f)=1}m^{−2s}; hence
  |O_N(ψ)| ≤ C_w(ε)|D|^{ε}N^{½}·max_u |L(½+iu,ψ)|(|u|+3)^{−2}, and trivially |O_N(ψ)| ≤ #{𝔟 prim : N𝔟 < 2N}.
- (c) With (H_δ), Halász–Montgomery over R ⊂ Ĝ (normalised block, |b_𝔞| ≤ d(𝔞)N^{−σ}, large values ≥ 1) gives
  |R| ≤ 2N^{2−2σ}|D|^{ε} as soon as N ≥ |D|^{(¼−δ)/(2σ−3/2)+ε}. This is below the large sieve threshold |D|^{½} only for
  σ > 1 − δ.

*Proof.* (a) If [𝔟] = [𝔟′], then 𝔟𝔟̄′ = (α) with α ∈ O prime to f (as in [D2] Lemma 2(i)); α = (x + y√D)/2, and
N(α) = N𝔟N𝔟′ < |D|/4 forces y = 0. So 𝔟𝔟̄′ = (n), n ∈ Z, 𝔟̄′ = (n/N𝔟)𝔟̄, and primitivity of 𝔟̄, 𝔟̄′ gives n = N𝔟, 𝔟′ = 𝔟.
(b) Every 𝔞 prime to f is uniquely m𝔟 with (m,f) = 1, 𝔟 primitive, and ψ(m𝔟) = ψ(𝔟). Mellin inversion, then the line
Re s = ½: L(s,ψ) is entire, ζ(1+2iu) ≠ 0 with 1/|ζ(1+2iu)| ≪ log(|u|+3), Π_{p|f}(1 − p^{−1−2iu})^{−1} ≤ d(f), and ŵ decays.
(c) Halász–Montgomery: Σ_{χ∈R}|P(χ)|² ≤ ‖b‖²(cN + |R| max_{ψ≠1}|O_N(ψ)|) (after the decomposition 𝔞 = m𝔟 of [C] Lemma 2,
with the m-weight as there); ‖b‖² ≤ N^{1−2σ+ε}. With (b) and (H_δ) the second term is ≤ |R|/2 once
N^{3/2−2σ}|D|^{¼−δ+2ε} ≤ 1/2. Finally (¼−δ)/(2σ−3/2) < ½ ⟺ σ > 1 − δ. □

For Dirichlet characters mod q, Pólya–Vinogradov (|Σ_{n≤N}χ(n)| ≪ √q log q) makes Halász–Montgomery effective from
N ≥ q^{1/(4σ−2)}, below the threshold q for every σ > ¾; that is the source of Montgomery's density 2(1−σ)/σ near σ = 1. Here
no analogue exists: (b) is the best uniform bound, the off-diagonal ψ = χχ̄′ is again a member of the family, and for
fundamental D every ψ ≠ 1 has conductor A(ψ*) = |D|, so Huxley's subdivision (into pieces whose off-diagonals have small
conductor) has nothing to subdivide. For D = f²d_K with |d_K| fixed, ψ on principal ideals is a character of (O_K/f)^*/((Z/f)^*O_K^*), so O_N(ψ) is a
combination of h(d_K) two-dimensional character sums over discs; completing it costs f ≍ h, so the two-dimensional Pólya–Vinogradov bound is again
of the size of the large sieve threshold, and Chang's Burgess bound in F_{p²} (GAFA 2009) saves only an unspecified p^{−δ}.

**Exceptional sets do not help either.** Let Bad_U := {ψ ≠ 1 : max_{|u|≤T₀}|L(½+iu,ψ)| ≥ U}. By (M_a) and a Sobolev step,
|Bad_U| ≤ C(ε)|D|^{a+ε}T₀^{B+1}U^{−4}. In Halász–Montgomery, pairs (χ,χ′) with χχ̄′ ∈ Bad_U cost at most N·|Bad_U| per χ,
the rest N^{½}U each; with U = N^{2σ−3/2}/2 this yields |R| ≤ N^{2−2σ+ε}(1 + |Bad_U|) ≪ N^{2−2σ+ε} + |D|^{a+ε}N^{8−10σ}.
This beats the large sieve bound |D|^{½}N^{1−2σ} iff N^{8σ−7} > |D|^{a−½}; for a = 1 that needs N > |D|^{½}.
So with the known fourth moment the exceptional-set device is exactly neutral below the threshold.

**Proposition 4 (what an improvement costs).** Let ½ < σ₀ < 1.
- (a) If (M_a) holds with a < 1/(3−2σ₀) (the condition of Theorem 2 at σ₀), then |L(½,χ)| ≤ C(ε)|D|^{1/(4(3−2σ₀))+ε} for
  every χ ≠ 1, i.e. (H_δ) at u = 0 with δ = (1−σ₀)/(2(3−2σ₀)). For σ₀ = 5/6 this is |D|^{3/16}: Burgess strength. Known:
  ¼ − 1/1889 (refereed), ¼ − 1/120 (claimed, arXiv 2607.04476).
- (b) Fed with the bound of Lemma 5(b), Halász–Montgomery beats the critical block N = |D|^{1/(2(3−2σ₀))} of Prop 1 iff
  max|L(½+iu,ψ)| ≤ |D|^{(4σ₀−3)/(4(3−2σ₀))−ε}; at σ₀ = 5/6 this is |D|^{1/16}, far beyond Weyl strength |D|^{1/6}.

*Proof.* (a) |L(½,χ)|⁴ is one term of the sum in (M_a). (b) The condition of Lemma 5(c) at N = |D|^{1/(2(3−2σ₀))} reads
N^{2σ₀−3/2} ≥ 2 max|L|·|D|^{ε}, and N^{2σ₀−3/2} = |D|^{(4σ₀−3)/(4(3−2σ₀))}. □

The only inputs that would move e(σ) in the middle range are therefore of the strength of Burgess-type subconvexity for all
class group L-functions or of a Lindelöf-type fourth moment over a thin family (size |D|^{½}, conductor |D|), the latter being
equivalent to an averaged form of the Michel–Venkatesh mixing conjecture for Heegner points (Σ_χ|L(½,χ)|⁴ is, by Hecke's formula
L(½,χ) ∝ |D|^{−1/4}Σ_C χ(C)E*(z_C,½), equal to h·Σ_C|Σ_{C₁}F(C₁)F̄(C₁C)|² with F(C) ∝ |D|^{−1/4}E*(z_C,½)); that conjecture is
known under GRH (Blomer–Brumley–Khayutin, arXiv 2212.06280) and, without rate, under extra hypotheses (Khayutin, Annals 2019).

## 5. Consequences

- **A₀ = 3/2 of [R2] Thm 2′ is unchanged.** Its binding range is σ ∈ [σ₁(A), σ₂(A)] ⊂ [0.80, 0.84] (A = 1.45: [0.810, 0.835]),
  where Theorems 2–3 with the known inputs give back e(σ). A per-D fourth moment lowers A₀ iff a < ¾: at A = 3/2 the per-D
  range must pass σ = 5/6, where e(5/6) = ¼ = (3/2)(1 − 5/6) and E_a(5/6) = (1+a)/7 < ¼ ⟺ a < ¾. qa1.py (4):
  A₀(a) = 1.5 for a ≥ 0.75; 1.4919 (a = 0.74); 1.4592 (0.70); 1.3772 (0.60); 1.2944 (0.50, matching [R2] hyp.py). By Prop 4(a),
  any such a gives Burgess-strength subconvexity |L(½,χ)| ≪ |D|^{3/16−η} for all ring class characters.
- **Conditional exponents of [R2] §5** (2/(3−2B♭) = 1.5997 under Liu's half-plane, x^{0.6251} for the ring class PNT) use
  σ ≤ B♭ = 0.8749 < σ_δ: unchanged.
- **Least primes, unconditional:** Theorem 3 keeps the |D|^{ε} loss, so it says nothing within O(1/log|D|) of σ = 1 (no
  log-free gain); the record |D|^{455} is untouched.
- **History check:** Prop 1(b) shows that [C]'s own Carlson argument already gives e(σ) once Y = |D|^{1/(2(3−2σ))} and the
  better of ⌊1/α⌋, ⌈1/α⌉ is used per block; [R2]'s long mollifier is a second route to the same exponent, not a stronger one.

## 6. Lead: one amplifier for the whole family

For fundamental D and χ² ≠ 1, θ_χ ∈ S₁(|D|, χ_D). At p inert in K, λ_{θ_χ}(p) = 0 and λ_{θ_χ}(p²) = λ(p)² − χ_D(p) = 1 for
every χ. For any π in F₀(|D|, χ_D) (real nebentypus), λ_π(p) ∈ iR at such p, so λ_π(p²) = 1 − |λ_π(p)|² ≤ 1: the amplifier
A(π) := Σ_{L<p≤2L, χ_D(p)=−1} λ_π(p²) is maximal, ≍ L/log L, simultaneously on all h dihedral forms of K (given that many inert
primes in (L,2L]). An amplified fourth moment Σ_{π∈F₀(|D|,χ_D)} ω_π|L(½,π)|⁴|A(π)|² ≪ |D|^{1+ε}L, of the kind behind the
reciprocity proofs of subconvexity (Blomer–Khan 2019; Yang, arXiv 2607.04476, Thm A), would give by positivity
Σ_χ|L(½,χ)|⁴ ≪ |D|^{1+ε}(log L)²/L, i.e. (M_a) with a = 1 − θ for L = |D|^{θ} — twice the saving of the sup × mean route
(a = 1 − 2δ with δ = θ/4). In the standard amplification count ¼ − 1/120 corresponds to θ = 1/30, which would improve e(σ)
for σ > 57/58 = 0.9828; reaching the middle range needs θ > ¼, i.e. amplification of Burgess strength (Prop 4). A second caveat:
the amplifier needs ≫ L/log L primes p ∈ (L,2L] with χ_D(p) = −1, which is not known unconditionally for every D when θ is small. Not verified: the admissible amplifier and the error terms of
Yang's Thm A were not checked, and ω_π ≫ |D|^{−ε} is ineffective (Siegel) for dihedral π.

## 7. Novelty and literature (2026-10-21)

- arXiv API, 22 successful queries (axq.out, axq2.out; rate-limited ones rerun): class group × fourth moment / moments / zero density / subconvexity /
  amplification / Lindelöf; mixing conjecture × Heegner / class group; ring class × moment / zero; dihedral × fourth moment.
  Relevant hits: Yang 2607.04476 (uniform subconvexity ¼ − 1/120 for Picard-group characters of orders, and ¼ − 1/80 at t = 0 for
  quadratic central character; full-family fourth moment q^{1+ε}, Thm E; read in TeX); Blomer–Brumley 2009.07093 and
  Blomer–Brumley–Khayutin 2212.06280 (moments of class group twists under GRH); Nordentoft 1903.03932 (hybrid subconvexity via
  Eisenstein sup norms); Vargas 2509.05473 (large values of L(½,χ), lower bounds); Assing–Humphries 2607.13518 (spectral fourth
  moments at Heegner points); Humphries–Khan 1904.05235 (L⁴ of dihedral Maaß forms). None has an unconditional per-D fourth moment
  below |D|^{1+ε} or a per-D density exponent below e(σ).
- Theorem 3's method is the classical "sup bound × mean value" argument (Montgomery, Topics in Multiplicative Number Theory,
  Ch. 12). Its statement for this family was not found; it is routine.

## 8. Dead ends (NOGO candidates)

- Per-D second-moment toolkit (large sieve, k-th powers, first moment, twisted second moment, any mollifier): optimum e(σ) (Prop 1).
- Known fourth moment a = 1 (large sieve; full GL₂ family by positivity): never helps (Thm 2).
- Fourth-moment route below σ ≈ 0.98: needs Burgess-type subconvexity (Prop 4(a)); Halász–Montgomery at the critical block needs
  |D|^{1/16} at σ = 5/6 (Prop 4(b)).
- Halász–Montgomery with the lattice: no Pólya–Vinogradov analogue; no subdivision (all conductors |D|); exceptional sets neutral.
- Moments off the line (Re s = σ₁ > ½): Lindelöf-quality 2k-th moments for σ₁ ≥ 1 − 1/(2k), but class II through them beats the
  first moment at ½ only if σ > 3/2.
- Mollifiers longer than |D|^{½}, or supported on primitive ideals: the first-moment class II does not depend on the length up to
  |D|^{½} and pays for more.
- Sub-family D = f²d_K, d_K fixed: even a Burgess-type bound for the Hecke characters of the fixed field K (δ = 1/16) gives gains
  only for σ > 13/14; 2D Pólya–Vinogradov is of threshold size.

## 9. Numerical checks (qa1.out)

- (1) Worst class I block over α ∈ [0.02, 1], k unrestricted: at α = 1/(3−2σ) with value e(σ) (σ = 0.55, 0.6, 0.7, 0.75, 5/6, 0.9,
  0.95, 0.99; agreement to 4–5 digits).
- (2) Grid optimum over (x, y, class II ∈ {first moment, Lemma 4}): 0.44444, 0.33333, 0.25000, 0.16667, 0.05660 at
  σ = 0.6, 0.75, 5/6, 0.9, 0.97; e(σ) the same.
- (3) σ_δ = 0.999957, 0.999470, 0.991525, 0.928571 for δ = 1/23041, 1/1889, 1/120, 1/16.
- (4) A₀(a) table of §5; (5) 1/(4(3 − 2·5/6)) = 0.1875 = 3/16.

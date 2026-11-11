# Map: Baluyot's proof of Conrey's slope-8/7 density estimate

Source: S. A. C. Baluyot, *On the zeros of Riemann's zeta-function*, PhD thesis, University of Rochester 2017, hdl 1802/32855 (URResearch itemId 33320, file 183460, 200 pp.), chapter 3, printed pp. 93–144 (printed page p = PDF page p + 9). Read from page images, not from pdftotext output.

## Statements

- Mollifier (1.2.6)–(1.2.7): M(s) = Σ_{n ≤ y^{1+ε}} μ(n) n^{−s} W(y/n)/W(y), y = T^θ, with
  W(ξ) = W(ξ; σ) = (1/2πi) ∫_(r) Γ(s+1) ζ(2σ+s)/ζ(2σ) · ζ(1+s) ξ^s ds (σ > ½, ξ > 0, r > 0).
  Termwise: W(ξ) = (ξ ζ(2σ))^{−1} Σ_b b^{1−2σ}/(e^{b/ξ} − 1) > 0; residues: W(ξ) = 1 + A ξ^{1−2σ} + Σ_k c_k ξ^{−k}, A = Γ(2−2σ)ζ(2−2σ)/ζ(2σ) → −1 as σ → ½. So M is a smoothed Selberg mollifier: its weight is W(y/n)/W(y) ≈ (1 − (y/n)^{1−2σ})/(1 − y^{1−2σ}), not Levinson's P(log(y/n)/log y).
- Thm 1.2.2 (p. 17): if 0 < θ < 4/7, there is ε_0 > 0 with (1/T)∫_T^{2T} |ζM(σ+it)|² dt = 1 + O(y^{1−2σ}), uniformly for ½ + 1/log y ≤ σ ≤ ½ + ε_0.
- Thm 1.2.1 (p. 16): if 0 < θ < 4/7, then N(σ, T) ≪_θ T^{1+θ(1−2σ)} log T uniformly for σ ≥ ½ and large T, where N(σ, T) = #{β > σ, 0 < γ ≤ T}. With 2θ = 8/7 − ε this is Conrey's Bull. AMS 20 (1989) Thm 2, as used by Lamzouri.
- Standing assumption (3.1.1): ε_0 < ¼ is small enough that |A| < e²/2 on the σ-range.

## Proof steps (inputs in brackets)

| # | Where | Content |
|---|---|---|
| L1 | Lemma 3.1.1, pp. 93–95 | W(ξ) = 1 + Aξ^{1−2σ} + O(ξ^{−3/4}) (shift to Re s = −¾); (i) W(y) = 1 + O(y^{1−2σ}), (ii) W(y) ≥ ½ + O(y^{−3/4}), using y^{1−2σ} ≤ e^{−2} and (3.1.1), (iii) W(y/h) ≪ e^{−h/2y} for h ≥ 4y, (iv) W(y/h) ≪ 1 [Stirling; ζ ≪ \|t\|², Titchmarsh (5.1.1); style of Soundararajan, Ann. Math. 2000, Lemma 2.1] |
| L2 | Lemma 3.1.2, pp. 95–96 | \|D(0, 2σ−1, 0, H/K)\| ≤ 2K² for the Estermann function D(s,α,β,H/K) = Σ m^{−s−α} n^{−s−β} e(mnH/K) [Hurwitz ζ; Euler–Maclaurin, Titchmarsh (2.1.2)–(2.1.4)] |
| L3 | Lemma 3.1.3, pp. 96–97 | expansion of S(w,α,β,H/K) = Σ m^{−α} n^{−β} e(mnH/K) e(mnw), Im w > 0: residues at s = 1−α, 1−β, 0 plus a Mellin integral of D(s,−α,−β,±H̄/K) [Conrey, Crelle 399 (1989), Lemma 5; Estermann's functional equation] |
| L4 | Lemma 3.1.4, pp. 97–101 | ∫_0^{∞e^{±iλ}} ∫_(c) (1+\|s\|)\|G_j(z,s)\| ds dz ≪ λ^{−η−5/2}(log Δ/Δ)^{η+7/2} for c = 2σ+η, ≪ λ^{2σ−7/2}(log Δ/Δ)^{η+5/2} for c = 2σ+η−1 [proof modelled on Balasubramanian–Conrey–Heath-Brown, Crelle 357 (1985), Lemma 5] |
| L5 | Lemma 3.1.5, p. 102 | incomplete Kloosterman sums: Σ_{r<b≤w,(b,vq)=1} e(nb̄/v) ≪ (vq)^ε v^{1/2}(n,v)^{1/2}(q + (w−r)/v) [Weil, via Hooley 1976, Lemma 4] |
| L6 | Lemma 3.1.6, p. 103 | bilinear Kloosterman fractions: Σ_b Σ_v \|Σ_a Σ_ℓ C(a,ℓ) e(ℓ\overline{ab}/v)\| ≪ (ABNV)^{1/2+ε}{(VB)^{1/2} + (N+A)^{1/4}[VB(N+A)(V+A²) + NA²B²]^{1/4}}, then expanded into 7 monomials A^a N^n (AB)^u V^v [Deshouillers–Iwaniec, Acta Arith. 43 (1984) 305–312, Lemma 1 with ϱ = 1; itself from DI, Invent. Math. 70 (1982), Thm 12 (Kuznetsov)] |
| L7 | Lemma 3.1.7, pp. 104–105 | (1/iΔ√π)∫_(c) e^{(s−s_0)²/Δ²} χ(1+β−s) ω^{−s} ds = ω^{−β}∫_0^∞ z^{s_0−β} e^{−Δ²(log z)²/4}(e(−ωz) + e(ωz)) dz/z [Conrey 1989, Lemma 3; Mellin–Parseval, Titchmarsh, Fourier integrals (2.1.17)] |
| L8 | Lemmas 3.1.8–3.1.9, pp. 106–111 | Gaussian integral; 𝒥 ≪ T^{1−2σ} + T^{−δ} for T/2 ≤ u ≤ 4T, Δ = T^{1−δ} |
| S2 | §3.2, Lemma 3.2.1, pp. 111–114 | ∫_σ^1 (N(ξ,2T) − N(ξ,T)) dξ ≤ (T/4π) log{(1/T)∫_T^{2T}\|ζM\|²} + O(log Ty) [Littlewood's lemma, Titchmarsh (9.9.1); Jensen's theorem for arg, Titchmarsh §9.4; M ≪ y^{1+ε}; Jensen's inequality] |
| S3 | §3.3, Lemma 3.3.1, pp. 114–119 | Gaussian smoothing w(t) of width Δ; g(u) as a contour integral; residue at s = 1; functional equation; L7 with ω = mnh/k, β = 2σ−1; rays at angle ±λ; L3 with w = H(z−1)/K ⇒ g(u) = 𝔐_1 + 𝔐_2 + 𝔇_1 + 𝔇_2 + 𝔈 (3.3.8) |
| S4 | §3.4, pp. 119–134 | Δ = T^{1−δ}, λ = 1/T, c = 2σ+η, η = 2ε (3.4.1), (3.4.39); dyadic split of B(s); Vaughan's identity μ = c_1+c_2+c_3 with cut y^{1/4} [Conrey 1989]; S_1 (type II) by L6: Lemma 3.4.1, L ≪ (ABV)^{x−2σ}N^{2σ−x−1}max{TN, ABV}T^{−1/2}y^{7/8+ε} for y ≤ T^{8/13}; 𝔈_1 ≪ T^{−1/2+7δ/2+ε}y^{7/8} via L4 (3.4.44); S_2 trivially: 𝔈_2 ≪ T^{−1+7δ/2+ε}y^{5/4}; S_3: Lemma 3.4.2 (α_0β_0 ≥ y^{1/4}, as S_1) and Lemma 3.4.3 (α_0β_0 ≤ y^{1/4}, type I by L5): 𝔈_3 ≪ T^{−1/2+7δ/2+ε}y^{7/8} + T^{−1+7δ/2+ε}y^{7/4}; total (3.4.61) |
| S5 | §3.5, pp. 134–140 | 𝔐_1 = ζ(2σ)Σ_{h,k} μ(h)μ(k)(h,k)^{2σ}(hk)^{−2σ}W(y/h)W(y/k)/W(y)² (𝓘 = 1 by the residue at z = 1) = 1/W(y)² + A y^{1−2σ}/W(y)² + 𝔐_12 + O(y^{3/4−2σ+ε}) = 1 + O(y^{1−2σ}) (Euler products F_1, F_1(0,0) = F_1(1−2σ,0) = 1); 𝔐_2 ≪ \|𝒥\|y^{1−2σ} ≪ y^{1−2σ} [thesis Lemma 2.2.1] |
| S6 | §3.6, pp. 140–143 | 𝔇_j ≪ e^{−½log²T}; g(u) = 1 + O(y^{1−2σ}) + O(T^{−1/2+7δ/2+ε}y^{7/8} + T^{−1+7δ/2+ε}y^{7/4}) (3.6.2); choose 7δ + 7θ/4 < 1, then ε_0 with 4θε_0 < δ and 7δ + (7/4 + 4ε_0)θ < 1 ⇒ g(u) = 1 + O(y^{1−2σ}) (3.6.8); integrate over u ⇒ (3.6.11); S2 ⇒ (3.6.12); extend to ½ ≤ σ via N(T) ≪ T log T; Selberg's σ-averaging; dyadic sum; Ingham's T^{3(1−σ)/(2−σ)}(log T)^5 for σ ≥ ½ + ε_0 ⇒ Thm 1.2.1 |
| S7 | §3.7, pp. 143–144 | lower bound with the inner weight v(t) ⇒ the two-sided Thm 1.2.2 (not needed for Thm 1.2.1) |

## Where θ < 4/7 comes from

Both binding error terms come from §3.4. Type II (DI, L6): the monomials (−2,2,4,4), (3,2,3,3) and the corrected (1,2,3,4) give T^{−2}y^{7/2}, and the fourth root is T^{−1/2}y^{7/8}. Type I (Weil, L5): qAV^{3/2} ≤ y^{7/4}. With the G-integral factor T^{−1+7δ/2} from L4 they give T^{−1/2}y^{7/8} and T^{−1}y^{7/4}. The second is the square of the first, so both are o(1) iff θ < 4/7. The side condition y ≤ T^{8/13} comes from the monomial (−1,4,3,4), and S_2 alone would allow θ < 4/5. DI's own Thm 4 (Acta Arith. 43, (5)) has the same term T^{1/2}M^{7/8}, hence M ≤ T^{4/7}.

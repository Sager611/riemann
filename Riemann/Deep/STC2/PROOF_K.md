# (K) by the minimum principle — proof (stc2, 20 Nov 2026)

Notation of Casaña-Giner (Zenodo 23866074): d > ½, B ∌ (I,1), U(b,t) = Σ_ι c_ι G_ι(b,t), G_ι(b,t) = Φ_X (d+b+it)^{−j} + Φ_X (d−b+it)^{−j}
(Φ_R = Re, Φ_I = −Im), f(b,t) = sgn(t)π − 2 arctan(t/(d+b)) − 2 arctan(t/(d−b)) + arctan(t/(2d+b)) + arctan(t/(2d−b)) (4.17).
(K): f ≤ U on 0 ≤ b ≤ ½, t ≠ 0.

**Reduction.** For τ ≥ 0 and s = ±1 put U_s(b,τ) = Σ_ι ε_s c_ι G_ι(b,τ) (ε_s = 1 if s = 1 or X = R, else −1) and
H_s(b,τ) = U_s(b,τ) − s f₊(b,τ), f₊(b,τ) = π − Σ_k α_k arctan(τ/y_k), (α_k) = (2,2,−1,−1), (y_k) = (d+b, d−b, 2d+b, 2d−b).
Since f(b,−τ) = −f(b,τ) and G_ι(b,−τ) = ε_X G_ι(b,τ): (K) at t = sτ (τ > 0) ⟺ H_s(b,τ) ≥ 0.

**Lemma.** H_s ≥ 0 on [−½,½] × (0,∞) if (a) U(b,0) − π ≥ 0 for b ∈ [0,½], (b) H_s(½,τ) ≥ 0 for τ ≥ 0.
Proof. (1) Harmonic: with ζ = b + iτ, (d+b+iτ)^{−j} = (d+ζ)^{−j} is holomorphic and (d−b+iτ)^{−j} = (d−ζ̄)^{−j} antiholomorphic on
|Re ζ| < d; arctan(τ/y_k) = arg(y_k + iτ) = Im log(c_k ± ζ) or −Im log(c_k − ζ) (c_k ∈ {d, 2d}, y_k > 0 for |b| < d). Real and imaginary
parts of holomorphic/antiholomorphic functions are harmonic, so H_s is harmonic on (−d,d) × ℝ ⊃ [−½,½] × [0,∞) (d > ½), and
continuous there (τ = 0 included: f₊ is the analytic form (4.17)). (2) Even in b: b ↦ −b swaps y_1 ↔ y_2, y_3 ↔ y_4 and the two
terms of G_ι. (3) Decay: |U_s| ≤ Σ|c_ι|·2/τ^j ≤ 2Σ|c_ι|/τ (τ ≥ 1), and f₊ = Σ_k α_k arctan(y_k/τ) for τ > 0 (Σα_k = 2), so
|f₊| ≤ Σ|α_k| y_k/τ; hence H_s → 0 as τ → ∞, uniformly in |b| ≤ ½. (4) For (b,τ) in the open half-strip and T > τ, the weak minimum
principle on the rectangle R_T = [−½,½] × [0,T] gives H_s(b,τ) ≥ min over ∂R_T, i.e. ≥ min(inf_{|b'|≤½} H_s(b',0), inf_{0≤τ'≤T} H_s(±½,τ'),
min_{|b'|≤½} H_s(b',T)). The first is U(b',0) − sπ ≥ 0 by (a) (f₊(b',0) = π; for s = −1 it is U + π > 0), the second ≥ 0 by (b) and (2),
the third → 0 as T → ∞ by (3). So H_s(b,τ) ≥ 0. □  (Strict: H_s ≢ 0 and the strong principle give H_s > 0 inside.)

**Certified inputs (certK2.py edge, 80-bit arb; log certK2_edge_all.log).** (a) U(b,0) − π on [0,½]; (b) H_s(½,τ) on [0,100], s = ±1;
(b') τ ≥ 100 at b = ½: Q_s(u) = t²U − 4d³u − K5|u|³ > 0 on 0 < su ≤ 1/100 (u = 1/t, t = sτ), which gives t²U > t²f by Lemma 4.10
(T1 = 100 ≥ 2d + ½; the tail form (4.19) of t²U). Each by 1-D bisection with v(m) − |v'(m)| r − |v''(box)| r²/2 > 0 (exact derivatives
from w ↦ Re Σ Ψ_k(w_k): H_b = Re Σσ_kΨ_k', H_τ = −Im ΣΨ_k', H_bb = −H_ττ = Re ΣΨ_k''; checked against mpmath finite differences).

**Independent check without the lemma (certK2.py box; log certK2_box_all.log).** 2-D bisection of [0,½] × [0,100] for s = ±1 with
H(m) − |H_b(m)| r_b − |H_τ(m)| r_τ − |H_bb(box)| (r_b² + r_τ²)/2 − |H_bτ(box)| r_b r_τ > 0 (Taylor, Lagrange remainder on the convex box), and
region B (|t| ≥ 100) as Q(b,u) > 0 on [0,½] × [−1/100, 1/100] with a second-order form in b and the natural extension in u.


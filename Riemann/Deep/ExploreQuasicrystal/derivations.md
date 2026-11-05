# explore5 derivations (Fourier quasicrystals / Lee–Yang → RH tree)

Notation: LY polynomial = P ∈ C[z_1..z_n], P ≠ 0 on D^n and on (C∖D̄)^n. ℓ ∈ R^n_{>0}, Q-independent. f(x) = P(e^{iℓ_1x}, …, e^{iℓ_nx}). μ_P = Σ_{f(x)=0} δ_x (real-rooted by LY: Im x > 0 ⇒ z ∈ D^n, Im x < 0 ⇒ z ∈ (C∖D̄)^n).

## D1. Kurasov–Sarnak formula and the axis coefficients
On Im x > 0, log P(z) = log P(0) + Σ_{m∈N^n∖0} c_m z^m converges (P ≠ 0 on D^n), so
f'/f(x) = Σ_m i⟨m,ℓ⟩ c_m e^{i⟨m,ℓ⟩x}; the argument principle on the strip |Im x| < ε gives
μ̂_P = (Σ_j d_j ℓ_j / 2π)·δ_0 + Σ_m (coefficient ∝ ⟨m,ℓ⟩c_m) at ±⟨m,ℓ⟩ (KS 2004.05678; d_j = degree in z_j).
ℓ Q-independent ⇒ distinct m give distinct frequencies, so μ̂_P determines every ⟨m,ℓ⟩c_m.

**Lemma (prime-power spectrum rigidity; elementary, ours).** If supp μ̂_P ∖ {0} ⊂ {±kℓ_j : k ≥ 1, j ≤ n}
(for ℓ_j = log p_j: the frequencies log p^k only, as in ζ's explicit formula), then μ_P is a finite union of
arithmetic progressions with differences 2π/ℓ_j.
Proof. c_m = 0 unless m = k e_j, so log P(z) − log P(0) = Σ_j g_j(z_j) on D^n, hence P = P(0)·Π_j Q_j(z_j)
(Q_j = exp g_j is the polynomial P(0,…,z_j,…,0)/P(0)). A product is LY iff every factor is: a root of Q_j inside
(outside) D gives a zero of P in D^n (in (C∖D̄)^n). So every root ω of Q_j is unimodular and
f = 0 ⟺ ∃j, ω: e^{iℓ_jx} = ω, i.e. x ∈ (arg ω + 2πZ)/ℓ_j. ∎
Consequence. ζ's explicit formula has (non-archimedean) spectrum on prime powers only; on any finite prime torus this
forces "Poisson" (periodic) zero sets. All non-periodic structure of ζ's zeros (GUE) must come from mixing primes,
which an FQ can only do through composite frequencies or reflected ones — ζ has neither; it has instead the
archimedean a.c. term (Gamma factor, density log(T/2π)/2π), which has no torus flow. DH: its −f'/f has
coefficients on all n (no Euler product), so the lemma does not even apply — the lemma uses the Euler product,
but it kills rather than proves.

## D2. Finite Lee–Yang toy (Euler product + exact conductor) — RH by Hermite–Biehler
E_P(s) = Π_{p≤P}(p^{s/2} − p^{−s/2}) = (P#)^{s/2}/ζ_P(s), F_P(s) = E_P(s) + E_P(1−s) (= 2 Re E_P on Re s = 1/2).
|p^{s/2} − p^{−s/2}|² = 4(sinh²(σL/2) + sin²(tL/2)), L = log p, is strictly larger at σ than at 1−σ iff σ > 1/2.
So |E_P(s)| > |E_P(1−s)| on Re s > 1/2 and F_P has all zeros on Re s = 1/2 (any P). In torus variables w_p = p^{1/2−s}:
F_P ∝ A(w) + B(w), A = Π(1 − a_p w_p), B = Π(w_p − a_p), a_p = p^{−1/2}: a multi-affine LY polynomial (|B| < |A| on D^n,
Blaschke). Its zero density is θ(P)/2π (constant), its prime-power lines are +Λ(n)/(2π√n) (the sign of 1/ζ: zeros
where arg E_P ≡ π/2 mod π, φ = tθ(P)/2 − arg ζ_P), and log(1 + B/A) adds composite and reflected lines at log(P#/m).
ζ's lines carry −Λ(n)/(2π√n). Zeta-sign version E = (P#)^{s/2}ζ_P(s): per prime the first-order HB margin at the line is
L(1 − 4a cos φ + 3a²)/(1 − 2a cos φ + a²), negative at cos φ = 1 iff a_p > 1/3 ⟺ p ≤ 7, and E^#(s) = E(1−s) has poles on
Re s = 1 (t = 2πk/log p): not HB. The LY mechanism produces the anti-zeta sign; the zeta sign needs poles (the Euler
factor (1 − p^{−1/2}z)^{−1} on each prime axis), i.e. no polynomial on the torus. (Gonek 0704.3448 builds RH-by-
construction models P_X(s) + χ(s)P_X(s̄); the analytic P_X(s) + χ(s)P_X(1−s) is his "wrong AFE": |P_X(1−s)| reaches
exp(X^{σ−1/2}/log X).)

## D3. Hermite–Biehler threshold of the differenced functions (Taylor 1945, Lagarias math/0601653)
F real entire of order 1, F(s) = F(1−s) (F = ξ or the completed DH function Φ), E_h(s) = F(s+h), h > 0.
E_h^#(s) := conj E_h(1 − s̄) = F(s−h). A_h = (E_h + E_h^#)/2, B_h = (E_h − E_h^#)/(2i).
(⇐) If F ≠ 0 on Re s > 1/2 + h: pair the zeros ρ, 1−ρ̄ (b = Re ρ − 1/2, |b| ≤ h, x = σ − 1/2 > 0, y = t − Im ρ):
|E_h/E_h^#|² factor = [((x−b+h)²+y²)((x+b+h)²+y²)] / [((x−b−h)²+y²)((x+b−h)²+y²)], and numerator − denominator
= 4xh·2(x² + h² − b² + y²) > 0. The exponential factors of the Hadamard product cancel. So E_h is Hermite–Biehler and
A_h, B_h have only real (critical-line), simple, interlacing zeros (de Branges / HB theorem).
(⇒) If F(ρ0) = 0 with Re ρ0 > 1/2 + h, E_h vanishes at ρ0 − h in Re s > 1/2: not HB, so A_h, B_h cannot both be
real-rooted with interlacing zeros (HB theorem for real entire functions of order 1 with real zeros).
Hence sup Re ρ − 1/2 = inf{h > 0 : A_h, B_h real-rooted and interlacing}: a zero-free half-plane in disguise (R47).
Corollary (routine): for ξ, h > B_2 − 1/2 = 0.37338196661320752 (Lean-checked half-plane, Bebel; 3/8 from 7/8)
suffices; Lagarias' unconditional range is |h| ≥ 1/2. A_h alone can be real-rooted slightly below the threshold
(t3: DH local threshold 0.2981 against Re ρ0 − 1/2 = 0.3085: the bubble {|E_h| < |E_h^#|} around ρ0 − h attaches to
the line once ρ0 − h is within ≈ 0.010 of it).

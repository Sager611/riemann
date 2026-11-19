# FodgeAudit — FULLPROOF: the chain of Appendix C (Zenodo 24325329 v40), re-derived

Source: `zfr_exact_stechkin_v40_2026-11-17.pdf`, `pdftotext -layout` lines (l.) and PDF pages (p.). Page starts: p.30 l.1972, p.31 l.2044, p.32 l.2134, p.33 l.2237, p.34 l.2307, p.35 l.2371, p.36 l.2436. Formulas read from page images (130 dpi), not from pdftotext.

## Setup (p.30, l.~2003–2010)
f(x) = (t/2π) log(N+x), W = π(N+L)³/t, M = ⌈cW^{1/3}⌉, g(x) = f(x+m) − f(x), 1 ≤ m < L; t > 0, N ≥ 1, L ≥ 2 integers, c > 0.

## Lemma C.2 (weighted Weyl differencing; p.30, l.1979–1995) — CORRECT
Statement: ξ_1..ξ_L ∈ C, M ≥ 1, w_j ≥ 0, W1 = Σw_j > 0, r(m) = Σ_{j=0}^{M−1−m} w_j w_{j+m}:
|Σξ_n|² ≤ ((L+M−1)/W1²)(r(0)Σ|ξ_n|² + 2Σ_{m=1}^{M−1} r(m)|Σ_{n=1}^{L−m} ξ_{n+m}ξ̄_n|).
- W1Σξ_n = Σ_{k=2−M}^{L}Σ_j w_jξ_{k+j}: each n ∈ [1,L] is hit by k = n−j ∈ [n−M+1, n] ⊂ [2−M, L] for every j: total weight W1. ✓
- Cauchy–Schwarz over the L+M−1 values of k. ✓ Expansion: for fixed (j, j′), n = k+j′ runs over [2−M+j′, L+j′] ⊇ [1, L] (0 ≤ j′ ≤ M−1), so the inner sum is the full Σ_n ξ_{n+j−j′}ξ̄_n. ✓ Σ_{j−j′=m} w_jw_{j′} = r(|m|). ✓
- Terms m and −m are complex conjugates; 2Re ≤ 2|·| needs r(m) ≥ 0, i.e. w_j ≥ 0. ✓ m ≥ L: inner sum empty. ✓
- w_j = 1: r(m)/W1² = (M−m)/M² = (1−m/M)/M ✓. Closed forms for (j+1)(M−j): exact for M ≤ 300 (R482) and re-proved symbolically here (`c2_sym.py`).

## Lemma C.3 (explicit B-process; statement p.30 l.2012–2022)
Hypotheses N ≥ 10⁴, L ≤ N/75, N′ = N − 1/2, µ′ = (2(N+L)^{3/2}/L)(N′^{−1/2} − (N′+L)^{−1/2}); for 1 ≤ m < L′ ≤ L:
|Σ_{n=0}^{L′−1−m} e(g(n))| ≤ µ′L(m/W)^{1/2} + 4.524(W/m)^{1/2} + 1.0004 + 2.24 m/W.

### Step 0, parameters (p.30) — CORRECT
a = −1/2, b = L′−m−1/2, λ = m/W, h = ((N+L)/N′)³, ε3 = 3/N′, ε4 = 12/N′².
- h ≤ (76/75 · 10⁴/(10⁴−½))³ = 1.04069… ≤ 1.0407 (computed in `consts.py`). ✓
- g″(x) = (t/2π)[(N+x)^{−2} − (N+x+m)^{−2}] = (t/π)∫_x^{x+m}(N+y)^{−3}dy ✓; y ∈ [a, b+m] = [−½, L′−½] ⇒ N+y ∈ [N′, N+L−½] ⊂ [N′, N+L] ✓ ⇒ λ ≤ g″ ≤ hλ ✓.
- g‴ = −(3t/π)∫(N+y)^{−4}, |g‴| ≤ ε3 g″ ✓; g⁗ = (12t/π)∫(N+y)^{−5}, |g⁗| ≤ ε4 g″ ✓. g″ > 0, g‴ < 0: g′ increases, g″ decreases ✓. (α, β < 0: g′ < 0.)

### Step 1, Poisson summation (p.30) — CORRECT
Euler: Σ_{a<n≤b}F(n) = ∫_a^bF + ∫_a^bρF′ − ρ(b)F(b) + ρ(a)F(a), ρ = {x}−½; ρ(a) = ρ(b) = 0 at half-integers; b ∉ Z so a<n<b; the integers in (a,b) are 0..L′−1−m ✓.
ρ = −Σ_{j≥1} sin(2πjx)/(πj), partial sums uniformly bounded ⇒ termwise integration (dominated convergence); by parts ∫sin(2πjx)F′ = −2πj∫cos(2πjx)F (sin(2πjx) = 0 at half-integers) ⇒ ∫ρF′ = Σ_{j≥1}(I_j + I_{−j}) ✓.

### Step 2, exceptional frequencies (p.30) — CORRECT
j_α = ⌊α+½⌋ ⇒ α ∈ [j_α−½, j_α+½), so j ≠ j_α ⇒ |j−α| ≥ ½ ✓ (same for β).
|I_j| ≤ 4(πλ)^{−1/2} for every j: {|ψ′_j| ≤ ρ0} is an interval of length ≤ 2ρ0/λ (ψ″_j = g″ ≥ λ); on each of ≤ 2 complementary intervals, ψ′_j monotone of one sign, |ψ′_j| ≥ ρ0: |∫e(ψ)| ≤ (1/2π)(1/|ψ′(c)| + 1/|ψ′(d)| + |1/ψ′(d) − 1/ψ′(c)|) ≤ 1/(πρ0) ✓. ρ0 = (λ/π)^{1/2}: 2ρ0/λ + 2/(πρ0) = 4(πλ)^{−1/2} ✓.

### Step 3, frequencies without a stationary point (p.31 top) — CORRECT
j ∉ E, j < α or j > β: ψ′_j ∈ [α−j, β−j] of one sign, |ψ′_j| ≥ d_j ≥ ½ ✓. Two partial integrations with B_j = e(ψ_j)/(2πiψ′_j):
I_j = B_j(b) − B_j(a) + R_j, R_j = −[e(ψ)ψ″/(4π²ψ′³)]_a^b + (1/4π²)∫e(ψ)(ψ″/ψ′³)′ ✓ ⇒ the displayed 4π²|R_j| ≤ g″(a)/|α−j|³ + g″(b)/|β−j|³ + ∫|(g″/ψ′³)′| ✓.
(g″/ψ′³)′ = g‴/ψ′³ − 3g″²/ψ′⁴; with dψ′ = g″dx: ∫|g‴|/|ψ′|³ ≤ ε3∫_{d}^{∞}v^{−3}dv = ε3/(2d²) ✓, 3∫g″²/ψ′⁴ ≤ 3hλ∫_d^∞v^{−4} = hλ/d³ ✓.
Sum: for j > β, |β−j| = d_j ≤ |α−j|, so the hλ-terms ≤ 3hλ/d_j³; d_j ∈ ½+Z≥0 at least termwise: Σd^{−3} ≤ Σ(k+½)^{−3} = 7ζ(3), Σd^{−2} ≤ π²/2 ⇒ (1/4π²)(21ζ(3)hλ + π²ε3/4) per side ✓, both sides 21ζ(3)hλ/(2π²) + ε3/8 ✓.

### Step 4, frequencies with a stationary point (p.31) — CORRECT
y = x−x_j; ψ′ = yP(y), ψ−ψ(x_j) = ½y²q(y), P, q ∈ [λ, hλ] (Taylor, integral remainder) ✓. u = y q^{1/2}: u du/dy = yP ⇒ du/dy = P/q^{1/2} > 0 ✓; a = q^{1/2}/P, a0 = g″(x_j)^{−1/2} ≤ λ^{−1/2} ✓.
∫_U^∞e(u²/2) = −e(U²/2)/(2πiU) + τ(U), τ = (1/2πi)∫_U^∞u^{−2}e(u²/2); second partial integration: |∫_U^∞u^{−2}e| ≤ U^{−3}/(2π) + (3/2π)∫_U^∞u^{−4} = U^{−3}/π ⇒ |τ| ≤ 1/(2π²U³) ✓.
a = a0 + u b1(u); the a0-part gives a0e(1/8) + a0[e(u_b²/2)/(2πiu_b) − e(u_a²/2)/(2πiu_a)] − a0(τ(u_b)+τ(|u_a|)); the u b1-part by parts gives [b1 e(u²/2)/(2πi)] − (1/2πi)∫e(u²/2)b1′. Boundary terms combine to e(u²/2)a(u)/(2πiu) ✓, and ψ′(x) = (dψ/du)(du/dx) = u/a(u) ✓ ⇒ B_j(b) − B_j(a) ✓. |R_j| bound ✓.
u_b ≥ (β−j)/(hλ)^{1/2}: u_b²/2 = ∫_{x_j}^bψ′dx = ∫_0^{β−j}v dv/ψ″ ≥ (β−j)²/(2hλ) ✓ (a naive y_b-argument gives only (β−j)/(hλ^{1/2}); the stated one is right). Sum over stationary j (β−j, j−α ∈ ½+Z≥0 termwise): (h^{3/2}λ/2π²)·14ζ(3) = 7ζ(3)h^{3/2}λ/π² ✓.
b1′(u) = (ua′(u) − a(u) + a0)/u² = u^{−2}∫_0^u s a″(s)ds ✓ ⇒ |b1′| ≤ ½ sup|a″| ✓. a″ (in u) = a(a_y² + a a_yy) ✓ (d/du = a d/dy).
|P′| ≤ ε3hλ∫s = ε3hλ/2, |q′| ≤ 2ε3hλ∫s(1−s) = ε3hλ/3, |P″| ≤ ε4hλ/3, |q″| ≤ 2ε4hλ∫s²(1−s) = ε4hλ/6 ✓.
a_y = q′/(2q^{1/2}P) − q^{1/2}P′/P² ⇒ |a_y| ≤ ε3λ^{−1/2}(h/6 + h^{3/2}/2) ✓.
a_yy = q″/(2q^{1/2}P) − q′²/(4q^{3/2}P) − q′P′/(q^{1/2}P²) − q^{1/2}P″/P² + 2q^{1/2}P′²/P³ (five terms) with bounds λ^{−1/2}·{ε4h/12, ε3²h²/36, ε3²h²/6, ε4h^{3/2}/3, ε3²h^{5/2}/2} ⇒ |a_yy| ≤ λ^{−1/2}[ε3²(7h²/36 + h^{5/2}/2) + ε4(h/12 + h^{3/2}/3)] ✓; |a| ≤ h^{1/2}λ^{−1/2} ⇒ |a″| ≤ h^{1/2}λ^{−3/2}Ξ with the stated Ξ ✓.
du/dx = 1/a = P/q^{1/2} ≤ hλ^{1/2} ✓; #stationary j ≤ β−α ≤ hλ(b−a) (j ∈ [α+½, β−½]) ✓; second part ≤ Σ_j (1/2π)·hλ^{1/2}(b−a)·½h^{1/2}λ^{−3/2}Ξ ≤ h^{5/2}Ξ(b−a)²/(4π) ✓, b−a = L′−m < L ✓; numeric value in `consts.py`.

### Step 5, main terms (p.31) — CORRECT
F(v) = g″(x)^{−1/2}, g′(x) = v: increasing ✓ (g″ decreasing). F(j) − ∫_{j−½}^{j+½}F ≤ ½(F(j) − F(j−½)) for increasing F; disjoint [j−½, j+½] ⊂ [α, β] ⇒ ΣF(j) ≤ ∫_α^βF + ½(F(β)−F(α)) ✓.
∫_α^βF = ∫_a^b g″^{1/2} ✓; g″ ≤ (t/π)m(N+x)^{−3} ⇒ ≤ (mt/π)^{1/2}·2(N′^{−1/2} − (N′+L)^{−1/2}) = µ′L(m/W)^{1/2} exactly ✓ (b < L−½).
F(β)−F(α) ≤ λ^{−1/2}(1 − h^{−1/2}) ≤ (h^{1/2}−1)λ^{−1/2} ✓ (stated form weaker, valid).

### Step 6, boundary terms (p.31) — CORRECT
e(−jb) = e(−ja) = (−1)^j at half-integers ✓; Σ_{k∈Z}(−1)^k/(ε−k) = π/sin πε (symmetric) ✓; ε = β−j_β ∈ [−½, ½); |π/sin πε − 1/ε| ≤ π−2 (odd, increasing in |ε|) ✓; j_α ≠ j_β ⇒ |β−j_α| ≥ ½ ⇒ |Σ_{j∉E}B_j(b)| ≤ (π−2+2)/(2π) = ½ ✓; same at a ✓. Shifting the centre of symmetric partial sums changes them by o(1) ✓.

### Step 7, collecting (p.31 bottom) — CORRECT
|Σ| ≤ 2·4(πλ)^{−1/2} (E) + µ′L(m/W)^{1/2} + ½(h^{1/2}−1)λ^{−1/2} + ½ + ½ + [21ζ(3)h/(2π²)]λ + ε3/8 + [7ζ(3)h^{3/2}/π²]λ + 0.00027 = the displayed line ✓. Constants: `consts.py` (256 bits).

## Lemma C.4 (p.32, l.2137–2177) — CORRECT
Statement: N ≥ 10⁴, L ≤ N/75, 1 ≤ L′ ≤ L: |Σ_{n<L′}e(f(n))|² ≤ (L/W^{1/3} + c)(αL + βW^{2/3}), α = 1/c + (8µ′/15)(c+W^{−1/3})^{1/2}, β = 8·4.524/(3√c) + 1.0004/W^{1/3} + 2.24(c+W^{−1/3})/(3W).
- C.2 with w_j = 1 (flat weights: the (1−m/M) factors), ξ_n = e(f(n−1)): inner sum Σ_{n=1}^{L′−m}e(f(n+m−1)−f(n−1)) = Σ_{n=0}^{L′−1−m}e(g(n)), the sum of C.3 with the same N, L ✓; empty for m ≥ L′ (C.3 bound ≥ 0 so extending to all m < M is harmless) ✓.
- Σ(1−m/M)m^{−1/2} ≤ ∫_0^M = (4/3)M^{1/2} (decreasing) ✓; Σ(1−m/M)m^{1/2} ≤ (4/15)M^{3/2} (concave, zero at 0 and M) ✓; Σ(1−m/M)m = (M²−1)/6 ✓; Σ(1−m/M) = (M−1)/2 ✓ ⇒ the displayed ratio bound (with 1.0004(M−1)/M ≤ 1.0004) ✓.
- L′+M−1 < L + cW^{1/3} ✓; L′/M ≤ L/(cW^{1/3}); M^{1/2} ≤ W^{1/6}(c+W^{−1/3})^{1/2}; M^{−1/2} ≤ (cW^{1/3})^{−1/2}; M ≤ W^{1/3}(c+W^{−1/3}) ⇒ bracket ≤ W^{−1/3}(αL + βW^{2/3}) term by term ✓. L′ = 1: 1 ≤ (L′+M−1)L′/M = 1 ✓.
- NB: the weighted C.2 ((j+1)(M−j), "−1.25%") is NOT used by C.4/C.1; only flat weights.

## Proposition C.1, proof (p.32–33, l.2178–2283)
c = 2.44, r = 80, ν = 0.95, n1 = ⌊√(t/2π)⌋, T1 = 2.9·10¹².
- RS input: |ζ(½+it)| ≤ 2|Σ_{n≤n1}n^{−1/2+it}| + 1.48t^{−1/4} + 0.127t^{−3/4}, t ≥ 200, from [7, Lemma 2.1] (= HPY 2207.02366, checked below) ✓; < 2 trivially at t ≥ T1 ✓.
- Blocks: N0 = ⌈νt^{1/3}⌉, N_{k+1} = N_k + L_k, L_k = ⌈N_k/r⌉ ≤ N_k(1/r + 1/N_k) ≤ u⁺N_k with u⁺ = 1/r + 1/(νT1^{1/3}) (N_k ≥ N0 ≥ νT1^{1/3}) ✓; last block starts ≤ n1, cut at n1, uncut end N_end ≤ (1+u⁺)n1 ✓.
- Head: 2Σ_{n<N0}n^{−1/2} ≤ 2(2√(N0−1) − 1) ✓; N0−1 < νt^{1/3} ✓ ⇒ ≤ 4√ν t^{1/6} − 2 ✓.
- Per block: 1/r ≤ u_k ≤ u⁺ = 0.0125738… < 1/75 ✓, N_k ≥ νT1^{1/3} = 13547.3 > 10⁴ ✓ ⇒ C.4 applies; x_k = W_k^{1/3} = π^{1/3}N_{k+1}t^{−1/3} ≥ y0 = π^{1/3}ν(1+1/r) ✓.
- Partial summation (n^{−1/2} decreasing): |Σ_block| ≤ N_k^{−1/2} max_{L′≤L_k}|Σ_{n′<L′}e(f(n′))| ✓ (covers the cut last block).
- (L/x + c)(αL + βx²) = (L²/x)(α + E), E = βx²/L + cαx/L + cβx³/L² ✓; L²/(Nx) = π^{−1/3}t^{1/3}u²/(1+u) ✓ ⇒ 2·block ≤ 2t^{1/6}π^{−1/6}(u/√(1+u))(α+E)^{1/2} ✓.
- µ′ = ((N+L)/(N′+L))^{3/2}φ(L/N′), φ(v) = 2(1+v)/(1+√(1+v)) — algebra re-derived ✓; φ increasing ✓; (N+L)/(N′+L) ≤ 1 + 1/(2N0′), L/N′ ≤ u⁺(N0′+½)/N0′ with N0′ = νT1^{1/3} − ½ ≤ N_k − ½ ✓ ⇒ µ′_k ≤ µ′⁺ ✓.
- α_k ≤ ᾱ(x_k) = 1/c + A(c+1/x_k)^{1/2}, A = 8µ′⁺/15, decreasing to α∞ = 1/c + A√c ✓; β_k ≤ β̄ (x_k ≥ y0) ✓.
- E_k ≤ (e1x_k + e0)t^{−1/3}: 1/L_k ≤ π^{1/3}(1+u⁺)r/(x_kt^{1/3}); βx²/L → β̄π^{1/3}(1+u⁺)r·x; cαx/L → cᾱ(y0)π^{1/3}(1+u⁺)r (ᾱ decreasing); cβx³/L² ≤ cβ̄π^{2/3}(1+u⁺)²r²x t^{−2/3} ≤ (…)T1^{−1/3}x t^{−1/3} ✓.
- (α_k+E_k)^{1/2} ≤ ᾱ^{1/2} + E_k/(2ᾱ^{1/2}) ≤ ᾱ(x_k)^{1/2} + E_k/(2α∞^{1/2}) ✓ (concavity; ᾱ ≥ α∞).
- u/√(1+u) ≤ κlog(1+u), κ = (2+u⁺)/(2√(1+u⁺)): log mean ≤ arithmetic mean, (2+u)/(2√(1+u)) increasing ✓.
- Σ_k (u_k/√(1+u_k))ᾱ(x_k)^{1/2} ≤ κΣ∫_{N_k}^{N_{k+1}}ᾱ(π^{1/3}Nt^{−1/3})^{1/2}dN/N (x_k ≥ π^{1/3}Nt^{−1/3} on the block, ᾱ decreasing) ✓ ≤ κ(α∞^{1/2}((1/6)log t + C′) + I∞), C′ = log((1+u⁺)/(√(2π)ν)), I∞ = ∫_{π^{1/3}ν}^∞(ᾱ^{1/2} − α∞^{1/2})dy/y: y-range ⊂ [π^{1/3}ν, π^{1/3}(1+u⁺)t^{1/6}/√(2π)] ✓.
- (p.33 top) Σu_kx_k = π^{1/3}t^{−1/3}Σ(1+u_k)L_k ≤ π^{1/3}t^{−1/3}(1+u⁺)²n1 (ΣL_k ≤ N_end ≤ (1+u⁺)n1) ✓; Σu_k ≤ κ″Σlog(1+u_k) = κ″log(N_end/N0) ≤ κ″((1/6)log t + C′), κ″ = u⁺/log(1+u⁺) (u/log(1+u) increasing) ✓.
- Collect: E-part ≤ (t^{1/6}/(π^{1/6}α∞^{1/2}))t^{−1/3}Σu_k(e1x_k + e0) (1/√(1+u_k) ≤ 1) ⇒ S's two terms ✓; ((1/6)log t + C′)t^{−1/6} decreasing for t ≥ T1 since (1/6)log T1 + C′ = 3.93 > 1 ✓. Head + RS: 4√νt^{1/6} − 2 + 2 ✓ ⇒ b1 = κα∞^{1/2}/(3π^{1/6}), b2 = 4√ν + (2κ/π^{1/6})(α∞^{1/2}C′ + I∞) ✓, S as printed ✓.
- SIGN NOTE: C′ = log((1+u⁺)/(√(2π)ν)) ≈ −0.855 < 0, so b2 decreases in α∞: an upper bound for α∞ plugged into b2 is the wrong direction. Correct only if the code evaluates b2 as one ball (checked in the code audit).
- I∞ tail: ᾱ(y) − α∞ = A((c+1/y)^{1/2} − c^{1/2}) ≤ A/(2√c y) ✓ (concavity); ᾱ^{1/2} − α∞^{1/2} ≤ (ᾱ−α∞)/(2α∞^{1/2}).
- Variant (p.33): ν = 1.15, T1 → e^{27.4}: N0 ≥ 1.15e^{9.133} = 10651 > 10⁴ ✓, u⁺ = 0.0125 + 1/10651 < 1/75 ✓; claimed 0.3080571, 2.428866, 485.31 (used in Lemma 8.9).

## Lemma C.5 (p.33, l.2284–2297) — CORRECT as a restatement
C.3's proof uses g only via: g ∈ C⁴, g″ non-increasing, λ ≤ g″ ≤ hλ, |g‴| ≤ ε3g″, |g⁗| ≤ ε4g″, a, b ∈ Z+½; N, L enter only via h, ε3, ε4, b−a ≤ L and ∫g″^{1/2} ≤ µ′L(m/W)^{1/2} ✓ (C⁴ needed: a_yy involves ψ⁗) ✓.

## Lemmas C.6 (l.2298–2334), C.8 (l.2344–2390) — CORRECT
- C.6: g̃ = g(x+k) − g(x); −g̃″ = (3t/π)∫_0^m∫_0^k(N+x+s1+s2)^{−4} ✓ (f⁗ = −(3t/π)(N+x)^{−4}); on [−½, L′−m−k−½] the argument lies in [N′, N+L] ✓ ⇒ C.5 for −g̃ with λ = mk/V, V = π(N+L)⁴/(3t), h = ((N+L)/N′)⁴, ε3 = 4/N′, ε4 = 20/N′² ✓; −g̃″ decreasing ✓; ∫(−g̃″)^{1/2} ≤ (3tmk/π)^{1/2}∫_{−½}^{L−½}(N+x)^{−2} = µ″L(mk/V)^{1/2} ✓ (µ″ = (N+L)²/(N′(N′+L))). Flat C.2 over k < K_m with C.4's four sums ⇒ |S_m|² ≤ (L+K_m−1)Q_m ✓; |S_m| ≤ L trivially ✓. Constants 4.52695, 1.00051, 2.27161 (N ≥ 5·10⁴) recomputed in `consts.py`.
- C.8: ĝ = g̃(x+j) − g̃(x), ĝ″ = (12t/π)∫∫∫(·)^{−5} ✓; λ = mkj/V3, V3 = π(N+L)⁵/(12t), h = ((N+L)/N′)⁵, ε3 = 5/N′, ε4 = 30/N′² ✓; ∫ĝ″^{1/2} ≤ (12tmkj/π)^{1/2}∫(N+x)^{−5/2} = µ3L(mkj/V3)^{1/2} ✓. Constants 4.53039, 1.00085, 2.30835 (N ≥ 2·10⁴) recomputed in `consts.py`.

## Proposition C.7 (l.2335–2343, proof l.2391–2447) — proof CORRECT; code `propC1p_v6.py` matches it
- Blocks N_k = ⌈ν_a sρ^k⌉, s = t^{1/3}, ρ = 1+1/80; block k occurs only if s ≥ 2πν_a²ρ^{2k} (N_k ≤ n1 ≤ s^{3/2}/√(2π)) ✓; y_k ∈ [ν_aρ^k, ν_aρ^k + 1/s_k], l_k ∈ [g0 ∓ 1/s_k], g0 = ν_aρ^k(ρ−1) ✓; N_k ≥ 1.001·2·10⁴ ✓; l_k/y_k ≤ 1/75 asserted ✓.
- Every term of B_k/s² non-increasing in s for fixed y, l and integers: (L+M−1)/s = l + (M−1)/s; ρ_m/s, σ_m/s, τ_{m,k}/s term by term (W fixed, V ∝ s, V3 ∝ s²: √V3/s constant); µ′, µ″, µ3 decrease in s (re-derived: µ3 ∝ (y−δ)^{−3/2} − (y+l−δ)^{−3/2}, δ = 1/(2s), increasing in δ) ✓ ⇒ evaluation at s_k with (y, l) in the balls is an upper bound ✓.
- Code: µ′ written as 2(N+L)^{3/2}/(√N′√(N′+L)(√N′+√(N′+L))) (= µ′, no cancellation) ✓; m-groups use increasing terms at the group's top, decreasing at its bottom ✓; group weight sums exact (flat closed form ✓; parabola via parab_forms) ✓; min(L, ·) with up(L) ✓; methods gated by lower bounds N ≥ 10⁴ / 5·10⁴ / 2·10⁴ ✓; term-by-term 2l/√y ✓; head 4√ν_a, RS 2/√s_a ✓; U rounded up ✓; any free integers (M, K, J; truncated last group reuses J) admissible ✓.
- Coverage: run_table_v6.py chains intervals by identical 4-decimal endpoint strings; prof_table_v6.txt: 421 contiguous intervals, 28.6900 → 87.2500, no gap or overlap ✓ (outer ends parsed as floats: a sliver < 10⁻¹⁴ at 28.69/87.25 is covered by C.1's e^{27.4} variant and C.1 itself).

## Cited B-process [4] (Dhiman–Kadiri–Quesada-Herrera, arXiv 2609.00537) — not an input
Fodge cites [4] only as another explicit B-process ("also given by"); C.3's proof is self-contained and does not use it. DKQH's Corollary (TeX l. 297–312, a, b ∈ Z+½, 0 < f′(b) < 1, f′ decreasing) has error 2.686/√λ2 + 2·3^{2/3}π^{−2/3}h2h3^{1/3}(b−a)λ3^{1/3} + (2/π)log(f′(a)−f′(b)) + D_f; Fodge instead sums the boundary terms in closed form (π/sin πε, |·| ≤ ½ per end), so no log term, at the price of 8/√π for the two exceptional frequencies. Different argument; no transcription to check.

## The chain (source lines) and what feeds N9/N7
- N9 record = Prop. C.1 (p.30 l.1975; proof p.32–33 l.2178–2283) ← HPY 2207.02366 Lemma 2.1 (TeX l.150–163, verbatim: t ≥ 200, n1 = ⌊√(t/2π)⌋, R(t) = 1.48t^{−1/4} + 0.127t^{−3/4}; its Gabcke input audited in RevAudit) + head bound + blocks [N_k, N_k+⌈N_k/80⌉) from ⌈0.95t^{1/3}⌉ ← Lemma C.4 (l.2137–2177) ← Lemma C.2 with w_j = 1 (l.1979–1995) + Lemma C.3 (l.2012–2136) + constants (`appC_sp_consts.py`, log `appC_sp_consts.log`).
- Variant from e^{27.4} (p.33, ν = 1.15; `appC_sp_consts_v2.py … e^{27.4}`, log `appC_sp_consts_v2_274.log`) → Lemma 8.9 / N7 below T1.
- N7 R₀ = 4.6892 (l.1316–1323, 1573–1578) ← Prop. C.1, its variant, Prop. C.7 (p.34–35) ← C.2 (flat or parabola weights) + C.3 / C.6 ← C.5 / C.8 ← C.5; code `propC1p_v6.py`, table `prof_table_v6.txt`. With C.1 alone: 4.7079.
- Every arrow re-derived in this file; numerical evidence in REPORT.md §2–§5.

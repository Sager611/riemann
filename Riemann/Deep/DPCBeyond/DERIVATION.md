# DPCBeyond derivation log

## D1 (13:55) Key idea: pointwise zero bound for R_T (replaces the 4th moment)

Notation: X = T^θ, L = log T, z = ρ − 1/2 − it − (shift), d_ρ = β − 1/2, d⁺ = max(d, 0), a_ρ = d_ρ L.
x_p(t) = ρ_η(log p/L) p^{−1/2−it}. c_n n^{−it} = (−1)^j f(Σ log p_i/L) ∏ x_{p_i}(t) for n = p_1⋯p_j squarefree.

D1.1 Explicit formula (Guinand–Weil, Perron shift; φ ∈ C_c^∞((0,∞))): with Φ(z) = ∫ φ(u) e^{zu} du,
  Σ_n Λ(n) φ(log n) n^{−s} = Φ(1 − s) − Σ_ρ Φ(ρ − s) − Σ_{k≥1} Φ(−2k − s)   (all zeros, multiplicity).
  Take φ(u) = ρ_η(u/L) χ(u/L)/u (χ smooth, = 1 on [0, θ], 0 beyond θ + η'): Λ(p)φ(log p) = ρ_η(log p/L).
  Q(τ) := Σ_n Λ(n)φ(log n) n^{−1/2−iτ} = Σ_p x_p(τ) + Q_pp(τ) (prime powers, Σ|coeff| = O_η(1)),
  Q(τ) = −Σ_ρ Φ(ρ − 1/2 − iτ) + O(T^{−10}) for τ ∈ [T/4, 4T] (pole term Φ(1/2 − iτ) ≪_A X^{1}(TL)^{−A}).
  NEEDS: ρ_η smooth (check DPC l.676–684).

D1.2 Free j-fold sums. S_j(t) := Σ_{p_1..p_j} f(Σ log p_i/L) ∏ x_{p_i}(t) = ∫ f̂(ξ) (Q − Q_pp)(t − 2πξ/L)^j dξ.
  Main part: Σ_{ρ_1..ρ_j} (−1)^j I_f(ρ_1..ρ_j; t), I_f := ∫_{R_+^j} f(Σu_i/L) ∏ φ(u_i) e^{z_i u_i} du, z_i = ρ_i − 1/2 − it.
  On supp: u_i ≥ ηL/2 > 0, Σu_i ≤ θL (f = 0 on [θ,∞)), so |e^{Σ z_i u_i}| = e^{Σ d_i u_i} ≤ X^{max_i d_i⁺}.
  IBP A times in each u_i with L|z_i| ≥ 1 (integrand vanishes smoothly at the ends; ∂^m φ ≪ L^{−m−1}):
    |I_f| ≤ C_{A,j} X^{max d_i⁺} ∏_i κ(L|z_i|),  κ(x) := min(1, x^{−A}),  |z_i| ≥ |γ_i − t|.
  Prime-power and distinctness (partition) corrections: blocks of size ≥ 2 are absolutely summable (Σ_p 1/p over
  (T^{η/2}, T^θ) = O_η(1)); for fixed block primes the free part is the same object with f shifted left (bound only improves).

D1.3 LEMMA Z (pointwise). With On(t) := Σ_ρ κ(L|γ_ρ − t|), W(t) := Σ_ρ X^{d_ρ⁺} κ(L|γ_ρ − t|), for t ∈ [T/2, 3T]:
  |R_T(t)| ≤ C (1 + W(t)) (1 + On(t))^{k−1},     (X^{max d⁺} ≤ Σ_i X^{d_i⁺}; Σ_tuples ≤ j W On^{j−1})
  and for 0 ≤ y ≤ 1/2 (s = 1/2 ± y + it changes Re z_i to d_i ∓ y):
  |R_T(t − iy)| ≤ C (1 + W_{−y}(t))(1 + On(t))^{k−1},  W_{−y} := Σ_ρ X^{(d_ρ − y)⁺} κ ≤ On + X^{−y} W,
  |R*(t − iy)| = |R_T(t + iy)| ≤ C X^{y} (1 + W(t))(1 + On(t))^{k−1}.
  C depends on f, η, θ, A, k only. No RH; uses only the explicit formula.

## D2 (14:00) Moments from Lemma Z

D2.1 Depth bins B_a = {ρ: γ ∈ [T/2, 5T/2], d_ρ L ∈ [a, a+1)}, |B_a| ≤ C TL e^{−κa} (density, uniform σ ≥ 1/2;
  κ = 8/7 − δ Conrey/Baluyot, or 1 − δ Jutila). W ≤ e^θ On + Σ_{a≥1} e^{θ(a+1)} F_a, F_a := Σ_{ρ∈B_a} κ(L|γ − t|).
D2.2 ‖On‖_{L^p[T,2T]} ≪_p T^{1/p}: On ≤ Σ_m (1+|m|)^{−A} n_m(t), n_m = zeros in [t + m/L, t + (m+1)/L) (Fujii, Step 5b).
D2.3 ‖F_a‖_p ≪ T^{1/p} e^{−κa(1−1/w)/p}: n^{(a)}_m ≤ n_m integer, (n^{(a)}_m)^p ≤ 1_{n^{(a)}_m ≥ 1} n_m^p, Hölder:
  ∫ ≤ meas{n^{(a)}_m ≥ 1}^{1−1/w} ‖n_m‖_{pw}^p ≤ (|B_a|/L)^{1−1/w} (CT)^{1/w}.
D2.4 Hence ‖W‖_p ≪ T^{1/p} iff Σ_a e^{θa − κa(1−1/w)/p} < ∞, i.e. pθ < κ (w large).
D2.5 COROLLARY (candidate (a)): ∫_T^{2T} |R_T|^q ≪ T whenever qθ < κ (Hölder: ‖(1+W)(1+On)^{k−1}‖_q ≤ ‖1+W‖_{qs}‖(1+On)^{k−1}‖_{qs'}).
  Conrey κ = 8/7: every θ < 4/7 has some q > 2. Jutila κ = 1: only θ < 1/2 (no gain over MV). So the gain is
  exactly Conrey's mollifier length 4/7 > 1/2, imported through the density theorem.
D2.6 Sanity: q = 4 gives θ < 2/7 (weaker than MV's 1/2); under RH (no off-line zeros) all q, θ < 1 (= DPC Lemma 5).
  Worst case = isolated off-line zeros, each giving |R_T| ~ X^{d}: threshold qθ < κ is the heuristic truth of the method.

## D3 (14:05) Step 7 with q > 2 only
  δ₁T ≤ ∫_{E_T} G_T |R_T|² |P(V)| ≤ ‖R_T‖_q² ‖P(V)‖_s meas(E_T)^{1 − 2/q − 1/s} (s large, 5c gives ‖P(V)‖_s ≪ T^{1/s}),
  so meas(E_T) ≥ c₂T with c₂ independent of ε. Needs only q > 2. OK for θ < κ/2.

## D4 (14:10) Step 6e replaced: hairpins via Lemma Z, threshold 2θ < κ

Off' := Σ_{ρ: d>0} X^{d} κ(L|γ − t|) (so W ≤ On + Off'). On the hairpin of ρ ∈ S_a (x = x_ρ, 0 < y < d_ρ):
  X^{y}(1 + W_{−y}) ≤ X^{y}(1 + On) + Off', hence
  |R_T R*|(x − iy) ≤ C [X^{y} + Off'(x)] (1 + Off'(x)) (1 + On(x))^{2k}.
∫_0^{d_ρ} X^{y} dy ≤ X^{d_ρ}/log X ≤ e^{θ(a+1)}/(θL); ∫_0^{d_ρ} dy ≤ (a+1)/L. With Z(x − iy) := (1 + |V_L| + 2Ñ_ρ)^D
(Hölder in y, exponents r, r'): hairpin(ρ) ≲ ((a+2)^C/L) [e^{θa} + Off'(x_ρ)] (1 + Off'(x_ρ)) Y_ρ,
  Y_ρ := (1 + On(x_ρ))^{2k} Z̄_ρ, Z̄_ρ := (avg_{0<y<Y_a} Z(x_ρ − iy)^{r'})^{1/r'}.
2 e^{θa} Off' ≤ e^{2θa} + Off'²  ⇒  Σ_{S_a} hairpin ≲ ((a+2)^C/L)[(1 + e^{2θa}) Σ_{S_a} Y_ρ + Σ_{S_a} Off'(x_ρ)² Y_ρ].
 (i) Σ_{S_a} Y_ρ ≤ |S_a|^{1−1/w} (Σ_{all ρ} Y_ρ^w)^{1/w} ≲ (a+2)^C σ_a^{1−1/w} TL   [discrete moments over all zeros:
     On via Lz DoubleSumZeros + Fujii; Z̄ via 6b (Poisson, Σ_ρ P_y(x_ρ − ·) ∈ L^{r'} with norm ≪ L T^{1/r'}) and 6d].
 (ii) Σ_{S_a} Off'² Y ≤ |S_a|^{1/v} (Σ_{all ρ} Off'(x_ρ)^{2u})^{1/u} (Σ_all Y^w)^{1/w}, 1/v + 1/u + 1/w = 1.
     Σ_ρ Off'(x_ρ)^{2u} ≪ TL needs 2uθ < κ (discrete analogue of D2.3: #{ρ: a B_{a'} zero within (|m|+1)/L of x_ρ} [radius fixed after referee B: m ∈ ℤ]
     ≤ Σ_{ρ'∈B_{a'}} Ñ_{ρ'} ≲ (TL e^{−κa'})^{1−1/w}(TL)^{1/w}). u > 1 possible iff θ < κ/2.
Total: Σ_a 𝓗_{T,a} ≲ T Σ_a (a+2)^C [e^{2θ(a+1)} σ_a^{1−1/w} + σ_a^{1/v}],  σ_a = min(2ε + o(1), C e^{−κa}).
 First sum → 0 as ε → 0 iff 2θ < κ(1 − 1/w) (split at A = c log(1/ε), c < (1−1/w)/(2θ)); second → 0 always.
CONCLUSION (candidate): Steps 5a, 6e, 7 hold for every θ < κ/2; with Conrey κ = 8/7 − δ: θ < 4/7, so φ > 1/(1 + 4/7) = 7/11.
 => λ₀ = 7/11 ≈ 0.6364 unconditionally?  MUST CHECK: (1) ρ_η smooth in DPC; (2) all other θ < 1/2 uses in PROOF.md;
 (3) Z̄ discrete moments; (4) Step 6c Poisson-tail costs T^{3θ/2}L = o(T)?  For θ > 2/3 fails, here θ < 4/7 fine.

## D5 (14:15) Checks
- DPC l.673–679: ρ_η(u) = ψ(u)ω(u/η), ψ, ω smooth; supp ρ_η ⊂ (η, b), b > θ; prime factors > T^η, k = ⌊θ/η⌋ (DPC's
  convention; PROOF.md writes η/2). f smooth, compact support on R, f = 0 on [θ, ∞) (DPC l.619). So φ(u) = ρ_η(u/L)/u ∈ C_c^∞((ηL, bL)): IBP valid.
- DPC's own RH route (l.665–670) bounds |e_r| ≪ (1 + |S_1|)^r: loses the length constraint (X^{rd} at an off-line zero);
  harmless under RH, fatal without. Lemma Z keeps f(Σu_i/L) inside the r-fold integral: X^{max d}, not X^{rd}. This is the new point.
- Shifted points: block (non-distinct) and prime-power terms for R*(t − iy): coefficient q^{−|B|/2 + |B|y} with q^{|B|} ≤ X/(rest):
  use the shifted-f bound (X/∏q^{|B|})^{y + max d⁺} for the free rest ⇒ total X^{y + max d⁺} Σ_q q^{−|B|/2}ρ^{|B|} = X^{y}·O(1)·W-structure. OK.
- Other θ < 1/2 uses in PROOF.md: only 5a, 6e (A₁, A₂), 6f, Step 7 Hölder. 6c Poisson tails T^{3θ/2}L = o(T) need θ < 2/3 only
  (and are not needed for R with Lemma Z). Step 0 construction needs 0 < φθ < 1 < φ(1+θ) < 2 and θ < 1.
- Hairpin y-integral for a = 0 (d_ρ tiny): use ∫_0^{d_ρ} Z^{r'} ≤ ∫_0^{Y_a} Z^{r'}, so Z̄_ρ := (Y_a^{−1}∫_0^{Y_a} Z^{r'} dy)^{1/r'} for every zero.
- Discrete moments over all zeros: Σ_ρ Z̄_ρ^m ≤ ∫ v^m · Σ_ρ Y_a^{−1} k_{Y_a}(· − x_ρ) ≤ ‖v^m‖_r ‖·‖_{r'} ≪ T^{1/r} · L T^{1/r'} = TL (6c with σ = 1).

## D6 (14:58) Review 1 and numerics
- Review 1 (logs/review1.txt): Lemma Z GAP — block prime q vs packet q^h bookkeeping at s = ½ − y + it; repaired (packets 𝔞 = ∏q^h,
  𝔞^{−1/2+y}(X/𝔞)^{y+D} = X^{y+D}𝔞^{−1/2−D}). Lemma M, 6e′, 7′: CORRECT (minor |m|+2, Fujii at comparable heights; applied).
- N1 (n1_mechanism.py, logs/n1.log): planted zero of depth a/L, j copies: I_j ≈ e^{θa}·(polynomial loss), lossy |I_1|^j ≈ e^{jθa}
  (a = 40: I_2 = 2.5e6, I_1² = 3.5e12, e^{θa} = 3.6e9). Confirms X^{max d}, not X^{jd}.
- N2 (n2_tails.py, logs/n2.log): T = 1e4, 1e6, θ = 0.45/0.55/0.65: E|R|²/S = 1.000, E|R|⁴/(E|R|²)² = 1.4–1.8, max |R|²/S ≤ 10.8,
  L²-tail mass beyond W = 8: ≤ 0.0008. No concentration at θ > ½ (numerically RH holds there, so this is the on-line regime of Lemma Z).
- N3 (n3_identity.py, logs/n3.log): the j = 1, 2 multivariate explicit formula with 482 real zeta zeros near t = 1000:
  |LHS − RHS| ≤ 4.5e−12 (j = 1), ≤ 9.8e−13 (j = 2). Signs, normalization and Fubini structure confirmed. Peak 22 MB.

## D7 (15:35) Close
- Review 3: 6e′ (all-zero moments) CORRECT; Theorem U Steps 0, 3, 7: no other θ < ½ use. Review 4: repaired Lemma Z CORRECT.
- Result: THEOREM U′, λ > max{2/(2 + κ), ½} = 7/11 (Conrey κ = 8/7). Wall for this method: θ < κ/2 (X^{2d} per isolated off-line zero in
  its hairpin, X^{qd} in the q-th moment; heuristically sharp, not proved — referee C). Lever: any uniform slope κ > 8/7 near σ = ½ with log T loss.


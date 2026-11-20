# hmean — derivations (checkable line by line)

Notation: ρ = ½ + iγ (RH where stated), τ = γ/2π, ℓ = log τ, L = log(T/2π), y = T^θ, B(s) = Σ_{k≤y} b(k)k^{-s}, b(k) = μ(k)P(log(y/k)/log y),
P(0) = 0, P(1) = 1, λ = log y·∫₀¹P = θL∫P, μ := λ/L = θ∫P, u = ½ + μ. S₁ = Σ_{0<γ≤T} ζ′(ρ)B(ρ), S₂ = Σ|ζ′(ρ)B(ρ)|², S_B = Σ ζ′(ρ)²B(ρ)².
H: S_B = (TL³/2π)(u² − 1/12) + o(TL³) for 0 < θ < ½.

## A. Exact identities (no RH)
A1. ζ(s) = χ(s)ζ(1−s) and χ(s)χ(1−s) = 1 give ζ′(s) = χ′(s)ζ(1−s) − χ(s)ζ′(1−s) and (χ′/χ)(s) = (χ′/χ)(1−s).
A2. At every nontrivial zero ζ(1−ρ) = 0, so ζ′(ρ) = −χ(ρ)ζ′(1−ρ); hence S_B = −Σ_γ χ(ρ)ζ′(ρ)ζ′(1−ρ)B(ρ)².
    Under RH ζ′(1−ρ) = conj ζ′(ρ), so S_B = −Σ_γ χ(ρ)|ζ′(ρ)|²B(ρ)²: a second moment twisted by χ(ρ)B(ρ)².
A3. Under RH, Re(ζ′(ρ)B(ρ)) = ½ζ′(ρ)ψ(ρ), ψ(s) = B(s) − χ(1−s)B(1−s), and Σ(Re Y)² = ½S₂ + ½Re S_B (explore55).

## B. Main term of H by the ratios recipe (re-derived; sympy: recipe_check.py, logs/recipe_check.log)
B1. Σ_γ G(ρ) = (1/2πi)∮(ζ′/ζ)G; right edge ≈ 0 (no n = 1 coefficient); on the left edge (ζ′/ζ)(s) = (χ′/χ)(s) − (ζ′/ζ)(1−s).
B2. G = ζ(s+α)ζ(s+β)B(s)², B(s) = (1/2πi)∫ y^z K(z)ζ(s+z)^{-1}dz, K(z) = Σ_j p_j j!/((log y)^j z^{j+1}) (P = Σ p_j x^j).
    Mollifier integrals used: (1/2πi)∫y^zK(z)dz = P(1) = 1; (1/2πi)∫y^zK(z)z^{-1}dz = λ; (1/2πi)∫y^zK(z)(z−α)z^{-1}dz = 1 − αλ.
B3. Recipe for ζ(s+α)ζ(s+β)ζ(1−s+δ)/(ζ(s+z₁)ζ(s+z₂)ζ(1−s+γ)), then ∂_δ at δ = γ → 0 (arithmetic factor = 1 at zero shifts; its
    derivatives are O(1/L) relative). Terms: χ′/χ: ℓ (no α, β); identity: (ζ′/ζ)(1+α) + (ζ′/ζ)(1+β) + 2λ; swaps α↔δ, β↔δ:
    D(α,β) = [β²g(α) − α²g(β)]/(αβ(β−α)), g(x) = e^{−xℓ}(1 − xλ)².
B4. Only D couples α and β: ∂_α∂_β D|₀ = −g₃ = ℓ³/6 + ℓ²λ + ℓλ² (sympy (i)). Integrating, S_B ~ (TL³/2π)(1/6 + μ + μ²) = (TL³/2π)(u² − 1/12).
B5. Same recipe, checks: Σζ′(ρ)B(ρ) ~ (T/2π)(L²/2 + Lλ) = BHB's S₁ (sympy (iv)); θ → 0 (B ≡ 1): Σζ′(ρ)² ~ (T/2π)L³/6, which is the
    leading term of Hughes–Pearce-Crump's conjecture (arXiv 2509.07792, conj:secmom; also conj:n1 with 1/Γ(k+2) at k = 2), obtained there
    by random-matrix theory and the ratios conjecture: an independent source for the θ = 0 constant.

## C. Modulus localization (new; recipe + numerics)
C1. Ẽ_K(s) := Σ_{n≤K}(ℓ − log n)n^{-s}. The approximate functional equation of ℓζ(s) + ζ′(s) with lengths K, τ/K, at a zero
    (ζ(ρ) = 0, ℓ + (χ′/χ)(ρ) = O(1/γ)), gives ζ′(ρ) = Ẽ_K(ρ) + χ(ρ)G_{τ/K}(1−ρ) + error, G_M(w) := Σ_{n≤M}(log n)n^{-w}; the error is
    small in the Riemann–Siegel range K ≍ √τ, not uniformly (sharp-cutoff endpoint terms; reviewer 1). Checked: at K = √τ the identity
    Σζ′B²(ζ′ − Ẽ_K) = Σχζ′B²G_{τ/K}(1−ρ) holds to 4 digits at 10¹² (logs/kappa_hi2.log). K = τ (G = 0) is used only as a reading:
    numerically Σζ′(ρ)Ẽ_τ(ρ)B(ρ)² = S_B to 4 digits at height 7·10⁴ (logs/kappa_low.log). Reading: S_B ≈ Σ_γ ζ′(ρ)A(ρ), A = Ẽ_τB² of
    length τy² = T^{1+2θ}, a first moment at zeros with a twist longer than T, for every θ ≥ 0.
C2. F(K) := Σ_γ ζ′(ρ)Ẽ_K(ρ)B(ρ)² keeps the moduli n ≤ K. Perron: Ẽ_K(s) = (1/2πi)∫_(c)(K^z/z)(ℓ + ∂_β)ζ(s+β)|_{β=z}dz, so
    F(K) = (1/2πi)∫(K^z/z)[(ℓ + ∂_β)∂_α S_B(α,β)]_{α=0,β=z}dz, S_B(α,β) := Σ_γ ζ(ρ+α)ζ(ρ+β)B(ρ)².
C3. From B3, the α-linear part at α = 0 (the ±1/α cancel) is (T/2π)Φ(β), Φ(β) = 1/β² + g₁/β + g₂ − g(β)/β² = −Σ_{j≥3}g_jβ^{j−2}
    (entire; sympy (ii)), g₁ = −(ℓ + 2λ), g₂ = ℓ²/2 + 2λℓ + λ².
C4. (ℓ + ∂)Φ = [−2/β³ + (ℓ − g₁)/β² + ℓg₁/β + ℓg₂] + e^{−βℓ}[2h/β³ − h′/β²], h = (1 − βλ)². With K = τ^κ the first bracket is
    a residue at z = 0 (contour to the left); the second carries (K/τ)^z and vanishes for κ < 1 (contour to the right).
C5. Result (sympy (iii)): F(τ^κ) ~ (T/2π)ℓ³ f(κ), f(κ) = 1/6 + μ + μ² + μ(1−κ)² + (1−κ)³/3 (0 < κ < 1; κ = 0 means Ẽ = ℓ, the
    limit κ → 0⁺; κ → 1⁻ gives H; residues are taken with λ/ℓ, replaced by μ = λ/L only at the end), f(1) = u² − 1/12,
    f(0) = ½ + 2μ + μ² (= ℓ·Σζ′(ρ)B(ρ)² by the same recipe; θ = 0: Fujii/Shanks ½).
C6. Reading: S_B − F(τ^κ) ~ −(T/2π)ℓ³[μ(1−κ)² + (1−κ)³/3]; at κ = ½ the bracket is 1/24 + μ/4, i.e. 25 % (θ = 0), 23 % (θ = ¼),
    22 % (θ = ½) of u² − 1/12 (P = x). The density −f′(κ) = (1−κ)² + 2μ(1−κ) is spread over all κ < 1.

## D. Why these routes do not prove H (θ = 0 included). Heuristic except D1 and A2; not a non-provability claim (reviewer 2)
D1. Status. θ = 0 (B ≡ 1) of H is Σ_γ ζ′(ρ)² ~ (T/2π)L³/6, conjectured (RH) by Hughes–Pearce-Crump (2509.07792, conj:secmom) and proved
    nowhere; every proved discrete moment of ζ^{(n)}(ρ) is a first moment (Fujii, Kaptan–Karabulut–Yıldırım, Hughes–Pearce-Crump 2106.03005,
    Durkan–Hughes–Pearce-Crump 2601.18025) or Hermitian (Gonek: Σ|ζ′(ρ)|²; BHB: S₂). So H at any θ would contain new mathematics beyond θ = 0.
D2. Contour (sketch; the O(T^{-A}) needs localization and contour estimates not written here): with smooth w on [T, 2T] and W(s) =
    w((s−½)/i), Σ_γ w(γ)ζ′(ρ)²B(ρ)² = −(1/2π)Σ_m Λ(m)m^{-1/2}∫ζ′(½+it)²B(½+it)²m^{it}w(t)dt + small (right edge and χ′/χ part: the
    Dirichlet series of ζ′²B² has no n = 1 term).
D3. χ² structure. ζ′(s)² = χ(s)²[Σ_a(log a − ℓ)a^{s−1}]² and ∫χ(½+it)²X^{it}w dt has its stationary point at t₀ = 2π√X with value
    √(πt₀)e^{i(2t₀+π/4)}w(t₀): the left edge becomes −(e^{iπ/4}/√2)Σ_h c_h h^{-3/4}Σ_{a,b,m}Λ(m)(log a−ℓ)(log b−ℓ)(abm)^{-1/4}w(t₀)e(2√(abm/h)),
    c_h = Σ_{k₁k₂=h}b(k₁)b(k₂). The phases are irrational unless abm/h is a square; exact squares give only O(T^{1/2}L³).
D4. Formal self-reference (stationary-phase heuristic). For m ≥ T^{1+ε}, Σ_m Λ(m)e(2√(abm/h))g(m) = ∫g − Σ_{ρ′}ĝ(ρ′); ∫g = O(T^{-A}) (no stationary point) and the ρ′-sum
    is stationary exactly at zeros with |γ′| ≈ t₀ ∈ [T, 2T], where it reassembles Σ_{γ′} χ(ρ′)²[Σ_{ab≤…}(log a−ℓ)(log b−ℓ)(ab)^{ρ′−1}]B(ρ′)²,
    the χ²-dual part of the approximate functional equation of ζ′(ρ′)², i.e. the original sum. RH (|ĝ| summed over zeros) only gives
    O(T^{3/2}L) here: the needed cancellation is among the zeros themselves. (Heuristic, stationary phase:) for fixed prime m ≤ T^{1−ε}, Poisson in (a, b) resonates only at
    ν | m (ν ∈ {1, m}), and those terms depend on α alone or β alone in the shifted version (they die under ∂_α∂_β): the coupled (swap) term,
    which carries H, needs the smooth part of Λ at composite m or the self-referential range.
D5. χ¹ form. By C1, S_B ≈ Σ_γ ζ′(ρ)A(ρ), A = Ẽ_τB² of length T^{1+2θ}. On the left edge ζ′(s) brings one χ, Gonek's lemma gives rational
    phases e(ma/k) with moduli k ≤ T^{1+2θ}: exactly BHB's S₁ problem with y replaced by T^{1+2θ}. BHB's range is y < T^{1/2}: their
    Lemma 2 (large moduli, via Heath-Brown's generalised Vaughan identity) gives M_{ν,3} ≪ y^{1/3}T^{5/6+ε}, which at y = T^{1+2θ} is
    T^{7/6+2θ/3+ε} ≫ TL³, and their O(yT^{1/2+ε}) term becomes T^{3/2+2θ+ε} (formal substitution; BHB's hypotheses not checked for this
    twist). The large-moduli estimate is where this route stops, for every θ ≥ 0.
D6. Not an error-term problem only. By C5–C6, the moduli beyond BHB's reach carry a fixed part of the main term itself: BHB can treat
    Ẽ_KB² only for Ky² < T^{1/2}, i.e. κ < ½ − 2θ (nothing for θ ≥ ¼), and F(τ^{½−2θ}) − S_B ~ (TL³/2π)[μ(½+2θ)² + (½+2θ)³/3]
    ≥ (TL³/2π)/24 when ∫P ≥ 0 (recipe; confirmed numerically, E2).
    The missing part is, by the approximate functional equation of ℓζ + ζ′ with lengths K, τ/K at a zero,
    S_B − F(K) ≈ Σ_γ χ(ρ)ζ′(ρ)B(ρ)²G_{τ/K}(1−ρ) = −Σ_γ χ(ρ)²ζ′(1−ρ)G_{τ/K}(1−ρ)B(ρ)², G_M(w) = Σ_{n≤M}(log n)n^{-w}: a χ²-object again.
D7. What would suffice: a first-moment theorem Σ_γ ζ′(ρ)A(ρ) = recipe + o(TL³) for divisor-bounded A of length T^{1+2θ} (θ > ¼ beats
    19/27), equivalently an asymptotic for the χ²-sums of D3/D6. Both cross the length-T barrier (DHPC's Σχ(ρ)X^ρ changes main term at
    X = T/2π); nothing in the cited literature crosses it for a holomorphic (non-Hermitian) moment. Accounting (θ = 0): the χ¹ piece of
    S₀ ≈ ΣẼ_{√τ}(ρ)Ẽ_τ(ρ) + Σχ(ρ)G_{√τ}(1−ρ)Ẽ_τ(ρ) is Σ_{n≤τ, n′≤√τ}(weights ≍ L²(nn′)^{-1/2})·Σ_γ χ(ρ)(n′/n)^ρ/n′; DHPC's error O(T^{1/2}log²T)
    per X, summed without cancellation over the ≍ T^{3/2} pairs, is O(T^{5/4}L⁴) ≫ TL³: cancellation over X is needed, not a better pointwise X-bound.
D8. Dead end, heuristic (avoid H by a lower bound): Σ(Im Y)² ≥ (Σ_γ Im Y_γ c_γ)²/Σc_γ² for c_γ = Re/Im D(ρ), D short (length T^κ,
    κ + θ < ½) needs only first moments and Landau–Gonek sums; but the projection of Y onto such D is its low-frequency part, coefficients
    ≈ (P′(1)/θ − P(1))Λ(N), N ≤ T^κ, which vanish at BHB's optimal P (P′(1) = 1 − θ, θ → ½); captured energy O(κ²)(1−2θ)²/θ² against the
    loss 3(1+θ)^{-4}κ ≈ 0.59κ from shortening θ to ½ − κ: no gain over 19/27 (heuristic; any c containing ζ′ or χ brings back S_B).

## E. Numerics (kappa.py; P = x, ∫P = ½; means over zeros of ζ′(ρ)·X·B(ρ)²/ℓ²; ζ′ by Riemann–Siegel)
E1. Zeros 50001–10⁵ (ℓ_rms = 9.11, 59 s, memcap peak 70 MB, logs/kappa_low.log): F(1) = exact S_B to 4 digits; F(κ) − F(1) vs f(κ) − f(1):
    θ = 0: κ = .5 .0436/.0417, .6 .0223/.0213, .7 .0096/.0090, .8 .0028/.0027; θ = ½ (μ = .257): .5 .1103/.1060, .6 .0644/.0625, .7 .0318/.0322.
E1b. Zeros 1981053–2001052 (Odlyzko zeros6 tail, ℓ_rms = 12.10, 320 s, peak 161 MB, logs/kappa_z6.log): F(κ) − F(1) vs f(κ) − f(1):
    θ = ½: κ = .3 .2390/.2369, .4 .1676/.1620, .5 .1069/.1042, .6 .0612/.0613, .7 .0289/.0315, .8 .0115/.0127; θ = ¼: .5 .0764/.0729, .6 .0437/.0413,
    .7 .0200/.0203; θ = 0: .5 .0439/.0417, .6 .0273/.0213; exact S_B/ℓ²: .5062 (θ = ½), .3350 (θ = ¼), .2047 (θ = 0).
E2. Zero indices 10¹²+1…+10⁴ (ℓ = 24.475, N_RS = 206393, 151 s, peak 242 MB, logs/kappa_hi2.log; κ ≤ ½ direct, κ ≥ ½ via the χ²-tail of C1/D6;
    log column 4 of tails_chi2 is the predicted F − S_B, the measured tail is S_B − F):
    | θ | exact S_B (f(1)) | F − S_B at κ = .1 .2 .3 .4 .5 (pred) | S_B − F at κ = .5 .6 .7 .8 .9 (pred, negated) |
    | ½ | .4847 (.4792) | .443 .341 .244 .162 .100 (.446 .331 .237 .162 .104) | −.1004 −.0673 −.0277 −.0115 −.0022 (−.1042 −.0613 −.0315 −.0127 −.0028) |
    | ¼ | .3161 (.3073) | .346 .259 .180 .111 .067 (.344 .251 .176 .117 .073) | −.0671 −.0413 −.0202 −.0087 −.0011 (−.0729 −.0413 −.0203 −.0077 −.0016) |
    | 0 | .1894 (.1667) | .237 .173 .128 .105 .055 (.243 .171 .114 .072 .042) | −.0553 −.0178 −.0053 −.0005 −.0018 (noisy: |Im| up to .04) |
    At κ = ½ the direct F and the χ²-tail agree to 4 digits (identity of C1 at K = √τ). κ = 0 deviates by O(1/L) (θ = ½: 1.009 vs 1.0625).
E3. Three-piece split at K = M = √τ (ζ′ = Ẽ + χG at a zero; S_B = P₀ + 2X₁ + X₂, P₀ = ΣẼ²B² (no χ), X₁ = ΣχẼG(1−ρ)B², X₂ = Σχ²G(1−ρ)²B²;
    per zero, units ℓ²; logs/kappa_low2.log, logs/kappa_hi3.log; the sum reproduces S_B to 4–5 digits):
    | ℓ | θ | P₀ | X₁ | X₂ (pure χ²) | S_B |
    | 9.11 | ½ | .911 | −.271 | .160 | .530 |  | 24.48 | ½ | .999 | −.414 | .313 | .485 |
    | 9.11 | ¼ | .589 | −.165 | .085 | .345 |  | 24.48 | ¼ | .564 | −.181 | .114 | .316 |
    | 9.11 | 0 | .288 | −.062 | .018 | .183 |  | 24.48 | 0 | .336 | −.091 | .036 | .189 (θ = 0 at 10¹²: |Im| ≤ .15, noisy) |
    | 12.10 | ½ | .898 | −.285 | .178 | .506 |  | 12.10 | ¼ | .577 | −.166 | .089 | .335 |  (12.10 | 0 | .307 | −.059 | .015 | .205)
    The pure χ² piece is not small and grows with ℓ (θ = ½: .160 → .178 → .313 at ℓ = 9.1, 12.1, 24.5), as do |X₁| and P₀, while their sum converges: a proof through this split
    needs more than an upper bound for the χ² piece (finite samples; no asymptotic law claimed). (At K·M = τ the recipe's double Perron integral loses its exponential factor and diverges
    formally, consistent with pieces larger than S_B; no prediction is claimed for them.)

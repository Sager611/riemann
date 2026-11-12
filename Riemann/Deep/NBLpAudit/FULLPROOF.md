# NBLp audit: line-by-line proofs of PROOF.md §3–§7

Referee: Orca agent `nblpaudit`, 12 Nov 2026. Object: `Riemann/Deep/NBLp/PROOF.md` (cited below as P).
Notation is P's (§0): `ρ_k(x) = {1/(kx)}` on (0,∞); `f_k = ρ_k − k^{-1}ρ_1` on (0,1) (= brief's `ρ_{1/k}`);
`χ = 1_(0,1)`, `χ₊ = 1_(0,1]`; `I_n = (1/(n+1), 1/n]`; `V`, `W_p`, `S_p`, `d_{N,p}`; `1/p = 1/2 + 1/r`;
`L²_σ = L²((0,∞), x^{2σ−1}dx)`; `γ(s) = ½s(s−1)π^{−s/2}Γ(s/2)`, `Γr_ε = γ(·+ε)/γ`; `P_N(s) = Σ_{k≤N} μ(k)(1−k/N)k^{−s}`.
Each step: the claim, the proof written out, and a verdict tag: ✓ (correct as written), ✓* (correct, P's text
needed the added detail given here), ✗ (false). No ✗ was found; the ✓* items are listed in §G.

## §A Proposition C and Proposition N (P §3)

**C(a).** `x ∈ I_n`, `y = 1/x ∈ [n, n+1)`. For integer `k ≥ 1`: `⌊y/k⌋ = ⌊⌊y⌋/k⌋` (standard: `⌊⌊y⌋/k⌋ ≤ y/k`, and if
`m ≤ y/k` for an integer `m` then `mk ≤ ⌊y⌋`), `= ⌊n/k⌋`. So `{y/k} − k^{−1}{y} = y/k − ⌊n/k⌋ − (y − n)/k = {n/k}`. ✓

**C(b).** Each `f_k` is constant on each `I_n`, so `V ⊂ S_p`. `F ↦ (F|I_n)_n` is an isometry `S_p → ℓ^p(w)`,
`w_n = 1/(n(n+1))`, onto. `S_p` closed: an L^p-convergent sequence has an a.e.-convergent subsequence, and a.e.
limits of a.e.-constants on `I_n` are a.e.-constant on `I_n`. Non-density: `g = 1_(1/2,3/4) − 1_(3/4,1)` ∈ L^q,
`∫_{I_1} g = 1/4 − 1/4 = 0`, `g = 0` off `I_1`, so `∫ gF = 0` for `F ∈ S_p`, but `∫ g 1_(1/2,3/4) = 1/4`.
Translates `g_n` supported in `I_n` (same construction) are linearly independent functionals vanishing on `S_p`,
so `codim S_p = ∞`. Hence no subset of `{f_k}` spans a dense subspace of `L^p(0,1)`, any `p ∈ [1,∞)`. ✓

**C(c).** `‖D_mF‖_p^p = ∫_0^{1/m}|F(mx)|^p dx = m^{−1}‖F‖_p^p`. For `x < 1/m` (so `mx < 1`):
`f_{km}(x) − k^{−1}f_m(x) = {1/(kmx)} − (km)^{−1}{1/x} − k^{−1}{1/(mx)} + (km)^{−1}{1/x} = f_k(mx)`.
For `x > 1/m`: `1/(mx) < 1`, `1/(kmx) < 1`, so `f_{km}(x) − k^{−1}f_m(x) = 1/(kmx) − 1/(kmx) + (km)^{−1}({1/x} − {1/x}) = 0 = D_m f_k(x)`.
(`k = 1`: `f_1 = 0` and `f_m − f_m = 0`.) `D_m` is bounded linear and maps `V` into `V`, hence `W_p` into `W_p`. ✓

**C(d).** `W_p ⊂ S_p` (closed). If `χ ∈ W_p`: `D_mχ = 1_(0,1/m) ∈ W_p` for all `m ≥ 1`, so
`1_{I_m} = D_mχ − D_{m+1}χ` (a.e.) `∈ W_p`; finite combinations of `1_{I_m}` are dense in `S_p ≅ ℓ^p(w)` (`p < ∞`);
so `W_p = S_p`. Conversely `W_p = S_p ∋ χ`. `d_{N,p}` is non-increasing and `χ ∈ W_p ⟺ inf_N d_{N,p} = 0`. ✓

**Proposition N** (`1 < p < ∞`; a zero `w_0`, `β_0 = Re w_0 > 1/p` ⟹ `χ ∉ W_p`). `q = p/(p−1)`;
`∫_0^1 |x^{w_0−1}|^q dx = ∫_0^1 x^{q(β_0−1)}dx < ∞ ⟺ q(β_0 − 1) > −1 ⟺ β_0 > 1 − 1/q = 1/p`. On `(1,∞)`:
`ρ_k(x) − k^{−1}ρ_1(x) = 1/(kx) − 1/(kx) = 0`. So, with Lemma T (needs `0 < β_0 < 1`, true for a non-trivial zero),
`∫_0^1 f_k x^{w_0−1}dx = ∫_0^∞ (ρ_k − k^{−1}ρ_1)x^{w_0−1}dx = −(ζ(w_0)/w_0)(k^{−w_0} − k^{−1}) = 0`. The bounded
functional `F ↦ ∫_0^1 F x^{w_0−1}` kills `V`, hence `W_p`, but gives `1/w_0 ≠ 0` at `χ`. ✓
(Trivial zeros are irrelevant: `β_0 > 1/p > 0`.) For `p > 2`, `w_0 = 1/2 + 14.1347…i` gives `χ ∉ W_p`. ✓

## §B Mellin–Plancherel toolkit (P §4)

Convention check: `ĥ(t) = ∫ h(u)e^{−iut}du`, `φ(u) = F(e^{−u})`. With `x = e^{−u}`, `dx = −x du`:
`∫_ℝ φ(u)e^{−σu}e^{−iut}du = ∫_0^∞ F(x)x^{σ+it}x^{−1}dx = ∫_0^∞ F(x)x^{s−1}dx`, and
`∫_ℝ|φ(u)e^{−σu}|²du = ∫_0^∞|F|²x^{2σ−1}dx`. So `𝓜_σ` is `√(2π)`×unitary `L²_σ → L²(ℝ)`. ✓
**Only Plancherel in L² is used anywhere in P.** No Hausdorff–Young, no L^p Mellin bound for `p ≠ 2`: every
L^p(0,1) norm is reached from an `L²_a` norm by Lemma M with `a < 1/p` (uses: §6 `a = σ_0`, A5 `a = σ'`,
A6 `a = σ''`; each `< 1/p`, checked at the place of use). ✓

**Lemma T** (`0 < σ < 1`). `ρ_k(x) = ρ_1(kx)`. `0 ≤ ρ_1 < 1` on `(0,1]`, `ρ_1(x) = 1/x` on `(1,∞)`. Hence
`∫_0^∞|ρ_1|x^{σ−1} ≤ 1/σ + 1/(1−σ)` and `‖ρ_1‖²_{L²_σ} ≤ ∫_0^1x^{2σ−1} + ∫_1^∞x^{2σ−3} = 1/(2σ) + 1/(2−2σ)`;
`x ↦ kx` gives `‖ρ_k‖_{L²_σ} = k^{−σ}‖ρ_1‖_{L²_σ}`. Mellin: `u = 1/(kx)`, `x = 1/(ku)`, `dx = −du/(ku²)`,
`x^{s−1} = (ku)^{1−s}`: `∫_0^∞{1/(kx)}x^{s−1}dx = k^{−s}∫_0^∞{u}u^{−s−1}du = −k^{−s}ζ(s)/s` by [C1]
(Titchmarsh (2.1.5), valid `0 < σ < 1`). `x^{σ−1}ρ_k ∈ L¹` makes the transform the pointwise integral. (b)
`∫_0^1x^{s−1}dx = 1/s`, `σ > 0`. ✓

**Lemma M** (`1 ≤ p < 2`, `a < 1/p`). Hölder with exponents `2` and `r` (`1/p = 1/2 + 1/r`) applied to
`|F| = |x^{a−1/2}F|·x^{1/2−a}`: `‖F‖_p ≤ ‖x^{a−1/2}F‖_{L²(0,1)}‖x^{1/2−a}‖_{L^r(0,1)}`, and
`‖x^{1/2−a}‖_r^r = ∫_0^1x^{r(1/2−a)}dx = 1/(1 + r(1/2 − a))`, finite iff `r(a − 1/2) < 1` iff `a < 1/2 + 1/r = 1/p`.
For `a = 1/p − ε/2`: `1 + r(1/2 − 1/p + ε/2) = 1 − 1 + rε/2 = rε/2`, constant `(2/(rε))^{1/r}`. Also
`‖x^{a−1/2}F‖_{L²(0,1)} ≤ ‖F‖_{L²_a}` (the `L²_a` norm integrates over `(0,∞)`). ✓
Sharpness remark: at `a = 1/p` the weight `x^{1/2−1/p}` lies only in weak-`L^r`, so `L²_{1/p}` controls only
the Lorentz norm `L^{p,2} ⊋ L^p` (`p < 2`); this is why P never works on the line `1/p` itself. ✓

**Lemma PW** (shift of the Mellin line). Hypotheses: `Ψ` analytic on a neighbourhood of `{a ≤ σ ≤ b}`,
`M_Ψ = sup_{a≤σ≤b}∫|Ψ(σ+it)|²dt < ∞`, `F ∈ L²_b`, `𝓜_bF = Ψ(b+i·)` a.e.
Proof. `ψ_σ := ` inverse Fourier transform of `Ψ(σ+i·) ∈ L²` (`σ = a, b`); `ψ_σ^{T}(u) := (2π)^{−1}∫_{−T}^{T}Ψ(σ+it)e^{iut}dt`.
(1) Tonelli: `∫_ℝ dT ∫_a^b(|Ψ(σ+iT)|² + |Ψ(σ−iT)|²)dσ = ∫_a^b dσ ∫_ℝ(...)dT ≤ 2(b−a)M_Ψ`, so the inner integral
`h(T)` is integrable on `[1,∞)`, hence `liminf_{T→∞} h(T) = 0`: pick `T_n ↑ ∞` with `h(T_n) → 0`.
(2) Cauchy on `[a,b]×[−T_n,T_n]` for `Ψ(s)e^{us}` (`u` fixed): `i∫_{−T_n}^{T_n}Ψ(b+it)e^{u(b+it)}dt − i∫_{−T_n}^{T_n}Ψ(a+it)e^{u(a+it)}dt`
equals the difference of the two horizontal integrals `∫_a^bΨ(σ±iT_n)e^{u(σ±iT_n)}dσ`, each bounded by
`e^{|u|max(|a|,|b|)}(b−a)^{1/2}h(T_n)^{1/2}` (Cauchy–Schwarz). So `e^{ub}ψ_b^{T_n}(u) − e^{ua}ψ_a^{T_n}(u) → 0` for every `u`.
(3) Plancherel [C7]: `ψ_σ^{T_n} → ψ_σ` in `L²(du)`; pass to a common subsequence converging a.e. for both `σ`.
Hence `e^{ua}ψ_a = e^{ub}ψ_b` a.e. Since `𝓜_bF = Ψ(b+i·)`, `ψ_b(u) = F(e^{−u})e^{−bu}` a.e. (inverse transform of
`𝓜_bF`), so `ψ_a(u) = F(e^{−u})e^{−au}` a.e.; `ψ_a ∈ L²` gives `F ∈ L²_a`, and `𝓜_aF = ψ̂_a = Ψ(a+i·)`. ✓
(Same conclusion for every `a' ∈ [a,b]`, applying the lemma on `[a',b]`.)

## §C Lemma Γ and Lemma H (P §5)

**Lemma Γ(a).** `log γ(w) = log(½) + log w + log(w−1) − (w/2)log π + log Γ(w/2)`, so
`γ'/γ(w) = 1/w + 1/(w−1) − ½log π + ½ψ(w/2)` (`ψ = Γ'/Γ`). ✓
Digamma. Binet [C3] for `Re z > 0`: `log Γ(z) = (z−½)log z − z + ½log 2π + ∫_0^∞ q(t)e^{−tz}t^{−1}dt`,
`q(t) = ½ − 1/t + 1/(e^t−1)`. Differentiation under `∫` is legitimate on `Re z ≥ c > 0` (integrand derivative
`−q(t)e^{−tz}`, dominated by `(t/12)e^{−ct}`), giving `ψ(z) = log z + (z−½)/z − 1 − ∫_0^∞q(t)e^{−tz}dt
= log z − 1/(2z) − ∫_0^∞q(t)e^{−tz}dt`. ✓
`q ∈ [0, t/12]`: `½ + 1/(e^t−1) = ½coth(t/2)`, and DLMF 4.36.3 at `z = t/2`:
`½coth(t/2) = 1/t + Σ_{k≥1} t/(t²/4 + k²π²) = 1/t + Σ_k 2t/(t² + 4k²π²)`; so `q(t) = Σ_k 2t/(t²+4k²π²) ≥ 0` and
`q(t) ≤ Σ_k 2t/(4k²π²) = (t/(2π²))(π²/6) = t/12`. ✓
For `Re z ≥ 1/4`: `|∫q e^{−tz}| ≤ ∫(t/12)e^{−t/4}dt = 16/12 = 4/3`; `|1/(2z)| ≤ 1/(2Re z) ≤ 2`; so
`Re ψ(z) ≤ log|z| + 2 + 4/3 ≤ log|z| + 4`. With `z = w/2` (`Re w ≥ 1/2`):
`½Re ψ(w/2) ≤ ½log|w| − ½log 2 + 2 ≤ ½log(1+|w|) + 2`. `Re(1/w) = Re w/|w|² ≤ 1/Re w ≤ 2`; `−½log π < 0`;
`Re 1/(w−1) = (Re w − 1)/|w−1|²` is `≤ 0` if `Re w < 1`, and `≤ 1/|w−1| ≤ 1/|Im w| ≤ 1` if `|Im w| ≥ 1`.
Total: `≤ 4 + ½log(1+|w|)` (`Re w < 1`), `≤ 5 + ½log(1+|w|)` (`|Im w| ≥ 1`). ✓
Numerical cross-check (§H, `numA.py`): on test grids `max(Re γ'/γ(w) − ½log(1+|w|)) = −0.92` (`Re w < 1`) and `−0.42`
(`|Im w| ≥ 1`, `Re w ≤ 3`), far below the claimed 4 and 5.

**Lemma Γ(b).** On `{Re w ≥ 1/2}` the only zero of `γ` is `w = 1` (`Γ` has no zeros; `w = 0` excluded;
`Γ(w/2)` has no poles there). Under the segment condition the segment `[s, s+ε]` avoids `w = 1`, so
`u ↦ log|γ(s+u)|` is `C¹` on `[0,ε]` with derivative `Re(γ'/γ)(s+u)`, and
`log|Γr_ε(s)| = ∫_0^ε Re(γ'/γ)(s+u)du ≤ ε(5 + ½log(1 + |s| + ε)) ≤ ε(5 + ½log(2+|s|))` (`ε ≤ 1`). ✓
If `Re(s+ε) < 1`: every `w = s+u` has `Re w < 1`, `|w| ≤ |Re w| + |t| ≤ 1 + |t|`, so
`2log|Γr_ε(s)| ≤ 2ε(4 + ½log(2+|t|))`, i.e. `|Γr_ε(s)|² ≤ e^{8ε}(2+|t|)^ε`. ✓

**Lemma H** (all zeros `β ≤ θ*`; `σ ≥ θ* − ε/2`, `σ > 0`, `s ≠ 1`; `ε > 0` real). If `ξ(s) = 0`: trivial.
Else `ξ(s+ε) ≠ 0` (`Re(s+ε) ≥ θ* + ε/2 > θ*`). [C2]: the Hadamard product converges absolutely, all factors at
`w ∉ {ρ}` are non-zero, so `log|ξ(w)| = Re A + Re(Bw) + Σ_ρ(log|1 − w/ρ| + Re(w/ρ))`, absolutely convergent.
Subtract (`w = s+ε` minus `w = s`): `Re(Bε) + Σ_ρ[log|(ρ−s−ε)/(ρ−s)| + ε Re(1/ρ)]`. Since `Σ_ρ Re(1/ρ)`
converges absolutely (`0 < Re(1/ρ) = β/|ρ|² ≤ |ρ|^{−2}`), split: `= ε(B + Σ_ρRe(1/ρ)) + Σ_ρ log|(ρ−s−ε)/(ρ−s)|`,
and `B = −Σ_ρ Re(1/ρ)` [C2] kills the first part. (`B` real; Davenport Ch. 12, `B = −½γ_E − 1 + ½log 4π`.)
Each term: `ρ − s` and `ρ − s − ε` have the same imaginary part `γ − t`, so
`|ρ−s−ε|² − |ρ−s|² = (β−σ−ε)² − (β−σ)² = ε(ε − 2(β−σ)) ≥ 0 ⟺ β ≤ σ + ε/2`, true as `β ≤ θ* ≤ σ + ε/2`. ✓
So `|ξ(s)| ≤ |ξ(s+ε)|`. `γ(s) ≠ 0` (`σ > 0`, `s ≠ 1`), `ζ = ξ/γ` off `s = 1`:
`|ζ(s)| ≤ |ζ(s+ε)||Γr_ε(s)|` (when `s + ε = 1` read the right side as `|ξ(1)|/|γ(s)|`). ✓
Equality case: if every zero had `β = θ*` and `σ = θ* − ε/2`, every factor is 1 (e.g. under RH, `θ* = 1/2`,
`|ξ(1/2 − ε/2 + it)| = |ξ(1/2 + ε/2 + it)|` by `ξ(s) = ξ(1−s)` and conjugation): the lemma is sharp, and the
line `σ' = θ* − ε/2` used in A4 is the leftmost line where it holds for an arbitrary zero set in `β ≤ θ*`. ✓

## §D Lemma A (Lindelöf-type bound in a zero-free half-plane)

Statement (P): `1/2 ≤ θ < 1`, `ζ ≠ 0` on `Re s > θ`, `0 < δ ≤ (1−θ)/2`, `η' > 0`; then for `σ ≥ θ + 2δ`:
`|ζ(s)|^{−1} ≤ C_A(2+|t|)^{η'}` (all `t`), `|ζ(s)| ≤ C_A(2+|t|)^{η'}` (`|t| ≥ 1`), `C_A = C_A(θ,δ,η')`.
Reduction `t ≥ 0`: `ζ(s̄) = conj ζ(s)`. ✓
*Compact part* `K = {θ+2δ ≤ σ ≤ 3, 0 ≤ t ≤ 6}`: `1/ζ` is continuous on `K` (`ζ` is
meromorphic with only the simple pole at 1, so `1/ζ` extends analytically across `s = 1` with value 0, and `ζ`
has no zeros in `K`: none off the real axis by hypothesis, none on `(θ,1)` (`ζ(σ) < 0` for `0 < σ < 1`), none on
`[1,3]`); `ζ` is continuous on `K ∩ {t ≥ 1}`. For `σ ≥ 3`: `|ζ|^{±1} ≤ ζ(3)` [C8]. Constants depend on `θ, δ`. ✓
*Branch.* `Ω = {Re s > θ} \ (θ,1]` is star-shaped w.r.t. `2` (a segment `[2,z]` with `Im z ≠ 0` meets ℝ only at `2`;
for real `z > 1` it stays in `(1,∞)`), hence simply connected. `ζ` is analytic and zero-free on `Ω`, so `log ζ`
exists on `Ω`; normalize it by `log ζ = Σ_pΣ_m p^{−ms}/m` on `Re s > 1 ⊂ Ω`. ✓
*(1) Borel–Carathéodory*, any `t ≥ 4`: `z_0 = 2+it`, `R = 2−θ−δ/2`, `r = 2−θ−δ`. Disc `|z−z_0| ≤ R`:
`Re z ≥ 2 − R = θ + δ/2 ≥ 1/2`, `Im z ≥ t − R ≥ t − 3/2 ≥ 5/2`, `Im z ≤ t + 3/2`. [C6] gives on it
`log|ζ(z)| ≤ log C_6 + ½log(t + 3/2) ≤ log C_6 + ½log(t+2) =: U` (`C_6 ≥ 1` so `U > 0`); `Re log ζ = log|ζ|`.
`|log ζ(z_0)| ≤ Σ_pΣ_m p^{−2m}/m = log ζ(2) = 0.4977 < 1`. [C4]: for `|z−z_0| ≤ r`,
`|log ζ(z)| ≤ (2r/(R−r))U + ((R+r)/(R−r))|log ζ(z_0)|`, `2r/(R−r) = 4(2−θ−δ)/δ ≤ 6/δ`,
`(R+r)/(R−r) = 2(4 − 2θ − 3δ/2)/δ ≤ 6/δ`. Every `σ + it` with `θ+δ ≤ σ ≤ 2` has `|σ+it − z_0| = 2−σ ≤ r`. For
`σ ≥ 2`: `|log ζ| ≤ log ζ(2) < 1`. So (★) `|log ζ(z)| ≤ (6/δ)(log C_6 + ½log(Im z+2) + 1) ≤ (C_7/δ)log(Im z + 2)`
for `Re z ≥ θ+δ`, `Im z ≥ 4`, `C_7 := 6(log C_6 + 1)/log 6 + 3` (absolute). ✓
*(2) Three circles* for `t ≥ 6`, `σ ∈ [θ+2δ, 1+δ]`: centre `z_1 = 1+2δ+it`, `r_1 = δ`, `r_2 = 1+2δ−σ ∈ [δ, 1−θ]`,
`r_3 = 1+δ−θ`. `r_3 ≤ (3/2)(1−θ) ≤ 3/4`; disc `|z−z_1| ≤ r_3` has `Re z ≥ θ+δ`, `Im z ∈ [t−3/4, t+3/4] ⊂ [4, ∞)`, so it
lies in `Ω` and (★) applies on `|z−z_1| = r_3`: `M(r_3) ≤ (C_7/δ)log(t+4) =: X`; enlarge `C_7` so `X ≥ 1`.
On `|z−z_1| = r_1`: `Re z ≥ 1+δ`, `M(r_1) ≤ log ζ(1+δ) ≤ log(1 + 1/δ)` [C8]. `r_2 < r_3` strictly; if `r_2 = r_1`
the point is on the inner circle. [C5]: `M(r_2) ≤ M(r_1)^{1−a}M(r_3)^a`, `a = log(r_2/r_1)/log(r_3/r_1) ∈ [0,1)`;
`a ≤ a_0 := log((1−θ)/δ)/log((1−θ+δ)/δ) < 1` (`(1−θ)/δ ≥ 2`); `M(r_1)^{1−a} ≤ max(1, log(1+1/δ))`,
`M(r_3)^a ≤ X^a ≤ X^{a_0}` (`X ≥ 1`). The target `σ+it` lies on `|z−z_1| = r_2`, so
`|log ζ(σ+it)| ≤ C_9(log(t+4))^{a_0}`, `C_9 = max(1, log(1+1/δ))(C_7/δ)^{a_0}`. ✓
*(3)* Young: `C_9y^{a_0} ≤ η'y + C_10(C_9, a_0, η')` (`a_0 < 1`). For `σ ≥ 1+δ`: `|log ζ| ≤ log ζ(1+δ)`. So for
`t ≥ 6`, `σ ≥ θ+2δ`: `|ζ|^{±1} ≤ e^{|log ζ|} ≤ e^{C_10 + log ζ(1+δ)}(t+4)^{η'} ≤ 2^{η'}e^{…}(2+t)^{η'}`. ✓
Dependence: `C_A` depends on `θ, δ, η'` and the absolute `C_6` only. ✓ (`C_6` is not made explicit in P; only its
existence is used. §H checks (★)'s shape and the final bound numerically for `t ≤ 10⁴`.)

## §E Lemma P (Riesz–Perron) and Theorem B (P §5–§6)

**Kernel.** For `y > 0`, `c > 0`: `(2πi)^{−1}∫_{(c)}y^w dw/(w(w+1)) = (1 − 1/y)_+`. `y > 1`: close left (`|y^w| ≤ y^c`
on the left arcs, integrand `O(|w|^{−2})`), residues `1` (at 0) and `y^{−1}/(−1)` (at −1). `y ≤ 1`: close right,
`|y^w| ≤ 1`, no poles. (`y = 1` included: absolutely convergent, closing right.) ✓
**Lemma P.** `c > max(0, 1−σ)`, `y = N/k`: `Σ_k μ(k)k^{−s}(1 − k/N)_+ = P_N(s)`; interchange of `Σ_k` and `∫` is
justified by `Σ_k k^{−σ−c}N^c∫|dw|/|w(w+1)| < ∞` (`σ + c > 1`); `Σ_kμ(k)k^{−s−w} = 1/ζ(s+w)` there. So
`P_N(s) = (2πi)^{−1}∫_{(c)}N^w dw/(ζ(s+w)w(w+1))`. In `−η ≤ Re w ≤ c`: `Re(s+w) ≥ σ−η ≥ θ+2δ > θ`, so `1/ζ(s+w)` is
analytic (value 0 at `s+w = 1`); `w = −1` lies outside (`η < 1`); only `w = 0` is a pole, residue `1/ζ(s)`
(which is `0` if `s = 1`: then `1/(ζ(1+w)w)` is regular at 0). Horizontal segments `Im w = ±T`: `|N^w| ≤ N^c`,
`|1/ζ(s+w)| ≤ C_A(2+|t|+T)^{η'}` (Lemma A where `Re(s+w) ≤ 3`, else `≤ ζ(3)`), `|w(w+1)| ≥ T²`, length `c+η`:
`→ 0`. So `P_N(s) − 1/ζ(s) = (2πi)^{−1}∫_{(−η)}…`; with `w = −η+iv`, `|N^w| = N^{−η}`,
`2 + |t+v| ≤ (2+|t|)(1+|v|)`, `|w||w+1| ≥ max(η,|v|)max(1−η,|v|)`:
`|P_N(s) − 1/ζ(s)| ≤ N^{−η}(2+|t|)^{η'}·(C_A/2π)∫_ℝ(1+|v|)^{η'}dv/(max(η,|v|)max(1−η,|v|))`, finite for `η' < 1`. ✓
Case `σ − η > 1`: `P_N(s) − 1/ζ(s) = −N^{−1}Σ_{k≤N}μ(k)k^{1−s} − Σ_{k>N}μ(k)k^{−s}`;
`|·| ≤ N^{−1}Σ_{k≤N}k^{−η} + Σ_{k>N}k^{−1−η} ≤ N^{−1}·N^{1−η}/(1−η) + N^{−η}/η` (integral comparison, both
monotone). No zero hypothesis used. ✓
Dependence: `C_P = C_P(θ, δ, η, η')` through `C_A(θ,δ,η')` and the `v`-integral. ✓

**Theorem B.** Given `η < 1/p − θ`: `σ_0 = (θ+η+1/p)/2`; `σ_0 − θ − η = (1/p − θ − η)/2 > 0` and `1/p − σ_0` is the
same positive number, so `σ_0 ∈ (θ+η, 1/p)`. `δ = min{(σ_0−η−θ)/2, (1−θ−η)/2, (1−θ)/2} > 0` (`η < 1/p − θ < 1 − θ`),
giving `σ_0 − η ≥ θ+2δ`, `1 − η ≥ θ+2δ`, `σ_0 ≥ θ+2δ`, `δ ≤ (1−θ)/2`. ✓
Algebra on `(0,1)`: for `c_k = μ(k)(1−k/N)`, `Σ_{k=1}^N c_kρ_k = Σ_{k=1}^N c_kf_k + (Σ_k c_k/k)ρ_1` and `Σ_kc_k/k = P_N(1)`,
`f_1 = 0`; so `F_N := χ + Σ_{k=2}^N c_kf_k = Φ_N − P_N(1)ρ_1` on `(0,1)`, `Φ_N := χ₊ + Σ_{k=1}^N c_kρ_k`. ✓
`‖F_N‖_p ≤ ‖Φ_N‖_{L^p(0,1)} + |P_N(1)|‖ρ_1‖_{L^p(0,1)}`, `0 ≤ ρ_1 < 1`; Lemma M (`a = σ_0 < 1/p`):
`‖Φ_N‖_{L^p(0,1)} ≤ (1 − r(σ_0 − 1/2))^{−1/r}‖Φ_N‖_{L²_{σ_0}}`. ✓
`Φ_N ∈ L²_{σ_0}` (finite sum; Lemma T, `0 < σ_0 < 1`), `𝓜_{σ_0}Φ_N = (1 − ζ(s)P_N(s))/s = ζ(s)(1/ζ(s) − P_N(s))/s`
(`ζ(s) ≠ 0`, `σ_0 > θ`). Lemma A (`η' = 1/8`, `σ_0 ≥ θ+2δ`): `|ζ| ≤ C_A(2+|t|)^{1/8}` (`|t| ≥ 1`), `≤ C_11` (`|t| ≤ 1`,
continuity on a compact set avoiding `s = 1`). Lemma P (`σ_0 − η ≥ θ+2δ`): `|1/ζ − P_N| ≤ C_PN^{−η}(2+|t|)^{1/8}`.
`‖Φ_N‖²_{L²_{σ_0}} = (2π)^{−1}∫|𝓜Φ_N|² ≤ (2π)^{−1}max(C_A,C_11)²C_P²N^{−2η}∫(2+|t|)^{1/2}(σ_0²+t²)^{−1}dt`
(integrand `~|t|^{−3/2}`). Lemma P at `s = 1` (`1 − η ≥ θ+2δ`, `1/ζ(1) = 0`): `|P_N(1)| ≤ 2^{1/8}C_PN^{−η}`. So
`‖F_N‖_p ≤ C(p,θ,η)N^{−η}` for all `N ≥ 1`. ✓ (This also reproves the open range of Theorem A: `p < 1/Θ`.)

## §F Theorem A, (i) ⟹ (ii) (P §7)

Hypothesis (i): no zero with `β > 1/p`; zeros ON `Re s = 1/p` are allowed (`Θ = 1/p` is the endpoint case).
**A0.** `ε_0 = min{1/2, (2/p−1)/2, 1−1/p}`, `0 < ε ≤ ε_0`, `σ' = 1/p − ε/2`. Then: `σ' ≥ 1/p − (2/p−1)/4 = 1/(2p) + 1/4 ≥ 1/2`
(`p ≤ 2`); `σ'+ε = 1/p + ε/2 ≤ 1/p + (1−1/p)/2 < 1`; `2σ' − 1 = (2/p − 1) − ε ≥ (2/p−1)/2 = (2−p)/(2p) > 0`;
`1 − ε − σ' = (1 − 1/p) − ε/2 > 0`; `ε ≤ 1/2`. ✓
**A1.** `φ_ε(j) = Σ_{d|j}μ(d)d^{−ε}`, `E_ε(y) = Σ_{j≤y}φ_ε(j) − y/ζ(1+ε)`. Fix `x > 0`, `y = 1/x`. For `k > y`, `{y/k} = y/k`,
so `Σ_kμ(k)k^{−ε}{y/k} = Σ_{k≤y}μ(k)k^{−ε}{y/k} + yΣ_{k>y}μ(k)k^{−1−ε}` (absolutely convergent)
`= yΣ_kμ(k)k^{−1−ε} − Σ_{k≤y}μ(k)k^{−ε}⌊y/k⌋ = y/ζ(1+ε) − Σ_{j≤y}φ_ε(j) = −E_ε(y)` (`Σ_{k≤y}μ(k)k^{−ε}#{m: km ≤ y}`
regrouped by `j = km`). So `Φ_ε := χ₊ − E_ε(1/·) = χ₊ + Σ_kμ(k)k^{−ε}ρ_k` pointwise. Checks: `x > 1`: both sides
`= y/ζ(1+ε)`; `x = 1`: both `= 1/ζ(1+ε)`. On `(0,1)`: `Φ_ε − 1 − {y}/ζ(1+ε) = −Σ_{j≤⌊y⌋}φ_ε(j) + ⌊y⌋/ζ(1+ε) = −E_ε(⌊y⌋) =: H_ε`,
constant on each `I_n`. Truncation: `c_k = μ(k)k^{−ε}(1−k/N)`, `Σ_kc_k/k = P_N(1+ε)`, so with `Φ_{ε,N} = χ₊ + Σ_{k≤N}c_kρ_k`,
`H_{ε,N} := Σ_{k=2}^Nc_kf_k = Φ_{ε,N} − 1 − P_N(1+ε){1/x}` on `(0,1)`; `H_{ε,N} ∈ V`. ✓
**A2.** `b ∈ (1−ε, 1)`. `Σ_k‖μ(k)k^{−ε}ρ_k‖_{L²_b} ≤ ‖ρ_1‖_{L²_b}Σ_kk^{−ε−b} < ∞` (Lemma T, `ε + b > 1`): the series
converges in `L²_b`; partial sums also converge pointwise (A1), so the `L²_b` sum is the pointwise sum a.e.;
`Φ_ε ∈ L²_b`. `𝓜_b` is continuous, so `𝓜_bΦ_ε = 1/s − (ζ(s)/s)Σ_kμ(k)k^{−ε−s}` in `L²(dt)`, the series converging
absolutely pointwise (`Re(s+ε) > 1`) to `1/ζ(s+ε)`: `𝓜_bΦ_ε = Ψ_ε(b+i·)`, `Ψ_ε = (1 − R_ε)/s`, `R_ε = ζ(s)/ζ(s+ε)`. ✓
**A3.** Analyticity: on `{1/p − ε < σ < 1}`, `ζ(s)` is analytic (`s ≠ 1`), `1/ζ(s+ε)` is analytic
(`Re(s+ε) > 1/p`: no zeros by (i); at `s+ε = 1` it vanishes), `1/s` analytic (`σ > 0`). This open strip contains
`σ' ≤ σ ≤ b`. Bounds for `σ' ≤ σ ≤ b`: Lemma H with `θ* = 1/p` (`σ ≥ θ* − ε/2`, `σ > 0`, `s ≠ 1`) gives
`|R_ε(s)| ≤ |Γr_ε(s)|` (if `ζ(s+ε) = ∞`, i.e. `s = 1−ε`, then `R_ε = 0`). `|t| ≥ 1`: Lemma Γ(b) (`σ ≥ 1/2`; segment
points have `|Im| ≥ 1`): `|Γr_ε| ≤ e^{5ε}(2+|s|)^{ε/2} ≤ e^{5/2}(3+|t|)^{1/4}` (`|s| ≤ 1+|t|`). `|t| ≤ 1`: `Γr_ε` is
continuous on the compact rectangle (`γ ≠ 0` there since `s ≠ 1`; `γ(s+ε)` finite), `≤ C_13(ε)`. So
`M_Ψ ≤ ∫_ℝ(1 + max(C_13, e^{5/2})(3+|t|)^{1/4})²(1/4+t²)^{−1}dt < ∞`. Lemma PW on `[a,b]`, any `a ∈ [σ',b]`:
`Φ_ε ∈ L²_a`, `𝓜_aΦ_ε = Ψ_ε(a+i·)`. ✓
**A4** (`B(ε) := ‖Φ_ε‖²_{L²_{σ'}} ≤ K_pε`). By A3 and Plancherel, `2πB(ε) = ∫_ℝ|1−R(σ'+it)|²|s|^{−2}dt`.
Pointwise `|1−R|² = 1 − 2Re R + |R|² ≤ 1 − 2Re R + |Γr_ε|²` (A3). All three pieces are integrable against `|s|^{−2}`:
`|R| ≤ |Γr_ε| ≤ max(C_13,e^{5/2})(3+|t|)^{1/4}`. So `2πB ≤ I_a − 2Re J + I_c` with
(a) `I_a = ∫dt/(σ'²+t²) = π/σ'`. ✓
(b) `J = ∫R(σ'+it)|s|^{−2}dt`. On `Re s = σ'`: `s̄ = 2σ' − s`, `|s|² = s(2σ'−s)`, `ds = i dt`, so `J = −i∫_{(σ')}f(s)ds`,
`f = R(s)/(s(2σ'−s))`. Singularities of `f` in `σ ≥ σ'`: `R` is meromorphic there with its only pole at `s = 1`
(simple, from `ζ(s)`; `ζ(s+ε) ≠ 0` for `Re(s+ε) ≥ 1/p + ε/2`; at `s = 1−ε` R has a zero, not a pole);
`1/(s(2σ'−s))` has the poles `0` (left of the line) and `2σ'` (`2σ' > 1`, so `2σ' ≠ 1`, simple). ✓
Rectangle `Q = [σ', X]×[−T,T]`, `X > 2σ'+1`, `T ≥ 1`, positively oriented. Residue theorem:
`∫_{Re s=X,↑} − ∫_{Re s=σ',↑} + (bottom, →) − (top, →) = 2πi(Res_1 f + Res_{2σ'} f)`.
Horizontal: on `Im s = ±T`, `σ' ≤ σ ≤ X`, `|s| ≥ T`, `|2σ'−s| ≥ T`, and Lemma H (valid for all `σ ≥ σ'`) + Lemma
Γ(b) (`|Im| = T ≥ 1`) give `|f| ≤ e^{5ε}(2+X+T)^{ε/2}T^{−2}`; times length `≤ X`: `→ 0` as `T → ∞`, `X` fixed. ✓
Vertical integrals converge absolutely (`|f| ≪ (3+|t|)^{1/4}|s|^{−2}` on `σ'`; `|f| ≤ ζ(X)²/((X−2σ')² + t²)` on `X`,
since `|ζ(s)| ≤ ζ(X)`, `|1/ζ(s+ε)| ≤ ζ(X+ε) ≤ ζ(X)`, `|s||2σ'−s| ≥ (X−2σ')² + t²` as `X > X−2σ' > 0`). So with
`T → ∞`: `∫_{(σ')}f = ∫_{(X)}f − 2πi(Res_1 + Res_{2σ'})`, and `|∫_{(X)}f| ≤ ζ(X)²π/(X−2σ') → 0` (`X → ∞`). ✓
Residues: `ζ(s) = 1/(s−1) + O(1)`, so `Res_1 R = 1/ζ(1+ε)` and `Res_1 f = 1/(ζ(1+ε)·1·(2σ'−1))`;
`f = −R(s)/(s(s−2σ'))`, `Res_{2σ'}f = −R(2σ')/(2σ')`. Hence
`J = −i·(−2πi)(Res_1 + Res_{2σ'}) = −2π(Res_1 + Res_{2σ'}) = 2π[R(2σ')/(2σ') − 1/((2σ'−1)ζ(1+ε))]`. ✓
(Real, as it must be: `R(s̄) = conj R(s)` makes `J = ∫Re R/|s|²`. Numerically confirmed, §H.)
(c) A0 gives `Re(s+ε) = σ'+ε < 1` on the line, so Lemma Γ(b) second part: `|Γr_ε|² ≤ e^{8ε}(2+|t|)^ε = e^x`,
`x = 8ε + εlog(2+|t|) ≥ 0`; `e^x − 1 ≤ xe^x`. So `I_c ≤ π/σ' + ε∫(8+log(2+|t|))e^{8ε}(2+|t|)^ε(σ'²+t²)^{−1}dt
≤ π/σ' + εe^4∫(8+log(2+|t|))(2+|t|)^{1/2}(1/4+t²)^{−1}dt = π/σ' + C_Γε` (`ε ≤ 1/2`, `σ' ≥ 1/2`). ✓
Sum: `2πB ≤ 2π/σ' − 2J + C_Γε = (2π/σ')(1 − R(2σ')) + 4π/((2σ'−1)ζ(1+ε)) + C_Γε`. ✓
`R(2σ') = ζ(2σ')/ζ(2σ'+ε) ≥ 1` (`ζ` decreasing on `(1,∞)`): first term `≤ 0`. `ζ(1+ε) ≥ 1/ε`, `2σ'−1 ≥ (2−p)/(2p)`:
`2πB ≤ 8πpε/(2−p) + C_Γε`, i.e. `B(ε) ≤ (4p/(2−p) + C_Γ/(2π))ε = K_pε`. ✓
`C_Γ`: recomputed (§H) `= 5993.05` (P: 5993.05…). ✓
*What A4 uses about the zeros:* only Lemma H on `σ ≥ σ'` (all `β ≤ 1/p`) and analyticity of `1/ζ(s+ε)` on
`Re(s+ε) ≥ 1/p + ε/2`. Nothing on `Re s = 1/p`; the cancellation of the main terms `2π/σ'` against
`2π R(2σ')/σ'` is exact (residue at `2σ'`), the remainder is the pole term `O(1/ζ(1+ε)) = O(ε)` and the gamma-
factor excess `O(ε)`. ✓
**A5.** On `(0,1)`: `χ + H_ε = Φ_ε − {1/x}/ζ(1+ε)` (A1). `‖{1/x}/ζ(1+ε)‖_p ≤ 1/ζ(1+ε) ≤ ε`. Lemma M (`a = σ' < 1/p`,
constant `(2/(rε))^{1/r}`) and A4: `‖Φ_ε‖_{L^p(0,1)} ≤ (2/(rε))^{1/r}B(ε)^{1/2} ≤ (2/r)^{1/r}K_p^{1/2}ε^{1/2−1/r}`, and
`1/2 − 1/r = 1/2 − (1/p − 1/2) = 1 − 1/p`. So `‖χ + H_ε‖_p ≤ A_pε^{1−1/p} + ε`. ✓
The exponent is positive exactly because `p > 1`; at `p = 1` the Hölder loss `ε^{−1/r}` equals the gain `ε^{1/2}`. ✓
**A6.** `σ'' = 1/p − ε/4 ∈ (σ', 1/p)`. `Φ_ε ∈ L²_{σ''}` with transform `Ψ_ε` (A3, `σ'' ∈ [σ', b]`); `Φ_{ε,N} ∈ L²_{σ''}` with
transform `1/s − (ζ(s)/s)Σ_{k≤N}μ(k)(1−k/N)k^{−s−ε} = (1 − ζ(s)P_N(s+ε))/s` (Lemma T). Difference:
`𝓜_{σ''}(Φ_ε − Φ_{ε,N}) = (ζ(s)/s)(P_N(s+ε) − 1/ζ(s+ε))`. ✓ On `Re s = σ''`:
* Lemma H with `ε/2` (`σ'' = 1/p − (ε/2)/2`): `|ζ(s)| ≤ |ζ(s+ε/2)||Γr_{ε/2}(s)|`. ✓
* Lemma A with `θ = 1/p ∈ [1/2,1)`, `δ = ε/8` (`≤ (1−1/p)/8 ≤ (1−θ)/2`), `η' = 1/8`: `Re(s+ε/2) = 1/p + ε/4 = θ+2δ`;
  `|ζ(s+ε/2)| ≤ C_A(2+|t|)^{1/8}` for `|t| ≥ 1`; `≤ C_14(ε)` for `|t| ≤ 1` (`1/p + ε/4 < 1`, compact, no pole). ✓
* Lemma Γ(b) with `ε/2` (`Re(s+ε/2) < 1`): `|Γr_{ε/2}| ≤ e^{5ε/2}(2+|s|)^{ε/4} ≤ e^{5/4}(3+|t|)^{1/8}` (P writes `e²`: weaker, fine). ✓
* Lemma P at the point `s+ε` (`Re = 1/p + 3ε/4`), shift `η = ε/2 ∈ (0,1)`, same `θ, δ`, `η' = 1/8`:
  `Re(s+ε) − η = 1/p + ε/4 = θ+2δ` ✓, so `|P_N(s+ε) − 1/ζ(s+ε)| ≤ C_P(p,ε)N^{−ε/2}(2+|t|)^{1/8}`. ✓
Product: `|𝓜(Φ_ε−Φ_{ε,N})|² ≤ C(p,ε)N^{−ε}(3+|t|)^{3/4}/(1/4+t²)`, integrable; `‖Φ_ε−Φ_{ε,N}‖²_{L²_{σ''}} ≤ C_15(p,ε)N^{−ε}`. ✓
Lemma P, case `σ−η > 1`, at `s = 1+ε`, `η = ε/2`: `|1/ζ(1+ε) − P_N(1+ε)| ≤ (1/(1−ε/2) + 2/ε)N^{−ε/2}`. ✓
`H_ε − H_{ε,N} = (Φ_ε − Φ_{ε,N}) − (1/ζ(1+ε) − P_N(1+ε)){1/x}` on `(0,1)` (A1); Lemma M (`a = σ''`, constant `(4/(rε))^{1/r}`):
`‖H_ε − H_{ε,N}‖_p ≤ C_17(p,ε)N^{−ε/2} → 0` (`N → ∞`, `ε` fixed). ✓
**A7.** Fix `ε`; take `N_0` with `C_17N_0^{−ε/2} ≤ ε`. For `N ≥ N_0`: `d_{N,p} ≤ ‖χ + H_{ε,N}‖_p ≤ ‖χ+H_ε‖_p + ‖H_ε−H_{ε,N}‖_p
≤ A_pε^{1−1/p} + 2ε` (`−H_{ε,N} ∈ span{f_2,…,f_N}`). `d_{N,p}` non-increasing, so `lim_N d_{N,p} ≤ A_pε^{1−1/p} + 2ε` for every
`ε ∈ (0,ε_0]`; hence `lim_N d_{N,p} = 0`. ✓
**Order of limits / the line `Re s = 1/p`.** Every constant (`C_13, C_14, C_15, C_17, C_P, C_A, N_0`) depends on
`(p, ε)` only; `N → ∞` is taken at fixed `ε > 0`, and every analytic input is used on lines `Re = σ', σ'', σ'+ε,
1/p+ε/4, 1/p+3ε/4` or to their right with `δ = ε/8 > 0` of room above `1/p` where zero-freeness is invoked
(Lemma A, Lemma P need `ζ ≠ 0` on `Re > θ = 1/p` only). The only use of information LEFT of `1/p` is Lemma H,
which needs `β ≤ 1/p` (not `< 1/p`). The `ε → 0` limit is taken only on the already-proved numerical bound
`A_pε^{1−1/p} + 2ε`. No uniformity in `ε` of any constant except `A_p, K_p` is needed, and those are explicit. ✓
Conclusion: Theorem A (i) ⟹ (ii) is proved, with the stated quantitative bound. Together with Prop. N and C(d):
Theorem A holds as stated. ✓

## §G Findings of the line-by-line check

No false step. Every ✓* is a gap of exposition, closed above in one or two lines; none changes a statement.
1. Lemma PW: P's "by Fubini there are `T_n → ∞` with `h_n → 0`" — needs `h ∈ L¹(1,∞) ⟹ liminf h = 0`; added. ✓*
2. Lemma A: `Ω` simply connected asserted; it is star-shaped w.r.t. `2`. ✓*
3. Lemma A (2): `M(r_3)^a ≤ X^{a_0}` needs `X = (C_7/δ)log(t+4) ≥ 1`; take `C_7 ≥ 1`. ✓*
4. Lemma A: `C_6` (Titchmarsh Thm 4.11 consequence) is never made explicit; only its existence is used, so
   `C_A`, `C_P` and all of Theorem A's `N_0(p,ε)` are ineffective. Theorem A's final bound `A_pε^{1−1/p} + 2ε` is
   explicit; the threshold `N_0` is not. (Matches P §9 "the bound is qualitative".) ✓
5. A1 writes "a function in `S`" for `S_p`. Typo. ✓*
6. A6 uses `e²` where `e^{5/4}` holds. Harmless. ✓
7. Lemma Γ(a)'s constants 4 and 5 are crude (true maxima `−0.92`, `−0.42` on grids); only `C_Γ` inherits them. ✓
8. Theorem B is stated for `N ≥ 2`; the proof gives `N ≥ 1`. ✓
9. Equivalence of formulations (not in P, needed to compare with Báez-Duarte's `L_p(0,∞)` wording): for
   `p > 1`, `χ₊ ∈ closure_{L^p(0,∞)} span{ρ_k : k ≥ 1}` ⟺ `χ ∈ W_p`. (⟹) on `(1,∞)`, `Σc_kρ_k − χ₊ = (Σc_k/k)/x`,
   so `|Σc_k/k|·‖1/x‖_{L^p(1,∞)} → 0`, and on `(0,1)` `‖χ − Σc_kf_k‖_p ≤ ‖χ − Σc_kρ_k‖_{L^p(0,1)} + |Σc_k/k|`.
   (⟸) `Σc_k(ρ_k − ρ_1/k)` is in the span, equals `Σc_kf_k` on `(0,1)` and `0` on `(1,∞)`. ✓ So Theorem A also
   answers Báez-Duarte's question in his own `L_p(0,∞)` formulation.
10. Load-bearing external inputs are only [C1]–[C8] (textbook) and, for Corollary 1, [C9]. [C9] and the
   `0.87338196661320752` half-plane are Lean-checked elsewhere (QRH port sorry-free, standard axioms; hpaudit
   11 Nov 2026, TREE.md N7); `0.87338196661320752 × 1.14497 = 0.9999962 < 1`, so `p ≤ 1.14497 ⟹ 1/p > 0.87338…`. ✓
11. Where the new content sits: Lemma H on the reflected line `σ' = 1/p − ε/2` (replacing Báez-Duarte's use of
   the functional equation on `1/2 − ε`, which only works for `p = 2`) plus the exact residue evaluation of `J`
   (A4), which turns `|R| ≤ |Γr_ε|` into `B(ε) = O(ε)`, beating Lemma M's loss `ε^{−1/r}` because `1/2 > 1/r`. The
   fixed-`ε` Möbius family `Σμ(k)k^{−ε}ρ_k` and "`n → ∞` first" are Báez-Duarte's 2003 scheme. ✓

## §I Cited results checked against sources (TeX via arXiv e-print; books from the standard editions)

* **Beurling 1955** (PNAS 41, 312–314; PMC528084, scan only, not fetched). Quoted verbatim in Dimitrov–Oliveira
  (1608.07887 TeX, Theorem A): with `C` = span of `{1/(θ_kx)}`, `θ_k ≥ 1`, `Σb_k/θ_k = 0`, "ζ does not vanish in
  `Re s > 1/p` iff `C^p = L^p(0,1)`", `p > 1`. This is the discovery's "Before" (`ρ_a`, `a = 1/θ ∈ (0,1]`). ✓
* **Báez-Duarte 2003** (math/0202141 TeX). Abstract/Theorem 1.1: RH ⟺ `χ ∈ closure_{L_2(0,∞)} span{ρ(1/(ax)) : a ∈ ℕ}`
  = [C10] after §G.9. §"Some comments": "One should extend it to the `L_p` case, that is, to the condition that ζ
  does not vanish in a half-plane `Re s > 1/p`." His Lemma (zratio) bounds `|ζ(½−ε+iτ)/ζ(½+ε+iτ)| ≤ C(1+|τ|)^ε` by
  the functional equation, unconditionally; P's Lemma H is the `1/p`-analogue and needs (i). ✓
* **Báez-Duarte 2005** (math/0505453 TeX), Remark after the General Strong Theorem: "One obviously expects a
  generalization to a sufficiency criterion for the non-vanishing of ζ in `Re s > 1/p` in terms of `L_p`
  approximations … we have neither worked out the necessity criterion in the `L_p`-case, nor … for the
  strengthened Nyman-Beurling Theorem." So "left open" is accurate (and already flagged in 2003). ✓
* **Dimitrov–Oliveira 2016** (1608.07887 TeX), Theorems 1–2: for Dirichlet `L`, `p ∈ (1,2]`: `L ≠ 0` on
  `Re s > 1/p` ⟺ `d_n(L,p) → 0`, `d_n² = inf ∫_0^∞|1_(0,1) − Σb_kκ(1/(kx))|²dx`, twisted
  `κ(x) = βx^α − Σ_{k≤x}χ(k)k^{1/2−1/p}`, `α = 3/2 − 1/p` (Thm 1); equivalently weighted `L²` on `Re s = 1/p` (Thm 2).
  Hilbert-space, twisted functions: not the `L^p(0,1)` statement with `{1/(kx)}`. Lemma M applied on a line
  `Re s = θ' ∈ (Θ, 1/p)` gives P's open range from it routinely (also needs `Σb_k/k → 0`, which the `L²(1,∞)` part
  supplies). The endpoint `p = 1/Θ` does not follow (§B, Lorentz remark). ✓
* **[C1]** Titchmarsh (2.1.5): `ζ(s) = −s∫_0^∞{x}x^{−s−1}dx`, `0 < σ < 1`; Báez-Duarte 2003 cites the same
  equation number for `−ζ(s)/s = ∫x^{s−1}ρ_1(x)dx`. ✓  **[C6]** Titchmarsh Thm 4.11 as quoted (uniform for
  `σ ≥ σ_0 > 0`, `|t| < 2πx/C`); the consequence `|ζ(σ+it)| ≪ |t|^{1/2}` (`σ ≥ 1/2`, `|t| ≥ 1`) follows as P says. ✓
* **[C2]** Davenport Ch. 12: product, order 1, `B = −Σ_ρRe(1/ρ)` (B real). ✓  **[C4], [C5]** Titchmarsh,
  *Theory of Functions* §5.5, §5.3, statements as used (B–C needs no sign condition on `U` in this form). ✓
* **[C3]** The `coth` expansion is DLMF 4.36.3 (checked on the DLMF page). The Binet formula
  `log Γ(z) = (z−½)log z − z + ½log 2π + ∫_0^∞(½ − 1/t + 1/(e^t−1))e^{−tz}t^{−1}dt` is on DLMF page §5.9(ii), but
  DLMF 5.9.1 is a different formula (`Γ(ν/μ)/(μz^{ν/μ}) = ∫exp(−zt^μ)t^{ν−1}dt`): wrong equation number, correct
  formula (Whittaker–Watson §12.31). Its consequence `Re ψ(z) ≤ log|z| + 4` is also checked numerically (§H). ✓*
* **[C9]** see §G.10. ✓

## §H Numerical falsification attempts (all single-threaded, under `memcap.sh 290`; peaks ≤ 25 MB)

Files: `numA.py`→`numA.txt`, `numB.py`→`numB_T4000.txt`, `numB_T1000_fine.txt`, `numB_T4000_fine.txt`,
`numC.py`→`numC.txt`, `lpnorm.c`→`lpnorm_M1e7.txt`. ζ by python-flint (arb, 64 bits), quadrature 8-point Gauss.
1. **`C_Γ`** = 5993.0534 (mpmath quad) = P's 5993.05. ✓
2. **Lemma Γ(a)** grid (`Re w ∈ [0.5,3)`, `Im w ≤ 200` and `10³…10⁵`): max excess `−0.92` / `−0.42` vs claimed 4 / 5. ✓
3. **Lemma H / Γ(b)** on the lines `σ'` for `(p,ε) = (1.1,0.1), (1.3,0.2), (8/7,0.05)`, `t ≤ 300`:
   `max|R|/|Γr_ε| = 0.998, 0.998, 0.999 ≤ 1`; `max|Γr_ε|²/(e^{8ε}(2+t)^ε) = 0.37, 0.14, 0.61 ≤ 1`. Sharpness with
   the true zeros (`θ* = 1/2`): `max|ξ(s)/ξ(s+ε)|` is `1.000000` on `σ = 1/2 − ε/2` (equality, functional equation),
   `≤ 1` at `+0.01`, and `1.22` (`ε = 0.2`), `2.33` (`ε = 0.05`) at `−0.01`, outside the lemma's range. ✓
4. **A4 residue identity** `J = 2π[R(2σ')/(2σ') − 1/((2σ'−1)ζ(1+ε))]`, `J` by quadrature on `[0,T]`:
   | (p, ε) | J numeric (T = 4000, fine panels) | formula | diff |
   |---|---|---|---|
   | (1.3, 0.2) | 3.43224623 | 3.43224623 | 2.6e−9 |
   | (1.1, 0.1) | 3.15289871 | 3.15289871 | 2.1e−9 |
   | (8/7, 0.05) | 3.43307983 | 3.43307983 | 9.2e−10 |
   A first run with width-0.5 panels gave `diff = −1.3e−6` at `(1.1, 0.1)`: the pole `s = 1` lies only `1 − σ' = 0.14`
   from the line and spoils Gauss on panels near `t = 0`; panels of width 0.01 on `[0,1]` remove it (table).
   The identity is confirmed to the quadrature error. ✓
5. **`B(ε)` and the A4 chain** (`2πB` numeric ≤ `π/σ' − 2J + I_c` ≤ final bound): `(1.3,0.2)`: 0.6647 ≤ 0.7399 ≤ 1201;
   `(1.1,0.1)`: 0.2726 ≤ 0.3139 ≤ 600; `(8/7,0.05)`: 0.0693 ≤ 0.0895 ≤ 300. `B/ε² = 2.6, 4.3, 4.4`: with the true
   zeros `B = O(ε²)`; the proof's `O(ε)` is the worst case over zero sets in `β ≤ 1/p`. ✓
6. **Lemma A shape** (`σ = 0.6, 0.75, 0.9`, `t ∈ [6, 6000]`): `max|log|ζ||/log(2+t)` over `[6,60], [60,600], [600,6000]`:
   `0.93, 0.46, 0.42` / `0.62, 0.30, 0.27` / `0.47, 0.25, 0.22`, decreasing with height; `max log|ζ|` grows sublinearly in
   `log t` (`σ = 0.75`: 1.07, 1.62, 1.89). Consistent with `|ζ|^{±1} ≪_{η'} t^{η'}`; `C_A` ineffective, so no sharper test. ✓
7. **Lemma P** (`numC`): `|P_N(s) − 1/ζ(s)|`, `N = 10²…10⁵`: at `t = 14.1347` local exponents `0.30, 0.30, 0.30`
   (`σ = 0.8`) and `0.44, 0.45, 0.45` (`σ = 0.95`) = `σ − 1/2`: the rate `N^{−η}`, `η < σ − θ`, holds and is sharp.
   Case `σ − η > 1` at `s = 1.1, 1.2`: errors `8e−6 … 2e−2` under the bounds `3.5 … 17`. ✓
8. **Theorem B decay** `‖F_N‖_p`, `F_N = 1 + Σ_{k≤N}μ(k)(1−k/N)ρ_{1/k}`, `M = 10⁷` breakpoints, exact integer sieve
   (direct-sum cross-check max deviation `2.5e−8`), tail by the mean of `|v|^p` on `(M/2, M]`:
   | N | 10 | 30 | 100 | 300 | 1000 | 3000 | 10000 |
   |---|---|---|---|---|---|---|---|
   | p = 1.1 | 0.345948 | 0.148676 | 0.059226 | 0.025883 | 0.015051 | 0.008404 | 0.004817 |
   | p = 1.3 | 0.349753 | 0.154557 | 0.068662 | 0.037184 | 0.025764 | 0.018951 | 0.012980 |
   | p = 2 | 0.364630 | 0.184879 | 0.132616 | 0.120989 | 0.120682 | 0.127894 | 0.122913 |
   Local exponents, `N ≥ 300`: `0.45, 0.53, 0.46` (`p = 1.1`; RH prediction `η < 0.409`), `0.31, 0.28, 0.31`
   (`p = 1.3`; `η < 0.269`); `p = 2` stalls (Báez-Duarte's `(log N)^{−1/2}`). Agreement with run1e9.txt: `N = 10, 100,
   1000` within its rigorous intervals (e.g. `0.0592261766 ∈ [0.0592261752, 0.0592262243]`). No decay slower than
   predicted; `p = 1.3` sits just above `1/p − 1/2`, as zero-driven terms `N^{ρ−1/p}` suggest. ✓

## §J What the endpoint adds (significance, for the DISCOVERIES gate)

Write `Θ = sup β`. Routine parts: Prop. N (`χ ∈ W_p ⟹ p ≤ 1/Θ`) and the open range (`p < 1/Θ ⟹ χ ∈ W_p`: Theorem B,
or Dimitrov–Oliveira Thm 2 on a line `θ' ∈ (Θ, 1/p)` + Lemma M). New part: `p = 1/Θ ⟹ χ ∈ W_p` (Theorem A at the
endpoint). Consequences:
* The endpoint matters only at the single exponent `p = 1/Θ`, which lies in `(1,2)` only if RH is false.
* Corollary 3 (RH ⟺ criterion for all `p ∈ (1,2)`) follows from the routine parts alone.
* Corollary 1 with the Lean-checked half-plane `Re s > 0.87338196661320752` (`1/B = 1.1449744…`): `p ≤ 8/7 = 1.142857`
  and `p ≤ 1.14497` are inside the open range, so they too follow from the routine parts. The endpoint adds,
  unconditionally, only `p = 1/0.87338196661320752 = 1.1449744…` itself (with 7/8 alone it would add `p = 8/7`).
* Corollary 2's exact equivalence `χ ∈ W_p ⟺ p ≤ 1/Θ` is where the endpoint is indispensable.
So the result is correct and answers Báez-Duarte's 2003/2005 question in exact form, but its new content is the
closing of a one-point gap; the discovery page's unconditional ranges should not be read as needing Theorem A.

## §K Side checks (outside the brief's §4–§7 scope)

* **NBLp2 (Dirichlet extension), changed steps only** (NBLp2/PROOF.md §5, §7): Lemma Γ_χ (`γ'/γ = ½log(q*/π) +
  ½ψ((w+a)/2)`, `Re z ≥ 1/4`, `|z| ≤ (1+|w|)/2`) ✓; Lemma E (`|1−z|/|1−zp^{−ε}| ≤ 1 + εlog p/(√p−1) ≤ p^{2.42ε}`,
  product `≤ q^{3ε}`) ✓; Lemma H_χ (only `Re B` enters since `ε` is real; zeros of `L(·,χ*)` in `Re > 0` are zeros of
  `L(·,χ)`) ✓; A4 for `δ = 0`: `J` complex, only `Re J` used, `2σ' ≥ u_0 = 1/2 + 1/p > 1`,
  `|1 − R(2σ')| ≤ ε|ζ'(u_0)|ζ(u_0)` (`|L(u) − L(u+ε)| ≤ Σn^{−u}(1−n^{−ε}) ≤ ε|ζ'(u)|`, `|1/L(u+ε)| ≤ ζ(u_0)`) ✓;
  (c) with `|R|² ≤ q^{7ε}e^{8ε}(2+|t|)^ε` ✓. No false step in the changed steps; the unchanged ones inherit §A–§F.
* **arXiv titles, 2026-11 listings to date** (math.NT 280, math.CA 109, math.FA 191 entries, 12 Nov): no title on
  Nyman–Beurling / Báez-Duarte / L^p closure criteria. (Novelty was not in this brief; title-level only.)

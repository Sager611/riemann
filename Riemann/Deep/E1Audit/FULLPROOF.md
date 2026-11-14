# E1Audit — line-by-line referee of Theorem E1 (santiaudit, 12 Nov 2026)
Referee: e1audit (Orca), 12 Nov 2026. Sources: SL = Santibáñez-Leal, Zenodo 23134787 v0.03; TeX = main.tex v0.03 from the CAOS source zip (Zenodo 23135354),
line numbers "l." refer to that TeX (sl.tex, 1167 lines); Wang = arXiv 2609.07918 v2 (zeros_short_intervals.tex). E1 = Riemann/Deep/SantiAudit/E1.md.
Equation map (PDF numbers = TeX labels): (10) eq:exactGaussian l. 363–369; (11) eq:representationerror l. 395–398; (12) eq:Estermann l. 433–438;
(13) eq:profileMellin l. 455–458; (14) eq:dualMellin l. 469–475; (15) eq:multipliernorm l. 485; profile A_{ε,j} l. 383–387; §3.6 Hankel l. 656–705;
§3.7 residues/compact windows l. 707–758; Thm 1.2 l. 158–176; Remark (hypotheses used) l. 336–346; Lemma 4.2 l. 823–888; Thm 1.3 proof l. 890–928; Thm 5.1 l. 932–963.

## 1. Inputs (10)–(14), re-derived independently
- (10)→profiles (l. 377–387). cos(2πxe^{2y/H}) = ½Σ_ε e(εx)·e^{ε2πix·2y/H}·exp[iφ_ε], φ_ε = ε2πx(e^{2y/H}−1−2y/H) REAL for real x, y. Taylor of exp(iφ) to order J−1:
  |remainder| ≤ |φ|^J/J!, |φ| ≤ 4πx(y²/H²)e^{2|y|/H} ⟹ kernel remainder ≪_J x^{−Re β}(x/H²)^J ✓ (l. 388). (2/√π)·½ = 1/√π and (iφ)^j/j! give exactly
  A_{ε,j}(x) = x^{−β}(ε2πix)^j/(j!√π)∫e^{−y²+(1−2β)y/H}e^{2i(T+ε2πx)y/H}(e^{2y/H}−1−2y/H)^j dy (l. 383–387) ✓; E1's Φ_j with w = (T−2πξ)/H is its ε = −1 case ✓.
- (11): per (h,k) the band x = nk/h ≍ T has ≍ Th/k integers n, outer factor 2π√(k/h)·b_hb_k/√(hk) ≍ 1/h, σ ≪ n^ε ⟹ Σ_{h,k≤y} T^{1+ε}/k ≪ T^{1+ν+ε}L ✓.
  Off-band: the y-phase 2Ty/H ± 2πxe^{2y/H} has derivative ≍ (T+x)/H on |y| ≤ √(A log T) unless x ≍ T/2π ⟹ O_A(x^{−Re β}(H/(T+x))^A) ✓.
- (12) checked three ways: (i) q = 1 reduces to ζ(s+a)ζ(s+b) = χ(s+a)χ(s+b)ζ(1−s−a)ζ(1−s−b) via 2 sin A sin B = cos(A−B) − cos(A+B) ✓;
  (ii) matches the standard Estermann FE (e.g. Bettin–Chandee–Radziwiłł Lemma: 2/k·(k/2π)^{2−2s−α−β}ΓΓ[cos(π(α−β)/2)D(1−s,h̄/k) − cos(π(s+(α+β)/2))D(1−s,−h̄/k)]) ✓;
  (iii) end-to-end numerical Voronoi identity with complex shifts (§4 below, voronoi_check.py). Residues (l. 444–447): the pole at s+a = 1 comes from
  u-sums with q | v (mean of e(Auv/q) over u is 1_{q|v}), giving q^{−1+a−b}ζ(1−a+b)F̃(1−a); symmetric term ✓.
- (14) from (12): u = 1−s−z/2 ⟹ 1−s−a = u−d/2, 1−s−b = u+d/2 (d = a−b), cos(π(s+z/2)) = −cos(πu) cancels the bracket's minus sign,
  q^{1−2s−z} = q^{2u−1}, (2π)^{2s+z−2} = (2π)^{−2u}, n^{s−1} = n^{−u−z/2}, F̃(s) = (q/p)^sÃ(s) (F(x) = A(xp/q)); ds = −du with reversed orientation keeps
  (2πi)^{−1}∫ upward. Product: (2/p)(q/p)^{−z/2}n^{−z/2}(pq/(4π²n))^u; times 2π/(gq) gives 4π/(gpq) ✓ (14) exactly. With a = α, b = −β: z = α−β, d = α+β ✓.
  Branch assignment: K-type (cos(πd/2)) ↔ e(Ān/q), J/Y-type ↔ e(−Ān/q), as in classical Voronoi (−2πY_0 with e(−ān/q), 4K_0 with e(ān/q)) ✓.

## 2. Steps 1–5 of E1 §2, re-derived
Step 1 (x-space form) ✓.
- DLMF 10.22.43: ∫t^{μ−1}J_d = 2^{μ−1}Γ((d+μ)/2)/Γ((d−μ)/2+1); reflection 1/Γ(1−(s−d)/2) = Γ((s−d)/2)sin(π(s−d)/2)/π gives E1's form ✓. 10.22.44 (Y_d) ✓.
  B̃(s) = (2^{s−1}/π)ΓΓ[sin(πd/2)sin(π(s−d)/2) − cos(πd/2)cos(π(s−d)/2)] = −(2^{s−1}/π)ΓΓcos(πs/2) ✓.
- s = 2u, x = 2√Y: Mellin transform of 𝒢(Y) := −πB(2√Y) is −π·2^{1−2u}·B̃(2u) = Γ(u−d/2)Γ(u+d/2)cos(πu) ✓; G(ξ) = 𝒢(4π²nξ/(pq)) has G̃(u) = (pq/(4π²n))^u ΓΓcos(πu) ✓.
- Absolute convergence of the FORWARD Mellin integral ∫|B(x)|x^{2σ−1}dx: at 0 needs 2σ > |Re d| (Y_d ≍ x^{−|Re d|}, Y_0 ≍ log x), at ∞ needs 2σ < 1/2 ✓.
  Caution (not an error): on Re u = σ < 1/4 the INVERSE integral (2πi)^{−1}∫ΓΓcos(πu)Y^{−u}du is not absolutely convergent (|ΓΓcos| ≍ |t|^{2σ−1}),
  so "𝒢(Y) = (2πi)^{−1}∫…" in E1 holds only conditionally; the proof never needs it: Parseval uses G ξ^{σ−1} ∈ L¹ and F̃(1−u) ∈ L¹(line) ✓.
- Contour: for each n the u-integrand of (14) is holomorphic on Re u > |Re d|/2 (Γ poles at u = ±d/2 − k; cos entire; Ã(1−u−z/2) entire, A ∈ C_c^∞(0,∞)),
  |ΓΓcos| ≪ (1+|t|)^{2σ−1} and Ã decays faster than any power uniformly on strips ⟹ horizontal segments vanish; the move from Re u > 1 is legitimate ✓.
  The dual series is then summed per n; absolute convergence of Σ_n afterwards (Step 5) justifies the per-n order ✓ (SL does the same by dyadic blocks, l. 476–478).
Step 2 (Hankel form) ✓. J = (H¹+H²)/2, Y = (H¹−H²)/(2i) ⟹ B = (−i/2)e^{iπd/2}H¹_d + (i/2)e^{−iπd/2}H²_d; with H^{1,2}_d = √(2/(πX))e^{±i(X−πd/2−π/4)}a_±
  (Watson 7.2(1)–(2)/DLMF 10.9 integral, Re d > −1/2, every X > 0) the d-phases cancel: B = √(2/(πX))(−i/2)[e^{i(X−π/4)}a_+ − e^{−i(X−π/4)}a_−] ✓ (= SL l. 677–681).
  |1 ± it/(2X)| ≥ 1, |arg| < π/2 ⟹ |(1 ± it/(2X))^{d−1/2}| ≤ e^{π|Im d|/2} (Re d < 1/2); with w = it/(2X): X∂_X(1+w)^c = −c(1+w)^c·w/(1+w), |w/(1+w)| ≤ 1,
  X∂_X[w/(1+w)] = −w/(1+w)+w²/(1+w)² ⟹ X^k∂_X^k a_± = O_k(1) uniformly on |d| ≤ 1/4 and ALL X > 0 (X ≥ 1 not needed) ✓.
Step 3 (profile bounds) ✓. ψ_j(y) = e^{−y²+(1−2β)y/H}[H²(e^{2y/H}−1−2y/H)]^j/(j!√π): since H²(e^{2y/H}−1−2y/H) = Σ_{m≥2}(2y)^mH^{2−m}/m! and H^{2−m} ≤ 1 (H ≥ 1),
  every derivative is ≤ C_k(1+y²)e^{2|y|} — ψ_j is Schwartz uniformly in H ≥ 1, |β| ≤ 1; Φ_j = H^{−2j}Ψ_j, Ψ_j = Fourier transform (frequency 2w) of ψ_j ✓.
  ∂_ξ on Ψ_j(w): −2π/H; on ξ^{j−β}, χ(ξ/T): ≪ 1/T ≤ 1/H ⟹ ∂^kA ≪_{j,k,K} (T/H²)^jH^{−k}E_K ✓ (|ξ^{−β}| = ξ^{−Re β} = O(1) for Re β ≪ 1/L).
  X = c√ξ: X' = X/(2ξ), X^{(k)} = c_kX/ξ^k, 1/X' = 2√ξ/c ⟹ (1/X')^{(k)} ≪ (T/X0)T^{−k}; Faà di Bruno with X^k∂_X^ka_± = O(1) gives ∂^k[X^{−1/2}a_±(X)] ≪ X0^{−1/2}T^{−k} ✓.
  X0 = 4π√(n/(2πpq))·√T = 2√(2π)√(nT/(pq)) ✓ (numerically 4π/√(2π) = 2√(2π)).
Step 4 (integration by parts) ✓. ∫Ge^{iX} = ∫G(−i/X')(e^{iX})' = i∫(G/X')'e^{iX}; induction: 𝒟^KG is a sum of ∂^mG·Π(1/X')^{(l)} with sizes ≪ (T/X0)^K H^{−m}T^{m−K} ≤ (T/(HX0))^K·|G|-scale ✓.
  T/(HX0) = (2√(2π))^{−1}√(N0/n), N0 = pqT/H² ✓; min{1,(c√(N0/n))^K} ≍_K (1+√(n/N0))^{−K} ✓; X0^{−1/2} ≍ (pq/(nT))^{1/4} ✓; constants uniform in n, p, q, shifts ✓.
Step 5 (summation) ✓. |σ_{−α,β}(n)| ≤ d(n)n^{|Re α|+|Re β|}; |n^{−z/2}| = n^{−Re z/2}; for n ≤ T^{O(1)} these are T^{O(1/L)} = O(1), beyond that (N0/n)^{K/2} absorbs them.
  N0 ≥ 1: Σ_{n≤N0}n^{−1/4+ε} + N0^{K/2}Σ_{n>N0}n^{−1/4+ε−K/2} ≪ N0^{3/4+ε} (K/2 > 3/4+ε). N0 < 1: ≪ N0^{K/2}ζ(K/2+1/4−ε) ≤ N0^{3/4} ✓.
  (1/(gpq))·H(pq/T)^{1/4}·(pq)^{3/4}(T/H²)^{3/4} = g^{−1}T^{1/2}H^{−1/2}: the pq powers cancel exactly ✓; Σ_{g≤y}g^{−1}(y/g)² ≤ ζ(3)y² ✓.
  Total ≪ T^{1/2+ε}H^{−1/2}y²; vs the main term ≍ H (SL l. 707–722: (H√π/g)[…], weights 1/(gpq), L^{O(1)}): relative T^{1/2+2ν−3θ/2+ε} ✓.
  Heuristic cross-check: Voronoi of a length-qH/p bump at qT/(2πp) modulo q has dual length q²X/Δ² = pqT/H² = N0 and term size (Δ/q)√q(nX)^{−1/4}: the same
  per-(g,p,q) total √(T/H)/g; at θ = 1 the range is the classical ν < 1/2 (Balasubramanian–Conrey–Heath-Brown), and ν < (3θ−1)/4 is Steuding's 2002 line ✓.
  Why SL's own §3.4–3.6 routes are weaker near 1/2 (not an error in E1): the L¹ multiplier (15) gives the per-n size H(N0/n)^σ; at σ = 1/4 this exceeds the
  x-space size H(pq/(nT))^{1/4} by exactly (T/H)^{1/2} (the stationary point of the u-integral at |t| ≍ X0 is lost); §3.6's Fourier Cauchy–Schwarz loses √(1+λ).

## 3. Other branches and every requirement at θ ↓ 1/2 (SL's own text)
| step (SL) | E1's claim | check against SL | at θ = 0.501, ν = 3/25 |
|---|---|---|---|
| ε = +1 primal (l. 464–465, 541–542) | ≪ (H/T)^A | Ψ_j at frequency (T+2πx)/H ≥ T/H: ≪ (H/T)^A; absolute primal sum ✓ | ✓ any θ < 1 |
| K-Bessel branch (l. 537–540) | O(H(y²/T)^σ) | ΓΓ exponentially decaying, move u right (no poles), Σ_n n^{−σ}, Σ_{p,q}(pq)^{σ−1} ≪ y^{2σ} ✓ | exponent −0.76σ ✓ |
| Estermann residues (l. 444–447, 709–728) | Young's diagonal, ν < 1/2 | j = 0 gives SL's (H√π/g)[…]; j ≥ 1 residues ≪ H(T/H²)^jL^4 = o(H/L) absolutely ✓; Lemma 2.2 needs only ν < 1/2 (Remark l. 336–338) ✓ | T^{1−2θ}L^4 → 0 ✓ |
| Taylor remainder (11) (l. 388–405) | J > (1+ν−θ+δ)/(2θ−1) | T^{1+ν+ε}(T/H²)^J ≪ HT^{−δ} ⟺ that J ✓ | at θ−η = 0.5005: J ≥ 621 |
| off-band kernel (l. 389–391) | (H/(T+x))^A | non-stationary y-phase ✓ | A > (1+ν−θ)/(1−θ) ≈ 1.24 ✓ |
| compact windows (l. 730–758) | width s0 = T^{θ−η} | w_D = φ_{s0}*Σ_{j<D}(−s0²/4)^jw^{(2j)}/j!; deficit ≤ u^D/D!; ‖w−w_D‖∞ ≪ (s0L/H)^{2D}; Σ|coeff| = O(H/s0) ✓ | η < 0.001, D ≥ 621 ✓ |
| shift derivatives (l. 723–725) | L^{O(1)} | Cauchy circles radius ≍ 1/L; all bounds uniform for |α|,|β| ≤ C/L ✓ | ✓ |
| (Emean) (l. 812–815) | Thm 1.2 for Q = 1, 1+x^j | both have Q(0) = 1, same P, R ✓ | ✓ |
| Lemma 4.2, Thm 1.3 (l. 823–928) | no θ condition | Riemann–von Mangoldt N(T,H) = (HL/2π)(1+o(1)) and Backlund–Jensen O(L) hold for H = T^θ, any θ > 0; only Thm 1.2 enters ✓ | ✓ |
| Thm 5.1 (l. 932–963) | Wang for 0 < λ < θ < 1 | Wang v2 TeX l. 229–240: Thm 2.2 for fixed 0 < λ < θ < 1, error O_g(H + T^λL²) = o(HL) ✓ | ✓ |
No step of Thm 1.2 (after replacing §3.4–3.5), Thm 1.3 or Thm 5.1 uses θ ≥ 0.527, ν < ν_A/ν_B or H ≥ T^{1/2+c} beyond θ − η > 1/2 (fixed η) ✓.
  Note (asymptotic only): the actual shifts are α = β = −R/L (+ derivative shifts), so |d| ≤ 1/4 needs L > 8R ≈ 71 (R = 8.85); harmless for a liminf.

## 4. Falsification (smallest failing cases) and numerics
Exponent table at the edge, θ = 0.501, ν = 3/25 (compact windows at θ' = θ − η = 0.5005):
| quantity | requirement | value |
|---|---|---|
| dual (new bound) | 1/2 + 2ν − 3θ'/2 < 0 | −0.01075 (−0.0115 at θ) ✓ |
| Taylor (11) | J(2θ'−1) > 1 + ν − θ' + δ | J ≥ 621 (fixed) ✓ |
| off-band / ε = +1 | A(1−θ') > 1 + ν − θ' | A ≥ 2 ✓ |
| K-branch | σ(2ν−1) < 0 | −0.76σ ✓ |
| j ≥ 1 residues | (T/s0²)L^4 = o(1/L) | T^{−0.001}L^5 → 0 only for L ≳ 7·10⁴ (asymptotic; SL's own remark l. 921–927 is of the same kind) ✓ |
| compact window | 2Dη > 1 + ν − θ + δ | D ≥ 621 ✓ |
| detector | κ(P,Q,R,ν) > k*(θ) | κ1 − k*(1/2) = 0.0078544 > 0, k* decreasing ✓ |
- Where it genuinely fails: ν ≥ (3θ−1)/4 (exponent ≥ 0; at equality the bound gives only T^ε·H) and θ ≤ 1/2 (T/H² ≥ 1: Taylor step and N0 ≥ pq). E1 claims neither ✓.
- d → 0: strip (|Re d|/2, 1/4) stays nonempty; Y_0 ≍ log x is integrable against x^{2σ−1}; a_±(X,d) continuous at d = 0 (Γ(1/2) ≠ 0) ✓.
- N0 < 1 (e.g. p = q = 1, all θ > 1/2): every n is in the tail; bound N0^{K/2} — no separate argument needed (SL's mod-1 convention l. 533–536 is unnecessary) ✓.
- Uniformity in j < J: constants C_j grow (ψ_j ≈ e^{−y²}(2y²)^j/j!, L¹ norm ≈ 2^j/√(πj)) but J is fixed for fixed θ ✓ — liminf only, as E1 states.
Numerics (all < 300 MB; peaks 17–62 MB):
- voronoi_check.py/.log (mpmath, dps 22): END-TO-END identity primal = Estermann residues + Σ_n dual (K-branch + −πB branch) for α = 0.03+0.02i, β = 0.05−0.01i,
  (p,q) = (1,3), A = complex Gaussian bump (ξc = 60, w = 6), 24 dual terms: |S| = 45.356107, relative error 1.7e−19 ✓. This checks (12), (14), Step 1's
  B = sin(πd/2)J_d + cos(πd/2)Y_d, the factor −π, the prefactor 4π/(gpq) and the branch/phase assignment at once (E1's own test compared only (14) with Step 1).
  The K-branch is ≈ e^{−56} here (not separately exercised; it is SL's unchanged bound). Near-coalescent residue case not run (time).
- perterm_check.py/.log (numpy; exact Hankel form, Gauss–Laguerre a_±, validated vs mpmath J_d, Y_d to ≤ 2.2e−2 worst case; real-d cross-check reproduces
  E1's e1_dual_size.log exactly: ratio 0.545, S/pred 0.612 at T = 25600, pq = 256): complex shifts α = 0.05+0.1i, β = −0.03+0.07i, j ∈ {0,1,2},
  T ∈ {1e5, 1e6}, θ ∈ {0.505, 0.55}, pq ∈ {6, 210}, N0 ∈ [1.9, 187]: max_{n≤N0}|𝒯_j(n)|/[(T/H²)^jHX0^{−1/2}] ∈ [0.19, 0.66];
  tail max at n ≥ 8N0 relative to that: 2.9e−6 … 1.1e−3 (j = 2 slowest: wider Ψ_2, j-dependent constant); S = Σd(n)|𝒯| over (T/H²)^jH(pq/T)^{1/4}N0^{3/4}:
  0.11 … 2.9 (grows with j, flat in T). No case exceeds the claimed scaling. (T = 1e7 rows cut by the 1200 s cap under background QoS.)
- constants_check.py/.log (independent route: Chebyshev-form Q, Gauss–Legendre 768/1536 nodes at 50/70 digits; Arb downstream): c = 3302.2869232044275656390124,
  κ1 = 0.08447790796974840320324474(±3e−27), c(1/2) = −0.1659638503649857462922256, k*(1/2) = 0.07662355506857570169601772,
  h(1/2;κ1) = 0.005535304238640537058439(±4e−25), dh/dθ(1/2) = 1.09746457895105; ν = 1249/10000, R = 173/20: κ = 0.08716133097499533713385541,
  h = 0.007429098196607165609897; h(θ;κ1) at θ = 0.501, 0.51, 0.527, 0.55, 0.6, 0.75: 0.0066312, 0.0163521, 0.0340454, 0.0566885, 0.1013780, 0.2069288 —
  all equal to E1 §4 to every printed digit ✓.

## 5. Reviewer (most delicate step: Steps 4–5, pasted inline with exact definitions)
VERDICT: CORRECT. IBP cost (T/X0)H^{−1} per step, uniform even for X0 < 1; min{1, c^{−K}(N0/n)^{K/2}} ≪_K (1+√(n/N0))^{−K}; n-sum needs K/2 > 3/4 + ε
(K ≥ 3 suffices), both regimes; pq cancellation exact; "ν < (3θ−1)/4 is where this bound gives a power saving, not a necessary condition for the dual to be small".

## 6. Verdict
Every step of E1 §2 and every "other branch" claim was re-derived from SL's TeX; no error found. Minor wording points (no effect): (i) the inverse-Mellin
identity for 𝒢 on Re u < 1/4 is only conditionally convergent — the proof uses Parseval, which is valid; (ii) Step 3's "X0 ≥ … → ∞ (Step 2 applies, X ≥ 1)":
the Hankel bounds hold for all X > 0, so ν < 1/2 is not needed there; (iii) E1 tested only j = 0 profiles — j = 1, 2 now tested; (iv) the result is purely
asymptotic (|d| ≤ 1/4 needs L > 8R; j ≥ 1 residues need T^{1−2θ}L^5 → 0).
- perterm_check2.py/.log (the rows cut from perterm_check): d = 0 EXACTLY (α = 0.01+0.02i, β = −α): ratios 0.507 (j = 0), 0.509 (j = 1) at T = 1e5, θ = 0.505,
  pq = 30, N0 = 26.7; 0.350 at T = 1e6, θ = 0.55, N0 = 7.54; tails ≤ 1.2e−4 at 8N0 ✓. N0 < 1 (p = q = 1, T = 1e6): |𝒯(n)|/(HX0(1)^{−1/2}) at n = 1, 2, 5:
  N0 = 0.251: 1.2e−3, 2.0e−6, 1.2e−14; N0 = 0.063: 9.5e−12, then roundoff floor (≈ 5e−15); N0 = 0.001: floor — super-polynomial tail decay as claimed ✓. Peak 57 MB.

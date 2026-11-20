# FULLPROOF — chain of the 0.3449 claim (two Dirichlet L-functions), with source lines

Sources (fetched 2026-11-20, latest versions): Wu arXiv 1802.09704v2 (TeX, single file: line numbers = file lines); Lamzouri arXiv 2610.26841v2 (IsolatedZeros2.tex; file line − 2 = TeX line, same numbering as DERIVATION). Ours: Riemann/Deep/IsoZero/DERIVATION.md (§k).

## A. Statement (Lamzouri Thm 1.4, l.248–249; ours: DERIVATION §8)
χ_1 ≠ χ_2 primitive mod q, χ_1(−1) = χ_2(−1) (same parity is a hypothesis, l.249, l.838), b_j real ≠ 0, θ_j = −arg(ε_{χ_j})/2 (l.845).
Λ_F(s) = (q/π)^{(s+a)/2}Γ((s+a)/2)F(s) = conj Λ_F(1 − s̄) (l.852–855): needs one common a, i.e. the same parity.
N_F(T) ~ (T/2π) log T (l.856–858, BH (3.12)); N(T,χ_j) = (T/2π)log(qT/2π) − T/2π + O(log T) (l.861–863).

## B. Reduction (Pro:BombieriHejhal, l.864–870; proof l.872–964, audited in IsoZero/AUDIT.md)
H_j(T) = {T < γ_n ≤ 2T: min neighbouring gap ≤ 2π/((log T)(log log T)^{1/4})}, with multiplicity (l.864–866).
#odd-multiplicity zeros of F(½+it), t ∈ [T,2T] ≥ (1 − ε)(T/2π) log T − |H_1(T)| − |H_2(T)| (l.868–870).
Use (ours, §8.1): a zero counted in I_1 of the block multiset has both ordinate gaps > 2πc/log T inside the block; for log log T > c^{−4} this exceeds the H-threshold, so I_1 ∩ H_j ⊂ {first, last zero of block}. |H_j| ≤ N_block − |I_1| + 2.
N_block = N(2T,χ) − N(T,χ) = (T/2π) log T + O(T). So |H_j| ≤ (1 − δ_j + o(1))(T/2π) log T with δ_j := liminf |I_1^{(j)}|/((T/2π)log T).
Proof of Thm 1.4 (l.1016–1029): dyadic summation of blocks (T/2^k, T/2^{k−1}], k ≥ 1: Σ (T/2^k/2π) log(T/2^k) = (T/2π) log T + O(T); with N_F(T) ~ (T/2π) log T the proportion is ≥ δ_1 + δ_2 − 1 − ε.
Order: ε fixed, T → ∞ (T ≥ max(T_0(ε), T_1(ε), T_2(ε))), then ε → 0. Blocks below T_0 are O(1).

## C. Per-character isolated proportion δ_j = p(c) (DERIVATION §2 on the block multiset)
C1. §2 (Hilbert-space inequality (E), (E′), §2.2–2.5) uses only: a finite conjugation-invariant multiset Z ⊂ ℂ, the kernel K = f̂, and ΣK(z−s)². Z_T^{(j)} (l.1018–1027) is conjugation-invariant: Λ(ρ,χ)=0 ⇒ Λ(1−ρ̄,χ)=0 (conjugate the FE Λ(s,χ)=ε_χΛ(1−s,χ̄)), same ordinate & multiplicity. ✓
C2. Pair correlation: Lem:GaPa (l.968–975; "a dyadic variant of their Lemma 3.1", l.967) with r = f*f, then La Lemma 3.2's argument removes the weight w(u) = 4/(4−u²): Lem:ConstructionKernelDirichlet (l.977–997): ΣK_ε² over the block = (C_{η_ε} + o_{ε,q}(1))(T/2π) log T. Constant: r(0) + 2∫_0^1 αr(α)dα = ∫f² + ∬|u−v|f(u)f(v) = C(f) (DERIVATION §1.2) ✓; log(qT) vs log T: same for fixed q ✓.
C3. Off-line term: hinge bound, §D below (replaces Lem:AnalogueConreyL, l.1010–1014, and LemmaWu, l.999–1008).
C4. Result: (1 − λ)|I_1| ≥ (2 − λ − C(f) − α_ρ B(ρ) − o(1)) N_block, i.e. δ_j ≥ p(c), the ζ number (all inputs χ-free).

## D. Hinge bound for fixed L(s,χ) (DERIVATION §8.2; replaces Lem:AnalogueConreyL)
D1. Wu Prop. 3.1 (1802.09704v2 l.388–413), verbatim content: χ primitive mod q, log q = o(log T), α = a/L_χ, β = b/L_χ, a, b ∈ ℂ, a, b ≪ 1; η > 0, Δ = T^{1−η}, y = T^θ, s_0 = ½ + iw, T ≤ w ≤ 2T;
 g(α,β,w) = (iΔ√π)^{−1}∫_{(1/2)} e^{(s−s_0)²/Δ²} L(s+α,χ)L(1−s+β,χ̄)B(s,χ)B_1(1−s,χ̄) ds, B_1(1−s,χ̄) = Σ_{n≤y} χ̄(n)ā(n)n^{−(1−s)} (l.392–396);
 "uniformly in a, b, and w" (l.398): g = Σ_{h,k≤y} χ_0(hk)(h,k)^{[α+β]}a(h)ā(k)h^{−1−β}k^{−1−α}[L(1−α−β,χ_0)(2πhk/(qw(h,k)²))^{α+β} + L(1+α+β,χ_0)] + O(T^{−ε_θ}) (l.400–401);
 (B1): ε_θ > 0 for θ < 4/7 when a(n) has the form (B) of Thm 1.1 (l.410), i.e. μ(n)(F_0 + F_1·(F_2*F_3))(n), F_i separable in 𝔽 = {F ≪ x^ε, F′ ≪ 1/x} (l.212–216; separability l.157–158).
 On Re s = ½: (s − s_0)² = −(t − w)², ds = i dt, so g = Gaussian mean (width Δ) at height w; with α = β real, L(½+α+it,χ)L(½+β−it,χ̄) = |L(½+α+it,χ)|² and B_1(½−it,χ̄) = conj B(½+it,χ). ✓
D2. Our mollifier a(n) = μ(n)P(log(y/n)/log y), P(0) = 0, P(1) = 1: F_0(x) = P(1 − log x/log y) = Σ_i c_i (log x/log y)^i, F_0(mn) splits binomially into products of (log x/log y)^i; on 1 ≤ x ≤ y (the only range used) each is ≤ 1 with derivative ≤ i/(x log y) ≪ 1/x. It is Wu's own Conrey part (l.265–266: a(n) = μ(n)(P_1(log(y/n)/log y) + P_2(·)Σ_{p|n}P(·)), case P_2 = 0), to which he applies (B). ✓
D3. Typos (authenticated in Wu's TeX):
 (i) exponent: l.400 (and l.202 Thm 1.1 (1.10), l.422 (x3.8)) print (h,k)^{α+β}; the diagonal of the twisted moment gives g^{1+α+β}h^{−1−β}k^{−1−α} (g = (h,k); re-derived: mh = nk, m = k′l, n = h′l); Wu's own α = β = 0 case (x3.4) l.405 has 1/[h,k] = (h,k)/(hk) and his Σ(α,β) l.1353 has (h,k)^{1+α+β}. Typo ✓.
 (ii) sign: I_R(Q,χ) is defined at V(½ + it + R/L_χ) (l.171) but evaluated at a = b = −R (l.203, l.426), and the Levinson mollifier ψ(s) = Σχ(n)a(n)n^{−s−R/L_χ} (l.262) is the left shift σ_0 = ½ − R/L_χ. Typo ✓. Irrelevant: Prop. 3.1 is used directly, uniform in the sign.
 (iii) further lower-order slips (new; harmless): (1.10) l.203 has (2^{1+α+β}−1)/(1+α+β) where ∫_T^{2T}(w/T)^{−α−β}dw/T = (2^{1−α−β}−1)/(1−α−β) (only this one reproduces (1.12)'s −1 + 2log2); (x3.4)/(1.12) have + c_q where L(1±z,χ_0) = (φ(q)/q)(±1/z + γ + c_q + O(z)) gives 2c_q (classical mean square of L(½+it,χ)); (1.12) has 2πHK for 2πhk. All O(1/log T) relative: no effect on leading constants.
D4. Shift: a = b = R_T := 2πρ L_χ/log T → 2πρ (real, bounded), so α = β = 2πρ/log T, σ⁺ = ½ + 2πρ/log T, ψ⁺(s) = Σ μ(n)χ(n)P[n]n^{α}n^{−s}, ψ⁺(σ⁺+it) = B(½+it,χ) = A_χ(t). g = Gaussian mean of |L(σ⁺+it,χ)A_χ(t)|². ✓


## E. Constants and certificate (DERIVATION §2.5, §3.7, §8.3; cert.py hinge_p, cert8.py)
E1. Main term at the shift, term by term (mainterm_sym.py, logs/mainterm_sym.log; sympy, general quartic P with P(0)=0, P(1)=1):
 Mellin residues → S_q(α,β) ~ (q/φ(q))(θL)^{−1}∫(P′+BθP)(P′+AθP) (Euler factor A_p(0) = 1 for p ∤ q; p | q removed → q/φ(q));
 Term 1 = L(1+α+β,χ_0)S_q(α,β) ~ ∫(P′+AθP)(P′+BθP)/((A+B)θ); Term 2 = (qw/2π)^{−α−β}L(1−α−β,χ_0)S_q(−β,−α) ~ −e^{−(A+B)}∫(P′−AθP)(P′−BθP)/((A+B)θ);
 q/φ(q) cancels in each term (the reviewer's "cancel at leading order" ✓); c(−R,−R) = c_1 (Conrey/BCY) ✓; c(R,R) = e^{−2R}c_1 = c_+ ✓; c(A,B) = e^{−(A+B)}c(−B,−A) ✓;
 P = x, R → 0: 1 + 1/θ (LemmaWu l.1004) ✓; inf_P: P = sinh(bx)/sinh b, ∫(P′+bP)² = b(1+coth b), c* = 1 + (1−e^{−2R})/(e^{2θR}−1) ✓.
E2. Numeric check of E1 with exact finite sums (mainterm_num.py, logs/mainterm_num.log; θ = ½, y = 10³…10⁶, q = 1 and 5, A = B ∈ {−1, 1, 2πρ}):
 M/c − 1 ≈ −4.5% (y = 10³) → −2.2% (y = 10⁶) (∝ 1/log y); M_5/M_1 − 1: −1.1% → −0.55% (→ 0); hinge case A = 2π·1.2036, sinh P: M − 1 = 3.1e−4 → 3.9e−4 (limit c* − 1 = 5.2e−4, from below).
 The finite expression satisfies M(A,A) = e^{−2A}M(−A,−A) exactly (FE symmetry of Wu's main term).
E3. Littlewood (DERIVATION §8.2): f = L(s,χ)ψ⁺(s) = Σχ(m)b_m m^{−s}, b_1 = P(1) = 1, b_m → Σ_{d|m}μ(d) = 0 (m ≥ 2 fixed), |b_m| ≤ d(m)max|P| m^{α}: |f(2+it) − 1| → 0 uniformly ✓;
 horizontal edges O(log T) (Re f(2+iT) > ½, Jensen) ✓; concavity: ∫_T^{2T} log|f(σ⁺+it)|dt ≤ (T/2) log((1/T)∫|f|²) ✓;
 2π∫_{σ⁺}^{2}ν(σ)dσ ≤ ∫log|f(σ⁺+it)| − ∫log|f(2+it)| + O(log T) ⇒ Σ_{T<γ≤2T, β>σ⁺}(β − σ⁺) ≤ (T/4π)log c_+ + o(T) ✓ (zeros of ψ⁺ only add to the left side).
 β − σ⁺ = (2π/log T)(y − ρ), y = (β − ½)log T/2π ⇒ Σ(y − ρ)^+ ≤ (B(ρ) + o(1))(T/2π)log T, B(ρ) = (1/4π)log c*(2πρ) (θ → 4/7) ✓.
E4. Unsmoothing (Wu l.425–428 "as in BCHB §3"): ∫_T^{2T}g(w)dw = ∫|f|²Φ, Φ = Gaussian-smoothed 1_{[T,2T]}, Φ = 1 − O(e^{−log²T}) on [T + Δlog T, 2T − Δlog T]; the two edge pieces (length T^{1−η}log T) are o(T) by Prop. 3.1 at a slightly smaller/larger T (y = T^θ = T_1^{θ′}, θ′ < 4/7). For an upper bound |f|² ≥ 0 suffices ✓.
E5. Certificate rerun (copy in .scratch/lcaudit/run, memcap): output identical to IsoZero/logs/cert8.log (logs/cert8_rerun.log): c = 10⁻⁴: p ≥ 0.672491254, 2p − 1 ≥ 0.344982508.
E6. Independent recomputation (indep_p.py, logs/indep_p.log; K by 400-node Gauss–Legendre and 40-digit tanh-sinh quadrature instead of closed forms; sup's by 1e−5 grids + golden section):
 C(f_0) = 1.32749929632059 (quadrature = closed form); c = 10⁻⁴: λ = 1.22382e−7 (β_T = −0.667 < 0, so λ_C = 2(1 − K(2c))), ρ_max = 1.2236206, ρ = 1.2036, α_ρ = 0.6702740537 at y* = 1.38857, B = 1.403822e−5:
 p = 0.672491254144, 2p − 1 = 0.344982508288 ✓ (cert: ≥ 0.672491254, ≥ 0.344982508). c = 10⁻³: 2p − 1 = 0.344974576351 ✓. ρ scan: ρ = 1.210 gives 0.672491327 (cert's ρ valid, not optimal; gain 7e−8).


## F. Falsification (zeros_F.py, logs/zeros_F.log; q = 5, Arb L-values, step 0.04, t ≤ 1000)
F1. Setting: χ_a = (5,2) (χ(2) = i), χ_b = (5,3) = χ̄_a, both odd; e^{iθ_j} = ε_j^{−1/2}; e^{iϑ(t)}e^{iθ_j}L_j(½+it) real to 5.7e−13 ✓. b = (1,1) is the Davenport–Heilbronn function: |F(0.808517 + 85.699348i)| = 1.25e−6 (known off-line zero; 6-digit input) vs 1.40 for b = (1,−1) ✓.
F2. Sign changes of Z_F in dyadic blocks vs the main term (T/2π)log(qT/2πe) (calibration: L_a alone 43/42.7, 99/99.2, 226/226.0, 507/507.1):
 (62.5,125]: DH 0.913, F(1,−1) 0.562; (125,250]: 0.938, 0.625; (250,500]: 0.956, 0.655; (500,1000]: 0.929, 0.682; (1000,2000] (zeros_F2.py): 0.899, 0.695 (L_a: 1124/1124.5). All ≥ 0.345: no conflict (a liminf claim cannot be refuted by finite data; this only checks the setting).
F3. Same parity is necessary. χ_a (odd) + χ_q = (5,4) (even): ϑ_odd − ϑ_even ≡ π/4 (mod π) (to 4 digits at t = 200…500, Stirling), so F(½+it) = e^{−iϑ_odd}(b_1Z_a + b_2e^{iπ/4+o(1)}Z_q) vanishes only at common zeros of L_a, L_q; numerically min_{t∈[10,500]} |F|/(|L_a| + |L_q|) = 0.383: no zero on the line. The FE Λ_F(s) = conj Λ_F(1 − s̄) (l.852) fails, and with it the reality of Z_F used in Pro:BombieriHejhal. Mod 8 has no admissible pair at all (its primitive characters χ_8 even, χ_{−8} odd). The hypothesis χ_1(−1) = χ_2(−1) is in Lamzouri's statement (l.249) and in ours ✓.
F4. Other hypotheses, each necessary for the FE/reality step: b_j real (complex b_j destroy reality on the line); χ_j primitive mod the same q (different conductors give different (q/π)^{s/2}); χ_1 ≠ χ_2 (else F = cL).

## G. Citation gap found: the "dyadic variant" of GaPa (2610.26841 l.967–975), and a repair
G1. Lamzouri: "The following is a dyadic variant of their Lemma 3.1" (l.967), stated for T < γ, γ′ ≤ 2T with log T/2π. GaPa's Lemma 3.1 (Glasnik Mat. 60(80) (2025) 207–227; accepted manuscript pGM7479.pdf, its l.788–806 in pdftotext) is for 0 < γ, γ′ ≤ T, normalisation log(qT)/2π, weight W(u) = 4/(4 − u²), uniform for 1 < q ≤ log T. The dyadic form is not proved in either paper; DERIVATION §8.2 uses it (via Lem:ConstructionKernelDirichlet). Plausible (same Montgomery method on [T,2T]), but not a citation-level fact.
G2. Repair without it (full range, only GaPa's actual statement):
 (a) Pair correlation on Z_T^{full} = {−i(ρ − ½)L′: 0 < γ ≤ T}, L′ = log(qT)/2π: GaPa Lemma 3.1 with r = f*f plus La Lemma 3.2's removal of W (verbatim as for ζ with BGST, which has the same W) gives ΣK_ε² = (C_{η_ε} + o(1))(T/2π)log(qT).
 (b) Full-range hinge from block hinges: blocks (T_k, 2T_k], T_k = T/2^k ≥ T^{1/2}, each with its own Wu mollifier y_k = T_k^θ and ρ_k := ρ log T_k/log T: Σ_block(y − ρ)^+ = (log T/log T_k)Σ_block(y^{(k)} − ρ_k)^+ ≤ (B(ρ_k) + o(1))(T_k/2π)log T. k ≤ K_0: B(ρ_k) → B(ρ); K_0 < k: (y^{(k)} − ρ_k)^+ ≤ (y^{(k)} − ρ/2)^+, ≤ B(ρ/2)·2^{−K_0}N; below T^{1/2}: O(T^{1/2}log²T). K_0 → ∞: Σ_{0<γ≤T}(y − ρ)^+ ≤ (B(ρ) + o(1))N(T,χ).
 (c) §2 on Z_T^{full} (literally the ζ argument): liminf |I_1^{full}|/N(T,χ) ≥ p(c).
 (d) Pro:BombieriHejhal on each block k ≤ K (T_k ≥ T^{1/2}), summed: Σ_k(1 − ε)(T_k/2π)log T_k = (1 − ε)(T/2π)log T − O(T). A zero in H_j(T_k) has a full-list neighbour gap ≤ 2π/((log T_k)(log log T_k)^{1/4}) ≤ 4π/((log T)(log log T^{1/2})^{1/4}) < 2πc/log(qT) (T large), so it is not in I_1^{full} (except the top zero): Σ_k|H_j(T_k)| ≤ N(T,χ_j) − |I_1^{full,(j)}| + O(1).
 Hence #odd on-line zeros of F in (T^{1/2}, T] ≥ (δ_1 + δ_2 − 1 − ε − o(1))(T/2π)log T with δ_j = p(c): the same 2p − 1. The same repair applies to Lamzouri's own 31%.


G3. Reviewer (logs/review.log): GAP, conditional on the isolation lemma: steps 1–3 correct; the isolation lemma (GaPa ⇒ |I_1|/N ≥ p(c)) "not supplied" — it was not pasted by design: it is DERIVATION §2 (audited) plus La 2609.02882 Lemma 3.2's proof, whose only input is BGST Lemma 5 (LamzAudit2 S3, re-derived), and GaPa Lemma 3.1 has BGST Lemma 5's exact form (same weight W, same class of r; error o(1) instead of O(1/√log T), which Lemma 3.2 does not need beyond o(1)). Correction adopted: BH's closed blocks [T_k, 2T_k] share endpoints; zeros at the K = O(log T) endpoints are counted at most twice: subtract O(log T). Residual: La's Lemma 3.2 proof is applied to L(s,χ) by inspection of LamzAudit2's re-derivation, not re-derived here.


G4. The pair-correlation half of the isolation lemma, made explicit for L(s,χ) (closes the reviewer's condition together with DERIVATION §2): La 2609.02882 Lemma 3.2's weight removal (LamzAudit2 P2: r_{δ,T} = Q_δ − Q_δ″/(4log²T), Q_δ = f_δ*f_δ) with log(qT) in place of log T:
 W(u) = 4/(4 − u²), u = ρ − ρ′, z = i u log(qT)/2π ⇒ 1/W = 1 − u²/4 = 1 + π²z²/log²(qT); FT of Q″ is −4π²z²Q̂, so r := Q − Q″/(4 log²(qT)) has r̂ = (1 + π²z²/log²(qT))K² = K²/W at every pair.
 GaPa Lemma 3.1 applied to Q and to Q″ separately (both C_c^∞((−1,1)), even, real): ΣK² = Σ r̂W = (T/2π)log(qT)[Q(0) + 2∫_0^1 αQ − (Q″(0) + 2∫_0^1 αQ″)/(4log²(qT)) + o(1)] = (C(f_δ) + o(1))(T/2π)log(qT), C(f) = ∫f² + ∬|u − v|ff. No ζ-specific input. ✓


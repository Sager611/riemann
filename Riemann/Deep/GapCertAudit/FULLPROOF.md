# gapcertaudit FULLPROOF: the chain behind Λ > 3.34 (= 167/50), with source lines

Sources (TeX e-prints, read 20 Nov 2026, deleted after use): BM = Bui–Milinovich arXiv 1410.3635v3 (main.tex, 758 lines);
BBLR = Bettin–Bui–Li–Radziwiłł arXiv 1609.02539 (1110 lines). "l." = line of that file.
Notation: L = log(T/2π), c = κπ, κ = 167/50, r = 2, ϑ = 2499/10000, y = T^ϑ, u = tL.
F(t) = e^{ivtL}G(α*;t), G(α;t) = ζ(½+α1+it)ζ(½+α2+it)M(½+α3+it), α* = (0, ic/L, ic/(2L)); real P, so
conj F = e^{−ivtL}G̃(β*;t), G̃(β;t) = ζ(½+β1−it)ζ(½+β2−it)M(½+β3−it), β* = conj α*.

## (a) Reduction (V) ⇒ "Σ_j a_jκ^{2j}c_j ≥ 0" under Λ ≤ κ (BM l.164-202 with Wirtinger replaced)
1. Λ ≤ κ: critical ordinates in [T, 2T] satisfy γ_{n+1} − γ_n ≤ (1+o(1))2c/L (2π/log γ = (2π/L)(1+O(1/L))).
   F vanishes at γ_n and γ_n − c/L; inside (γ_n, γ_{n+1}) of length d ∈ (c/L, 2c(1+o(1))/L] the point γ_{n+1} − c/L
   splits it into pieces c/L and d − c/L. So consecutive real zeros t̃_n of F in [T, 2T] are ≤ (1+o(1))c/L apart:
   BM l.180-183 (the mollifier only adds zeros). Checked.
2. Fix ε > 0, λ = (1+ε)κ/L; for T ≥ T_0(ε) every gap ℓ_n = t̃_{n+1} − t̃_n ≤ λπ, so the rescaled gap g_n = ℓ_n/λ ∈ (0, π]
   (zero-length gaps: skip). f(u) = Re F(t̃_n + λu) and Im F(t̃_n + λu) are real-analytic, vanish at u = 0, g_n: (V) applies.
   d/du = λ d/dt: ‖f^{(j)}‖²_u = λ^{2j−1}∫|·^{(j)}|²dt, [f′f″]_u = λ³[·], [f‴f″]_u = λ⁵[·]. Sum Re + Im (real-bilinear):
   Σ_j a_jλ^{2j}∫_{t̃_n}^{t̃_{n+1}}|F^{(j)}|² + φ12λ⁴[B12]_{t̃_n}^{t̃_{n+1}} + φ32λ⁶[B32] ≥ 0 (after multiplying by λ),
   B12 = Re(F′ F̄″), B32 = Re(F‴ F̄″): single smooth functions of t, so Σ_n telescopes to B(t̃_N) − B(t̃_1). Checked.
3. End terms: t̃_1 − T, 2T − t̃_N ≤ 2c/L (BM l.194: ≪ 1). Weyl ζ(½+it) ≪ t^{1/6+ε}, Cauchy on circles of radius 1/L for
   the derivatives (≪ L^k t^{1/6+ε}), |M^{(k)}(½+it+ic/(2L))| ≪ (log y)^k Σ_{h≤y}|d_r(h)P[h]|h^{−½} ≪ T^{ϑ/2+ε}:
   |F^{(k)}|² ≪ T^{2/3+ϑ+ε}, ϑ < 1/4 ⇒ ≪ T^{11/12+ε}. So the end pieces ∫_T^{t̃_1}, ∫_{t̃_N}^{2T} and the brackets
   B(t̃_1), B(t̃_N) (with λ⁴, λ⁶ ≤ 1) are O(T^{11/12+ε}) = o(T L^{(r+2)²}). Checked (BM l.194 states the same bound).
4. With ∫_T^{2T}|F^{(j)}|² = c_jKTL^{(r+2)²+2j}(1+o(1)) (item (b)), λ^{2j}L^{2j} = ((1+ε)κ)^{2j}:
   Σ_j a_j((1+ε)κ)^{2j}c_j ≥ −o(1); T → ∞ then ε → 0 (polynomial in ε): Σ_j a_jκ^{2j}c_j ≥ 0. Sign: H := −Σ a_jκ^{2j}c_j/K'
   > 0 contradicts it. a_3 < 0 is harmless (only continuity in ε is used). Checked. Verdict (a): HOLDS.

## (b) Mean values ∫|F^{(j)}|², j = 3, 4 (BM Lemma main l.370-384, l.393-434; own derivation)
1. Exact identity: d/dt G(α;t) = iΣ_k∂_{α_k}G, so F^{(j)} = (iL)^j e^{ivtL}D_α^jG|_{α*}, D_α = v + L⁻¹Σ_{k≤3}∂_{α_k};
   conj F^{(j)} = (−iL)^j e^{−ivtL}D_β^jG̃|_{β*}. Hence ∫w(t/T)|F^{(j)}|² = L^{2j}D_α^jD_β^j I_{α,β}(M)|_{(α*,β*)},
   I as in BM l.364-367 (w under the integral, compact support: differentiation under ∫ is fine). BM l.395-397 is j = 1.
2. Main term: BM l.373, K' = A_{r+2}(log y)^{r²+4r}L⁴ŵ(0)/(2(r²−1)!((r−1)!)⁴) is shift-free; c(α,β) (l.377-381) depends
   on the shifts only through one exponential exp(Σ_kα_kℓ_k + Σ_kβ_km_k) (y- and (Ty^{−x1−x3}), (Ty^{−x2−x4})-powers).
   Own computation: Σ_kℓ_k = −[(log y)(x+x1+x2+x3+x4) + (log T − (x1+x3)log y)t1 + (log T − (x2+x4)log y)t2]
   (the (α2−α1)-terms cancel in Σ∂), the same for Σ_k m_k (the (β2−β1)-terms cancel). So L⁻¹Σ_k∂_{α_k} and L⁻¹Σ_k∂_{β_k}
   both multiply the integrand by −S′, S′ = (log T/L)S, S = ϑ(x+Σx_i) + (1−ϑ(x1+x3))t1 + (1−ϑ(x2+x4))t2, independent of
   the shifts: D_α^jD_β^j c = ∫(v − S′)^{2j}(integrand), exactly, any j (= BM l.422, l.431). (v − S′)^{2j} = (v − S)^{2j}
   + O_j(1/L) uniformly (S ≤ 3 on the domain): the change is O(TL^{(r+2)²+2j−1}). Checked.
3. Error: E = I − Main is holomorphic on the polydisc P = {|α_k − α*_k|, |β_k − β*_k| ≤ 1/L} (I: ζ's pole at ½ − it is
   far for t ∈ [T, 2T]; Main: entire) and |E| ≤ C T L^{(r+2)²−1} on P (Lemma main "uniformly for α_j, β_j ≪ L⁻¹",
   l.383, with constant (c + 1)/L; BBLR l.104 allows complex shifts ≪ 1/log T). Cauchy on P: |∂^nE| ≤ n!L^{|n|}sup_P|E|,
   so |(L⁻¹Σ∂_α)^a(L⁻¹Σ∂_β)^bE| ≤ a!b!C(a+2,2)C(b+2,2)·sup|E|, and D_α^jD_β^jE = O_j(TL^{(r+2)²−1}); times L^{2j}:
   O_j(TL^{(r+2)²+2j−1}) (= BM l.417 "each derivative adds a factor O(L)", which is j-free). Checked.
4. Mollifier coefficients with complex α3 on P: |h^{−α3}| ≤ y^{|Re α3|} = e^{O(ϑ log T/L)} = O(1), so a_h ≪ h^ε (BBLR l.108).
5. Unsmoothing: BM's w is supported in [1, 2], w^{(j)} ≪ T^ε (l.290). Lower bound: w_− ≤ 1_{[1,2]}, ŵ_−(0) ≥ 1 − η.
   Upper bound: a finite smooth partition 1_{[T,2T]} ≤ Σ_i w_i(t/T_i), T_i ∈ [T/2, 2T], each w_i supported in [1, 2],
   Σ_iT_iŵ_i(0) ≤ (1+η)T; F keeps its own L, Lemma main at scale T_i changes log T_i/L = 1 + O(1/L) and K' by
   (1 + O(1/L)). Hence ∫_T^{2T}|F^{(j)}|² = c_jK T L^{(r+2)²+2j}(1 + O(η) + o(1)), η → 0: the leading term (all (a) needs).
   BM's own h1 (l.200) is an unsmoothed ratio computed from the smoothed theorem: the same step. Verdict (b): HOLDS.

## (c) Lemma 2, sufficiency direction (own derivation; ExploreShape REPORT "Lemma 1–2")
1. f ∈ H⁴(0,g), f(0) = f(g) = 0, b_k sine coefficients, ω_k = kπ/g, A0 = f″(0), A1 = f″(g), d_k = (2/g)((−1)^kA1 − A0).
   Integrating by parts (boundary terms contain only f(0), f(g), sin(0), sin(ω_kg)): f″ ↔ −ω²b (sine); f‴ ↔ d_k − ω³b_k
   (cosine), constant term (A1 − A0)/g; f⁗ ↔ ω⁴b_k − ω_kd_k (sine). Parseval: ‖f‴‖² = (A1−A0)²/g + (g/2)Σ(d_k − ω³b_k)².
2. Decay: e_k := ω⁴b_k − ω_kd_k ∈ ℓ², so ω_kb_k = d_k/ω² + e_k/ω³ and d_k − ω³b_k = −e_k/ω_k are in ℓ¹: the cosine series
   of f′, f‴ converge absolutely, f′(0) = Σω_kb_k, f′(g) = Σ(−1)^kω_kb_k, f‴(0) = (A1−A0)/g + Σ(d_k − ω³b_k), f‴(g) likewise.
   Hence [f′f″]_0^g = (g/2)Σω_kb_kd_k and [f‴f″]_0^g = (A1−A0)²/g + (g/2)Σ(d_k − ω³b_k)d_k.
3. Collecting (each series converges as a combined series; p b², γd² diverge separately when a4 ≠ 0, so the termwise
   grouping must be kept): Q = (g/2)Σ_k[p(ω_k)b_k² + β(ω_k)b_kd_k + γ(ω_k)d_k²] + (a3 + φ32)(A1−A0)²/g with
   p = Σa_jω^{2j}, β = φ12ω − (2a3+φ32)ω³ − 2a4ω⁵, γ = a3 + φ32 + a4ω². Re-derived; equals the page.
4. If p(ω_k²) > 0 ∀k: each bracket ≥ ρ(x_k)d_k², ρ = γ − β²/(4p) = N/(4p), x = ω²; Σρ_kd_k² converges (ρ = O(1/x)).
   d_k² = (4/g²)(A1 + A0)² (k odd), (4/g²)(A1 − A0)² (k even): Q ≥ (2/g)[(A0+A1)²S_odd + (A1−A0)²(S_even + (a3+φ32)/2)].
   With φ32 = −a3 the last term vanishes. So p > 0 on [1, ∞) (x_k ≥ (π/g)² ≥ 1 for g ≤ π), S_odd ≥ 0, S_even ≥ 0 ⇒ Q ≥ 0.
   N = 4pγ − xq², q = φ12 − (2a3+φ32)x − 2a4x²: own expansion gives x⁵-coefficient 4a4² − 4a4² = 0 and x⁴-coefficient
   4a4(a3 + (a3+φ32) − (2a3+φ32)) = 0 identically; n0 = 4a0(a3+φ32) = 0, n1 = 4a0a4 − φ12², n2 = 4a1a4 + 2φ12a3,
   n3 = 4a2a4 − a3² + 4φ12a4 (φ32 = −a3). These are cert4.py's n (l.31-34) with c = a3+φ32, e = 2a3+φ32. Verdict (c): HOLDS
   (only sufficiency is used; the "iff" is not needed).

## (d) Main-term pipeline certJ.py vs GapBar (gtw.gamma_arrays_tw, gcert.h_series/uhat; audited by GapAudit)
- h_series_K = gcert.h_series with the closed recursion m_j = (e^z − j m_{j−1})/z for any k (GapBar: k ≤ 2) — same
  Taylor series of h(y) = ∫_0^{1−y}t^ke^{zt}dt (own check: h′ = −e^z(1−y)^ke^{−zy}). uhat_K = B1·Hank(hh)·B3ᵀ, entrywise
  identical to gcert.uhat's double sum. Twisted phases (±c/2) and the four signed terms: verbatim gtw's.
- τ-weights: g_m = Σ_{k+l=m}C(m,k)conv2(u_k,w_l) is (τ1+τ2)^m; for m ≤ 2 it is gtw's combos [(0,0,1)], [(1,0),(0,1)],
  [(2,0,1),(1,1,2),(0,2,1)]. S-kernels H_n = Σ_jC(n,j)sh^j(g_{n−j}): n = 0,1,2 give gtw's a, b = sh(a)+g1,
  cc = sh²(a) + 2sh(g1) + g2. ⟨S^m⟩ = Σ_iC(m,i)Q(H_{m−i}, X^i) is the binomial expansion of S = X + (R + τ1 + τ2) with
  X ↔ ϑx (Legendre multiplication T), sh ↔ R = ϑ(x1+..+x4): all multiplications commute, so exact for every m.
- Windows: entries M, N < D of conv2/sh/uhat use only indices < D (h_series has 2D terms): exact; X^q from T^q in the
  window nf = n2 + q0 + K + 4 ≥ n2 − 1 + q0 + K + 1: exact (as GapBar). Hermitian part: real φ, true Q real ⇒ Re of the
  symmetrised array is exact; imok only tests symmetry of Im (relevant for complex φ, as GapAudit noted).
- Lemma T (truncation), own derivation. h_k(y) = ∫_0^{1−y}t^ke^{zt}dt, z ∈ {0, ±ic}: for q ≥ 1,
  |h^{(q)}(0)| = |d^{q−1}/ds^{q−1}(s^ke^{zs})|_{s=1} ≤ Σ_iC(q−1,i)k^i|z|^{q−1−i} = (k+|z|)^{q−1}; |h(0)| ≤ 1/(k+1). Leibniz with
  phases |α1|, |α3| ≤ c/2: |u_k[a][b]| ≤ (c/2 + k + c + 1)^{a+b} ≤ β^{a+b}, β = 25 ≥ 3c/2 + 9 = 24.74. Fibre weights
  wr[a][b] = (a+r−1)!(b+r−1)!/(a!b!); Chu–Vandermonde Σ_{a1+a2=M}C(a1+r−1,r−1)C(a2+r−1,r−1) = C(M+2r−1,2r−1); four terms,
  Σ_kC(m,k) = 2^m: |g_m[M][N]| ≤ 4·2^m((r−1)!)⁴C(M+2r−1,2r−1)C(N+2r−1,2r−1)β^{M+N}. sh: C(M+2r−2,2r−1)(M+2r−1) =
  M·C(M+2r−1,2r−1), so each sh multiplies the bound by ≤ (M + N) (β ≥ 1); Σ_jC(n,j)2^{n−j}(M+N)^j = (2+M+N)^n ≤
  (2+M+N+4r)^n (the code's, weaker). ‖I^M‖ ≤ ϑ^M/M! (Volterra), ‖X^q‖ ≤ ϑ^q on [0,ϑ]: Cauchy–Schwarz gives the inner-product
  bound. Convexity (2+4r+M+N)^n ≤ 2^{n−1}((1+2r+M)^n + (1+2r+N)^n) and the region {M ≥ D or N ≥ D} ⊂ {M ≥ D} ∪ {M < D ≤ N}
  give certJ l.157-163 exactly; the cut at 700 (ratio βϑ/M < 1/40, remainder ≪ 10⁻¹⁰⁰⁰) is covered by the added 10⁻³⁰⁰.
  errs[m] = Σ_iC(m,i)tails[m−i]ϑ^{q0+i}‖φ‖², added to the radii before forming c_j and H. Verdict (d): HOLDS.

## (2) Code audit (cert4.py, certJ.py; run copies in .scratch/gapcertaudit/run/ES = git HEAD, GapBar symlinked read-only; deleted after use)
- cert4.py: a from the float witness by limit_denominator(10¹²), then exact Fractions (a0 += 2·10⁻⁵, φ12 += 2·10⁻⁶,
  φ32 := −a3); arb only from fmpq (exact) and arb.pi(). (1) p(1) > 0 exact; fmpq_poly.complex_roots (Arb-rigorous isolation,
  real roots with exact 0 imaginary part): real roots −5.441, 0.96554 < 1, complex pair 6.689 ± 5.143i ⇒ p > 0 on [1, ∞).
  (2) n3 > 0 exact; X_N := upper end of the largest real root of N (1.2856683898) ⇒ N > 0, ρ > 0 for x > X_N.
  (3) cover starts at lo = (π/√X_N).lower() ≤ π/√X_N (checked: yes, the cell list starts at gstar.lower()) and ends at
  π.upper() ≥ π; for g < lo every x_k > X_N (termwise positive). (4) Each cell is one arb ball G, x_k = (kπ/G)² encloses
  all g in the cell; k ≤ 61, dropped k ≥ 62 have x_k ≥ 62²(π/π.upper())² > X_N (positive): truncation is a lower bound.
  The comment "x_k >= 4" at l.66 is stale (harmless). S_even + (a3+φ32)/2 = S_even since a3 + φ32 = 0 exactly.
- certJ.py: θ, κ, v, A_n parsed by GC.dec (exact decimal → fmpq); a parsed as exact Fractions; the only floats choose the
  integer β (= 25, any β ≥ 24.74 is valid) and the precision (382 bits). a in amp_d64_k3.34.json = cert4's a (string-equal,
  checked); φ12, φ32 do not enter H (they telescope). Decision ok = H > 0 and c0 > 0 and imok and θ < 1/4 on Arb balls.

## (3) Falsification runs (one process at a time, memcap 290 MB; logs in logs/)
| run | result | time, peak |
|---|---|---|
| cert4.py 200 61 (copy) | JSON identical to the agent's; CERTIFIED True; S_odd(π) ≥ 1.0791e-5 (k ≤ 61), min cell lower bounds S_odd 4.0095e-6 (g ≈ π), S_even 0.016050 | 0.4 s, 35 MB RSS |
| certJ.py amp_d64_k3.34 (D = 72, 382 bits) | identical to the agent's JSON (all fields but timing): H = 9.3142330897e-35 ± 1.92e-46, H/c0 = 0.00136298995295 ± 5.25e-15 | 195 s, 167 MB memcap, 211 MB RSS |
| certJ.py same amplifier, D = 64 (361 bits): other truncation and precision | H = 9.31e-35 ± 5.38e-38, H/c0 = 0.001363 ± 1.78e-7, CERTIFIED True; every s_m, c_j ball overlaps D = 72's; mid differences (c4: 3.8e-6 rel.) < D = 64 radii (c4: 6.1e-6 rel.) | 204 s, 104 MB memcap, 163 MB RSS |
| certJ.py test amplifier (d = 8, κ = 13/4, Wirtinger a, D = 60) | identical to the agent's test JSON | 115 s, 81 MB |
| indep8.py (own, BM TeX only; extends GapAudit's indep7) same d = 8 test, n = 16 / 20 Gauss points | ⟨S^m⟩/⟨1⟩, m ≤ 8: n = 16 and n = 20 agree to ≤ 1.4·10⁻¹⁰; n = 20 vs certJ at D = 72: ≤ 5.5·10⁻¹⁴ for every m ≤ 8 (float limit); every D = 60 ball (relative radii 4·10⁻¹¹, 4·10⁻⁹, 3·10⁻⁷ at m = 6, 7, 8) contains the value | 17 s / 54 s, 61 / 126 MB |
| certJ.py test amplifier at D = 72 | radii ≤ 7·10⁻¹⁴ relative for all m ≤ 8; CERTIFIED (Wirtinger form, H/c0 = 0.0153629584984787) | 175 s, 159 MB memcap, 222 MB RSS |
| hiprec.py: D = 60 test with +160 bits (510) | every s_m ball unchanged ⇒ the D = 60 radii are the Lemma T tail bound, not arithmetic error | 125 s |
| majcheck.py: entrywise Lemma T on the certificate's arrays (κ = 3.34, D = 72, n ≤ 8) | max \|H_n[M][N]\|/majorant = 0.49 (n = 0, at (0,0)), ≤ 4.7·10⁻⁸ (n = 8); with the sharper (2+M+N)^n ≤ 0.49; for max(M, N) ≥ 64: ≤ 1.9·10⁻⁴⁶ | 174 s, 162 MB |
- Reading: the pipeline (τ-weights, sh-kernels, binomials) is right for every m ≤ 8 (an error would show at O(1)).
  Arb prints only certain digits, so a printed midpoint is rounded to its radius; the apparent "D = 60 vs D = 72"
  differences (1.4·10⁻⁷ at m = 8) are that rounding. Lemma T holds on every computed entry and is very conservative
  where it is used (10⁻⁴⁶ at the truncation edge): the tails are bounds, not estimates.
- FINDING (documentation): the page and R494 say crosscheck.py "matches the j = 3, 4 moments (0.0088906238, 0.020059772)".
  crosscheck.py has K = 4 (l.52): those numbers are ⟨(v−S)^3⟩/⟨1⟩, ⟨(v−S)^4⟩/⟨1⟩, i.e. only c_0, c_1, c_2 were
  cross-checked; c_3 = ⟨(v−S)^6⟩, c_4 = ⟨(v−S)^8⟩ had no independent check before indep8.py (now: agree, see table).

## (3b) Per-gap checks by another basis, and the shifts (vbasis.py, shifts.py; logs/vbasis*.log, logs/shifts.log)
- Basis sines sin(kπx/g), k ≤ K, plus x(g−x)P_m(2x/g−1), m < MP; Gram matrices of f^{(j)} and of the symmetrised brackets by
  192-node Gauss–Legendre at 50 digits; λ_min of the pencil (Cholesky + eigsy). MP = 2 (two cubics only) converges like 1/K
  (λ_min = 5.6·10⁻³, 3.1·10⁻³ at K = 6, 12): the minimiser's sine tail b_k ~ d_k/ω³ needs the polynomial jets.
- MP = 8, g = π: λ_min(cert) = 1.817086256·10⁻⁵ at K = 20, 30 and at MP = 12 (10 digits); orch_vnum's 1.82·10⁻⁵ (pure
  polynomials, degree ≤ 27). Other g (K = 12): 0.02407, 0.08503, 0.24289, 1.1020, 7.2603, 92.476, 3406.1 at g = 3.1, 3.0,
  2.7707, 2, 1, 0.3, 0.05 (orch_vnum: 0.0241, 0.0850, 0.243, 1.10, 7.26, 92.5 at the common g).
- λ_min(noshift) = −1.829137·10⁻⁶ = λ_min(cert) − 2·10⁻⁵; λ_min(sdp) = +1.92·10⁻⁹ at g = π; shifts.py: sdp has n3 = −7.8·10⁻¹⁰
  and S_odd(10⁻³) = −7.3·10⁻¹⁶, S_odd(10⁻⁴) = −9.4·10⁻¹⁷ (fails for tiny gaps); noshift has S_odd(π) = −1.15·10⁻⁶.

# FULLPROOF — the RH-free Breuillard–Varjú chain, with source lines (audit bvirraudit, 14 Nov 2026)

BV = arXiv 1810.13360, irreducibility-via-large-primes.tex. D = Riemann/Deep/BVIrr/DERIVATION.md. Each step: source, what is used, audit check.

## A. Inputs
- A1 (CT) every number field K: ζ_K(s) ≠ 0 for Re s > 31/32, s ≠ 1 (HPAudit D: ∃ continuation f, all K : Type, constant 31/32 absolute;
  no n_K, no exceptional-zero clause). With the functional equation (D §0 FE): nontrivial zeros have 1/32 ≤ β ≤ 31/32, so |ρ| ≥ 1/32.
  Checked: trivial zeros 0 (order r1+r2−1), −2j (r1+r2), −(2j−1) (r2), each ≤ n_K (Γ_R, Γ_C pole orders).
- A2 (EF) BV Thm 10 (l.704–743): Σ_N Λ_K(N) g(log N) = G(1) − Σ_ρ G(ρ), g ∈ C² compactly supported in (0,∞), G(s) = ∫ e^{su} g.
  Checked: BV's ĝ(−iρ) = G(ρ); A(n) log n summed over n^m = N is Λ_K(N); residues at 1 and at all zeros; contour to −∞ legitimate
  (|G(σ+it)| ≤ e^{σ u_0}‖g''‖₁/|s|², u_0 = min supp g > 0).
- A3 (ZC) Lagarias–Odlyzko 1977 Lemma 5.4 (BV l.892–897): #{ρ: |γ − t| ≤ 1} ≤ C_0(log Δ_K + n_K log(|t|+2)).

## B. Prop 9′ (replaces BV Prop 9, l.650–700; D §1)
BV Prop 9: Σ_p A(p) log p h_X(log p) = 1 + O(X² log Δ e^{−X/2}), h_X = 2e^{−X} 1_{(X−log 2, X]}; RH enters only at l.670–675
(Grenié–Molteni |ψ_K(x) − x| ≤ C√x(log x log Δ + d log² x)); the proper-power step (l.689–698) is RH-free.
Prop 9′: |ψ_K(e^X) − ψ_K(e^X/2) − e^X/2| ≤ C_2 x^{31/32}(1+X)(log Δ_K + n_K(1+X)), X ≥ 2. Audit, step by step (D §1 Steps 1–8):
- S1 kernel φ_η = I_{η/4}^{*4}: cubic B-spline, C², supp [−η/2, η/2], Laplace transform (sinh z/z)^4, z = sη/8. Correct.
- S2 |sinh z| ≤ cosh Re z (|sinh z|² = sinh²x + sin²y), sinh z/z = ∫_0^1 cosh(tz)dt ⇒ |Φ_η(s)| ≤ e^{|σ|η/2} min(1, (8/(η|s|))^4). Correct.
- S3 sandwich g_− ≤ 1_{(a,b]} ≤ g_+ with g_± = 1_{J±} * φ_η, J± = [a ∓ η/2, b ± η/2]; supp ⊂ [a−η, b+η] ⊂ (1,∞) for X ≥ 2
  (a − η ≥ 1.057); uses Λ_K ≥ 0. Correct.
- S5 |G_±(1) − x/2| ≤ ηx (0.82ηx·(1+η²/64) + (x/2)η²/64). Correct.
- S6 |G_±(ρ)| ≤ 2.6 x^{31/32}|ρ|^{−1} min(1,(8/(η|ρ|))^4); Z_j ≤ 2C_0 ℓ_K(j+1); j=0 term uses |ρ| ≥ 1/32; tail uses
  (j+3)/(J+3) ≤ j/J, Σ_{j>J} j^{−5} ≤ 1/(4J⁴), Σ j^{−4} ≤ 1/(3J³), n_K/3 ≤ 0.24 ℓ_K. Total 5.2C_0 x^{31/32}(34+log(9/η))(log Δ + n log(12/η)). Correct.
- S7 trivial zeros ≤ 0.95 n(1 + Σ e^{−j}) ≤ 1.6 n_K. Correct.
- S8 η = x^{−1/32}/4: ηx = x^{31/32}/4; logs X/32 + log 36, X/32 + log 48. Correct; C_2 absolute, uniform in K.
- Remark 1.2: Λ_{Q(√D)} = Λ(1 + χ_D) checked for split / inert / ramified p (all j). Correct.
- (9′) normalisation (D §2): prime powers p^e (e ≥ 2) in (x/2, x] number ≤ √x + (X/log 2)x^{1/3} ≤ 2X√x (X ≥ 2); Λ_K(p^e) ≤ n_K X.
  So Σ_p A_K(p) log p h_X(log p) = 1 + O(e^{−X/32}(1+X)(log Δ_K + n_K(1+X))). K = Q: w − 1 = O(X²e^{−X/32}). Correct.

## C. Prop 19′ (replaces BV Prop 19, l.1125–1137, proof l.1198–1231; D §3)
BV: Σ_O Σ_p A_{K_O}(p) log p h_X = |Ω^m/G| + O(d^m X² · 2m d^{m+1/10} e^{−X/2}) (l.1201–1204), using Lemma 21 (l.1154–1196).
- (a) n_{K_O} = |O| ≤ d^m (index of the stabiliser); ≤ d^m orbits; log Δ_{K_O} ≤ 2m d^m log(Hd) ≤ 4m d^{m+1/10} (log d ≤ d^{1/10}
  needs d ≥ 10^{16}); 4m d^{m+1/10} ≤ 4d^m X; |E_1| ≤ 9C_3 d^{2m} X² e^{−X/32} (C_3 the constant of (9′)). Correct. Lemma 21's chain re-derived: Landau
  M ≤ ‖P‖₂, Mahler |Δ_Q| ≤ n^n M^{2n−2}, (d+1)^{d−1} ≤ d^d (d ≥ 3), invariant order of primitive Q (disc = Δ_Q), Toyama by induction.
- (b) discrepancy primes (RH-free): log|Res(P̃,R)| ≤ 10mX(d log 2 + d^{1/10} + log d) + d ≤ 11mdX, from
  |Res(A,B)| ≤ 2^{deg A deg B} M(B)^{deg A} M(A)^{deg B} (‖B‖₁ ≤ 2^{deg B}M(B)); ≤ 23 m d e^{X/10} primes; |E_2| ≤ 46 m d^{m+1} X e^{−0.9X}. Correct.
  Side note: BV's Dubickas–Konyagin count of exceptional R (l.1046–1048, l.2041, ≤ e^{X/10}) is RH-free; any e^{X/5} would also do.
- (c) margin: X/160 ≥ 2m log d + 2 log X + log(18C_3); with X ≥ m d^{1/10} needs d^{1/10} ≥ 640 log d (d ≥ 10^{48.6}; recomputed
  71,620 ≤ 72,444 at 10^{48.6}); inside condiX (X > 100 m² d^{1/10}) needs d^{1/10} ≥ 6.4 log d (d ≥ 10^{25.8}; recomputed). Uniform in m.
  Only "< 1/3" is needed downstream, so these thresholds are cosmetic.

## D. Theorem 2′ (BV Thm 2, statement l.283–310, proof §8.1 l.2411–2470; D §4)
BV's RH calls in the proof: l.2433 (part 1: Prop 19 for K = Q(a)), l.2464–2466 (part 2: Prop 19 for K = Q(a_1..a_m)), l.2467 (|w−1| = O(X²e^{−X/2}) under RH for ζ_Q). No others (D's TeX census, spot-checked here).
- Window condiX (l.1980–1983): 100 m² max{κ^{−1}, τ^{−1}, d^{1/10}} < X < (κτ/C_0) d/(m log md)³, κ = τ/100.
- Part (1): m = 1, X = (√c_0/10) τ√d/log d (l.2430). Recomputed: lower end X > 10⁴/τ ⟺ dτ⁴ > (10^{10}/c_0)(log d)²; X > 100 d^{1/10}
  follows (X ≫ C^{1/4} d^{1/4}(log d)^{−1/2}); upper end ⟺ τ√d > 10 C_0 √c_0 (log d)², implied for d large. Cor cor:exp-root-2 exponent
  c_0τ²d/(100X(log d)²) = X exactly, so failure ≤ 2e^{−X/6}. Prop 19′ (X ≥ d^{1/10}, d ≥ C ≥ D_0) ⇒ ||Ω/G| − 1| ≤ 2/3 + O(e^{−X/40}) < 1.
  Rest RH-free (l.2443–2451). Union bound ≤ e^{−c′τ√d/log d}. Correct.
- Part (2): m = ⌈30(log d)²⌉, X² = c_0τ²d/(100 m (log md)²) (l.2463), X ≍ τ√d/(log d)². Lower end X > 10⁴m²/τ ⟺ dτ⁴ > (10^{10}/c_0) m⁵(log md)²
  ≍ (log d)^{12}; upper end ⟺ τ√d > 10C_0√c_0 m^{5/2}(log md)² ≍ (log d)^7, implied by τ√d ≥ C^{1/4}d^{1/4}(log d)³. Prop 19′ uniform in m.
  |w − 1| = O(X²e^{−X/32}) from (9′) with K = Q (CT covers ζ). B_m ≤ 2^{m²} ≤ e^{X/100} since X > 10⁴m². ⇒ |Ω^m/G| = B_m ⇒ m-transitive
  (lm:perm, needs deg P̃ ≥ m) ⇒ Alt by Wielandt (m ≥ 6 log d; BV cite it l.2453). Correct.
- BV's own statement "d ≥ C/τ⁴" (l.288): at dτ⁴ = C the chosen X is below the window's lower end (ratio → 0); the narrowing to
  dτ⁴ ≥ C(log d)² / C(log d)^{12} is right and is RH-independent. (Other X choices might recover parts of BV's range; not needed.)

## E. Descendants (where GRH/ERH enters in the source TeX, and the replacement)
- VX 2511.04240. RH only in Prop pr:number-distinct-factors (l.361–375; proof l.416–458, error O(d²X² log(dM) e^{−X/2}) via Prop
  pr:prime-number-theorem l.384–394 = BV Prop 9), §6 (l.972–1031: X = 5 log d, l > 5/ε, (eq:weights) PIT for Q), and Cor cr:pd via
  GS Prop 9.1 (l.1038–1050, RH for Dirichlet L). Prop pr:number-roots-random-poly (l.821–830) needs only q > d, l ≥ ε^{−1} log q/log d,
  d, q large in terms of l; VX say so (l.832–841); its proof (l.879–968) uses q only via q^{5/(l+1)} ≤ d/(10 l³ log d) and q^{1/(l+1)} ≤ d^{ε}. So X = 81 log d, l > 81/ε: errors (log d)² d^{−81/32} and d²X² log d · d^{−81/32}
  = O(d^{−0.53}(log d)³). Proper-power step (l.1004–1027) needs Σ ≥ 1/2 per realization: holds. Thm 1.1 unconditional. Correct.
- VX Cor 1.3 (statement l.220–231, proof l.1050–1086): needs Σ_{p ≤ f(d), pattern ω} log p ≤ C f(d) 2^{−π(d)} AND Thm 1.1.
  Indicator 2^{−π(d)} Σ_{m | Π r} ω_m (m/p); (m/p) = χ_D(p), D = disc Q(√m), |D| ≤ 4m, log m ≤ θ(d); ψ(x,χ_D) = ψ_{Q(√D)}(x) − ψ(x).
  Error O(x^{T}(d + log x)²) with T = 31/32 (CT) or T = B_2 + 2.2·10^{−254} (Bebel, Lean). Range f ≥ 2^{32π(d)}d^{65} resp. f ≥ 2^{aπ(d)},
  a > 1/(1 − T*) = 7.8978 (a = 7.9 leaves a(1−T*) − 1 = 2.8·10^{−4}; d_0 large but finite). Correct.
  Dependency: Cor 1.3 uses Thm 1.1, hence CT (fields Q(α), deg ≤ d, not abelian). HPAudit's map ("if HCT fails ... the Legendre
  corollary survives with HB") is wrong: without CT the Legendre corollary falls with Thm 1.1.
- Eberhard 2008.01223 (Thm 'main', l.154–171; ERH footnote l.156). RH only in Prop prime-ideal-theorem (l.604–612, = BV Prop 9) used in
  Prop local-global-general (l.614–664) / local-global (l.693–700): error O(m n^{2m+1} X² e^{−X/2} log(Hn)); PIT for Q at l.735, l.819.
  m = 1 (l.722–744): E R_φ(p) = 1 + O(1/p + p e^{−cn}) for every p > 2H+1 (no upper limit), so X = 32cn/33 balances e^{X−cn} and e^{−X/32}:
  E|Ω/G| = 1 + O(n⁵ log(Hn) e^{−cn/33}). m ≥ 2 (l.788–842): X = cn^{1/2}/m^{5/2} unchanged; (9′) extra n^{2m}X² term subsumed; m fixed.
  Thm 1.3 unconditional, same rates (e^{−cn}, e^{−cn^{1/2}}), smaller c. Correct.
- FJSS 2106.04049 (Thm 1.1 l.157–160). ERH only via Prop prop:ERH (l.481–484, = BV Prop 9) and Eb Prop 3.5 (l.486–489) at l.512–521;
  PNT with Riemann error for Q at l.514. Local input valid for 2 < p ≤ exp(c n^{1/4}); X = c n^{1/4} unchanged:
  E|Ω/G| = 1 + O(n³X² e^{−X/32} log n) = 1 + O(e^{−c n^{1/4}/32 + O(log n)}) (FJSS's printed "exp(cn^{−1/4})", l.519, l.521, is a sign typo).
  Simple spectrum (prop:simple) RH-free. Thm 1.1 unconditional, probability ≥ 1 − 2e^{−c′n^{1/4}}. Correct.
- BGS 2502.17218 (Thm 1 l.357–359; Thm 2 l.821–823; minors l.685). RH only via (eq:bv) (l.936–947, from BV Prop 9 + Lemma 21; their
  right side C k n^k log(Hn) log²x/x^{1/2} undercounts orbits, n^{2k} is right), used at l.1149 (x = exp((n/k⁸)^{1/(2C)}), k ≤ n^c),
  l.1241 (x = e^{c′n}, k = 1, 6) and l.2043. Height: P(log H ≥ e^{εn}) ≤ C_ε e^{−(b/2)εn} for every ε > 0 (cor:height l.1049–1057),
  so ε < c′/64 is allowed; error C k n^{2k}(e^{εn} + log n)(c′n)² e^{−c′n/32} ≤ e^{−c′n/64}. Thms 1, 2 and the minors corollary
  unconditional, probability ≥ 1 − Ce^{−cn}, smaller c. Correct. (eq:rh) covers only Q(a) but k = 6 needs Q(a_1..a_6): moot with CT.
- BBK 2609.40123 (Thm thm:main l.165–171 unconditional: H_2(μ) > 12, μ finitely supported, |ω_j| ≤ H ⇒ P(disc R square) ≤ C(μ)/log n;
  remark l.189: general μ "if we further assume GRH ... follows from [BV, Thm 2]"). Model monic, ω_0 own law (conditioning on ω_0 ≠ 0
  allowed: BV allow A_0, A_d other laws). Thm 2′(2) + P(R(±1) = 0) = O(n^{−1/2}) + Σ_{k≥3} P(Φ_k | R) = O(1/n) (BV l.2494–2513) +
  Gal ⊆ A_n ⇔ disc square: P(R irreducible, Gal = S_n | ω_0 ≠ 0) ≥ 1 − C′(μ)/log n for every finitely supported μ with H_2 > 12. Correct.
- BSHKP 2308.04878 (l.113): f uniform among degree-n ±1 polynomials (all n+1 coefficients ±1); "known to hold under GRH"; they prove limsup = 1.
  Our model is theirs. BV Cor 3 (l.333–338) with N = 2, U_2 = {0, ±1} (U_N = 0 and roots of unity of DEGREE < N):
  n even: f(±1) odd sums, f(0) = a_0 ≠ 0 ⇒ 1 − O(1/n); n odd: 1 − 2·C(n+1,(n+1)/2)2^{−n−1} + O(1/n) = 1 − 2√(2/(πn)) + O(1/n). Gal ⊇ A_n
  via Thm 2′(2) (τ = 0.49 fixed). Correct.
- BV Cor 3 as worded on the discovery page ("roots of unity of order ≤ N", no 0): FALSE. (i) μ uniform on {0,1}, N = 2: P(x | P) = 1/2,
  so P(irreducible) ≤ 1/2, while the page's formula gives 1 − √(2/(πd)) + O(1/d). (ii) μ = ±1, N = 3, n ≡ 3 mod 4: the omitted
  Φ_4 term P(f(i) = 0) = p(⌊n/2⌋+1) p((n+1)/2) ~ 4/(πn) ≫ n^{−3/2} (exact: ratio 40 at n = 999, 403 at n = 99999; logs/cor3_misquote.txt).
  BV's own U_N (REPORT.md and DERIVATION.md of BVIrr use it) is right; only the page's paraphrase must change.

## F. Parameter corners and reruns
- Smallest d: Theorem 2′ needs d ≥ D_0 (Prop 19′(c): 10^{48.55} with only X ≥ m d^{1/10}; 10^{25.81} with condiX; recomputed in
  logs/kernel_check.txt). Absorbed: dτ⁴ ≥ C(log d)² with C ≥ D_0 forces d ≥ D_0 (τ < 1). With BV's constants set to 1, τ = ½, the true
  margin is positive from 10^{15} (part 1) and condiX is nonempty from 10^{50} (part 2) (params_rerun.txt): "d ≥ 10^{26}" is sufficient, not needed.
- Window lower end: X > 10⁴m²/τ forces exactly the narrowed size hypotheses (C above); small τ: at τ = (C/d)^{1/4} X falls below the window.
- 31/32: Lean statement ∀ K ∀ s (31/32 < Re s → s ≠ 1 → f s ≠ 0) — no n_K, no Siegel-zero clause; a real zero in (31/32, 1) is excluded.
- Pointwise small x: only VX (X = 5 log d, primes ≈ d⁵) uses the RH bound at x polynomial in d; fixed by X = 81 log d, which VX's local
  Prop allows (q > d, l ≥ ε^{−1} log q/log d). All other uses have log x a power of d or n.
- Reruns on copies under memcap 290 MB: params.py identical to params_output.txt (peak 7 MB footprint); littlewood.gp n = 2..15 counts
  identical (peak 7 MB); kernel_check.py (Steps 2, 3, 5, 7 numerically; 80,000 samples, 0 violations).

RESULT: Breuillard–Varjú's Theorem 2 holds without the Riemann hypothesis (Cureton–Thorner's Lean-checked strip Re s > 31/32 suffices), both parts, same rates up to absolute constants: for independent coefficients bounded by e^{d^{1/10}}, the middle ones i.i.d. with any law μ with ‖μ‖₂² < 1 − τ, and dτ⁴ ≥ C(log d)^{12} (every fixed τ, d large), P = Φ·P̃ with P̃ irreducible (Φ cyclotomic times x^m, deg Φ ≤ (C/τ)√d) with probability ≥ 1 − e^{−cτ√d/log d}, and Gal(P̃) ⊇ Alt(deg P̃) with probability ≥ 1 − e^{−cτ√d/(log d)²}; consequently random ±1 (Littlewood) polynomials are irreducible with probability → 1, and the RH/ERH-conditional theorems of Varjú–Xu (its Legendre-symbol corollary in a smaller range), Bary-Soroker–Garzoni–Sodin, Blachar–Breuillard–Kozma (general μ), Eberhard and Ferber–Jain–Sah–Sawhney hold unconditionally.

# BV without RH (node N11)

Agent bvirr (Orca), 2026-11-06, 09:58–10:55 UTC. Proofs: DERIVATION.md. Search: NOVELTY.md.

## Claim
1. **Prop 9′** (DERIVATION §1–2). For every number field K and X ≥ 2:
   |ψ_K(e^X) − ψ_K(e^X/2) − e^X/2| ≤ C e^{31X/32}(1 + X)(log Δ_K + n_K(1 + X)), hence
   Σ_p A_K(p) log(p) h_X(log p) = 1 + O(e^{−X/32}(1 + X)(log Δ_K + n_K(1 + X))) — BV's Prop 9 with e^{−X/32} for e^{−X/2}.
   From Cureton–Thorner's Lean-checked strip alone (any strip Re s > 1 − δ gives x^{1−δ}).
2. **Prop 19′** (§3). Σ_p B_P(p)^m log(p) h_X(log p) = |Ω^m/G| + O(e^{−X/40}) for X ≥ m d^{1/10}, d ≥ D_0, no hypothesis.
3. **Theorem 2′** (§4). BV Thm 2 as in the RESULT line, without hypothesis on zeros: both parts hold for dτ⁴ ≥ C(log d)^{12}
   (part (1) already for dτ⁴ ≥ C(log d)²), in particular for every fixed τ and d ≥ d_0(τ), which covers all corollaries below.
4. **Corollaries** (§5–6):
   - BV Cor 3 for every μ with finite variance, not a point mass: P(P_d irreducible) = 1 − P(∃ω ∈ U_N : P_d(ω) = 0) + O_{μ,N}(d^{−N/2}).
   - Littlewood polynomials f = Σ_{i≤n} ±x^i: P(f irreducible) = 1 − O(1/n) for n even and 1 − 2√(2/(πn)) + O(1/n) for n odd;
     so lim P = 1 (BSHKP 2308.04878 proved limsup = 1); Gal(f) ⊇ A_n with probability → 1.
   - Gal(P̃) ⊇ Alt(deg P̃) for every law as in Thm 2; so Gal(P) ⊇ A_d with probability 1 − P(P has a root in U_N) − O(d^{−N/2})
     for every finite-variance law (→ 1 for the 0/1 and ±1 models; unconditionally new for both).
   - BBK 2609.40123: for every bounded μ with H_2(μ) > 12, P(R irreducible with Gal(R) = S_n | ω_0 ≠ 0) ≥ 1 − C(μ)/log n.
   - VX 2511.04240 Thm 1.1: random multiplicative ±1 coefficients: P(P_d(x)/x irreducible) ≥ 1 − C(ε)d^{−1/2+ε}. VX Cor 1.3
     (Legendre-symbol polynomials) for f(d) > 2^{7.9π(d)} (Lean-checked Dirichlet half-plane) or 2^{32π(d)}d^{65} (CT alone);
     VX's GRH range f(d) > 2^{2π(d)}d⁴ would need strip width ½.
   - BGS 2502.17218 Theorem 1 (tridiagonal: irreducible, Gal ∈ {S_n, A_n}, probability ≥ 1 − Ce^{−cn}), Theorem 2, and their
     eigenvector-minor corollary (l.685).
   - Eberhard 2008.01223 Thm 1.3 (random matrix: irreducible with probability 1 − O(e^{−cn}), Gal ≥ A_n with 1 − O(e^{−cn^{1/2}})),
     FJSS 2106.04049 Thm 1.1 (random symmetric ±1 matrix: irreducible with probability ≥ 1 − 2e^{−cn^{1/4}}), 2304.12075 Cor 1.2.
   - BV Thms 5, 6: hypotheses hold automatically for large d; Thm 7 superseded (every ζ_K is zero-free in |1 − s| < 1/32).

## Before
BV (Acta Math. 223, 2019) proved all of the above under RH for Dedekind zeta functions. Unconditionally: BKK (Invent. Math. 2023)
for laws uniform mod 4 primes (e.g. ≥ 35 consecutive integers); BSHKP (2023) limsup = 1 for ±1; Kogler 2609.26771 (Sep 2026, Lean)
the 0/1 model, P(reducible) = √(2/(πn)) + O(1/n), by twin factorizations, general laws announced. BV note (l.548–551) that the
classical zero-free region would need primes of size e^{Cd log d}, beyond the random-walk window (condiX, l.1980–1983).

## Method (BV step, TeX line, change)
1. Prop 9 (l.650–698): RH used only for |ψ_K(x) − x| ≤ C√x(log x log Δ + d log² x) (Grenié–Molteni, l.670–675). Replaced by
   Prop 9′: BV's own explicit formula (Thm 10, l.704–743) applied to a smoothed upper and lower bound of 1_{(X−log 2, X]}
   (4-fold box kernel of width η = e^{−X/32}/4), Lagarias–Odlyzko's zero count (BV l.892–897), and the fact that the strip
   plus the functional equation puts every nontrivial zero in 1/32 ≤ β ≤ 31/32 (so |ρ| ≥ 1/32). Error e^{−X/32} for e^{−X/2}.
2. Prop 19 (statement l.1125–1137, proof l.1198–1231): orbit fields K_O have degree ≤ d^m, log Δ ≤ 2m d^m log(Hd) (Lemma 21);
   error d^{2m}X² e^{−X/32} ≤ e^{−X/40} once X ≥ m d^{1/10} and d ≥ D_0 (d ≥ 10^{48.6} suffices for the margin; within BV's
   window, d ≥ 10^{25.8}). Discrepancy primes unchanged (RH-free).
3. Proof of Thm 2(1) (§8.1, l.2419–2451): same κ = τ/100, m = 1, X = τ(c_0 d/100)^{1/2}/log d (l.2430). The RH step l.2433–2441
   now uses Prop 19′: ||Ω/G| − 1| ≤ 2/3 + O(e^{−X/40}) < 1. Mahler-measure, proper-power and cyclotomic-degree steps unchanged.
4. Proof of Thm 2(2) (l.2453–2470): BV's m = ⌈30(log d)²⌉ and X ≍ τd^{1/2}/(log d)² (l.2463); m ≥ 6 log d, so Wielandt's bound
   (cited at l.2453) gives Alt; their "(30 log d)²" for Bochert–Jordan there is a harmless misquote. Prop 19′ at l.2464–2466;
   |w − 1| = O(X²e^{−X/32}) from Prop 9′ for Q at l.2467 (was RH for ζ_Q); B_m ≤ 2^{m²} ≤ e^{X/100} < e^{X/32}.
5. Cor 3, Cor 4, Thm 1 (proofs l.2494–2525): only through Thm 2(1). Thms 5–6: (log d)^α/d < 1/32.
6. Descendants: VX X = 5 log d → 81 log d and l > 5/ε → 81/ε (l.974–976; error d^{2−81/32} = d^{−0.53}); VX Cor 1.3 via
   ψ(x, χ_D) = ψ_{Q(√D)}(x) − ψ(x) and Prop 9′ (GS Prop 9.1 with x^{31/32}); BGS x = e^{c′n} (l.1241) with ε < c′/64 in
   Cor height (l.1049–1057); Eberhard X = cn/2 → 32cn/33 (l.739), m ≥ 2 unchanged (l.811); FJSS X = cn^{1/4} unchanged (l.518);
   BBK l.189 and BSHKP via Thm 2′ and Cor 3.
7. Sanity: params.py tabulates X, the margins and the smallest strip width that would suffice (δ_min ≈ 10^{−5} at d = 10^{20},
   part 1); littlewood.gp counts all monic ±1 polynomials of degree n ≤ 20 (littlewood_output.txt): the reducible ones with a root
   ±1 match the exact count, the rest stay below 2.7/n (n ≤ 20).

## Checks: every RH/GRH/ERH occurrence (rg 'GRH|ERH|Riemann hypothesis|RH' on the TeX)
- **BV 1810.13360.** Statements: l.181–184 (abstract), l.228 (Thm 1), l.293 and l.301 (Thm 2 parts 1–2), l.333 (Cor 3), l.349
  (Cor 4); zero hypotheses near 1: Thm 5 (l.378–390), Thm 6 (l.419–432), Thm 7 conclusion (l.438–449); remarks l.243, l.274,
  l.360–375, l.481–486, l.506, l.550, l.573, l.646, l.700. Uses: Prop 9 (l.659, l.670) → Prop 19 (l.1128, l.1199) → §8.1
  (l.2433, l.2464, l.2467); Cors 3–4 via Thm 2 (l.2517–2521); Prop 13/20 (mild hypotheses, l.904–1024, l.1138–1150, l.1233–1262)
  → Thms 5–7 (l.2532–2645; l.2565 and l.2623 use only unconditional facts about ζ_Q and Deuring–Heilbronn). Other analytic
  inputs, all unconditional: PNT in Lemma 26 (l.1522), PNT for w → 1 (l.2419), PNT in progressions of modulus ≤ d^{1/5}
  (l.2294, l.2314). Slips (harmless): Prop 9's error vanishes for K = Q as printed (log Δ = 0; the input has d log² x);
  l.2453 quotes Bochert–Jordan as (30 log d)², above their m = ⌈30(log d)²⌉ (l.2463); Wielandt's 6 log d, cited in the
  same sentence, covers m (misquote only).
- **BSHKP 2308.04878.** l.113 (abstract: "known to hold under GRH"), l.128 (BV). No use; their theorem is unconditional.
- **BBK 2609.40123.** l.160, l.214 (BV under GRH), l.163 (their result is GRH-free), l.189 (the general-μ remark: the only
  GRH-conditional statement). No use in proofs.
- **VX 2511.04240.** Statements l.112, l.168 (Thm 1.1), l.223 (Cor 1.3), l.210; uses: Prop 2.2 (l.384–394) in Prop 2.1
  (l.361–375, l.416–458) and §6 (l.984–1000); GS Prop 9.1 under RH for Dirichlet L (l.1038–1050) in Cor 1.3 (l.1051–1090);
  remarks l.138, l.250. All replaced (DERIVATION §6.3).
- **BGS 2502.17218.** (eq:rh) l.318–322; statements l.229, l.358 (Theorem 1), l.823 (Theorem 2), l.685 (minors); uses only
  through (eq:bv) l.936–947 at l.1149, l.1241, l.2043; mentions l.444, l.555, l.592. All replaced (§6.4). Scope slip: (eq:rh)
  covers only Q(a), but (eq:bv) for k = 6 needs Q(a_1, …, a_6). Moot now.
- **Also checked:** Eberhard 2008.01223 l.113, l.156, l.215, l.321, l.576, l.606, l.617, l.695 (uses l.733–743, l.800–830);
  FJSS 2106.04049 l.144, l.151, l.153, l.158, l.481–489 (use l.514–521); 2304.12075 l.103, l.139–149, l.237, l.653 (only via
  Eberhard and FJSS). Hokken–Koukoulopoulos 2510.18857 and 2404.13559 only mention BV.
- **Reviewer**, pass 2 (DERIVATION §4 against BV l.2411–2470): every RH use replaced and every inequality holds; it flagged
  BV's m below the quoted Bochert–Jordan threshold — a misquote: Wielandt's 6 log d, cited by BV, suffices.
- **Reviewer**, pass 3 (orchestrator's independent review of §3–4): §3(b) used 80 m d² e^{X/10} ≤ d⁴e^{X/10}, false for
  m > d²/80; repaired keeping m: |E_2| ≤ 46 m d^{m+1} X e^{−0.9X} ≤ ½e^{−X/40} uniformly in m. Rest confirmed.
- **Reviewer**, pass 1 (DERIVATION §1–3): no false step in Prop 9′ or §2; it flagged (i) the Dirichlet remark (the sandwich
  needs nonnegative coefficients) — replaced by the quadratic-field identity Λ_{Q(√D)} = Λ(1 + χ_D), which is all §6.3 uses;
  (ii) BV's discriminant and exceptional-prime bounds imported in §3 — now derived there (Mahler, Toyama, invariant order,
  Dubickas–Konyagin as stated by BV l.2041, a resultant bound).

## Novelty
NOVELTY.md: arXiv API (33 queries), arXiv listings 2026-10/11 (math.NT/PR/CO/GR), Hexagon (all 716 works), Cureton–Thorner's own
TeX, the GXL evidence API (3 claims), Semantic Scholar citations of BV, Eberhard and FJSS (13 citing works since 2025-09 read in
TeX), Zenodo, HAL, the three papers that use the strip (2610.18491, 2610.37179, 2611.06705), our files. Nothing states any of the
claims; three October 2026 papers (2610.11047, 2609.15813, 2608.27624) still call Eberhard's and FJSS's theorems ERH-conditional. Known: only the 0/1 irreducibility (BV Thm 1, Cor 4), proved unconditionally by Kogler 2609.26771 by other means. Risk
medium (Kogler's announced general-measure sequel; BV's authors helped Kogler).

## Caveat
- Rests on Cureton–Thorner's 31/32 (Lean-checked, R399; paper an unrefereed Hexagon preprint). Any fixed strip would do; the
  number 31/32 only changes constants (and VX Cor 1.3's range).
- Classical inputs cited, not re-proved: the explicit formula (BV Thm 10), Lagarias–Odlyzko's zero count, Mahler's and Toyama's
  discriminant bounds, the invariant order of a non-monic polynomial, Dubickas–Konyagin's count.
- RH-free parts of BV and of the descendants (random walks, Konyagin's estimates, Green–Tao input, SL₂ expansion, local
  problems for matrices) are used as published; we checked only that they assume nothing about zeros.
- Constants effective, not explicit: our margin needs d ≥ 10^{26} inside BV's window (10^{48.6} without it), on top of BV's
  own unspecified constants.
- Size hypothesis (RH-independent, BV's proof needs it too): BV's "d ≥ C/τ⁴" (l.288) and "τ⁴d large" (l.2430, l.2463) understate
  what their window (condiX) needs: its lower end 10⁴m²/τ (κ = τ/100) requires dτ⁴ ≥ C(log d)² in part (1) and
  dτ⁴ ≥ C(log d)^{12} in part (2). Theorem 2′ is stated with these; for fixed τ (all corollaries here) nothing changes.
- 2304.12075 Cor 1.3 (G(n, p) graphs): its transfer from FJSS's model not checked.
- Method is a substitution; the novelty is the statements. Filing suggestion: DISCOVERIES for the Littlewood limit, the
  general-law Galois statement and the unconditional matrix theorems; COROLLARIES if routine substitutions are filed there.

RH relevance: N11 (GRH families: a uniform zero-free strip removes GRH from the Breuillard–Varjú method and its descendants).

## Files
Riemann/Deep/BVIrr/: REPORT.md, DERIVATION.md, NOVELTY.md, params.py + params_output.txt (parameter bookkeeping),
littlewood.gp + littlewood_check.py + littlewood_output.txt (exhaustive Littlewood counts, n ≤ 20, checked against the exact
±1-root count).

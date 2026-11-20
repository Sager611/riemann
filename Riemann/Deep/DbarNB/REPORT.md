THEOREM: Under RH alone, lim_{λ→0} D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|² for BBLS's continuous Nyman–Beurling distance (dilations θ ∈ [λ,1]): Lemmas R″, R′ proved, reduction (3v) written out without RH, no hypothesis on ζ′(ρ); reviewers CORRECT after fixes, SIGNIFICANT (node N5)

# DbarNB — dbarnb (Orca, 20 Nov 2026, 10:02–10:55 UTC)

## Statement
- D(λ) = dist_{L²(0,∞)}(1_(0,1), span{ {θ/x} : λ ≤ θ ≤ 1 }), L = log(1/λ). Theorem: RH ⇒ L·D(λ)² → Σ_ρ m_ρ²/|ρ|²
  (= 2 + γ − log 4π if all zeros are simple). Upper bound new; lower bound is Burnol's (math/0103058, unconditional).
- Lemma R″ (RH): sup_L L·U_{Z(L⁶)}(L − √L) < ∞, hence D(λ)² ≪ 1/log(1/λ).
- Lemma R′ (RH; form needed for the limit): for each fixed H, L·U_{Z(L⁶)}(L−√L) ≤ (1+o(1))[L·U_{Z(H)}(L−√L) + ε(H)], ε(H) → 0.
- Reduction (3v), now unconditional: D(λ)² ≤ (1 + L^{−1/2})·U_{Z(T₂)}(L − √L) + O(L^{−3/2}), T₂ ∈ [L⁶, L⁶+4].
- Where RH is used: only in R″/R′, twice: (a) all ordinates real, so the dual exponential sums have real frequencies
  (Turán–Nazarov with no growth factor); (b) S(t) ≪ log t/log log t, so a cluster at scale L^{−3/4} at height γ has
  M ≪ log γ/log log γ zeros and the Turán constants e^{O(M)} cost only |γ|^{o(1)}, summable against 1/|γ|².

## Proof (DERIVATION.md, §1–7; one paragraph each)
1. Dual form (exact): U_Z(ℓ) = sup_{p ∈ span{u^k e^{iγu}}} |p(0)|²/N_ℓ(p), N_ℓ = ∫_0^ℓ(|p′|² + ¼|p|²) + ½|p(0)|² + ½|p(ℓ)|² (OU RKHS).
   So R″ is a weighted Nikolskii inequality for exponential sums with frequencies at the zeros.
2. Grouping: components of "gap < s₀ = L^{−3/4}"; RH bounds each group by M ≪ log γ/log log γ.
3. Local: Lemma A (complex extension of a group's envelope to radius 1/W, W = max(width, KM²/ℓ), via a multiplier and
   Turán–Nazarov tails; Cauchy estimates); Lemma B: |p_g(0)|² ≤ e^{O(M)}(ℓγ_g²)⁻¹∫|p_g′|²χ (γ² gain by reverse Bernstein).
4. Almost orthogonality (Lemma D): Σ_g‖p_g′‖²_χ ≤ (1 + e^{−cL^{1/8}})‖p′‖²_χ, from Lemma Λ: a group's χ-windowed spectrum leaks
   only e^{o(log²L) − c(dℓ)^{1/2}} at distance d (Newton–Hermite expansion + Gevrey-2 decay), and gaps are ≥ s₀ (s₀ℓ ≍ L^{1/4}).
5. R″ by Cauchy–Schwarz over groups: L·U ≤ (1+o(1))Σ_g e^{O(M_g)}/γ_g² < ∞. R′: fixed low part vs high groups, cross terms
   in the OU norm O(L^{−3/4+o(1)}) (Lemma E); with Claim F (explore59) the Corollary lim L·U_{Z(L⁶)}(L−√L) = Σ m²/|ρ|².
6. Lemma C = Titchmarsh Thm 9.7 only (Cartan curves dropped): rectangles between good heights tile the zeros above T₂ with no
   omission or double count; on the line ζR_n = (2πi)⁻¹∮H(w)B(s−w)ζ[s,w]/ζ(w)dw (divided difference, the referee's repair),
   so ‖ζR_high‖² ≤ e^{−cL²}; R_high is a limit of rectangle partial sums, the same for four tilings (two widths × two height
   families; finite contour sums, no bound on 1/ζ′(ρ)), which give the growth of G; (3v) follows, without RH (complex nodes below T₂).
7. Theorem = 6(v) + Corollary + Burnol's liminf.

## Evidence (numerics are falsification only; logs/ in this folder; double precision, peaks via memcap)
- Dual form (dbarnb_falsify.py dual; logs/dual.log): the OU-RKHS Gram ⟨e^{iγ_k u}, e^{iγ_j u}⟩_N equals explore59's kernel matrix
  K(γ_j, γ_k) to 1e-16 (80 zeros + two triple clusters), so L·U agrees to 10 digits: 0.0403698458 (L = 5), 0.0382192552 (20), 0.0377548398 (80).
- RH-allowed worst-case clusters (orth/clusters modes; first 300 zeros ±γ, every 8th zero above n = 20 replaced by an m-cluster of
  spacing h/L, i.e. unresolved at every L): ratio U_all/Σ_g U_g (cost vs decoupled sum) for L = 10, 20, 40, 80, 160:
  m=4, h=0.3: 1.0242, 1.0119, 1.00002, 0.99999, 0.99982 (L·U_all ≈ 0.0634–0.0650, bounded);
  m=5, h=0.3: 1.0636, 1.0144, 1.0005, 1.0012, 1.0001;  m=6, h=1.0: 1.1035, 1.0210, 1.0010, 1.0010, 0.9998.
  m=6, h=0.1: 1.057, 1.008, 0.988, 1.0002, 0.998 (roundoff-limited: single-cluster Gram condition ~1e14; not used).
- Block almost-orthogonality in the OU norm, sup_p Σ_g N(p_g)/N(p) (generalized eigenvalue): m=4,h=0.3: 524, 4.63, 1.43, 1.16, 1.08;
  m=6,h=1.0: (not PD at L=10), 1145, 1.92, 1.32, 1.13 → 1, as Lemma D predicts for the χ-windowed norm; the cost ratio is ≈ 1 even
  where this block constant is large (L = 10, 20): the proof bounds the cost, not a Riesz constant (orchestrator's R508 point (3)).
- Orchestrator's R508 checks (40 digits): single coalescent jet is the worst m-cluster, C/m² ≤ 1 (m ≤ 6); pairs decouple.
- Memory: peaks ≤ 120 MB RSS (/usr/bin/time -l), memcap ≤ 100 MB.

## Checks
- Turán–Nazarov statement (absolute c, exponent = #terms − 1, factor e^{|I|max|Re λ|} = 1 for real frequencies) as quoted in
  Friedland–Yomdin arXiv:1107.0039 Thm 1.1 (scout). Erdélyi arXiv:1602.02315 abstract: |f(0)| ≤ (1+ε_n)3n‖f e^{−9nt/2}‖_{L²[0,1]}
  for n real frequencies — a polynomial local constant (not needed; consistent with R508).
- RH inputs: S(t) ≪ log t/log log t (Titchmarsh 14.13) and real ordinates. Unconditional inputs: Titchmarsh 9.7, convexity bounds,
  functional equation, Gevrey-2 Fourier decay, Paley–Wiener–Schwartz, Burnol's liminf, explore59 (3i), (3ii), §7, Claim F.


## Review (three reviewers, one per decisive block, formulas inline; gpt-6-astra, xhigh)
- Block 1, interpolant and cost bound (§1 dual form, §3 grouping, Lemmas A, B, Theorem R″ given Lemma D): CORRECT, ROUTINE.
  Wording points adopted: jet-target vector in the duality; minimum attained (finite-measure representer ν_p); explicit counting
  in §3(iii) (reflect, start at the smallest |γ|, t₀ off the ordinates); charge each group to its largest-|γ| ordinate.
- Block 2, tail domination (Lemmas Λ, D, E, Theorem R′, Corollary): GAP, minor — Lemma D needs windows with gaps ≥ s₀, but §3
  gave enlarged containing windows; fixed (windows = convex hulls, midpoint, half-length). All other steps checked (Newton–Hermite
  coefficients incl. jets, Gevrey prefactors, cells, cross terms, Lemma E, R′, Corollary). ROUTINE given the inputs.
- Block 3, final assembly (§6, §7): GAP, fixed — (O2) second tiling omitted zeros between T₂ and t′_{n₀}, and regrouping across
  tilings needs finite contour sums (no absolute convergence over single zeros); (O3) ζ[s,w] has the pole of ζ at 1 (identity for
  s ≠ 1, R_high(1) via I_n(1)); (O4) maximum principle only for the entire G on the trivial-zero circles; (O5) Z(T₂) = complex nodes
  for the unconditional reduction. Checked: rectangles isolate each high zero once, residue signs, integrability on the line,
  growth, (v), the theorem with Burnol's liminf; RH used only in the interpolation asymptotics. Verdict SIGNIFICANT for the
  RH-only sharp continuous asymptotic; contour repairs routine.
- Second pass on the fixes (O1–O6 incl. block 1 wording): all CLOSED, no new error; verdict CORRECT for revised §3, Lemma D, §6.

## Novelty (adversarial; scouts' reports in .scratch/dbarnb/*_scout.md)
- arXiv listings this and last month (lst.py 'Nyman|Beurling|B.ez-Duarte|Burnol|Balazard', math.NT, math.CA) and Zenodo
  "Nyman-Beurling": no RH-alone upper bound for D(λ) or d_N of order 1/log. GXL (2 calls): 0 supports; hits are Burnol's lower
  bound (math/0103058), BCF 1211.5191 (RH + moment condition on 1/ζ′(ρ), discrete d_N), Dimitrov–Oliveira 1608.07887
  ((log log n)^{5/2+ε}/√log n, quoting Balazard–de Roton), and a structural preprint (preprints.org 202506.0772, abstract only).
- Balazard–de Roton arXiv:0812.1689 introduction (TeX): under RH only d_N² ≪ (log log N)^{−2/3} (Báez-Duarte) and their
  (log log N)^{5/2+δ}(log N)^{−1/2}; Burnol's d_N² ≥ (C+o(1))/log N; nothing on the continuous D(λ).
- Repo: TREE/NOGO list the limit and the 1/log order as explore59/R506 IDEA with open Lemmas R′, R″.
- Burnol (TeX of math/0103058, 0112254, 0105120, 0203120, 0208121, 0407443): only the liminf lower bound; 0112254 calls the
  matching lim "conjecture" ("The authors of [Bal3] conjecture that equality holds (also with lim in place of liminf) when one
  counts the zeros with their multiplicities"); no upper bound, no O(1/log(1/λ)) anywhere.
- Not read (not on arXiv): BBLS 2000 (Adv. Math. 149) and 2005 (Ramanujan J.), Balazard–Saias "Notes IV" (Adv. Math. 2004);
  citing sentences found (BdR, BCF, Dimitrov–Oliveira, Burnol) attribute to them only lower bounds, the conjecture and the
  discrete RH bounds above.

## Caveat
- Ineffective: L₀ is astronomically large (Turán constants e^{O(M)} against L^{3/4}); no rate for L·D² − Σ m²/|ρ|².
- Black boxes: Nazarov's lemma (absolute constant), Titchmarsh 9.7 and 14.13, Burnol's liminf; explore59's accepted parts
  (3i) window, (3ii) trivial zeros, §7 closure, Claim F (referee: correct). Identification of Burnol's D(λ) with the
  unconstrained distance used here is explore59's (REPAIR (3i), accepted); for Nyman's constrained variant the upper bound holds
  too (subtract G(1)·1, REPAIR (3i)), the lower bound then follows since constrained ≥ unconstrained.
- R′ is proved for fixed H (enough for the limit), not uniformly in H ≤ L⁶ as the brief's form asks.
- RH-blind mechanism (explore59 §4): the same proof gives the analogous limit for any function whose zeros are real in z with
  S-type cluster bounds; it does not touch d_N (explore38's integer-dilation penalty) or RH itself.

RH relevance: N5 (continuous Nyman–Beurling rate: exact constant under RH alone, replacing (log L)^{5/2+ε}/√L; no move toward RH).

## Files (Riemann/Deep/DbarNB/)
DERIVATION.md (§1 dual form, §2 inputs, §3 grouping, §4 Lemmas A, B, Λ, §5 Lemma D, R″, Lemma E, R′, Corollary, §6 Lemma C
and (3v), §7 theorem); dbarnb_falsify.py (modes dual, clusters, orth); logs/ (dual, orth_m4_h0.3, orth_m6_h1.0, clusters_m5_h0.3,
clusters_m6_h0.1 and .err with memcap peaks); dbar_*.py (explore59's scripts, copied, unchanged); REPORT.md.


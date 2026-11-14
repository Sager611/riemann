IDEA: Rudnick–Sarnak triple correlations (vanishing bispectrum on |ξ₁| + |ξ₂| + |ξ₃| < 2) carry multiplicity information that pair correlation lacks: in the exact CUE(N) analogue (Diaconis–Shahshahani moments as constraints) they cut the maximal Σm²/N by ≈ 0.022 at every N = 3–10, a pivot sum-of-squares 3-point SDP keeps ≈ 0.018 of it (pointwise 3-point LPs lose it like 1/N), pointing to ≈ 0.70 simple zeros under RH against the pair method's 0.6792 (node N8)

# ExploreSchurWeyl: Diaconis–Shahshahani exact moments → N8 simple zeros (explore50, 14 Nov 2026)

## 1. Transfer
- Tool: Diaconis–Shahshahani/Diaconis–Evans exact moments. For U Haar in U(N), p_k = Tr U^k, partitions λ, μ with |λ| = |μ| ≤ N:
  E p_λ conj(p_μ) = δ_{λμ} z_λ, z_λ = Π_j j^{a_j} a_j! (character orthogonality, Schur–Weyl). The all-index n-level correlation sums of the
  eigenvalue multiset are exactly products of traces, Σ_{j_1..j_n} e^{i Σ k_i θ_{j_i}} = Π Tr U^{k_i}; so "Haar moments for Σ_{k_i>0} k_i ≤ N"
  is the exact finite-N form of the Rudnick–Sarnak range Σ|ξ_i| < 2 (mean spacing 1, ξ = k/N), and of Hughes–Rudnick mock-Gaussianity.
- Node N8, quantity: the multiplicity sum N*(T) = Σ_{0<γ≤T} m_ρ = Σ_ρ m_ρ² (distinct ρ) and the simple share N_s/N ≥ 2 − N*/N.
  Under RH, pair correlation only: N* ≤ 4/3 N (Montgomery), 1.3275 (Montgomery–Taylor, optimal for f̂ ⊂ [−1, 1], CCLM), 1.3208 N
  (Chirre–Gonçalves–de Laat 1810.08843 Thm 1, Cohn–Elkies class f ≥ 0, f̂ ≤ 0 off [−1, 1]), so N_s ≥ 0.6792 N. Best under RH by other
  means: 19/27 = 0.7037 (Bui–Heath-Brown 2013, mollifier). CGdL's TeX has no triple/n-level/k-point input (rg 'Rudnick|triple|n-level|k-point': 0).
- Question: does the RS-range n-level information (in particular E p_a p_b conj(p_{a+b}) = 0, a + b ≤ N: the vanishing bispectrum)
  lower the maximal N*/N compatible with the data?

## 2. Finite-N test problem (exact; Schur–Weyl moments as constraints)
Random rotation- and reflection-invariant multisets of N points on Z_M (z = e^{2πix/M}). LP over the law (one variable per dihedral orbit):
maximise E[Σ_ρ m_ρ²]/N (also min E[#simple]/N, min E[#distinct]/N) subject to the Haar moments for a constraint set:
P = pair only (λ = μ = (k), k ≤ N: E|p_k|² = k); P3 = P + triples (E p_a p_b p̄_{a+b} = 0, a + b ≤ N); L4 = all with ℓ(λ) + ℓ(μ) ≤ 4;
DS = all |λ| = |μ| ≤ N. Also "relax" = the finite-N CGdL-type pair-measure LP (ν ≥ 0 on Z_M, ν̂ ≥ 0, ν̂(k) = k for k ≤ N), relax_lp.py.
Scripts: cue_lp.py, relax_lp.py; logs/.

## 3. Results (max E[Σm²]/N; simple share min = 2 − this, all extremals have m ∈ {1, 2})
Lattice LP (cue_lp.py, exact enumeration of dihedral orbits on Z_M; logs/cue_lp_small.log, cue_lp_N5.log):
| N | M | P | P3 | L4 | DS |
|---|---|---|---|---|---|
| 2 | 32 | 1.098477 | = | = | = |
| 3 | 96 | 1.151941 | 1.131801 | = P3 | = P3 |
| 4 | 48 | 1.182863 | 1.162829 | 1.155959 | = L4 |
| 5 | 40 | 1.204795 | 1.180307 | 1.173856 | 1.171684 |
Continuous positions (colgen.py: column generation, pricing by batch Adam over positions for each multiplicity pattern; the master
value is a feasible law, so a lower bound; master + max reduced cost found is a heuristic upper bound; logs/colgen_N*_*.log):
| N | relax (relax_lp.py, M = 32N) | P | P3 | L4 |
|---|---|---|---|---|
| 2 | 1.099579 | 1.1 (exact: q ≤ 1/10, below) | = | = |
| 3 | 1.153470 | 1.152294 | 1.132180 | (= P3 on lattice) |
| 4 | 1.185818 | 1.184376 | 1.164423 | 1.158253 |
| 5 | 1.207820 | – | 1.185970 | ≥ 1.179616 (ub ≈ 1.179682, time cap) |
| 6 | 1.223696 | – | 1.201701 | – |
| 8, 12, 20, 40 | 1.245011, 1.268036, 1.288035, 1.303676 | | | |
- Smallest case N = 2, exact: the constraints reduce to E cos Δ = −1/2, E cos 2Δ = 0 for the gap Δ; 3 × 3 Toeplitz positivity gives
  P(double) ≤ 1/10, so max E Σm²/N = 11/10 for P and for DS alike (the 4th-moment constraints are implied). No triple effect at N = 2.
- N ≥ 3: the triple (bispectrum) constraints cut the maximum by relax − P3 = 0.0213, 0.0214, 0.0218, 0.0220 at N = 3, 4, 5, 6 (almost
  constant in N), 4-level constraints by a further 0.006; pair-only realisable P sits only 0.0012–0.0014 below the pair relaxation.
- Mechanism at N = 3: for a configuration {a, a, b} (one double point), Re p_1² p̄_2 = 3(1 + 2cos Δ)² ≥ 0, zero only at the equispaced
  gap Δ = 2π/3; so every double point pushes the bispectrum E p_1² p̄_2 up and must be paid for by simple configurations with negative
  bispectrum, which the pair constraints do not require. The extremal laws (N = 3, M = 48, logs/cue_lp_N3_support.log): pair-only mixes
  the near-rigid {0, 13, 26}/48 (weight 0.647) and {0, 12, 25} (0.126) with the double points {0, 0, 24} (0.134), {0, 0, 13} (0.093);
  with triples, clustered simple configurations {0, 9, 18}, {0, 9, 19} (0.098) enter and the double-point weight falls from 0.227 to 0.197.


## 4. The transferred claim (precise)
- Finite N (exact constraints from Schur–Weyl): for a rotation-invariant random N-multiset on the circle with E|p_k|² = k (1 ≤ k ≤ N),
  E[Σ_ρ m_ρ²]/N ≤ v_P(N); adding E p_a p_b p̄_{a+b} = 0 (a + b ≤ N) lowers the maximum to v_P3(N) < v_P(N), and all DS-range moments
  to v_DS(N) ≤ v_P3(N). Numbers in §3: v_P − v_P3 ≈ 0.020 for N = 3–6, v_P3 − v_DS ≈ 0.006–0.009 (N = 4, 5). For N = 2 nothing changes.
- ζ (the IDEA, not proved): under RH, Montgomery's F(α) on |α| ≤ 1 plus Rudnick–Sarnak's triple correlation for f̂ supported in
  |ξ₁| + |ξ₂| + |ξ₃| < 2, i.e. max|ξ_j| < 1 for three frequencies summing to 0 (vanishing mixed third cumulants of linear statistics,
  the mock-Gaussian statement; Hughes–Rudnick's single-function case needs supp f̂ ⊂ [−2/3, 2/3])
  imply N*(T) ≤ (c₃ + o(1))N(T) with c₃ < c₂ := the pair-correlation LP optimum (≤ 1.3208, CGdL; our finite-N relaxation extrapolates
  to ≈ 1.320). Proof route: a 3-point certificate, i.e. f₂ (f̂₂ ⊂ [−1, 1]) and f₃ (f̂₃ in the RS region) with the pivot inequality
  Σ_k f₂(x_j − x_k) + Σ_{k,l} f₃(x_j − x_k, x_j − x_l) + c ≥ m_j at every point of every multiset; the finite-N configuration-LP duals are
  the circle versions of exactly this (polynomials in p_λ p̄_μ dominating Σ m² on all configurations). Checking it pointwise in (x, y)
  (an LP class) loses the gain as N grows; checking it as pointwise terms plus a pivot sum of squares Σ_i (c_i + Σ_k u_i(x_j − x_k))²
  (an SDP class) keeps it (§8).
- Why triples see multiplicity: a double point contributes m³ = 8 to the diagonal of the all-index triple sum against m² = 4 in the pair
  sum, while the RS-range bispectrum must vanish; e.g. for {a, a, b}, Re p₁²p̄₂ = 3(1 + 2cos Δ)² ≥ 0, so double points force a positive
  bispectrum that simple configurations must cancel at a pair-correlation cost.


## 5. Prior art and novelty check
- Alpöge–Furman 2608.13637 (TeX l. 589): for their Gram-matrix certificate, tr G^k is evaluable by the diagonal method "exactly in the
  Rudnick–Sarnak range X^k ≤ T^{2−ε}; at X ≍ T this allows only k = 1. Thus, unconditionally, higher moments add nothing"; conditionally
  HL*(4) (Hardy–Littlewood-type, bandwidth one, k ≤ 4) would give 13/18. Their bandwidth-one ceiling (Lean ceiling_law256, explicit
  256-periodic law) is 0.6818287 for pair-only certificates. They keep one bandwidth for all k; the RS-range triple data at reduced
  bandwidth (pairs on [−1, 1], triples on |ξ₁| + |ξ₂| + |ξ₃| < 2) used here are not considered there.
- CGdL 1810.08843: pair only (rg of the TeX: no Rudnick/triple/n-level/k-point). Rodgers 1502.05658: RS/mock-Gaussian moments give tail
  bounds for zero counts in boxes of height 1/log T (hence Σ m_ρ^k ≪_k N under RH), not a better simple-zero proportion.
  Lagarias–Rodgers 1905.12123 / 1907.03391: RS-range n-level data cannot exclude the Alternative Hypothesis (simple, lattice-supported
  mimics); they say nothing on multiplicities. Devine (Zenodo 22066689, 22287432), Santibáñez-Leal 24227401, Wang 2609.24167 (withdrawn):
  "three-point"/"consecutive-triple" Gram arguments built on pair correlation only (realisability), not RS triple correlations.
- Searches with no relevant hit: arXiv API abs:"triple correlation" AND zeta AND zeros; "n-level correlation" AND "simple zeros";
  Rudnick AND "simple zeros"; "higher correlation" AND "simple zeros"; correlation AND "multiplicities of zeros"; "three-point" AND zeta
  AND zeros; bispectrum AND zeta; monthly listings 2026-10/11 (7 categories, title regex); Zenodo '"triple correlation" zeta',
  '"simple zeros" "triple correlation"', '"n-level correlation" multiplicity zeros'; GXL evidence API (1 claim: 248 passages, the only
  "refute" is Alpöge–Furman's sentence above).


## 6. Scalable 3-point relaxation (the finite form of a 3-point Cohn–Elkies certificate)
- tri_relax.py / tri_sym.py: LP over the pair measure ν₂ on Z_M and the all-index triple measure ν₃(x, y) = E Σ_{j,k,l} 1[x_j − x_k = x,
  x_j − x_l = y] on Z_M², with ν₂, ν₃ ≥ 0, ν̂₂ ≥ 0, the Haar values of ν̂₂ (|k| ≤ N) and of ν̂₃(a, b) = E p_{a+b} p̄_a p̄_b on the DS range
  (|a| + |b| + |a+b|)/2 ≤ N, the marginal Σ_y ν₃(x, y) = N ν₂(x), and the multiplicity links for m ∈ {1, 2} on each coincidence line
  (ν₂ ≤ ν₃|_line ≤ 2ν₂, Σ_line = N ν₂(0), ν₃(0, 0) = 3ν₂(0) − 2N from m(m − 1)(m − 2) = 0). Its LP dual is a lattice/periodic 3-point
  certificate (a function of x, one of (x, y), line multipliers) — the object a ζ proof would need in the continuum.
- Checks: at N = 3 it equals the exact configuration LP on the same lattice (M = 24: 1.126526 = 1.126526; M = 48: 1.131282 = 1.131282),
  so the 3-point relaxation is exact there; at N = 4, M = 48 it gives 1.167122 against the exact 1.162829 (pair mode 1.184225 against
  1.182863; tri_sym.py, symmetric: 1.167103), i.e. it keeps 0.0171 of the exact 0.0200 triple gain. Coarser lattices lose more
  (N = 6, M = 36: 0.0127; N = 8, M = 32: 0.0060), and even at M = 12N the LP gain decays like 1/N (§8).


## 7. Adversarial test and method ceiling
- The LP optima are explicit adversaries: rotation-invariant N-point laws with the CUE density and exactly the CUE moments in the whole
  DS (= Rudnick–Sarnak) range, yet with double points: P(configuration has a double point) = 0.198, 0.308, 0.383, 0.406, 0.496 for P3 at
  N = 3, 4, 5, 6, 8 (colgen logs). So no argument that uses only RS-range correlations can prove more than 2 − v_DS(N) simple share at
  size N; this ceiling covers every admissible test function (all (λ, μ) in the range, both signs, full support), since it is the LP
  optimum over all laws (lattice LP exact up to HiGHS float tolerance; continuum colgen converged to ub − val < 2·10⁻⁶; its pricing is
  heuristic (random-restart Adam) and at N = 10 it offered no triple-point patterns, so v_P3(10) is a lower bound).
- The gain is not implied by pair agreement: the same LP with pair constraints only reaches v_P(N), and the pair-only realisable laws sit
  within 0.0014 of the pair relaxation (the CGdL-type bound), so the 0.020–0.022 drop is information carried by the triple (and higher)
  correlations, i.e. it beats the pair-correlation bound at every N ≥ 3 tested, which is what the brief's adversarial criterion asks.
- It does not separate ζ from the adversaries (by construction no RS-range statement can), and it is conditional on RH like Montgomery's
  and Rudnick–Sarnak's inputs (we did not re-read RS 1996 for its exact hypotheses; positivity arguments need real ordinates anyway).


## 8. Results at scale
Exact finite-N optimum (colgen, continuous positions) against the pair relaxation:
| N | relax (pair LP, M = 32N) | v_P3 (exact, colgen) | gain relax − v_P3 | P(double point in a configuration) |
|---|---|---|---|---|
| 3 | 1.153470 | 1.132180 | 0.02129 | 0.198 |
| 4 | 1.185818 | 1.164423 | 0.02140 | 0.308 |
| 5 | 1.207820 | 1.185970 | 0.02185 | 0.383 |
| 6 | 1.223696 | 1.201701 | 0.02200 | 0.406 |
| 8 | 1.245011 | 1.222849 | 0.02216 | 0.496 |
| 10 | 1.258515 | 1.236538 | 0.02198 | 0.581 |
The gain is flat at ≈ 0.022 for N = 3–10. relax(N) → ≈ 1.3199 (fit R − a/(N + b) to N = 12, 20, 40), consistent with CGdL's 1.3208.
3-point relaxations at M = 12N (LP, tri_sym.py) and M = 8N (SDP, tri_sdp.py, pivot positivity [[N, ν₂ᵀ], [ν₂, ν₃]] ⪰ 0, from
E Σ_j (c + Σ_k u(x_j − x_k))² ≥ 0), gain = pair mode − triple mode at the same (N, M):
| N | LP gain (M = 12N) | SDP gain (M = 8N) | exact gain on the same lattice |
|---|---|---|---|
| 3 | exact (M = 24, 48) | exact (M = 24) | – |
| 4 | 0.01712 (M = 48), 0.01726 (M = 64) | 0.02085 (M = 32); 0.02010 (M = 48) | 0.02192 (M = 32), 0.02003 (M = 48) |
| 6 | 0.01229 | 0.02144; 0.02052 (M = 36) | 0.02282 (M = 36: P 1.210697, P3 1.187874) |
| 8 | 0.00777 | 0.01749 | – |
| 10 | 0.00548 | 0.01727 | – |
| 12 | 0.00452 | 0.01821 | – |
| 14 | – | 0.01954 (tri 1.251312, pair 1.270852) | – |
| 16 | 0.00359 | – (stopped at the time budget) | – |
- Same-lattice check N = 6, M = 36 (cue_lp.py, 63,133 orbits): the SDP keeps 0.0205 of the exact 0.0228 (90%); its triple value
  1.191163 is 0.0033 above the exact 1.187874. Allowing any multiplicities (only m(m − 1)(m − 2) ≥ 0, ν₃|_line ≥ ν₂, no m ≤ 2) leaves
  every SDP value unchanged (N = 4, M = 48: 1.163123 / 1.183224; N = 8, M = 64: 1.220750 / 1.238241), so m ≤ 2 is not an assumption.
- The LP relaxation (pointwise ν₃ ≥ 0 only) keeps only ≈ 0.056/N of the triple gain: a pointwise-positive 3-point Cohn–Elkies
  certificate extrapolates to nothing in the limit (float data N ≤ 16, not a proof). The pivot SDP gain does not decay (0.0175–0.0195 for N = 8–14, ≈ 80–90% of the exact 0.022). So a ζ
  proof must use sum-of-squares (pivot-PSD) 3-point certificates, not pointwise-positive kernels.


## 9. Verdict
- New bound, rate or barrier for ζ: none proved. What is shown (float LP/SDP, finite N, the exact CUE analogue of RH + Montgomery +
  Rudnick–Sarnak): (a) RS-range triple correlations carry multiplicity information that pair correlation lacks — the maximal Σm²/N
  compatible with the data drops by ≈ 0.022 at every N = 3–10 (4-level data: ≈ 0.006 more, all n: ≈ 0.002 more at N = 5);
  (b) barrier for the obvious certificate class: pointwise-positive 3-point kernels (the LP relaxation, i.e. Cohn–Elkies with ν₃ ≥ 0
  and the multiplicity lines) keep only ≈ 0.056/N of it, so they give nothing in the limit; (c) pivot sum-of-squares 3-point
  certificates (SDP) keep ≈ 0.0175–0.0195 for N = 8–14, without decay.
- What it would give: if (c) survives N → ∞, RH + Montgomery + RS triples give N* ≤ (1.3199 − 0.017 + o(1))N, i.e. ≈ 0.697 simple
  zeros (0.702 with the full exact gain, ≈ 0.70–0.71 with all RS levels), against 0.6792 for the pair method (CGdL) and 19/27 = 0.7037
  (Bui–Heath-Brown, RH, mollifiers). So a method gain of ≈ 2–3 points, at best level with the RH record; not a record by itself.
- Conditional on RH (as both inputs are). Not RH-relevant beyond N8: multiplicity bounds do not constrain real parts (NOGO R18).
- Next step (Orca-sized): the continuum pivot-SOS 3-point bound for ζ — choose f₂ (f̂₂ ⊂ [−1, 1]), f₃ (f̂₃ ⊂ RS simplex), a PSD kernel
  Q(x, y) and line multipliers with f₃ + line terms − Q ≥ 0 pointwise and the pair conditions; first reproduce the ≈ 0.017 gain in the
  scaled limit (N = 20–40 with a lower-memory SDP), then certify.

## Files
- cue_lp.py (exact lattice configuration LP), colgen.py (continuous column generation), relax_lp.py (pair relaxation),
  tri_relax.py / tri_sym.py (3-point LP), tri_sdp.py (3-point pivot SDP, cvxpy + Clarabel/SCS), ds_mc.py (Haar Monte Carlo check of the
  DS moments: ≤ 2.7 standard errors inside |λ| ≤ N for N = 3, 4; clear failures at |λ| = N + 1, e.g. E|p_4|² = 3.00 for N = 3);
  queue1.sh, queue2.sh. Logs in logs/ (memcap peaks ≤ 128 MB for every finished run; runs above 290 MB were killed and not used).


## Suggested NOGO line
- explore50 (Orca, exploration: Schur–Weyl duality / Diaconis–Shahshahani exact moments, Hughes–Rudnick mock-Gaussianity, Rudnick–Sarnak n-level → N8) IDEA: in the exact CUE(N) analogue (laws of N-multisets with the Haar moments E p_λ p̄_μ = δz_λ for |λ| = |μ| ≤ N, the RS range) the triple moments E p_a p_b p̄_{a+b} = 0 cut max E Σm²/N by 0.0213–0.0222 below the pair relaxation at N = 3–10 (4-level ≈ 0.006 more); pointwise 3-point LP certificates keep only ≈ 0.056/N of it, pivot-SOS 3-point SDPs 0.017–0.020 (N = 8–14, no decay); if this persists, RH + Montgomery + RS triples give ≈ 0.70 simple zeros (pair method 0.6792, Bui–Heath-Brown 19/27); nothing proved for ζ; Alpöge–Furman's "higher moments add nothing" is single-bandwidth [ExploreSchurWeyl](Riemann/Deep/ExploreSchurWeyl/REPORT.md)

RH relevance: N8 (proportion of simple zeros under RH; a method gain for the pair-correlation approach, RH itself untouched).


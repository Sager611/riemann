# grhharv5 working notes (20 Nov 2026)

Harvest: axq.py brief query → 40 ids (< 100, no split); 21 NEW + R504's 4 A-candidates downloaded (.scratch/grhharv5, dl.sh). Other SEEN ids (R504's B/C/criteria, 2511.23381, 2509.24119, 2508.15715, 2508.14888, 2508.11125) skipped. 2512.23302, 2512.01588, 2510.05501 appear only in Harvest4/5 HARVEST.md notes (not REPORT.md): re-triaged.

## 2510.14167 (ω*(n) maximal order; SPD.tex)
- GRH use: only π(x;d,1) = Li(x)/φ(d) + O(√x log x) ≫ x/(φ(d) log x) for d ∈ 𝒟, d = x^{1/2−ε+o(1)} (l.203–205, 248–252); no other (l.296 Harman for unconditional).
- Paper's own general principle (l.372–383): if π(x;d,1) ≫ x/(φ(d)log x) for d | k up to x^α then ω*(n) ≥ exp((max_t f_α(t) + o(1)) log n/log log n), f_α(t) = (t log t − (t−α)log(t−α) − α log α)/(t+1−α) (l.320). Max: f_α* = −log(1 − τ), τ ∈ (0,1) with 1 − τ = τ^α (derivation: d/dτ of αH(τ)/(α + (1−α)τ), τ = α/t; check α = 1/2 → log φ; α = 0.4736 → 0.4669 as the paper).
- θ₁₄₁ pointwise: ψ(x;d,1) = x/φ(d) + O(x^θ log²(dx)) ⇒ α = 1 − θ₁₄₁ − ε = 0.12665…: f* ≈ 0.2107 (to compute) < paper's unconditional 0.6736 log 2 = 0.4669 (Harman Prop l.298–312, smooth d < x^{0.4736}). Density-exponent route: worst A(σ) on [1/2, θ₁₄₁] is 12/5 at σ = 3/4 (Huxley), range x^{5/12} < 0.4736: half-plane does not help. → A-candidate FAILS (no gain).

## 2507.05862 (norm-Euclidean cyclic cubic fields)
- Thm T:big (l.267–275) is UNCONDITIONAL ("Our paper completely solves the problem", l.267; abstract l.198): GRH only in McGown's prior result (l.264–266, 372–376 overview). R504's A-candidate moot → class 0.

## 2510.02309 (effective Brauer–Siegel for Artin L; effective-bounds-artin-arxiv-v2.tex)
- Prop grh-bound (l.262–271): GRH for ζ_K ⇒ |κ(χ)| ≪ (loglog eD_k)^{χ̃(1)−χ(1)}(loglog q(χ))^{χ(1)−⟨χ,1⟩}, ≫ (loglog eD_k)^{μ(χ̃)−μ(χ)}(loglog q(χ))^{μ(χ)−⟨χ,1⟩}; ≪,≫ depend on [k:Q], |G|, χ(1). Unconditional: Thm cond (l.131–145), log in place of loglog, times ε(χ) (exceptional factors).
- Proof (l.1010–1048) uses GRH ONLY via LO Chebotarev π_C(x) = |C|/|G| Li(x) + O(x^{1/2} log(D_K x)), x ≥ (log D_K)² (l.1012–1015) → (pnt-grh) l.1016–1019 → Mertens (mertens-grh) l.1021–1024 (y ≥ log D_k; primes in (log D_k, (log D_k)²] discarded trivially) → tail product l.1027–1030 with {1 + O(y^{−1/2} log(D_K y))} → short Euler product l.1031–1034 (T ≥ (log D_K)²) → prop:residue/prop:small-primes with T(ψ) ≥ log q(ψ).
- Replacement: FSW generaldelta δ = 1/32 for K/k (needs zero-free Re s > 31/32 for ζ_K only; valid all x > 2): error c₁(|C|/|G|) x^{31/32}(32 log D_K + n_K log x). Then tail-product error ∫_y^∞ t^{−33/32} log(D_K t) dt ≪ y^{−1/32} log(D_K y): small for y ≥ (C log D_K)^{32}. Mertens: PIT for k with x^{31/32} error for N𝔭 > (C log D_k)^{32}; (log D_k, (C log D_k)^{32}] contributes ≤ n_k Σ 1/p ≤ n_k log 33 + o(1) = O_{[k:Q]}(1) (as l.1024 and lem:mertens l.659–668 do for their ranges). Same trivial discard (Euler factors in [(1+1/N𝔭)^{−ψ(1)}, (1−1/N𝔭)^{−ψ(1)}]) brings T(ψ) back to log q(ψ). No exceptional zero: β ≤ 31/32 < 1 − 1/(4 log D_K) once D_K > e^8. Subfields K_ψ, k: half-plane holds for every number field.
- ⇒ Prop grh-bound holds unconditionally given 31/32 (constants change by factors depending on χ(1), [k:Q]). R504's 'only from x ≥ D_K^c' was wrong: FSW's error is power-saving from x ≥ (C log D_K)^{32}.

## 2511.13659 (SIVP for module lattices)
- ERH "for the L-function of every Hecke character of K of trivial modulus" (introduction.tex l.216, 438; self-reducibility-bulk.tex l.18, 579) enters through BDPW20 Thm 3.16 (l.592): equidistribution on X_1 = Pic_K^0 (l.539, 856: χ over all unitary Hecke characters = characters of the Arakelov class group, infinity part H/Log O_K^×, dim r1+r2−1); and |𝒫(B)| ≥ B/(2 log B) for B ≥ (12 log|Δ_K| + 8d + 28)^4 (preliminaries.tex l.984–989, BDPW20 Lemma A.3).
- Characters with nontrivial torus part are Größencharaktere of infinite order: not factors of any Dedekind zeta; no half-plane available (Harvest4 l.51 agrees; R504's 'divide ζ_H' holds only for class-group characters). Covered case: unit rank 0 (K = Q, imaginary quadratic), where Pic^0 = Cl_K and ζ_H covers all; but then n = rd = 2r is fixed (r fixed in Thm main l.216) and SIVP in fixed dimension is poly-time anyway: no gain. → class B.

## 2512.24147 (large quadratic character sums revisited; main.tex)
- Thm 1.1 (l.239–244): GRH, exp((log X)^{1/2+ε}) < x < X^{1/2}: max_{X<|d|≤2X, d∈F} Σ_{n≤|d|/x} χ_d(n) ≥ √(X/x) exp((1+o(1))√(log(√X/x) log_3(√X/x)/log_2(√X/x))).
- GRH only via Lemma 2.2 (l.265–276, DM Lemma 1): Σ_{|d|≤X, d∈F} χ_d(n) = main·1_{n=□} + O(X^{1/2+ε} f(n_0) g(n_1)), f = exp((log n_0)^{1−ε}), g = Σ_{d|n_1} μ²(d)/d^{1/2+ε}; used in M_1, M_2 (l.316–329) with resonator N = X^{1/2−δ}/x (l.310), resonator primes ≤ (log N)^{1+o(1)} (l.293, BT).
- Half-plane: on Re s ≥ θ+δ, |log L(s,ψ)| ≪ (log q(|t|+2))^{1−η} (Borel–Carathéodory + three circles, as GRHHarv4 A3) ⇒ error X^{θ+δ+ε} exp(C(log n_0)^{1−η}) g_θ(n_1); log n_0 ≤ (1+o(1))P_+(n) ≤ (log N)^{1+o(1)} ⇒ X^{o(1)}. Then N = X^{1−θ−δ}/x, x < X^{1−θ−δ}: error/main X^{θ−1+ε}Nx² = X^{−δ+ε}x as l.329. Candidate A: same theorem with √X/x → X^{1−θ₁₄₁}/x, range x < X^{1−θ₁₄₁−δ}.

## Dong et al. cluster (all GRH only via DM Lemma 1 = Darbar–Maiti Math. Ann. 2025): 2509.07651 (DZ; abstract misspells 'Riemnn', query missed it), 2508.16967 (abstract says 'conditional' only), 2509.20192, 2510.09005, 2512.24147
- Polynomial resonators (support ≤ y = X^{1/2−δ}·(...)): 2509.07651 Thms 1.1–1.2, 2508.16967 Thm 1.1, 2509.20192 Thm 1.2, 2510.09005 Thm 1.2. Half-plane: y ≤ X^{1−θ−δ}, constants √(1/2) → √(1−θ₁₄₁) = 0.355885…; but the unconditional GS Lemma 4.1 (error X^{1/2}n^{1/4} log n, cited by the papers, e.g. 2509.20192 l.268) supports y ≤ X^{1/3−δ}/x (short sums: M_2 error X^{1/2}y^{3/2}x^{5/2}Σr² vs main X x Σr²) or X^{1/6−δ}x^{−5/6} (Pólya-expanded sums: Σ_{|k|≤z}|a_k|k^{1/4} ≪ z^{1/4}, z ≍ (Xx)^{1/2}): constants √(1/3), √(1/6) > √(1−θ₁₄₁). No gain.
- GCD resonators (BT extremal sets, |M| = N, primes ≤ (log N)^{1+o(1)}, elements super-polynomial): 2509.07651 Thm 1.3, 2509.20192 Thm 1.3, 2512.24147 Thm 1.1 (log_3 factor). GS n^{1/4} fatal; DM's exp((log n_0)^{1−ε}) ≤ exp(P_+(n)^{1−ε}) = X^{o(1)} essential. Half-plane lemma DM_θ: error X^{θ+δ} exp(C(log n_0)^{1−η}) g_θ(n_1) → N = X^{1−θ−δ}/x. TRANSFERS.
- Erratum (2512.24147 l.322–329): M_2 error written X^{1/2+ε}N²x² (Σ_{k,ℓ≤x} 1); true Σ_{|k|≤z}|a_k| ≪ log z gives X^{1/2+ε}N²(log z)²; as written the ratio X^{−δ+ε}x beats the main term x^{−1}exp(…) only for x ≤ X^{δ/2}; l.363–370 'x/N', 'x exp' are typos for x^{−1}/N·…, x^{−1}exp. Harmless after the repair.

## 2510.14167 constant (prachar.py, mpmath): α = 1/2: 0.481211825059603 = log φ ✓; α = 0.4736: 0.466944945316581 (paper 0.4669, /log2 = 0.67366) ✓; α = 1 − θ₁₄₁: 0.210472675714423 (/log 2 = 0.30365) < 0.4669. FAILS (no gain).

## Reviewers (round 1, 3 parallel 'reviewer' subagents, verdicts in .scratch/grhharv5/review1.txt)
- A1 2510.02309: CORRECT, ROUTINE (tail control at (C log D_K)^{32}; no other GRH use; discard valid for moduli both ways; small D absorbed).
- A2 cluster: GAP, ROUTINE: (a) take δ, κ → 0 slowly for 1+o(1); (b) '4^{ω(n)} = X^{o(1)}' unsupported (ω(n) can be ≍ (log N)^{1+o(1)}/log_2 N): fix log G(n_1) ≤ 3Σ_{p|n_1}p^{−θ} ≪ (log 3n_1)^{1−θ} = o(log X) as log n_1 ≤ (log X)^{1+o(1)}; (c) App. 2 needs δ < κ; (d) whether a polynomial-support GCD construction + unconditional GS Lemma 4.1 recovers log_3 results: unresolved.
- A3/A4: (1) 2510.14167 CORRECT (0.6736 log 2 = 0.46690394; paper's f_θ(u) = 0.466945); (2) polynomial resonators CORRECT (no gain), ROUTINE.

## Extra novelty queries (arXiv API)
- (abs "Artin L-function(s)") AND (zero-free | quasi | Littlewood | residue), Oct 2025 – Nov 2026: only 2510.02309 itself and 2611.20592 (SEEN, Deep RH).
- abs "resonance method" AND (quadratic | fundamental discriminants), Sep 2025 – Nov 2026, and abs Darbar | (quadratic Dirichlet AND Omega): GRH or L-value papers only (2608.15773, 2607.25403, 2607.20408, 2603.21256, 2609.26786, 2511.14012); no unconditional log_3 character-sum result.

## Raw reviewer replies
3/3 succeeded

=== Task 1: reviewer ===
1. OK — The tail is \(O(y^{-1/32}\log(D_Ky))\), hence \(<1/2\) at \(y=(C\log(eD_K))^{32}\) for sufficiently large \(C\). This controls Euler-product tails; it does **not** give relative Chebotarev asymptotics at that cutoff.
2. OK — Under the supplied proof audit, no further GRH use remains: enlarge cutoffs, apply conductor–discriminant, then discard intermediate primes. Small discriminants/conductors require the usual regularization or absorption into constants.
3. OK — Interpret the Euler-factor inequalities as bounds on **moduli**. The discarded interval has bounded reciprocal-prime mass, so its product modulus is bounded above and away from zero; both bounds survive.
4. OK — **ROUTINE** transfer given Q1 and H: fixed positive zero-free width gives polylogarithmic cutoffs, whose logarithms preserve the claimed loglog exponents. Q1’s own significance is separate.
VERDICT: CORRECT, ROUTINE

=== Task 2: reviewer ===
1. **OK** — DM_theta follows from the stated zero-free hypothesis: three-circles gives the subpower logarithmic conductor bound; Perron absorbs the height dependence using part of δ. Isolate the square/principal contribution.

2. **GAP** — The \((\log z)^2\) correction and error exponents work. Two missing justifications: recover coefficient \(1+o(1)\) by taking δ, κ→0 sufficiently slowly; replace the unsupported \(4^{\omega(n)}=X^{o(1)}\) inference with \(\log G(n_1)\ll(\log(3n_1))^{1-\theta}=o(\log X)\). Application 2 also requires δ<κ.

3. **GAP** — The supplied \(\exp((\log N)^{1+o(1)})\) support bound neither guarantees polynomial support nor rules it out. Thus it does not justify using Granville–Soundararajan, but also does not establish that an unconditional polynomial-support construction cannot recover these results. That comparison remains unresolved from the inline evidence.

4. **OK** — **ROUTINE relative to Q2**: a standard zero-free-half-plane character-sum estimate followed by resonator bookkeeping. Q2’s own significance is outside this review; independent novelty of the applications is not established.

VERDICT: GAP, ROUTINE

=== Task 3: reviewer ===
(1) OK — stationarity gives \(1-\tau=\tau^\alpha\), hence \(C=-\log(1-\tau)\); \(C(1-\theta)\approx0.210473<0.4669\).
Numerical correction: \(0.6736\log2=0.46690394\), not \(0.466945\); comparison unaffected.
Since \(3/4<\theta\), the stated \(12/5\) density bottleneck survives. This addresses that density estimate—not every possible density method.

(2) OK — error/main ratios are \(X^{-1/2}y^{3/2}x^{3/2}\) and \(X^{-1/4}y^{3/2}x^{5/4}\), giving the claimed unconditional cutoffs.
The displayed constants require \(x=X^{o(1)}\). Otherwise retain the \(x\)-losses: the analogous half-plane Pólya estimate requires \(y\ll X^{1-\theta}/[x(\log z)^2]\), still weaker than the unconditional cutoff because \(\theta>5/6\).
These are routine comparisons of the supplied methods, not universal impossibility results.

VERDICT: (1) CORRECT, ROUTINE; (2) CORRECT, ROUTINE.
1. OK — Fixed-δ estimates permit diagonalization; choose κ = 3δ → 0 jointly, making all error exponents and implicit-constant contributions o(δ).
2. OK — The stated support bound makes both losses X^{o(1)} for fixed η; diagonalization must also accommodate η(δ).
3. GAP — δ < κ/2 alone does not control X^{δ−κ+o(1)}x throughout the stated range. Correction (5), not this inequality alone, repairs the problem.
4. GAP — A polynomial-support construction giving no log₃ gain is not an upper bound. The argument neither excludes weighted use of GS nor proves that polynomial support cannot attain the gain; the literature search does not establish nonexistence.
5. OK — ∑_{0<|k|≤z}|a_k| ≤ 4H_{⌊z⌋}. Thus the stated termwise error becomes X^{θ+δ+o(1)}N²(log(2z))²; dividing by XN gives x⁻¹X^{δ−κ+o(1)}, sufficient with the joint slow choice.

VERDICT: GAP, ROUTINE
## Reviewer round 3 (A2, narrow)
Yes—assuming the displayed estimates have uniform \(o(1)\) throughout the admissible range.
Write \(F=(\log X)^{-c}\exp((2+o(1))\sqrt{\log N\,\log_3N/\log_2N})\).
Direct relative error: \(\dfrac{xX^{\delta-\kappa+o(1)}}{xF}\le X^{-\kappa/2+o(1)}(\log X)^c=o(1)\).
Pólya relative error: \(\dfrac{x^{-1}X^{\delta-\kappa+o(1)}}{x^{-1}F}\), exactly the same bound.
Both comparisons resolve the objection: \(x\) cancels; even growth of \(F\) is unnecessary—a logarithmic loss is harmless.
VERDICT: CORRECT, ROUTINE

IDEA: build unconditional N8 certificates from the NPA level-2 moment matrix of Rudnick–Sarnak-range zero words, not from word lists: in the exact CUE(N) analogue with off-line pairs of every depth, sparse deep pairs fake the top triples for free (triples alone keep 1–8 % of their gain, N = 4–6), while the full level-2 set (triples plus their quartic diagonal partners) keeps 89–96 % (N = 4–6; N = 4, angles Z₄₈: level 2 ≥ 0.842987 for every law with off-line pairs of any depth, partly interval-checked, against 0.817137 at bandwidth one) (node N8)

# ExploreNPA — NPA hierarchy (moment matrices of words, SOS duals) → N8 simple zeros on the line (explore52, 19 Nov 2026, from 20:16 UTC)

Code: [offline_lp.py](offline_lp.py) (direct LP), [npa_cg.py](npa_cg.py) (column generation, same model; `-rays` = infinitely deep pairs),
[verify_n4.py](verify_n4.py) (50-digit primal and arb dual checks), [colgen_open.py](colgen_open.py) (explore50's continuous-position
column generation plus the open-top sets), queues [run_n6.sh](run_n6.sh), [run_inside.sh](run_inside.sh), [run_colgen.sh](run_colgen.sh); logs in [logs/](logs/).
Every run single-threaded under `.scratch/memcap.sh 290`; finished runs peaked ≤ 256 MB maxrss (≤ 230 MB by memcap); runs above the cap
(the first N = 6 L4 off-line attempt, N = 7 L4, N = 6–7 at M = 28–36) were killed and are reported as such.

## 1. Transfer (bridges measured before choosing)
| bridge | measurement | verdict |
|---|---|---|
| N1: NPA lift of Weil positivity (translation/dilation words of a test function) | translations commute, so every word moment matrix is Toeplitz in one variable: level k = Weil positivity on a larger span (Bochner/Krein), no new constraint; dilations are not symmetries of Weil's distribution, so their Gram matrices carry no state structure (= NOGO R8.9, R40, R458) | collapses at level 1 |
| N13: self-testing / rigidity of the CCM ground state | the rigidity statement is classical moment theory (Carathéodory–Fejér–Pisarenko for Toeplitz forms, Krein for the continuous case: the extension at the bottom of the spectrum is unique and atomic); simple even ground states also occur for Davenport–Heilbronn at its failing window (TREE N13) | no RH information |
| N8 small gaps (explore34/40) | the commutative NPA limit is exact local realizability, already in explore40 (float floors 0.744–0.764, certified 0.768); 7/11 DPP cap | saturated |
| N8 simple zeros on the line (zeta23/TEN, record 0.6736393) | zeta23's rank–trace inequality is NPA level 1: the 2×2 Hankel matrix of (d, tr G̃, tr G̃²) (AF 2608.13637 §7(d): Chebyshev–Markov–Stieltjes, m = 1); AF §7(e): single-bandwidth higher moments add nothing (they need X^k ≤ T^{2−ε}, and dimension d = λN caps the method at λ); bandwidth-one ceiling 0.6818287 (AF's Lean law); explore50: RS-range triples help, but only for on-line (RH-world) laws | CHOSEN |

Tool (NPA, Navascués–Pironio–Acín 0803.4290): outer approximations by moment matrices Γ_{u,v} = ⟨u*v⟩ over words of length ≤ k;
level k adds the words of length k; PSD of Γ couples odd words to even ones (|⟨u*v⟩|² ≤ ⟨u*u⟩⟨v*v⟩); dual = SOS certificate.
Transfer: words = holomorphic power sums of the zero configuration. For an off-line zero γ_ρ is complex, so a pair {ρ, 1 − ρ̄}
enters p_k with weight 2cosh(k s) instead of 2 (s ∝ (β − 1/2) log T): these are zeta23's huge Gabor entries |φ̂(x + iy)| ≤ L e^{L|y|/2}.
Level 1 = pair words E|p_k|² (bandwidth one, Montgomery/BGSTB, unconditional in the complex-zero frame); level 2 adds the
Rudnick–Sarnak-range words of length 3 and 4: E p_a p_b conj p_{a+b} (vanishing bispectrum) and E p_a p_b conj(p_c p_d), a + b = c + d
(Hughes–Rudnick mock-Gaussian moments). Load-bearing input (not re-derived here): AF §7(e) and the Rudnick–Sarnak/Hughes–Rudnick prime-side
computations evaluate every such word exactly in the RS range Σ|ξ| < 2 by the diagonal method, with no RH (zero side = holomorphic sum over all
zeros by the explicit formula).

## 2. Statement and finite test problem
Exact finite-N analogue (extends explore50's cue_lp.py to the unconditional world): laws of N-point configurations on Z_M made of on-line
points (count 1) and off-line pairs (count 2, depth s from a grid up to 5–7, weight 2cosh(ks) in p_k, p_{−k} = conj p_k by ρ ↦ 1 − ρ̄);
constraints = Diaconis–Shahshahani values of the chosen words, |λ| = |μ| ≤ N (the RS range); objective = min E[#simple on-line]/N.
Sets: P = pairs (level 1 = bandwidth one); P3 = P + triples; Q4 = P + all length-4 words with ℓ(λ) = ℓ(μ) = 2; L4 = all words with
ℓ(λ) + ℓ(μ) ≤ 4 (the level-2 data); world "on" = no pairs (explore50, RH world), "off" = pairs allowed (unconditional world).
- Claim A (barrier). Triples alone (P3) do not survive off-line zeros: sparse deep off-line pairs, with phases making |p_{a+b}| small while
  p_a p_b is large, move E p_a p_b conj p_{a+b} at vanishing pair cost; the P3 gain over P tends to 0 as the depth and M grow.
- Claim B (idea). The level-2 data (L4) survive: the quartic words price the same move (|p_a p_b conj p_{a+b}| ≤ (|p_a p_b|² + |p_{a+b}|²)/2,
  the NPA PSD inequality); off-line pairs of every depth remove only 4–11 % of the RH-world gain (N = 4–6).
- Claim C (NPA closure). A free fake needs a constrained off-diagonal word whose diagonal partner is unconstrained (the top triples in P3,
  the (3,1) words in L4); diagonal-closed sets (P, Q4, L4c) admit none, and words strictly inside the RS range are immune (N = 4–6);
  what remains are priced trade-offs (L4 −4 to −11 %, L4c −5 to −13 % of the gain).

## 3. Tests (min E[#simple on-line]/N; logs/)
a) Smallest case N = 3 (M = 96, logs/off_N3_M96.log): P 0.848059, P3 0.868199 in both worlds (explore50's 1.151941/1.131801 reproduced as 2 − max Σm²/N).
b) N = 4, M = 48, depths s ≤ 5 on a 12-point grid, plus "rays" (logs/off_N4_M48_*.log, cg_N4_M48_rays.log):
   rays = infinitely deep pairs at infinitesimal weight; they reach only the top words |λ| = |μ| = N, with directions from
   c_k = 1 + e^{−kΔ}w^k (two pairs, depth gap Δ, relative angle w; Δ = ∞: one pair), 7201 directions (npa_cg.py -rays).
| words | RH world (on) | off-line pairs, s ≤ 5 | + rays | gain over P kept off-line |
|---|---|---|---|---|
| P (level 1) | 0.817137 | 0.817137 | 0.817137 | – |
| P3 = P + triples | 0.837171 | 0.818792 | 0.817382 | 1.2 % |
| Q4 = P + length-4 words | 0.825142 | 0.825142 | 0.825142 | 100 % |
| L4 = all words ℓ(λ)+ℓ(μ) ≤ 4 (level 2) | 0.844041 | 0.842991 | 0.842988 | 96.1 % |
   Depth sweep (P3 off / L4 off): s ≤ 0.6: 0.831990 / 0.843882; ≤ 1.5: 0.819366 / 0.843149; ≤ 3: 0.818808 / 0.842999; ≤ 5: 0.818792 / 0.842991.
   Angle grid (off, s ≤ 7, no rays; cg_N4_Mdep.log), M = 24 / 48 / 96: P 0.819055 / 0.817137 / 0.816262, P3 0.824783 / 0.818792 / 0.817407,
   L4 0.853095 / 0.842991 / 0.841740: the P3 gain → 0 with finer angles, the L4 gain stays 0.0255–0.034.
c) N = 5, M = 40 (cg_N5_M40.log, cg_N5_M40_rays.log): P 0.795205; P3 on 0.819693, off 0.802105, + rays 0.797131 (kept 7.9 %);
   Q4 0.799171 both; L4 on 0.826144, off 0.824675, + rays 0.824016 (kept 93.1 %).
d) N = 6, M = 24, 8 depths + 27,937 rays incl. three-pair ones (cg_N6_M24.log): P 0.806701; P3 on 0.825733, off + rays 0.807557 (kept 4.5 %);
   L4 on 0.838170, off + rays 0.834591 (kept 88.6 %; rerun with 12,385 rays, three-pair rays on a 24-angle grid, cg_N6_M24_L4off.log,
   230 MB; the first attempt with 27,937 rays was killed by the 290 MB cap).
e) Certificates (verify_n4.py; logs/verify_*.log):
   - primal, P3 off, N = 4, M = 48, s ≤ 5: the optimal law re-solved on its 9-column support at 50 digits: value 0.818792412597, residual
     1.2·10⁻⁵³, all weights ≥ 0; the faking column is ONE configuration, two off-line pairs of depth 5 at angles 0 and 7/48, weight 1.9·10⁻¹⁸.
   - dual, L4 off + rays, N = 4, M = 48: y from the LP; slack g = ns/N − y₀ − Σ y_r f_r ≥ −9·10⁻¹⁶ on all on-line configurations (mpmath
     recheck); on every "one pair + two on-line points" configuration g ≥ −10⁻⁹ for all depths s ∈ [0, 12] (h(u) = u⁸g, u = e^{−s}, arb
     bisection with the depth-uniform acceptance lo(h) ≥ −10⁻⁹u⁸, 4.3·10⁶ boxes), and for s → ∞ g ~ c e^{8s} with c ≥ 0.0676 > 0 for every
     shape (the remainder for s > 12 is bounded by hand, not interval-checked); on two-pair configurations g ≥ −9·10⁻¹⁶ (float scan of
     (s₁, s₂) ∈ [0, 15]², plus the rays for s → ∞; not interval-checked). Hence, up to these two hand/float steps, every law on Z₄₈ with
     off-line pairs of any depth matching the level-2 words has E[#simple on-line]/N ≥ 0.8429878883, while level 1 admits 0.817137
     (explicit on-line law): a level-2 gain of 0.02585 in the unconditional analogue.
   - the same dual for L4 WITHOUT rays fails at depth 15 (two pairs, angle 9/48: g = −6·10⁴⁶ = −2.3·10⁻⁶ relative): finite depth grids
     alone are not enough; the rays are what make the certificate depth-uniform.
   - The same NPA reading holds for every N tested: with rays, the off-line optimum for P3 equals the optimum for P plus the triples strictly
     inside the RS range (N = 6: P3 off + rays 0.807557 vs P3o 0.807541); the rays remove exactly the top (|λ| = N) triples.

f) Where the gain sits (logs/cg_opentop.log, cg_inside_N6N7.log). "o" sets keep words of length ≥ 3 only strictly inside the RS range
   (|λ| ≤ N − 1; all pair words kept). These words are immune to off-line pairs (on = off to 6 digits at N = 4, 5, with every depth and rays):
| N, M | P | P3o | P3 (on) | L4o | L4 (on) | inside share of the level-2 gain (L4o − P)/(L4 − P) |
|---|---|---|---|---|---|---|
| 4, 48 | 0.817137 | 0.817275 | 0.837171 | 0.817885 | 0.844041 | 2.8 % |
| 5, 40 | 0.795205 | 0.796848 | 0.819693 | 0.800269 | 0.826144 | 16.4 % |
| 6, 24 | 0.806701 | 0.807541 | 0.825733 | 0.814086 | 0.838170 | 23.5 % |
   The top level |λ| = N is Σ|ξ| = 2, the boundary of the Rudnick–Sarnak range (explore50's dictionary ξ = k/N): almost all of the small-N
   gain sits there, and there the triples are fakeable while the quartics protect them. The inside share grows with N.


g) More checks (logs/cg_N6_M24_opentop_off.log, cg_L4c.log, cg_inside_N6N7.log, colgen_N8_*.log):
   - N = 6, M = 24, off-line world (8 depths ≤ 5, clusters of up to three pairs): P3o 0.807541 and L4o 0.814086, equal to the RH world.
   - L4c = level 2 closed under diagonal partners (L4 without the (3,1) words): N = 4 on 0.841363, off + rays 0.840035; N = 5 on 0.822794,
     off + rays 0.819318: no free fake, but a priced loss of 5.5 % and 12.6 % of the gain (L4 loses less: its (3,1) words still constrain).
   - N = 7, M = 21 (on-line): P 0.824515, P3o 0.825753, P3 0.849848, L4o 0.835332 (L4 killed by the cap): inside share of the triple gain
     4.9 % (N = 5–7: 6.7, 4.4, 4.9 %), while the inside level-2 gain grows (L4o − P: 0.0008, 0.0051, 0.0074, 0.0108 for N = 4–7).
   - N = 8, continuous positions (colgen_open.py = explore50's column generation plus the open-top sets; colgen_N8_P.log, colgen_N8_P3o.log;
     max E[Σm²]/N, converged to reduced cost < 10⁻⁷): P 1.242943, P3o 1.240929, P3 1.222849 (explore50): the inside share of the triple
     gain is 10.0 %; the triple information still sits mostly on the top level.

## 4. Why (NPA mechanism; proofs in the finite analogue)
- Rays fake top triples (Claim A). Two pairs of equal depth s → ∞ with relative angle w, weight ε: they add ε e^{2Ns}v to the top words,
  v = Re(Π_{k∈λ} c_k Π_{k∈μ} conj c_k)/√(z_λ z_μ), c_k = 1 + w^k, and o(1) to all lower words. For a top triple (a, b; a + b = N) the pair cost is
  ε e^{2Ns}|c_N|²/N and the shift ε e^{2Ns}Re(c_a c_b conj c_N)/√(z_{(a,b)}N); letting w^N → −1 (c_N → 0, c_a c_b ≠ 0) the ratio
  cost/shift is O(|c_N|) → 0. This is
  cglunc's Prop O mechanism (deep pairs break configuration-wise certificates; CGLUnc) moved from CGL kernels to triple words; it needs
  K ≍ N T^{−2δ} off-line zeros at depth δ, allowed by every known density slope (8/7 < 2), as in Prop O.
- Quartics price them (Claim B). The level-2 moment matrix Γ (rows/columns 1, p_k, p_a p_b) is PSD for every law and a ray adds the rank-one
  PSD block t·u u* with u = (c_N, c_a c_b): |Re(c_a c_b conj c_N)| ≤ |c_a c_b||c_N|, so a shift Δ of the triple entry raises the diagonal
  entries E|p_N|² and E|p_a p_b|² by amounts whose geometric mean is ≥ |Δ|. Constraining the quartic diagonal (Q4/L4) turns the free move
  into a priced one; Q4e (pairs + diagonal quartics only) gains nothing itself (0.817216 at N = 4): the information is in the triples, the
  quartics only defend it. Inside words (|λ| ≤ N − 1) are immune because any deep configuration that moves them by Δ moves the positive
  top pair word E|p_N|² by ≫ Δ.
- Closure (Claim C). Any configuration contributes the rank-one PSD matrix u u*, u = (p_λ)_λ, to the moment matrix; deep configurations
  move diagonal words only upward. So a move is free only through an off-diagonal word whose diagonal partner is not constrained: the top
  triples of P3 (partner (a,b)|(a,b) missing) and the (3,1) words of L4 (sextic partner missing). Diagonal-closed sets (P, Q4, L4c) have
  no free move; their off-line losses are trade-offs that the LP prices (0 % for P and Q4, 5–13 % for L4c).
- Clusters. J deep pairs give c_k = Σ_j e^{−kΔ_j}w_j^k, a recurrence of order J, which can vanish on J − 1 top levels; so in larger
  analogues free fakes could reach levels n ≥ N − J + 2 for non-closed word sets (not seen at N ≤ 6 on the lattices tried); closure
  blocks them at every level.

## 5. ζ dictionary and status
- Off-line zero ρ = 1/2 + δ + iγ: a factor of frequency ξ ∈ (0, 1] gets weight T^{ξδ} (zeta23's |φ̂(x + iy)| ≤ L e^{L|y|/2}); a word of total
  weight σ = Σξ gets ≍ T^{2σδ}. Sparse deep zeros therefore act only near the top σ → 1 of the RS range, like the rays; by Conrey's slope
  8/7 words of weight σ < 4/7 are untouched by off-line zeros (K T^{2σδ}/N ≤ T^{(2σ − 8/7)δ}).
- Transferred claim (conditional on the load-bearing input of §1, not proved here): an unconditional N8 certificate built from RS-range
  word statistics of ζ's zeros (holomorphic, complex γ) must carry the quartic words next to every triple word it uses near the top of
  the range; triples-only certificates (explore50's P3 route) are defeated by sparse deep off-line zeros, exactly as CGL kernels are (Prop O).
- What is NOT shown: (i) the continuum limit N → ∞ (the analogue's top level is ζ's RS boundary Σ|ξ| = 2, reachable only as a limit of
  inside words, where faking costs are finite, not zero; at N ≤ 7 most of the gain sits on that level, the inside share of the level-2
  gain growing from 3 % to ≈ 24 % for N = 4–6); (ii) any bound for ζ: no certificate in zeta23's Gabor-matrix frame was built (it needs a
  configuration-wise degree-4 matrix inequality, an Orca-size job); (iii) the unconditional evaluation of RS-range triple/quartic words in
  zeta23's frame (AF §7(e) asserts it for tr G̃^k; the mixed-frequency words are the same diagonal computation, not re-derived here).
- Comparison with the node: record 0.6736393 (arena TEN), TEN ceiling 153/227, AF's bandwidth-one ceiling 0.6818287, which is this
  analogue's P problem at N = 256 (their explicit 256-periodic law with near-CUE form factor). P falls with N (0.848, 0.817, 0.795 at
  N = 3–5; 0.757057 at N = 8, continuous positions, logs/colgen_N8_P.log) towards it; in the RH world triples add 0.019–0.025 (N = 4–8)
  and level 2 adds 0.027–0.032 (N = 4–6); 89–96 % of the level-2 gain survives off-line pairs at N = 4–6. Level 2 is the first input class beyond bandwidth one that the analogue
  says can survive off-line zeros; explore50's triples-only route needs RH (Claim A).
- Next step (if pursued): the continuum version of the off-line LP (stationary laws, Montgomery form factor on [0, 1], RS-range triple and
  quartic form factors on Σ|ξ| ≤ 2 − ε, off-line pairs with weights T^{ξδ}) to see whether a level-2 gain survives ε → 0 and N → ∞.

## 6. Adversarial
- Model check: the certified level-2 bound (≥ 0.8429873 at N = 4) is met by CUE(4) itself (all zeros simple, share 1); the extremal laws
  have m ∈ {1, 2} plus, off-line, only rays (infinitesimal weight): the adversary never uses off-line zeros at positive density.
- Davenport–Heilbronn: the inputs are the RS-range word values, which come from ζ's Euler product through the diagonal method (prime powers
  only, Λ weights); DH's explicit formula has no such prime side and its words are not the CUE/GUE values, so nothing here applies to DH or
  to an arbitrary positive-definite distribution (level 1 alone holds for every configuration, |p_k|² ≥ 0, and is the part that gains nothing).
- The barrier half (Claim A) is about methods, not ζ: it says what a triples-only certificate cannot do unconditionally, and it is the same
  sparse-deep-pair mechanism as cglunc's Prop O (CGLUnc), here for words of length 3.
- Discretisation: lattice angles (M = 24–96) and a depth grid; rays remove the depth cut; continuous angles only via the M sweep
  (N = 4) and the continuous-position column generation of §3(g).

## 7. Novelty
- arXiv math.NT listings Oct and Nov 2026 (show=2000; title scan for simple zeros, n-level, triple correlation, pair correlation, moment
  matrix, semidefinite): nothing on higher correlations for simple-zero proportions; 2609.02882 (Lamzouri, 67.25 % simple on the line)
  and 2501.14545 (Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh, box hypothesis) are pair-only; arXiv search for NPA/moment-matrix
  hierarchies applied to zeta zeros: none (0803.4290 itself only).
- Our files: explore50 (ExploreSchurWeyl: triples/L4 gains, on-line laws only, "nothing proved for ζ"); AF §7(d)–(f) (rank–trace as
  Chebyshev–Markov–Stieltjes, single-bandwidth "higher moments add nothing", HL*(4) ⇒ 13/18); cglunc Prop O (deep pairs kill CGL kernels);
  Yang–Yang (Zenodo 21975237, GUE moments, REFUTATIONS). The off-line extension, the rays, the quartic protection and the boundary
  concentration are not in any of them.
- Zenodo (API queries: "triple correlation" zeros zeta simple; "n-level" "simple zeros"; "Rudnick-Sarnak" "simple zeros"; moment matrix
  off-line; higher correlations critical line proportion; NPA hierarchy Riemann zeta): 24227401 (a sharp ratio inside Wang's withdrawn
  three-point kernel framework, conditional), 23489033 (67.33 %, ten-zero LP windows) and unrelated physics items; nothing with RS-range
  words, off-line laws or moment-matrix closure.

## Suggested NOGO line
- explore52 (Orca, exploration: quantum information / NPA hierarchy (moment matrices of operator words, SOS duals, Tsirelson, self-testing) → N8, N1, N13) IDEA: N1 lifts collapse at level 1 (commuting translations give Weil positivity again), N13 self-testing is Krein/Pisarenko (no RH information); N8: in the exact CUE(N) analogue with off-line pairs (the unconditional world) sparse deep pairs fake the top Rudnick–Sarnak triples for free, so triples alone (explore50's route) keep 1–8 % of their gain without RH (certified primal, N = 4), while the level-2 moment matrix (triples with their quartic diagonal partners) keeps 89–96 % (N = 4–6; N = 4: ≥ 0.842987 for every law and depth, partly interval-checked, vs 0.817137 at bandwidth one); free fakes need a word whose diagonal partner is unconstrained (NPA closure); at small N the gain sits mostly on the RS boundary Σ|ξ| = 2; no bound for ζ [ExploreNPA](Riemann/Deep/ExploreNPA/REPORT.md)

RH relevance: N8 (unconditional proportion of simple zeros on the line: which correlation data can pass AF's bandwidth-one ceiling 0.6818 without RH).

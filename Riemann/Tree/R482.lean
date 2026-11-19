/-!
# Round 482 — node N9 (moments, Lindelöf, subconvexity), round 16, [check]/[falsify] of Fodge's explicit bound |ζ(1/2+it)| ≤ 0.3080532 t^{1/6} log t + 2.405897 t^{1/6} + 506.6 for t ≥ 2.9·10^12 (Zenodo 24325329): his explicit B-process (Lemma C.3) survives 26,500 exact-sum tests, Lemma C.2's closed forms are exact, the bookkeeping and constants reproduce; an external record above 2.9·10^12, ours below; free probability, Borodin–Okounkov and Grothendieck's inequality transfer no tool

* Node: N9 (16th round; R474 before). Round type: [check]/[falsify] of a fresh certified Zenodo claim that beats our node record (lessons 147, 201),
  from the 14–19 Nov gap scan: arXiv monthly listings math.NT/CV/SP above 2611.17255 (206 entries, 72 keyword hits, none moving a node record;
  nearest N9: 2611.24820 sixth/twelfth moments in prime-power moduli, 2611.21551 complex moments of ζ′ under Lindelöf, 2611.19383 Ω± for S_n(t)
  under RH); Zenodo from 14 Nov with paging (60 records); HAL since 8 Nov (Oukil v23, Bonelli); QRT PRs 142/143 merged; arena leader Bakken/Bei Yi.
* Claim: Fodge, 'Exact Stechkin weights and improved zero-free regions for the Riemann zeta-function and Dirichlet L-functions'
  (Zenodo 24325329, 17 Nov 2026), Prop. C.1: |ζ(1/2+it)| ≤ 0.3080532 t^{1/6} log t + 2.405897 t^{1/6} + 506.6 for t ≥ 2.9·10^12, by Revers's
  blocks (2602.05614) of length ⌈N/80⌉ with an explicit B-process (Lemma C.3) in place of the second-derivative estimate.
* [falsify] Lemma C.3 (|Σ_{n=0}^{L′−1−m} e(g(n))| ≤ µ′L(m/W)^{1/2} + 4.524(W/m)^{1/2} + 1.0004 + 2.24 m/W for N ≥ 10^4, L ≤ N/75):
  22,000 random and 4,521 adversarial cases (stationary points within 10⁻⁴ of an integer or half-integer at an endpoint, m/W ∈ [0.5, 2],
  L′ = L, m = 1): max ratio 0.5469; the stationary-phase regime (N = 10^4, every m, 5000 t per m): max 0.3589; 4.524 = 8/√π + (√h − 1)/2,
  h ≤ 1.0407: HOLDS (scouts c3rand, c3stat).
* [check] Lemma C.2 (weighted Weyl differencing): W1 = M(M+1)(M+2)/6 and 30 r(m) = (M−m)(M−m+1)(M−m+2)(M² + 3Mm + 2M + m² + 3m + 2)
  exact for all M ≤ 300 (integers); the inequality on constant, quadratic, n log n and random sequences: max ratio 1 (equality for constant
  sequences, as Cauchy–Schwarz allows): HOLDS (orchestrator; scout c2c4 timed out).
* [check] Prop. C.1 numbers: 0.5397, 0.3787, 0.3518, 0.3321 times t^{1/6} log t at log t = 28.6957, 41.3378, 56, 100; it beats our 0.5601
  (t ≥ 10^10) from its start 2.9·10^12 and our 0.482061 (t ≥ 8.97·10^17) from 2.5·10^13; inputs HPY 2207.02366 Lemma 2.1, Revers 2602.05614
  (refereed by our revaudit), Dhiman–Kadiri–Quesada-Herrera 2609.00537 (explicit B-process, cited); no newer explicit bound below 0.3081
  (axq.py): HOLDS (c1num).
* [check] Prop. C.1 bookkeeping at t = 2.9·10^12: N0 = 13548, n1 = 679374, 315 blocks, the last cut at n1; head 463.566 < 4√0.95 t^{1/6} − 2;
  Riemann–Siegel remainder < 2; b1, b2, S as his interval log: HOLDS (c1head). The orchestrator reran appC_sp_consts.py 2.44 80 0.95
  0.3080532 2.405897 506.6 on a copy: log identical (10 s, 56 MB).
* Verdict: Prop. C.1 survives; it is the explicit N9 record for t ≥ 2.9·10^12 (unrefereed; the B-process proof, l.2021–2136, not refereed
  line by line) → TREE N9; ours stays best for 3 ≤ t < 2.9·10^12. The same paper's N7 claims (classical shape 4.6892 for ζ through Prop. C.1,
  4.9774 for Dirichlet L with 3 ≤ q ≤ 400 000) and its Dedekind companion (Zenodo 24354828) → the next N7 round.
* Exploration: free probability (R-transform additivity of zero blocks fails: union variance 320.9 against free 641.7), Borodin–Okounkov
  (RH does not make the Li Toeplitz matrix positive semidefinite: λ₀ = 0, det T₂ = −0.0080; R446 had strong Szegő), Grothendieck's inequality
  (true for every matrix, so it cannot separate RH from off-line Weil matrices): no transfer.
-/

namespace Fayed.R482
end Fayed.R482

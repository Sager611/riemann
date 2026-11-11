/-!
# Round 427 — node N4 (Robin/Lagarias), round 12: Gracia Llorente's 2-adic Robin theorem (Zenodo 24044336: v_2(n) ≤ 73 from the 7/8 half-plane) authenticated: statement, exact identity and 2-slot lemma hold, it inherits Chourasiya–Simonič's unsupported b_4 (harmless); float model reproduced (K = 75 with our range and B_2); prior art: Axler v_2(n) ≤ 20 (Thm 1.4), an unreviewed 25-free claim; biomechanics, cosmic rays and avalanches give no route

* Node: N4 Robin/Lagarias (twelfth round; R422 before). Round type: [check] a candidate's source before an Orca launch
  (robin2adic, harvest6 candidate 1).
* [check] Zenodo 24044336 (A. Gracia Llorente, 7 Nov 2026, PDF): Theorem 1 as quoted (hypotheses: zeta != 0 for Re s > 7/8 and
  cited explicit estimates; conclusion: Robin for every n > 5040 with v_2(n) <= 73); Proposition 2 is an exact identity
  log(sigma(n)/(e^gamma n log log n)) = E(x) - B(n) - D(n) - H(U) with B, D, H >= 0; Lemma 3: D_a(2) > 2^(-a-1): HOLDS.
* [check] Inputs: Broadbent-Kadiri-Lumley-Ng-Wilk Chebyshev bounds, Platt-Trudgian RH to 3e12, Bellotti-Wong zero count,
  Chourasiya-Simonic zero density (Table 1, all rows with b_4 = 0.54, e.g. sigma in [27/32, 28/32], B_1 = 21.77): our
  REFUTATIONS entry finds b_4 = 0.54 unsupported (0.5501 needed, B_1 up to 1.6% larger): FLAW in an input, harmless for his
  theorem (his E at L = 106 is 1.227e-23, more than 4 times below 2^(-74)).
* Computation: the harvest's float model (Harvest6/k_est.py) reproduced: 1.227e-23 at L = 106 (his value); 7.793e-24 at
  L = 108.58 (our Robin range 10^(10^46.79)) with beta <= B_2; 7.916e-24 with B_1 raised 1.6%; hence K = 75, K = 76 from
  L >= 110: HOLDS. This is far beyond R294's Morrill-Platt endpoint, where the half-plane was too weak; no contradiction.
* Prior art: Axler (arXiv 2110.13478) Thm 1.2: every 21-free n >= 5041; Thm 1.4: every n >= 5041 with v_2(n) <= 20 (before
  Hertlein, v_2(n) <= 19); no later unconditional valuation result (arXiv, Zenodo); zeta.teal-sea.com (robin_tfree,
  Corollary 1) claims every 25-free n unconditionally from E(x) <= 2.4807e-8, a candidate pending review (it reports that a
  25-free claim of Morrill-Platt was withdrawn). TREE N4 corrected (theorem numbers; 21-free and v_2 <= 20 are different
  families). Candidate robin2adic (v_2(n) <= 75 unconditionally) stays open; brief updated.
* Exploration: biomechanics (Hill force-velocity passivity, positive-real tissue models), cosmic-ray transport (leaky box and
  Fokker-Planck positive semigroups, power laws) and avalanche science (critical branching, depinning): no prime structure: no
  route.
-/

namespace Fayed.R427
end Fayed.R427

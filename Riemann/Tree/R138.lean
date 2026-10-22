/-!
# Round 138 — audits of six conjecture claims + second check of the Blaschke gap

Scouts (8/8 ran):
* 2610.26686 (Erdős–Szemerédi sunflower conjecture, all k: |F| ≤ 3(2 − 1/(16k³))^n): none found by two
  scouts (interpolation lemma l. 405–532, dimension count Lemma 3.7 l. 538–605, block transfer l. 612–685;
  brute force k = 3 small D, n). For k = 3 this is Naslund–Sawin 2017; for k ≥ 4 the naive slice-rank tester
  gives no saving (Σ_(s ≤ (1−2/k)n) C(n,s) ≥ 2^(n−1)), so the new step needs a deep audit → Orca.
* 2610.28828 (counterexample to Gale's conjecture in R³): verified. Exact rational check of the 101 points:
  2030 edges at distance > 0.88928·D, no proper 4-colouring (exhaustive DSATUR, 18401 branches), so every
  4-partition has a part of diameter > 0.88928·D > β₄(ball) = √((3+√3)/6)·D = 0.8880738·D.
* 2610.29584 (R(5,10;5,8) = 21): not reproducible as published: no explicit order-20 graphs or
  certificates, the UNSAT runs are only asserted, no code URL.
* 2610.24442 (strong Erdős–Frankl–Füredi cover-free conjecture, T(r) ≥ (r+1)²): none found; brute force
  r = 1, n ≤ 3 and r = 2, n ≤ 4; 20000 random weighted checks of κ(H)² ≤ Σ 1/p_x.
* 2610.25868 (Terao freeness, 15–16 lines): none found; 22,544,280 matroids, code and certificates public.
* 2610.26232 (Baranov–Bessonov–Kapustin conjecture): none found.
* 2609.23822 (Blaschke conjecture, adapted complex structures): the gap is confirmed by a second scout.
  l. 716–723 build the Stein-factor cover whose fibres are the components of 𝒦_z and only a section of it;
  a section of a finite étale cover does not force one-point fibres, so "𝒦_z is connected, hence
  𝒞_z = Z_p" (l. 735) is unsupported and the VMRT recognition step is unjustified. → reviewer.
-/

namespace Fayed.R138
end Fayed.R138

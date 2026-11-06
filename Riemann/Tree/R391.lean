/-!
# Round 391 — node N7 (zero-free regions and zero density), round 10: B_2 + 2.2e-254 stands (tracker: #119 merged 0.8749570194201, weaker; hexagonmath, Zenodo, HAL, arXiv nothing new); the Analytic Number Theory Exponent Database has no 2026 half-plane (last commit 26 Sep); Liu's B_2 recomputed exactly from the Hexagon TeX; isoperimetric rearrangement, neural networks (a Nyman-Beurling restatement) and SAT/SOS proof complexity give no route

* Node: N7 zero-free regions and zero-density constants (tenth round; R381 before).
* [check] QuasiRiemannTracker (GitHub API): #119 merged 4 Nov (0.8749570194201, weaker), #114 open (0.96875), #117 closed; nothing
  later beats Bebel's Lean-checked B_2 + 2.2e-254: HOLDS.
* [check] off arXiv: hexagonmath.org's newest item 2611.00127 (Hexagon numbering) is a math.OC paper on Hurwitz interval polynomials
  (seen R348); Zenodo ('"zero-free region"', '"zero density"', '"quasi-Riemann"', '"zero-free half-plane"') 23935840, 23963781,
  23971353, 23933662: no new region or density bound; HAL: none: HOLDS.
* [check] arXiv since 28 Oct: 2611.06416, 2611.06964, 2611.05657, 2611.01771, 2610.40510 (statistical-mechanics zero-freeness) improve
  none of B_2 + 2.2e-254, R_0 = 4.896, C = 53.989: HOLDS.
* [check] the Analytic Number Theory Exponent Database (Tao-Trudgian-Yang, github.com/teorth/expdb; first check here): last commits
  26 Sep 2026; zero_density.tex (l. 276, 774-776, 975) has the Guth-Maynard / Ivic / 2025 bounds and zero_free_region.tex (l. 84-137)
  only classical and Vinogradov-Korobov regions; no 2026 half-plane or better exponent than our N7 line: HOLDS.
* Adversarial: Liu's contact constant from the Hexagon 2611.00082 v1 TeX (eq. vtwo:B2, l. 86-89; contact factorization l. 3901-3916):
  1278 B^2 - 2331 B + 1061 = 0 has discriminant 9729 = 9 * 1081, smaller root (777 - sqrt 1081)/852 = 0.87338196661320751993,
  matching Hexagon2's audit: HOLDS.
* Exploration: isoperimetric and rearrangement inequalities (they control positive local kernels; Weil's form is signed, and positivity
  on a cone of symmetric decreasing functions would already be RH-equivalent; R352), neural networks (2309.09171 recasts
  Nyman-Beurling approximants as one-layer networks; the zero-free region it gives is the classical Cauchy-Schwarz bound
  |rho| >= (2 Re rho - 1)^{1/2}/delta, which with delta = 0.073 at N = 5120 excludes only |rho| < 14: KNOWN, N5), SAT/SMT and SOS
  proof complexity (finite Boolean or PSD certificates, no analytic tail control; R345): no route.
-/

namespace Fayed.R391
end Fayed.R391

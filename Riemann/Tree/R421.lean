/-!
# Round 421 — node N8 (critical-line proportion and gaps), round 11: Lamzouri's unconditional small-gap constant 0.7892 (arXiv 2611.14486; was 1 - 2.36e-6) authenticated and recorded with the large-gap values; Yin-Zhao-Hao-Sun's RH proof (2611.12438) refuted; candidate: optimise Lamzouri's minorant and inputs; population genetics, array processing and kinetic theory give no route

* Node: N8 critical-line proportion and gaps (eleventh round; R418 before). Round type: [check] new papers from the listings after
  the outage; gaps between zeros were not on the N8 line yet.
* [check] arXiv 2611.14486 (Y. Lamzouri, 10 Nov 2026), Thm 1.1 (TeX l.195-197): D(lambda) > 0 for every lambda >= 0.7892,
  unconditionally (a positive proportion of gaps at most lambda times the mean spacing); before Zhao 1 - 2.36e-6 (arXiv 2603.17334,
  l.190-192) and Simonic-Trudgian-Turnage-Butterbaugh 1 - 10^{-3e13}; mu <= 1/2 (l.184); no large-gap theorem: RECORD (theirs),
  added to the N8 line.
* [check] Inputs: asymptotic N(T), Conrey's density (1989), Fujii's moments of S(t+h) - S(t), the BGST pair-correlation mean
  O_f((log T)^{-1/2})N(T); no 2025-2026 unrefereed input; the final inequality at lambda = 0.7892, theta = 4/7 - 1e-5 recomputed by
  quadrature: 1.00010014764 > 1: HOLDS.
* [check] Large gaps: best unconditional Lambda > 2.64 (Bui-Hall-Subira-Jorge 2511.14415, Thm gapthm, P = 1), 3.18 under RH (1410.3635);
  no 2025-2026 paper or Zenodo item beats them (our gaplambda agent's target is not pre-empted): HOLDS.
* Harvest: 2611.14486's remark (l.199) 'no serious attempt to obtain the best constant': the Selberg minorant of [-1, 1] was chosen for
  easy verification, a degree-162 Gaussian polynomial already gives below 0.7864; mu <= 1/2 and lambda below 0.7892 stay open:
  candidate (Orca): optimise the minorant and substitute the 2026 inputs (simple zeros on the line).
* Refutation: arXiv 2611.12438 (Yin, Zhao, Hao, Sun, 9 Nov 2026) claims RH (l.157-160); its Prop. 3.1 (l.1218-1232) infers
  |zeta(a+ib)| != |zeta(1-a-ib)| from |chi(a+ib)| != 1, which fails when both vanish, as at every zero: REFUTATIONS.
* Exploration: population genetics (Ewens, Poisson-Dirichlet: value distribution only, already R370), array signal processing
  (Pisarenko, root-MUSIC: finite Toeplitz certificates, already R413) and kinetic theory (H-theorem, Bobylev: no zero mechanism,
  already R389): no route.
-/

namespace Fayed.R421
end Fayed.R421

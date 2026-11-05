/-!
# Round 347 — node N2 (de Bruijn-Newman), round 6: records stand; Hakimoglu-Brown's 11/12 table obsolete under the Lean-checked half-plane; density alone cannot beat 2a^2

* Node: N2 de Bruijn-Newman constant (sixth round; R337 before).
* [check] Zenodo and arXiv since R337: Lambda <= 0.107559925318 (ours, DBNB2) stands; Raymond's Zenodo 23752103 is still v1.0
  (0.107964910079); Planat-Raymond 2611.06245 (0.156821, classical strip) and Hakimoglu-Brown 2610.12803 (0.13963 from the
  7/8 half-plane) are weaker; Zenodo 23964801 concerns Davenport-Heilbronn: HOLDS.
* Adversarial: Zenodo 23931764 proves only Lambda_chi <= 9/32 for real chi (from 7/8; ours is 0.278828 with B); newer deposits
  are surveys or say they are not proofs: HOLDS.
* Harvest: Hakimoglu-Brown's uncertified sigma = 11/12 comparison table (TeX l.141-147; about 79 rows, 300 rational entries) would
  bound Lambda from a zero-free half-plane weaker than the Lean-checked Re s > 0.87338196661320752, so its bound cannot beat
  the record: not worth computing.
* Wild: zero-density estimates cannot push the strip bound below 2a^2 = 0.27883 without RH to a height: one pair of zeros at the
  edge (z = +-2ai) already needs heat time 2a^2 (de Bruijn 1950 Thm 13), and N(sigma, T) bounds allow one zero: KNOWN.
* Exploration: Stein's method (no Euler-product mechanism beyond Berry-Esseen CLTs for log zeta), positive geometry (canonical
  forms are signed sums; Weil's prime terms oppose the archimedean term, 2610.27903), quantum computing (2303.04602 samples
  xi zeros with a six-qubit Fourier transform; no mechanism, classical computation is far ahead): no idea.
-/

namespace Fayed.R347
end Fayed.R347

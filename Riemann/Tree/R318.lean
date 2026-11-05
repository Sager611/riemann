/-!
# Round 318 — node N8 (critical-line proportion), third round: 67.354% stands; the paper's open ground-state question is our n8opt run; transfer to Dirichlet families untried

* Node: N8 critical-line proportion (third round).
* State of the art: 2610.08965 (at least 67.354% simple and on the line) stands. Newer papers are weaker: 2610.15329 (6.78% on
  the line), 2610.26841 (65% isolated); 2610.23467 (simple or on the line) was audited earlier. Zenodo: nothing newer.
* Open question: 2610.08965, section 'How far the method goes' (TeX l.443-447): for a fixed window the best constant is a
  ground-state problem. Our Orca agent n8opt answers it for the paper's window: the 8-point certificate is LP-saturated
  (re-optimised data <= 0.673540, proved ceiling 0.673566); more needs a new window or longer frames.
* Families: Conrey-Iwaniec-Soundararajan (more than 56% of zeros of Dirichlet L-functions on the line, on average) has not been
  revisited with 2610.08965's detector: transfer untried (Orca candidate).
* Adversarial: no new Zenodo or arXiv claim above 67.354% (KNOWN items only).
* Exploration: Kac-Rice counts on random models of the detected function give heuristics, not ceilings (weak idea);
  random-matrix mollifier ceilings do not apply because 2610.08965's detector is not a mollifier (no idea). The xi' proportion
  and Taylor-model scouts timed out.
-/

namespace Fayed.R318
end Fayed.R318

/-!
# Round 145 — counterexample hunts in freshly posed conjectures

Scouts (8 launched, 7 finished):
* FOUND: Sun's circular-permutation conjecture for primitive elements (as stated in arXiv 2610.06289, Conj. 1.1;
  Z.-W. Sun 2013, Conj. 3.8) fails for q = 11, a₀ = 2. In G_(11,2) (x ∼ y iff x ≠ y and 2 + xy is a primitive root
  mod 11, i.e. in {2, 6, 7, 8}) the vertices 2 and 9 both have neighbourhood {3, 8}, forcing the 4-cycle 2–3–9–8 in any
  Hamilton cycle on 10 vertices. Exhaustive search agrees; the only failure with q > 7 among primes q ≤ 61.
  [note](Riemann/Deep/SunPrim/NOTE.md), script Riemann/Deep/SunPrim/verify.py. Pending: Sun's original wording,
  prime powers, novelty, reviewer.
* FOUND: the open problem of arXiv 2610.22248 (l. 1351–1354; over an arbitrary field) fails over Q for the
  biregular path P3: a lex Gröbner basis of the full isomorphism system (9 matrix entries + det inverse) forces
  e³ = 4; isomorphic over C (256t⁶ − 1 in the basis). Similarly K_(1,3) (α³ = 9, diagonal normalization).
  Pending: independent check, the singular non-biregular graphs over C, novelty, reviewer.
* None: 2609.28594 power-of-two cycles with at most one low-degree vertex (85,254 graphs) and the rainbow version
  (180 clique colourings); 2610.03729 "every polyhedral graph is schurian" (64 graphs); 2610.05776 Wythoff rows
  (n ≤ 60, 2000 entries); 2610.25661 factorial separation (n ≤ 200). Necklace Hamiltonicity timed out.
-/

namespace Fayed.R145
end Fayed.R145

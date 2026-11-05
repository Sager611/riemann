/-!
# Round 250 — Sun's 2013 conjectures: three more fail in tiny cases (3.3, 3.4(ii), 3.5(i)); candidates under novelty check

* Counterexample candidates (scouts; orchestrator re-checked by brute force, .scratch/sun38/c35.py), Z.-W. Sun, arXiv 1309.1679:
  - Conj. 3.3 (distinct consecutive differences from a fixed first element when n ∤ |G|): G = Z/6, elements 0, 1, 3, 5 with
    a₁ = 0 — any ordering puts the three odd elements consecutively, and their two differences lie in {2, 4} and coincide.
  - Conj. 3.4(ii) (|G| odd, n ∤ |G|: a cyclic ordering with distinct differences): G = Z/5, A = {0, 1, 2}; a 3-term
    progression always gives cyclic differences {−d, −d, 2d}.
  - Conj. 3.5(i) (|G| ≢ 0 mod 3, n > 3: cyclic ordering with distinct a_i + 2a_{i+1}): fails for every 5-subset of Z/7 and for
    some 4- and 5-subsets of Z/8.
  Together with Conj. 3.10 (R249). Orca "sun38" checks the latest arXiv and 2019 journal statements and prior notes.
* Other Sun conjectures tested in their smallest cases hold: 3.1, 3.7, 3.8, 3.9, 3.11–3.13, 3.15–3.18, 4.1, 4.2, 4.4–4.10,
  4.12, 4.13; 4.3 and 4.11 not testable. A scout's "counterexample" to 3.14 at n = 4 was wrong: 0-2-3-1-4 has sums
  2, 5, 4, 5, 4 = (p+1)/6 for Sophie Germain p = 11, 29, 23.
-/

namespace Fayed.R250
end Fayed.R250

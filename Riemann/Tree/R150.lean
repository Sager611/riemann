/-!
# Round 150 — two more Sun counterexamples (Conj. 3.12(ii), OEIS A228728); OEIS sweeps; three Orca agents finish

* Sun, arXiv 1309.1679, Conj. 3.12(ii) (adjacent differences and products distinct, except {±s, ±t}): fails for
  {2, 3, 4, 6} (2·6 = 3·4), {−8, −2, 1, 4}, {−4, −1, 2, 8}, {−6, −4, −3, −2} and {±1/3, ±1, ±3}. 3.12(i) survives 25,776 sets.
* OEIS A228728 (Sun: "unable to find a counterexample"): fails for {0, 1, 2, 3, 4, 5, 8} (cyclic arrangement
  (0,3,4,2,8,1,5) has distinct distances; none with ends 0 and 8). Verified: Riemann/Deep/SunPrim/verify4.py.
* OEIS sweeps ("In general", "primitive root", "any field"): A228762 survives (odd |G| ≤ 11, Z/3 × Z/3); A185645(ii)
  (= Conj. 3.3) also fails in Z/3 × Z/3 (n = 4) and (Z/2)³ (n = 5). Three sweeps timed out.
* Orca smale: arXiv 2610.21472 (Smale's mean value constant ≥ √(8081/6656)) is sound (exact Schur–Cohn re-check of the
  pair (9/16, 3/8)); the earlier disproof is GPT-6 Astra's Lean proof (github.com/tadamcz/mean-value-problem), also
  non-explicit; no explicit P with S(P) > 1 (best S = (d−1)/d up to d = 48; known pairs need degree ≳ 10⁴).
* Orca alldisc: candidate DISCOVERY: A > 1 for almost all discriminants Δ (all orders, fundamental or not), hence every
  primitive form of discriminant Δ represents a prime ≤ 2|Δ|^A (pending review).
* Orca stlempfix: the STLemP gap is repaired (N = max{1, 2^⌈log₂(P/s₀)⌉}); thresholds unchanged (pending re-review).
-/

namespace Fayed.R150
end Fayed.R150

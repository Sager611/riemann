/-!
# Round 449 — node N2 (de Bruijn–Newman), round 12, [check] of Missiaen's certified Davenport–Heilbronn bounds (Zenodo 24155142 v15.2.0): Λ_DH > 28825920237/125000000000 = 0.230607361896 holds (our rerun of the full 512-bit suite passes; tail bounds, Interval-Newton step and the Φ_DH ↔ ξ_DH normalization audited) and is new; Λ_DH ≤ 2(σ* − 1/2)² = 0.769699258304 by de Bruijn's strip theorem, the strict '<' by a qualitative Ki–Kim–Lee localization; σ* = 1.1203624981833 is Bombieri–Ghosh's prime-angle root; so no argument that also applies to Davenport–Heilbronn can prove Λ ≤ 0.2306; Renegar derivatives, Sturmian words and Haemers/Lovász bounds transfer no tool

* Node: N2 (12th round; R431 before). Round type: [check] of a new unrefereed certified claim next to the node
  (Zenodo 24155142, T. Missiaen, v15.2.0, 12 Nov 2026; TeX, verify script and README in the record's zip).
* [check] Certificate (orchestrator): the full verification suite (python-flint, 320/384/512-bit balls) reran: ALL PASS in
  100.6 s, peak 28 MB; z_240(t*) = 480.690231042507472707… − 9.535301365553739·10⁻⁷ i, N(B) ⊂ int B at radius 10⁻³²,
  Im N(B) < 0, so H_{t*}^DH has a non-real zero at t* = 0.230607361896 and Λ_DH > t* (de Bruijn monotonicity).
* [check] Code (script l.433-512, TeX Lemma l.574-592): E1 (n ≥ 34 truncation of Φ_DH), E2 (u ≥ 2.3 tail; the worst
  case k = 4, |y| = 10⁻⁶ is ≤ 6·10⁻¹³⁴⁷ of the bound), M4 and the Interval-Newton step are sound; acb.integral is rigorous
  for the entire integrand (the ignored 'analytic' flag only affects speed): HOLDS.
* [check] Normalization: Φ_DH(0.15) = Φ_DH(−0.15) to 10⁻³⁵; |D_DH(ρ_240)| = 9·10⁻¹¹ against 1.6·10⁻² at ρ_240 + 0.01;
  2∫_0^{2.3} Φ_DH(u)cos(zu)du / [(5/π)^{(s+1)/2} Γ((s+1)/2) D_DH(s)] = 1 at z = 10 and 20 (s = 1/2 + iz/2): the
  heat flow is Polymath15's convention, so Λ_DH > 0.2306 > 0.10756 ≥ Λ_ζ (our N2 record): HOLDS.
* [check] σ* = 1.12036249818332506180… (script's 320-bit enclosure; an independent crude sieve enclosure
  [1.1124, 1.1335] contains it). No zero has Re ρ ≥ σ*: only the phase bound |arg((1 − ix)/(1 + ix))| = 2|Re arctan x|
  ≤ 2 arctan|x| for |x| < 1 is needed, x = p^{−s} for p ≡ 2, 3 (mod 5) (primes ≡ ±1 cancel in L(s, χ̄)/L(s, χ)), so a
  scout's 'GAP' about torus points with non-unit moduli is no gap (re-derived); equality only at the aligned torus point,
  excluded for real γ (TeX l.377-380).
* [check] Upper bound: de Bruijn 1950 Thms 12-13 apply with Δ = 2σ* − 1, so Λ_DH ≤ 0.769699258304; the strict '<'
  (Cor. l.548-560) localizes à la Ki–Kim–Lee with the finitely many low zeros handled existentially (not explicit).
* [known] σ* and the prime-angle identity are Bombieri–Ghosh's (Russian Math. Surveys 66 (2011)); Dobner 2005.05142
  proves a generalized Newman conjecture, not bounds for DH; no earlier bound on Λ_DH found: the lower bound is new.
* Consequence (node line): an argument that also applies to D_DH (functional equation, Γ-factor, periodic Dirichlet
  coefficients, no Euler product) cannot prove Λ ≤ 0.2306; our 0.10756 uses the Euler product via the zero-free half-plane.
* Exploration (tool transfers, tested): Renegar derivatives of Gårding hyperbolic polynomials (→ N6: they only descend
  from already hyperbolic Jensen polynomials), Sturmian complexity and Fine–Wilf (→ N7: zero-occupancy words carry no
  analytic content), Haemers rank and Lovász theta bounds (→ N1: Weil positivity supplies no canonical graph): no transfer.
-/

namespace Fayed.R449
end Fayed.R449

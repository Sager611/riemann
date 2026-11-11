/-!
# Round 423 — node N6 (Jensen hyperbolicity), round 11: the quoted inputs of today's explicit Kim–Lee wedge authenticated (Kim–Lee Thm 3 = Obreschkoff's sector theorem, Riemann's integral and the factor 32, DLMF 10.21.15, 10.16.9, 10.21.3, the Lean-checked half-plane) and its constants recomputed; Bui–Milinovich's large-gap 3.18 corrected to unconditional; volcanology, sediment transport and crystal growth give no route

* Node: N6 Jensen hyperbolicity (eleventh round; R413 before). Round type: authenticate the quoted inputs of a fresh discovery
  (kimeff, discoveries/jensen-kim-lee-wedge.md, Riemann/Deep/KimEff/DERIVATION.md (I1)-(I4)).
* Kim–Lee arXiv 2105.05386 Thm 3 (TeX l.168-189): S(delta) = {|Im z| <= delta |z|}; real P, Q, Z(P) in S(delta), Q hyperbolic,
  deg Q <= delta^(-2) imply P(D)Q hyperbolic; their corollary: J(P; d) hyperbolic for d <= delta^(-2), no condition P(0) != 0:
  exactly as quoted: HOLDS. Their Thm 1 (l.132-136): real entire transcendental f of order rho < 2 with zeros in a strip has
  N(f; d) = O(d^(c/2)) for every c > rho, ineffective; for xi (rho = 1) this is our quoted 'n >> d^(1/2+eps)': HOLDS.
* Riemann's integral xi(1/2 + z) = 4 int_0^oo phi(t) cosh(z t) dt, phi(t) = sum_k pi k^2 e^(5t/2) (2 pi k^2 e^(2t) - 3) exp(-pi k^2 e^(2t)),
  and gamma(m) = 32 m! Mt(m)/(2m)!: independent mpmath check (dps 30) at z = 0, 0.5, 2 and m = 0..3, relative errors <= 1.8e-31: HOLDS.
* DLMF (page text read by the orchestrator; the scout timed out): 10.21.15 the product formula for J_nu with nu >= 0; 10.16.9
  J_nu(z) = (z/2)^nu 0F1(-; nu+1; -z^2/4)/Gamma(nu+1); 10.21.3 nu <= j'_(nu,1) < y_(nu,1) < y'_(nu,1) < j_(nu,1) < ...: as quoted: HOLDS.
* Half-plane: a = 0.37338196661320752 is Bebel's Lean-checked Re rho <= 0.87338196661320752 (QuasiRiemannTracker, standard axioms)
  with rho -> 1 - conj(rho); the exact B_2 = (777 - sqrt 1081)/852 is smaller by 7.1e-20 (no visible change); the parabola step
  (Im w)^2 <= 4a^2(a^2 - Re w) re-derived with sympy (difference 4(a - b)(a + b)(a^2 + t^2) >= 0 for |b| <= a): HOLDS.
* Computation: constants recomputed independently: n_* = 998,616,377,771; D_n >= 6.45 n^2/log^2 n (3.59 for a = 1/2) at
  n = n_*, 1e13, 1e20, 1e50, 1e200; D_n > D_all from n = 4.69741725e13 (D_n/D_all = 1.000000022 at 4.6974173e13); the pairs
  (c, C) = (0.55, 0.87160) and (1, 1.1853e-10) hold with smallest ratio 1.130: HOLDS.
* Side (gaplambda's flag): Bui–Milinovich's Lambda > 3.18 (1410.3635v3 l.53, l.70) is unconditional for consecutive critical
  zeros (RH only to pass to all zeros); R421 and the N8 line said 'under RH': corrected. Under RH: lambda > 3.63846723
  (Durkan–Pearce–Crump 2610.28527).
* Exploration: volcanology (eruption-repose renewal, chamber passivity), sediment transport (Exner nonlocal bedload as a positive
  jump generator with jumps log p) and crystal growth (Wulff convexity gamma + gamma'' >= 0): each scout 'IDEA' only re-encodes
  Weil's explicit formula as the field's positivity condition; nothing in the field proves it, and the dissipation of a positive
  jump generator is positive for every jump law and has no archimedean term: no route.
-/

namespace Fayed.R423
end Fayed.R423

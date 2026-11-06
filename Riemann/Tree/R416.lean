/-!
# Round 416 — node N12 (Function-field and Selberg-class analogues), round 12: the quoted inputs of today's two N12 discoveries authenticated: Dousselin's printed constants, setting, odd-character remark and the range varkappa < 1/8 (TeX l.614-655, 710-723, 1534-1538), his Lemma 12 and eqintero giving C5' = 2/sqrt(pi varkappa (1 - theta)), Bauer's unconditional twisted second moment; the classical and GRH discriminant limits and our A_c(15/32) = 23.38789 recomputed, Odlyzko's 1976 Table 2 below our 31/32 column at n = 8, 20, 100, 1000, the kernel lemma's Fourier transform to 1e-49; gravitational lensing, cardiac excitable media and enzyme kinetics give no route

* Node: N12 Function-field and Selberg-class analogues (twelfth round; R406 before). Round type: authentication of the quoted
  inputs of today's discoveries selbergn4 (discoveries/selberg-linear-combinations.md) and odlyzko31
  (discoveries/discriminant-bounds-strip.md).
* [check] Dousselin 2311.10285: Theorems mainth/application (TeX l.614-640) print the range, the table (7.38, 4.91, 3.68, 2.94,
  1.46e-9), the asymptotic bound and the remark that odd characters are handled similarly; l.651-655 the setting (distinct
  primitive characters, real nonzero c_j); l.710-723 the hypothesis and l.1534-1538 the error exponents T^{1/2+4 varkappa},
  T^{1/3+5 varkappa} that force varkappa < 1/8, so our varkappa = 0.1249999 is admissible: HOLDS.
* [check] Lemma 12 (l.1107-1112) and eqintero (l.1640-1645): exact main term -Pi_{qd}(v) Delta_2(vL)/(2 sqrt(pi v) log xi) with
  Delta_2(X)/sqrt X -> -4, giving C5' = 2/sqrt(pi varkappa (1 - theta)); Bauer's twisted second moment (Theorem l.1191-1205) is
  unconditional for h, k <= T^omega, omega < 1/2, which covers xi = T^varkappa: HOLDS.
* [check] Discriminant limits (mpmath, two quadratures): A_c(1/2) = 22.3816160955 = 4 pi e^gamma, A_c(0) = 44.7632321910,
  A_c(15/32) = 23.3878942309; totally real 60.8395403239 and 215.3325159535; consistent with Odlyzko's survey (JTNB 2 (1990),
  doi 10.5802/jtnb.22): HOLDS.
* [check] Odlyzko's 1976 Table 2 (his web page): our 31/32 column exceeds it at n = 8, 20, 100, 1000 by 0.015, 0.065, 0.306,
  0.676 (all signatures) and 0.062, 0.344, 1.576, 3.341 (totally real); Diaz y Diaz's 1980 tables: no public source found: HOLDS.
* [check] Lemma A: the Fourier transform of cosh(ax)/cosh(wx) agrees with the closed form to 1.9e-49 at four test points and is
  positive; with a positive-definite S the zero terms Re Phi(rho) >= 0 for |Re rho - 1/2| <= w: HOLDS.
* Exploration: gravitational lensing (Khavinson-Neumann image counts need finite rational degree; no arithmetic positivity),
  cardiac excitable media (phase-singularity winding counts census all zeros and cannot single out the critical line) and enzyme
  kinetics (King-Altman positivity does not give real-rooted binding polynomials): no route.
-/

namespace Fayed.R416
end Fayed.R416

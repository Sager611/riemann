/-!
# Round 307 — node N10 (explicit prime bounds), third round: best-bound table by range (Chirre-Helfgott to 1.6*10^83, then Fiori-Jaskari); 113.67 rechecked; Mertens and pi - li records added

* Node: N10 explicit psi and pi - li bounds (third round).
* Computation (relative error eps(x) = |psi(x) - x|/x): Chirre-Helfgott (2512.15709) gives 1.05*10^-12 for
  10^40 <= x; Fiori-Jaskari (2609.23222, 7/8) wins from x* = 1.614*10^83 (log10 x* = 83.208), e.g. 8.3*10^-15 at 10^100
  and 2.6*10^-27 at 10^200; the classical Fiori-Kadiri-Swidinsky bound is far weaker throughout (0.198 at 10^83).
  Already implicit in R286/R295 and NOGO: KNOWN.
* Adversarial: Chirre-Helfgott's constant: C_T = log^2(T/2pi)/(2pi) - log(T/2pi)/(6pi) = 113.66888055 < 113.67 at
  T = 3*10^12 + 1 + pi/3 (Theorem mainthmB -> Cor. psiexpl): HOLDS.
* State of the art (Mertens): |M(x)| < sqrt(x) exhaustively for x <= 10^16 (Hurst); isolated M(10^24) = 7,189,337,839
  and M(10^25) = -258,560,632,948 (2607.07566); limsup M(x)/sqrt(x) > 1.826054, liminf < -1.837625 (Hurst 2018);
  unconditional M(x) << x log x exp(-eta (log x)^{3/5} (log log x)^{-1/5}) explicit forms in 2208.06141. Added to TREE.md N10.
* Wild (pi - li): |pi(x) - li(x)| < sqrt(x) log x/(8 pi) for 2657 <= x <= 1.101*10^26 (Johnston 2109.02249 Cor 3.3);
  Chirre-Helfgott prove psi bounds only. A two-sided bound does not give pi(x) < li(x); the verified sign range stays
  as published (Buthe 2018: 10^19); a crossing exists below 1.3982*10^316. No new claim to refute.
* Exploration: S-matrix bootstrap positivity does not transfer (off-line cross terms are sign-indefinite; Gram/SDP
  certificates only cover a chosen span) (N1, none); Lieb-Thirring/CLR counting fails at zeros (-Z''/Z has
  non-integrable (m - m^2)/(z - r)^2 singularities) (N7/N13, none); S-fraction coefficients of sum p_{k+1} z^k are
  positive for truncated zero sets automatically (finite positive measures), so they carry no information (N3/N6, none).
-/

namespace Fayed.R307
end Fayed.R307

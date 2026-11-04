/-!
# Round 214 — a second Zenodo proof of RH breaks (a non-log-concave 'PF∞' function); leads for long searches; hunts clean

* Refutation: Bonelli, Zenodo 10.5281/zenodo.23788542 (30 Oct 2026), Theorem 68.1: the claimed PF∞ function
  Λ₀(u) = e^{u/2}ϑ₃(0, e^{−πe^{2u}}) is even with Λ₀(u) ≥ e^{|u|/2}, hence not log-concave: log Λ₀(0) = 0.0829,
  log Λ₀(±1) = 0.5000 (.scratch check, mpmath). PF∞ of the moments is equivalent to RH anyway. Logged in REFUTATIONS.md.
* Lead (needs an Orca audit): openai/math "A sharp Fourier certificate for planar circle packing" (claims a sharp
  Cohn–Elkies function in dimension 2, open since 2003): F = e^{−πhs}(p₁ + K₂), s = b|x|², b = √3/2, h = 2/5, with
  Fourier-measure kernels K_i and tabulated coefficients (build/sections/fourier.tex, signs.tex, certificates.tex).
* Leads for long computer searches: 2610.28547 (BCC crystallization), 2610.29593 (filling-area bound 5.42828),
  2610.21263 (binary/spherical code gaps), 2610.39116, 2611.03763, 2611.01064 (lists in .scratch, not kept).
* No counterexample: 2611.04923 (x ≤ 10⁵, 16 cases), 2611.05060 (Erdős–Gyárfás, cyclic Cayley graphs of orders 7–11),
  Erdős problem 9 (n ≤ 15, uninformative). One scout timed out (automated graph conjectures).
-/

namespace Fayed.R214
end Fayed.R214

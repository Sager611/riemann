/-!
# Round 302 — node N13 (Hilbert-Polya), second round: no candidate gained a property; every 2024-2026 operator's missing step is equivalent to RH; sign obstruction confirmed

* Node: N13 Hilbert-Polya and spectral approaches (second round).
* State of the art: Yakaboylu (2309.00405, 2408.15135 v17): a non-symmetric half-line operator has point spectrum
  i(1/2 - rho) for simple zeros; the decisive positivity W > 0 is equivalent to RH, and only then is the derived h
  self-adjoint. Connes-Consani-Moscovici 2511.22755 Thm finmain: finite self-adjoint operators approximate low zeros,
  convergence unproved; 2310.18423 stays semilocal. Bender-Brody-Muller 2211.01899: self-adjointness missing. No
  candidate has a new proved property.
* Adversarial: a scout reports that Yakaboylu's similarity lemma (v17) uses omega(t) = t e^t/(1 + e^t)^2, so
  omega(T)^{-1} and V = omega(T)^{-2} are unbounded and not boundedly invertible on L^2(dt), hence the similarity does
  not give self-adjointness on the stated space (the first-zero eigenfunction is in L^2). Not verified by us and the
  paper makes no unconditional RH claim: recorded as a possible gap only. The Zenodo/arXiv claim scout timed out.
* Wild: the trace-formula sign obstruction is known (Connes 9811068 absorption spectrum; Berry-Keating/Sierra
  0805.4079; Bolte-Keating 0912.3183): positive orbit amplitudes cannot give Weil's negative prime terms; 2511.18309
  (spectral-shift jumps) is formal beyond finite truncations. Inverse spectral theory: Wu-Sprung-type potentials match
  finitely many zeros (0811.1389); no unconditional infinite construction.
* Exploration: Lapidus-Maier fractal strings and spectral-operator invertibility are proved RH-equivalences (1501.05362),
  no new mechanism (none). Super-resolution: matrix-pencil fits recover damped modes (toy: error < 0.002 at 40 dB), so a
  'no damped mode in a band' certificate from psi data needs explicit model/noise bounds (N7/N10, idea). Ruelle
  resonances: Dolgopyat-type cancellation for twisted transfer operators gives pole-free strips (1811.05367 Thm 5.1);
  Mayer's operator gives the modular Selberg zeta (1008.4229): twisting Mayer operators by characters and proving
  non-integrability for a prime roof is a concrete zero-free-strip program (N7/N12, idea, nothing for zeta itself).
-/

namespace Fayed.R302
end Fayed.R302

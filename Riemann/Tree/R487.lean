/-!
# Round 487 — node N13 (Hilbert–Pólya and spectral approaches), round 15, REFUTATIONS of fresh claimed proofs: Prodromov v10.3 (an assumed spectrum–zero bijection; no Lean file), Witt ×2 (an index argument without a spectral correspondence; a 'confinement' divergence that is ζ's pole at s = 1), AL-Bahli (assumes von Koch's bound); Oukil v23 changes only an O(τ⁻²) term, so the R476 refutation stands; Gaussian correlation, Skorokhod embedding and minimax transfer no tool

* Node: N13 (15th round; R465 before). Round type: REFUTATIONS of fresh claimed proofs (Zenodo since 14 Nov, HAL since 8 Nov;
  texts pre-extracted with curl's default agent, lessons 228, 263, 274).
* [falsify] Prodromov v10.3 (Zenodo 24373726, 19 Nov; concept 22148892): Theorem 3.4 / Corollary 3.5 give only a distributional
  resolvent convergence; the bijection between the spectrum of the Krein-space operator and the zeros of ξ is assumed; no .lean file
  in the record despite the 'Lean 4 CIC DAG': FLAW (scout prodromov; REFUTATIONS).
* [falsify] Witt (Zenodo 24291648): §6 l.218–236, an 'orphaned zero' declared topologically forbidden; only real eigenvalues of a
  Hecke-type operator are shown (l.202–203), no spectrum–zero correspondence: FLAW (scout witt1).
* [falsify] Witt (Zenodo 24292077): Theorem 6.9, eq. (49), V(σ₀) ∼ −2(ζ′/ζ)(σ₀ + 1/2) → +∞ as σ₀ → 1/2+ is the pole of ζ at
  s = 1 (orchestrator read l.1212–1220): FLAW (scout witt2).
* [falsify] AL-Bahli (Zenodo 24284330, 16 Nov): Definition 2.1.2 (x = the (π(x) − 1)-th term of a progression with free x₁, d) is
  not supplied by Dirichlet's theorem, and (3.8) assumes von Koch's bound (RH): FLAW (scout albahli).
* [check] Oukil, HAL hal-05552822 v23 (14 Nov) against v22 (refuted R476): the only change (diff of pdftotext) adds
  (π/τ)Ψ[η♯](1)∫₀¹ sin(2πv)g_r(2π/τ, v)dv to ℑ I(s); g_r(t, v) = t(1 − 2v) + O(t²) as t → 0 (expansion, orchestrator), so the term is
  O(τ⁻²) and the limit passage 'by Lebesgue's theorem' is v22's verbatim: the R476 refutation applies unchanged (the scout read
  the repaired formula as closing the gap but did not test the Lebesgue step; REFUTATIONS line extended).
* [check] Bonelli, HAL hal-05865650 (16 Nov, 'Kaleidoscopic Wave Mechanics of the Riemann Zeta Function'): the HAL record has
  no file (document URL 404, API files_s empty): not testable.
* Exploration: Royen's Gaussian correlation inequality (needs the Weil form already PSD: circular), Skorokhod embedding / Doob's
  maximal inequality (centred zero counting is not a martingale), von Neumann's minimax (simplex duality loses the negative
  directions of a signed quadratic form; R363 tried game theory): no transfer.
* Scouts: 8, all answered.
-/

namespace Fayed.R487
end Fayed.R487

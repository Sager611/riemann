/-!
# Round 476 — node N5 (Nyman–Beurling / Báez-Duarte), round 15, the unseen HAL claimed proofs (Nov 2026): Robart's ℓ²-summability of Σ_{k≤n} μ(k)/k is false (an H² argument at the first zeta zero; exact sieve to 10⁷), Oukil's inequality (4) contradicts his own hypothesis H2, Agélas v2's Proposition 4.1 fails at s = 3/4 for the η coefficients, Arbai's determinant identity contradicts itself and his Part III conclusion rests on a proposition proved only above a height: four entries in REFUTATIONS; N5 records stand; Beurling–Malliavin, Fock-space sampling and Hörmander's ∂̄ estimates transfer no tool

* Node: N5 (15th round; R461 before). Round type: the unseen HAL claimed-proof notes (queued since 12 Nov; HAL search on
  'Riemann hypothesis' shows nothing newer than 8 Nov); PDFs fetched with curl's default agent (HAL's bot check answers
  browser agents), pdftotext into .scratch/r476/txt; Agélas v3 is under embargo, v2 used.
* [falsify] Robart, HAL hal-05850770 v1: Theorem 2 (Σ γ_n² < ∞, γ_n = Σ_{k≤n} μ(k)/k) is false unconditionally:
  Σ_{n≤N} γ_n² = ∫_1^{N+1} m², ∫_1^∞ m(x)x^{−s−1}dx = 1/(sζ(s+1)) would be H² on Re s > −1/2, impossible at
  1/2 + 14.1347i; sieve: S(10⁴) = 1.52994, S(10⁶) = 1.66344, S(10⁷) = 1.72904 (slope ≈ 0.029 per log N).
  Broken step: S(u, r) = Σ_{k≤u} k^{3/2}μ(k)r^k called bounded on [1, ∞) × [0, 1]: FLAW (REFUTATIONS).
* [falsify] Oukil, arXiv 2603.15716 v18 / HAL v22: τ·Im I(σ + iτ) → Ψ[η♯](1) along every τ → ∞ (integration by parts,
  Riemann–Lebesgue; G′ ∈ L¹ by (1)); (4) would bound it by |C(r)| and so contradicts H2; for η = {t}: Ψ[η♯](1) =
  0.0792909824, H2 margin 6.1·10⁻⁶ at r = 0.1 (closed form K(x) = 1/(2x) − 1/(x−1) + ζ(x)/x); grid check of τ·Im I
  at r = 0.1, 0.3, τ ≤ 400: values around 0.079 with ζ-driven fluctuations: FLAW (REFUTATIONS).
* [falsify] Agélas, HAL hal-05060598 v2: Proposition 4.1's (30) fails for a_n = (−1)^{n+1} at s = 3/4 (difference
  7.69 … 6188 for p_k = 97 … 99991, ratio to the bound 6.0 … 27509); the scout's Davenport–Heilbronn test satisfies
  Proposition 4.2's hypothesis (bounded partial sums) but P_k(s, a_DH) is a sum of two Euler products, so no
  contradiction there; the decisive failure is (30)/(64): FLAW (REFUTATIONS).
* [falsify] Arbai Part I (hal-05852361): Theorem 3.5 (3.10) with Ω non-vanishing makes det_Fred(I − L_s) infinite at
  s = ρ/2 and contradicts its own next sentence; Z_even undefined; κ(G_N) ≥ e^{c√N} asserted without proof or source.
* [falsify] Arbai Parts II–III (hal-05852389, hal-05855807): Part III Theorem 5.6(iii) claims RH; (5.29)–(5.30) for all
  τ₀ use Proposition 5.5 (proved for |τ₀| ≥ T*(δ, θ)) and an unquantified distant-pair bound: FLAW (REFUTATIONS,
  one entry for the program).
* [explore] Beurling–Malliavin completeness radius → N5: finite zero sets have BM density 0 (Gram of the first four
  zeros at T = 30.42 has rank 4); the infinite set needs spacing input beyond N(T): no new bound.
* [explore] Seip–Wallstén Fock-space sampling → N1: the planar density of the zeros is 0: no transfer.
* [explore] Hörmander's L² ∂̄ estimates → N1: no canonical plurisubharmonic weight or ∂̄-closed datum for Weil's
  kernel: no transfer.
-/

namespace Fayed.R476
end Fayed.R476

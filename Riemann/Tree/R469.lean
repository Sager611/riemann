/-!
# Round 469 — node N7 (zero-free regions and zero-density constants), round 17, the unseen Zenodo N7 records (no date filter, 225 records, 11 unseen on topic): Murphy's Lean-kernel-checked Λ ≤ 0.2811856084 for ζ (Liu's half-plane composed with de Bruijn's strip theorem, standard axioms) HOLDS (N2 line); Rafiq's rerun of Bellotti's log-free density code reproduces 31 of 38 table rows, the other 7 printed constants are larger than the code's, and it does not touch our factor-2 and signed-sum flaws (REFUTATIONS line updated: published in JNT 2025); Güvenç's three-wall barrier is explore4's NOGO line; Missiaen's Paper VII gives no record; no N7 record moves; Fekete points, Tverberg and Markov–Bernstein transfer no tool

* Node: N7 (17th round; R466 before). Round type: [check] of the unseen Zenodo records next to the node (lessons 216, 227): Zenodo
  without a date filter on "zero-free region", "zero free region", "zero-density", "zero density", "density hypothesis",
  "quasi-Riemann", "Vinogradov-Korobov", "de la Vallee Poussin", "zero-free half-plane", "Siegel zero", "Guth-Maynard", "Guth
  Maynard", "large values", "explicit zero-free", "zero-free regions", "log-free", "Ingham", "Mossinghoff", "Korobov": 225 records,
  11 unseen on topic; texts pre-extracted (lesson 228).
* [check] P. W. Murphy, Zenodo 23324360 v1.1.0 (12 Oct 2026), 'A kernel-checked bound Λ ≤ 9/32 for the de Bruijn–Newman constant from
  OpenAI's quasi-Riemann hypothesis; v1.1: Λ ≤ 0.2811857 from Liu's improved half-plane': ArdaDBNLiu.lean l.23–26
  dbn_real_zeros_of_liu_unconditional with H t z := ∫₀^∞ e^{tu²}Φ(u)cos(zu) du (challenge file l.30–35) and t ≥ 2(θ − 1/2)²,
  θ = (1507 − 2√921)/1653, value 0.2811856084; axioms logged: propext, Classical.choice, Quot.sound; no sorry in the solution: HOLDS.
  With PR 141's θ₁₄₁ = 0.8733458561067185768 the formula gives 0.27877425654 (orchestrator, Decimal; = COROLLARIES R466's Λ_χ value);
  9/32 at θ = 7/8. N2 line updated (a Lean-kernel-checked Λ for ζ; the record Λ ≤ 0.107559925318 is not Lean-checked).
* [check] A. Rafiq, Zenodo 23131177 (4 Oct 2026), verification note on Bellotti, 'An explicit log-free zero density estimate for the
  Riemann zeta-function' (arXiv 2405.12545; J. Number Theory 269 (2025) 37–77): the supplementary code with Table 4's parameters
  reproduces 31 of 38 Table 2 rows; 7 rows differ by −0.056% to −0.829% (code below the printed C: exp(50)–exp(60) 128806.85 vs
  128461.4; exp(170.2)–exp(500) 1.12·10⁶ vs 1.1107·10⁶; ...), so the printed constants are the larger ones; a prose typo ((u + w)⁴ for
  (v + w)⁴, confirmed by the author), an argument-order slip (fracd2 calls contrd2(u,v,w,x), signature (u,x,v,w); small effect) and
  double-precision scipy quadratures under mp.dps = 10. Scout FLAW read in the note by the orchestrator: a reproducibility note, not a
  disproof; it does not address the derivation flaws of our REFUTATIONS line (factor 2 in c₄, signed sum for |G(1/2 + it)|, R160),
  which stay open; REFUTATIONS line updated with the publication and Rafiq's findings.
* [check] G. Güvenç, Zenodo 21764256 (3 Aug 2026), 'Three-Wall Barrier Map for Zero-Density Estimates at the Guth–Maynard Frontier':
  wall W2 caps the exponent at σ = 7/10 (Thm 4, §VII): KNOWN, explore4's NOGO line (30/13 pinned at (σ, τ) = (7/10, 6/5)).
  A. Rafiq, Zenodo 22268018 (3 Sep 2026), review note on Turnage-Butterbaugh's survey (arXiv 2607.04632): flags its range for
  Montgomery's density-hypothesis theorem for ζ (σ ≥ 9/10, not 3/4; §4, l.440–450); a survey misstatement, not checked here.
* [check] T. Missiaen, Zenodo 24141101 v7.1.1 (11 Nov 2026), Paper VII (mirror Poisson–Hadamard kernel, quasiconcave minorant,
  'rigorous Trudgian Stieltjes bounds'; Thm main_p7 l.65–82): no explicit zero-free region for ζ; its R_eff > 5.558691 (l.80) is
  weaker than R₀ = 4.896 and than the Lean-checked half-plane; its verification script needs matplotlib (not run): no record.
* Not checked: T. Marechal, Zenodo 19128923 (Mar 2026), 'transverse curvature' zero-free framework (its zero-free conclusion is
  conditional on a 'Strong Spacing Conjecture'; scout gave no answer). Not read beyond titles: claimed proofs 16242777, Hodge-via-
  Korobov 14928133, Kreĭn-string 20114086, comparative sign changes 23597138, mixed zero measures 23861147. N7 records unchanged.
* Exploration: Fekete points and transfinite diameter (normalized-zero Vandermonde means V₂, …, V₆ = 0.3276, …, 0.2373 only obey the
  tautological V_n ≤ 1/4 type bounds), Tverberg's theorem (affine partitions of zero points carry no constraint on Re ρ), Chebyshev–
  Zolotarev/Markov–Bernstein (a polynomial surrogate of Hardy's Z has approximation error above the zero-free margin: −2.756 on
  [0, 14]): no transfer.
* Scouts: 8, 7 answered (marechal empty); the decisive values (Λ formula at θ₁₄₁; the direction of Rafiq's 7 rows) checked by the
  orchestrator.
-/

namespace Fayed.R469
end Fayed.R469

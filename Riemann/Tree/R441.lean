/-!
# Round 441 — node N8 (critical-line zeros, Gram-point statistics), round 14, audit of Jerby's Gram-point program after our counterexample to his repulsion bound: the claimed proof of Montgomery's pair correlation conjecture under RH (arXiv 2511.18275, Thm 9.2) treats the deterministic sequence PC_N(f; 1̄) as a stationary strongly mixing process because covariances decay (l. 4134–4141; also S′(t), l. 3888–3897) and bounds local zero counts by the average spacing (l. 4016–4024): broken; the 'corrected Gram's law ⟺ RH' of 2406.06548 (Thm B) proves only RH ⇒ the law (l. 458–465); Rosser-rule violations are frequent (Gourdon: 32.06 per million to 10^13; Brent–van de Lune–te Riele–Winter list zero-free length-2 blocks), so our counterexamples sit in a known block type; 2610.21952 states 11 laws, 6 finite (a confined-pair liminf law, l. 653–663); no prior test of μ > 4; multigrid, Allan/Hadamard variance and Erlang formulas transfer no tool

* Node: N8 (14th round; R425 before). Round type: falsify the other claims of a program whose central conjecture we just refuted (jerbyfals:
  μ = 3.968 at n = 1,565,734,083 and 3.945 at n = 2,565,007,238; DISCOVERIES).
* [falsify] 2511.18275 (e-print of 22 Dec 2025, Disc-PCC-15.tex), Thm 9.2: under RH, for every Schwartz f, the pair correlation of the zeros
  of Z follows the sine-kernel law. Proof, l. 4134–4146: 'Cov(PC_N(f; 1̄), PC_{N'}(f; 1̄)) → 0 ... Hence the sequence {PC_N(f; 1̄)} forms a
  stationary, strongly-mixing process ... By the ergodic theorem for mixing sequences'; l. 3888–3897 read Cov_T(S′(t), S′(t+u)) = O(1/log T)
  as strong mixing; l. 4016–4024 bound the neighbours of a zero within 2πR/log N by O(R) from the average spacing (under RH the count is
  R + O(log N/log log N)). Covariance decay of a deterministic sequence gives neither stationarity nor mixing, so the transfer from the
  randomized model to ζ is unproved: FLAW (lines read by the orchestrator; reviewer CORRECT, SIGNIFICANT; OpenAlex: cited only by the
  author's two later papers) → REFUTATIONS.
* [falsify] 2406.06548, Theorem B (mainthm:B2, l. 452–456: for any n, RH ⟺ (−1)^n Δ_n(1̄) > 0 with Δ_n(ā) := Z_{N(g_n)}(g_n(ā); ā)): the proof
  (l. 458–465) uses Theorem ESHO (RH ⟺ non-colliding paths from 0̄ to 1̄) and Δ_n(0̄) = (−1)^n, so it argues only RH ⇒ the corrected law;
  the converse is not argued (reviewer CORRECT, SIGNIFICANT) → REFUTATIONS. The paper does not claim that the second-order form is positive
  (l. 665, 723); at n = 1,565,734,083 it is (−1)^n Z(g_n) + (Z′(g_n)/log(g_n/2π))² = −2.5171.
* [check] Rosser's rule: first violation at n = 13,999,525 (t ≈ 5.35·10⁶; Trudgian 2011 §7.3 prints 13,999,825); Gourdon (2004): 320,624,341
  violations up to 10^13 (32.06 per million); Brent–van de Lune–te Riele–Winter (1982, Table 6): 43 zero-free length-2 blocks (pattern OO) in
  [g_156,800,000, g_200,000,000). Our counterexamples sit in this known block type.
* [harvest] 2610.21952 (Clouds.tex): 11 stated laws or conjectures, 6 finite: the repulsion bound (l. 313–326, refuted), a confined-pair law
  liminf (2π/log(g_n/2π))·μ(g_n) = 2 over isolated bad points with a confined pair (l. 653–663), the trunk-slope law (l. 1831–1851), the
  lower-scale law (l. 2002–2041) and two Davenport–Heilbronn numerics (l. 762–800, 840–855).
* [check] adversarial novelty: no work besides Jerby's tests μ at Gram points (arXiv API, lst.py, Zenodo 'Gram point', OpenAlex); 2610.21952
  restates the bound as open.
* Exploration (tool transfers, tested): multigrid local Fourier analysis (a damped-Richardson smoothing factor is no ζ multiplier and exceeds
  1 at the first zero), Hadamard/Allan variance (generic squared-difference positivity; Allan variance tried R422), Erlang B/C monotonicity (no
  intrinsic load or server count): no transfer.
-/

namespace Fayed.R441
end Fayed.R441

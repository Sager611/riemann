/-!
# Round 468 — node N1 (Weil positivity), round 11, the unseen Zenodo N1 records (no date filter, 110 records, 32 unseen on topic): De Jesus's certified two-window bracket agrees with ours (c = 5: intersection [8.2734, 9.2810]·10⁻¹⁸); Al Mohder's negative Connes–van Suijlekom Galerkin matrices for Dirichlet L are a truncation effect by the paper's own account; Mittermeier's conductor-23 architecture and Ursachi's criterion hold as stated; Zenodo 22062591's certificate is incomplete (raw Arb data missing); Franchi Viceré's claimed proof breaks at Theorem 2 (REFUTATIONS); Chekmachev's c ≤ 47 stands; Freud weights, circle packing and Edmonds' polytope transfer no tool

* Node: N1 (11th round; R445 before). Round type: [falsify]/[check] of the unseen Zenodo records next to the node (lessons 216, 227):
  Zenodo without a date filter on "Weil positivity", "Weil quadratic form", "Weil criterion", "Weil explicit formula", "explicit
  formula" "positivity", "Connes" "positivity", "Weil functional", "Weil distribution", "Bombieri" "Weil", "Yoshida" "explicit
  formula": 110 records, 32 unseen on topic; 12 texts pre-extracted (lesson 228); Spark's 24193168 (13 Nov) has no files.
* [check] E. De Jesus, Zenodo 24171237 (12 Nov 2026), 'Certified Weil Positivity in Two Windows': certified two-sided bounds for the
  spectral infimum, L = 0.8: [1.43300048827856, 1.674]·10⁻¹⁷; L = log √5 = 0.80471895621705018: [8.27343868254787, 9.475]·10⁻¹⁸.
  Orchestrator: our WeilSimple bracket at c = 5 (CCM's QW_λ, REPORT l.226) is [6.1089, 9.28104]·10⁻¹⁸: consistent, intersection
  [8.2734, 9.2810]·10⁻¹⁸ (A401172's a(5) = 17 either way); below the record: HOLDS.
* [check] Zenodo 22105737 (Aug 2026, positivity up to a = (1/2) log 5): complete certificate; its boxed λ_min(S₊) =
  5.35928198745683500559·10⁻¹⁸ (TeX l.1411–1428) is the Schur minimum of a proof majorant, 'not claimed to equal the spectral bottom':
  no conflict with ours. Zenodo 22062591 (up to support radius (1/2) log 3): the raw Arb JSON is missing by the record's own notice
  (l.7–14, 50–52) and its Schur margins are no eigenvalue bound: the certificate is incomplete (below the record; no proof claim).
* [falsify] A. Al Mohder, Zenodo 21914853 (13 Aug 2026, 'Ground-state negativity of truncated Connes–van Suijlekom forms for Dirichlet
  L-functions', Thm 9, Cor 10): negative even-sector CvS Galerkin matrices (N = 20, T = 800, c = 7–37; χ₅, χ₋₃: −0.05128, −0.00297 at
  c = 7; −0.67349, −0.79166 at c = 29); the paper says this is consistent with GRH (the truncated CvS form is not a compression of
  Weil's form; no N → ∞ or c → ∞ claim): HOLDS, no GRH conflict.
* [falsify] R. A. Mittermeier, Zenodo 24186091 (13 Nov 2026, 'A Conductor-23 Operator Architecture for Weil Positivity'): no RH
  claim; Theorem 44 (zero-floor obstruction) rules out a coercive sign transfer; eq. (37) η(z)η(23z) = ½(Θ_{Q₀} − Θ_{Q₁}),
  Q₀ = x² + xy + 6y², Q₁ = 2x² + xy + 3y², checked to n = 500 (0 mismatches): HOLDS.
* [falsify] A. Ursachi, Zenodo 24192507 (13 Nov 2026, 'Local Weil Positivity and Signed Prime-Deficit Criteria'): Thm 6.1 claims
  RH ⟺ W(X) ≥ −C eventually; the 'prime-square density' null is ∫e^{−b/2}h_x(b) db = 0, not a termwise cancellation; the smallest
  instance (bump K, h = 1/3) reproduced by integration by parts (∫F₀ ≈ −2.7·10⁻³⁴, null ≈ −7.6·10⁻³⁶); nothing beyond the record:
  HOLDS as stated.
* [falsify] C. Franchi Viceré, Zenodo 19546495 v6.5 (12 Apr 2026), claimed RH proof by 'semilocal spectral descent': breaks at
  Theorem 2 (l.123–145): the identity QW_λ^N(g ∗ ḡ♯) = Σ_k g̃(ρ_k)g̃(1 − ρ̄_k) over the spectral zeros of the truncated operator
  σ_S(N, λ) is not in the cited Connes 1999 Thm 4 or CCM 2511.22755 §5; with Theorem 1 (real spectrum, fine) it would make every
  truncated Weil matrix positive semidefinite, which is Weil positivity itself (circular); the S-limit of Theorem 4 ('Euclid
  guarantees the exhaustion') is harmless, since the local terms stabilize for fixed g (its own Theorem 3). Scout FLAW, checked by
  the orchestrator in the text: REFUTATIONS.
* Withdrawn: A. Fernandez, Zenodo 23991459 v2 (5 Nov 2026), withdraws the RH-proof claim of his Gaussian Weil–Toeplitz certificate
  (Zenodo 17567141): REFUTATIONS (status line).
* Not read beyond titles: 24 further unseen records (frameworks; 23440529, an audit of Desogus's Three Gates, already in REFUTATIONS;
  Davenport–Heilbronn negativity 23873865; Shoeib's Fourier-multiplier obstruction 21455507; claimed proofs 21184595, 21218854).
  Record (Chekmachev, c ≤ 47) stands.
* Exploration: orthogonal polynomials for Freud weights (Christoffel extremal principle → N3: the Hankel matrix (λ_{i+j+1})_{i,j<2}
  of Li coefficients has determinant −0.0037321667 < 0, so RH does not make them a moment sequence), Koebe–Andreev–Thurston circle
  packing (no canonical triangulation from zeros; Rodin–Sullivan needs domain data), Edmonds' matching polytope (inequalities only
  after arbitrary zero-gap encodings): no transfer.
* Scouts: 8, all answered; the decisive comparison (c = 5 brackets) and the Franchi Viceré step checked by the orchestrator.
-/

namespace Fayed.R468
end Fayed.R468

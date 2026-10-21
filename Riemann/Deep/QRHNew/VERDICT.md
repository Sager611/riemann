KNOWN: Rubinstein–Restrepo Borrero, Hexagon 2610.00220v1 (13 Oct 2026); Hu, Hexagon 2610.00216v1 (13 Oct 2026)

- **Claim checked.** L_K(s,ψ) ≠ 0 for Re s > 11/12 for unitary ψ of K = Q(√−3) with angular type m ≠ 0. This gives
  17/12 for CM elliptic curves and (w−1)/2 + 11/12 for CM newforms.
- **R–RB.** Part I: OpenAI's 11/12 argument with ∂_z̄²θ gives L(E_d,s) ≠ 0 for Re s > 17/12, E_d: y² = x³ + d. This is
  our G1 for j = 0, by the same route. Part II: 7/8 for ᾱη (every finite-order η), so 11/8 for every E with CM by K.
- **Hu, Thm 1.2 and §9.** 7/8 for every unitary Hecke character of every number field, infinite order included (the only
  pole is 1 + it for N^{it}); CM motives at w/2 + 7/8. This implies all of Theorem G and G1–G3, with 7/8.
- **Repairs (NOTE.md §2–3), each re-derived against paper2's TeX:**
  - The power is k = 1 + m, since (4.15) carries ν̄(z₁)ν(z₂). Conjugating when m < 0 gives k ≥ 2, so no radial or
    negative-k case.
  - The multiplier must be conj κ(g_{1,h}); without it B_{p,1} ≠ χ_p³ and the large sieve fails. Also ᾱ(c_h)^{2k}.
  - The dual Mellin factor is (i^k/4)(2π)^{−2s}, not (27/(2π)²)^s; the kernel V♯_k and |C_k| ≤ 1/81 follow.
  - The principal-case pole is at 1 + it, handled meromorphically as in paper2 §3.
  - Only m ≠ 0 went beyond paper2: finite-order characters (dihedral weight 1) and N^{it} twists were already covered.
- **Closure re-checked.** (4.15), (5.3), (5.8), (7.1), (7.3) and Lemma 7.2 keep the same k. **Matches R–RB at k = 2:**
  (8.4) factor (cv)⁻⁴, kernel Γ(4/3±t)Γ(5/3±t), (2π)^{−2s}, conj κ in (B.16).
- **Numerics.** `angular_check.py` (16 s, 22 MB) checks Lemma A for k = 1–4, that the constant term is killed, the
  Mellin constants and the one-mode functional equation, each to ≤ 1e−16.
- **Search.** arXiv HTML search (the API was rate-limited), the Hexagon API, OpenAI's 749 preprints, Liu ×2,
  Xu–Sang, Cureton–Thorner, FSW, the web and the evidence API. Only the two Hexagon papers state the result; the rest
  treat finite order only.
- **Filing.** Not DISCOVERIES. Not REFUTATIONS (the flaws were in our own draft). Suggested NOGO line: "QRHNew: known
  (Hu; R–RB, Hexagon 13 Oct)". The B♭ record hunt stays negative. What remains is a corrected general-k write-up of the
  11/12 route, weaker than Hu's result.

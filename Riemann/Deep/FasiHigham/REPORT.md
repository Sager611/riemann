THEOREM: the m = 3 case of the Fasi–Higham positivity conjecture holds for every k: c_n(k,3) > 0 for all k, n ≥ 0 (short proof); no diagonal counterexample: c_n(m,m) > 0 for ALL n ≥ 0 for every odd m ≤ 201 (rigorous, computer-assisted).

Context: arXiv 2610.39051 (read in TeX) refutes the conjecture for k ≫ m and proves k = 0, m = 1 and
eventual positivity; it leaves the diagonal k = m and every m ≥ 3 open. Notation, proofs, data: NOTE.md.

New, rigorous
1. m = 3, every k (NOTE §2): with N = k+1, (k+1)(k+2)(k+3) q = y³ + 3Ny + 2N, y = N − z. The positive
   zero dominates (|z±|² − z₀² = 3N(1−d) > 0) and 2|R±| < R₀ reduces to one exact polynomial identity
   plus g(d) = −3d³ + 12d² − 30d + 12 < 0 on [1/2, 2/3). Tight: 1 − 2|R±|/R₀ ≈ 1/(3(k+1)). Not in the
   paper; no other source found (arXiv search, evidence API).
2. Diagonal: c_n(m,m) > 0 for all n, odd m ≤ 201 (Arb balls to n₁ ≈ 0.73 m², then a certified pole
   bound). The earlier check (m ≤ 61, n ≤ 12m) never reached n₁; verify.py re-proves m ≤ 101 with mpmath.
3. Columns: for every odd m ≤ 101, every k < k_min(m) gives c_n(k,m) > 0 for all n; this covers all
   k ≤ 2m + 11, hence the whole near-diagonal, and makes the k_min table below exact.
4. First failures k_min(m) for every odd m ≤ 101, with the first negative index (NOTE §5):
   k_min(m) = 2m + c, c ∈ {12, …, 18} (c = 13 for 35 ≤ m ≤ 71); first negative at n = k + m − (2…6).
   k_min(m) − m = m + c > 0 and grows: it never changes sign. (k_min(9) = 42, below the paper's 45.)

Understanding (heuristic, fits all data)
- For n ≤ k + m, c_n = [z^n] e^z/p_{k,m}(z) (Padé condition); for odd m, 1/q = Σ_ℓ e^z z^{(k+m+1)ℓ}E^ℓ/p^{ℓ+1}
  with E having positive coefficients. A saddle-point competition predicts failures iff
  k/(k+m) > σ*(n/(k+m)), σ* = 0.90, 0.82, 0.72, 0.68 at 0.6, 0.8, 0.95, 0.99 and → 2/3 as n → k+m.
  This explains k_min ≈ 2m + O(1) and n ≈ k + m; checked against 4 independent cases.
- Diagonal: the oscillating part is below e^{−0.34m} of the main term for n ≤ 2m; beyond, the closest
  approach is an O(1) dip at n ≈ 2.17m (relative dip 0.73–0.83 for 5 ≤ m ≤ 197, no trend). Partial
  fractions mislead here (residues ~10¹⁸ at m = 61 cancel to c₀ = 1).
NO RECORD for the diagonal itself: a proof for all odd m needs uniform-in-m error terms for that saddle
analysis and an explicit pole tail (NOTE §4). All evidence says the diagonal conjecture is true.
Files: Riemann/Deep/FasiHigham/NOTE.md, verify.py (stdlib + mpmath, ≈ 3 min, prints ALL OK),
arb_cert.py (python-flint; reruns the m ≤ 201 and column certificates).

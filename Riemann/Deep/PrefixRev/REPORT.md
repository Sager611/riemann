DISCOVERY: For 2 ≤ k < m < n with m+k < n, ⟨r_n,r_m,r_k⟩ = S_n iff n is odd, m+k = n−1, gcd(m+1,k+1) = 1, (k+1)∤m, and not (n ≡ 1 mod 4 and 4 | k).

What it settles. Conjecture 4 of Blanco et al. (arXiv 2511.16959) says m+k < n forces a proper subgroup.
arXiv 2610.30892 proved this except on the line n odd, m+k = n−1. On that line the literal conjecture
is false (smallest counterexamples (7,4,2), (11,6,4), (13,10,2); 318 with n ≤ 101; infinitely many,
e.g. (n, n−3, 2) for n ≡ 1 mod 6). The theorem above answers the whole line exactly.

Full structure on the line n = m+k+1 (α = k+1, β = m+1, g = gcd(α,β)):
- orbits are the classes x ≡ ±c (mod g), so ⌊g/2⌋+1 orbits (any n);
- n odd, g = 1, α | n: imprimitive, blocks = intervals of length α, G ≅ D_{2α} ≀ D_{2n/α};
- n odd, g = 1, α ∤ n: primitive and ⊇ A_n; A_n iff n ≡ 1 mod 4 and 4 | k, else S_n.

Proof idea. P = r_k r_n r_m is the first-return map of x ↦ x+α on Z/(n+1) minus 0, hence an n-cycle
when g = 1. In the coordinates where P is the shift, r_n is a reflection and r_k a reflection
restricted to a Beatty set; block systems must be residue classes mod d, and a lap-counting argument
shows one survives iff α | n. Then Jones's classification of primitive groups with an n-cycle (CFSG)
plus a short fixed-point bound for collineations (r_k fixes more than n/2 points) gives A_n.

Verification (verify.py, 31 s, 62 MB, exit 0): every claim for all n ≤ 101; G ⊇ A_n there is also
proved without CFSG via explicit Jordan p-cycles; sympy group orders for n ≤ 25 and the wreath orders
for n ≤ 51; the fixed-point bound brute-forced on 9 small projective spaces.

Novelty check: arXiv search, full math.GR and math.CO listings Sep–Oct 2026, GXL evidence API:
nothing beyond the two papers above. Caveat: the counterexamples alone are minor (2610.30892 flags the
wording, and Theorem 1 of 2511.16959 already gives (n, n−3, 2), n ≡ 1 mod 6); the new part is the
exact classification of the boundary line.

Files: NOTE.md (statements, proofs, tables), verify.py.

DISCOVERY: every primitive positive definite binary quadratic form of discriminant D represents a prime p ≪_ε |D|^{1.7498} (2B♭+ε; 7/4+ε from Hu's 7/8). The published record is |D|^694 (Thorner–Zaman 2017). New tool: an unconditional zero-density bound |D|^{2σ(1−σ)+ε} for ring class L-functions. Minor: least prime ideal in a ray class of an imaginary quadratic field ≪ (|d_K|N𝔪)^{3+ε} (was |d_K|^694·N𝔪^521); Piltz exponents α_k ≤ 0.8749 for all k (new for k ≥ 59).

# Combos deep dive: report (2026-10-20)
Proofs: [NOTE.md](NOTE.md). New results assume only Liu's claimed theorem (every finite-order Hecke L-function is zero-free for Re s > B♭ = 0.874887) or Hu's 7/8. The density theorems are unconditional.

## Main result (brief item 1)
- Thm 2: p(Q) ≤ C(ε)|D|^{2B♭+ε}, with an effective constant. For x ≥ C(ε)|D|^{2B♭+ε}, Q represents ≫ x/(h(D)log x) primes in (x/2, x].
- Cor 3: π_C(x) ~ Li(x)/h(D) for every ring class C, uniformly for |D| ≤ x^{0.5715−ε}. Before, an asymptotic was known only for |D| = x^{o(1)}. Cor 4: the least prime of the form x² + ny² is ≪ n^{1.7498}.
- Previous bounds: |D|^694 (Thorner–Zaman 2017, best published); |D|^{20/3} on average over D (Ditchen); |D|^A with A unspecified (Hexagon 2610.00178 §54, Oct 2026; Gaudet 2024); |D|^{3.996} from Liu or Hu with a trivial zero count (stated nowhere). GRH gives |D|^{1+ε}; for x² + ny² the least prime is trivially ≥ n.
- How it works: the ring class family has h ≈ |D|^{1/2} characters, each of conductor ≤ |D|.
  - Every ideal is m·(primitive ideal), and χ(m𝔞) = χ(𝔞). This gives a large sieve as good as for Dirichlet characters, with constant N + |D|^{1/2} (Lemma 2).
  - The approximate functional equation then gives a mollified second moment ≪ |D|^{1/2+ε} (Lemma 4), and a Carlson-type zero detector gives Σ_χ N(σ) ≪ |D|^{2σ(1−σ)+ε} (Thm 1).
  - The half-plane stops σ at B♭, so x ≥ |D|^{2B♭} suffices.
- The exponent equals 2B, so every improvement of the half-plane lowers it; B = 1/2 would give GRH's exponent 1. In Linnik's problem, by contrast, the density peak at σ = 5/7 makes B irrelevant.

## Other brief items
2. Ray classes: for imaginary quadratic K, (|d_K|N𝔪)^{3+ε} by the same method (Thm 8); for example, the least Gaussian prime ≡ α mod q has norm ≪ q^{6+ε}. For general K, only the substitution (h_𝔪 log)^{7.99}.
3. Goldbach: KNOWN. Li–Liu (arXiv 2610.26756, 19 Oct) prove E(N) ≪ N^{1/2}(log N)^100 under any half-plane with Θ < 11/12. This supersedes Pintz–Althöfer's 7/12, which was limited by the 19/24 minor-arc exponent.
4. Least non-residue and least prime residue: (log q)^8 is Hexagon 178 Lemma 1.6; Liu gives 7.993, a substitution. The exponent 1/δ is a barrier.
5. Short intervals, almost all intervals, Bombieri–Vinogradov: no change. The density maximum 30/13 sits at σ = 7/10 < 7/8, and zero methods cap Bombieri–Vinogradov at 1/2.
6. Class numbers: Hexagon 178 already has the explicit h(D) > 0.655√|D|/(πe log|D|) and an inexplicit √|D|/log log|D|. An explicit loglog constant is routine.
7. Piltz divisor problem: α_k ≤ B♭ for all k (Thm 9), a two-line corollary of Lindelöf right of B♭. It beats the database records (Bellotti 2023) for k ≥ 59, and every bound stated with an explicit range from k = 40.

## Significance (for the reviewer)
- Thm 2 is the real item: 694 → 1.75, well below the 4 that the 2026 inputs give directly. The method is classical (Carlson density plus the large sieve); the family-specific large sieve (Lemma 2) is the key step. No density theorem for ring class L-functions with a small exponent was found.
- Thms 8 and 9 are minor and routine. Do not list items 3–6 in DISCOVERIES.md.

## Novelty checks (NOTE §8)
arXiv API (least prime and quadratic form, ring class, class group zero density, quasi-Riemann, divisor problem, Piltz); Hexagon site search (7 terms); GXL evidence API (2 claims). Read Hexagon 178 §§1, 19, 54 and Thorner–Zaman 2017/2019; grepped FSW, Hu, Liu and the OpenAI preprints (749 titles).

## Leads
Gaussian primes in narrow sectors (Hu covers Hecke characters of infinite order; needs a Huxley-type density for {λ^m}); indefinite forms.

Suggested DISCOVERIES line: "Least prime represented by a binary quadratic form of discriminant D is ≪ |D|^{1.7498} (record |D|^694), from a new ring class zero-density bound and Liu's half-plane ([report](Riemann/Deep/Combos/REPORT.md))".

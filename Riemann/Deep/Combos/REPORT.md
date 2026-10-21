DISCOVERY (conditional on Liu's claimed half-plane (L)): every primitive positive definite binary quadratic form of discriminant D represents a prime p ≪_ε |D|^{2B♭+ε}, so p ≪ |D|^{1.7498} (7/4+ε from Hu's 7/8). The unconditional record is |D|^455 (Zaman's PhD thesis 2017; |D|^694 in Thorner–Zaman 2017). New tool: an unconditional zero-density bound Σ_χ N(σ;t) ≪ |D|^{2σ(1−σ)+ε}T³ per unit height window for ring class L-functions; no published bound implies it ([NOVELTY.md](NOVELTY.md)). Minor, also under (L): least prime ideal in a ray class of an imaginary quadratic field ≪ (|d_K|N𝔪)^{3+ε} (was |d_K|^694·N𝔪^521); Piltz exponents α_k ≤ 0.8749 for all k (new for k ≥ 59).

# Combos deep dive: report (2026-10-20)
Proofs: [NOTE.md](NOTE.md). New results assume only Liu's claimed theorem (every finite-order Hecke L-function is zero-free for Re s > B♭ = 0.874887) or Hu's 7/8. The density theorems are unconditional.

## Main result (brief item 1)
- Thm 2 (under (L)): p(Q) ≤ C(ε)|D|^{2B♭+ε}, with an effective constant. For x ≥ C(ε)|D|^{2B♭+ε}, Q represents ≫ x/(h(D)log x) primes in (x/2, x].
- Cor 3 (under (L)): π_C(x) ~ Li(x)/h(D) for every ring class C, uniformly for |D| ≤ x^{0.5715−ε}. Before, an asymptotic was known only for |D| = x^{o(1)}. Cor 4 (under (L)): the least prime of the form x² + ny² is ≪ n^{1.7498}.
- Previous bounds: |D|^455 (Zaman's PhD thesis 2017, Cor 1.4.1, best unconditional); |D|^694 (Thorner–Zaman 2017, best in a journal); |D|^{20/3} on average over D (Ditchen); |D|^A with A unspecified (Hexagon 2610.00178 §54, Oct 2026; Gaudet 2024); from Liu or Hu, |D|^{3.996} with a trivial zero count and |D|^{2.40} (fundamental D) with Pascadi–Thorner's family density (both stated nowhere). GRH gives |D|^{1+ε}; for x² + ny² the least prime is trivially ≥ n.
- How it works: the ring class family has h ≈ |D|^{1/2} characters, each of conductor ≤ |D|.
  - Every ideal is m·(primitive ideal), and χ(m𝔞) = χ(𝔞). This gives a weighted large sieve with constant N + |D|^{1/2}: the weight m(𝔞)^{1+2ε} sits on |c_𝔞|² and costs only N^{2ε} for divisor-bounded coefficients (Lemma 2). Published large sieves for this family have constant N + |D|.
  - The approximate functional equation then gives a mollified second moment ≪ |D|^{1/2+ε} (Lemma 4), and a Carlson-type zero detector gives Σ_χ N(σ;t) ≪ |D|^{2σ(1−σ)+ε}T³ per unit window, so T⁴ up to height T (Thm 1; only the |D|-aspect is claimed).
  - The half-plane stops σ at B♭, so x ≥ |D|^{2B♭+ε} suffices.
- The exponent is 2B + ε, so every improvement of the half-plane lowers it; B = 1/2 would give GRH's exponent 1 (+ε). In Linnik's problem, by contrast, the density peak at σ = 5/7 makes B irrelevant.

## Other brief items
2. Ray classes: for imaginary quadratic K, (|d_K|N𝔪)^{3+ε} by the same method (Thm 8); for example, the least Gaussian prime ≡ α mod q has norm ≪ q^{6+ε}. For general K, only the substitution (h_𝔪 log)^{7.99}.
3. Goldbach: KNOWN. Li–Liu (arXiv 2610.26756, 19 Oct) prove E(N) ≪ N^{1/2}(log N)^100 under any half-plane with Θ < 11/12. This supersedes Pintz–Althöfer's 7/12, which was limited by the 19/24 minor-arc exponent.
4. Least non-residue and least prime residue: (log q)^8 is Hexagon 178 Lemma 1.6; Liu gives 7.993, a substitution. The exponent 1/δ is a barrier.
5. Short intervals, almost all intervals, Bombieri–Vinogradov: no change. The density maximum 30/13 sits at σ = 7/10 < 7/8, and zero methods cap Bombieri–Vinogradov at 1/2.
6. Class numbers: Hexagon 178 already has the explicit h(D) > 0.655√|D|/(πe log|D|) and an inexplicit √|D|/log log|D|. An explicit loglog constant is routine.
7. Piltz divisor problem: α_k ≤ B♭ for all k (Thm 9), a two-line corollary of Lindelöf right of B♭. It beats the database records (Bellotti 2023) for k ≥ 59, and every bound stated with an explicit range from k = 40.

## Significance (for the reviewer)
- Thm 2 is the real item: 455 → 1.75 under (L), below the 4.0 (trivial count) and 2.40 (Pascadi–Thorner, fundamental D) that the 2026 inputs give with published densities. The method is classical (Carlson density plus the large sieve); the family-specific large sieve (Lemma 2) is the key step. No published density for ring class L-functions implies Thm 1 ([NOVELTY.md](NOVELTY.md)).
- Thms 8 and 9 are minor and routine. Do not list items 3–6 in DISCOVERIES.md.

## Novelty checks (NOTE §8, [NOVELTY.md](NOVELTY.md))
arXiv API (least prime and quadratic form, ring class, class group zero density, quasi-Riemann, divisor problem, Piltz); Hexagon site search (7 terms); GXL evidence API (2 claims). Read Hexagon 178 §§1, 19, 54 and Thorner–Zaman 2017/2019; grepped FSW, Hu, Liu and the OpenAI preprints (749 titles).

## Leads
Gaussian primes in narrow sectors (Hu covers Hecke characters of infinite order; needs a Huxley-type density for {λ^m}); indefinite forms.

Suggested DISCOVERIES line: "Assuming Liu's claimed half-plane, the least prime represented by a binary quadratic form of discriminant D is ≪ |D|^{1.7498} (unconditional record |D|^455), via a new unconditional ring class zero-density bound ([report](Riemann/Deep/Combos/REPORT.md))".

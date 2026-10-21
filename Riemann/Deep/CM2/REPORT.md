DISCOVERY: the almost-all threshold for least primes in CM fields drops from 3/2 to (3+√69)/10 = 1.1307: for a fixed totally real field F and every A > (3+√69)/10, all but O(X^{1−c}) CM fields K/F with prime (or δ₀-rough) relative discriminant of norm ≤ X have prime ideals of degree one in every ideal class in every interval (x, x + x·D^{−κ}], x ≥ D^A (D = N𝔡_{K/F}); so every class has a prime of norm ≤ 2d_K^A. Tool: the first family large sieve over pairs (K, χ ∈ Ĉl_K) for CM fields over F at the family size X^{3/2}.

# CM2 report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Brief: `.scratch/orca/cm2-brief.md`. Paper proofs, not refereed, no Lean. Unconditional. F fixed, n = [F:Q].

## Theorems (complete proofs in the NOTE)
- **Thm 1 (family large sieve).** For the class I sequences b of the F-version of the RCApps3 detector (squarefree ideals of F, mollified, pole killed):
  Σ_{K CM, X/2<N𝔡_{K/F}≤X} Σ_{χ∈Ĉl_K} |Σ_𝔫 b(𝔫)λ_χ(𝔫)|² ≪ X^{O(δ)}(X^{3/2} + N³/X)N^{1−2σ}, N ≤ X^{1−4δ}. The sum of the per-field sieves (CMQuartic) is X·N·N^{1−2σ}.
- **Thm 2 (density over the family).** Σ_K #{zeros of L(s,χ), χ ≠ 1, σ ≤ β ≤ 9/10, |γ − t| ≤ 1} ≪ X^{(13−8σ)/(7−2σ)+O(δ)} (19/16 at σ = 5/6; was 5/4).
- **Thm 3 (least primes).** CMQuartic Theorem B for every A > A₀ = (3+√69)/10 (was A > 3/2). Same crossover as for F = Q: σ_c = 5 − √69/2.

## How it works (what is new beyond F = Q)
- Positivity over all O_F-orders O_𝔣 = O_F + 𝔣O_K of all CM K/F, with weight ω = |μ(O)|·[N(O^×):U²]; then a dual identity (Lemma 4): the family sum becomes a
  count of t ∈ 𝔩 (an ideal lattice in Rⁿ) and conductors 𝔶, summed over the unit classes u ∈ U/U². Class group of F: split b by classes, reference primes S₀.
- Poisson in t over an ideal lattice on a hyperbolic region: the dual lattice is an ideal lattice, so only the product of the widths matters (Lemmas 7, 9).
- Local densities at 𝔭 | 2 are bounded functions of μ mod squares (Lemma 10). Mellin in N𝔶 with ζ_F, stopping at Re z = δ/2.
- The unit sum turns the generator-dependent symbols [ν_𝔫/𝔮] into ray class characters (Lemma 11).
- Secondary terms: the classical large sieve for Hecke characters suffices, because the modulus range shrinks with κ (true for F = Q too:
  Heath-Brown's quadratic large sieve is not needed in RCApps3). Main terms: pole killed by (Z); Jutila-type mean value over F (Lemma 14).

## Checks (.scratch/cm2/, all pass)
- Dual identity exact: 14 CM fields (over Q(√5), Q(√3), Q(√2), Q(√10), Q(√15); incl. Q(ζ₅), Q(ζ₈), Q(ζ₁₂), class number 2, U⁺ ≠ U²) and 7 non-maximal orders.
- Lemma 1 (orders, 2-adic criterion) against pari's rnfdisc: 371 prime tests, 249 above 2, no mismatch.
- Poisson main term against brute force (F = Q(√5)): ratios 0.93, 1.015, 1.001.
- Family sums over CM fields of Q(√5): pole-killed sequences stay at the family size up to N = X; pole-carrying ones grow like N/√X.

## Novelty
- 49 arXiv API queries, 3 GXL claims: nothing implies Theorems 1–3. Closest: Miller–Peckner 1003.5336 (1-level density over CM fields, GRH),
  Fouvry–Iwaniec 2003 (F = Q), Goldmakher–Louvel 1112.1642 (not needed), Thorner–Zaman, An (fixed families). Previous best for CM fields: our CMQuartic (A > 3/2).

## Not done / open
- A₀ → 1 needs Thm 1 without N³/X (as for F = Q). All CM K/F (small genus conductors) and F varying: as in CMQuartic §7. Sectors, non-CM K/F: not attempted.

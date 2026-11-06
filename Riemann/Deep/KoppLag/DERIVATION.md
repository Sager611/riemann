# DERIVATION: unit-generated real quadratic orders with 2-torsion class group (Kopp–Lagarias Problem 1.7)

Agent kopplag, 2026-11-06. Files: `cert.py` (certificate), `forms.py` (exact form arithmetic), `checks.py` + `pari_checks.gp` (independent checks), `logs/`.

## 0. Setting and statement

N ∈ {+1, −1}; n ≥ 3 if N = +1, n ≥ 1 if N = −1; Δ = Δ_n = n² − 4N (> 0, not a square). O = O_Δ = ℤ[ω], ω = (Δ + √Δ)/2.
ε = ε_n = (n + √Δ)/2 ∈ O^×, ε' = (n − √Δ)/2, εε' = N. Δ = f²Δ₀ with Δ₀ fundamental, χ = (Δ₀/·) primitive mod Δ₀ ≥ 5.
Cl(O) = invertible fractional ideals / principal ideals (wide); Cl⁺(O): modulo principal ideals with totally positive generator.
These are exactly the real unit-generated orders (KL Prop. 2.3); Δ = 5 occurs twice (n = 3, N = +1 and n = 1, N = −1).

**Theorem.** Assume (L): |L(1, χ)| ≥ 0.11/max(3, log log q) for every primitive Dirichlet character χ mod q ≥ 3
(Littlewood78, `discoveries/littlewood-half-plane.md`, from the Lean-checked half-plane Re s > 0.8734).
(a) Cl(O_Δ) = Cl(O_Δ)[2] holds for exactly the 53 discriminants Δ = n² − 4 and the 37 discriminants Δ = n² + 4 of KL Tables 2A/2B
(89 distinct orders; largest Δ = 68640 = 262² − 4).
(b) Cl⁺(O_Δ) = Cl⁺(O_Δ)[2] (one class per genus) holds for exactly those, minus Δ ∈ {221, 725, 1517, 3965, 7565}: 84 orders.
(c) Without (L), and without GRH: (a) and (b) hold for every Δ < 10¹⁸ (elementary: Lemma A, genus theory not even needed).
(d) (L) may be replaced by any bound |L(1,χ)| ≥ c/max(3, log log q) with c ≥ 0.037.

## 1. Facts used

F1 (elements, norms). O = {(x + y√Δ)/2 : x, y ∈ ℤ, x ≡ yΔ (mod 2)}; α = u + vω has x = 2u + vΔ, y = v. Since Δ ≡ n² ≡ n (mod 2), x ≡ ny (mod 2).
Call α *primitive* if α ∉ kO for every integer k ≥ 2; equivalently gcd(u, v) = 1. Multiplication by α on the lattice O has determinant N(α), so [O : αO] = |N(α)|.
O ∩ ℚ = ℤ.

F2 (class number formula for orders; Neukirch, ANT I.12.12). h(O) = h(O_K)/[O_K^× : O^×] · #(O_K/fO_K)^×/#(O/fO_K)^×,
with #(O_K/fO_K)^× = f² ∏_{p|f}(1 − 1/p)(1 − χ(p)/p) and (O/fO_K)^× ≅ (ℤ/fℤ)^× (since f√Δ₀ ≡ 0 mod fO_K; KL l.840–850),
so the last quotient is f ∏_{p|f}(1 − χ(p)/p). Both unit groups are {±η^k}, hence [O_K^× : O^×] = R(O)/R(O_K), R = log of the fundamental unit.
With Dirichlet's h(O_K)R(O_K) = (√Δ₀/2) L(1, χ) (so L(1, χ) > 0, and (L) bounds L(1, χ) itself):
  h(O) R(O) = (√Δ/2) · L(1, χ) · ∏_{p|f} (1 − χ(p)/p).            (F2)
Checked: `logs/pari_checks_output.txt` (2080 orders, PARI lfun at 38 digits, |Q − 1| ≤ 1.2·10⁻³⁸ with the exact h from forms.py).

F3 (genus theory; Halter-Koch, Quadratic Irrationals 2013, Thm 6.5.2 = KL Prop. prop:genusorder, l.1543–1550).
|Cl⁺(O)/Cl⁺(O)²| = |Cl⁺(O)[2]| = 2^{μ(Δ)−1}, μ as in KL l.1527–1537; μ(Δ) − 1 ≤ ω(Δ) (μ = ω + 1 only when 32 | Δ).
Checked: `logs/checks_output.txt` C4, all 1998 orders with n ≤ 1000 (count of narrow classes with C² = 1 computed exactly).

F4 = (L).

F5 (forms; used only in the computation, §4). Proper equivalence classes of primitive forms of discriminant Δ ↔ Cl⁺(O)
via 𝔞 = ℤα₁ + ℤα₂ (positively oriented) ↦ N(xα₁ + yα₂)/N(𝔞) (Halter-Koch Thm 6.4.2; Cohen, CCANT Thm 5.2.9). Reduced forms:
|√Δ − 2|a|| < b < √Δ, i.e. 0 < b < √Δ and √Δ − b < 2|a| < √Δ + b; every proper class is one ρ-cycle of reduced forms (Cohen §5.6, Def. 5.6.2, 5.6.4, Prop. 5.6.6;
Buchmann–Vollmer, Binary Quadratic Forms, §6.10).

## 2. Lemma A (elementary)

**Lemma A.** Let n ≥ 1 (so ε > 1), Δ = n² − 4N ≥ 100 and α ∈ O_Δ primitive with m := |N(α)| < √Δ/2. Then m = 1.

Proof. (2.1 normalization) For a unit η = ±ε^j, |αη/(αη)'| = |α/α'| · (ε/|ε'|)^j = |α/α'| ε^{2j} (|ε'| = 1/ε); αη is primitive iff α is, and |N| is unchanged.
So we may assume 1 ≤ |α/α'| < ε² and α > 0. Then α² = |α/α'|·m gives √m ≤ α < ε√m, and |α'| = m/α ≤ √m.
Write α = (x + y√Δ)/2 as in F1. Then y√Δ = α − α' ≥ α − |α'| ≥ 0.
If y = 0, α = x/2 ∈ O ∩ ℚ = ℤ and primitivity forces α = ±1, m = 1. Let y ≥ 1. Then x = α + α' ≥ α − |α'| ≥ 0 and
  y√Δ = α − α' ≤ α + |α'| < (ε + 1)√m,  so  y < Y := (ε + 1)√m/√Δ;                  (A1)
  x = y√Δ + 2α' ≥ y√Δ − 2√m,  so  x + ny ≥ y(n + √Δ) − 2√m = 2yε − 2√m.          (A2)
Put z := x − ny (even, by F1). From 4N(α) = x² − Δy² = x² − n²y² + 4Ny²: z(x + ny) = 4N(α) − 4Ny², so
  |z| (x + ny) ≤ 4m + 4y².                                                       (A3)
(2.2 the inequality) Claim: φ(y) := yε − y² − m − √m > 0 for 1 ≤ y ≤ Y. Then 2yε − 2√m > 2(m + y²) > 0, and (A2), (A3) give |z| < 2,
so z = 0 (z even). Then x = ny, 4N(α) = 4Ny², m = y², and α = (ny + y√Δ)/2 = yε ∈ yO; primitivity forces y = 1, m = 1.
φ is concave, so it suffices to check y = 1 and y = Y (if Y < 1 there is no y ≥ 1). We use, for Δ ≥ 100:
ε ≥ (√(Δ − 4) + √Δ)/2 ≥ √Δ − 2/√Δ (as n² = Δ + 4N ≥ Δ − 4 and √(Δ−4) ≥ √Δ − 4/√Δ for Δ ≥ 4), ε < √Δ + 1, √m < Δ^{1/4}/√2.
  φ(1) = ε − 1 − m − √m > √Δ/2 − 1 − 2/√Δ − Δ^{1/4}/√2  (= 1.5639… at Δ = 100, increasing in Δ).
  φ(Y)/√m = (ε + 1)ε/√Δ − √m((ε + 1)²/Δ + 1) − 1 > √Δ − 4/√Δ − (Δ^{1/4}/√2)((1 + 2/√Δ)² + 1) − 1  (= 3.144… at Δ = 100, increasing),
using (ε+1)ε/√Δ ≥ ε²/√Δ ≥ (√Δ − 2/√Δ)²/√Δ ≥ √Δ − 4/√Δ and (ε+1)²/Δ < (1 + 2/√Δ)². ∎
Checks: C8 evaluates both lower bounds (arb) on Δ = 10^{t/100}, t = 200..4000: minima 1.5639 and 3.1440 at Δ = 100.
C1 searches the whole normalized window (y ≤ (ε+1)(Δ/4)^{1/4}/√Δ + 2, |x − y√Δ| ≤ 4) for n ≤ 3000 and 80 random n ≤ 10⁹: no counterexample.
C2: the principal cycle is {(1, n−2, 2−n), (2−n, n−2, 1)} (N = +1, n ≥ 4) resp. {(1, n, −1), (−1, n, 1)} (N = −1) for all n ≤ 20000,
the forms-theoretic form of Lemma A (Lagrange: |m| < √Δ/2 properly represented ⇒ m is a first coefficient on the cycle).
Sharpness: √Δ/2 cannot be replaced by √Δ: for n = p², N = −1, the split prime p has p² = n < √Δ and 𝔭² principal (C3b: n = 25, 49, 121, …).

## 3. The bound D₀

**3.1 Corollary B (split-prime criterion).** Let Δ = n² − 4N ≥ 100 with Cl(O) = Cl(O)[2]. Then no prime p with 4p⁴ < Δ (⇔ p² < √Δ/2) splits in O,
i.e. the Kronecker symbol (Δ/p) ∈ {0, −1} for all such p.
Proof. Let (Δ/p) = 1. Then p ∤ Δ, so p ∤ f, O is p-maximal and pO = 𝔭𝔭̄ with 𝔭 ≠ 𝔭̄ invertible of norm p
(𝔭 = pℤ + ((b + √Δ)/2)ℤ, b² ≡ Δ mod 4p, p ∤ b; 𝔭̄ = 𝔭 would need b ≡ −b mod 2p, i.e. p | b). 𝔭² is invertible of norm p² and primitive:
𝔭² ⊆ kO with k ≥ 2 forces k² | p², k = p, 𝔭² ⊆ pO = 𝔭𝔭̄, so 𝔭 ⊆ 𝔭̄ (cancel the invertible 𝔭), and equal norms give 𝔭 = 𝔭̄, impossible.
If [𝔭]² = 1 in Cl(O) then 𝔭² = αO with α ∈ 𝔭² ⊆ O; α is primitive (else 𝔭² ⊆ kO) and |N(α)| = [O : 𝔭²] = p² (F1).
For p² < √Δ/2, Lemma A gives p² = 1, a contradiction. So [𝔭] has order > 2. ∎
(For our families 2 never splits: n odd gives Δ ≡ 5 mod 8; n even gives 4 | Δ.)
Check C3: 58836 pairs (n, p) (all n ≤ 2000, 60 random n ≤ 10⁹): the reduced cycle of 𝔭∘𝔭 (Cohen Alg. 5.4.7) never contains a form with |a| = 1.

**3.2 Lower bound for h(O).** By (F2), with R(O) ≤ log ε_n (ε_n ∈ O^×, ε_n > 1, so the fundamental unit is ≤ ε_n):
- 2R(O) ≤ 2 log ε_n < log(Δ + 4): N = +1: ε_n < n = √(Δ+4); N = −1: ε_n < √Δ;
- L(1, χ) ≥ 0.11/max(3, log log Δ₀) ≥ 0.11/max(3, log log Δ) by (L) (Δ₀ ≥ 5, Δ₀ ≤ Δ);
- ∏_{p|f}(1 − χ(p)/p) ≥ ∏_{p|f}(1 − 1/p) = φ(f)/f ≥ 1/P_{ω(f)} ≥ 1/P_{ω(Δ)}, P_k := ∏_{i≤k} p_i/(p_i − 1)
  (the i-th prime of f is ≥ p_i and p/(p − 1) decreases; primes of f divide Δ).
Hence  h(O) > g(Δ)/P_{ω(Δ)},  g(Δ) := 0.11 √Δ / (log(Δ + 4) · max(3, log log Δ)).          (B1)

**3.3 g is increasing on [16, ∞).** d/dΔ log g = 1/(2Δ) − 1/((Δ+4) log(Δ+4)) − [log log Δ > 3]/(Δ log Δ log log Δ) > 1/(2Δ) − 1/(2.99Δ) − 1/(60Δ) > 0
for Δ ≥ 16 (log(Δ+4) ≥ 2.99; when log log Δ > 3, log Δ > 20); g is continuous at the switch.

**3.4 D₀.** If Cl = Cl[2] then h = |Cl/Cl²| ≤ |Cl⁺/Cl⁺²| = 2^{μ−1} ≤ 2^k, k := ω(Δ) (F3; Cl⁺ ↠ Cl). With (B1): g(Δ) < T_k := 2^k P_k, and Δ ≥ Π_k := p₁⋯p_k.
For each k, cert.py finds by bisection (arb, 200 bits) an integer D_k⁺ with g(D_k⁺) > T_k certified; by 3.3, every Δ ≥ D_k⁺ has g(Δ) > T_k.
So a 2-torsion order with ω(Δ) = k has Π_k ≤ Δ < D_k⁺; k is *excluded* if Π_k ≥ D_k⁺.
Computed (logs/cert_output.txt): k ≤ 14 possible, D₁₄⁺ = 20485358601780401 ≈ 2.05·10¹⁶ the largest; k = 15..60 excluded.
Induction for k ≥ 15: with p = p_{k+1} ≥ 53 and X = log Π_k ≥ log Π₁₅ = 40.96, Euclid's p ≤ Π_k + 1 gives log p ≤ X + 1, so
log(Π_{k+1} + 4) ≤ log p + log(Π_k + 4) and
  g(Π_{k+1})/g(Π_k) ≥ √p · X/(2X + 1) · log X/log(2X + 1) ≥ 3.0220 > 2·53/52 ≥ 2p/(p − 1) = T_{k+1}/T_k
(both factors increase with X and p; the value at p = 53, X = log Π₁₅ is computed in cert.py). With the base g(Π₁₅) > T₁₅ (k = 15 excluded)
this gives g(Π_k) > T_k for every k ≥ 15, so every k ≥ 15 is excluded.
**Result: every unit-generated real quadratic order with Cl(O) = Cl(O)[2] has Δ < D₀ = 20485358601780401.**
(The brief's route via Robin's ω(n) ≤ 1.3841 log n/log log n gives a cruder D₀; the k-split is sharper and needs no extra input.)

**3.5 Robustness.** Replacing 0.11 by c, cert.py bisects the smallest c (3 decimals) with D₀(c) ≤ 10¹⁸ (the sieve range): c = 0.037.
So a constant three times weaker than Littlewood78's suffices. For c = 0.037 the degree k = 15 is no longer excluded (its D₁₅⁺(c) ≤ 10¹⁸ enters D₀(c)) and the induction starts at k = 16; cert.py picks the base K for each c (the smallest K ≥ 15 with every k ∈ [K, 60] excluded), and the step constant computed at k = 15 bounds every later step (orchestrator reviewer note).

## 4. The computation (cert.py)

**4.1 Exact class groups for n ≤ 1000** (both families, 1998 orders; this covers every Δ < 100, where Lemma A is not claimed).
forms.py enumerates all reduced primitive forms (b ≤ ⌊√Δ⌋, b ≡ Δ mod 2, a | (Δ − b²)/4 with ⌊√Δ⌋ < 2|a| + b and 2|a| − b ≤ ⌊√Δ⌋;
these integer tests are equivalent to the real ones because √Δ is irrational), splits them into ρ-cycles (h⁺ = number of cycles,
each ρ-image is asserted to be reduced and cycles disjoint), and sets h = h⁺ or h⁺/2 according as (−1, b₀, ·) lies on the principal cycle.
For each class C = [F]: narrow test C² = 1, wide test C² ∈ {1, J}, J = [(−1, b₀, ·)] = kernel of Cl⁺ → Cl (the class of √Δ·O).
C² is computed by composition (Cohen Alg. 5.4.7) and reduction (ρ until reduced); the result is located by its cycle id.
Every test is repeated without composition: C⁻¹ = [(a, −b, c)] = [(c, b, a)] and J·C⁻¹ = [(−c, b, −a)]
(J·[(a, b, c)] = [(−a, b, −c)]: for 𝔟 = √Δ·𝔞 the basis (√Δα₁, −√Δα₂) is positively oriented and N(x√Δα₁ − y√Δα₂)/N(𝔟) = −f_𝔞(x, −y));
(c, b, a) and (−c, b, −a) are reduced when (a, b, c) is (2|c| = (√Δ − b)(√Δ + b)/(2|a|) lies in (√Δ − b, √Δ + b)). The two tests agree (asserted) on all 2080 orders.
Group structures (KL Tables 2A/2B columns) from the numbers of classes of order dividing 2^k.

**4.2 Sieve, n > 1000, Δ_n < 10¹⁸** (n ≤ 10⁹): by Cor. B a 2-torsion order must pass "(Δ/p) ≠ 1 for every prime p with 4p⁴ < Δ".
cert.py: chunks of 4·10⁶ values of n; per chunk the primes with 4p⁴ < Δ_{min} (a subset of each member's primes, so the chunk filter keeps a superset),
residue tables A_p = {r mod p : ((r² − 4N)/p) ≠ 1}; then every chunk survivor is re-tested exactly with its own full prime list.
Full survivors: 58 (N = +1, largest n = 5182) and 24 (N = −1, largest n = 5968), all with n ≤ 5968, i.e. Δ < 3.6·10⁷.
Check C5: an independent brute-force test (Jacobi-symbol code) over all n ∈ [1001, 3·10⁶] gives exactly the same survivor sets;
C9: an independent C program (`sievecheck.c`, residue tables by squaring, 128-bit bounds) over the full range n ∈ [1001, 10⁹] gives the same 58 + 24 survivors
(`logs/sievecheck_output.txt`).
**4.3 Survivors**: exact class groups (4.1) show none is 2-torsion (logs/survivors.txt: 82 orders, h between 21 and 192, e.g. n = 1022, N = +1: Cl ≅ ℤ/16 × (ℤ/2)²).
PARI quadclassunit (GRH-conditional, used only as a cross-check) agrees on h and on 2-torsion for all 2080 orders.

## 5. Proof of the Theorem
Let Cl(O_Δ) = Cl(O_Δ)[2]. By 3.4, Δ < D₀ < 10¹⁸. If n ≤ 1000, 4.1 decides directly. If n > 1000, then Δ > 10⁶ ≥ 100, so by Cor. B n is a sieve survivor (4.2),
and 4.3 shows no survivor is 2-torsion. Hence the 2-torsion orders are those found in 4.1: 53 + 37 entries (Δ = 5 twice), 89 orders, identical with KL Tables 2A/2B
including Cl, Cl⁺ and the conductors. (b): Cl⁺ = Cl⁺[2] implies Cl = Cl[2] (Cl⁺ ↠ Cl), and 4.1 computes Cl⁺ exactly. (c): §3 is not needed for Δ < 10¹⁸. (d): 3.5. ∎

## 6. Remarks
- KL's text (l.1609) says "a total of 86" orders; their tables contain 55 + 37 = 92 entries (with −4, −3 and Δ = 5 twice), i.e. 89 distinct real orders. Count erratum only.
- KL l.1743: "The complete classification ... seems out of reach at this time, in parallel with the idoneal number problem being unsolved."
  The obstruction is the same (an effective lower bound for L(1, χ)); (L) removes it in both problems.
- The list has 60 maximal orders (f = 1) and 29 non-maximal ones; one class per genus: 84 (56 maximal).
- For the maximal orders the GRH-conditional lists of Louboutin–Mollin–Williams (1993, "Class groups of exponent two in real quadratic fields":
  "with one possible exceptional value ruled out by the GRH") become unconditional for the n² ± 4 families.

NO DISCOVERY: OpenAI "Primitive roots for every admissible integer base" (4 Oct 2026): no false step found. Kummer-field part is uniform and correct; the result rests on the unrefereed marked Type II estimate of companion [21] ("The Poisson–Dirichlet law for prime predecessors", 24 Sep 2026, Thm 3.1).

Source: openai/math preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026 (92 pp). Text: `.scratch/primroot/p.txt`; [21]: `.scratch/primroot/pd.txt`. No Lean proof of either (CONTENTS.md entries 011, 029 have no Lean link; Wikipedia: "unreviewed, no Lean").

## How the argument goes
1. Prop 2.1 (§§10–12): ≫ x/L² primes p ∈ (x,2x), p ≡ u (mod M), p−1 = c·r·Q, c ∈ {2,4}, Q > x^0.9 prime, all primes of r in (e^{L^0.1}, e^{L^0.3}). Class u (Lemma 2.4) forces (a/p) = −1, so the index i_p(a) is odd.
2. An odd prime q | i_p(a) is either Q (then p | ∏_{j≤2x^0.1}(a^j−1): O(x^0.2/L) primes lost) or q | r, q ≤ e^{L^0.3}; then p splits completely in K_q = Q(μ_q, a^{1/q}).
3. Prop 2.2: #split ≪ x/(q(q−1)L) + x^{1−δ0}; union over q is o(x/L²). Hooley's medium range never occurs: Kummer degrees ≤ e^{2L^0.3} = x^{o(1)}, never x^{1/2}.

## Kummer uniformity and zero-free input: checked, sound
- Lemma 9.1: n_q = q(q−1), D_q ≤ q^{q(q−2)}(q^q|a|^{q−1})^{q−1}, log D_q + n_q ≪ n_q log 2q.
- Descent: K̃ = K_q(μ12) is abelian over K_q and over F_q = Q(μ_{12q}); ζ_K̃ = ∏ L_{F_q}(s,η) (Thm 1.2 applies), ζ_K̃/ζ_{K_q} entire.
- Lemma 9.2: Θ ≪_f x + x^{1−δ}(log D + n), constant independent of K (functional equation on Re s = −1/2, Hasanalizade–Shen–Wong uniform zero count, residue at s = 0 is O(n)). Each split p gives n_q ideals, so the degree cancels: error x^{1−δ0} log(2q)/L. Only an upper bound is used; no Chebotarev main term.
- Uniform enough: Thm 1.2 gives each cyclotomic F ⊇ μ12 a zero-free half-plane of fixed width 10^-6, exactly what Prop 2.2 needs. Liu's 0.874887 would also suffice. OpenAI's 7/8 alone would not (ζ, Dirichlet L, Hecke over Q(√−3) only; Prop 2.2 needs Hecke L over Q(μ_{12q}) for every q).
- Width traced through §§4, 7: M0 = 1033/1000, I_Z ≪ Z^{933/2000+ε} = Z^{7/15−1/6000+ε} (exact). Exponents are numbers; F enters constants only (incl. d_F in Lemma 7.1). §§3–8 and App. A (cubic theta) not line-checked.

## Where the weight sits
- Prop 2.1 alone gives ≫ x/L² primes with P⁺(p−1) > p^0.89: far past the unconditional record 0.677 (Baker–Harman). The full Poisson–Dirichlet law for p−1 was known only under Elliott–Halberstam (Bharadwaj–Rodgers). So Thm 1.1 is exactly as solid as [21, Thm 3.1].
- Lemma 10.2 is not a black-box citation: it argues the proof of [21, Thm 3.1] works for other coefficients. Checked against [21]; the extension holds. [21] §5 (major term) uses only Lemma 5.2 (= hypothesis (10.6)) plus divisor mean squares; §3.3, §4, §4.9 use only sup/ℓ² norms and roughness.
- [21] core spot-checked, no error: Lemma 3.3 (Kloosterman 4th moment, det-2 map); (4.27) lattice box from goodness test (i); (4.33) clean edge ‖T(ϑ)‖²/Y ≪ 1/d + Y^-0.1 log Y + d log d/Y with d > L^{A0}; (4.36) phase separation forces a unique omission tuple (mass M^-K); (4.44)–(4.55) rank expansion.
- Trace bound (4.1) is not contradicted by backtracking walks (symmetrization costs M^-K per step). Character-mod-R test coefficients cannot break the coefficient-independent replacement (4.2): each costs ≤ 1/R ≤ L^-A0.

## Elementary steps verified (`.scratch/primroot/table.py`)
- Lemma 2.4 table: (d0/p) = −1 for every prime p < 2·10^5 in each listed class; a = 2, class 5 mod 8 too.
- B–V level (κ−ε)(1+ε)/4 > 0; sieve margin 0.9 − 0.2 − 1/4 > 1/4; Buchstab identity (11.6); Lemma 12.5 densities (3j−2)/j², (j−1)/j²; exponent a/2 + h0/6 − 1 = −8/15.
- Cosmetic only: Thm 10.1 misquotes [21]'s title.

## Residual risk (unverified; any failure kills Thm 1.1)
- [21] Prop 4.1 signed-memory trace bound not checked line by line; [21] Lemma 3.4 error budget unchecked; Thm 1.2 §§5–8 unchecked unless Liu is accepted.
- Next: line-check [21] §4.6 (reversal table) and §4.8 (ghost and birth bookkeeping).

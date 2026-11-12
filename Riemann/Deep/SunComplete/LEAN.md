FORMALIZED: Theorem A (Sun–Hou Conjecture 3.5(i): every finite field with q > 7 elements has a numbering a_1, …, a_q whose cyclic adjacent sums are all primitive) as `SunLean.theoremA_of` and `SunLean.theoremA_of_mid2`, from the named hypotheses `ChvatalErdos` (Chvátal–Erdős 1972), `Lemma3A` (NOTE Lemma 3, Cohen's sieve with Gauss sums) and `Lemma6Amid2` (NOTE Lemma 6 only for 30,140 ≤ q < 118,380,036,096 with ω(q − 1) ≤ 11 and φ(q − 1) < 7·2^ω(q−1)·(⌊√q⌋ + 1); Lemma 6 itself is proved for every larger q, every q with ω(q − 1) ≥ 12 and every q with φ(q − 1) ≥ 7·2^ω(q−1)·(⌊√q⌋ + 1): `lemma6A_big2`, `csa_noSieve_phi`), the certificate hypotheses `Prop7A` (verify.py per-q sieve data, 7 < q < 30,140) and `CertA` (certA.json, q ∈ E_A), and `FFacts F` (three Mathlib finite-field facts: `card_prim`, in raw form, compiled in FFCoreC.lean; `primePow`, `two_eq_zero` and the repackaging as `FFacts F` are proved in FFInst.lean, which does not compile under the 300 MB cap); Proposition 4 and Lemma 5 are proved in full; kernel-checked by `decide`: q = 3, 5, 7 fail for every field, q = 11, 13 work for every field, q = 4, 9 work in explicit models.

# sunlean — Lean 4 formalization of Theorem A (Sun–Hou Conjecture 3.5(i))

2026-11-12, agent sunlean. Side result: discoveries/sun-hou.md, Theorem A of NOTE.md. Lean v4.35.0-rc2, Mathlib v4.35.0-rc2. Compiled modules are sorry-free, use only the standard axioms, and contain no `native_decide`, `skipKernelTC` or `axiom` declarations.

## Files (`Riemann/Deep/SunComplete/Lean/`)

| file | content | compiled (memcap peak) |
|---|---|---|
| Graph.lean | `HasHamCycle`, `IndepSet`, `ConnDel` (G − S nonempty and connected), `ChvatalErdos` (statement), `finRotate` lemmas | yes |
| Defs.lean | `IsPrim`, `SunHou`, `sumGraph` (= H_q), `sunHou_iff`, `not_isPrim_zero`, `two_nbrs`, `primSet` | yes |
| Arith.lean | the ℚ inequalities behind Prop 4 (i), (ii), the F-analysis and Lemma 5 | yes |
| Prop4.lean | `L3Data`, `indep_bound`, `alpha_le`, `cut_bound`, `ham_of_bounds`, **`prop4`**, **`lemma5`** | yes |
| Sieve.lean | sieve data (θ, δ, Δ, W, τ, E), `pairCount`, `Lemma3A`, `CritA`, `CSA`, `EA`, `Prop7A`, `Lemma6A`, `CertA`, `FFacts` | yes |
| CaseA1.lean | τ > 0, E ≥ 0; e(I, I) ≤ ε·|I|, e(A, B) = 0, every degree of H_q ≥ φ(q − 1) − ε | yes |
| TheoremA.lean | `l3dataA`, `ham_of_critA` (Prop 4 for H_q), `ham_of_CSA` (Lemma 5 for H_q), **`theoremA_of`** | yes |
| L6Arith.lean | the ℚ inequalities of the no-sieve case of Lemma 6 (`noSieve_core`, `noSieve_arith`) | yes |
| L6NT.lean | Euler product φ(m) = m·∏(1 − 1/p); (ω + 1)·∏(1 − 1/p) ≥ 1; (ω + 1)! ≤ ∏p ≤ m | yes |
| L6.lean | `noSieve`, `csa_noSieve`, `bound_fact`, **`lemma6A_big`**, `Lemma6Amid`, **`lemma6A_of_mid`** | yes |
| L6b.lean | `oddProd_le_prod` (k distinct primes: product ≥ 2·3·5⋯(2k − 1)), `bound_odd`, **`lemma6A_big2`**, `Lemma6Amid2`, **`lemma6A_of_mid2`** | yes |
| TheoremA2.lean | **`theoremA_of_mid`**, **`theoremA_of_mid2`** (Theorem A with `Lemma6Amid`, `Lemma6Amid2` in place of `Lemma6A`) | yes |
| Small.lean | transport along ring isomorphisms, `isPrim_zmod_iff`, **`not_sunHou_small`** (q = 3, 5, 7) | yes |
| Cycles.lean | the field `F4`, **`sunHou_F4`**; **`sunHou_zmod11`**, **`sunHou_zmod13`**, `sunHou_11_13` | yes |
| Cycles9.lean | the field `F9` = F_3[x]/(x² + x + 2), **`sunHou_F9`** (certA.json cycle) | yes |
| CheckA.lean, CheckS.lean | `#print axioms` | yes |
| FFCoreC.lean | `card_prim` from Mathlib (light imports) | yes |
| FFCoreP.lean | `primePow`, `two_eq_zero` from `Mathlib.FieldTheory.Finite.Basic` | no: killed at 309 MB |
| FFInst.lean | `ffacts : FFacts F` for every finite field | no: killed at 312 MB |

Peaks per module: see Build. FFInst.lean and FFCoreP.lean are the uncompiled instance files the brief allows: they name the Mathlib theorems (`FiniteField.isPrimePow_card`, `FiniteField.even_card_iff_char_two`, `ringChar.spec`, `IsCyclic.card_orderOf_eq_totient`, `Fintype.card_units`, `orderOf_units`); `theoremA_of` keeps `FFacts F` as an explicit hypothesis.
## Main statements (verbatim)

```lean
def IsPrim (x : F) : Prop := orderOf x = Fintype.card F - 1
def SunHou : Prop := ∃ a : Fin (Fintype.card F) ≃ F, ∀ i, IsPrim (a i + a (finRotate _ i))
def sumGraph : SimpleGraph F where
  Adj x y := x ≠ y ∧ IsPrim (x + y)
def HasHamCycle (G : SimpleGraph V) : Prop :=
  ∃ a : Fin (Fintype.card V) ≃ V, ∀ i, G.Adj (a i) (a (finRotate _ i))
theorem sunHou_iff (h2 : 2 ≤ Fintype.card F) : SunHou F ↔ HasHamCycle (sumGraph F)

theorem theoremA_of (hCE : ChvatalErdos) (hL3 : Lemma3A) (hP7 : Prop7A) (hL6 : Lemma6A)
    (hCert : CertA) (F : Type) [Field F] [Fintype F] (hF : FFacts F)
    (h7 : 7 < Fintype.card F) : SunHou F
theorem theoremA_of_mid (hCE : ChvatalErdos) (hL3 : Lemma3A) (hP7 : Prop7A)
    (hL6 : Lemma6Amid) (hCert : CertA) (F : Type) [Field F] [Fintype F] (hF : FFacts F)
    (h7 : 7 < Fintype.card F) : SunHou F
theorem lemma6A_big (q : ℕ) (hq : 3 ≤ q)
    (h : 53876069761024 ≤ q ∨ 16 ≤ (primeDivisors (q - 1)).card) :
    ∃ (d : SieveData q) (s : ℚ), CSA q d s
theorem lemma6A_of_mid (h : Lemma6Amid) : Lemma6A
theorem theoremA_of_mid2 (hCE : ChvatalErdos) (hL3 : Lemma3A) (hP7 : Prop7A)
    (hL6 : Lemma6Amid2) (hCert : CertA) (F : Type) [Field F] [Fintype F] (hF : FFacts F)
    (h7 : 7 < Fintype.card F) : SunHou F
theorem lemma6A_big2 (q : ℕ) (hq : 3 ≤ q)
    (h : 118380036096 ≤ q ∨ 12 ≤ (primeDivisors (q - 1)).card) :
    ∃ (d : SieveData q) (s : ℚ), CSA q d s
theorem lemma6A_of_mid2 (h : Lemma6Amid2) : Lemma6A

theorem prop4 (hCE : ChvatalErdos) (D : L3Data G) {δ₀ : ℚ} (hdeg : ∀ v, δ₀ ≤ (deg G v : ℚ))
    (hn : 3 ≤ Fintype.card V) (hNE : D.E < D.Nt) (C : ℕ)
    (hC : (Fintype.card V : ℚ) * (D.ε / D.τ + D.E) / (D.Nt + D.E) < C + 2)
    (ha1 : 1 ≤ δ₀ - C + 1)
    (hF : (Fintype.card V : ℚ) / 2 < δ₀ - C + 1 ∨
      (0 < Fq D.Nt D.E (Fintype.card V) C (δ₀ - C + 1) ∧
        0 < Fq D.Nt D.E (Fintype.card V) C ((Fintype.card V : ℚ) / 2))) :
    HasHamCycle G
theorem lemma5 (hCE : ChvatalErdos) (D : L3Data G) {δ₀ : ℚ} (hdeg : ∀ v, δ₀ ≤ (deg G v : ℚ))
    (hn : 3 ≤ Fintype.card V) (hNt : 0 ≤ D.Nt) {A' Dd : ℚ}
    (hA' : (Fintype.card V : ℚ) * (D.ε / D.τ + D.E) / (D.Nt + D.E) ≤ A')
    (hD : Dd ≤ δ₀ + 1) (hDA : A' < Dd) (h2A : 2 * A' < Fintype.card V)
    (hCS : D.E ^ 2 * Dd < (Dd - A') * D.Nt ^ 2 * (1 - 2 * A' / Fintype.card V)) :
    HasHamCycle G
def Fq (Nt E n C a : ℚ) : ℚ := Nt ^ 2 * a * (n - a - C) - E ^ 2 * (a + C) * (n - a)

theorem not_sunHou_small (F : Type) [Field F] [Fintype F]
    (h : Fintype.card F = 3 ∨ Fintype.card F = 5 ∨ Fintype.card F = 7) : ¬ SunHou F
theorem sunHou_F4 : SunHou F4
theorem sunHou_F9 : SunHou F9
theorem sunHou_11_13 (F : Type) [Field F] [Fintype F]
    (h : Fintype.card F = 11 ∨ Fintype.card F = 13) : SunHou F
```

## Hypotheses (named and certificate) and their sources

```lean
def ChvatalErdos : Prop :=
  ∀ (V : Type) [Fintype V] [DecidableEq V] (G : SimpleGraph V) (k : ℕ),
    3 ≤ Fintype.card V →
    (∀ I : Finset V, IndepSet G I → I.card ≤ k) →
    (∀ S : Finset V, S.card < k → ConnDel G S) →
    HasHamCycle G
def Lemma3A : Prop :=
  ∀ (F : Type) [Field F] [Fintype F] [DecidableEq F] (d : SieveData (Fintype.card F)) (s : ℚ),
    0 < d.δ → 0 ≤ s → (Fintype.card F : ℚ) ≤ s ^ 2 → ∀ X Y : Finset F,
      d.τ * (((Fintype.card F : ℚ) - 1) * X.card * Y.card) -
          (Fintype.card F : ℚ) * (pairCount X Y : ℚ) ≤ 0 ∨
      (d.τ * (((Fintype.card F : ℚ) - 1) * X.card * Y.card) -
          (Fintype.card F : ℚ) * (pairCount X Y : ℚ)) ^ 2 ≤
        (d.τ * d.E s) ^ 2 * ((X.card * ((Fintype.card F : ℚ) - X.card)) *
          (Y.card * ((Fintype.card F : ℚ) - Y.card)))
def Prop7A : Prop :=
  ∀ q : ℕ, PrimePow q → 7 < q → q < 30140 → q ∉ EA →
    ∃ (d : SieveData q) (s : ℚ) (C : ℕ), CritA q d s C
def Lemma6A : Prop := ∀ q : ℕ, PrimePow q → 30140 ≤ q → ∃ (d : SieveData q) (s : ℚ), CSA q d s
def Lemma6Amid : Prop :=
  ∀ q : ℕ, PrimePow q → 30140 ≤ q → q < 53876069761024 → (primeDivisors (q - 1)).card ≤ 15 →
    phi (q - 1) < 7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) →
    ∃ (d : SieveData q) (s : ℚ), CSA q d s
def Lemma6Amid2 : Prop :=
  ∀ q : ℕ, PrimePow q → 30140 ≤ q → q < 118380036096 → (primeDivisors (q - 1)).card ≤ 11 →
    phi (q - 1) < 7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) →
    ∃ (d : SieveData q) (s : ℚ), CSA q d s
def CertA : Prop :=
  ∀ (F : Type) [Field F] [Fintype F], Fintype.card F ∈ EA → 7 < Fintype.card F → SunHou F
structure FFacts (F : Type) [Field F] [Fintype F] : Prop where
  primePow : PrimePow (Fintype.card F)
  card_prim : (primSet F).card = phi (Fintype.card F - 1)
  two_eq_zero : Fintype.card F % 2 = 0 → (2 : F) = 0
```

- `CritA q d s C`: `0 < δ`, `0 ≤ s`, `q ≤ s²`, `E(s) < q − 1`, `α* < C + 2`, `1 ≤ δ₀ − C + 1`, and `q/2 < δ₀ − C + 1` or `F(δ₀ − C + 1) > 0 ∧ F(q/2) > 0`, with `α* = q(ε/τ + E)/(q − 1 + E)`, `δ₀ = φ(q − 1) − ε`, `ε = [q odd]` (Sieve.lean). This is verify.py `crit_exact('A', …)` with `s = sqrt_ub(q)` and `C = ⌊α*⌋ − 1`.
- `CSA q d s`: `0 < δ`, `0 ≤ s`, `q ≤ s²` and some `A' ≥ α*`, `D ≤ δ₀ + 1` with `A' < D`, `2A' < q`, `E²D < (D − A')(q − 1)²(1 − 2A'/q)` (Lemma 5's (CS)).
- Sieve data: `core ∪ siev = primeDivisors (q − 1)`, disjoint; `θ = ∏_core (1 − 1/p)`, `δ = 1 − Σ_siev 1/p`, `Δ = (s − 1)/δ + 2` (`s ≥ 1` sieving primes) or `1`, `W = 2^|core|`, `τ = θδ`, `E(s) = (WΔ − 1)s + 1` (NOTE §2). `pairCount X Y = #{(x, y) ∈ X × Y : x + y ∈ P}`.

| hypothesis | kind | source | what would discharge it |
|---|---|---|---|
| `ChvatalErdos` | named theorem | V. Chvátal, P. Erdős, Discrete Math. 2 (1972) 111–113 | a Lean proof of Chvátal–Erdős (longest-cycle argument; not in Mathlib) |
| `Lemma3A` | named theorem | NOTE Lemmas 1–3 (Cohen's prime sieve, Gauss sums of absolute value √q, Parseval, Cauchy–Schwarz) | the sieve and Lemmas 1–3 in Lean (Mathlib has `gaussSum` and `gaussSum_mul_gaussSum_eq_card`, i.e. |g|² = q, but not the sieve) |
| `Lemma6Amid2` (replaces `Lemma6A` in `theoremA_of_mid2`; `Lemma6Amid`, the same with q < 53,876,069,761,024 and ω ≤ 15, in `theoremA_of_mid`) | named theorem, the rest of Lemma 6 (q − 1 very smooth: φ(q − 1) < 7·2^ω(⌊√q⌋ + 1)) | NOTE Lemma 6 for ω(q − 1) ≤ 11 (cores of the t = 0, 1, 2 smallest primes; table points q₀ ≤ 30,138 for r ≤ 6); verify.py `check_lemma6` | the sieve case in Lean: j-th prime factor ≥ j-th prime, θ ≥ α_r, δ ≥ δ_{t,r}, the table for r ≤ 15 and monotonicity in q |
| `Prop7A` | certificate | verify.py `check_prop4` / `crit_best('A', q)` (all 2^ω(q−1) partitions, `Fraction`) for the 3,256 prime powers 8 ≤ q < 30,140 outside E_A | per-q kernel checks (too large for `decide` here) |
| `CertA` | certificate | certA.json (She x for 62 primes, explicit cycles for q = 9, 13, 16, 25, 27, 29, 49, 61, 121, 169, 361, 841, 1681); verify.py `check_certificates` (gen_certs.py generates) | Lemma 8 in Lean + per-field `decide` + uniqueness of F_q; done here: q = 11, 13 (every field), q = 9 (one model) |
| `FFacts F` | Mathlib facts | `FiniteField.isPrimePow_card`; `IsCyclic.card_orderOf_eq_totient` on `Fˣ` with `Fintype.card_units`, `orderOf_units`; `FiniteField.even_card_iff_char_two`, `ringChar.spec` | FFCoreC.lean compiles `card_prim` (raw form, 283 MB); FFCoreP.lean and FFInst.lean are over the cap |

`Lemma3A` is Lemma 3 in squared form: with `L = τÑ|X||Y| − n·e(X,Y)` and `Z = |X|(n−|X|)|Y|(n−|Y|) ≥ 0`, `e ≥ (τ/n)(Ñ|X||Y| − E√Z)` ⟺ `L ≤ τE√Z` ⟺ `L ≤ 0 ∨ L² ≤ (τE)²Z` (as `τE ≥ 0`). It is stated for every rational `s ≥ √q`: Lemma 3 holds with any larger `E`, and `WΔ − 1 ≥ 0` makes `E(s)` increasing in `s`.

## What is proved, and how

- **Proposition 4 and Lemma 5, in full, for any finite graph** (Prop4.lean, Arith.lean). Input `L3Data G`: a pair count `e` with Lemma 3 in squared form, `e(I,I) ≤ ε|I|` on independent sets, `e(A,B) = 0` across an empty cut; and a minimum-degree bound `δ₀`.
  - (i) `indep_bound`, `alpha_le`: |I| ≤ α* = n(ε/τ + E)/(Ñ + E). For `X = Y = I` the root is exact, `√Z = s(n − s)`, so the squared form gives `L ≤ τE·s(n − s)` (`indep_arith`).
  - (ii) `cut_bound`: an empty cut (A, B) of G − S with |A| ≤ |B| gives `δ₀ ≤ a − 1 + c` (a vertex of A has all neighbours in (A ∖ v) ∪ S), `2a + c ≤ n` and `Φ(a, c) = Ñ²a(n − a − c) − E²(a + c)(n − a) ≤ 0` (`cut_arith`: divide `(τÑab)² ≤ (τE)²a(b + c)b(a + c)` by `τ²ab`).
  - `ham_of_bounds`: Chvátal–Erdős with `k` once `α ≤ k` and `Φ(a, c) > 0` for every admissible `(a, c)` with `c < k`.
  - `prop4`: `k = C + 1`, any natural `C` with `α* < C + 2` (NOTE: `C = ⌊α*⌋ − 1`); `C + 1 ≤ n` from `a₁ ≥ 1` and `deg ≤ n − 1`; `Φ(a, c) − F(a) = (C − c)(Ñ²a + E²(n − a)) ≥ 0`; F concave: `F(x) − F(y) = (Ñ² − E²)(x − y)(n − C − x − y)`, so `F(a) ≥ F(a₁)` or `F(a) ≥ F(n/2)` on `[a₁, n/2]` (`prop4_arith`).
  - `lemma5`: `k = ⌈A'⌉₊`; `n·Φ(a, c) − (n − a)[a(Ñ²(n − 2c) − E²n) − E²nc] = Ñ²ac(n − 2a) ≥ 0`, and (CS) makes the bracket positive for `a ≥ D − c`, `c < A'` (`lemma5_arith`).
- **Case (A)** (CaseA1.lean, TheoremA.lean): `e = pairCount`, `n = q`, `Ñ = q − 1`, `ε = [q odd]`. Off-diagonal pairs of an independent set are non-edges; for even q, `x + x = 2x = 0` is not primitive (`FFacts.two_eq_zero`). Every degree is `≥ φ(q − 1) − ε`: `y ↦ x + y` is a bijection and `|P| = φ(q − 1)` (`FFacts.card_prim`). τ > 0 (primes ≥ 2, δ > 0), E ≥ 0 (WΔ ≥ 1).
- **`theoremA_of`**: q ∈ E_A → `CertA`; q < 30,140 → `Prop7A` → `ham_of_critA`; q ≥ 30,140 → `Lemma6A` → `ham_of_CSA`; then `sunHou_iff`. `theoremA_of_mid` is `theoremA_of` composed with `lemma6A_of_mid`.
- **Lemma 6 outside a finite middle range** (L6Arith, L6NT, L6). No sieve: core = all r prime divisors of q − 1 (δ = Δ = 1, W = 2^r, τ = θ ≥ 1/(r + 1), as the prime divisors are r distinct integers ≥ 2), s = ⌊√q⌋ + 1, A' = 1/τ + E, D = φ(q − 1) = (q − 1)θ (Mathlib's Euler product). With K = (r + 1)2^r·s, (CS) holds once 7K ≤ q − 1, since then K²(n + 1) < n(n − 2K)(n − 4K), n = q − 1 (`noSieve_core`); and 7K ≤ q − 1 follows from q ≥ 196(r + 1)²4^r (`csa_noSieve`). For r ≥ 16, q − 1 ≥ ∏p ≥ (r + 1)! ≥ 196(r + 1)²4^r (`fact_le_prod`, `bound_fact`); for r ≤ 15 the bound holds once q ≥ 196·16²·4^15 = 53,876,069,761,024. Hence `lemma6A_big`. With r₁ = 1/θ in place of r + 1 the same algebra gives (CS) whenever φ(q − 1) ≥ 7·2^r·(⌊√q⌋ + 1) (`csa_noSieve_phi`). Together: `lemma6A_of_mid : Lemma6Amid → Lemma6A`. Sharper (L6b): k distinct primes have product ≥ 2·3·5⋯(2k − 1), as at most one is even (`max_prime_bound`, `oddProd_le_prod`), which is ≥ 196(r + 1)²4^r from r = 12 on (`bound_odd`); so `lemma6A_big2` (q ≥ 196·12²·4^11 = 118,380,036,096 or ω(q − 1) ≥ 12) and `lemma6A_of_mid2 : Lemma6Amid2 → Lemma6A`. This replaces NOTE's tail argument (r ≥ 61) by a cruder one valid from r ≥ 16 (r ≥ 12 with L6b).
- **Sharpness and small fields, by kernel `decide`** (Small.lean, Cycles.lean, Cycles9.lean): a Sun–Hou numbering gives every element two distinct neighbours (`two_nbrs`), but in `ZMod 3, 5, 7` the vertices `1, 1, 5` have at most one; `ZMod.ringEquivOfPrime` and `sunHou_of_ringEquiv` carry this to every field with 3, 5 or 7 elements. `F4` and `F9` are explicit fields (axioms by `decide`; F9's associativity and distributivity by `ring` on `ZMod 3` coordinates); `(0, ω, 1, ω²)` and certA.json's q = 9 cycle `[3, 0, 4, 5, 2, 6, 7, 8, 1]` work. `ZMod 11` (She's Lemma 8 with x = 6: `(4, 2, 6, 0, 8, 9, 10, 7, 1, 5, 3)`) and `ZMod 13` (certA.json cycle) work, hence every field with 11 or 13 elements.
- Primitivity tests used by `decide`: `isPrim_zmod_iff`, `F4.isPrim_iff`, `F9.isPrim_iff`, all `orderOf_eq_iff` (order exactly q − 1).

## Not proved (left as hypotheses)

- `Lemma6Amid2`: Lemma 6 for the prime powers 30,140 ≤ q < 118,380,036,096 with ω(q − 1) ≤ 11 and φ(q − 1) < 7·2^ω(q−1)·(⌊√q⌋ + 1), i.e. q − 1 very smooth. There the sieve is needed (without it (CS) needs q of order 10^6 already for r = 6). NOTE's proof: the core is the t smallest prime divisors, the j-th prime divisor is at least the j-th prime, so θ ≥ α_r, τ ≥ α_tδ_{t,r} and q − 1 ≥ P_r; (CS) at the tabulated points (r ≤ 11 here) in exact arithmetic; monotonicity in q. Not done in this budget. Lemma 5, which turns (CS) into a Hamilton cycle, is proved.
- `ChvatalErdos` (not in Mathlib), `Lemma3A` (Gauss sums and the sieve), `Prop7A` (3,256 per-q checks), `CertA` (62 She primes and 13 explicit cycles up to q = 1681; only q = 9, 11, 13 are kernel-checked here, and transport to every field of order q = p^k, k ≥ 2, needs the uniqueness of F_q).
- `FFacts.primePow`, `FFacts.two_eq_zero`: Mathlib proves them (FFInst.lean), but importing `Mathlib.FieldTheory.Finite.Basic` exceeds the cap.
- `sunHou_F4`, `sunHou_F9` concern one model each; transport to every field with 4 or 9 elements needs `FiniteField.ringEquivOfCardEq` (GaloisField, not imported).

## Fidelity

The conjecture as printed (NOTE.md §0, Theorem A; Z.-W. Sun, Nanjing Univ. J. Math. Biquarterly 36 (2019), Conjecture 3.5(i), with Q.-H. Hou): *for every prime power q > 7, the elements of F_q can be numbered a_1, …, a_q so that all q sums a_1 + a_2, …, a_{q−1} + a_q, a_q + a_1 are primitive*; false for q = 3, 5, 7, true for q = 4. She, arXiv:2609.19114 (e-print TeX, lines 111–117): *"Write P_q for the set of primitive elements of F_q, that is, the generators of the multiplicative group F_q^* = F_q∖{0}. Define the graph H_q on F_q by xy ∈ E(H_q) ⟺ x ≠ y and x + y ∈ P_q."* Line 131: *"Sun and Hou conjectured that these three graphs are Hamiltonian, apart from a finite list of small fields."* Theorem 1.2(1), line 136: *"If q > 18 888 871 is a prime power, then H_q is Hamiltonian."*

| concept | Lean | agreement |
|---|---|---|
| F_q | any `F : Type` with `Field F`, `Fintype F`; `q = Fintype.card F` | every finite field; q is a prime power by Mathlib (`FFacts.primePow`) |
| primitive | `IsPrim x := orderOf x = q − 1` | F_q^* is cyclic of order q − 1, so order q − 1 ⟺ generator; `orderOf 0 = 0` excludes 0 (`not_isPrim_zero`) |
| numbering a_1, …, a_q | `a : Fin q ≃ F`, a_{i+1} = `a i` | a bijection onto F_q |
| cyclic sums | `a i + a (finRotate q i)`, `finRotate q i = i + 1 mod q` | includes a_q + a_1 |
| H_q | `sumGraph F`: `x ≠ y ∧ IsPrim (x + y)` | She's H_q verbatim |
| Hamiltonian | `HasHamCycle` | `sunHou_iff`: equals the Sun–Hou property for q ≥ 2; degenerate only for n ≤ 2 |
| κ ≥ α ⇒ Hamiltonian | `ChvatalErdos` with an integer k | implied by Chvátal–Erdős: k < n gives κ ≥ k ≥ α; k = n forces G complete; k > n is vacuous (S = univ) |
| Lemma 3 | `Lemma3A` | exact squared form; quantified over rational s ≥ √q with E(s) ≥ E(√q), so implied by NOTE's Lemma 3 |
| Prop 4 | `prop4` | any C with α* < C + 2 (NOTE's C = ⌊α*⌋ − 1 is one); otherwise verbatim (a₁, F, the two alternatives) |
| Lemma 5 / (CS) | `lemma5`, `CSA` | A' ≥ α*, D ≤ δ₀ + 1 free (NOTE's A' = ε/τ + E, D = δ₀ + 1 and Lemma 6's worst-case values are instances) |
| Prop 7 | `Prop7A` | verify.py `crit_exact('A', …)` condition for each q ∉ E_A, 7 < q < 30,140 |
| E_A | `EA` | the 79 values of NOTE §4 / verify.py `EA`, copied |
| sharpness | `not_sunHou_small`, `sunHou_F4` | q = 3, 5, 7 fail for every field; q = 4 works in the model `F4` |

**Reviewer** (statement fidelity, one `reviewer` subagent; the definitions, main statements, hypotheses and the conjecture pasted inline): verdict **ISSUES**, all documentation-level, "no identified mismatch affecting the conjecture for q > 7". Points: (1) definitions faithful; (2) `HasHamCycle` degenerates for n ≤ 2 (harmless: n ≥ 3 throughout) — now stated in its docstring; (3) `ConnDel` and `ChvatalErdos` correct; the k > n case is vacuous — now stated in the docstring; (4) `Lemma3A` is exactly Lemma 3 with E(s), not stronger; (5) `Lemma6A`'s docstring named NOTE's core choice that the definition does not impose — docstring corrected; `CertA` correctly quantifies over every field of each cardinality; `FFacts` satisfiable; (6) `theoremA_of` is conditional on the listed hypotheses; `sunHou_F4` covers one model of F_4. All addressed above; no statement changed.

## Build

Clean rebuild, one module at a time, plain `lean -j1` (no lake process), oleans to `.scratch/sunlean/olean` (build output, deleted after the run), every compile under the watchdog (`Riemann/Deep/SunComplete/Lean/audit/verify.sh`; log `verify.log`, peaks `peaks.log`):

```
LP=.scratch/sunlean/olean:$(cat Riemann/Deep/DPCUncond/Lean/audit/leanpath.txt)      # lake's LEAN_PATH
for m in Graph Arith L6Arith Prop4 Defs Sieve L6NT L6 L6b CaseA1 TheoremA TheoremA2 Small Cycles Cycles9 CheckA CheckS FFCoreC; do
  d=Riemann/Deep/SunComplete/Lean; mkdir -p .scratch/sunlean/olean/$d
  LEAN_PATH=$LP OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 /usr/bin/time -l \
    ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 -o .scratch/sunlean/olean/$d/$m.olean $d/$m.lean
done
```

Clean rebuild 2026-11-12 (second run; module end times 19:31:51–19:44:13 UTC), all rc = 0; memcap peak = top MEM, the cap's measure; footprint = `/usr/bin/time -l`. Max peak 283 MB:

| module | s | memcap peak | footprint |
|---|---|---|---|
| Graph | 138 | 176 MB | 222 MiB |
| Arith | 43 | 247 MB | 287 MiB |
| Prop4 | 43 | 253 MB | 263 MiB |
| Defs | 58 | 239 MB | 252 MiB |
| Sieve | 52 | 259 MB | 269 MiB |
| CaseA1 | 58 | 253 MB | 269 MiB |
| TheoremA | 160 | 250 MB | 262 MiB |
| Small | 33 | 242 MB | 258 MiB |
| Cycles | 62 | 254 MB | 267 MiB |
| Cycles9 | 101 | 259 MB | 272 MiB |
| CheckA | 50 | 252 MB | 260 MiB |
| CheckS | 39 | 213 MB | 250 MiB |
| FFCoreC | 40 | 283 MB | 295 MiB |

Follow-up (Lemma 6), 2026-11-12 20:02–20:07 UTC: dependencies rebuilt (Arith 286 MB, Graph 211, Defs 233, Sieve 257, Prop4 240, CaseA1 252, TheoremA 255), then the new modules, all rc = 0 (memcap peak / footprint): L6Arith 272 MB / 278 MiB, L6NT 256 MB / 272 MiB, L6 277 MB / 295 MiB (final version, with `csa_noSieve_phi`), L6b 257 MB / 279 MiB, TheoremA2 (final, both theorems) 226 MB / 264 MiB.

Killed by the watchdog (not part of the build): a first, unsplit Prop4 with `nlinarith` (294 MB), FFInst (312 MB), FFCore and FFCoreP (309 MB each; `Mathlib.FieldTheory.Finite.Basic` alone sits near the cap), a first Cycles9 with 729-case `decide` field axioms (316 MB).

## Static scan

`rg -w 'sorry|admit|native_decide|skipKernelTC' Riemann/Deep/SunComplete/Lean/`: no match. `rg '^\s*axiom\s' Riemann/Deep/SunComplete/Lean/`: no match. 20 files, 1,978 lines.

## `#print axioms`

All 25 checked declarations depend on exactly `[propext, Classical.choice, Quot.sound]` (CheckA.lean, CheckS.lean, FFCoreC.lean, L6.lean, L6b.lean, TheoremA2.lean):

```
'SunLean.theoremA_of' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.L3Data.prop4' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.L3Data.lemma5' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.L3Data.ham_of_bounds' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.ham_of_critA' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.ham_of_CSA' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.sunHou_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.deg_ge' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.pairCount_diag' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.not_sunHou_small' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.not_sunHou_zmod7' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.sunHou_F4' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.sunHou_F9' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.sunHou_zmod11' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.sunHou_zmod13' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.sunHou_11_13' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.sunHou_of_ringEquiv' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.FFCore.card_prim' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.lemma6A_big' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.csa_noSieve_phi' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.lemma6A_of_mid' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.theoremA_of_mid' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.lemma6A_big2' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.lemma6A_of_mid2' depends on axioms: [propext, Classical.choice, Quot.sound]
'SunLean.theoremA_of_mid2' depends on axioms: [propext, Classical.choice, Quot.sound]
```

## Repository state

The repo-root `.lake` (a real directory: the Mathlib cache of the Lean agents) was left as it is: no `lake build`, no `lake exe cache get`, the placeholder not restored (the orchestrator restores it). No commits.

## Entry for FORMALIZATIONS.md

- Sun–Hou Conjecture 3.5(i) (Theorem A, discoveries/sun-hou.md): Lean `theoremA_of_mid2` from Chvátal–Erdős, Lemma 3, Lemma 6 on 30,140 ≤ q < 1.2·10^11 with ω(q−1) ≤ 11 and φ(q−1) < 7·2^ω(⌊√q⌋+1) (named), Prop 7 data and certA (certificates), three Mathlib facts; Prop 4, Lemma 5 and the rest of Lemma 6 proved; q = 3, 5, 7 fail, q = 4, 9, 11, 13 work (`decide`) — Riemann/Deep/SunComplete/LEAN.md

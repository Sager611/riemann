# serrenumK — DERIVATION (every replaced step and constant)

Source: arXiv 2611.19556 (TeX Finalv2.tex, fetched 19 Nov 2026; line numbers below). ρ := log N_{K/Q}𝔪,
𝔪 := rad(2O_K𝔑_{E₁}𝔑_{E₂}) (Thm A) or rad(2O_K𝔑_E) (Thm B). n_K = [K:Q], d_K = |disc K|. log = natural log.

## 0. Inputs
- **chebexp (A)** (Riemann/Deep/ChebExp/REPORT.md, audit Riemann/Deep/ChebExpAudit/REPORT.md: PASSED). For every Galois
  L/K with L ≠ Q and every conjugacy class C of Gal(L/K) there is a prime 𝔭 of K, unramified in L, with N𝔭 = p prime,
  Frob_𝔭 = C and N𝔭 ≤ (c(Λ)(Λ + log 3))^32, Λ = log Δ_L, c = 0.484 (all Λ), 0.29 (Λ ≥ 10), 0.151 (Λ ≥ 100),
  0.108 (Λ ≥ 1000), 0.094 (Λ ≥ 10⁴). Only analytic input: ζ_L ≠ 0 on Re s > 31/32 (Cureton–Thorner, Lean-checked).
- Everything else is the paper's algebra (Mayle–Wang's 2-adic deviation group, Serre–Tate, Lemmas 4.x), re-checked in §2.

## 1. GRH audit of 2611.19556 (rg -n -i 'GRH|Riemann|Bach|Lagarias|Odlyzko|conditional|assum'; every hit read)
| line | what | GRH? |
|---|---|---|
| 137, 142 | intro (Serre, Mayle–Wang use Bach–Sorenson) | text only |
| 164–170 | Cor. 'total theorem' (K = Q(√D), formally real, abs(D) < 50) | "Assume GRH" (via Thm 5.4 and via Banwait–Derickx, §5) |
| 178–200, 916–945 | Thm 'Eff iso for K' | "Assume GRH": only through Lemma 'bach norm prime' |
| 359–375 | Thm 'bachthm' = Bach–Sorenson Thm 5.1: N𝔭 ≤ (4 log d_L + 2.5n + 5)², degree-1 𝔭 | **the GRH input** |
| 398–460 | Lemma 'bach norm prime' (L̃ = L·K(√α)) | only through 'bachthm' applied to L̃/K |
| 976–1056 | proof of Thm 'Eff iso for K': (4 log d + 2.5n + 5)² and Bach–Sorenson Table 1 rows (1.796, 1.89, 6.3), (1.336, 3.45, 45.7) | **GRH constants** (replaced) |
| 1173–1206 | Thm 5.4 ('5.4mainthm') | only through 'Eff iso for K' |
| 233–342 | Dickson, Borel (Prop. 'borel iso'), exceptional, C_ns, Thm 'bound theorem' (ℓ > 163) | no GRH written, **but** l.290 cites Banwait–Derickx Thm 1.9, whose Type 2 step needs GRH (§5) |
| 379–396, 535–660, 823–912 | Lemma 'prime ramify extn', 'bound disc', Serre–Tate, 'norm bound' | unconditional (algebra) |
Other hits ('assume' at 269, 276, 327, 610, 857, 1122): wlog phrases in algebraic proofs.

## 2. The paper's algebra kept (re-derived)
- **[L:K] ≤ 192** (l.950, Mayle–Wang proof of Thm 2): L = the field cut out by the 2-adic deviation group of
  (ρ_{E₁,2^∞}, ρ_{E₂,2^∞}); C ⊂ Gal(L/K) conjugation-closed, nonempty, with tr ρ₁ ≠ tr ρ₂ on C (Faltings + non-isogeny).
  Pure group theory, valid over every K. L/K unramified outside S := {𝔮 | 𝔪}.
- **Lemma 'norm bound'** (l.901–911): α ∈ O_K, v_𝔮(α) = 1 for 𝔮 ∈ S, log|N α| ≤ 2n_K log n_K + (n_K/2) log d_K + 2n_K ρ.
  Check: CRT gives α₀ ≡ a_𝔮 mod 𝔮², a_𝔮 ∈ 𝔮∖𝔮²; reduce coordinates mod u = (N𝔪)² in a basis w_j with
  ∏_σ Σ_j |σ w_j| ≤ n_K^{2n_K} d_K^{n_K/2} (Minkowski II: λ_n ≤ √d_K since λ₁ ≥ 1, and a basis with
  ‖w_j‖ ≤ ((j+1)/2)λ_j ≤ n_K √d_K); then |σα| ≤ u Σ_j |σw_j|. Correct.
- **Lemma 'prime ramify extn'** (l.379–396): v_𝔮(α) odd ⇒ 𝔮 ramifies in M = K(√α) (also for 𝔮 | 2). Correct.
- **Lemma 'bound disc'** (l.535, Hensel–Serre): log d_L ≤ [L:K] log d_K + ([L:K] − 1) log N(rad 𝔡_{L/K}) + [L:K] n_K log[L:K].
  Check: v_𝔓(𝔇) ≤ e − 1 + v_𝔓(e) gives log N𝔡_{L/K} ≤ (m−1)Σ log N𝔮 + m n_K Σ_p v_p(e_p) log p ≤ (m−1)ρ + m n_K log m.
  Increasing in m = [L:K], so m ≤ 192 gives (4.4): log d_L ≤ 192 log d_K + 191ρ + 192 n_K log 192.
- **(4.5)–(4.6)**: 𝔡_{L̃/L} | 4α (basis 1, √α) ⇒ d_{L̃} ≤ d_L² N_{L/Q}(4α) = d_L² (4^{n_K}|Nα|)^{[L:K]} (also if √α ∈ L).
- **R**: log d_{L̃} ≤ 2(192 log d_K + 191ρ + 192 n_K log 192) + 192(n_K log 4 + 2n_K log n_K + (n_K/2) log d_K + 2n_K ρ)
  = (382 + 384n_K)ρ + (384 + 96n_K) log d_K + 384 n_K log(384 n_K) =: R(K, ρ)   [384 n_K log 192 + 192 n_K log 4 = 384 n_K log 384].
  Recomputed in cert.py (exact Fraction coefficients; 384² = 4·192² for the logs) — matches l.1084–1088.
- [L̃:K] ≤ 384, [L̃:Q] ≤ 384 n_K (not needed below: chebexp's bound has no degree term).

## 3. The replaced steps (GRH → chebexp)
- **Envelope (Lemma 0).** F(Λ) := c(Λ)(Λ + log 3) jumps down at 10, 100, 1000, 10⁴, so use G(R) := sup_{0≤Λ≤R} F(Λ)
  (nondecreasing). Piecewise: Λ < 10: F ≤ 0.484(10 + log 3) < 5.374; [10,100): ≤ 0.29(100 + log 3) < 29.32;
  [100,1000): ≤ 0.151(1000 + log 3) < 151.17; [1000, 10⁴): ≤ 0.108(min(R,10⁴) + log 3) ≤ 1080.12; [10⁴, R]: 0.094(R + log 3).
  Hence **G(R) ≤ 0.108(R + log 3) for R ≥ 1400** (0.108·1401.09 = 151.32 > 151.17) and **G(R) ≤ 0.094(R + log 3) for
  R ≥ 11490** (0.094·11491.09 = 1080.16 > 1080.12). Also G(R) ≤ ĝ(R) := 0.108R + 152.3 for every R ≥ 0. (cert.py, Arb.)
- **Lemma 1 (replaces Lemma 'bach norm prime', l.398–460).** Same hypotheses, no GRH: there is 𝔭 ∉ S, unramified in L,
  N𝔭 = p prime, Frob_{𝔭,L/K} ⊂ C and N𝔭 ≤ G(log d_{L̃})^32. Proof: L̃ = LM is Galois over K (paper), L̃ ≠ Q
  ([L̃:Q] ≥ [M:Q] = 2n_K). Pick σ ∈ Gal(L̃/K) with σ|_L ∈ C, C̃ its class. chebexp for (L̃/K, C̃): 𝔭 unramified in L̃,
  so 𝔭 ∉ S (each 𝔮 ∈ S ramifies in M ⊂ L̃) and 𝔭 is unramified in L; Frob_{𝔭,L/K} = C̃|_L ⊂ C (C conjugation-closed).
  N𝔭 ≤ (c(Λ)(Λ + log 3))^32 ≤ G(log d_{L̃})^32 ∎. (Bach–Sorenson's 2.5n term disappears: chebexp has no degree term.)
- **Theorem A (replaces Thm 'Eff iso for K', l.916–1100).** K Galois, E₁, E₂/K non-CM, not K-isogenous. The paper's proof
  verbatim with Lemma 1: 𝔭 ∤ 2𝔑₁𝔑₂, a_𝔭(E₁) ≠ a_𝔭(E₂) (Mayle–Wang Prop. 12(3)), N𝔭 = p prime and
  N𝔭 ≤ G(R)^32, R = R(K, ρ) of §2. Since R ≥ 384 n_K log(384 n_K) ≥ 384 log 384 > 2285.0 > 1400:
    **N𝔭 ≤ (a ρ + b(K))^32, a = 0.108(382 + 384n_K) = 41.256 + 41.472 n_K,
      b(K) = 0.108((384 + 96n_K) log d_K + 384 n_K log(384n_K) + log 3) = (41.472 + 10.368n_K) log d_K + 41.472 n_K log(384n_K) + 0.1187;**
    and, when R(K, ρ) ≥ 11490, N𝔭 ≤ (a′ρ + b′(K))^32 with 0.094 in place of 0.108 (a′ = 35.908 + 36.096 n_K).
  K = Q: (82.729ρ + 246.904)^32 [(72.005ρ + 214.898)^32 for ρ ≥ 12.017]: serrenum's shape (17ρ + 100)^32, constants worse
  because the paper's L̃ doubles the degree and log|Nα| ≤ 2n_Kρ + … (§6).
  GRH (paper): (C₁ρ + C₂)², C₁ = 1.336(382 + 384n_K) = 510.352 + 513.024n_K — same structure, 1.336 → 0.108, square → 32nd power.
- **Typo in the paper's C₂.** The proof (l.1094–1100) gives C₂ = (513.024 + 128.256n_K) log d_K + 513.024 n_K log(384n_K)
  + 1324.8n_K + 45.7; the statements (l.153–156, 939–942) write 513.024(1 + (n_K/4) log d_K + n_K log(384n_K)) + …,
  i.e. 513.024 where 513.024 log d_K is proved: the stated C₂ is too small by 513.024(log d_K − 1) when d_K > e
  (e.g. 554 for Q(√47)). Minor; our b(K) follows the proof.

## 4. Theorem B (replaces Thm 5.4, l.1173–1206) — the non-split Cartan step
K Galois, E/K non-CM, ρ_E := log N rad(2O_K𝔑_E), ℓ prime with ℓ > 2G(R(K, ρ_E))^16 (≤ 2(aρ_E + b(K))^16). Claim:
if G_E(ℓ) ⊂ C_ns^+(ℓ) then G_E(ℓ) ⊂ C_ns(ℓ). (The paper states "not contained in C_ns^+" but its proof, l.1181, only treats
"in C_ns^+, not in C_ns"; the full statement holds for formally real K, where complex conjugation excludes C_ns.)
- ε: G_K → C_ns^+/C_ns = {±1}, nontrivial. **ε unramified at every 𝔭 ∤ 𝔑_E**, also 𝔭 | ℓ (paper's Lemma 'character
  properties'(i) cites Mayle–Wang Lemma 15 over Q; over K, with e = e(𝔭/ℓ) ≤ n_K): supersingular at 𝔭 | ℓ ⇒ ρ̄(I_𝔭) cyclic
  of order (ℓ²−1)/gcd(e, ℓ²−1) ≥ (ℓ²−1)/e (Raynaud, e < ℓ − 1); a cyclic subgroup of C_ns^+ not in C_ns has order | 2(ℓ−1)
  (g ∈ C_ns^+∖C_ns ⇒ g² scalar), and (ℓ²−1)/e > 2(ℓ−1) for ℓ > 2e − 1, so ρ̄(I_𝔭) ⊂ C_ns. Ordinary at 𝔭 | ℓ: excluded below.
- ε = χ_D (K(√D)/K unramified outside 𝔑_E); E_D := twist; E_D non-CM, not K-isogenous to E (Lemma 'E and ED');
  rad 𝔑_{E_D} | rad 𝔑_E (ρ_{E_D,ℓ'} = ρ_{E,ℓ'} ⊗ ε unramified where E is good and ε unramified), so 𝔪(E, E_D) = rad(2O_K𝔑_E).
- Theorem A for (E, E_D): degree-one 𝔭₀ ∤ 2𝔑_E, N𝔭₀ ≤ G(R)^32, a_{𝔭₀}(E) ≠ a_{𝔭₀}(E_D) = ε(Frob_{𝔭₀}) a_{𝔭₀}(E).
  Hence ε(Frob_{𝔭₀}) = −1 and a_{𝔭₀}(E) ≠ 0.
- **𝔭₀ ∤ ℓ (step the paper skips, as Mayle–Wang do over Q).** If 𝔭₀ | ℓ then N𝔭₀ = ℓ, 0 < |a_{𝔭₀}| ≤ 2√ℓ < ℓ, so
  E is ordinary at 𝔭₀; then ρ̄(I_{𝔭₀}) ∋ g = (u *; 0 1) with u ranging over χ_ℓ(I_{𝔭₀}), a cyclic group of order
  (ℓ−1)/gcd(e, ℓ−1) ≥ 3 for ℓ > 2e + 1; take u ≠ ±1: g has distinct eigenvalues u, 1 ∈ F_ℓ, so g ∉ C_ns, and tr g = u + 1 ≠ 0,
  so g ∉ C_ns^+∖C_ns (trace 0 there): contradiction with G_E(ℓ) ⊂ C_ns^+. (ℓ > 2n_K + 1 is implied: 2G^16 ≥ 2·151^16.)
- 𝔭₀ ∤ ℓ𝔑_E and ε(Frob) = −1 ⇒ Frob_{𝔭₀} ∈ C_ns^+∖C_ns ⇒ a_{𝔭₀} ≡ tr = 0 mod ℓ ⇒ ℓ ≤ |a_{𝔭₀}| ≤ 2√N𝔭₀ ≤ 2G(R)^16. ∎
  Numbers: ℓ ≤ 2(aρ_E + b(K))^16 with a, b(K) of §3 (table in REPORT; e.g. Q(√5): 2(124.201ρ_E + 651.302)^16).

## 5. Corollary (replaces Cor. 'total theorem', l.164–170): a second GRH input, Banwait–Derickx Type 2
- The paper's Thm 'bound theorem' (l.338–341, K = Q(√D), 0 < D < 50, ℓ > 163 ⇒ surjective or in a Cartan normalizer but
  not in the Cartan) is stated without GRH, but its Borel step (l.290) cites Banwait–Derickx (ANT 19 (2025);
  arXiv 2203.06009v3) Thm 1.9 ('thm:more_quadratic', TeX l.352), proved with their Algorithm 'combined', whose Type 2 part
  needs GRH: TeX l.1457 "This is the only result in our work which requires the Generalised Riemann Hypothesis"
  (Prop. 'type_2_bound': p ≤ (8d log(12p) + 16 log Δ_k + 10d + 6)⁴, via Larson–Vaintrob Cor. 6.3 = Bach–Sorenson).
  Unconditionally BD give: every K-isogeny prime p > 163 is of Momose Type 2 (Type 3 needs K ⊃ a Hilbert class field
  of an imaginary quadratic field: impossible for real K) and satisfies their Condition CC (Prop. 'cond:CC', l.1427, no GRH).
- **Lemma 2 (Type 2 primes, unconditional).** K = Q(√D) real, d_K ≤ 188. No K-isogeny prime of Momose Type 2 has
  log p ≥ y* = 162.0 (p ≥ 10^70.356). Proof. F := Q(√D, √−p), biquadratic (K real), Δ_F = |d_K d(−p) d(−pD)| ≤ 16p²d_K²
  (conductor–discriminant). chebexp (base Q, class {1}): q₁ splits completely in F, q₁ ≤ G(R_F)^32,
  R_F = 2 log p + 2 log d_K + log 16. L₂ := F(√q₁): q₁ ramifies in Q(√q₁) ⊂ L₂ but not in F, so [L₂:F] = 2 and
  log Δ_{L₂} ≤ 2 log Δ_F + 4 log(4q₁) (as (4.5)); chebexp (class {1}): q₂ ≠ q₁ splits completely in L₂ ⊃ F,
  q₂ ≤ G(R₂)^32. cert.py: for every y = log p ≥ 162.0, q₁, q₂ < p/4 (7352 intervals of width ¼ on [162, 2000], then a
  monotone tail with G(R) ≤ 0.108R + 152.3 ≤ R). Each q_i splits in K (a degree-1 prime 𝔮 | q_i, f = 1 odd) and in
  Q(√−p), q_i < p/4. If both had q_i² + q_i + 1 ≡ 0 mod p they would be the two primitive cube roots of 1 mod p
  (distinct residues in (0, p)), so p | q₁ + q₂ + 1 ∈ (0, p/2 + 1): impossible. So one q_i meets CC's (1)–(3) and splits
  in Q(√−p), contradicting CC. ∎
- **Corollary C (unconditional).** K = Q(√D), D squarefree, 1 < D < 50; E/K non-CM; ρ̄_{E,ℓ}(G_K) ⊄ C_s^+(ℓ). Then ρ̄_{E,ℓ}
  is surjective for every ℓ > max(e^{162.0}, 2G(R(K, ρ_E))^16), in particular for ℓ > max(10^70.36, 2(124.201ρ_E + b(K))^16).
  Proof: paper's case analysis with Borel ⇒ ℓ ≤ 163 or log ℓ < 162 (Lemma 2); exceptional (ℓ > 7, l.292–307) and C_ns
  (complex conjugation, l.309–330) kept from the paper (no GRH; not re-derived); C_s^+ excluded by hypothesis;
  C_ns^+∖C_ns: Theorem B. Type 2 dominates while 2G^16 < 10^70.36, i.e. N rad(2O_K𝔑_E) < about 10^93.

## 6. Remarks
- The constants are ~6× (a) worse than serrenum's because the paper's L̃ doubles [L:K] and uses log|Nα| ≤ 2n_Kρ + …;
  serrenum's excluded-prime term (n_L‖k‖ log m/x inside chebexp's master inequality) would avoid L̃ and give the
  coefficient 0.1·191 ≈ 19 per unit ρ for every K; not done here (needs new chains).
- K = Q sanity: Theorem A gives (82.729ρ + 246.904)^32, serrenum (17ρ + 100)^32: same shape, as required.

## 7. Constants → cert.py (rg -n -F '<constant>' cert.py)
| constant | meaning | cert.py line(s) |
|---|---|---|
| 191, 192 (4.4); 2, 1/2, 2 (4.2); 192·(n log 4 + log|Nα|) (4.6) | R's ingredients, exact Fractions | 32–34 |
| 382 + 384n, 384 + 96n, 384n log(384n) | R(K, ρ) and its check | 36–40 (exact), 71–75 (Arb), 131 |
| 0.484, 0.29, 0.151, 0.108, 0.094 | chebexp table (decimal Arb balls) | 17 |
| 5.374, 29.32, 151.17, 1080.12 | envelope piece caps | 45 |
| 1400, 11490 | envelope thresholds | 46, 47 (asserted), 75, 78 |
| 152.3 (ĝ), 171, 0.892 | tail envelope ĝ(R) = 0.108R + 152.3 ≤ R for R ≥ 171 | 48, 112 |
| 2 log p + 2 log d_K + log 16 | log Δ_F (Lemma 2) | 97 |
| 2R_F + 4 log 4 + 4 log q₁ | log Δ_{F(√q₁)} (Lemma 2) | 99 |
| y* = 162.0, width ¼, Y₀ = 2000, tail t1–t4 | Type 2 chain | 101–115 |
| 510.352, 513.024, 128.256, 1324.8, 45.7 | the paper's GRH constants (comparison only; C₂ as proved) | 91 |
Each a, b(K), a₁, b₁(K) is rounded up to 3 decimals and certified > its exact Arb value. It is also checked against
0.108(R + log 3) at ρ = 1, 2, which is affine in ρ, to within 10⁻⁴⁰ (l.62–80).

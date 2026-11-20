# serrenumK follow-up — DERIVATION2: Momose Type 2 primes over every number field; Larson–Vaintrob without GRH

Inputs: chebexp (A) (every Galois L/K, L ≠ Q, every class: unramified degree-one 𝔭, N𝔭 ≤ G(log Δ_L)^32, envelope G as in
DERIVATION §3); Banwait–Derickx (BD, arXiv 2203.06009v3) Prop. 'cond:CC' (TeX l.1427–1450, no GRH), Def. Momose Type 2
(l.858: λ^{12} = χ_p^6), Thms 'main-brief', 'type_1_brief' (l.251–263: "Let k be a number field", nonzero integers MMIB(k),
TypeOneBound(k)), Type 3 needs k ⊇ HCF of an imaginary quadratic field (l.1806); LV14 (arXiv 1110.0255) Thm 1 (l.186–212).

## 1. Lemma 2 for every k (replaces BD Prop. 'type_2_bound' = LV Cor. 'avo', the only GRH step of BD, l.1457)
k any number field, M its Galois closure over Q, n_M = [M:Q], Δ_M = |disc M|, y = log p.
- F := M(√−p), Galois over Q. If √−p ∈ M, F = M; else [F:M] = 2 and 𝔡_{F/M} | (−4p), so
  log Δ_F ≤ 2 log Δ_M + n_M log(4p) =: R_F(y). chebexp over Q, class {1}: q₁ splits completely in F, q₁ ≤ G(R_F)^32.
- F₂ := F(√q₁): q₁ is unramified in F but ramified in Q(√q₁), so [F₂:F] = 2, and 𝔡_{F₂/F} | 4q₁, [F:Q] ≤ 2n_M:
  log Δ_{F₂} ≤ 2R_F + 2n_M log(4q₁) =: R₂. chebexp, class {1}: q₂ splits completely in F₂, q₂ ≠ q₁ (unramified), q₂ ≤ G(R₂)^32.
- q_i splits completely in M ⊃ k, so every prime 𝔮 | q_i of k has f = 1 (odd), and q_i splits in Q(√−p).
- cert2.py: q₁, q₂ < p/4 for every y ≥ y*(k) (table). Then at most one of q₁, q₂ has q² + q + 1 ≡ 0 mod p (two such
  would be ω, ω² mod p, so p | q₁ + q₂ + 1 < p/2 + 1). That q_i satisfies CC's (1)–(3) and splits in Q(√−p): CC fails.
- **Lemma 2(k).** No elliptic curve over k has a k-rational p-isogeny of Momose Type 2 with log p ≥ y*(k). No hypothesis on k.
- **Lemma 2′ (closed form, every k).** y*(k) ≤ Y₀(k) := max(3000 n_M, 2 log Δ_M + n_M log 4). Proof for y ≥ Y₀:
  R_F ≤ (n_M + 1)y ≤ 2n_M y (y ≥ 2 log Δ_M + n_M log 4); G(R) ≤ ĝ(R) = 0.108R + 152.3 ≤ R for R ≥ 171 (R_F ≥ n_M y ≥ 3000);
  so log q₁ ≤ 32 log(2n_M y); R₂ ≤ 4n_M y + 2n_M log 4 + 64n_M log(2n_M y) ≤ 5n_M y iff y ≥ 2 log 4 + 64 log(2n_M y);
  log q₂ ≤ 32 log(5n_M y) < y − log 4. Both inequalities hold at y = 3000n_M for every n_M ≥ 1
  (3000n − 2.78 − 64 log 6000 − 128 log n > 0, 3000n − 1.39 − 32 log 15000 − 64 log n > 0 at n = 1, increasing in n) and
  their left minus right sides are increasing in y (derivatives 1 − 64/y, 1 − 32/y > 0). ∎
- Earlier Lemma 2 (REPORT.md, real quadratic) is the case M = k = Q(√D): R_F = 2y + 2 log d_K + log 16, same R₂.

## 2. Isogeny primes without GRH (LV14 Cor. 'isogeny')
- **Theorem 3.** Let k be a number field not containing the Hilbert class field of an imaginary quadratic field. Then every
  prime p such that some E/k has a k-rational p-isogeny satisfies p | MMIB(k)·TypeOneBound(k) or log p < y*(k) ≤ Y₀(k).
  So IsogPrimeDeg(k) is finite, with an explicit bound computable from k, unconditionally.
  Proof: BD's isogeny character λ is (a) not of Momose Type 1, 2, 3 ⇒ p | MMIB(k) (Thm 'main-brief'; covers signature
  Type 2 not Momose Type 2 via Cor. 'type_2_not_momose', l.887–910); (b) Momose Type 1 ⇒ p | TypeOneBound(k) ('type_1_brief');
  (c) Momose Type 2 ⇒ log p < y*(k) (Lemma 2(k)); (d) Momose Type 3 ⇒ k ⊇ HCF (excluded). BD state (l.1457) that only
  their Type 2 bound needs GRH; (a), (b) hold for every k, Galois or not (their non-Galois versions of Momose's Lemmas 1, 2:
  Props. 'momose_1_non_galois_md' l.596, 'momose_lemma_2_non_galois' l.719, replacing Momose's flawed passage to the
  Galois closure, l.236–247). ∎ The 'only if' half of LV's Cor. (CM curves) is unconditional anyway.
- **Corollary 4 (LV14 Thm 1 with case (2) removed).** For every K and every prime ℓ ∉ S_K with log ℓ ≥ y*(K), alternative (2)
  of LV Thm 1 (ψ^{12} = χ_ℓ^6, "GRH fails for K[√−ℓ]") cannot occur. Proof: LV state that in case (2) ρ̄_{E,ℓ} is reducible
  over F_ℓ (l.205–206), so E has a K-rational ℓ-isogeny with character λ, ψ ∈ {λ, χ_ℓλ^{−1}}, and either way λ^{12} = χ_ℓ^6:
  Momose Type 2, impossible by Lemma 2(K). ∎ This is Harvest4 item 5's goal.

## 3. Consequences checked in the citing papers (only what was read)
- **Allen–Genao 2511.23381.** GRH appears only in the abstract (l.50), the convention (l.75), Thm
  'UniformBds_OnPrimeLevel_OfAbelianDivisionFields' (l.77: "if GRH is true and F has no rational CM, then p is contained in
  a finite set S_F from [LV14, Theorem 1]"), the quoted LV Thm 1 (l.265–280) and its use (l.283).
  - Their l.75 definition: "a number field is said to have rational CM if it contains an imaginary quadratic field".
  - **Case (1) excluded by LV's own conclusion**, not merely by "a CM field". LV Thm 1 (1110.0255 final.tex l.192–196),
    verbatim: "There exists a CM elliptic curve $E'$, which is defined over $K$ and whose CM-field is contained in $K$,
    with an $ell$-adic weight $1$ associated character whose mod-$ell$ reduction $psi'$ satisfies $psi^{12} =
    (psi')^{12}$".
  - So case (1) produces a CM elliptic curve over F all of whose endomorphisms are defined over F. Its CM field is an
    imaginary quadratic field inside F, so F has rational CM in Allen–Genao's sense (l.75). It also has rationally defined
    CM in Genao's sense (2405.05507 l.313); then F ⊇ K(j(E')), a ring class field, hence the Hilbert class field of K.
  - For F without rational CM, case (1) is impossible for p ∉ S_F.
  - **Case (2)** forces an F-rational Momose Type 2 isogeny (Corollary 4), impossible for log p ≥ y*(F).
  - **Consequence (unconditional).** If F has no rational CM, F(E[p])/F is abelian, and p satisfies the hypotheses of l.77,
    then p ∈ S_F ∪ {p : log p < y*(F)}. This is an enlarged exceptional set, not necessarily the original S_F: the GRH
    statement's set S_F is replaced by S_F ∪ {p : log p < y*(F)}, with y*(F) ≤ max(3000n_M, 2 log Δ_M + n_M log 4).
- **Genao 2405.05507.** Inside the paper's proofs GRH is used once (l.507): LV Thm 1 for the F-rational ℓ-isogeny ⟨P⟩
  (proof of Thm 'ModEllRep_ofEllCurve_With1_EllIsogeny' a), to discard case (2). Corollary 4 discards it for log ℓ ≥ y*(F),
  so that step holds with c(F) replaced by max(c(F), ⌈e^{y*(F)}⌉) (c(F), l.182: ℓ ≥ max(74, 15[F:Q]+2), unramified, ∉ S_F).
  But the GRH theorems (l.145, 151, 166, 172) also use Thm 'BaseChangeImprovement' (l.150, "under GRH"), described as an
  improvement of [Gen24, Thm 2] and cited in the proofs (l.405, 407, 445, and once in §477–580); its GRH input lies in [Gen24],
  **not checked**: no unconditional claim for 2405.05507 beyond the l.507 step.
- **LV's effective S_K.** LV's Thm 'effectivethm' (product of ℓ ∈ S_K) is stated under GRH (l.3056), but LV Remark 'besteff'
  (l.3068) says the bound ∏ℓ ≤ b_char·b_rat is proven independently of GRH, and the remark at l.3044–3049 that an
  unconditional Chebotarev makes it unconditional; for Theorem 3 we do not need S_K (BD's integers replace it).

## 4. Comparison with the planned route (Harvest4/HARVEST.md item 5)
Planned: chebexp in LV's own step (Cor. 'avo' on K(√±ℓ, √±3), trivial class, "not over 3"), then ℓ ≤ (1 + √N v)^4
(LV l.2660–2684), i.e. log ℓ ≲ 64 log G(·): exponent 64, plus the unproved claim that a split prime avoids 3 (fails when
3 ramifies in K with even index, e.g. K(√−3)/K unramified above 3). Here: BD's Condition CC needs only q < p/4, so
exponent 32, any q (no condition at 3), at the price of a second prime. Float comparison (cert2.py, BD's bound
log Δ_{E'} ≤ 4 log Δ_k + 2d log(12ℓ)): planned route y ≈ 321–440 vs y*(k) = 162–218 on BD's table fields.

# Liu, arXiv 2610.26835, general-field step: working notes (2026-10-21)

Brief: `.scratch/orca/liuhecke-brief.md`. Verdict: [REPORT.md](REPORT.md), NO REFUTATION.
Sources read as TeX or text:
- Liu `number_fields_v2.tex` (arXiv e-print);
- Cureton–Thorner `DedekindZeta.tex` (Hexagon 2610.00200v2);
- Hu `weil_seven_eighths_revised.tex` (Hexagon 2610.00216v1);
- OpenAI primitive roots, PDF to text (cited as [Artin] in Liu);
- Blomer–Goldmakher–Louvel, arXiv 1112.1650v1 TeX;
- Goldmakher–Louvel, arXiv 1112.1642v2 TeX.

Line numbers refer to Liu's TeX.

## 1. Structure of the claim
- Thm vi:thm:main (l. 10626): every finite-order L_K(s,χ), for every K, is nonzero for Re s > B♭.
- Proof: F = K(ζ12); Thm vi:prop:transfer (l. 10740) gives β_*(F) ≤ B♭; Prop vi:prop:descent (l. 10875) descends to K.
- Thm vi:prop:transfer says the proofs of Part II (cyclotomic F ⊇ μ12) "apply verbatim" to any F ⊇ μ12. Its proof is a list of the facts used, with their sources.

## 2. Does Part II use more than μ12 ⊂ F?
- Search of l. 1585–10423 with rg:
  - "cyclotomic" and Q(ζ_m) appear only in statements and set-up lines, e.g. l. 1775, 4782, 7149, 7285, 8124, 8802, 9850, and the remark at l. 3078 ("no uniformity in the field");
  - Galois, Frobenius, automorphism and embedding appear only in local Kummer and Artin-map uses, and in the fixed embedding ȷ of μ6;
  - nothing uses CM structure, Z[ζ_m], ramification above m, or splitting laws.
- Cyclotomic structure is used only in the Part II descent (§10424ff, Gal(F/Q) abelian, metabelian scope) and in the class 𝒞. Part III replaces both.
- The cited source [Artin] is written the same way. §3 opens with "Fix a cyclotomic field F containing μ12", and §5 with "F is the fixed cyclotomic field … In particular, μ12 ⊂ F". All its Neukirch citations are general.
- Cureton–Thorner (Prop. "field-scope" and §§ "Arithmetic choices over a general field", "Global realization") check every construction for F ⊇ μ12 that need not be Galois or CM, and formalize the KP cubic rank-two constructions in Lean.

## 3. Ingredients of Part II's base case and their status for general F
Key: (a) proved in this paper; (b) in the literature for every F ⊇ μ12; (c) missing or false (none found).

| Ingredient (locator) | Status | Remark |
|---|---|---|
| S, R = O_{F,S} with trivial class group, allowed generators, H_S isotropic complement, sextic reciprocity, unit characters Υ_S (Lemmas v:lem:allowed-reciprocity l. 1884, v:lem:unit-products) | (b) | [Artin] Lemmas 3.1–3.3; CT26 Lemma "generators". The dimensions 2(\|S_f\|+r_F) and \|S_f\|+r_F were re-derived. The pairing is alternating since −1 = ζ12⁶. |
| Gauss corrections Λ = Λ_{Q(ω)}∘N, G_0 (Lemma v:lem:theta-characters l. 2324) | (b) | [Artin] Lemma 3.4. Davenport–Hasse prime by prime, γ2(p)³ = −Λ(p). Weil-index product formula. At each complex place the angular type is ±1. No normality needed (CT26 Lemma "corrections"). |
| Exceptional theta rep, exterior Whittaker table W_e, recurrence W_{e+3} = q⁻¹W_e (Prop v:prop:checked-inputs l. 2023) | (b)+(a) | KP over any global field ⊇ μ3, c = 2; tame uniqueness Cor I.3.6 (2(c+1) = 6 ≡ 0 mod 3). The table is re-derived in the paper. The S-functional is a finite sum of local products, so wild uniqueness is not used. Liu's citation of KP II.2.5 for wild places is therefore unneeded; CT26 cites Cai–Fan–She. |
| Cube in every completion implies cube (Lemma vi:lem:mu12(5)) | (b) | Grunwald–Wang has no exception for odd exponent. |
| Radial functional equation, derivative vectors of angular type ±1 (Prop v:prop:radial l. 2884) | (a)+(b) | [Artin] §5 (Bessel K_{1/3}), extended in the paper to arbitrary compact radial tests. Gamma multiplier raised to the power r_F. |
| Orbit Poisson, compact dual support, norm-sheet kernels (Prop v:prop:orbit-poisson l. 3290, Lemma v:lem:compact-dual l. 3394, §v-detail-8) | (a) | Self-dual lattice R ⊂ F_S with covolume 1. The support bound x > C_Φ ⇒ 𝒦 = 0 is correct. The kernel log power d_F = r_F + \|S_f\| − 1 is the dimension of the norm sheet. |
| Quadratic large sieve (Lemma v:lem:primary-kummer-transport l. 2934) | (b) | GL Thm 1.1: "Let 𝔠 be an integral ideal of a number field k … ≪ (MN)^ε(M+N)". The family conditions of Def. 1 are checked in the lemma. |
| Cubic large sieve (same lemma) | (b) | BGL Thm 1.3 (thm:thm3): K ⊇ μ_n, n ≥ 3, squarefree ideals prime to S, (MN)^ε(M+N+(MN)^{2/3}). BGL App. II gives the Kummer construction. |
| Inverse moment M1, unmarked (Thm v:thm:unmarked-inverse l. 4689) → seed β_* ≤ 11/12 (Thm v:thm:unmarked-seed l. 4781); marked (Thm v:thm:marked-inverse-new l. 5197) | (a) | Modelled on OpenAI §17, which is Lean-checked over Q(√−3). Hu has a parallel version (Thm "can-inverse"). CT26 deliberately avoided this general-field Möbius mean square. |
| Product moment M2 (Thm v:thm:M2-assembled l. 6971): 13/18 ≤ κ ≤ 1, z = 0 at arbitrary bounded lengths, hypothesis β_*(F) ≤ (1+κ)/2 when κ < 1 | (a) | Modelled on Liu 2610.12234 Prop 28.1, which extends OpenAI L18.1 below its printed κ ≥ 3/4. The hypothesis is a uniform zero-free bound for the same F, discharged by the contradiction with κ = 2β_*(F) − 1. |
| Hecke functional equation and AFE, prime ideal theorem in ray classes, log control of L (Lemma v:lem:refined-afe, Lemma vi:lem:mu12(2)) | (b) | Standard (Hecke; Iwaniec–Kowalski Thm 5.3). Γ(s)^{r_F} because F is totally complex. |
| Detector, row count R_κ, nonzero principal correction (Lemma v:lem:general-field-detector l. 8123; §§v-detail-28/29) | (a) | Small primes are put in S so the principal Euler factor stays nonzero (CT26 does the same). |
| Averaged numerator (Prop v:prop:refined-numerator l. 9158) | (a) | Hölder plus amplified M2 (Xu–Sang). Lengths m ≤ 1/2 + r_F/100 are covered: H ≤ U^{1/100} because τ_η < d_min/10⁴ with d_min = 1/100 (l. 9576); for m ≥ 1/2, D^{−θ} ≤ U^{−θ/2}. |
| Scalar system, tangent family, quartic for B♭ (Lemma v:lem:refined-family; Thm v:thm:refined-limit l. 9849) | does not depend on F | Reproduced earlier to 10⁻¹² (Riemann/Deep/QRHImprove). |
| Order of quantifiers (Lemma v:lem:arithmetic-independent-mesh l. 7144) | (a) | Mesh and margins are numerical and chosen before F. Constants, height orders and Z-thresholds may depend on F, S and η. |

## 4. Descent and scope (Prop vi:prop:descent)
- The three factorizations are standard class field theory and hold for full primitive L-functions, including characters ramified at real places.
- Nonprincipal factors are entire, so no pole cancels a zero.
- At s = 1, the ratio of residues is positive.
- An imprimitive L-function differs by finite Euler factors, whose zeros lie on Re s = 0.
- Equivalent to the brief's L_F(s,χ∘N) = ∏_ψ L_K(s,χψ) (Hu Lemma "descent", same proof).
- No hidden height threshold:
  - β_*(F) is a supremum over all zeros;
  - the Mellin continuation of e^{(s−5/6)²}H(s)/L^S(s,η) excludes zeros at every height;
  - constants depending on K only move thresholds in the averaging scale Z.
- The boundary B♭ is the same for every K and χ.
- Corollary vi:cor:consequences:
  - (1) and (3) are immediate;
  - (2) correctly limits itself to direct sums of monomial representations. General Artin L-functions need Brauer induction plus pole bookkeeping, done in Hu §9.

## 5. Comparison with the other general-field papers
- **Cureton–Thorner**, Hexagon 2610.00200v2:
  - result: 31/32 for ζ_K and all finite-order Hecke L over every K, Lean-checked (AxiomMath/DedekindZetaZeroFree);
  - route: the [Artin] comparison (theta upper bound plus Poisson principal term), which needs no general-field Möbius mean square;
  - they say the 7/8 refinements "would impose further requirements beyond the cubic transformation and quadratic sieve".
- **Hu**, Hexagon 2610.00216v1, 13 Oct:
  - result: 7/8 for all unitary Hecke characters, of finite or infinite order, over every number field; Weil and Artin L-functions follow by Brauer induction;
  - route: F(ζ12) plus abelian descent;
  - App. A tabulates the field-sensitive inputs (generators, Gauss calibration, orbit Poisson, theta reflection, sieve identification, the two inverse transforms, ordinary reflection, fourth-moment recursion). This matches Liu's list;
  - the recursions of OpenAI §§17–18 are transferred by citing their locations, not reproduced. Not Lean-checked.
- **Rubinstein–Restrepo Borrero**, Hexagon 2610.00220v1: Q(√−3) only (y² = x³ + d; angular characters ᾱη). Not a general-field result.
- The brief's "arXiv 2610.00216/00220" are Hexagon identifiers, and Hu is not a special-field paper.

## 6. What remains unverified
These steps do not depend on F; a failure in them would already break Q(√−3).
- M2 for 13/18 ≤ κ < 3/4, which the B♭ step needs (κ♭ = 2B♭ − 1 = 0.749773).
- The M2 clause without prime sums at arbitrary lengths.
- The Xu–Sang averaged numerator.
- The general-field M1/M2/count derivations of Part II (~5,000 lines). They were spot-read here (statements, Lemma v:lem:primary-kummer-transport, orbit Poisson, compact dual, the M2 induction outline, Prop v:prop:refined-numerator), not line-checked.

Lean status:
- OpenAI 7/8 (ζ and Dirichlet L-functions): Lean-checked; our port builds with no sorry (Riemann/Deep/QRH/REPORT.md).
- Liu's Lean: scalar and selected arithmetic lemmas only (Remark v:rem:analytic-inputs).
- CT26's Lean: descent and KP constructions over F ⊇ μ12.

No scripts were run, so there is no CPU or memory cost.

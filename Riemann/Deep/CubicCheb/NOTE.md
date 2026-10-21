# CubicCheb: Chebotarev for almost all cubic fields from the ring class results

Date 2026-10-21. Brief: `.scratch/orca/cubiccheb-brief.md`. Scratch: `.scratch/cubiccheb/` (`dict.gp` → `dict.out`
dictionary check; `least.gp` → `least.out` least primes; `ax.py`, `ax2.py` → `q3.out` arXiv API; `gxl.py` → `gxl_c1.out`,
`gxl_c2.out` GXL evidence API; sources read as text: `ptbw.txt`, `lotz.txt`, `cloz.txt`, `p_1606.06103.txt` (EPW),
`p_2405.08383.txt` (Lemke Oliver–Smith), `p_1601.03012.txt`, `p_2004.14410.txt`). [C], [D2], [R], [R2], [R3], [R4] as in
[RCApps4](../RCApps4/NOTE.md). Paper proofs, not refereed, no Lean, no 2026 half-planes.

## Verdict

| Item | Result | Status |
|---|---|---|
| entry 2 ⇒ cubic fields (Theorem B) | A > A₀ = (3+√69)/10: for all but O(X^{1−c}) fundamental D, every cubic field of discriminant D has primes of each S₃ class in every (x, x+x\|D\|^{−κ}], x ≥ \|D\|^A | routine corollary; gives nothing for almost all *fields* (c < 1.4·10^{−4} < 1/3); beaten by EPW, PTBW, LOTZ; for 3-cycles even by the every-field bound \|d_K\|^{3/4+ε} |
| strongest version with our tools (Theorem A) | any A > 0: all but O(X^ε) fundamental D, asymptotic count of each class in (x, x+x\|D\|^{−κ}]; all but O(X^{1/3+ε}) cubic fields with fundamental discriminant | routine: the PTBW mechanism with [R] Prop 1 as input; also follows from LOTZ's zero-free region plus classical facts; the family large sieve plays no role |
| non-fundamental d_K (27% of cubic fields) | not covered by entry 2 or [R] Prop 1 (Remark 2) | not claimed |
| overall | | **NO DISCOVERY** |

## 0. The dictionary

D fundamental, F = Q(√D), G = Cl(F) (wide class group, also for D > 0), h = |G|.

**Lemma 1** (Hasse). (i) Cubic fields K with d_K = D, up to isomorphism, correspond to subgroups H ≤ G of index 3, i.e. to pairs
{ψ, ψ̄} of characters of G of order 3 (H = ker ψ); there are (|G[3]| − 1)/2 of them. Each is an S₃-field, and its Galois closure
K̃ is the class field of F for H. (ii) For a prime p ∤ D: Frob_p(K̃/Q) is a transposition iff (D/p) = −1. If p = 𝔭𝔭̄ splits in F,
Frob_p = 1 iff ψ(𝔭) = 1, and Frob_p is a 3-cycle iff ψ(𝔭) ≠ 1.

*Proof.* (i) For a non-Galois cubic K with resolvent F, d_K = d_F f², where fO_F is the conductor of K̃/F (Hasse). d_K = D = d_F
forces f = 1, so K̃/F is unramified at finite primes, and at infinite ones because [K̃:F] = 3: K̃ lies in the Hilbert class field
and corresponds to an H of index 3. Conversely the class field L_H is Galois over Q, because Gal(F/Q) acts on G by inversion
(𝔞𝔞̄ = (N𝔞)) and so fixes H, with group (G/H) ⋊ C₂ ≅ S₃; its cubic subfields have discriminant d_F·1² = D. (ii) An inert p has
Frob_p ∉ Gal(K̃/F) = A₃. For split p, Frob_p ∈ Gal(K̃/F) is the Artin symbol of 𝔭, i.e. the image of [𝔭] in G/H. □

*Check* (`dict.gp`): 14 fundamental D (−23, −31, −59, −83, −107, −283, −3299, −4027, −11651, 229, 257, 316, 1957, 32009; 3-rank 2
for four of them), all 26 cubic fields K (cubic subfield of `bnrclassfield(F, ker ψ)`): `nfdisc(K) = D` every time; for every
prime 3 ≤ p ≤ 3000 with p ∤ D·disc(f_K), the factorisation type of f_K mod p is the one predicted by (ii): 11,135 tests, 0 mismatches.

*Counts.* By Davenport–Heilbronn (mean of |G[3]| is 2 for D < 0 and 4/3 for D > 0; N₃^∓(X) ~ X/(4ζ(3)), X/(12ζ(3))), cubic fields
with fundamental discriminant are a proportion 6ζ(3)/π² = ζ(3)/ζ(2) ≈ 0.7308 of all cubic fields, in each signature. The rest have
d_K = d_F f², f > 1, and correspond to primitive cubic ring class characters of conductor f (Remark 2).

## 1. Theorem A: the strongest statement our tools give

**Theorem A.** Let A > 0, 0 < ε ≤ 1, η := ε/100, κ := η·min{A, 1}/4. There is an effective X₀(A, ε) such that for X ≥ X₀ all
fundamental D with |D| ≤ X outside a set E(X), #E(X) ≤ X^ε, satisfy: for every cubic field K with d_K = D, every conjugacy
class C ⊆ Gal(K̃/Q) ≅ S₃, every x ≥ |D|^A and every H ∈ [x|D|^{−κ}, x],

  #{p ∈ (x, x+H] : p ∤ D, Frob_p ∈ C} = (|C|/6)·(H/log x)·(1 + O_{A,ε}(|D|^{−κ/2} + exp(−c(log x)^{1/4}))).

In particular every class has a prime ≤ 2|D|^A. Moreover the cubic fields with fundamental discriminant |d_K| ≤ X and d_K ∈ E(X) number ≪_ε X^{1/3+2ε}.

*Proof.* **Exceptional sets** (dyadic, D ∈ 𝔇(X) ∪ 𝔇⁺(X)); T_X := X^{3κ} ≤ X^{3η/4}.
- E₁(X): some χ ∈ Ĝ_c(D) (complex characters, so every cubic ψ) has a zero with β ≥ 1 − η, |γ| ≤ T_X. By [R] Prop 1 (D < 0)
  and [R2] Prop 2 (D > 0, narrow class characters; ψ∘(G⁺ → G) has the same L-function), with ε₀ := η:
  #E₁(X) ≤ Σ_D N^c_D(1 − η, T_X) ≤ C(η)(XT_X)^{10η/(1+2η)+η} ≤ C(η)X^{12η}.
- E₂(X): L(s, χ_D) has a zero with β ≥ 1 − η, |γ| ≤ T_X. χ_D is primitive of conductor |D| ≤ X and T_X ≤ X, so by (I4)
  (Thorner–Zaman, arXiv 2208.11123, Thm 1.2; η ≤ 1/40): #E₂(X) ≤ 10^{88}(10^{421}X^{99})^{η}.
Summing over dyadic X: #E(X) ≤ C X^{99η}log X ≤ X^{100η} = X^ε for X ≥ X₀.

**Explicit formula.** Fix D ∉ E with |D| large, K ↔ ψ, x, H as stated; H₁ := H|D|^{−κ}, T₁ := |D|^{3κ} ≤ T_X. Take smooth
w₋ ≤ 1_{(x, x+H]} ≤ w₊, supported in [x − H₁, x + H + H₁], with |w^{(k)}| ≤ C_kH₁^{−k}; W(s) := ∫w(u)u^{s−1}du, W(1) = H + O(H₁),
|W(ρ)| ≤ C_k H x^{β−1}(4x/(H₁|ρ|))^k and x/H₁ ≤ |D|^{2κ}. For g ∈ G, orthogonality on G/H gives

  Σ_{𝔭 ∤ D, deg 1, ψ(𝔭) = ψ(g)} log N𝔭·w(N𝔭) = (1/3)[ψ_w(1) + ψ_w(χ_D) + 2 Re(ψ̄(g)S_ψ(w))] + O((Hx^{−1/2} + log x)log x),

S_ψ(w) := Σ_𝔞 Λ_F(𝔞)ψ(𝔞)w(N𝔞), ψ_w(χ) := Σ_n Λ(n)χ(n)w(n) (ζ_F = ζ·L(·,χ_D)). By (I6) each of ψ_w(1), ψ_w(χ_D), S_ψ(w) equals
δ·W(1) − Σ_ρ W(ρ) + O((H/x + x^{−1/2})log(|D|x)) (δ = 1 only for ζ).
- Zeros with |γ| > T₁, all three L-functions (≤ C log(|D|(|k|+3)) zeros per unit window, (Z)): with k ≥ 3 + 2/κ,
  Σ ≤ C_k H Σ_{m≥T₁} log(|D|m)(4|D|^{2κ}/m)^k ≤ H|D|^{−1}.
- ζ, |γ| ≤ T₁ ≤ x^{3η/4}: Vinogradov–Korobov and Ingham's N(σ,T) ≪ T^{3(1−σ)}log⁵T give Σ x^{β−1} ≪ exp(−c(log x)^{1/4}).
- χ_D and ψ, |γ| ≤ T₁: D ∉ E₂ ∪ E₁, so β < 1 − η, and there are ≤ CT₁log(|D|T₁) such zeros; their total is
  ≤ CHT₁log(|D|T₁)x^{−η} ≤ Hx^{3κ/A − η}log² x ≤ Hx^{−η/4}log² x ≤ H|D|^{−κ/2} (3κ/A ≤ 3η/4, Aη/4 ≥ κ).
So the left side is W(1)/3 + O(H(|D|^{−κ/2} + e^{−c(log x)^{1/4}})). Each split p gives 𝔭, 𝔭̄ with ψ(𝔭̄) = ψ̄(𝔭): Frob_p = 1 gives
two ideals with ψ = 1, a 3-cycle gives one ideal in each nontrivial coset (Lemma 1). Hence Σ_{Frob_p = 1}log p·w(p) = W(1)/6 + …,
Σ_{3-cycles} = 2W(1)/6 + …, and Σ_{(D/p)=−1} log p·w(p) = ½(ψ_w(1) − ψ_w(χ_D)) + … = 3W(1)/6 + …. Sandwich with w₋, w₊ and
log p = log x + O(H/x) give the formula (H ≥ x^{1−η/4}, so the prime-power and remainder terms are of lower order).

**Fields.** #{K: d_K ∈ E(X)} ≤ Σ_{D∈E(X)}(|G_D[3]| − 1)/2 ≤ #E(X)·max_{|D|≤X}|G_D[3]| ≪_ε X^{ε}·X^{1/3+ε} (Ellenberg–Venkatesh,
IMRN 2007: |Cl(D)[3]| ≪_ε |D|^{1/3+ε}). □

*Where the 1/3 comes from.* The ψ-part alone fails for ≤ #E₁(X) ≤ X^{12η} *fields* (Prop 1 counts characters, i.e. fields). Only
χ_D, one L-function shared by all cubic fields of discriminant D, forces the multiplicity bound; this is PTBW's τ* = 1/3 for S₃.

## 2. Theorem B: what entry 2 itself gives

**Theorem B.** For A > A₀ = (3+√69)/10 let c, κ, δ be as in [R3] Theorem 3 (D < 0) and [R4] Theorem 3⁺ (D > 0). For all
fundamental D, |D| ≤ X, outside a set of size ≪ X^{1−c}, every cubic field K with d_K = D satisfies, for x ≥ |D|^A,
H ∈ [x|D|^{−κ}, x]: #{p ∈ (x, x+H]: Frob_p = 1} ≥ ν_D H/(96 log 2x), #{3-cycles} ≥ ν_D H/(48 log 2x),
#{transpositions} ≥ H/(5 log 2x), with ν_D ≥ c'|D|^{−δ/2}/log|D| as in [R] Theorem 2.

*Proof.* Sum [R] (2.1) (narrow version for D > 0) over the h/3 (h⁺/3) classes of a coset of H and halve (Lemma 1(ii)). For
transpositions use ½(ψ_w(1) − ψ_w(χ_D)) as in Theorem A: D lies outside [R]'s E₂, so the zeros of χ_D (d = D ∈ 𝔉_l) with
|γ| ≤ T₁ have β ≤ 0.995 and cost ≤ H|D|^{−0.005A+κ+ε₁+o(1)} ([R] §2.4(d)); ζ as above. □

*Fields.* The cubic fields over the exceptional D are ≤ X^{1−c}·X^{1/3+ε}, which is o(X) only if c > 1/3. Here c ≤ 0.505δ and
δ ≤ A/4400, so c < 1.4·10^{−4} for A ≤ 1.2 (and c = 10^{−4} for A ≥ 1.7, [R]). The power-saving exceptional set of entry 2,
counted in discriminants, therefore says nothing about almost all cubic fields (3 | h(D) for a positive proportion of D, but the
exceptional D could carry large 3-torsion; Cauchy–Schwarz would need Σ_{|D|≤X}|Cl(D)[3]|² ≪ X^{1+θ} with θ < c ≈ 10^{−4}, while
Ellenberg–Venkatesh give θ = 1/3). Theorem A supersedes Theorem B in every
respect. The reason: a cubic field sees 3 of the h(D) class group characters, so a zero-free box near 1 for one complex character
suffices, and log-free family densities give that for all but X^{O(η)} fields. Entry 2's family large sieve is needed only for all
h(D) classes at once (that is entry 2 itself).

## 3. Remarks

1. *Large A, fields below X^{1/3}* (sketch, not used). For A ≥ 171 the χ_D-part needs no exceptional set: (I4) with
   Deuring–Heilbronn (the N* bound) and (I5) give ≥ ν·W(1)/8 split primes with ν ≫ |D|^{−1/2}log^{−2}|D|; the ψ-part must then
   beat ν, which needs Aη > 1/2, leaving ≤ X^{6/A+ε} exceptional fields. Same order of strength as LOTZ's region plus Linnik.
2. *Non-fundamental d_K = d_F f².* ψ is a primitive cubic ring class character of conductor f of O_{d_K}. Entry 2 and [R] Prop 1
   are proved for fundamental D only ([R] §6: multiplicity of π_χ over D = f²d_K). Near σ = 1 any log-free GL₂ density works, e.g.
   Humphries–Thorner (Compos. 2024) with the convexity mean value, for the dihedral newforms of level |d_K| ≤ X attached to these ψ;
   χ_{d_K} = χ_{d_F}·1_{(·,f)=1} then carries the multiplicity #{K with resolvent F, |d_K| ≤ X} ≪ X^{1/2+ε}, giving all but O(X^{1/2+ε})
   cubic fields. Sketch only; not checked in detail, not claimed.
3. *Numerics* (`least.gp`): all 2716 fundamental D with 3 | h(D), 5 ≤ |D| ≤ 20000 (5564 cubic fields). Max over fields of
   log P/log|D| for the least prime in the class: identity median 0.31, 99% 0.70, max 1.30 (D = −23, P = 59); 3-cycles median
   0.16, 99% 0.34, max 0.39; transpositions median 0.12, 99% 0.36, max 0.52 (D = −231). Consistent with (log|D|)^{O(1)}, far below 1.13.

## 4. Literature (2026-10-21; arXiv API 8 queries `q3.out`, GXL 5 claims, sources read)

| Source | Family | Exceptions (count ≤ X) | Least prime / range | Short intervals |
|---|---|---|---|---|
| Ellenberg–Pierce–Wood, ANT 11 (2017), Thm 2.1 + Rem 2.1 | all cubic fields | ≪ X^{1−δ+ε}, any δ ≤ 2/25 | ≫ X^δ/log X primes ≤ (X/2)^δ of each splitting type | no |
| Pierce–Turnage-Butterbaugh–Wood, Invent. 219 (2020), Thm 3.3(2) | cubic, tame part of d_K squarefree (⊃ fundamental) | ≪ X^{1/3+ε} | asymptotic for x ≥ κ₁exp(κ₂(log log D^{κ₃})^{5/3}(log log log)^{1/3}) | no |
| Thorner–Zaman, arXiv 1909.01338 | removes the strong Artin hypothesis from PTBW | | | |
| Lemke Oliver–Thorner–Zaman, Invent. 235 (2024), Thm 1.1, Cor 3.8 | all cubic fields, 2-dim part | ≪ X^ε | Σ_{p≤x}tr ρ(Frob_p) ≪ x e^{−c√log x} for x ≥ (log D)^{O(1/ε)}; zero-free strip of width ≍ε for heights ≤ D^{O(1)} | routine from the strip |
| Lemke Oliver–Smith, arXiv 2405.08383 | every G | bad Galois fields ≤ Δ^{ε(1+o(1))}; fields over a bad subfield (here the quadratic resolvent) uncontrolled (their Rem 1.11) | H ≥ (log Δ)^{(2+2/ε)[K:F]} | no |
| every field: Ge–Milinovich–Pollack (Math. Z. 2019); Cho–Lemke Oliver–Zaman (arXiv 2512.24963, Thm 2.5); Cho–Kim (PAMS 2021); Burgess | every cubic field | none | for fundamental d_K (D_{K̃} = \|d_K\|³): identity \|d_K\|^{3/2+ε}; non-identity \|d_K\|^{3/4+ε}; 3-cycles \|d_K\|^{9/16+ε} if d_K < 0; transpositions \|d_F\|^{1/(4√e)+ε} | – |
| 2026 half-planes (unrefereed: Liu B♭, OpenAI 7/8); Frei–Sofos–Widmer arXiv 2610.18491; [ArtinHalfPlane](../ArtinHalfPlane/NOTE.md) | every field | none | (C log Δ_{K̃})^{7.99} | – |
| ours, Theorem B (entry 2) | fundamental d_K | O(X^{1−c}) discriminants; field count not controlled | \|d_K\|^{1.1307+ε} | yes |
| ours, Theorem A | fundamental d_K | O(X^ε) discriminants; O(X^{1/3+ε}) fields | \|d_K\|^{A}, any A > 0 | yes (asymptotic) |

2025–26 hits: only Cho–Lemke Oliver–Zaman (Dec 2025, every field) and Frei–Sofos–Widmer (Oct 2026, conditional on the half-planes).
Other hits, not stronger: Cho–Kim arXiv 1601.03012 (averages of least primes over S₃, S₄ fields), An arXiv 2004.14410 (log-free
density, almost-all Chebotarev for a degree-n family), Thorner–Zaman ANT 2017/2019, Kadiri–Ng–Wong PAMS 2019 (exponent 16, all G).
GXL found no almost-all Chebotarev in short intervals; none was found elsewhere either. The short-interval form of Theorem A may
be unstated, but it is a routine consequence of LOTZ's zero-free strip and the classical quadratic-character facts.

**Decision: NO DISCOVERY.** The cubic corollary of entry 2 is weaker than published almost-all bounds in exponent (1.13 against
D^{o(1)} or 0.08), in exceptional set (discriminants, not fields), and for 3-cycles even weaker than the every-field bound. The
strongest version our tools give (Theorem A) is routine and matches PTBW's 1/3 in the field count.

## 5. Dead ends (NOGO candidates)
- Cubic fields from the every-class theorem: needs the exceptional discriminants to carry little 3-torsion (c > 1/3 with
  Ellenberg–Venkatesh); entry 2 has c ≈ 10^{−4}. Use the three characters directly (Theorem A).
- Below X^{1/3} exceptional *fields* at small A: blocked by χ_D, which is shared by all cubic fields of discriminant D; needs
  moments of |Cl(D)[3]| beyond the first, or a cubic-field-weighted large sieve for χ_{d_K} (Davenport–Heilbronn with local
  conditions), not attempted.

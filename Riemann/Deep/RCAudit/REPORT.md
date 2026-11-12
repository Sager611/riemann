AUDIT PASSED: Combos/NOTE.md §1 is correct as stated — Theorem 1 (Σ_χ N(σ;t,χ) ≪_ε |D|^{2σ(1−σ)+ε}T³ for 0.6 ≤ σ ≤ 1, all D < 0), Lemmas 1–5, and, under (L), Theorem 2 and Corollary 3; only the discovery file drops the σ-range, and its claim is false as written for σ < 1/2.

# RCAudit: referee of the ring class zero-density theorem (Combos NOTE §0–1), 2026-11-12

Brief: .scratch/orca/rcaudit-brief.md. Long derivations: [DERIVATION.md](DERIVATION.md). Scripts and outputs: [checks/](checks/).

## Verdict table

| Step | NOTE lines | Inputs | Verdict |
|---|---|---|---|
| (F1) a ↦ a∩O, and a ⊂ mO_K ⟺ a∩O ⊂ mO for (m,f) = 1 | 44–46 | Cox Prop 7.20, 7.22 | correct (mO_K∩O = mO re-derived) |
| (F2) a_Q, N(xa+y(−b+√D)/2) = aQ(x,−y), b ⊂ mO ⟺ m \| gcd(x,y) | 47–51 | Cox Thm 7.7 | correct |
| (F3) h(D) ≤ 3\|D\|^{1/2}log²(3\|D\|) | 52–53 | Cox Thm 7.24, L(1,χ) ≤ log\|d\|+2 | correct |
| (F4) zeros per unit window ≤ c₀log(A(\|t\|+3)²) | 54–55 | IK Prop 5.7 | correct, absolute c₀ |
| (F5) AFE, \|ε\| = 1, \|V_s(y)\| ≤ C_B min(1,(y/(\|u\|+2))^{−B}) | 56–59 | IK Thm 5.3 (G = e^{u²}, X = 1), Prop 5.4 | correct, uniform in u |
| (F6) explicit formula, \|R\| ≤ C(log A+2) | 60–63 | ArtinHalfPlane Lemma 1 | correct as used: the display's proof does not use that lemma's zero-free hypothesis, and (L) supplies it in Thm 2 |
| Lemma 1 primitive lattice points | 67–75 | — | correct; numerics below |
| Lemma 2 weighted large sieve, constant (30/ε)log²(3\|D\|)(N+\|D\|^{1/2}) | 76–97 | L1, F1–F3 | correct; one unstated line: primitivity transfers O_K ↔ O (every m with b ⊂ mO_K is prime to f) |
| Lemma 3 coefficient sums | 98–107 | d_k(mO_K) ≤ d_k(m)², d_kd_l ≤ d_{kl} | correct |
| Lemma 4 mollified second moment ≪ \|D\|^{1/2+ε}X^{1+ε}(\|u\|+3)^{1+ε} | 108–133 | F5, L2, L3 | correct (c_b χ-free inside a conductor group; extension to all of Ĝ is an upper bound; χ̄(a) = χ(ā)) |
| Lemma 5 zero detector | 134–149 | contour shift, \|Γ(v+iu)\| ≤ C₃e^{−\|u\|} on −1/2 ≤ v ≤ −1/10 | correct; C₃ is why σ ≥ 0.6 |
| Thm 1(b) | 160 | F3, F4 | correct |
| Thm 1(a), class II | 166–169 | L4, L5 | correct; cosmetic: L4 at ε₀ gives \|D\|^{1/2+2ε₀+ε₀²}, not \|D\|^{1/2+2ε₀} (take L4 with ε₀/2) |
| Thm 1(a), class I (tail, dyadic j, removing β, k-th power, Hölder, Sobolev, L2+L3) | 171–202 | L2, L3, L5 | correct; exponent recomputed: (σ+ε₀)(2−2σ+5ε₀) ≤ 2σ(1−σ)+6ε₀, ×J ⇒ +10ε₀ |
| Thm 1 conclusion, T³ per window, T⁴ to height T | 151–158, 203 | — | correct for 0.6 ≤ σ ≤ 1 |
| Thm 2 (least prime p(Q) ≪ \|D\|^{2B♭+ε}) | 205–252 | Thm 1, F6, (L) | correct under (L); (Hu) variant 7/4+ε correct |
| Cor 3 (ring class PNT, \|D\| ≤ x^{1/(2B♭)−ε}) | 253–281 | Thm 1, F6, row count | correct under (L) |

## Findings (none breaks a proof)

1. **Statement range (discovery file).** discoveries/ring-class-zero-density.md states the bound with no σ-range. NOTE Thm 1(a) is for 0.6 ≤ σ ≤ 1 only (l.153). For σ < 1/2 the range-free claim is false: zeros are symmetric under ρ ↦ 1−ρ̄, so Σ_{|k|≤T}N_D(σ;k) ≥ ½Σ_χN(T,χ*) ≫ hT log|D| ≫ |D|^{1/2−ε}T, while |D|^{2σ(1−σ)+ε}T⁴ is a power of |D| smaller (σ = 0: |D|^ε). On 1/2 < σ < 0.6 it is unproved in NOTE but a routine extension (DERIVATION §6). Fix: add "for 0.6 ≤ σ ≤ 1" to the Claim.
2. **Cosmetic.** Class II ε-bookkeeping (row above); (F6) citation should say "first display, whose proof does not use the zero-free hypothesis"; Lemma 2's proof should state the primitivity transfer.
3. **Internal supersession (not audited).** RCApps2 NOTE Prop 1(a) claims Σ_χN(σ;t,χ*) ≪ (|D|T)^{2(1−σ)/(3−2σ)+ε}, 0.6 ≤ σ ≤ 1, for the same family. That exponent is smaller than 2σ(1−σ) for every σ in (1/2,1) (σ(3−2σ) > 1); e.g. σ = 3/4: 0.333 vs 0.375. So Combos Thm 1 is no longer the project's best density for this family; the N7/N11 DISCOVERIES line rests on the weaker version.

## Numerical falsification tests (all passed)

- **Lemma 1.** Every reduced form of every D ∈ [−10⁴, −3] (5000 D; h(D) = qfbclassno(D) for all) at N ∈ {⌊|D|^{1/2}⌋, |D|, 4|D|}, plus 1244 (D,N) with N up to |D|² (|D| ≤ 400) and |D|^{1.5} (141 D ≤ 10⁴): max P_Q(N)/bound = 0.62 (D = −3, N = 1).
- **Lemma 2, exact sharp constant.** The best constant is K = h·max_C Σ_{a∈C, Na≤N, (a,f)=1} m(a)^{−1−2ε} (equality in Cauchy–Schwarz at c_a = m(a)^{−1−2ε} on the extremal class), computed exactly from primitive representations by each reduced form (norms prime to f, divided by |O^*|, all m prime to f). Same 5000 D × 3 N and 1244 (D,N): max K/claim = 1.5·10⁻³ (ε = ½), 3.1·10⁻⁴ (ε = 0.1), 3.2·10⁻⁵ (ε = 0.01). Growth: K/((N+|D|^{1/2})log(3|D|)) ≤ 0.33, 0.354, 0.363 (ε = ½, 0.1, 0.01; max at D = −479, −311) in every test; K/(N+|D|^{1/2}) rises 0.67 → 1.8 → 2.6 → 3.3 over the decades of |D| (one log, from h ≍ |D|^{1/2}L(1,χ)), so log² in the claim is generous. The Remark (weight needed) is confirmed: the unweighted constant at N = |D|^{1/2} has median 0.45 → 0.64 → 1.12 over the decades 10–10⁴ (≈ |D|^{1/4} growth).
- **Lemma 4.** S := Σ_{χ≠1}|L(½+iu,χ)M_X(½+iu,χ)|² with PARI 2.19 (primitive Hecke L-functions via bnrinit/bnrconductor/bnrmap, E_χ restored, ring class characters = ray class characters mod f trivial on (Z/f)^*, counts checked against qfbclassno): 123 fundamental and 85 non-fundamental D (f ≤ 60, |D| ≤ 10⁴), X ∈ {1,5,20}, u ∈ {0,10,30}. R := S/(|D|^{1/2}X(|u|+3)): max 3.72 (X = 1, u = 0, D = −9455), 2.0 non-fundamental (D = −5831), 0.12 at X = 20, 0.25 at u = 30. The X = 1, u = 0 case was recomputed exactly by Chowla–Selberg, S = hΣ_C ζ(½,C)² − ζ_K(½)², ζ(½,C) = (1/2)a^{−1/2}[2γ + 2log(√|D|/(8πa)) + 8Σ d(n)cos(πnb/a)K₀(πn√|D|/a)]; it matches PARI to 14 digits (D = −23, −119, −1031, −9455, −9719). Over 182 fundamental D up to 10⁹, R/(L(1,χ_D)²log³|D|) rises to 0.0025 at 10⁴ and stays in 0.0029–0.0050 on 10⁵–10⁹ (a 20% drift over 10⁶–10⁹, the size of a lower-order log term): the excess over |D|^{1/2} is L(1,χ)²log³|D|, no detectable power of |D|. The formula also proves it: |ζ(½,C)| ≤ a^{−1/2}(γ + |log(√|D|/(8πa))| + 0.02) and Σ_C a_C^{−1} ≪ log²|D|.
- Max process memory 76 MB (l12.py), 57 MB (epstein.py), 14 MB (gp); all single-threaded under alarms.

## Literature check

- **Pascadi–Thorner 2508.14888** (v2, 13 Sep 2025; no later version). Thm 1.3 (thm:ZDE, read in the TeX): Σ_{π∈S}N_π(σ,T) ≪ (D_F^{−n(n+1)}Q^{2n+1}T^{[F:Q]n(n+1)+4}|S|⁴)^{(1−σ)/(3−2σ)+o(1)}, with 𝔠_π = D_F^n N𝔮_π Π(3+|μ|^{d(v)}). For n = 1, F = K, S = Ĝ: (|d_K|^{−2}Q³T⁸h⁴)^{…}, Q ≍ |D|, h ≍ |D|^{1/2}: |D|^{3(1−σ)/(3−2σ)} (fundamental D) up to |D|^{5(1−σ)/(3−2σ)}. Theorem 1 is smaller for every σ: 2σ(3−2σ) < 3 ⟺ 4σ²−6σ+3 > 0. σ = 3/4: 0.375 vs 0.50–0.83; σ = 0.9: 0.18 vs 0.25–0.42.
- **Kowalski–Michel 2002, Weiss/Thorner–Zaman, An, Humphries–Thorner** (NOVELTY2 table, not re-read): exponents ≥ 2.75 at σ = 3/4; log-free ones win only within O(1/log|D|) of σ = 1.
- **New since the 20–21 Oct checks.** arXiv listings Oct + Nov 2026, math.NT/CA/CV/SP/PR/FA/CO, two title regexes (zero density, class group, ring class, Hecke, binary/quadratic forms, large sieve, density, dihedral, imaginary quadratic, theta series, ideal class, least prime, Linnik, Epstein): 82 title hits, none a density estimate for class group, ring class or imaginary-quadratic Hecke families. Abstracts read: 2610.04045 (de Faveri, large sieve for n-th order Hecke characters, n ≥ 3), 2611.00924 (Chang–Wang, half-plane B = 0.87446, no density), 2610.15989 (Ye, short-strip densities for tempered GL(n)), 2610.32470 (Chambers, Shintani lifts). Zenodo, three queries ("zero density" "class group"; ring class / class group L-function / class group characters; "binary quadratic form" "least prime"): nothing relevant.
- Verdict: still new, and still stronger in |D| than every published bound for every fixed σ ∈ [0.6, 1) (on (1/2, 1) after the routine extension); method classical: new but routine.

## Dependents (iv)

The uses match: RCApps NOTE (I1) quotes Thm 1 with 0.6 ≤ σ ≤ 1, unit windows, (b) for all σ; Prop 1 uses (I2) = Lemma 4 at X = 1, uniform in u, integrated over |v| ≤ T; Thm 2 §2.4(b) uses (a) on [0.6, σ₀] with Σ_{|k|≤T₁}(|k|+3)³ ≪ T₁⁴ and (b) on [½, 0.6]: all match the statements (σ-range, per-window t-dependence, uniformity in D). RCApps3/4/5 and AllDisc use Lemmas 2, 3, 5 and Thm 1(b); not refereed here.

## Reviewer

Two reviewer runs (decisive arguments pasted inline with NOTE quotes; run 52a559a7 and a follow-up).
- Run 1 on Lemma 2, Lemma 4 and the Theorem 1 exponent (condensed excerpt): S1 (Lemmas 1–2) CORRECT; S2 GAP "missing ε-reparameterization"; S3 GAP "k need not be defined for small |D| (D = −1000, ε = 0.1, j = 1)", "Class II gives |D|^{1/2+2ε₀+ε₀²}", "needs a per-character zero count"; overall GAP/INCOMPLETE, ROUTINE, "repairable completion gaps, not a counterexample".
- Each point is handled by a NOTE line the excerpt omitted: l.130–132 "≤ C(ε)|D|^{15ε}P^{1+3ε} ≤ C(ε)|D|^{1/2+19ε}(X(|u|+3))^{1+5ε} (ε ≤ 1/3). Rename ε."; l.160 "For (a) we may assume |D| ≥ D₀(ε); otherwise use (b)." with l.187 "Since 2^{j−1} > X/4 ≥ |D|^{ε₀/2}"; l.164 "N_D(σ;t) ≤ c₀log(|D|T²)(|𝒳| + 1)" via (F4). Run 2 on these exact lines: G1–G4 all CLOSED; overall CORRECT, ROUTINE (G4: Lemma 4 at δ = ε₀/(1+ε₀) gives exactly |D|^{1/2+2ε₀}).

RH relevance: N7, N11. Theorem 1 is the family zero-density tool behind the least-prime, Sato–Tate and CM least-prime lines; it moves no RH node by itself (records only), but its audit removes an unrefereed load-bearing input.

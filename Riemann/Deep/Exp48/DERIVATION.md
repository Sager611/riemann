# exp48 DERIVATION: the exponent-4 and exponent-8 lists of Elsenhans–Klüners–Nicolae are complete

Agent exp48, 2026-11-06. Every numbered constant below is certified in `cert.py` (output `logs/cert_output.txt`), checked
by `scan.c` (unconditional class-group scan), `conductor.py` or `spot.gp`; the check is named in brackets.

## 0. Notation and inputs

- D < 0 a discriminant, d = |D|, χ = χ_D, h = h(D) the form class number, Cl(D) the class group of primitive positive
  definite forms of discriminant D (= Pic of the order), E(D) its exponent, ω(D) the number of distinct primes dividing D,
  d_k = p_1 ⋯ p_k the product of the first k primes.
- (L78) Littlewood78 Cor. 2 (`Riemann/Deep/Littlewood78/REPORT.md` l.1 and l.35; `discoveries/littlewood-half-plane.md`;
  paper proof, two reviews, not formalized): for every fundamental D < −4,
  h(D) ≥ L78(d) := 0.0350·√d / max(3, log log d).
  It comes from |L(1, χ)| ≥ 0.11/max(3, log log q) for every primitive χ mod q ≥ 3 and h = √d·L(1, χ_D)/π (d > 4, w = 2);
  0.11/π = 0.035014 ≥ 0.0350 [cert (0a)]. Input: the Lean-checked half-plane L(s, χ) ≠ 0 for Re s > B = 0.87338196661320752
  (TREE.md N7). No GRH. Fallback with OpenAI's 7/8: 0.105/π ≥ 0.0334 [cert (5a)].
- EKN = Elsenhans, Klüners, Nicolae, arXiv 1803.02056 (only version, TeX `class4.tex` of 6 Mar 2018; Acta Arith. 193 (2020)
  217–233). Line numbers refer to `class4.tex`.
  - Thm 1 (l.149–175): exactly 9, 56, 17, 203, 27, 432, 33, 778 fields of exponent 1, …, 8 with |D| ≤ 3.1·10^20; largest |D|
    163, 5460, 4027, 435435, 37363, 5761140, 118843, 430950520; the nine with |D| > 10^7 are listed (l.172–173).
  - Thm 2 (l.176–183): "Assuming the non-existence of Siegel-zeros, our computation found all fields with exponents 2,4,8.
    Without assuming any unproven statement we can conclude that there is at most one missing field with exponent 2,4 or 8."
  - Lemma 1 (l.214–217): p split, p^c < |D|/4 ⇒ the order of a prime above p exceeds c.
  - Lemma Siegel1 (l.306–310, Tatuzawa), Lemma Siegel (l.312–321): if L(s, χ) ≠ 0 on [1 − 1/(4 log|D|), 1) and |D| > e^{11.2},
    then h > 0.655/(πe)·√|D|/log|D| (0.655/(πe) = 0.0767003 [cert (0c)]); Lemma Siegel_one (l.323–336); Lemma Chen (l.339–352):
    for A ≥ 10^6 and all |D| ≥ A with at most one exception, h ≥ min(√|D|/(7.732π log|D|), m·1.5·10^6·√|D|/(πe^m log|D|)),
    m = log|D|/log A — inequality (b_ch); 1/(7.732π) = 0.0411679 [cert (0b)].
  - Lemma lem6 (l.375–387): E(D) = 2^r ⇒ h ≤ 2^{r(ω−1)}; Lemma lem7 (l.389–399): |D| ≥ d_ω.
  - Thm upp_bound (l.401–409, proof l.410–468) and Example (l.474–476): N_2 = 11, N_4 = 24, N_8 = 58, d_N ≤ 2.01·10^11, 2.38·10^34, 3.17·10^110.
  - §6 (l.484–731): Rédei-matrix search, for each k = ω(D) ≤ N_c, of all D with |D| > 10^6 (l.712–715); Algorithm algo Step 1
    (l.695–696): B_ℓ such that min(b_ch) > 2^ℓ c^{k−1−ℓ} for all |D| > B_ℓ; B_s is used in Check Step 6 (l.641) and NextTuple
    Step 2 (l.676). l.705–710: "missing examples have the property that the bounds of Lemma Chen are wrong and this can happen
    at most one time." Running time 60 h with Lemma Chen, 17 h with the Siegel bound (l.718–728).
  - §7 (l.734–890): direct search without L-functions (smallest split prime ≤ 197: 268 fields with E | 4 and 778 with E = 8,
    largest 430950520; sieve to 3.14·10^20 for no split prime ≤ 193: 1,002,279 candidates, none with E ≤ 100).
- ARSX = Aiken, Ross, Shvydkoy, Xue, arXiv 2605.26020 (TeX `boundary_CM_pts.tex`): Thm all-CMpts-on-bdry (l.231–250);
  Thm exp-tables (l.261–270); Lemma property-of-f (l.603–627: L(f) | 12·|O_K^×|·E(D), from Fan–Pollack arXiv 2503.07801
  Prop. 2.3); Lemma f-divides (l.631–637); algorithm (l.645–665); l.666–668: "this is only reason we need such conditions; the
  rest of the algorithm is based on the unconditional results about conductors". Tables (l.676–800): E = 1: 13 discriminants
  (largest |D| 163), E = 2: 88 (7392), E = 4: 485 (887040), E = 8: 2329 (1723802080); data and Sage code at
  github.com/erick-ross/research-data, `class-groups-small-exp` (`find_fund_discs` brute-forces fundamental discriminants up
  to EKN's largest |D|, so it relies on EKN's completeness).

## 1. Where a Siegel-type input enters (exhaustive, from the two TeX sources)

| # | place | input | role |
|---|---|---|---|
| S1 | EKN Thm upp_bound (l.401–468) | Lemma Siegel (conditional) or Lemma Siegel_one with A = d_N (one exception) | caps ω(D): no D with \|D\| ≥ d_{N_c} and E(D) = c |
| S2 | EKN §6 Algorithm algo Step 1 (l.695–696), used at l.641, l.676 | (b_ch) of Lemma Chen (60 h run) or Lemma Siegel (17 h run) | caps \|D\| for each k = ω(D) and each lower bound s of the Rédei rank: \|D\| ≤ B_s |
| S3 | EKN Thm Bach_Bound (l.251–260), Table l.264–280 | ERH | only the ERH statements (exponents 3, 5, 6, 7); not used here |
| S4 | ARSX algorithm l.649–651, l.666–668 | EKN's list of fundamental D_0 with E(D_0) \| E | the only conditional step of ARSX |

Not Siegel-type: EKN Lemma 1, genus theory (lem6), lem7, Rédei theory (Thm redei, Sueyoshi's rank bound), Check Steps 1–5
and 7, §7 and the small fields (tables, Magma, l.494–499); ARSX's conductor lemmas. EKN's single possible missing field is
the exception of Lemma Chen in S2 (l.705–710); the unconditional half of Thm upp_bound has its own Tatuzawa exception
(Lemma Siegel_one at A = d_58), and EKN do not argue that the two coincide. Both disappear below. Side remark: with Chen's
constant in place of 0.655/(πe), the three conditions of Thm upp_bound give N = 12, 24, 59 rather than 11, 24, 58 [cert (4)];
EKN's sentence at l.602–604 ("Using Lemma Chen … we get 11, 24, and 58") uses the Siegel constant. Immaterial here.

## 2. Replacing each Siegel-type input by L78

**(R1) for S2.** For every d ≥ 35: L78(d) > √d/(7.732π log d) ≥ min(b_ch)(d) [cert (1a)–(1c)].
Proof. If log log d ≤ 3 the claim is 0.0350·log d > 3/(7.732π) = 0.123504; the left side increases in d and at d = 35 equals
0.124437 (margin 9.34·10^−4; at d = 34 it fails). If u = log d ≥ e^3, the claim is 0.0350·u/log u > 1/(7.732π); u/log u
increases for u ≥ e, and at u = e^3 the left side is 0.0350·e^3/3 = 0.2343 > 0.0412. The second term of the minimum is never
needed. For d ≥ 10^6 the ratio L78/Chen1 is at least 3.915 (at d = 10^6) [cert (1d)].
Consequence: every fundamental D with |D| > 10^6 satisfies (b_ch); EKN's §6 search (Chen version) misses no field.

**(R2) for the faster run and for Thm upp_bound's first half.** For d ≥ e^{6.6}: L78(d) > 0.655/(πe)·√d/log d
[cert (2a)–(2c): 0.0350·6.6 = 0.231 > 3·0.0767003 = 0.230101; 0.0350·e^3/3 > 0.0767003]. So the conclusion of Lemma Siegel
holds for every |D| > e^{11.2} with no hypothesis on zeros. (Its hypothesis also follows directly from the half-plane:
1 − 1/(4 log|D|) > B once log|D| > 1/(4(1 − B)) = 1.974 [cert (2d)]; this is the route of OpenAI's 5 Oct note for exponent 2, but it
rests on EKN's reading of Tatuzawa's Lemma 9, which (R2) does not need.)

**(R3) for S1.** Lemma A: if E(D) | 2^r (r ≥ 1) then h(D) ≤ 2^{r(ω(D)−1)} — lem6 with a_r ≥ 0 allowed: Cl(D) ≅
⊕_{j ≤ r} (Z/2^j)^{a_j}, Σ a_j = ω − 1 (genus theory, fundamental D), h = 2^{Σ j a_j}.
L78 is increasing in d [cert (3'): on log log d < 3 it is 0.0350√d/3; for u = log d > e^3, d/du(u/2 − log log u) =
1/2 − 1/(u log u) > 0]. With |D| ≥ d_ω (lem7): L78(d_ω) ≤ L78(|D|) ≤ h ≤ c^{ω−1}. [cert (3)]: L78(d_k) > c^{k−1} for every
k ≥ K_c with K_2 = 10, K_4 = 21, K_8 = 53, checked exactly for k ≤ 399; for k ≥ 399 by induction: p_{k+1} < d_k (Euclid,
[cert (3')]) gives log d_{k+1} < 2 log d_k, so max(3, log log d_{k+1})/max(3, log log d_k) < 1 + log 2/3, and
L78(d_{k+1})/L78(d_k) = √p_{k+1}·(that ratio)^{−1} > √p_{K_c+1}/(1 + log 2/3) > c (p_{K_c+1} = 31, 79, 251).
Hence E(D) | 8 ⇒ ω(D) ≤ 52, E(D) | 4 ⇒ ω(D) ≤ 20, E(D) | 2 ⇒ ω(D) ≤ 9: inside EKN's search ranges k ≤ 58, 24, 11.
The bound fails at k = K_c − 1 [cert (3)], so K_c is the best this argument gives. With the 7/8 fallback the same K_c
(10, 21, 53) [cert (5c)] and (R1) holds for d ≥ 41 [cert (5b)].

## 3. Theorem A (fundamental discriminants)

Prior proof (found in the novelty pass, NOVELTY.md): M. Jaswal, Zenodo 23251444 (9 Oct 2026; v5 23547020, 21 Oct),
§5.1, Thm 5.2, under (H) = no real zero in (7/8, 1) for real primitive χ, by route (R2) and Tatuzawa. Theorem A below is
his Theorem 5.2 with (H) supplied by the half-plane; our proof replaces Tatuzawa's localization by (R1)/(R3).

Theorem A. Let D < 0 be a fundamental discriminant. Then E(D) = 4 if and only if D is one of EKN's 203 discriminants
(largest |D| = 435435), and E(D) = 8 if and only if D is one of EKN's 778 (largest |D| = 430950520). With E(D) | 8 there are
exactly 9 + 56 + 203 + 778 = 1046 fields. Hypotheses: the half-plane (through L78) and the correctness of EKN's §6
computation for |D| > X1 (X1 = our scan bound, §4).

Proof. D = −3, −4 have h = 1. Let D < −4, E(D) | 8. If |D| ≤ X1 the scan (§4) decides it. If |D| > X1 ≥ 10^6:
ω(D) ≤ 52 by (R3), and D satisfies (b_ch) by (R1); EKN's §6 search with c = 8 and k = ω(D) ≤ 58 returns every D with
|D| > 10^6, k prime factors, exponent dividing 8 and (b_ch) (l.705–710: a missed field violates (b_ch)). So D is on EKN's
list, whose largest |D| is 430950520 < X1 (§4): a contradiction. Hence every such D has |D| ≤ X1 and is one of the scan's
1046, which are EKN's. The exponent-4 statement is the case E(D) | 4 of the same search (or the
c = 4 search with ω ≤ 20). □
The exponent-2 statement (56 fields, largest 5460) follows in the same way; for D ≡ 0 mod 4 it is Euler's idoneal-number
list, completed by OpenAI's 5 Oct note through EKN Thm 2 (paper2.pdf §1, pdftotext l.89–94). Not claimed here.

## 4. Independent unconditional scan of all discriminants with |D| ≤ X1 (scan.c)

For every D < 0 with D ≡ 0, 1 mod 4 and 3 ≤ |D| ≤ X1 (fundamental or not):
1. For the split primes q (q ∤ D, (D/q) = 1) in increasing order: if 4q^8 < |D| the class of the prime form (q, b, ·) has
   order > 8 (Lemma 1; the proof holds for any order since q ∤ conductor: 𝔮^8 = (α), α = (x + y√D)/2, y = 0 would force
   𝔮^8 = (q^4) = 𝔮^4 𝔮̄^4); reject. Otherwise compute f_q^8 by three exact squarings (Cohen Alg. 5.4.7 + reduction, 128-bit);
   if it is not principal, reject (an explicit element of order not dividing 8). At most 6 split primes are tested.
2. Survivors: enumerate every reduced primitive form (|b| ≤ a ≤ c, b ≥ 0 if |b| = a or a = c, gcd(a, b, c) = 1);
   h = their number (exact). E | 8 iff every f^8 is principal; n2 = #{f: f^2 = 1}, n4 = #{f: f^4 = 1} fix the structure.
No step uses an L-function or GRH. X1 = 10^10 (four runs: 3..10^7, 10^7..4·10^8, 4·10^8..2·10^9, 2·10^9..10^10;
`logs/scan_stats.txt`, list `logs/scan_E_divides_8_upto_1e10.txt`). Results: up to X1 the discriminants with E(D) = 1, 2, 4, 8 number
13, 88, 485, 2329 (largest |D| 163, 7392, 887040, 1723802080); the fundamental ones 9, 56, 203, 778 (largest 163, 5460,
435435, 430950520). This equals ARSX's Tables 1, 2, 4, 8 discriminant by discriminant and EKN's counts and maxima
(`cmp.py` → `logs/cmp_output.txt`). No discriminant with E(D) | 8 lies in (2·10^9, 10^10]. Survivors of step 1 that failed the
full check: 3 (all between 10^7 and 4·10^8). The 37 odd-or-(−4) entries with E | 2 equal the table of ARSX Thm
all-CMpts-on-bdry. (Jaswal's Prop. 5.1 did the fundamental part of this up to 2·10^9; see NOVELTY.md.)

## 5. Theorem B (ARSX's tables over all discriminants; boundary CM points)

Theorem B. ARSX's Tables 2, 4 and 8 (88, 485, 2329 discriminants) are complete without any hypothesis on zeros (given the
half-plane through L78 and EKN's §6 computation as in Theorem A); Thm exp-tables items (2) and (3) hold unconditionally.
Proof. D = D_0 f^2, E(D) | E ∈ {2, 4, 8}. Cl(D) ↠ Cl(D_0) (Cox Cor. 7.17), so E(D_0) | E and D_0 is one of the 65, 268,
1046 fundamental discriminants of Theorem A (all with |D_0| ≤ 430950520). Lemma property-of-f (algebraic) gives L(f) | 12uE,
u = |O_{K}^×| ∈ {6, 4, 2}, L(f) = lcm_{p^k ‖ f} p^{k−1}(p − (D_0/p)); this leaves finitely many f. `conductor.py` enumerates them by depth-first search over prime powers p^r with
L(p^r) | 12uE and cuts a branch as soon as a divisor f' of f is certified bad: |D_0 f'^2| ≤ X1 and D_0 f'^2 is not on the
scan list with exponent dividing E (then Cl(D_0 f^2) ↠ Cl(D_0 f'^2) has exponent not dividing E). With X1 = 10^10 the leaves reached
are 101, 586, 3879 candidates for E = 2, 4, 8 (f = 1 included; without cuts D_0 = −3 alone has 786,432 candidates for
E = 8, matching ARSX's README "0.8M"). For |D_0 f^2| ≤ X1 the scan decides. The 964 leaves with |D_0 f^2| > X1 (all for
E = 8) are excluded by a divisor f' | f: 433 through the scan list as above, 531 (D_0 large, so |D_0 f'^2| > X1 already
for small f') by an exact test on D_0 f'^2: Lemma 1 (291 cases) or a prime form f_q with f_q^8 non-principal (240 cases).
No candidate is left (`logs/conductor_output.txt`; with X1 = 2·10^9 the same holds: 4025 (= 2 + 4023) leaves beyond X1, 1936 decided
by the exact test, `logs/conductor_output_X1_2e9.txt`). □
(`conductor.py` self-test: its exclusion test never fires on the 3602 known pairs (D, E) with E(D) | E, and it detects an
order-8 element in all 2329 exponent-8 discriminants.)
Corollary. ARSX Thm all-CMpts-on-bdry (l.231–250) holds unconditionally: its proof (l.528–555) needs only the E | 2 tables.
That part rests on the exponent-2 classification (OpenAI's 5 Oct note / EKN Thm 2 for E = 2) plus ARSX's conductor step; it is
a routine consequence and is listed here only for completeness.
Class-number route (brief item 3): h(D_0 f^2) = h(D_0)·f·∏_{p|f}(1 − (D_0/p)/p)/[O_K^× : O^×] ≥ φ(f)/3 together with
h(D) ≤ 8^{μ−1}, μ ≤ ω(D) + 1 (genus characters), also bounds f, but the bound (φ(f) ≤ 3·8^{ω(D_0)+ω(f)}) is far weaker than
ARSX's divisibility lemma and is not needed.

## 6. Checks

- Siegel-type inputs: the table in §1 is every occurrence of Lemmas Siegel1/Siegel/Siegel_one/Chen and of B_ℓ, B_s in EKN
  (`rg -n 'Siegel|Chen|B_' class4.tex`) and of "Siegel"/"conditional" in ARSX.
- d = 3, 4 (units 6, 4): h = 1, E = 1; L78 is applied only for D < −4, where w = 2 and h = √d L(1, χ)/π.
- D ≡ 0 mod 4 and mod 8: L78 holds for every fundamental discriminant (χ_D primitive of conductor |D|); EKN's search runs the
  four 2-parts {1, −4, −8, 8} (Algorithm algo Steps 2–5); the scan treats every D ≡ 0, 1 mod 4.
- Small d (below 10^6): EKN used tables; the scan covers them unconditionally (exact h by counting forms).
- No GRH: L78 uses the half-plane only; the scan and conductor.py use exact form arithmetic only; PARI's quadclassunit
  (GRH-conditional in general) appears only in spot.gp, where bnfcertify certifies the fundamental fields.

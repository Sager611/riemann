# NOVELTY — BV without RH (agent bvirr, 2026-11-06, 10:15–10:55 UTC)

Claims checked: (A) Breuillard–Varjú Thm 2 (both parts) without RH, for every coefficient law; (B) BV Cor 3 for every law;
(C) Littlewood ±1 polynomials irreducible with probability → 1 (rate 1 − O(1/n), n even; 1 − 2√(2/(πn)) + O(1/n), n odd);
(D) Gal ⊇ A_n with high probability for every law (0/1 and ±1 included); (E) Blachar–Breuillard–Kozma's Sym(n) for every μ
with H_2(μ) > 12; (F) Varjú–Xu Thm 1.1, and Cor 1.3 for f(d) > 2^{7.9π(d)}; (G) Bary-Soroker–Garzoni–Sodin Theorems 1–2;
(H) Eberhard 2008.01223 Thm 1.3, Ferber–Jain–Sah–Sawhney 2106.04049 Thm 1.1, 2304.12075 Cor 1.2, all without ERH.

## Searches
- arXiv API (`.scratch/axq.py`, 33 queries): au:Breuillard, au:Varju, au:Kogler, au:"Bary-Soroker", au:Kozma, au:Poonen,
  au:Hokken, au:Garzoni, au:Sodin, au:Eberhard, au:Sawhney (2025-06 onwards); abs: "random polynomial(s)" AND irreducible;
  Odlyzko AND Poonen; Littlewood AND irreducible AND polynomial; "Galois group" AND "random polynomial"; "characteristic
  polynomial" AND irreducible (AND random / "Galois group" / "extended Riemann hypothesis"); "Riemann hypothesis" AND
  irreducible AND random; "Riemann hypothesis" AND "Galois group" AND random; Cureton AND Thorner; "Dedekind zeta" AND
  "zero-free" AND polynomial; and, without date limit from 2023-06: Littlewood AND irreducible AND polynomial, "Littlewood
  polynomials" AND irreducible, "random polynomial" AND irreducible AND coefficients, ti:irreducibility AND ti:random.
- arXiv monthly listings 2026-10 and 2026-11 (`.scratch/lst.py`), math.NT, math.PR, math.CO, math.GR (the math.GR 2026-11 fetch
  failed), title regex `random polynomial|irreducib|Galois group|Littlewood polynomial|Odlyzko|random matri.*polynomial|
  characteristic polynomial`: 33 hits, none on irreducibility of random polynomials or random matrices (closest: 2610.22764
  trinomials with small coefficients, 2610.10614 Laguerre polynomials).
- Hexagon API: all 716 works (2609.00001–2611.00127), titles and abstracts, regex `random polynom|irreducib|Breuillard|Varj|
  Odlyzko|Littlewood polynom|Galois group|Bernoulli convol|Weinberger|Konyagin|Bary-Soroker|characteristic polynomial`:
  38 hits, none relevant (2611.00087 Laguerre polynomials; 2610.00483 real roots of heat-evolved random polynomials).
- Cureton–Thorner's own TeX (Hexagon 2610.00200v2 source, l.150–180): no mention of Breuillard–Varjú, random polynomials or
  irreducibility; their stated consequences are the prime ideal theorem with an x^θ error and the Hecke corollary.
- Zenodo API, 9 queries ("Breuillard", "Odlyzko-Poonen", Odlyzko AND irreducible, "Littlewood polynomials" AND irreducible,
  "Cureton" AND "Thorner", "quasi-Riemann" AND polynomial, irreducibility AND "random polynomial", "Galois group" AND
  "random polynomials", "random polynomials"): nothing relevant. HAL API, 4 queries: nothing.
- Papers that use Cureton–Thorner, read in TeX: Frei–Sofos–Widmer 2610.18491 (QRH Chebotarev; class groups), Baldwin–Sivek
  2611.06705 (Kuperberg's coNP), 2610.37179 (sparse irreducibility testing in PP through a Weinberger-type root count with CT
  and FSW; it notes that any fixed δ > 0 suffices, l.1073, but treats deterministic sparse polynomials and does not cite BV).
  None mentions BV or random polynomials.
- GXL Evidence API (3 claims: unconditional irreducibility of i.i.d. random polynomials; Littlewood polynomials irreducible
  with probability → 1; a uniform zero-free half-plane makes Breuillard–Varjú unconditional): supporting passages only from
  BV, BSHKP, BKK, Kogler, FJSS and older work; nothing states the claims without RH beyond Kogler's 0/1 case.
- Semantic Scholar citations: BV (53), Eberhard (31), FJSS (24). All 13 citing works dated 2025-09 or later were read in TeX
  (2609.09354, 2608.24586, 2608.06576, 2601.16295, 2510.26191, 2609.15813, 2603.15856, 2511.08333, 2608.27624, 2601.09704,
  2510.16225, 2608.22297, 2610.22764, besides the ones above): none removes the hypothesis. Three say in October 2026 that it is
  still needed: 2610.11047 (8 Oct, dgs.tex l.164–166: irreducibility for random symmetric matrices "known only ... under the
  extended Riemann hypothesis, by FJSS"), 2609.15813 (l.360: Eberhard and FJSS "both assuming the extended Riemann
  hypothesis"), 2608.27624 (l.1392: FJSS "conditional on the extended Riemann hypothesis").
- Our files: `rg -i 'Breuillard|1810\.13360|Odlyzko.Poonen|Littlewood polynomial'` outside .scratch finds only
  Harvest4/HARVEST.md (the source of this brief) and R386 (flat Littlewood polynomials, unrelated). No NOVELTY file in
  Riemann/Deep covers BV.

## Closest prior work and what it leaves open
| work | date | proves | leaves open |
|---|---|---|---|
| BV 1810.13360 (Acta Math. 223, 2019) | 2018 | Thm 2, Cor 3, Cor 4 under RH for Dedekind zeta | unconditional versions |
| Bary-Soroker–Koukoulopoulos–Kozma 2007.14567 (Invent. Math. 233, 2023) | 2020 | unconditional for μ uniform mod 4 primes (e.g. ≥ 35 consecutive integers), Alt included | 0/1, ±1, general μ |
| Bary-Soroker–Hokken–Kozma–Poonen 2308.04878 | 2023 | Littlewood: limsup = 1 | lim = 1 |
| Kogler 2609.26771 (current source dated 15 Oct 2026; Lean) | Sep 2026 | 0/1 model unconditional: P(reducible) = √(2/(πn)) + O(1/n); P = ΦI with probability 1 − Ce^{−cn/(log n)^4} | other laws ("future work", introduction l.60); Galois groups |
| Blachar–Breuillard–Kozma 2609.40123 | 30 Sep 2026 | P(disc square) ≤ C/log n for H_2(μ) > 12; Sym(n) for uniform μ | general μ needs GRH (l.189) |
| Varjú–Xu 2511.04240; BGS 2502.17218; Eberhard 2008.01223; FJSS 2106.04049; 2304.12075 | 2020–2025 | under RH/ERH | unconditional |
| 2402.06643 | Jan 2024 | BKK's rate improved to 1 − O(n^{−1/2}) under near-uniformity mod 12 primes | 0/1, ±1, general μ |
| Hokken–Koukoulopoulos 2510.18857 | Oct 2025 | reciprocal model, unconditional by BKK's method | (no RH used; not a descendant) |
| 2511.10359 | Nov 2025 | completely multiplicative ±1: limsup = 1 | lim (VX Thm 1.1 covers VX's two models) |

## Verdict
New as far as we can see on 2026-11-06, except the 0/1 irreducibility statements (BV Thm 1 and Cor 4), which Kogler proved
unconditionally in September 2026 by other means. The method is a substitution of Cureton–Thorner's strip into BV's proof; the
novelty is in the statements. Pre-emption risk: medium. Kogler announces general measures; Breuillard and Varjú, who helped
Kogler, may notice the strip; 2610.37179 and Baldwin–Sivek apply the same strip to Weinberger-type counts.

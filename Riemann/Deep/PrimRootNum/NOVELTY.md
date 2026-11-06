# Novelty checks (primroot, 2026-11-06)

Question: is any numeric bound g*(p) ≤ (c log p)^E (or any explicit polylogarithmic or p^{o(1)} bound) for the least prime primitive
root, valid for every prime, already published or posted? Answer: no. Details:

- **arXiv monthly listings** (lst.py, this and last month: 2026-11, 2026-10). math.NT, regex 'primitive root|least prime primitive|Artin':
  7 title hits, none relevant (2610.08569 Johnston–Naik, least prime with a GIVEN primitive root, a different question; 2610.05914
  conflict-avoiding codes; Artin–Springer, Artin–Davenport, Artin L at s = 1, Artin characters, Artin–Schreier). math.NT + math.CO,
  'primitive root|Vinogradov|nonresidue|non-residue': 3 hits, none relevant (2610.05515 Vinogradov subsystem).
- **Zenodo** (lst.py --zenodo): '"least prime primitive root"' 7 records, '"prime primitive root"' 11 records, 'primitive root explicit bound'
  (noise). The relevant ones are the October notes by "dr vengy", all conditional on OpenAI's 7/8 half-plane (QRH), all read (description
  + pdftotext, rg for constants; PDFs deleted):
  - 23341593 (13 Oct): Martin Prop. 3 ⇒ g*(p) ≪ (log p)^{24}; "No numerical value of C is asserted."
  - 23352868 (14 Oct): character sums ⇒ g*(p) ≪_ε p^ε; Martin route ⇒ (log p)^{24}; no constants.
  - 23394711 (15 Oct), 23423890 (16 Oct, with a Lean "Level1" file of the substitution): same deduction, ≪ only; 23423890 quotes
    Martin on Math.SE agreeing with the deduction.
  - 23405184 (15 Oct): N_p(p − 1) ~ (φ(p−1)/(p−1)) π(p − 1); "no explicit cutoff".
  - 23621868 (24 Oct): g*(p) ≪ (log p)^{24} = o(p); "The implied constant in (3) is not made explicit here."
  - 23623384 (24 Oct): explicit g*(p) < p for p ≥ 10^{41} (then 5·10^{38}, 1.4·10^{33}) under QRH: plain Vinogradov sieve at x = p with
    Robin's bound for 2^{ω}, explicit smoothed character sums; no polylogarithmic numeric bound. Our Corollary E(i) (p ≥ 10^{41.032})
    reproduces this threshold by a different route (now unconditional given Bebel's half-plane, which implies 7/8); not claimed as new.
  - Also listed: 23421906, 23420395 (Erdős 985 for special families, elementary), 23442654 (non-residues), 23116412 (unrelated).
- **Martin 1997** (arXiv math/9807104, TeX read): Prop. 3 g*(q) ≪_σ (ω(φ(q))² log ω · c₀^{−1} log q)^{1/(1−σ)} under a zero-free rectangle;
  the shifted-sieve lemma uses Iwaniec's Lemma 2 asymptotically (O-constants). No numbers.
- **GXL evidence API** (2 calls): (1) "explicit unconditional upper bound for the least prime primitive root for every prime":
  supports only for g(p) (p^{5/8} all p > 3, 2206.14193; p^{0.68}, 1803.11061; g^□(p) < p^{0.88}, 1701.05980; √p − 2 for p > 1.77·10^{54},
  2608.20431) and GRH (ĝ(p) ≤ ((8/5)(2^{ω(p−1)} − 1) log p)² for p ≥ 10⁹ and ĝ(p) < √p − 2 for p > 2791, McGown–Treviño–Trudgian
  1508.05182); refutes: "It is not known whether ĝ(p) < p for all p, or even for all sufficiently large p" (Cohen–Trudgian 1602.02440),
  Martin 1998 ("no uniform upper bounds for g*(p) that improve upon … the least prime in an arithmetic progression").
  (2) "a zero-free half-plane implies an explicit polylog bound for g*(p)": 3 supports, none explicit (GRH 1508.05182; Fan–Pollack;
  Martin's almost-prime paper).
  Carella's arXiv notes (1806.01150, 1709.01172, 1405.0161) claim unconditional g*(p) = O((log p)^{1+ε}) or p^{5/log log p}: unrefereed,
  no constants, and contradicted by the field's statements above (2016–2024 papers still give p^{1/4+ε} as the best unconditional bound
  for g(p)); not treated as prior art.
- **Hexagon** (hexagonmath.org search 'primitive root', 'primitive roots', 'Vinogradov', 'Artin'): Artin-type positive lower density
  (every nonsquare a ≠ −1, elliptic curves with cyclic reduction), Pólya–Vinogradov with log log q (2610.00256); nothing on g*(p).
- **OpenAI** "Primitive roots for every admissible integer base" (4 Oct; our audit Riemann/Deep/PrimRoot/REPORT.md): Artin-type
  (a fixed base is a primitive root for infinitely many p), a different question. Riemann/Deep/ArtinHalfPlane/REPORT.md: Linnik 7/3,
  Chebotarev; no g*(p).
- **Our repo** (rg 'primitive root' over *.md; DISCOVERIES, COROLLARIES, TREE, NOGO; ExplProg/UnitCover NOVELTY files): no numeric
  g*(p) bound; ExplProg's Theorem 2 gives only g*(p) ≤ (p log p)^{7.898} (least prime in the class of a primitive root < p).
- **Verdict.** New: (a) the first numeric bound for g*(p) valid for every prime that is polylogarithmic, (1.956 log p)^{3.4/(1−B)} ≤ (1.956 log p)^{26.86};
  (b) (0.814 log p)^{3/(1−B)} ≤ (0.814 log p)^{23.70} for every p ≥ 5 with ω(p − 1) ≤ 10⁴; (c) the per-ω bounds (κ_ω log p)^{7.898};
  (d) g*(p) < √p − 2 for p ≥ 10^{97.547} without GRH. Not new: the asymptotic shape (Martin 1997 + the half-plane, Zenodo notes);
  g*(p) < p for p > 10^{41} (Zenodo 23623384, QRH-conditional).

# DetFactor: novelty scan (6 Nov 2026, 11:30–13:00 UTC)

Verdict: the three goals are KNOWN (OpenAI 30 Sep and 4 Oct 2026; Hexagon 2610.00178, 11 Oct; Zenodo 23205034, 6 Oct).
Not found anywhere: that OpenAI's factoring theorem needs only Cureton–Thorner's Lean-checked strip, a numeric bound for
its auxiliary primes, or its exponent with that input (Prop. C, Thm B: routine).

## Pre-empting sources (read in TeX)
- **OpenAI, "The Quasi-Riemann Hypothesis" (github.com/openai/math, preprints/The-Quasi-Riemann-Hypothesis-September-30-2026,
  build/paper.tex)**, Corollary cor:nonresidues-square-roots (l. 320–362): n(p) ≤ C(log p)^A (A = 32 via BIMS weak GRH with
  ε = 7/16) and a deterministic polynomial-time square-root algorithm (Tonelli–Shanks; "its stopping rule does not require
  knowing C"). → task (1) KNOWN for square roots.
- **OpenAI, "Deterministic Polynomial Factorization over Prime Fields" (4 Oct 2026, CONTENTS.md entry 142; build/main.tex,
  sections 00–10)**: Theorem thm:factorization (00-introduction l. 27–38): complete factorization of degree-n polynomials over
  F_p in O(((n + 1)⌈log₂ p⌉)^{10^{12}}) deterministic bit operations, no GRH. Analytic input: companion "Primitive roots for
  every admissible integer base" (4 Oct, entry 029) Theorem 1.2, Hecke strip of width 10^{−6} over cyclotomic fields ⊇ μ_12
  (10-analytic l. 9–18); used only in Section 10 (rg of all sections, DERIVATION §3). No Lean link in CONTENTS for 029 or 142;
  OpenAI's Lean comparator challenges have Hecke 7/8 only for ℚ(√−3) (lean/ComparatorChallenges/HeckeSevenEighths.lean, sorry
  in the challenge statement as usual). Its table section builds q-primary generators of F_p(ζ_q) (task 2 for these fields).
  → task (3) KNOWN (claimed; much stronger than unconditional Rónyai/Evdokimov).
- **Hexagon 2610.00178 (11 Oct 2026), "Unconditional implications of the Quasi-Riemann Hypothesis"** (quasi_rh.tex):
  principal-conclusions guide l. 184–195; Section sec:roots (l. 10430–10470): all roots of X^n − a in F_p in deterministic
  time poly(n, log p) (Adleman–Manders–Miller without ERH); Theorem thm:even-splitting (l. 10605, Gao; Rónyai) and Theorem
  thm:three-divisible-splitting (l. 10735, Rónyai's case r = 3), and Lemma lem:eisenstein-cubic-nonresidue (l. 10745–10815):
  a non-cube a + bω of F_p[ω] (p ≡ 2 mod 3) with a² − ab + b² ≤ C(log 3p)^{32}, C effective, not computed (precursor of our
  Thm E). → r-th roots, Rónyai r = 2, 3 and the r = 3 case of small F_{p²} non-residues KNOWN.
- OpenAI factoring preprint, Lemma ff:unitroot (02-finite-fields.tex l. 106–): q-primary generators of extension
  components by the eigenvector equation v^{Q^{d/q}} = ζv (q odd) and root extraction in k[X]/H — the composite-field
  step we first planned to derive (not claimed).
- **Zenodo 23205034 (Kriebel, 6 Oct 2026)**: consequences of OpenAI's family 142 over explicitly represented finite fields
  (dense factorization, prescribed-degree irreducible construction via Berlekamp, Shoup, Rai); states that the prime-field
  theorem "requires family029 general cyclotomic Hecke zero-freeness". → F_{p^k} extension KNOWN.
- Cureton–Thorner, Hexagon 2610.00200v2 (12 Oct): generalises OpenAI's Theorem 1.2 to ζ_K for every K (31/32, Lean); its
  TeX has no algorithmic application (rg 'algorithm|factor|deterministic|square root': only abelian factorization).

## Other checks (nothing pre-empting the routine additions)
- Zenodo (quoted phrases; zen.py): "least quadratic nonresidue" (23442654, 23442990, 23702717: exponents 8, 16 + ε,
  constants not computed; 10836569 old; 23423890 primitive roots); "non-residue"/"nonresidue" "quasi-Riemann" (25 records,
  incl. 23645699: binomials X^d − ℓ irreducible with ℓ ≪ (log p)^{24} under QRH, no constant); "Tonelli-Shanks" (3,
  unrelated); "square roots modulo" (4, unrelated); "deterministic square root" (0); "deterministic factorization" (18,
  only 23205034 relevant); "family142" (only 23205034); "Deterministic Polynomial Factorization" (1, unrelated);
  "auxiliary prime" (3, unrelated); "Cureton" "factorization" (no relevant hit). Descriptions of 23442654, 23702717,
  23442990: no algorithmic statement.
- arXiv listings 2026-10 and 2026-11 (lst.py 'non.?residue|square root|deterministic|factor(ing|ization) of
  polynomials|finite field|polynomial factori' math.NT cs.CC cs.DS cs.SC math.AG; cs.SC 2026-10 by direct listing, 37
  entries): relevant only 2610.32508 (Spaenlehauer, 23 Oct: deterministic root finding with elliptic curves; unconditional
  online phase given a preprocessing datum; cites Rónyai's GRH algorithm, not the half-planes). 2611.02280 (primality),
  2610.27168, 2611.00647: unrelated.
- Papers using the 31/32 strip: 2610.18491 (FSW, Chebotarev), 2610.37179 (sparse irreducibility over ℤ in PP; abstract read),
  2611.06705 (Baldwin–Sivek, 3-manifold problems in coNP; abstract read; harvest4 read the TeX: no factoring). Pathak
  2405.13316 (least character non-residue under zero-free strips/rectangles, constants not numeric; abstract read).
- Hexagon API, all 716 works (titles and abstracts; hxall.py regex 'square root|non.?residue|deterministic … factor|
  Tonelli|Berlekamp|Rónyai|Evdokimov|finite field … algorithm|primitive root|r-th root|Cureton'): algorithmic hits are about
  integers (2611.00115 square roots modulo composites under a #SAT algorithm, 2610.00600/00711 integer factoring, 2610.00697
  primality), none on polynomials over finite fields; search pages 'deterministic', 'square roots', 'finite field',
  'non-residue', 'nonresidue', 'factoring polynomials'.
- HAL (API, 2025–2026): '"polynomial factorization" deterministic', '"factoring polynomials" deterministic',
  '"quadratic non-residue"', '"square roots" deterministic "finite field"', '"quasi-Riemann"', 'Evdokimov factoring',
  title searches: nothing relevant.
- Faithful restatements of the classical ERH steps: Ivanyos–Karpinski–Saxena 0804.1974 (l. 94–97, 100–105, 111–128,
  465–466, 485), Ivanyos–Karpinski–Rónyai–Saxena 0811.3165 (l. 150–167), Bhargava–Ivanyos–Mittal–Saxena 1702.00558
  (l. 727–812). Huang 1991, Rónyai 1988/1992, Evdokimov 1994 and Buchmann–Shoup 1996 themselves were not read (not on
  arXiv; the AMS link returned a login page).
- Our files: Riemann/Deep/OAINum/REPORT.md (R255, 5 Nov 2026): audit of the factoring preprint, no flaw; it records that
  the unconditional theorem rests on the companion's strip (δ = 10^{−6}); no Cureton–Thorner substitution, no numeric
  auxiliary primes, no exponent. Harvest4/HARVEST.md (candidate 7 and the routine line flagged this; "not checked beyond arXiv"), ExplProg
  NOVELTY (OpenAI 30 Sep Cor 1.2 cited for the non-residue bound, not for the algorithm), COROLLARIES (no factoring entry).

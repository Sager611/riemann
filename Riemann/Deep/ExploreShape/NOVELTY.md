# explore53 NOVELTY (19 Nov 2026)
Question: has anyone used per-gap quadratic inequalities with derivatives of order ≥ 2 (only f vanishing at the
gap ends, telescoping jet terms) in the large-gap (Λ) problem, or beaten Bui–Milinovich's Wirtinger step this way?

## Searches
- `python3 .scratch/axq.py` (arXiv API): 'abs:Wirtinger AND abs:zeta' (1 hit: Saker 0906.5458),
  'abs:Wirtinger AND abs:zeros AND abs:gaps' (+1511.04007, band-limited frames), 'all:"large gaps" AND all:zeta AND
  all:consecutive' (12 hits: 2610.28527, 2304.05181, 1704.05834, 1305.4068, 1112.6038, 1109.3855, 1101.3197,
  1002.2695, 1001.0494, 0906.5458, 0903.4007, 0510530), 'abs:"higher order" AND abs:Wirtinger' (6 hits, none on ζ;
  0412269 Böttcher–Widom: higher-order Wirtinger–Sobolev with all derivatives vanishing at the ends), 'ti:Wirtinger'.
- web_search(approach "arxiv"): four queries (higher order Wirtinger gaps zeta; Bui Milinovich Wirtinger mollifier;
  …): only BM 1410.3635, 0903.4007, 0910.2052, 2208.02359, 2610.31744; one query hit HTTP 429.
- `python3 .scratch/lst.py 'Wirtinger|large gap|gaps between (consecutive )?zeros|zero spacing|consecutive zeros|higher derivative' --zenodo '"large gaps" zeta'`:
  2026-10/11 listings (7 cats): 2611.14486, 2610.14127, 2610.17349, 2610.22927, 2610.26038, 2610.28527, 2610.31744,
  2610.36293, 2611.21544: none uses higher-derivative per-gap inequalities; Zenodo: no large-gap item.
- Repo: `rg -i "higher[- ]order wirtinger|third[- ]order|Opial|Beesack|F'''"` over NOGO, TREE, Tree rounds, Deep reports:
  nothing on this.

## Read (TeX)
- Bui–Hall–Subirà Jorge 2511.14415 §§ "The Wirtinger inequality and large gaps", "Hall's first/second approach",
  "Hall's Wirtinger type inequality" (l.431-520): the only known method uses Wirtinger-type inequalities in f and f'
  only: L² Wirtinger (Hall 2005, BM, Bredberg), Beesack's L^{2k} form, Hall's L⁴ form (3λ0(v)∫|f|⁴ ≤ …∫|f'|⁴ +
  6v…∫|ff'|²). No second or higher derivative.
- Saker 0906.5458 (l.360-420) on Hall's "Generalized Wirtinger inequalities, random matrix theory, and the zeros of
  the Riemann zeta-function" (JNT 97, 2002): ∫H(y'/y)y^{2k} ≥ (2k−1)L∫y^{2k}, H even convex: first derivative,
  RMT-conditional mixed moments Z^{2k−2h}Z'^{2h}. Saker 0906.5458 / 1001.0494: Opial-type, first derivative.
- Bredberg 1101.3197: L² Wirtinger with amplifier.

## Verdict
No prior use found of order ≥ 3 per-gap inequalities (with the boundary curvature f''(0), f''(g) and telescoping
terms) in the large-gap problem. Classical ingredients (sine expansions, Markov–Lukács duality, telescoping /
S-procedure) are standard; the transfer and the numbers are ours. Not yet a discovery: floats only (no Arb
certificate of c_j, j ≥ 2), and the mean values c_j, j ≥ 2, need an audit line.

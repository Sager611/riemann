# exp48 NOVELTY (2026-11-06)

Verdict: KNOWN. The completeness of Elsenhans–Klüners–Nicolae's exponent-2, -4, -8 lists given a zero-free interval
(7/8, 1) for real characters was stated and proved by M. Jaswal on Zenodo before us. The all-discriminant extension
(2605.26020's tables) and the boundary-CM corollary are not stated in any source found, but they follow in one line from
Jaswal's theorem and 2605.26020 l.666–668 (routine).

## Prior statement found

- M. Jaswal, "Imaginary quadratic fields of class number at most 3000" (Zenodo concept 23251443; unrefereed):
  v1 23251444 (9 Oct 2026), v2 23269461 (9 Oct), v3 23402528 and v4 23403019 (15 Oct), v5 23547020 (21 Oct).
  Every version has §5.1 "Class groups of exponent dividing 8": "We show that, under (H), the lists of imaginary quadratic
  fields whose class group has exponent 2, 4 or 8 given by Elsenhans, Klüners and Nicolae are complete"
  (v5, pdftotext l.956–1035: Prop. 5.1 unconditional — exactly 1046 fundamental D with |D| < 2·10^9 and
  E(D) | 8, counts 9, 56, 203, 778, largest 163, 5460, 435435, 430950520; Thm 5.2 under (H) — exactly 56, 203, 778 fields).
  (H): L(σ, χ) ≠ 0 for 7/8 < σ < 1 and every real primitive χ; it is the real case of the 7/8 half-plane, hence implied by
  the Lean-checked half-plane of TREE.md N7. Jaswal's proof: (H) gives the hypothesis of EKN's Lemma Siegel for
  |D| > e^{11.2}, Tatuzawa's lemma excludes the exception of Lemma Siegel_one, and 0.0767 > 0.0412 covers Lemma Chen — our
  route (R2) in DERIVATION.md. Our route (R1)/(R3) through Littlewood78 differs only in avoiding Tatuzawa's localization.
  Jaswal mentions 2605.26020 only as extending EKN "to all negative discriminants ... under the same hypothesis on Siegel
  zeros" (l.1033–1035) and does not state that its tables become complete.
- OpenAI, "The Quasi-Riemann Hypothesis" (5 Oct 2026, github.com/openai/math preprints, paper2.pdf; pdftotext
  l.89–94): cites EKN Thm 2 and states only exponent dividing two (Euler's 65 idoneal numbers).

## Checked, no statement of the result

- arXiv TeX (`rg` for Elsenhans, Klüners, Weinberger, idoneal, 1803.02056, "class group(s) of exponent", "exponent of the
  class group"): 2610.21099 (Kerr), 2610.30476 (Guo), 2610.15483, 2610.31402, 2610.21511, 2610.18491 (Frei–Sofos–Widmer),
  2610.26835 (Liu), 2005.12512: no hit.
- OpenAI QRH 30 Sep (paper.pdf) and "Uniform exclusion of Landau–Siegel zeros" (1 Oct): no hit.
- Zenodo 23607492 (Dr Vengy, 23 Oct): E(D) ≥ (δ + o(1)) log|D|/log log|D| under a half-plane; cites EKN only for Lemma 1;
  no lists.
- arXiv search (web_search, arxiv approach; 6 queries: exponent 4 unconditional; EKN + Siegel; exponent dividing 8
  complete list; quasi-RH + exponents; idoneal + half-plane; Weinberger + Siegel): nothing beyond EKN, 2610.18491 and
  2610.26835 (`logs/novelty_websearch.txt`).
- arXiv monthly listings math.NT and math.GM, 2026-10 and 2026-11, titles matching 'class groups? of (small )?exponent|idoneal|
  imaginary quadratic|Siegel zero|Landau-Siegel|exponent of (the )?class', and Zenodo '"class group" AND exponent'
  (`python3 .scratch/lst.py ...`, output `logs/novelty_lst.txt`): 9 arXiv titles (2610.15483, 2610.21099, 2610.21511, 2610.30476,
  2610.31402 read in TeX above; 2611.06002, 2610.10452, 2610.11210, 2610.34556 off topic) and 7 Zenodo records (23547020 =
  Jaswal v5, 23607492 above, 5 off topic).
- arXiv API, full archive (`.scratch/axq.py`): abs "class group" AND exponent AND "imaginary quadratic" (1809.05750,
  divisibility); exponent AND idoneal (none); "Siegel zero" AND "class group" AND exponent (EKN only).
- HAL API: '"class group" AND exponent AND "imaginary quadratic"' and text search: 0 records; latest "class group" records
  off topic.
- Our files: `rg -i '23547020|Jaswal'` finds Jaswal only as a class-number source (discoveries/class-number-5000.md,
  -10000.md, Harvest2 l.67, Combos/NOVELTY2.md: "Jaswal 23403019 (class numbers ≤ 3000 assuming 7/8)"); no file records
  his §5.1. Harvest3 listed Zenodo 23547020 among the sources to check for candidate 1 but did not read its §5.1.

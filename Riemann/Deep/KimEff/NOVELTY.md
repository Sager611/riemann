# Novelty check (kimeff, 6 Nov 2026)

Question: is an explicit (effective) hyperbolicity region of Kim–Lee shape n >= C d^c, c < 3/2, for J^{d,n} of Xi known?
Answer: none found. The best explicit large-n region before this work is our two-thirds wedge (savvak); Kim–Lee's
d^{1/2+eps} shape is ineffective.

## Sources scanned

- arXiv monthly listings (show=2000), 2026-11 and 2026-10, math.NT, math.CV, math.CA, titles matching
  'Jensen|hyperbolic|derivatives|Laguerre|Polya' (.scratch/lst.py): 33 title hits, none relevant (2610.38870 is an asymptotic
  Jensen–Pólya criterion for GRH at shift 0; 2610.35783 Turán expressions for Stirling-type polynomials; 2610.37931
  generalised hyperbolic polynomials, combinatorial).
- arXiv API (sorted by date): abs:"nonreal zeros" AND derivatives (4 hits: 2609.37164 heat flow, 1702.00626, 0902.0413
  Hawaii conjecture, 0705.1854); "non-real zeros"+derivatives+entire (math/0510502, Bergweiler–Eremenko–Langley, Pólya's
  conjecture, qualitative); "Jensen polynomials" AND effective (2211.05072, 2210.08617, 1910.01227 GORTTW: c e^{d/2}, c not
  computed); "Jensen polynomials" AND Riemann (9 hits incl. 2105.05386 Kim–Lee, 2008.07206 Farmer, 2007.13582 O'Sullivan,
  1902.07321 GORZ, 2605.31356 free probability); "Jensen polynomials" AND Bessel (1405.1638, Turán stable regions);
  "Jensen polynomials" AND explicit (2605.29991, 2306.07459, 1904.12727: other sequences); repeated differentiation of xi;
  derivatives of xi + zeros (0803.0425 pair correlation of xi' zeros etc.). None gives an explicit Jensen region of this shape.
- OpenAlex: Kim–Lee W3162564323 (0 citers). Kim PAMS 1996 W1528064814 (13 citers: Craven–Csordas composition theorems,
  Farmer–Rhoades 2005, Ki–Kim 2000 Fourier–Pólya, de Bruijn's question 2003, Farmer 2022, Laguerre–Pólya multiplicities 1999,
  zeros of derivatives of xi 2019, critical points 2021, Kim–Lee 2021, 2024 entire-function classes): none effective for
  Xi_0. Kim JKMS 1996 W775911138 (0 citers). Farmer–Rhoades "Differentiation evens out zero spacings" (37 citers: random
  polynomials, xi' gaps, Farmer 2022, ...): qualitative. Ki's "The Riemann Xi-function under repeated differentiation"
  (J. Number Theory 120, 2006; known to us from memory: A_m Xi^{(2m)}(C_m z) -> cos z on compacta) is about the other family
  (derivatives in z) and is a limit statement, not an explicit region.
- Zenodo '"Jensen polynomials"' (60 records, newest 2026-11-05): Chung 23987495 and Chen 22951677 (RH claims), Fagiri
  23287787 (J^{2,0} in Lean), Savva 22293642 (two-thirds wedge, audited by us), others unrelated. No explicit Kim–Lee shape.
- HAL '"Jensen polynomials"': Planat hal-05675741 (fitted constant, n >= C0 d^4), hal-03505288 (unrelated).
- GXL evidence API: HTTP 429 (daily quota spent) on both claims; fell back to the searches above.
- Our files: TREE N6, Riemann/Deep/Jensen9/SOTA.md (no explicit N(d) beyond the all-n range except the wedges below),
  Riemann/Deep/SavvaK (n^2 log(n+2) >= 24 d^3), Riemann/Deep/HollandK (n^3 log^2(n+2) >= 7000 d^5), R313/R393/R413 rounds
  (large-n thresholds ineffective; Kim's theorems the blocking step).

## Kim's 1996 papers

Not obtained: the AMS PDF (doi 10.1090/s0002-9939-96-03083-3; Unpaywall lists only the publisher URL) returned the AMS
"Verifying your connection" bot page to curl; the koreascience search page did not download. Our route does not use them:
we prove the needed zero location for Xi_0^{(n)} directly (DERIVATION Sec. 8). So we cannot say which steps of Kim's
general proofs are non-constructive; the general theorem (arbitrary real entire f of order < 2) stays ineffective as far
as we know.

## Verdict

New as far as found: an explicit Kim–Lee-shape region (d <= 6.45 n^2/log^2 n for n >= 9.99e11), proved by a
Bessel-mixture form of Pólya's integral plus Obreschkoff. The ingredients (Pólya's kernel, Obreschkoff, Gauss–Lucas,
Bessel zeros) are classical; the mixture identity Theta^{(n)} = gamma(n) E_n[0F1(; n+1/2; t^2 w/4)] is elementary and may
appear in the literature in some form, but no explicit zero-free region or Jensen region derived from it was found.

## Refresh (11 Nov 2026, after a 4.5-day pause)

- arXiv listings 2026-11 (through 2611.13932) and 2026-10, math.NT/CV/CA, wider title regex (adds 'entire function|real zeros|
  xi-function|Riemann hypothesis'): 23 hits not seen before, none relevant (2611.10726 Weil-form two-interval criterion,
  2611.12438 zeros of zeta and zeta', 2610.04972 Liu–Mao log-concavity, already known; the rest unrelated).
- Zenodo '"Jensen polynomials"' (60 records): new 24101654 (10 Nov; continuous double-Turán inequalities for the xi kernel by
  the Csordas–Dimitrov method, with a barrier beyond the first derivative shift; no Jensen region) and 24104077 (RH exposition).
- Verdict unchanged.

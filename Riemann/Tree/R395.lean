/-!
# Round 395 — node N8 (critical-line proportion), round 9: the quoted inputs of today's discoveries authenticated in the TeX (Conrey's Theorem 2 via Bui-Conrey-Young, Lamzouri's Thm 1.1, Thm 1.4, Lem:Conrey and Pro:BombieriHejhal, Wu's Prop. 3.1 case (B1), BMOR's zero count); two scout 'flaws' were paraphrase mismatches, fixed in the pages; the GXL evidence API finds no pre-emption; Miller test with prime bases up to (log n)^7.8978 (routine, COROLLARIES); Herglotz-Toeplitz positivity, mixing times and auction theory give no route

* Node: N8 critical-line proportion (ninth round; R392 before), used to authenticate the quoted inputs of today's discoveries
  (isolated zeros and linear combinations, N8; numeric Linnik, N10/N11).
* Conrey's Crelle 1989 Theorem 2 as stated by Bui-Conrey-Young (1002.4127, TeX l. 338-346): sigma_0 = 1/2 - R/log T with R > 0
  (left of the line), P(0) = 0, P(1) = 1, theta < 4/7, and c_1 for Q = 1 exactly as quoted in IsoZero DERIVATION 3.1: HOLDS.
* Lamzouri 2610.26841: Thm 1.1 (0.6552 for a non-explicit c, 0.50875 at c = 1/6; l. 211-220), Lem:Conrey (0.0806 + o(1), an
  announcement of Conrey; l. 757-767), Thm 1.4 (31%) and Pro:BombieriHejhal (distinct primitive characters to the same modulus, a
  fixed small c > 0; l. 247-252, 868-870): HOLDS. Our pages said 'any fixed c': corrected to 'a fixed small c' (a smaller c only
  increases the isolated proportion, so nothing changes).
* Wu 1802.09704 Prop. 3.1 (l. 390-413): uniform in complex shifts a, b << 1 of either sign; theta < 4/7 only in case (B1),
  a(n) = mu(n)(F_0 + F_1 (F_2 * F_3))(n) with F_i separable in {F << x^eps, F' << 1/x}, which covers mu(n) P(log(y/n)/log y)
  (DERIVATION 8 says so). The (h,k)^{alpha+beta} of l. 400 is a typo (l. 1353 and the alpha = beta = 0 case have 1 + alpha + beta);
  Wu's (1.6) (+R) and his evaluation a = b = -R disagree in sign, immaterial here because the proposition is uniform in the sign:
  HOLDS (the scout's FLAW compared a shortened paraphrase).
* Bennett-Martin-O'Bryant-Rechnitzer 2005.02989 Thm 1.1 (l. 163-176): conductor q > 1, T >= 5/7, l > 1.567, constants 0.22737,
  2 and -0.5, N(T, chi) counting 0 < beta < 1, |gamma| <= T, as used by explprog: HOLDS.
* Novelty (GXL evidence API): 'more than 34% of the zeros of linear combinations of two Dirichlet L-functions on the line,
  unconditionally': 0 supporting, 1 refuting passage (2610.26841's 31%); 'an explicit numerical polynomial bound for the least prime
  in progressions, valid for every modulus': 1 supporting, 9 refuting (1802.00085's exponential bounds, 1510.07711): no pre-emption.
* Corollary (routine, COROLLARIES): every odd composite n >= 9 has a prime strong witness p <= (log n)^{7.8978}, from explprog's
  non-residue bound and the proper-subgroup lemma for strong liars (under GRH: 2 log^2 n, Bach 1990).
* Exploration: Herglotz and Caratheodory-Toeplitz positivity of the Li coefficients (2 x 2 minors force lambda_n <= lambda_0,
  impossible since lambda_n ~ (n/2) log n; OPUC already R333), mixing times and cutoff for walks on (Z/q)^* (small-prime
  generation gives hitting, not spectral gaps; 1909.09053 uses GRH averages) and auction theory (Prekopa and order-statistic
  closure preserve density log-concavity, not coefficient total positivity or Jensen hyperbolicity): no route.
-/

namespace Fayed.R395
end Fayed.R395

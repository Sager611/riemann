/-!
# Round 414 — node N11 (GRH families, Siegel zeros, class numbers), round 14: the inputs of today's unit-generated-orders discovery authenticated in the TeX (Kopp-Lagarias Prop. 2.3 with its n ranges, the Delta < 10^10 table captions, Hua's ineffective finiteness, the genus count 2^{mu-1} with mu - 1 <= omega, the class number formula for orders with ratio 1 at Delta = 68640); prior art: Biro and Byeon-Kim-Lee settle class number one only, and Louboutin-Mollin-Williams's exponent-two list for maximal orders leaves one possible exception (complete under GRH); no pre-emption; radio interferometry, tensegrity and pharmacokinetics give no route

* Node: N11 GRH families, Siegel zeros and class numbers (fourteenth round; R411 before). Round type: authentication of the
  quoted inputs of a fresh discovery (kopplag, discoveries/unit-generated-orders.md), prior art, pre-emption.
* [check] Kopp-Lagarias 2512.11311 v3: Prop. 2.3 (TeX l.543-553) gives the unit-generated real quadratic orders as
  Delta = n^2 - 4 (n >= 3) and n^2 + 4 (n >= 1); Tables 2A/2B are complete for Delta < 10^10 (captions l.1668, 1704); Thm 1.5 is
  ineffective via Hua (l.307); no unconditional list is claimed: HOLDS.
* [check] Genus theory (KL Prop. 3.4 = Halter-Koch Thm 6.5.2, l.1525-1553): |Cl^+[2]| = 2^{mu-1} and mu - 1 <= omega(Delta),
  stated at l.1551-1553; recomputed from the definition for all 631 Delta = n^2 -+ 4 < 10^5: equality only for 32 | Delta: HOLDS.
* [check] Class number formula for orders (KL eq. cnorder l.828, eq. BS00 l.887-889) agrees with our (F2); at Delta = 68640
  (f = 2, Delta_0 = 17160, h(O) = 16) both sides are 89.0932789684... (PARI, ratio 1.000000000000000; field h = 8, R = 5.5683...): HOLDS.
* Prior art: Biro 2003 (Yokoi n^2 + 4, Chowla 4n^2 + 1) and Byeon-Kim-Lee 1508.05644 (n^2 - 4) settle class number one only (KL
  Thms 4.1, 4.4, 4.6); Louboutin-Mollin-Williams 1993 give the exponent-two list for maximal orders of extended Richaud-Degert type
  'with one possible exceptional value ruled out by the GRH' (the GXL passage in NOVELTY.md; KL l.1711: complete under GRH; the
  scout's attribution of the exception to the class-number-one paper is contradicted by that passage): no unconditional
  exponent-two classification, maximal or not: NOT KNOWN.
* [check] Pre-emption (arXiv 'unit-generated', au:Kopp AND au:Lagarias, OpenAlex citations 0, Zenodo, HAL 0): no settlement and
  no part II: HOLDS.
* Exploration: radio interferometry (positive-definite visibilities, closure invariants), tensegrity (PSD stress matrices,
  Gortler-Healy-Thurston rigidity) and pharmacokinetics (reversible compartments realise Jacobi matrices with real spectra) give
  finite positivity or real-spectrum certificates only, with no encoding of the Euler product and no passage from finite
  windows to all test functions: no route.
-/

namespace Fayed.R414
end Fayed.R414

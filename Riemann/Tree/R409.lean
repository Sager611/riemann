/-!
# Round 409 — node N11 (GRH families), round 12: the quoted inputs of today's two N11 results authenticated in the TeX: Kogler covers only 0/1 polynomials (general laws announced as future work), Bary-Soroker-Hokken-Kozma-Poonen prove only limsup = 1 for +-1 polynomials, Eberhard's and Ferber-Jain-Sah-Sawhney's ERH enters only through Breuillard-Varju-type prime sums (two scout FLAWs rejected on re-derivation), Cohen-Trudgian and McGown-Trevino-Trudgian quoted exactly; superconductivity, viscoelasticity and dynamo theory give no route

* Node: N11 GRH families and Siegel zeros (twelfth round; R408 before). Round type: authenticate the quoted inputs of today's
  N11 results (least prime primitive root, discoveries/prime-primitive-root.md; random-polynomial irreducibility without RH, bvirr).
* Kogler 2609.26771 thm:main (intro l. 4, 16-24): 0/1 coefficients only, P(irreducible) -> 1 with reducibility sqrt(2/(pi n)) + O(1/n),
  unconditional; l. 60: the method 'generalizes to arbitrary coefficient measures', future work; no +-1 or Galois statement: HOLDS.
* Bary-Soroker-Hokken-Kozma-Poonen 2308.04878 Thm thm:short (l. 136): limsup P(+-1 polynomial irreducible) = 1; the full limit is
  conjectured (l. 126) and 'known to hold under GRH' (l. 113, via Breuillard-Varju); no 2024-2026 proof of the limit found: HOLDS.
* Eberhard 2008.01223: ERH only in BV's Prop 9 (l. 604-608) and the local-global Prop (l. 695-700), used at X = cn/2 (part 1) and
  X = c n^{1/2}/m^{5/2} (part 2). The scout's FLAW read the premise (Cureton-Thorner's strip) as unproved; re-derived: with error
  n^5 log(Hn) e^{-X/32} against e^{X - cn}, X = 32cn/33 balances both (e^{-cn/33}), and part 2's m n^{2m+1} X^2 e^{-X/32} -> 0: HOLDS.
* Ferber-Jain-Sah-Sawhney 2106.04049 (l. 481-487, 514-521): the scout's FLAW claimed log Delta_K = O(n^2 log n) overwhelms e^{-X/32}
  at X = c n^{1/4}; re-derived: the error n^3 X^2 e^{-X/32} log n still tends to 0 (exponential in n^{1/4} beats n^3), and l. 514's
  PNT for Q with Riemann error is replaced by the 7/8 error: HOLDS (constants only).
* Cohen-Trudgian 1602.02440 (l. 40): 'It is not known whether g^(p) < p for all p, or even for all sufficiently large p', g^ the
  least prime primitive root; McGown-Trevino-Trudgian 1508.05182 Thms 1-2 (l. 45, 72-78): the GRH bounds exactly as quoted: HOLDS.
* Exploration: superconductivity (the BCS/Eliashberg gap equation is a top-eigenvalue problem for a positive compact kernel),
  linear viscoelasticity (completely monotone relaxation moduli and Kramers-Kronig give Stieltjes positivity) and dynamo theory
  (kinematic induction operators are dissipative non-self-adjoint semigroups): no arithmetic kernel, no route.
-/

namespace Fayed.R409
end Fayed.R409

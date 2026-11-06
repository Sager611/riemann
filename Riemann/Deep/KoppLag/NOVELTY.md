# NOVELTY: Kopp–Lagarias Problem 1.7 (unit-generated real quadratic orders with 2-torsion class group)

Scan by agent kopplag, 2026-11-06, 12:25–13:05 UTC. Raw outputs (Zenodo, GXL) in `logs/novelty_scan.txt`; hits summarized here.

## Statement checked
The complete list of unit-generated real quadratic orders O = ℤ[ε] (Δ = n² ∓ 4) with Cl(O) = Cl(O)[2], and with one class per genus
(Cl⁺(O) = Cl⁺(O)[2]); proof that the Kopp–Lagarias tables (Δ < 10¹⁰) are complete.

## Source problem
- arXiv 2512.11311 v3 (Kopp–Lagarias, 21 Apr 2026), TeX read: Thm 1.5 (l.303–306) finiteness "ineffective; it depends on the ineffective result of Hua" (l.307–308);
  Problem 1.7 (l.329–336); tables "complete for Δ < 10¹⁰" and "conjecturally complete" (abstract l.87–90, l.145, l.342, l.428, Table captions l.1667, l.1703);
  l.1743: "The complete classification of unit-generated orders with 2-torsion class group seems out of reach at this time, in parallel with
  the idoneal number problem being unsolved." Their code (github.com/gskopp/UnitGeneratedOrders, functions.sage + notebooks): Sage class numbers
  for n ≤ 10⁵, compared with the 2-rank.
- Count erratum found: text l.1609 says 86 orders; Tables 2A/2B list 55 + 37 = 92 entries = 89 distinct real orders (+ Δ = −4, −3).

## Searches (all negative for a solution)
- arXiv monthly listings math.NT 2026-10 and 2026-11 (`.scratch/lst.py`, titles matching unit.generated|Kopp|Lagarias|real quadratic orders|class group|
  one class per genus|idoneal|Richaud|Degert): 8 title hits, none on unit-generated orders or 2-torsion class groups of real orders.
  Nearest: 2610.29347 (21 Oct 2026; d − x² factorizations for Q(√p), Q(√2p): different problem), 2610.18491 (FSW, quasi-RH applications: ℓ-torsion, small generators; no real-quadratic exponent-2 list).
- arXiv API: au:Kopp AND au:Lagarias (3 papers: 2512.11311, 2407.08048, 2212.09177; no sequel "II" posted); abs:"unit-generated" (only 2512.11311 in math.NT);
  abs:"one class per genus" (math/0401334 only); abs:"Richaud-Degert" (2512.11311, 1812.02488, 1501.06051, 1209.4958); abs:"exponent two" AND "real quadratic" (0);
  abs:"class group" AND exponent AND "real quadratic" (2607.11926, 2211.12279: unrelated); abs:idoneal (10 hits, none on real quadratic orders).
- Semantic Scholar citations of 2512.11311: 2 (2407.08048, its predecessor; 2610.16530, bounded continued fractions in Q(√127): unrelated).
- Zenodo (quoted phrases): "unit-generated orders" (23333043, 21986780), "unit-generated" (+21484481), "one class per genus" (22743167), "Kopp-Lagarias",
  "Kopp and Lagarias", "2512.11311", "unit-generated order", "unit generated order" (0 each). 23333043 and 21484481 are 1.1 kB AI summaries
  (E8 Intelligence Research) with no theorem; 22743167 is about imaginary one-class-per-genus discriminants (Kronecker limit formulas).
- HAL API: "unit-generated orders", "one class per genus" "real quadratic", "Kopp" "Lagarias" "unit-generated": 0 hits.
- GXL Evidence API (2 claims): the only supporting passages are KL's own tables ("complete for Δ < 10¹⁰", "conjecturally complete");
  Louboutin–Mollin–Williams 1993 ("Class groups of exponent two in real quadratic fields") is returned as neutral:
  "determining all real quadratic fields of ERD-type having class group of exponent ≤ 2 (with one possible exceptional value ruled out by the GRH)".
- Our files: `rg -i 'idoneal|unit-generated|Kopp|one class per genus'` over the repo: only Harvest4 (this candidate), ImqClass2/novelty.txt (title seen),
  Exp48 (imaginary exponents 2, 4, 8: Jaswal Zenodo 23251444/23547020), Combos/ArtinHalfPlane notes (OpenAI QRH 5 Oct: idoneal numbers, imaginary only).

## Prior results and what is new
| result | scope | status |
|---|---|---|
| Kopp–Lagarias 2026 | unit-generated orders, Δ < 10¹⁰ | computed list; completeness open (Problem 1.7) |
| Louboutin–Mollin–Williams 1993 | maximal orders of ERD type | complete with one possible exception; none under GRH |
| Biró 2003, Byeon–Kim–Lee | class number 1, maximal n² ± 4 | unconditional (h = 1 only) |
| OpenAI QRH note (5 Oct 2026), Jaswal (Zenodo) | imaginary: idoneal numbers, exponents 2/4/8 | complete given the half-plane |
| **this work** | all unit-generated real orders, wide and narrow | complete given the half-plane; elementary and unconditional for Δ < 10¹⁸ |

New: (1) the solution of Problem 1.7 (both parts), (2) the elementary split-prime criterion (Lemma A, Cor. B) for unit-generated orders, which alone
extends KL's verified range from 10¹⁰ to 10¹⁸ without GRH or the half-plane, (3) removal of LMW's possible exception for the n² ± 4 maximal orders,
(4) the count erratum. Lemma A is a self-contained version of the classical principle (Lagrange) that integers |m| < √Δ/2 represented by the principal form
lie on its (here very short) cycle, a standard tool in Richaud–Degert class-number-one criteria (Yokoi, Mollin–Williams, Louboutin); its use with genus
theory for non-maximal orders and an explicit L(1, χ) bound is not in any source found.

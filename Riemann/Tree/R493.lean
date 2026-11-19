/-!
# Round 493 — node N11 (GRH families and Siegel zeros), round 19, harvest of the GRH-conditional arXiv papers of 15 Sep – 12 Nov 2026 not yet in our files (16 ids, 15 read): four use GRH only through a power-saving Chebotarev or character-sum bound, so the Lean-checked half-planes would remove it (2610.25793 generalized Euclidean pure orders of odd degree, 2609.30573 Koblitz–Zywina primes with κ < 1/128, 2610.23660 signs of Alladi–Johnson constants, 2611.10489 divisor moments: candidates for one write-up agent); 2610.01837 already applies the quasi-RH itself; seven need the full critical line

* Node: N11 (19th round; R492 before). Round type: harvest of GRH-conditional statements made unconditional by our strips (lesson 96),
  continuing R492 backwards: arXiv API (abs: GRH / generalized / generalised / extended Riemann hypothesis, submitted 15 Sep – 12 Nov
  2026) gave 16 ids not in our files; 15 TeX sources pre-extracted in .scratch/r493 (2609.37434, a random-model paper, skipped).
  Load ≈ 9: all eight scouts answered.
* HALF-PLANE SUFFICES (candidates; one write-up agent with a reviewer per item, COROLLARIES if confirmed):
  - 2610.25793 ('Orders in the Euclidean Court', 19 Oct): Theorem maingeneralGRH (l.93: for fixed d ≥ 2 finitely many generalized
    Euclidean pure orders Z[m^{1/d}]; GRH needed only for odd d) uses GRH only through Prop. chebotarev (l.292–297, Serre (20_R):
    |π_C(x) − |C|/|G| Li(x)| ≤ A|C| x^{1/2} log([E:Q]Px)) in Lemma generalGRH (l.199–202, 299–313: two primes q₁ < q₂ ≤ |m|^{1/6},
    q ∤ m, with m not a p-th power mod q, from π_C(|m|^{1/6}) ≥ 2 + ω(p|m|)). Orchestrator read l.93, 199–202, 291–314: with
    Cureton–Thorner's zero-free Re s > 31/32 for ζ_E, E = Q(m^{1/p}, ζ_p), the count has error O_d(x^{31/32} log|mx|), enough at
    x = |m|^{1/6} (or (log|m|)^{33}); no other GRH use (rg). The scout's x = (log|m|)³ arithmetic dropped the log|pmx| factor.
  - 2609.30573 (refined Koblitz conjecture, congruence obstructions): Theorem Quanticep (l.264–270: for every 0 < κ < 1/8 infinitely many
    good primes p with the stated size condition on #E_p(F_p), under GRH for Dedekind zeta functions via Theorem ECDTA, l.1378): 31/32
    gives κ < 1/128 (scout; not yet checked).
  - 2610.23660 (Alladi–Johnson series, μ(n)ω(n)/n over P⁻(n) ≡ l mod k): the GRH count of classes with C(k, l) < 0 (φ(k)(1 − o(1)),
    l.381–510) uses GRH through a character prime-sum bound; any θ < 1 gives the same conclusion (scout; not yet checked).
  - 2611.10489 (local vanishing density and divisor moments): Chebotarev/Sato–Tate error terms used only as a power saving (l.1113;
    scout: same result on a smaller range; not yet checked).
* KNOWN: 2610.01837 (cubic field counts, error X^{2/3−θ}) states the quasi-RH application itself (Theorem quasiapplthm, l.191–220,
  with τ = 1/32 from the Lean-checked Dedekind half-plane): nothing for us.
* NEEDS THE FULL CRITICAL LINE: 2609.16807 and 2609.23390 (central values of quadratic Hecke L), 2610.24958 (L^{(ℓ)}(1/2)), 2609.18569
  (Euler–Kronecker extremes: square-root cubic character sums), 2610.38884 (Serre's conjectures over Q(√13): GRH Odlyzko discriminant
  bounds), 2610.23201 (S(t, π) on GL₃). NO GRH RESULT: 2610.40267, 2609.17390, 2610.10177, 2609.28315 (GRH only in cited prior work).
* Exploration (tool transfers, tested): Coulomb/log gases with Serfaty's renormalized energy → N1 (a quadratic log-gas equilibrium does
  not constrain ζ zeros: the first gap violates its Euler–Lagrange equation), actuarial ruin theory → N1 (Lundberg's inequality needs a
  positive claim law; the Weil prime kernel is signed), optimal recovery (Micchelli–Rivlin) → N13 (no convex class of candidate
  functions or sampling operator for ζ): no transfer.
-/

namespace Fayed.R493
end Fayed.R493

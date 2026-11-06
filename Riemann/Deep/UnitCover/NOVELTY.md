# Novelty check (2026-11-06, 06:59–07:25 UTC)

Claim checked: for every Q ≥ 3, an explicit-size set of O(log Q) primes of size ≤ (c' log Q)^{7.8978} covering (ℤ/qℤ)^× by subset
products for EVERY q ≤ Q (no exceptional moduli), with numeric constants; and a leading constant below 2/log 2 for an all-moduli
statement.

## Sources read
- **Klurman–Shparlinski–Teräväinen, arXiv 2609.25393** (v1 only, 21 Sep 2026; TeX read): Thm 1.1 (TeX l.109–122): (log Q)^{1+ε} primes
  ≤ (log Q)^{C_*(A+1)}, all but O(Q(log Q)^{−A}) moduli, C_* not computed; GRH sketch l.101–105 uses all primes ≤ (log q)^B (≫ log Q
  primes); l.129: the exceptional set "is unavoidable with present technology" without Vinogradov's least-non-residue conjecture.
- **Hughes, arXiv 2610.04678** (v1 only, 3 Oct 2026; TeX read): Thm main (l.71–79): (1/log 2 + o(1)) log Q primes, all but O_A(Q ℓ^{−A})
  moduli, three exceptional sets, the third (random failures, l.382–395) independent of zeros; Thm structured (l.87–95):
  2 log Q/log 2 + O_A(log Q/log log Q) primes for every q ∉ 𝓔_an coprime to ∏𝒫, Q sufficiently large; Remark l.119–124:
  "We do not track effectivity". With the half-plane every character is good ((ZF), l.262, σ < 1 − B), 𝓔_an = ∅, so his structured
  theorem becomes all-moduli (coprime form) with 2/log 2 + o(1) and no numbers. Not stated by Hughes.
- **Frei–Sofos–Widmer, arXiv 2610.18491 v2** (14 Oct; TeX searched with rg for subset product, cover, Regev, Boolean, Klurman,
  Hughes, Shparlinski): no covering statement; its "small generators" are of number fields.
- **Zenodo 23581412** (Tolbert, "Deterministic Prime Bases for Regev's Algorithm", 22 Oct, v3): relation lattices for almost all
  inputs; explicitly does not assume that small primes generate the unit group; different statement (cite only).
- Zenodo 23599288 (Na, products of distinct integers, Ang's conjecture): different problem.

## Searches
- arXiv listings 2026-10 and 2026-11, math.NT and math.CO (.scratch/lst.py, regexes 'unit group|subset product|small primes|cover|
  Boolean|Regev|generat' and 'unit group|subset product|small primes|cover|Boolean|Regev|generat|multiplicative group|quasi-Riemann|
  zero-free half|half-plane'): relevant hits only 2610.04678 (Hughes), 2610.18491 (FSW), and half-plane papers 2610.12234,
  2610.15291, 2610.26835, 2611.00924 (abstracts read: zero-free regions only, no covering application); 2610.25228 (covering
  𝔽_q by translates of multiplicative subgroups: different problem).
- Zenodo API, quoted phrases: "unit group" (115 hits, none on topic), "subset products" (3: off topic), +"subset product" +primes (0),
  +"unit group" +"small primes" (1: Tolbert 23581412), +"subset products" +"unit group" (0), +Regev +"small primes" (1: lattice
  enumeration 23018079, off topic), +"quasi-Riemann" +"subset" (1: off topic), +"quasi-Riemann" +"unit group" (0),
  +"half-plane" +"subset products" (0), +"Erdos-Renyi" +"unit group" (0), "cover the unit group" (0), "prime alphabet" (off topic).
- HAL API: "subset products" AND primes (0), "unit group" AND "small primes" (0), Regev AND "small primes" (0), "quasi-Riemann" (2, old).
- GXL evidence API, 2 claims: (1) "For every Q there is a set of O(log Q) small primes whose subset products cover the unit group
  (ℤ/qℤ)^* for every modulus q up to Q, with explicit constants": 254 passages, supports only Hughes's optimality line, refutes =
  KST l.129 and Hughes's 𝓔_an (exceptional sets); (2) "A zero-free half-plane for all Dirichlet L-functions implies that subset
  products of few small primes represent every unit modulo every q": 244 passages, nothing on covering (products of three primes
  0708.1562; 1705.06087; 2609.37434 Granville model).
- Our repo (rg over all files except this folder and .scratch): only Harvest2 candidate 1 (this brief) and Harvest1's list entry
  of 2610.04678; no NOVELTY file mentions unit-group covering.

## Verdict
No source states a covering of every q ≤ Q by O(log Q) small primes with numbers, nor an all-moduli leading constant below
2/log 2. Known ingredients: Erdős–Rényi doubling and Hughes's first-moment covering lemma (lem:cover, TeX l.132–177; k = 1 here), the explicit-formula machinery of
explprog. New here: the empty exceptional set made explicit (Lemma D with numbers), the k-th moment refinement of the collision
step with the enumeration of closed relation sets (k ≤ 4), all constants, the small-Q exact tables. Risk: KST or Hughes v2 could
add a half-plane remark (both v1 on 6 Nov 2026).

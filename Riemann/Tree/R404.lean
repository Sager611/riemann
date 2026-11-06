/-!
# Round 404 — node N11 (GRH families and Siegel zeros), round 9: the quoted inputs of today's covering discovery authenticated in the TeX (Klurman-Shparlinski-Teravainen Thm 1.1, l. 101-105 and l. 129; Hughes Thms 1.1-1.2, l. 121, the (ZF) condition and E_an, l. 262-325; Frei-Sofos-Widmer prove no covering); its Proposition 3 rechecked independently; KST's effectivity remark (l. 252) is superseded; plasma physics, rogue waves and jamming give no route

* Node: N11 GRH families and Siegel zeros (ninth round; R399 before). Round type: authenticate the inputs of a fresh discovery
  (discoveries/unit-group-cover.md, unitcover).
* Input authenticated: Klurman-Shparlinski-Teravainen 2609.25393 (Generation.tex): Thm 1.1 (l. 109-122): (log Q)^{1+eps} primes
  <= (log Q)^{C_*(A+1)}, an absolute C_* >= 1 not computed, covering all but O_{eps,A}(Q (log Q)^{-A}) moduli; l. 129 'unavoidable
  with present technology'; l. 101-105 the GRH sketch (B > 4, P = (log q)^B): HOLDS.
* Input authenticated: Hughes 2610.04678: Thms 1.1-1.2 (l. 71-96), the effectivity remark (l. 121) and the zero-free condition (ZF)
  defining E_an (l. 262-325) as quoted; 1 - B = 0.1266 suffices; Hughes states no all-moduli consequence with numbers: HOLDS.
* Pre-emption check: Frei-Sofos-Widmer 2610.18491v2 prove small splitting primes in Chebotarev form, class-group bounds and small
  generators of number fields (constants effective but unspecified); no covering of (Z/qZ)^x by subset products and nothing
  simultaneous in q: HOLDS.
* Computation: Proposition 3 rechecked independently (numpy subset-product closure): the first 11, 13, 15, 17 primes cover every
  q <= 100, 300, 1000, 2000, and the 16 least primes above 2000 cover every q <= 2000; 17 is the least m for q <= 2000 (4.2 s): HOLDS.
* Harvest: KST l. 252 (Thm 1.1 could be made effective with eta = 1/2, keeping the exceptional-zero term) is superseded by the
  all-moduli numeric statement; Hughes l. 116-117 (is 2/log 2 optimal) and l. 120-123 (effectivity) are handled in UnitCover's
  REPORT: no new candidate.
* Exploration: plasma physics (the Penrose criterion is a Nyquist, i.e. argument-principle, test whose sufficient condition is a
  Hermite-Biehler-type sign condition coming from a positive equilibrium; for xi that sign condition is RH-equivalent, de Branges),
  rogue waves and NLS (the defocusing Zakharov-Shabat operator is self-adjoint, but realizing xi's zeros by such a Dirac or canonical
  system is the de Branges-Lagarias inverse problem, solvable only if RH holds) and jamming (Maxwell-Calladine's R^T R is positive by
  construction; an arithmetic R is the Hilbert-Polya problem itself): no route.
-/

namespace Fayed.R404
end Fayed.R404

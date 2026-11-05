/-!
# Round 362 — node N4 (Robin/Lagarias), round 6: records stand; no GA2 number in [5041, 10^(10^46.79)] and no extraordinary number in (4, 10^(10^46.79)] (routine); the odd-Robin gap of 2008.04787 harvested

* Node: N4 Robin/Lagarias (sixth round; R353 before).
* [check] state of the art: the Oct-Nov listings and Zenodo give 2610.04579, 2610.13864 (unrelated) and no new proof claim: HOLDS.
* Computation (routine, COROLLARIES): Caveney-Nicolas-Sondow (1112.6010, eq. GA2>ega, TeX l. 325-327) prove that every GA2 number
  N has G(N) = sigma(N)/(N log log N) >= e^gamma; with our Robin range 5041 <= n <= 10^(10^46.79) no GA2 number lies in
  [5041, 10^(10^46.79)], and as no GA1 number lies in [5, 5040], no extraordinary number lies in (4, 10^(10^46.79)] (RH iff 4 is
  the only one); before 10^8576 (their Cor. lowerbound).
* Adversarial: the newest Zenodo Robin-based RH claim found (F. Vega, 7018145, ICRDM 2022 slides): its Lemma 2 restates Robin's
  Proposition correctly, so a scout's FLAW there is not one; the first false step was not located this round.
* Harvest: 2008.04787 asks (l. 85) whether f_{2,1}(n) exceeds its limsup infinitely often and (l. 400) whether the odd Robin
  inequality holds on the finite gap 277093800961598968947975375 < n < c_0 (candidate: close it with RH to 3e12 and explicit
  theta bounds); 1905.09327 (l. 178) asks to improve Ramanujan's constant c_1 > 1.393.
* Wild (Robin-type criteria for other L-functions): no sigma-over-ideals or progression criterion equivalent to GRH found
  (1010.3239: the Dedekind psi function, for zeta only; 1012.0671: t-free integers): NO IDEA.
* Exploration: GCD graphs and the Duffin-Schaeffer compression (1907.04593, 2404.15123 control overlaps, not GCD sums; 2201.12968
  is already sharp for log-GCD sums), quantum error correction (shadow and quantum MacWilliams constraints stay finite LPs;
  Duursma-type zeta analogues are classical, 1709.03396), integral geometry (Helgason-Ludwig moments constrain Radon data, not
  the kernel of Xi; R345): no new route.
-/

namespace Fayed.R362
end Fayed.R362

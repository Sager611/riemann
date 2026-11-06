\\ body of pari_checks.gp (appended by checks.py after the data matrix R)
\\ per order: quadclassunit(D).no = exact h; wide 2-torsion iff all cyc entries = 2 (PARI: GRH-conditional, cross-check only);
\\ quadunit(D) = (n + sqrt D)/2 (trace n, norm N) for D != 5; class number formula h(O) log eps = (sqrt D/2) L(1,chi_D0) prod_{p|f} (1 - (D0/p)/p);
\\ L(1, chi_D0) >= 0.11/max(3, log log D0) (Littlewood78) as a sanity check.
bad = 0; worstQ = 0.; minLratio = 10^9; cnt = 0;
{
for (i = 1, #R~,
  my(n = R[i,1], N = R[i,2], D = R[i,3], f = R[i,4], D0 = R[i,5], h = R[i,6], w2 = R[i,7]);
  my(q = quadclassunit(D));
  if (q.no != h, bad++; print("quadclassunit h mismatch ", [n,N,D,q.no,h]));
  my(isw2 = (#q.cyc == 0) || (vecmax(q.cyc) == 2));
  if (isw2 != w2, bad++; print("2-torsion mismatch ", [n,N,D,q.cyc,w2]));
  my(u = quadunit(D));
  if (D != 5 && (trace(u) != n || norm(u) != N), bad++; print("unit mismatch ", [n,N,D,u]));
  my(eps = if (D == 5, (1 + sqrt(5))/2, (n + sqrt(D))/2));
  my(L = lfun(lfuncreate(D0), 1));
  my(F = factor(f), E = 1.);
  for (j = 1, #F~, E *= 1 - kronecker(D0, F[j,1])/F[j,1]);
  my(Q = h * log(eps) / (sqrt(D)/2 * L * E));
  if (abs(Q - 1) > 1e-25, bad++; print("class number formula mismatch ", [n,N,D,Q]));
  worstQ = max(worstQ, abs(Q - 1));
  minLratio = min(minLratio, L / (0.11 / max(3, log(log(D0)))));
  cnt++;
);
}
print("PARI checks: ", cnt, " orders; mismatches ", bad, "; max |Q-1| = ", worstQ, "; min L(1,chi_D0) / (0.11/max(3,loglog D0)) = ", minLratio);

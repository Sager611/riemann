\\ exp48 spot checks (PARI/GP): class group structure of the largest listed discriminants.
\\ quadclassunit is GRH-conditional in general; for the fundamental ones bnfcertify makes the result unconditional.
\\ The scan (scan.c) obtained the same h and 2-/4-torsion counts by counting reduced forms (unconditional).
L = [-5460, -435435, -11148180, -12517428, -15337315, -15898740, -17168515, -28663635, -29493555, -31078723, -430950520, -887040, -494689888, -1723802080];
{
for (i = 1, #L,
  D = L[i];
  q = quadclassunit(D);
  cyc = q.cyc; h = q.no; ex = if(#cyc, cyc[1], 1);
  fd = isfundamental(D);
  om = omega(D);
  cert = "-";
  if (fd, bnf = bnfinit(x^2 - D, 1); cert = bnfcertify(bnf); if (bnf.cyc != cyc, error("bnf/quadclassunit mismatch ", D)));
  printf("D=%d fundamental=%d omega=%d h=%d exponent=%d cyc=%s bnfcertify=%s factor=%s\n", D, fd, om, h, ex, cyc, cert, factor(D));
);
}
quit;

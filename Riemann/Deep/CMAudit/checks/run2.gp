read("Riemann/Deep/CMAudit/checks/drv.gp");
setrand(2); OUT = "Riemann/Deep/CMAudit/checks/out_q2.txt"; fac = [1/8, 1/4, 1/2, 1, 2, 4, 8, 16];
write(OUT, "BOX Q(sqrt2) [trials, max count/bound] = ", boxtest(y^2 - 2, 300));
write(OUT, "BOX Q(sqrt10) [trials, max count/bound] = ", boxtest(y^2 - 10, 200));
{ my(K = mkK(y^2 - 10, findalpha(y^2 - 10, 300, "p"))); write(OUT, "EIG Q(sqrt10) [h, #ideals, top eigenvalue, h*max_C sum w^-1, rel diff] = ", eigcheck(K, 200)); }
foreach ([300, 3000, 3*10^4, 3*10^5, 1500000], T, runfield(OUT, y^2 - 2, T, "p", fac));
runfield(OUT, y^2 - 2, 3*10^4, "c", fac);
foreach ([1000, 10^4, 6*10^4], T, runfield(OUT, y^2 - 10, T, "p", fac));
write(OUT, "DONE");

read("Riemann/Deep/CMAudit/checks/drv.gp");
setrand(5); OUT = "Riemann/Deep/CMAudit/checks/out_q5.txt"; fac = [1/8, 1/4, 1/2, 1, 2, 4, 8, 16];
write(OUT, "BOX Q(sqrt5) [trials, max count/bound] = ", boxtest(y^2 - y - 1, 300));
{ my(K = mkK(y^2 - y - 1, findalpha(y^2 - y - 1, 2000, "p"))); write(OUT, "EIG [h, #ideals, top eigenvalue, h*max_C sum w^-1, rel diff] = ", eigcheck(K, 150)); }
foreach ([300, 1000, 3000, 10^4, 3*10^4, 10^5, 3*10^5, 10^6, 4*10^6], T, runfield(OUT, y^2 - y - 1, T, "p", fac));
foreach ([10^4, 10^5], T, runfield(OUT, y^2 - y - 1, T, "c", fac));
foreach ([10^4, 10^5], T, runfield(OUT, y^2 - y - 1, T, "e", fac));
write(OUT, "DONE");

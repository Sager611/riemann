import cypari2; gp=cypari2.Pari(); gp.allocatemem(200*10**6); gp.set_real_precision(30)
gp("L3=lfundiv(lfuncreate(x^4+8*x+12),lfuncreate(1)); Lr=lfundiv(lfuncreate(x^3-x-1),lfuncreate(1));")
print("P6b lfuncheckfeq (bits of agreement, very negative = OK): L(rho3) A4:", gp("lfuncheckfeq(L3)"), " L(rho) S3:", gp("lfuncheckfeq(Lr)"))
print("P6b |L(rho3, 1/2+i*1.81487947273257758)| =", gp("abs(lfun(L3,1/2+I*1.81487947273257758490668747514))"), " |L'| =", gp("abs(lfun(L3,1/2+I*1.81487947273257758490668747514,1))"))
print("P6b conductor of L(rho3):", gp("lfunparams(L3)[1]"), " disc F4:", gp("nfdisc(x^4+8*x+12)"), " Galois group:", gp("polgalois(x^4+8*x+12)[4]"))
K=gp("polredabs(nfsplitting(x^4+8*x+12))"); print("P6b K12 =", K, " degree", gp("poldegree(%s)"%K), " disc K12 =", gp("factor(nfdisc(%s))"%K))


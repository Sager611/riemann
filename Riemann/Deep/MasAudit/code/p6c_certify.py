import cypari2; gp=cypari2.Pari(); gp.allocatemem(200*10**6); gp.set_real_precision(40)
gp("L3=lfundiv(lfuncreate(x^4+8*x+12),lfuncreate(1)); Lr=lfundiv(lfuncreate(x^3-x-1),lfuncreate(1));")
print("P6c root number L(rho3):", gp("lfunrootres(L3)[3]"), " L(rho2,S3):", gp("lfunrootres(Lr)[3]"))
print("P6c Hardy Z of L(rho3) at t=1.80, 1.83:", gp("[lfunhardy(L3,1.80), lfunhardy(L3,1.83)]"), "-> sign change = zero on Re s = 1/2 in (1.80,1.83)")
print("P6c Hardy Z of L(rho2) at t=5.10, 5.13:", gp("[lfunhardy(Lr,5.10), lfunhardy(Lr,5.13)]"))


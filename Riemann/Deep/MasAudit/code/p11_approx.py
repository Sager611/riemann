# P11: the §1.4.5 approximation phi_z(t) = t^{-z} 1_{t>1}: M(phi_z)(rho-1/2) = 1/(z-(rho-1/2)).
# Spectral side over zeta zeros: absolute sum diverges ~ log^2 T, symmetric sum converges.
import cypari2; gp=cypari2.Pari(); gp.allocatemem(100*10**6); gp.set_real_precision(30)
gp("Z=lfunzeros(lfuncreate(1),1000); z=3/2;")
for T in [100,300,1000]:
    print("T=%5d  #zeros=%5d  sum_abs=%.4f  sum_sym=%.8f" % (T, int(gp(f"#select(t->t<{T},Z)")),
          float(gp(f"sum(k=1,#Z,if(Z[k]<{T},2/abs(z-I*Z[k]),0))")), float(gp(f"sum(k=1,#Z,if(Z[k]<{T},2*real(1/(z-I*Z[k])),0))"))))

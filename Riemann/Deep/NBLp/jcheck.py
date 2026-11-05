# /// script
# dependencies = ["mpmath"]
# ///
# Check the exact residue formula of A4(b): J = int R(s')/|s|^2 dt on Re s = s' = 1/p - eps/2, R = zeta(s)/zeta(s+eps)
#   J = 2 pi [ R(2s')/(2s') - 1/((2s'-1) zeta(1+eps)) ]   (any zero location with Re rho <= 1/p; here true zeta)
# and compute B(eps) = (1/2pi) int |1-R|^2/|s|^2 dt (truncated) vs the bound K_p eps.
import mpmath as mp, time
mp.mp.dps = 15
t0=time.time()
for (p,eps) in [(mp.mpf(8)/7, mp.mpf('0.1')), (mp.mpf(8)/7, mp.mpf('0.05')), (mp.mpf('1.5'), mp.mpf('0.1'))]:
    sp = 1/p - eps/2
    R = lambda t: mp.zeta(mp.mpc(sp,t))/mp.zeta(mp.mpc(sp+eps,t))
    T = 300
    pts = [0]+[k*2.5 for k in range(1,int(T/2.5)+1)]
    Jnum = 2*mp.quad(lambda t: mp.re(R(t))/(sp**2+t**2), pts)
    # tail estimate beyond T: R ~ 1 + O(eps log t), so tail ~ 2*(1/T)
    Jtail = 2/mp.mpf(T)
    Jform = 2*mp.pi*( mp.zeta(2*sp)/mp.zeta(2*sp+eps)/(2*sp) - 1/((2*sp-1)*mp.zeta(1+eps)) )
    Bnum = mp.quad(lambda t: abs(1-R(t))**2/(sp**2+t**2), pts)/mp.pi   # (1/2pi)*2*int_0^T
    Kp = mp.mpf('5993.0534')/(2*mp.pi) + 4*p/(2-p)
    print("p=%s eps=%s sigma'=%s: J_num(0..T)+tail~%s  J_formula=%s | B_num(trunc)=%s  K_p*eps=%s  B/eps=%s  [%.0fs]" % (
        mp.nstr(p,5), mp.nstr(eps,3), mp.nstr(sp,5), mp.nstr(Jnum+Jtail,8), mp.nstr(Jform,8), mp.nstr(Bnum,6), mp.nstr(Kp*eps,6), mp.nstr(Bnum/eps,5), time.time()-t0), flush=True)

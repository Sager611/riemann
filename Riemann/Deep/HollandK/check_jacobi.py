
import mpmath as mp, random
mp.mp.dps = 50
def roots_q(U,V,d):
    # coefficients of 2F1(-d,U;V;y/U) = sum_j (-1)^j C(d,j) (U)_j/(U^j (V)_j) y^j
    co = [(-1)**j*mp.binomial(d,j)*mp.rf(U,j)/(U**j*mp.rf(V,j)) for j in range(d+1)]
    r = mp.polyroots(co[::-1], maxsteps=400, extraprec=400)
    return sorted([mp.re(x) for x in r]), max(abs(mp.im(x)) for x in r)
def jac(U,V,d):
    H = U-d-1; al = V-1; be = U-V-d
    Mt = mp.matrix(d,d)
    for k in range(d):
        Mt[k,k] = U*(H*(V+2*k)+2*k*(k+1))/((H+2*k)*(H+2*k+2))
    for k in range(1,d):
        o = U/(2*k+H)*mp.sqrt(k*(k+al)*(k+be)*(k+H)/((2*k+H-1)*(2*k+H+1)))
        Mt[k-1,k] = o; Mt[k,k-1] = o
    ev = mp.eigsy(Mt)[0]
    return sorted([ev[i] for i in range(d)]), Mt
random.seed(1)
for trial in range(6):
    d = random.randint(3,9)
    V = mp.mpf(random.randint(32*d, 80*d)) + mp.mpf(random.random())
    U = V + d + mp.mpf(random.randint(0, 500))
    r1, im = roots_q(U,V,d)
    r2, Mt = jac(U,V,d)
    err = max(abs(a-b) for a,b in zip(r1,r2))
    diag = max(abs(Mt[k,k]-V) for k in range(d)); off = max(abs(Mt[k-1,k]) for k in range(1,d))
    print('d',d,'U',mp.nstr(U,8),'V',mp.nstr(V,8),'root err',mp.nstr(err,3),'imag',mp.nstr(im,3),
          ' |diag-V|/d',mp.nstr(diag/d,4),' off/sqrt(Vd)',mp.nstr(off/mp.sqrt(V*d),4),
          ' maxdev/sqrt(Vd)', mp.nstr(max(abs(x-V) for x in r1)/mp.sqrt(V*d),4))

# check the exact diagonal-difference formula (jacobi-diagonal-difference)
for trial in range(5):
    d = random.randint(3,9); V = mp.mpf(random.randint(32*d, 80*d)) + mp.mpf(random.random()); U = V + d + mp.mpf(random.randint(0,500))
    H = U-d-1
    worst = 0
    for k in range(d):
        dg = U*(H*(V+2*k)+2*k*(k+1))/((H+2*k)*(H+2*k+2))
        fm = (2*k*(k+H+1)*(U-2*V) + V*H*(d-1))/((H+2*k)*(H+2*k+2))
        worst = max(worst, abs(dg - V - fm))
    print('diag-difference formula max error', mp.nstr(worst, 3))

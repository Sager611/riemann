
import mpmath as mp
mp.mp.dps = 30
def kap(Qc, R, th):
    Q = lambda y: mp.polyval(Qc[::-1], y)
    dQc = [k*Qc[k] for k in range(1, len(Qc))]
    dQ = lambda y: mp.polyval(dQc[::-1], y)
    W0 = mp.quad(lambda y: mp.e**(2*R*y)*Q(y)**2, [0, 1]); W1 = mp.quad(lambda y: mp.e**(2*R*y)*(R*Q(y) + dQ(y))**2, [0, 1])
    c = (1 + mp.e**(2*R)*Q(1)**2)/2 + W0/th + th*W1/3
    return 1 - mp.log(c)/R
for th in ['2/3', '1/2', '1/4', '1/6', '1/8']:
    t = mp.mpf(eval(th)); R = mp.sqrt(mp.mpf(3)/5)/t
    print(f"theta={th}: Levinson Q=1-y, P=x, R=sqrt(3/5)/theta={mp.nstr(R,6)}: kappa={mp.nstr(kap([1,-1], R, t), 8)}")

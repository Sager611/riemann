import mpmath as mp
mp.mp.dps = 20
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1); chi = [0, 1, k, -k, -1]
DH = lambda s: mp.dirichlet(s, chi)
tot = 0
for g in [mp.mpc(0.8085, 85.6993), mp.mpc(0.6508, 114.1633), mp.mpc(0.5744, 166.4793), mp.mpc(0.7243, 176.7025)]:
    try:
        r = mp.findroot(DH, g); w = 2*(r.real - 0.5)/r.imag**2; tot += w
        print('zero', mp.nstr(r, 10), '|DH|', mp.nstr(abs(DH(r)), 3), 'Carleman weight 2(beta-1/2)/gamma^2 =', mp.nstr(w, 5))
    except Exception as e:
        print('fail', g, e)
print('sum', mp.nstr(tot, 5))

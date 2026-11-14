from mpmath import mp, mpf, sqrt, findroot
mp.dps = 40
B2 = (777 - sqrt(1081)) / 852
P = lambda t: 3888 * mpf(19)**5 * (15*t - 13)**4 * (1 - t)**5 - 7 * (7 - 6*t)**4 * (9*t - 7)**5
th = findroot(P, mpf('0.8733458561'))
def Bsig(s):
    r = sqrt(964*s**2 - 1156*s + 1081)
    return (103528*s**2 + 7316*s*r + 51710*s - 3823*r - 51257) / (12*(8720*s**2 + 640*s*r + 4912*s - 356*r - 5179))
s = findroot(lambda s: Bsig(s) - th, mpf('0.0003'))
w = 2*(6*th - 7)/(9*(12*th - 11)); d = (45*w - 14)/(6*(9*w + 1))
print('theta_inf =', mp.nstr(th, 30)); print('B_2 - theta_inf (lambda_eq) =', mp.nstr(B2 - th, 12)); print('sigma_eq =', mp.nstr(s, 12))
print('contact width w* =', mp.nstr(w, 15), ' delta* =', mp.nstr(d, 15), ' count 5/9+delta/6 =', mp.nstr(mpf(5)/9 + d/6, 15))
print('theta check (14+25d)/(18+24d) =', mp.nstr((14 + 25*d)/(18 + 24*d), 30))

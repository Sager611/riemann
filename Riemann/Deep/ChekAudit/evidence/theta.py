import mpmath as mp
mp.mp.dps = 60
def Phi(u):
    s = mp.mpf(0)
    for n in range(1, 30):
        s += (2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2) - 3*mp.pi*n**2*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n**2*mp.e**(2*u))
    return 2*s
def Theta(c):
    c = mp.mpf(c); a = mp.log(c)/2; kap = 2*mp.pi*c - mp.mpf(9)/2
    return (2*mp.pi*Phi(a)/kap)**4
paper = {2:'1.870055251e-3',3:'3.726005778e-8',5:'4.863756675e-18',10:'7.153054191e-44',16:'1.296812318e-75',22:'1.126704465e-107',23:'4.897398523e-113',25:'9.005302038e-124',31:'4.712163115e-156',45:'1.183307824e-231',47:'1.786944521e-242',71:'1.461291502e-372'}
for c,v in paper.items():
    t = Theta(c)
    print(c, mp.nstr(t, 12), v, 'rel', mp.nstr(t/mp.mpf(v)-1, 3))
for c in [mp.e**mp.mpf('3.1'), mp.e**mp.mpf('3.2'), mp.mpf(5), mp.mpf(16)]:
    print('c=', mp.nstr(c, 8), 'Theta=', mp.nstr(Theta(c), 6))

import random
def T(r,z): return (14-18*r-12*z)/(30*(3*r-2))
mx=0
for _ in range(2000):
    r=random.uniform(0.67,0.7778); z=random.uniform(-0.01,0.02)
    zp=(z-T(r,z))/r
    dT=(36*z-6)/(30*(3*r-2)**2) + (-2/(5*(3*r-2)))*zp
    w=(33*r-14)/48 - z
    claim=-16*w/(25*r*(3*r-2)**2)
    wp=33/48 - zp; claimw=(3*r-8/5)/(r*(3*r-2))*w
    mx=max(mx,abs(dT-claim)/(1+abs(dT)),abs(wp-claimw)/(1+abs(wp)))
print("max rel err of both identities:",mx, " w(7/9)=",(33*7/9-14)/48, 35/144)

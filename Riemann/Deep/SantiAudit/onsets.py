# Onsets: least theta with h(theta;kappa)>0 for each certified detector, under three moment ranges.
from mpmath import mp, mpf, cot, sqrt, findroot
mp.dps = 40
det = [("199/10000","0.01427299118820739777907682971367218570370"),("299/10000","0.02144830107160461829756580252369722582753"),
       ("339/10000","0.02431764197988143277823702075131968426128"),("349/10000","0.02503497665228033293076056772439978810048"),
       ("399/10000","0.02862164964681450058004672694436292391142"),("229/5000","0.03285392367255858089368599993099127340889"),
       ("499/10000","0.03579499548791817334727167601167804381953"),("999/10000","0.07166400518009452424270000076357507802514")]
c = lambda t: 2 - t/2 - cot(t/sqrt(2))/sqrt(2)
def theta_h(k):
    cstar = 2 - 2/(1-k)
    return findroot(lambda t: c(t) - cstar, mpf('0.52'))
print("theta0 (c=0):", findroot(c, mpf('0.55')))
print("c(1/2) =", c(mpf(1)/2), " k*(1/2) =", 1-2/(2-c(mpf(1)/2)))
for nu_s, k_s in det:
    a,b = map(int, nu_s.split('/')); nu = mpf(a)/b; k = mpf(k_s)
    th = theta_h(k)
    sec2 = max(th, nu + mpf(1)/2)          # Sec 2 only: nu < theta - 1/2
    nuA  = max(th, (1+nu)/2, (6*nu+2)/5)   # paper's nuA
    triv = max(th, (4*nu+1)/3, mpf(1)/2)   # Steuding / trivial dual: nu < (3theta-1)/4
    print(nu_s, "kappa/nu=%.6f" % (k/nu), "theta_h=%.7f" % th, "sec2=%.7f" % sec2, "nuA=%.7f" % nuA, "triv=%.7f" % triv)
# optimal-detector heuristic kappa = 0.71733 nu
for name, rng in [("sec2", lambda t: t-mpf(1)/2), ("nuA", lambda t: min(2*t-1,(5*t-2)/6)), ("triv", lambda t: (3*t-1)/4)]:
    f = lambda t: mpf('0.71733')*rng(t) - (1-2/(2-c(t)))
    lo, hi = mpf('0.45'), mpf('0.56')
    # scan
    import mpmath
    ts = [lo + (hi-lo)*i/2000 for i in range(2001)]
    pos = [t for t in ts if f(t) > 0]
    print(name, "heuristic onset (kappa=0.71733*nu_max):", ("%.5f" % pos[0]) if pos else None)

# 2510.14167: constant max_t f_a(t), f_a(t) = (t log t - (t-a) log(t-a) - a log a)/(t+1-a) (paper l.320),
# closed form: max = -log(1 - tau), tau in (0,1) with 1 - tau = tau^a. a = 1/2 (GRH), 0.4736 (Harman), 1 - theta141.
from mpmath import mp, mpf, log, findroot, sqrt
mp.dps = 30
theta = mpf('0.87334585610671857680')
def f(a, t): return (t*log(t) - (t-a)*log(t-a) - a*log(a))/(t+1-a)
for a in [mpf(1)/2, mpf('0.4736'), 1-theta]:
    tau = findroot(lambda x: 1 - x - x**a, mpf('0.3'))
    c = -log(1-tau); u = a/tau
    print('a=%s tau=%s u=%s c=%s f(u)=%s c/log2=%s' % (mp.nstr(a,12), mp.nstr(tau,12), mp.nstr(u,12), mp.nstr(c,15), mp.nstr(f(a,u),15), mp.nstr(c/log(2),12)))
print('log golden =', mp.nstr(log((1+sqrt(5))/2),15), ' sqrt(1-theta)=', mp.nstr(sqrt(1-theta),20))

# e1rate: on theta in [9/11, 1], is c(theta) > h(theta; k) and c > k for k up to 0.4090 (Conrey's 4/7 one-piece ceiling 0.40896)?
import mpmath as mp
mp.mp.dps = 20
s2 = mp.sqrt(2); c = lambda t: 2-t/2-mp.cot(t/s2)/s2
h = lambda t, k: (3+k-mp.sqrt((1-k)*(9-k-8*c(t))))/4
k = mp.mpf('0.4090'); worst = None
for i in range(0, 1001):
    t = mp.mpf(9)/11 + (1-mp.mpf(9)/11)*i/1000
    m = c(t) - max(h(t, k), (c(t)+2*k)/3)
    if worst is None or m < worst[0]: worst = (m, t)
print("min over theta in [9/11,1] of c - max(h(theta;0.409), (c+2*0.409)/3) =", mp.nstr(worst[0], 8), "at theta =", mp.nstr(worst[1], 6))

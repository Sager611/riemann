# bvirraudit: numeric falsification of Prop 9' kernel facts (D sec. 1, Steps 2, 3, 5, 7). Pure python.
import cmath, math, random
random.seed(1)
def Phi(s, eta):                 # Laplace transform of phi_eta = I_{eta/4}^{*4}: (sinh z / z)^4, z = s eta / 8
    z = s * eta / 8
    return (cmath.sinh(z) / z) ** 4 if abs(z) > 1e-12 else 1.0
worst = 0.0; n = 0
for eta in [0.25, 1e-2, 1e-4, 1e-8]:
    for _ in range(20000):
        sig = random.uniform(-40, 40); t = 10 ** random.uniform(-3, 9) * random.choice([-1, 1])
        s = complex(sig, t)
        bound = math.exp(abs(sig) * eta / 2) * min(1.0, (8 / (eta * abs(s))) ** 4)
        worst = max(worst, abs(Phi(s, eta)) / bound); n += 1
print("Step 2: max |Phi_eta(s)| / bound over %d samples = %.6f (must be <= 1)" % (n, worst))
# Step 3 sandwich: g_+ = 1_{J+} * phi_eta >= 1 on [a,b]; g_- <= 1_{(a,b)}; check on a grid via exact B-spline CDF
def bspline_cdf(u, eta):          # P(sum of 4 iid U[-eta/8, eta/8] <= u), exact (Irwin-Hall)
    w = eta / 4; y = (u + eta / 2) / w   # in units of box width, shifted to [0, 4]
    if y <= 0: return 0.0
    if y >= 4: return 1.0
    return sum((-1) ** k * math.comb(4, k) * (y - k) ** 4 for k in range(int(y) + 1)) / 24
def g(u, lo, hi, eta): return bspline_cdf(u - lo, eta) - bspline_cdf(u - hi, eta)
bad = 0
for X in [2.0, 5.0, 40.0]:
    a, b = X - math.log(2), X
    for eta in [0.25, 0.01]:
        for i in range(4001):
            u = a - eta + (b - a + 2 * eta) * i / 4000
            gp_ = g(u, a - eta / 2, b + eta / 2, eta); gm = g(u, a + eta / 2, b - eta / 2, eta)
            ind = 1.0 if a < u <= b else 0.0
            if not (gm <= ind + 1e-12 and ind <= gp_ + 1e-12 and gm >= -1e-12 and gp_ <= 1 + 1e-12): bad += 1
            if (gm > 1e-12) and not (a < u < b): bad += 1
print("Step 3: sandwich violations on grid =", bad)
# Step 5: |G_pm(1) - x/2| <= eta x
w5 = 0.0
for X in [2, 3, 10, 50]:
    x = math.exp(X)
    for eta in [0.25, 0.1, 1e-3, 1e-6]:
        P1 = (math.sinh(eta / 8) / (eta / 8)) ** 4
        for sgn in [1, -1]:
            ap, bp = X - math.log(2) - sgn * eta / 2, X + sgn * eta / 2
            G1 = (math.exp(bp) - math.exp(ap)) * P1
            w5 = max(w5, abs(G1 - x / 2) / (eta * x))
print("Step 5: max |G(1) - x/2|/(eta x) = %.4f (must be <= 1)" % w5)
# Step 7: trivial zeros: 0.95 * (1 + sum e^{-j*1.05}) <= 1.6 ; a - eta >= 1.05 at X = 2, eta = 1/4
print("Step 7: a - eta at X=2, eta=1/4: %.4f; 0.95*(1+sum e^-j) = %.4f" % (2 - math.log(2) - 0.25, 0.95 * (1 + 1 / (math.e - 1))))
# Thresholds: d^(1/10) >= a log d
for a in [640, 8, 6.4]:
    k = 1.0
    while 10 ** (k / 10) < a * k * math.log(10): k += 0.01
    print("d^(1/10) >= %g log d  <=>  log10 d >= %.2f" % (a, k))
print("1/(1 - T*) with T* = 0.87338196661320752:", 1 / (1 - 0.87338196661320752), " margin at a=7.9:", 7.9 * (1 - 0.87338196661320752) - 1)

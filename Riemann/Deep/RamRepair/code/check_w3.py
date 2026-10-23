
# check_w3.py -- formal proof that w(u) = sec^4(theta) * int_0^{D-u} g(v) g(v+u) dv on [0, D], g(v) = cos(k v - theta) - cos(theta),
# k = tan(theta), D = 2 theta/k.  Everything is written with X = e^{i theta}, Y = e^{i k u}; theta, u, k appear as further symbols,
# with k = tan(theta) = -i (X^2 - 1)/(X^2 + 1) substituted; the difference must vanish as a rational function in (X, Y, theta, u).
import sympy as sp
X, Y, th, u = sp.symbols('X Y theta u')
I = sp.I
k = -I*(X**2 - 1)/(X**2 + 1)
cth = (X + 1/X)/2; sth = (X - 1/X)/(2*I)
# int_0^{D-u} e^{i n (k v - theta)} dv = (e^{i n (theta - k u)} - e^{-i n theta})/(i n k)   (k D = 2 theta)
def E(n):  # e^{i n (theta - k u)} = X^n Y^{-n}
    return X**n*Y**(-n)
def intE(n):
    if n == 0: return 2*th/k - u                     # length D - u
    return (E(n) - X**(-n))/(I*n*k)
# g(v) g(v+u) with g = (e^{i(kv-th)} + e^{-i(kv-th)})/2 - c ; shift by u multiplies e^{i n (k v - th)} by Y^n
terms = {}   # coefficient of e^{i n (k v - th)}
def add(n, coef): terms[n] = terms.get(n, 0) + coef
# (A1 + A_{-1})/2 - c times (Y A1 + Y^{-1} A_{-1})/2 - c, where A_n = e^{i n (kv - th)}
gv = {1: sp.Rational(1, 2), -1: sp.Rational(1, 2), 0: -cth}
gvu = {1: Y/2, -1: 1/(2*Y), 0: -cth}
for n1, c1 in gv.items():
    for n2, c2 in gvu.items():
        add(n1 + n2, c1*c2)
h = sum(coef*intE(n) for n, coef in terms.items())
s2 = 1/cth**2
cosku = (Y + 1/Y)/2
def sin_(a_exp_i):  # sin of angle whose e^{i angle} is given
    return (a_exp_i - 1/a_exp_i)/(2*I)
sin2th = sin_(X**2)
w = s2*(s2*(th*cth/sth - u/2)*cosku + 2*th*cth/sth - u + sin_(X**2/Y)/sin2th - 2*(1 + sin_(X/Y)/sth))
diff = sp.together(sp.expand(h*s2**2 - w))
num = sp.expand(sp.numer(diff))
print('numerator of (sec^4 h - w) after expansion:', num)
print('IDENTITY PROVED' if num == 0 else 'NOT ZERO')

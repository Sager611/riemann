# splice.py (bothachain): close the step c_i -> c_w with CCMFMM's certified single window c_w (35 or 40; ../CCMFMM/logs/int.jsonl) instead of new
# certificates. CCMFMM certified at c_w: lambda_1^odd(c_w) >= L_odd and lambda_2^even(c_w) >= L_2 (both rigorous, nu > kappa eta), U(c_w) >= lambda_1^even(c_w).
# If c_i < c_w and U(c_i) < min(L_odd, L_2), monotonicity gives for every c in [c_i, c_w]:
#   lambda_1^even(c) <= U(c_i) < min(L_odd, L_2) <= min(lambda_1^odd, lambda_2^even)(c_w) <= min(lambda_1^odd, lambda_2^even)(c),
# the same step inequality as chain_ccm.py with eps := min(L_odd, L_2). Appends the splice record; the chain resumes from c_w with U(c_w).
# usage: splice.py <chain.jsonl> <c_w>
import sys, json, math
from fractions import Fraction as F
def lg(s): s = str(s).strip('[').split(' ')[0]; m, _, e = s.partition('e'); return math.log10(float(m)) + int(e)
out, cw = sys.argv[1], sys.argv[2]
W = {d['c']: d for d in map(json.loads, open(__file__.rsplit('/', 1)[0] + '/../CCMFMM/logs/int.jsonl'))}[cw]
R = [d for d in map(json.loads, open(out)) if d.get('ok')]; last = R[-1]
ci, Ui = last['c'], last['U_up']
assert W['even_simple'] and W['defl']['rigorous'] and W['odd']['rigorous'], 'window not certified'
assert float(W['defl']['nu_cert']) > W['defl']['kappa_eta'] and float(W['odd']['nu_cert']) > W['odd']['kappa_eta']
Ls = [W['defl']['eps'], W['odd']['eps']]; eps = min(Ls, key=lg)
assert F(ci) < F(cw), 'chain already past the window'
ok = lg(Ui) < lg(eps)
rec = dict(c=cw, cprev=ci, Uprev=Ui, eps=eps, NG=W['NG'], bits=W['bits'], K2=W['K2'], ok=ok, U_up=W['U_up'], U=W['U'], lemmaG=True,
           splice='CCMFMM logs/int.jsonl c = %s: L_2 = %s, L_odd = %s certified (nu %s, %s; kappa eta %.1e, %.1e), U from its gstate NG = %d' % (
               cw, W['defl']['eps'], W['odd']['eps'], W['defl']['nu_cert'], W['odd']['nu_cert'], W['defl']['kappa_eta'], W['odd']['kappa_eta'], W['NG']),
           defl=W['defl'], odd=W['odd'], g_rss=W.get('g_rss'), g_sec=W.get('g_sec'), sec=0)
print(json.dumps(dict(c=cw, cprev=ci, Uprev=Ui, eps=eps, ratio=10**(lg(eps) - lg(Ui)), ok=ok)))
if ok:
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')

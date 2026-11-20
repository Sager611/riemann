R='/Users/adriansager/playground/riemann/'
T=open(R+'TREE.md',encoding='utf-8').read()
o="they hold given the positive selector, none claims RH, R507, R518)"
n=("they hold given the positive selector, none claims RH, R507, R518; version 5 (24400244, 20 Nov 15:10 UTC) adds §§55–61 on short and long excursions at ordinary primes (Lemmas 55.1, 56.1, 57.2, Props. 56.2, 57.1, 58.1, 59.1, 61.1, Cor. 57.3, Thm 60.1, a cofinal RH-equivalent through sharp short-interval prime sums with s = (14/13)(1 − Θ), Guth–Maynard input): all hold as checked, none claims RH, R528)")
assert T.count(o)==1
open(R+'TREE.md','w',encoding='utf-8').write(T.replace(o,n,1)); print('ok')

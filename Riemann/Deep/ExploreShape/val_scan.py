import sys, json, ast; sys.path.insert(0,'.')
from validate import minform
res = json.load(open('logs/scanJ_r2_d24.json'))
for line in res:
    for J in (4,5,6):
        a = line['a%d'%J]; ph = {ast.literal_eval(k): v for k,v in line['phi%d'%J].items()}
        out=[]
        for N,smin in ((16,1e-6),(30,1e-6),(30,1e-3)):
            ev,g = minform(a,ph,J,N,300,smin=smin)
            out.append('N%d smin%.0e: %.3g @g=%.3g' % (N,smin,ev,g))
        print('kappa %.2f J%d  ' % (line['kappa'],J) + ' | '.join(out), flush=True)

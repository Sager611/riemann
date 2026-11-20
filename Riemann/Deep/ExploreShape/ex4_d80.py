import sys, json, ast; sys.path.insert(0,'.')
from exact4 import check4
for line in json.load(open('logs/scanJ_r2_d80.json')):
    a=line['a4']; ph={ast.literal_eval(k):v for k,v in line['phi4'].items()}
    print('d80', line['kappa'], 'h4=%.6f'%line['h4'], 'h5=%.6f'%line['h5'], check4(a, ph[(1,2)], ph[(3,2)], ng=300, KM=20000), flush=True)

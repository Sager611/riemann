
# mkcells_hp141.py -- pbpaudit: one cell per k = 3..6 closing [stated decimal L0, first double of HP141's theta141 scan)
import json, math
from decimal import Decimal
d = json.load(open('Riemann/Deep/HP141/Ingham3/logs/cert_result.json'))['main']
cells = []
for k in (6, 5, 4, 3):
    L0 = d[str(k)]['L0']; dec = Decimal(str(L0)); dbl = Decimal(L0)
    La = math.nextafter(L0, 0.0)
    assert Decimal(La) < dec               # the cell [La, L0] contains the decimal, so [decimal, L0) is covered
    cells.append(dict(k=k, La=La, Lb=L0, tag='hp141k%d' % k, decimal=str(dec), double=str(dbl), sliver=str(max(dbl - dec, Decimal(0)))))
    print(k, 'decimal', dec, 'double', dbl, 'sliver', max(dbl - dec, Decimal(0)), 'cell La', Decimal(La))
json.dump(cells, open('Riemann/Deep/Ingham3Audit/logs/cells_hp141.json', 'w'), indent=1)


import math, numpy as np, sys
sys.path.insert(0,'Riemann/Deep/ExploreLK')
from lk_cone import base_moments, column, target
PI=math.pi
for dl in [0.3,0.5,0.7,1.0]:
    xs=np.array([0.0,dl]); m,d=base_moments(xs,2.0,2)
    a=math.cos(PI*dl/2)
    print(dl, 'd/m1', d/m[1], '2sin^2', 2*math.sin(PI*dl/2)**2, 'm1/m0', m[1]/m[0], (PI**2/2)/(a*a+0.5))
# sine check: zeros at 0,1 period 2
m,d=base_moments(np.array([0.0,1.0]),2.0,2); print('sine', m/m[0], d/m[0], 'expect', [1,PI**2,PI**4], 2*PI**2)
print('target', target(2))

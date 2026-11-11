import numpy as np, math
d=np.loadtxt('zeros_6670_6750.txt'); n=d[:,0].astype(int); z=2*d[:,1]
i=int(np.where(n==6709)[0][0]); zm,zp=z[i],z[i+1]; m=np.ones(len(z),bool); m[[i,i+1]]=False
D=math.log(zm/(4*math.pi))/(4*math.pi)
g=np.sum(1/(z[m]-zm)**2+1/(z[m]-zp)**2) + 2*D*(1/(zm-z[0]+0.9)+1/(z[-1]-zm+0.9))  # window + smooth tail
Dl=zp-zm; x=Dl*Dl*g
lam=((1-5*x/4)**0.8-1)/(8*g)
print('CSV (z-units): Delta=%.6f g=%.5f Delta^2 g=%.6f lambda_CSV=%.6e' % (Dl,g,x,lam))

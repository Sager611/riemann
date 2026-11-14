import mpmath as mpm
import numpy as np

from mpmath import mpf, nstr, fprod, fsum, power, factorial, fmul, fdiv, exp, sqrt, pi, log, findroot, zeta, euler
from numpy import inf
#from scipy.optimize import differential_evolution, NonlinearConstraint
from scipy.integrate import quad
from scipy.special import polygamma, factorial
from scipy.optimize import optimize, minimize, NonlinearConstraint, differential_evolution
from sympy import Symbol, solveset, S, erf, div, solve

##############

mpm.mp.dps = 10 #50

# ---- BellottiFix (2026-11-14): corrected copy; every change is marked FIX ----
import os, json
MODE = os.environ.get('BFIX_MODE', 'fix')  # FIX modes: fix = c3 x2 + rigorous |G(1/2+it)| (c4 as printed, its 2.4 checked);
                                            # r162 = only c4 x2 (earlier reading); both = c3 x2, c4 x2, rigorous G; c3only; c5only
K3 = 2 if MODE in ('fix', 'both', 'c3only') else 1   # FIX (4.11): representatives on both sides of gamma_r
K4 = 2 if MODE in ('r162', 'both') else 1            # FIX (4.13) variant: doubling on top of the printed 2.4
RIGG = MODE in ('fix', 'both', 'c5only')             # FIX l.933: |G(1/2+it)| <= c^2 P W/log^2 W, P = prod_p(1+(sqrt p-1)/p^2)
G16P = float(os.environ.get('BFIX_G', '9.7572'))     # FIX: >= c^2 P (sqfree.py: c^2 P <= 9.7571; crude 16P <= 23.5604, prodP.py)

##############

def latex_float(num, we):
    strng = format(float(num),we)
    if "E" in strng:
        a, b = strng.split("E")
        return "$" + str(float(a)) + "\cdot 10^{" + str(int(b)) + "}" + "$"


t = 3*power(10,12)
tupp= exp(29)
u=2.1781287
x=7.0798370
v=4.4744028
w=0.4808273
sigma=0.9927
T=3*pow(10,12)
t2=exp(46.2)
t3=exp(170.2)
t4=exp(481958)
err=8.066394225e+8390

#number of zeros in circle

r=power(2,0.5)*(1-sigma)

#constant b2

def constna(sigma):
    return fdiv(r*2*4.7908,4)+fdiv(2*r*log(1+fdiv(1,r))*(fdiv(4,pi)-1),1+2*r)+2*r*2.6908+2+fdiv(64*power(r,2),power(1+2*r,2))

CN2=constna(sigma)

#constant b1

def constinfrontint(sigma):
    return fdiv(fdiv(2*power(2,0.5)*log(t+power(2,0.5)*(1-sigma)),4)+fdiv(32*power(r,2),power(1+2*r,2))*log(t+power(2,0.5)*(1-sigma))-fdiv(16*power(r,2),power(1+2*r,2)),log(t))

CN1=constinfrontint(sigma)

#contribution square of divisor function VW

def varX(v,w):
    
    return power(t,v+w)

#constant d_1,1

def contrd2(u, x, v, w):
    
    return fdiv(varX(v,w)*power(log(varX(v,w)),3)*power(pi,-2)+0.746*varX(v,w)*power(log(varX(v,w)),2)+0.825*varX(v,w)*log(varX(v,w))+0.462*varX(v,w)+9.73*power(varX(v,w),3/4)*log(varX(v,w))+0.73*power(varX(v,w),0.5),varX(v,w)*power(log(varX(v,w)),3))

#contribution fraction of square of divisor function d_2,1

def fracd2(u, x, v, w):
    
    return fdiv(contrd2(u,v,w,x)*power(log(varX(v,w)),3)+fdiv(contrd2(u,v,w,x)*power(log(varX(v,w)),4),4),power(log(varX(v,w)),4))

#W^2

def varXw(w):
    
    return power(t,2*w)

#constant d_1,2

def contrd2w(u, x, v, w):
    
    return fdiv(varXw(w)*power(log(varXw(w)),3)*power(pi,-2)+0.746*varXw(w)*power(log(varXw(w)),2)+0.825*varXw(w)*log(varXw(w))+0.462*varXw(w)+9.73*power(varXw(w),3/4)*log(varXw(w))+0.73*power(varXw(w),0.5),varXw(w)*power(log(varXw(w)),3))

#contribution fraction of square of divisor function d_2,1

def fracd2w(u, x, v, w):
    
    return fdiv(contrd2w(u,v,w,x)*power(log(varXw(w)),3)+fdiv(contrd2w(u,v,w,x)*power(log(varXw(w)),4),4),power(log(varXw(w)),4))


#contribution Psi2 - constant d_3

def Psisquare(u, x, v, w):
    
    return fdiv((fdiv(4,power(v-u,2)*power(log(t),2)))*(1+1.0061*(u*log(t)+1.333+fdiv(11,power(u*log(t)*power(t,u),0.5)))),log(t))

#contribution PSiTheta

def PsiTheta(u, x, v, w):
    
    return fdiv((fdiv(2,(v-u)*w*power(log(t),2)))*(1+1.0061*(w*log(t)+1.333+fdiv(11,power(w*log(t)*power(t,w),0.5)))),log(t))

#estimate for G(1) constant d_5

def thsquare(u,x,v,w):
    return fdiv(1.0061*(w*log(t)+1.333+fdiv(3.95,power(t,w/2))), w*log(t))

#zero free regions

def nu1(t):
    return fdiv(1,5.558691*log(t))

def J(t):
    return fdiv(log(t),6)+log(log(t))+log(0.618)

def R(t):
    return fdiv(J(t)+ 0.685 +0.155 *log(log(t)),log(t)*(0.04962 -fdiv(0.0196,J(t)+1.15)))

def nu2(t):
    
    return fdiv(1, R(t)*log(t))

def nu3(t):
    return fdiv(log(log(t)),21.233 *log(t))

def nu4(t):
    return fdiv(1,53.989 *power(log(t),2/3)*power(log(log(t)),1/3))

#estimate for Lemma 4.1 first factor d_4

def Boundfirstfactor(u,x,v,w):
    return 3.09*fdiv((u+v)*(1.301*(1+power(fdiv(v,u),2))+1.084*(fdiv(v,u)+1)-0.116),(v-u)*(fdiv(v,u)-1))+fdiv(Psisquare(u, x, v, w),exp(power(t,u+v-x))*log(t))+fdiv(Psisquare(u, x, v, w),log(t))*(2+(x-(u+v))*log(t))

#stirling near zero
integralfromgamma = lambda y: np.power(y,2-2*sigma)*exp(-y)
integralfromgamma = quad(integralfromgamma, 1, np.inf)[0]

print('Intgamma:', integralfromgamma)

def BoundfromGamma1(u,x,v,w):
    return (0.632121+integralfromgamma)


#estimate for Lemma 4.1 when diff <delta

# FIX (4.9): X^y/y is increasing only for y >= 1/log X; Re z >= 2 nu(T) and a = 2 nu(T) x log T is < 1 in rows 26-38, so a term
# capped by X^{2delta}/(2delta) may need the factor FA = e^{a-1}/a; capped terms: s = r and the two nearest per side (NOTE 1)
def nuT0():
    if t3 <= t and t < t4: return nu3(t)
    elif t2 <= t and t < t3: return nu2(t)
    elif 0 < t and t < t2: return nu1(t)
    return nu4(t)
aMono = 2*x*log(t)*nuT0()
FA = exp(aMono-1)/aMono if aMono < 1 else mpf(1)

def constc2(u,x,v,w):
    return (1+5*(FA-1) if K3 == 2 else 1)*fdiv(thsquare(u,x,v,w),w)*BoundfromGamma1(u,x,v,w)*Boundfirstfactor(u,x,v,w)

def Creqs1(u, x, v, w):
    return constc2(u,x,v,w)*(CN1+fdiv(CN2*power(nu1(t),-1),log(t)))


def Creqs2(u, x, v, w):
 return constc2(u,x,v,w)*(CN1+fdiv(CN2*power(nu2(t),-1),log(t)))

def Creqs3(u, x, v, w):
 return constc2(u,x,v,w)*(CN1+fdiv(CN2*power(nu3(t),-1),log(t)))

def Creqs4(u, x, v, w):
 return constc2(u,x,v,w)*(CN1+fdiv(CN2*power(nu4(t),-1),log(t)))

#estimate for Lemma 4.1 when r different from s and imaginary part > delta and less than 1

def numbterm(u,x,v,w):
    if t3 <= t and t < t4:
        return log(power(nu3(tupp),-1))
    elif t2 <= t and t < t3:
        return log(power(nu2(tupp),-1))
    elif 0 <= t and t < t2:
        return log(power(nu1(tupp),-1))
    elif t >= t4:
        return log(power(nu4(tupp),-1))

def constc3(u,x,v,w):
    return K3*fdiv(2*thsquare(u,x,v,w),w)*BoundfromGamma1(u,x,v,w)*Boundfirstfactor(u,x,v,w)*numbterm(u,x,v,w)*1.12

def Creqs11(u, x, v, w):
    return constc3(u,x,v,w)*(CN1+fdiv(CN2*power(nu1(t),-1),log(t)))


def Creqs21(u, x, v, w):
 return constc3(u,x,v,w)*(CN1+fdiv(CN2*power(nu2(t),-1),log(t)))

def Creqs31(u, x, v, w):
    return constc3(u,x,v,w)*(CN1+fdiv(CN2*power(nu3(t),-1),log(t)))

def Creqs41(u, x, v, w):
    return constc3(u,x,v,w)*(CN1+fdiv(CN2*power(nu4(t),-1),log(t)))

#estimate for Lemma 4.1 when r different from s and imaginary part greater than 1

integralfromgammastir = lambda y: np.power(y,1.5-2*sigma)*exp(fdiv(-pi*y,2))*exp(fdiv(1,6*y))
integralfromgammastir = quad(integralfromgammastir, 1, np.inf)[0]

def BoundfromGamma2(u,x,v,w):
    
    return power(2*pi,0.5)*integralfromgammastir

def constc4(u,x,v,w):
    return K4*fdiv(2*thsquare(u,x,v,w),w)*BoundfromGamma2(u,x,v,w)*Boundfirstfactor(u,x,v,w)*2.4

def Creqs12(u, x, v, w):
    return constc4(u,x,v,w)*(CN1+fdiv(CN2*power(nu1(t),-1),log(t)))


def Creqs22(u, x, v, w):
    return constc4(u,x,v,w)*(CN1+fdiv(CN2*power(nu2(t),-1),log(t)))

def Creqs32(u, x, v, w):
    return constc4(u,x,v,w)*(CN1+fdiv(CN2*power(nu3(t),-1),log(t)))

def Creqs42(u, x, v, w):
    return constc4(u,x,v,w)*(CN1+fdiv(CN2*power(nu4(t),-1),log(t)))

#Contribution integrals in Lemma 4.1

#integral from -infty to -3-a J1

integralJ1= lambda y: np.power(2*y,27/164)*exp(fdiv(-pi*y,2))*power(y,1-2*sigma)
integralJ1 = quad(integralJ1, 3, np.inf)[0]

def BoundJ1(u,x,v,w):
    
    return power(2*pi,0.5)*exp(fdiv(1,6*(2*sigma-1.5)))*66.7*integralJ1

 
def BoundJ2(u,x,v,w):
     return 1.461*4.05206

#integral from  a to infty  J4

integralJ4= lambda y: np.power(2*y,27/164)*exp(fdiv(-pi*y,2))*power(y,1-2*sigma)
integralJ4 = quad(integralJ4, 3, np.inf)[0]

def BoundJ4(u,x,v,w):
    
    return power(2*pi,0.5)*exp(fdiv(1,6*(2*sigma-1.5)))*66.7*integralJ4

#integral from  3-a to a  J3 division in J31, J32, J33

integralJ33= lambda y: np.power(fdiv(y,3)+1,27/164)*exp(fdiv(-pi*y,2))*power(y,1-2*sigma)
integralJ33 = quad(integralJ33, 1, np.inf)[0]

def BoundJ33(u,x,v,w):
    
    return power(2*pi,0.5)*exp(fdiv(1,6*(2*sigma-1.5)))*66.7*integralJ33


integralJ31= lambda y: np.power(fdiv(y,3)+1,27/164)*exp(fdiv(-pi*y,2))*power(y,1-2*sigma)
integralJ31 = quad(integralJ31, 1, np.inf)[0]

def BoundJ31(u,x,v,w):
    
    return power(2*pi,0.5)*exp(fdiv(1,6*(2*sigma-1.5)))*66.7*integralJ31

#J32 used with sigma 0.985

def BoundJ32(u,x,v,w):
    return 66.7* 1.0486*3.62642
 

#contribution for the single integral  of zeta and gamma in Lemma 4.1

def Boundintzetagamma(u,x,v,w):
    
    return BoundJ1(u,x,v,w)+BoundJ2(u,x,v,w)+BoundJ4(u,x,v,w)+(BoundJ31(u,x,v,w)+BoundJ32(u,x,v,w)+BoundJ33(u,x,v,w))*power(t,fdiv(27,164))

#total contribution of integral term for fixed r and s

def factorbeforeint(u,x,v,w):
    Gfac = fdiv(G16P, w*log(t)) if RIGG else thsquare(u,x,v,w)   # FIX: replaces d_5 in the bound for |G(1/2+it)| log W / W
    return fdiv(Boundfirstfactor(u,x,v,w)*Gfac,2*pi)*fdiv(power(t,w),w*log(t))*(power(t,x*(1.5-2*sigma))+power(fdiv(power(t,u),power(log(t),2)),1.5-2*sigma))


#constant c5

def condition2(u,x,v,w):
    
    return (factorbeforeint(u,x,v,w)*Boundintzetagamma(u,x,v,w))
 
 #Lower Bound for zero detector
 #error term given by I2 and I3

integralerrorterm = lambda y: np.power(1+fdiv(y,0.447*log(t)),27/164)*exp(fdiv(-pi*y,2))*power(y,-sigma)
integral = quad(integralerrorterm, 1, np.inf)[0]

I2I3=2*(122.092+exp(1/6)*power(2*pi,0.5)*66.7*integral)


#condition 1 (lower bound, Lemma 3.4)
def firstterm(u, x, v, w):
    return exp(-fdiv(pow(log(t),2),pow(t,u)))


def S1(u, x, v, w):
    return fdiv( fracd2(u, x, v, w)*pow(u+w,4)*pow(2*pi,0.5)*exp(fdiv(1,6*0.447 *log(t)))*pow(log(t),4)*pow(t,x*(1-sigma)),exp(fdiv(0.447*pi*log(t),2))*pow(0.447*log(t),sigma-1+0.5))


def SU(u, x, v, w):
    return fdiv(fracd2(u, x, v, w)*I2I3,2*pi)*pow(v+w,4)*pow(t,fdiv(v+w,2)+x*(0.5-sigma)+27/164)*pow(log(t),4)


def cond1(u, x, v, w):
    return exp(-fdiv(pow(log(t),2),pow(t,u)))-fdiv(fracd2(u, x, v, w)*pow(u+w,4)*pow(2*pi,0.5)*exp(fdiv(1,6*log(t)))*pow(log(t),4)*pow(t,x*(1-sigma)),exp(fdiv(0.447*pi*log(t),2))*pow(0.447*log(t),sigma-1+0.5))-fdiv(fracd2(u, x, v, w)*I2I3,2*pi)*pow(v+w,4)*pow(t,(v+w)/2+x*(0.5-sigma)+27/164)*pow(log(t),4)-fdiv(1,pow(10,30))



#Final constant in front of J^2
 
def condition4(u, x, v, w):
    
    return pow(cond1(u, x, v, w),2)-condition2(u, x, v, w)



#constant in front of J before dividing

#
def C1(u, x, v, w):
    return (Creqs1(u, x, v, w)+Creqs11(u, x, v, w)+Creqs12(u, x, v, w))

def C2(u, x, v, w):
    return (Creqs2(u, x, v, w)+Creqs21(u, x, v, w)+Creqs22(u, x, v, w))

def C3(u, x, v, w):
    return (Creqs3(u, x, v, w)+Creqs31(u, x, v, w)+Creqs32(u, x, v, w))

def C4(u, x, v, w):
    return (Creqs4(u, x, v, w)+Creqs41(u, x, v, w)+Creqs42(u, x, v, w))

def FinalC(u,x,v,w):
    if t3 <= t and t < t4:
        return fdiv(C3(u, x, v, w),condition4(u, x, v, w))
    elif t2 <= t and t < t3:
        return fdiv(C2(u, x, v, w),condition4(u, x, v, w))
    elif 0 < t and t < t2:
        return fdiv(C1(u, x, v, w),condition4(u, x, v, w))
    elif t >= t4:
        return fdiv(C4(u, x, v, w),condition4(u, x, v, w))

#Constants output

print('Constant C:', FinalC(u,x,v,w))
print('Constant B:', 2*x)
print('d_1,1:', contrd2(u, x, v, w))
print('d_1,2:', contrd2w(u, x, v, w))
print('d_2,1:', fracd2(u, x, v, w))
print('d_2,2:', fracd2w(u, x, v, w))
print('d_3:',Psisquare(u, x, v, w))
print('d_5:',thsquare(u,x,v,w))
print('d_4:',Boundfirstfactor(u,x,v,w))
print('c_1:',cond1(u, x, v, w))
print('c_2:',constc2(u,x,v,w))
print('c_3:',constc3(u,x,v,w))
print('c_4:',constc4(u,x,v,w))
print('c_5:',condition2(u,x,v,w))
print('b2:', CN2)
print('b1:', CN1)

# FIX checks: (4.13) needs 2.4 >= 2 + 4 delta g(1)/I (both sides + two boundary terms per side), delta <= 1 - alpha0;
# (4.11) needs H_floor(1/nu) <= 1.12 log(1/nu), nu = nu_i(T1) (the per-side sum is H_floor(1/delta)/delta, delta > nu)
F4need = 2 + fdiv(4*(1-sigma)*exp(fdiv(1,6)-fdiv(pi,2)), integralfromgammastir)
nuT1 = exp(-numbterm(u,x,v,w))
Hn = mpm.harmonic(mpm.floor(fdiv(1,nuT1)))
print('extra:', json.dumps({'mode': MODE, 'F4need': float(F4need), 'ok4': bool(F4need <= 2.4), 'invnuT1': float(1/nuT1),
      'H': float(Hn), 'Hbound': float(1.12*log(1/nuT1)), 'ok3': bool(Hn <= 1.12*log(1/nuT1)),
      'C1': float(C1(u,x,v,w)), 'cond4': float(condition4(u,x,v,w)), 'Iint': float(integralfromgammastir), 'aMono': float(aMono), 'FA': float(FA)}))


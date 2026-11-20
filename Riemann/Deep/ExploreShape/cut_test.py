import sys, json; sys.path.insert(0,'.')
from sdpcut import sdp_cut
mom=json.load(open('logs/mom_r2_d24_k3.32.json')); kap=mom['kappa']; cen=mom['central_about_mean']
for J in (4,5):
    m=[cen[2*j]*kap**(2*j)/cen[0] for j in range(J+1)]
    val,a,ph,mn,rd=sdp_cut(m,J,verbose=True)
    print('J',J,'value',val,'minval',mn,'rounds',rd,'a',a,ph, flush=True)

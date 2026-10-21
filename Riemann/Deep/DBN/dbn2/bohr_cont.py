
import subprocess, shutil
def run(t,y2,load,save,it=300):
    out=subprocess.run(['nice','-n','10','./bohr','690988',repr(t),repr(y2),str(it),'1','0',load,save],capture_output=True,text=True).stdout
    d=dict(kv.split('=') for kv in out.split() if '=' in kv)
    return float(d['bestPhi'])
res={}
for y2 in [0.001,0.002,0.005,0.01,0.02,0.035,0.05]:
    t=0.128; load='-'; last_neg=None
    while t<0.16:
        p=run(t,y2,load,'th_tmp.bin')
        shutil.copy('th_tmp.bin','th_prev.bin'); load='th_prev.bin'
        if p<0: last_neg=(t,p)
        else:
            print('y2',y2,'last Phi<0 at t=',last_neg,'first Phi>=0 at t=',round(t,4),'Phi',p,'Lambda_obstr>=',None if last_neg is None else round(last_neg[0]+y2/2,5),flush=True)
            break
        t=round(t+0.001,4)

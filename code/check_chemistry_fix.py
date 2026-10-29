#!/usr/bin/env python
"""Run the corrected (v1.1) generator at the v1.0.1 adopted constants and seeds and
print the six-linkage distribution, phenolic ends, unreacted side chains, branch
coefficient and fragment count per wood; then exercise every diagnostic switch and
check determinism. Produced results/chemistry_fix_at_v101_constants.log (2026-09-14).
    conda run -n lignin python code/check_chemistry_fix.py
"""
import sys, time, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import mini_forward_generator as mfg, calibrate_full as cf
P={'k_bo4':8.0,'k_b5':2.4904,'k_b1':0.2715,'k_bb':16.7825,'k_55':0.7108,'k_4o5':0.1737,'n_nuclei':4.716,'k_bo4_nuc':0.9213}
K=cf.ALL_KEYS
print("%-9s %6s %6s %6s %6s %6s %6s | %5s %5s %6s %6s %5s" % ("wood","bO4","b5","bb","b1","5-5","4O5","RMSE","ends","CA/100","branch","frags"))
for sg,name in [(0.0,"softwood"),(1.0,"poplar"),(3.0,"high-S")]:
    agg={k:0 for k in K}; ends=[];ca=[];br=[];fr=[]
    for r in range(30):
        g=mfg.MiniForwardKMC(sg,n_target=120,params=P,seed=7000+r+int(sg*100)); s=g.run()
        for k in K: agg[k]+=g.bonds[k]
        ends.append(s["n_phenolic_ends"]); ca.append(s["unreacted_sidechain_per100"]); br.append(s["branch_coeff"]); fr.append(s["n_frags"])
        assert s["n_units"]==120 and not g.wasted
    f=cf.full_fracs(agg); rmse=np.sqrt(np.mean([(100*(f[k]-cf.REFERENCE[sg][k]))**2 for k in K]))
    print("%-9s %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f | %5.2f %5.1f %6.1f %6.3f %5.2f" % ((name,)+tuple(100*f[k] for k in K)+(rmse,np.mean(ends),np.mean(ca),np.mean(br),np.mean(fr))))
a=mfg.MiniForwardKMC(0.0,120,P,seed=5).run(); b=mfg.MiniForwardKMC(0.0,120,P,seed=5).run(); assert a==b; print("deterministic: OK")
for extra in [{"k_ao4":0.5},{"oligo_bb":0.3},{"oo_route":1},{"oo_route":2},{"k_nuc":0.5},{"nuclei_poisson":1},{"end_restricted":0},{"nucleation_on":0},{"pct_h":0.1,"frac_pca":0.3,"frac_tricin":0.15},{"k_55_nuc":1.0}]:
    p=dict(P); p.update(extra); s=mfg.MiniForwardKMC(0.8,120,p,seed=3).run(); print("  ok", extra, "n_units",s["n_units"],"frags",s["n_frags"],"bo4 %.1f"%(100*s["frac"]["bo4"]))
t=time.time(); mfg.MiniForwardKMC(0.0,480,P,seed=1).run(); print("480 units: %.2f s"%(time.time()-t))

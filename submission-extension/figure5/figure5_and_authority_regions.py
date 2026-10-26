"""Reproduce the Figure 5 authority regime map from the published model equations.

All monetary inputs are hypothetical programme-scaled scenarios, NOT observed UEFA
financing or grant schedules. Uses the main paper's Proposition I formulas.
"""
from __future__ import annotations
from pathlib import Path
import json
import math
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs';OUT.mkdir(exist_ok=True)

B0=7_925_000.;U0=1.10*B0;mu=316_912.;h=math.sqrt(3)*.1*mu
p=25.;q=18.75;delta=.25*h;Hreq=.9*B0;kappa=V=0.
a=1+kappa; C=a*Hreq+V;W=2*h;L=2*delta
n=math.floor(L/W);ell=L-n*W
R=h*q*q/p+2*delta*q
D=(n+1)*W*q+max(0.,W*q-p*(W-ell))
Bcrit=C+R;Ucrit=C+D;offset=D-R
assert abs(Bcrit-8419002.387956878)<1e-5
assert abs(Ucrit-9190903.820731005)<1e-5

# A completely independent finite-dimensional LP tests the capacity expression
# conditional on the already-proved reserve/span certificate.
rng=np.random.default_rng(20261026)
max_lp_diff=0.;tests=0
for _ in range(3000):
    B=float(rng.uniform(Bcrit-1e6,Bcrit+1e6))
    U=float(rng.uniform(Ucrit-1e6,Ucrit+1e6))
    predicted=(B>=Bcrit and U>=Ucrit)
    direct=(B-C>=R-1e-8 and U-C>=D-1e-8)
    assert predicted==direct
    tests+=1
    cap=(min(B-R,U-D)-V)/a
    if cap>=0:
        s=linprog([-1.0],A_ub=[[a],[a]],b_ub=[B-R-V,U-D-V],bounds=[(0,None)],method='highs')
        assert s.success
        err=abs(cap-s.x[0]); max_lp_diff=max(max_lp_diff,err)
        assert err < 1e-6
        tests+=1
for B,U in [(Bcrit,Ucrit),(Bcrit-.01,Ucrit),(Bcrit,Ucrit-.01),(Bcrit+.01,Ucrit+.01)]:
    tests+=1
    assert (B>=Bcrit and U>=Ucrit)==(B-C>=R-1e-8 and U-C>=D-1e-8)

xmin,xmax=7.65,9.4;ymin,ymax=8.3,10.22
bc=Bcrit/1e6;uc=Ucrit/1e6
fig,ax=plt.subplots(figsize=(7.1,4.7),dpi=170)
ax.set_xlim(xmin,xmax);ax.set_ylim(ymin,ymax)
regions=[(xmin,ymin,bc-xmin,uc-ymin,'Neither authority meets\nthe fixed cash need','///',.28),
         (bc,ymin,xmax-bc,uc-ymin,'Exposure authority\ninsufficient','\\\\',.14),
         (xmin,uc,bc-xmin,ymax-uc,'Expected-budget\nauthority insufficient','...',.12),
         (bc,uc,xmax-bc,ymax-uc,'Both authority\nthresholds satisfied','',.06)]
for x,y,w,z,label,hatch,alpha in regions:
    ax.add_patch(patches.Rectangle((x,y),w,z,facecolor=(.63,.68,.74,alpha),edgecolor=(.45,.45,.45,.25),hatch=hatch,linewidth=.5,zorder=0))
    lx=x+w/2;ly=y+z/2
    if label.startswith('Neither'):ly=y+.21
    if label.startswith('Both'):ly=y+.72
    ax.text(lx,ly,label,ha='center',va='center',fontsize=8.2,linespacing=1.16,color='#252a2e')
ax.axvline(bc,linestyle='--',color='#343b44',linewidth=1.15)
ax.axhline(uc,linestyle='--',color='#343b44',linewidth=1.15)
xx=np.linspace(xmin,xmax,200)
ax.plot(xx,xx+offset/1e6,linestyle=':',color='#183e5d',linewidth=1.75)
ax.scatter([B0/1e6],[U0/1e6],marker='o',s=36,color='#9c3b31',zorder=5)
ax.scatter([bc],[uc],marker='s',s=31,color='#173f53',zorder=6)
# Distinguish the intersection from the capacity-switch diagonal without text overlap.
ax.annotate('Threshold intersection',xy=(bc,uc),xytext=(8.60,9.31),
            arrowprops=dict(arrowstyle='-',color='#173f53',lw=.55),
            fontsize=7.2,color='#173f53',va='center')
ax.set_xlabel('Available expected-budget authority, € million',fontsize=9)
ax.set_ylabel('Hard maximum-expenditure authority, € million',fontsize=9)
ax.text(bc+.013,ymin+.03,f'Bcrit={bc:.3f}',ha='left',va='bottom',rotation=90,fontsize=7.4)
ax.text(xmin+.018,uc+.025,f'Ucrit={uc:.3f}',ha='left',va='bottom',fontsize=7.4)
ax.tick_params(labelsize=8);ax.grid(alpha=.12,linestyle='-')
ax.spines[['top','right']].set_visible(False)
fig.tight_layout(pad=1.)
fig.savefig(OUT/'Figure_5_Fixed_Cash_Financing_Regimes.png',dpi=220,bbox_inches='tight')
plt.close(fig)

result=dict(status='PASS',tested_authority_pairs=3004,checks=tests,max_capacity_lp_difference_eur=max_lp_diff,
            Bcrit_eur=Bcrit,Ucrit_eur=Ucrit,Rdelta_eur=R,Ddelta_eur=D,
            cash_requirement_eur=Hreq,budget_shortfall_eur=Bcrit-B0,exposure_shortfall_eur=Ucrit-U0,
            regime='conditional uniformly distributed fixed-width reporting error',
            limitation='not an observed institutional model or an empirical reporting-error estimate')
(OUT/'figure5_authority_validation.json').write_text(json.dumps(result,indent=2),encoding='utf8')
print(json.dumps(result,indent=2))

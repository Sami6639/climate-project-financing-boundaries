"""Deterministic reproduction of the bounded-expenditure analysis.
Run from any directory: python /path/to/code/reproduce.py
No empirical reporting errors or recipient cash observations are simulated.
"""
from pathlib import Path
import csv, json, math, platform, time, hashlib
import numpy as np
import scipy
from scipy.optimize import linprog, brentq
from scipy.special import ndtr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
START=time.perf_counter()
P=json.loads((ROOT/'data'/'parameters.json').read_text())
B=P['budget_eur'];mu=P['reference_tonnes'];p=P['rate_cap_eur_per_tonne'];q=P['target_eur_per_tonne']
phi=lambda z:np.exp(-np.asarray(z)**2/2)/np.sqrt(2*np.pi)
g=lambda z:phi(z)+np.asarray(z)*ndtr(z)
def uniform_max(mean,cap,floor=0,h=1,p=1,ceiling=False):
 if mean<floor or cap<floor or (not ceiling and mean>cap):return float('nan')
 vals=[p,math.sqrt(max(0,p*(mean-floor)/h)),(cap-floor)/(2*h)]
 if not ceiling:vals.append(math.sqrt(max(0,p*(cap-mean)/h)))
 return min(vals)
def gaussian_max(mean,cap,floor=0,sigma=1,p=1,ceiling=False):
 R=cap-floor
 if mean<floor or R<0 or (not ceiling and mean>cap):return float('nan')
 if ceiling:mean=min(mean,(cap+floor)/2)
 if mean==floor or mean==cap or R==0:return 0.
 k=R/(p*sigma)
 z=brentq(lambda z:p*sigma*(g(z)-g(z-k))-(mean-floor),-40,40+k,xtol=1e-13)
 return float(p*(ndtr(z)-ndtr(z-k)))
def lp_max(mean,cap,n,distribution='uniform',ceiling=False):
 # Baseline plus unrestricted piecewise-constant slopes; no ramp imposed.
 x=np.linspace(-1,1,n+1) if distribution=='uniform' else np.linspace(-12,12,n+1)
 dx=np.diff(x)
 if distribution=='uniform':
  weights=dx/2; mean_weights=((1-x[:-1])**2-(1-x[1:])**2)/4
 else:
  weights=np.diff(ndtr(x));mean_weights=g(-x[:-1])-g(-x[1:])
 objective=np.r_[0,-weights];cost=np.r_[1,mean_weights];exposure=np.r_[1,dx]
 if ceiling:r=linprog(objective,A_ub=np.vstack([exposure,cost]),b_ub=[cap,mean],bounds=[(0,cap)]+[(0,1)]*n,method='highs')
 else:r=linprog(objective,A_ub=exposure[None,:],b_ub=[cap],A_eq=cost[None,:],b_eq=[mean],bounds=[(0,cap)]+[(0,1)]*n,method='highs')
 if not r.success:raise RuntimeError(r.message)
 return -float(r.fun)
def csvout(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
h=math.sqrt(3)*mu*P['relative_sd'];reserve=h*q*q/p;span=2*h*q;H=.9*B;U=1.1*B
central={'uniform_halfwidth_tonnes':h,'mean_reserve_eur':reserve,'payment_span_eur':span,'advance_90pct_eur':H,'minimum_cap_for_90pct_eur':H+span,'cap_110pct_eur':U,'max_advance_eur':min(B-reserve,U-span),'max_advance_percent':100*min(B-reserve,U-span)/B,'max_obligation_at_90pct':uniform_max(B,U,H,h,p),'exact_cap_ratio_gate':1+reserve/B,'no_cap_advance_percent':100*(B-reserve)/B,'span_sd_threshold':(U-H)/(2*math.sqrt(3)*mu*q),'mean_sd_threshold':(B-H)*p/(math.sqrt(3)*mu*q*q),'upper_mean_sd_threshold':(U-B)*p/(math.sqrt(3)*mu*q*q)}
rows=[]
for a in [1.05,1.1,1.2]:
 for k,V in [(0,0),(.02,25000),(.05,100000)]:
  hm=(min(B-reserve,a*B-span)-V)/(1+k)
  rows.append({'cap_ratio':a,'loading':k,'fixed_cost_eur':V,'exact_advance_eur':hm if B<=a*B-reserve else 'Infeasible','ceiling_advance_eur':hm})
csvout('table3_advances.csv',rows)
rows=[]
for r in [.05,.075,.07700141167815729,.10,.10266854890420972,.125,.15,.20]:
 hh=math.sqrt(3)*mu*r;c=hh*q*q/p;s=2*hh*q;hm=min(B-c,U-s)
 rows.append({'relative_sd':r,'mean_reserve_eur':c,'span_eur':s,'mean_test_90pct':B-H>=c-1e-7,'span_test_90pct':U-H>=s-1e-7,'exact_regime_exists':U-B>=c-1e-7,'max_ceiling_advance_percent':100*hm/B})
csvout('table5_dispersion.csv',rows)
rows=[]
for typ,mean,cap in [('Lower mean',.1,2),('Upper mean',1.9,2),('Peak span',.375,.75)]:
 for ceil in [False,True]:rows.append({'case':typ,'regime':'Ceiling' if ceil else 'Exact','lower_margin':mean-.25,'upper_margin':'Not required' if ceil else cap-mean-.25,'span_margin':cap-1,'maximum_obligation':uniform_max(mean,cap,ceiling=ceil)})
csvout('table2_constraint_failures.csv',rows)
# Verification grid, deterministically enumerated in the saved CSV.
ratios=[0,.005,.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99,.995,1]
checks=[]
for n in [16,32,64,128,256,1024]:
 for cap in [.25,.5,.75,1,1.5,2]:
  for a in ratios:
   m=a*cap;formula=uniform_max(m,cap);lp=lp_max(m,cap,n)
   checks.append({'n':n,'cap':cap,'mean':m,'formula':formula,'lp':lp,'formula_minus_lp':formula-lp})
csvout('uniform_lp_audit.csv',checks)
ceiling=[]
for cap in [.5,1,2,4]:
 for a in [0,.01,.1,.25,.5,.75,1,1.25]:
  m=a*cap;f=uniform_max(m,cap,ceiling=True);l=lp_max(m,cap,1024,ceiling=True)
  ceiling.append({'cap':cap,'budget':m,'formula':f,'lp':l,'formula_minus_lp':f-l})
csvout('ceiling_lp_audit.csv',ceiling)
normal=[]
for cap in [.25,.5,1,2,4,8]:
 for a in [.01,.1,.25,.5,.9,.99]:
  m=a*cap;f=gaussian_max(m,cap);l=lp_max(m,cap,2048,'normal')
  normal.append({'cap':cap,'mean':m,'formula':f,'lp':l,'formula_minus_lp':f-l})
csvout('gaussian_lp_audit.csv',normal)
# Exact piecewise-linear integration of constructed uniform schedules.
construct=[]
for cap in [.5,1,2,4]:
 for a in [.05,.25,.5,.75,.95]:
  mean=a*cap;target=.8*uniform_max(mean,cap);D=2*target;w=D;c=target**2
  weight=(mean-c)/(cap-2*c)
  ramp=lambda y: (1-weight)*np.clip(y-1+w,0,w)+weight*(cap-D+np.clip(y+1,0,w))
  knots=np.unique([-1,1,1-w,-1+w]);knots=knots[(knots>=-1)&(knots<=1)]
  avg=float(np.sum(np.diff(knots)*(ramp(knots[:-1])+ramp(knots[1:]))/4))
  derivative=float((ramp(1)-ramp(-1))/2)
  construct.append({'cap':cap,'mean':mean,'target':target,'constructed_mean':avg,'mean_error':abs(avg-mean),'marginal_error':abs(derivative-target)})
csvout('constructive_audit.csv',construct)
summary={'uniform_lp_cases':len(checks),'uniform_max_gaps_by_n':{str(n):max(x['formula_minus_lp'] for x in checks if x['n']==n) for n in [16,32,64,128,256,1024]},'ceiling_cases':len(ceiling),'ceiling_max_gap':max(x['formula_minus_lp'] for x in ceiling),'gaussian_cases':len(normal),'gaussian_max_gap':max(x['formula_minus_lp'] for x in normal),'constructive_cases':len(construct),'constructed_mean_max_error':max(x['mean_error'] for x in construct),'constructed_marginal_max_error':max(x['marginal_error'] for x in construct)}
assert all(x['formula_minus_lp']>=-1e-7 for x in checks+ceiling+normal)
assert summary['uniform_max_gaps_by_n']['1024']<.001
assert summary['gaussian_max_gap']<.0001
assert summary['constructed_mean_max_error']<1e-10
# Figure series (all are scenarios, never fitted data).
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.25,'savefig.dpi':220})
blue='#17618a';orange='#d57b1d';red='#aa3030'
x=np.linspace(0,2.5,501);ce=np.array([uniform_max(v,2,ceiling=True) for v in x]);ex=np.array([uniform_max(v,2) for v in x])
fig,ax=plt.subplots(figsize=(8.5,4.6));ax.plot(x,ce,color=orange,label='Available expected budget');ax.plot(x,ex,'--',color=blue,label='Exact expected spending');ax.axvspan(2,2.5,color='grey',alpha=.12);ax.set(xlabel='Required or available expected expenditure',ylabel='Maximum local marginal obligation',ylim=(0,1.15));ax.legend();fig.tight_layout();fig.savefig(OUT/'figure1.png');plt.close(fig)
csvout('figure1_series.csv',[{'mean_or_budget':a,'exact_max':b,'ceiling_max':c} for a,b,c in zip(x,ex,ce)])
x=np.linspace(.9,1.4,301);fig,ax=plt.subplots(figsize=(8.5,4.8));f2=[]
for k,V,col,label in [(0,0,blue,'No donor costs'),(.05,100000,orange,'5% loading; EUR 100,000 fixed cost')]:
 y=(np.minimum(B-reserve,x*B-span)-V)/(1+k)/B*100
 ax.plot(x,y,color=col,label=label);ax.axhline((B-reserve-V)/(1+k)/B*100,color=col,ls=':')
 f2.extend({'cap_ratio':xx,'loading':k,'fixed_cost_eur':V,'ceiling_advance_percent':yy} for xx,yy in zip(x,y))
ax.axvspan(.9,central['exact_cap_ratio_gate'],color='grey',alpha=.12);ax.axvline(1.1,color='grey',ls='--');ax.scatter([1.1],[central['max_advance_percent']],color=blue,zorder=3);ax.scatter([1.1],[90],color=red,marker='x',s=65,zorder=4)
ax.annotate('84.03% maximum advance',xy=(1.1,central['max_advance_percent']),xytext=(1.17,78),arrowprops={'arrowstyle':'-','color':blue},color=blue)
ax.annotate('90% advance infeasible at U/B = 1.10',xy=(1.1,90),xytext=(1.13,97),arrowprops={'arrowstyle':'-','color':red},color=red)
ax.text(.92,95,'Exact parity\ninfeasible',color='dimgray');ax.set(xlabel='Hard donor-expenditure cap / budget scale',ylabel='Maximum early project cash (% of budget)',ylim=(58,102));ax.legend(loc='lower right');fig.tight_layout();fig.savefig(OUT/'figure2.png');plt.close(fig);csvout('figure2_series.csv',f2)
fig,ax=plt.subplots(figsize=(8.5,4.7));x=[0,1,2];ax.plot(x,[48.25,0,0],'-o',color=blue,label='Grant before supplier bill');ax.plot(x,[4.825,-43.425,0],'--o',color=red,label='Grant after supplier bill');ax.axhline(0,color='grey');ax.set_xticks(x,['Initial finance','Full supplier bill','Final transfer']);ax.set(ylabel='Cash balance before bridge finance (EUR thousands)',xlabel='Hypothetical event order; actual dates are unobserved',ylim=(-52,70));ax.text(.02,.95,'Reported cost EUR 48,250; grant EUR 43,425\nAssumed liquid own cash EUR 4,825; no credit',transform=ax.transAxes,va='top',fontsize=10);ax.legend(loc='lower left');fig.tight_layout();fig.savefig(OUT/'figure3.png');plt.close(fig)
csvout('figure3_series.csv',[{'event':a,'advance_cash_eur':b,'late_grant_cash_eur':c} for a,b,c in zip(['Initial finance','Full supplier bill','Final transfer'],[48250,0,0],[4825,-43425,0])])
fig,axs=plt.subplots(1,3,figsize=(10.5,4),sharey=True);f4=[]
for ax,R in zip(axs,[1,2,4]):
 x=np.linspace(0,1,201);uu=[uniform_max(a*R,R,h=math.sqrt(3)) for a in x];nn=[gaussian_max(a*R,R) for a in x]
 ax.plot(x,uu,color=blue,label='Uniform');ax.plot(x,nn,'--',color=orange,label='Gaussian');ax.set(title=f'Cap range / (rate × SD) = {R}',xlabel='Mean position within payment bounds',ylim=(0,1.04));ax.axvline(.5,color='grey',ls=':');f4.extend({'cap_to_rate_sd':R,'mean_fraction':a,'uniform_max':b,'gaussian_max':c} for a,b,c in zip(x,uu,nn))
axs[0].set_ylabel('Maximum obligation / rate cap');axs[1].legend(loc='lower center');fig.tight_layout();fig.savefig(OUT/'figure4.png');plt.close(fig);csvout('figure4_series.csv',f4)
(OUT/'central_results.json').write_text(json.dumps(central,indent=2))
(OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2))
runtime={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'seconds':time.perf_counter()-START,'stochastic_sampling':False,'seed':'Not applicable: no random draws'}
(OUT/'runtime.json').write_text(json.dumps(runtime,indent=2))
print(json.dumps({'central':central,'audit':summary,'runtime':runtime},indent=2))

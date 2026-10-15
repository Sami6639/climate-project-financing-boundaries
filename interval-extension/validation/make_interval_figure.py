import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from stress_formula import theory,OUT
h=np.sqrt(3)*.1*316912.; p=25.; q=18.75;B=7925000.;U=1.1*B
x=np.linspace(0,2,1601)
y=[];old=[]
for a in x:
 D,R=theory(a,q/p)
 H=min(B-p*h*R,U-p*h*D)
 Ho=min(B-q*h*(1+2*a),U-2*q*h*(1+a))
 y.append(100*H/B if H>=0 else np.nan)
 old.append(100*Ho/B if Ho>=0 else np.nan)
y=np.array(y);old=np.array(old)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(6.5,3.6),layout='constrained')
ax.fill_between(x,old,y,color='#17828C',alpha=.10,zorder=1)
ax.plot(x,y,color='#125A78',lw=2.1,label='Sharp interval-wide frontier',zorder=3)
ax.plot(x,old,color='#A4622B',lw=1.8,ls='--',label='Slope-q sufficient construction',zorder=2)
ax.set(xlim=(0,2),ylim=(28,90),xlabel=r'Interval half-width relative to reporting error ($\delta/h$)',ylabel='Maximum secured advance (% of budget)')
ax.set_xticks(np.arange(0,2.01,.25)); ax.set_yticks([30,40,50,60,70,80,90]);ax.grid(axis='y',alpha=.2)
ax.legend(loc='upper right',frameon=False)
fig.savefig(OUT/'Figure_3_Interval_Frontier.png',dpi=300)
fig.savefig(OUT/'Figure_3_Interval_Frontier.pdf')
np.savetxt(OUT/'figure_3_data.csv',np.c_[x,y,old],delimiter=',',header='delta_over_h,sharp_advance_percent,sufficient_advance_percent',comments='')

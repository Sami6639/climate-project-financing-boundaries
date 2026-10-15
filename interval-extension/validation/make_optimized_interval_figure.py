"""Independent audit cross-check: matched global, optimal one-ramp, conservative comparison."""
from pathlib import Path
import os,tempfile,csv,json
os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='interval-mpl-'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ramp_interval_checks import global_frontier,best_ramp,fixed_ramp
OUT=Path(__file__).resolve().parent

def run():
    h=np.sqrt(3)*.1*316912.;p=25.;q=18.75;B=7925000.;U=1.1*B
    alpha=np.unique(np.r_[np.linspace(0,2,801),.1,.25,.5,.75,1.,1.25,1.5])
    rows=[]
    for a in alpha:
        D,R=global_frontier(h,p,q,h*a);best=best_ramp(h,p,q,h*a,B,U);simple=fixed_ramp(h,p,q,h*a,q)
        sharp=min(B-R,U-D);conservative=min(B-simple['reserve'],U-simple['span'])
        rows.append(dict(delta_over_h=float(a),global_advance_eur=sharp,optimal_one_ramp_advance_eur=best['c'],conservative_slope_q_advance_eur=conservative,global_minus_optimal_eur=sharp-best['c'],global_minus_conservative_eur=sharp-conservative,optimal_slope_eur_per_tco2e=best['r']))
    with (OUT/'audit_crosscheck_data.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
    x=np.array([v['delta_over_h'] for v in rows]);sharp=np.array([v['global_advance_eur'] for v in rows]);best=np.array([v['optimal_one_ramp_advance_eur'] for v in rows]);simple=np.array([v['conservative_slope_q_advance_eur'] for v in rows])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,1,figsize=(7.1,6.4),sharex=True,gridspec_kw={'height_ratios':[2.1,1.]})
    axes[0].plot(x,sharp/1e6,label='Sharp class frontier',color='#174f78',lw=2.6,zorder=3)
    axes[0].plot(x,best/1e6,label='Optimal single ramp',color='#ce7326',lw=2.,ls='--',zorder=4)
    axes[0].plot(x,simple/1e6,label='Conservative slope-q ramp',color='#66865b',lw=1.8,ls=':',zorder=2)
    axes[0].set_ylabel('Maximum advance (€ million)');axes[0].legend(loc='upper right',frameon=False,fontsize=9)
    axes[1].plot(x,(sharp-best)/1e6,color='#ce7326',lw=2.,label='Sharp minus optimal single ramp')
    axes[1].plot(x,(sharp-simple)/1e6,color='#66865b',lw=1.8,ls=':',label='Sharp minus conservative ramp')
    axes[1].set_ylabel('Capacity difference (€ million)');axes[1].set_xlabel(r'Interval half-width / reporting half-width, $\delta/h$');axes[1].legend(loc='upper right',frameon=False,fontsize=8.5)
    for ax in axes:
        ax.grid(axis='y',color='#d8dee4',lw=.6);ax.set_xlim(0,2);ax.axvline(.25,color='#888888',lw=.8,ls='--',alpha=.7)
    axes[1].set_ylim(bottom=-.025)
    fig.tight_layout()
    fig.savefig(OUT/'audit_crosscheck_Three_Way_Frontier.png',dpi=300,bbox_inches='tight')
    fig.savefig(OUT/'audit_crosscheck_Three_Way_Frontier.pdf',bbox_inches='tight')
    chosen=[r for r in rows if any(abs(r['delta_over_h']-x)<1e-12 for x in [0,.1,.25,.5,.75,1,1.5,2])]
    (OUT/'audit_crosscheck_selected_values.json').write_text(json.dumps(dict(inputs=dict(h=h,p=p,q=q,B=B,U=U,kappa=0,V=0),rows=chosen),indent=2))
    print(json.dumps(chosen,indent=2));return chosen
if __name__=='__main__':run()

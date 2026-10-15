#!/usr/bin/env python3
"""Reproduce the central three-way frontier figure from exact formulas."""
from pathlib import Path
import argparse
import csv
import math
import os
import tempfile
_cache=Path(tempfile.gettempdir())/"interval_finance_plot_cache"
_cache.mkdir(parents=True,exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",str(_cache/"matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME",str(_cache))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from single_ramp_optimum import solve_single_ramp, class_frontier, conservative_frontier


def make_figure(out_dir):
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    mu=316912.;h=math.sqrt(3)*.1*mu;B=7925000.
    rows=[]
    for i in range(901):
        lam=i/400
        x=dict(h=h,delta=lam*h,p=25.,q=18.75,B=B,U=1.1*B,kappa=0.,V=0.,mu=mu)
        r=solve_single_ramp(**x);other={k:v for k,v in x.items() if k!='mu'}
        g=class_frontier(**other);c=conservative_frontier(**other)
        rg=g['H_algebraic']-r['H_algebraic'];cg=g['H_algebraic']-c['H_algebraic']
        if min(rg,cg)<-1e-7:raise RuntimeError('Frontier ordering violated')
        rows.append(dict(delta_over_h=lam,global_H=g['H_algebraic'],ramp_H=r['H_algebraic'],
                         conservative_H=c['H_algebraic'],global_minus_ramp=max(0,rg),
                         global_minus_conservative=max(0,cg),optimal_r=r['r']))
    with (out/'Figure_3_three_way_data.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':10,
                         'axes.titlesize':10,'legend.fontsize':8.5,'axes.spines.top':False,
                         'axes.spines.right':False,'savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,1,figsize=(7.2,6.5),sharex=True,
                          gridspec_kw={'height_ratios':[1.25,1]},layout='constrained')
    xs=[x['delta_over_h'] for x in rows]
    colors=['#17365D','#C15522','#777777']
    for key,label,color,style in [('global_H','Full-class frontier',colors[0],'-'),
                                 ('ramp_H','Best single ramp',colors[1],'--'),
                                 ('conservative_H','Slope-q construction',colors[2],':')]:
        axes[0].plot(xs,[100*x[key]/B for x in rows],label=label,color=color,ls=style,lw=2.0)
    axes[0].set_ylabel('Maximum advance (% of B)')
    axes[0].set_title('(a) Capacity under the same interval-wide requirements',loc='left',pad=9)
    axes[0].legend(loc='upper right',frameon=False)
    axes[0].set_ylim(20,89)
    axes[0].annotate('The best ramp attains the full frontier\nthrough δ/h = 0.25',
                     xy=(.18,83.8),xytext=(.08,58),fontsize=8.5,
                     arrowprops={'arrowstyle':'-','color':'#555555','lw':.8})
    axes[1].plot(xs,[x['global_minus_ramp']/1000 for x in rows],color=colors[1],lw=2,
                 label='Full class minus best ramp')
    axes[1].plot(xs,[x['global_minus_conservative']/1000 for x in rows],color=colors[2],lw=1.8,ls=':',
                 label='Full class minus slope-q construction')
    axes[1].set_ylabel('Advance difference (€ thousand)')
    axes[1].set_xlabel('Coverage half-width, δ/h')
    axes[1].set_title('(b) Comparator choice changes the measured advantage',loc='left',pad=9)
    axes[1].legend(loc='upper right',frameon=False)
    axes[1].set_ylim(-15,570)
    axes[1].axhline(0,color='#888888',lw=.7)
    for ax in axes:
        ax.grid(axis='y',color='#DADDE1',lw=.6)
        ax.set_xlim(0,2.25)
        ax.set_axisbelow(True)
        for boundary in [.25,1.,2.]:ax.axvline(boundary,color='#CED3D8',lw=.8,ls='--',zorder=0)
    axes[1].set_xticks([0,.25,.5,.75,1,1.25,1.5,1.75,2,2.25])
    fig.savefig(out/'Figure_3_Three_Way_Frontier.png',dpi=300)
    fig.savefig(out/'Figure_3_Three_Way_Frontier.pdf')
    plt.close(fig)
    print(f'PASS: saved {len(rows)} exact-formula points and PNG/PDF under {out}')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',default=str(Path(__file__).parent))
    make_figure(p.parse_args().output_dir)

"""Regenerate publication figures, numerical tables, and illustrated trace data."""
import json
import os
from pathlib import Path

MAN = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(MAN/'.mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from molecular_layout import projection, view_basis

E = MAN/'evidence'
F = MAN/'figs'
F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
    'axes.labelsize':9, 'axes.titlesize':10, 'axes.titleweight':'bold',
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,
    'svg.fonttype':'none','savefig.dpi':240})
GREEN='#007F68'; BLUE='#3269A8'; ORANGE='#CC6B26'; INK='#24363C'; GRAY='#ABB6BA'

def read(n): return json.loads((E/n).read_text())
def save(fig,n):
    for ext in ['pdf','svg','png']: fig.savefig(F/f'{n}.{ext}',bbox_inches='tight',facecolor='white')
    plt.close(fig)

seed=read('seed_comparison.json'); slap=read('slap_sweep.json'); methods=seed['methods']; N=1851
assert all(methods[k]['golden_cases']==N for k in ['seeds1','seeds2','seeds3','seeds10'])
assert len(seed['common_case_indices'])==1807
table=[]
for k,label in [('seeds1','1 seed order'),('seeds2','2 seed orders'),('seeds3','3 seed orders'),('seeds10','10 seed orders')]:
    d=methods[k]; r=d['golden_outcomes']
    table.append(f"{label} & {r['recovered']:,} & {d['golden_recovery_percent']:.2f} & {r['not_recovered']} & {r['unknown']} & {d['common_mean_cpu_seconds']:.2f} \\\\")
(MAN/'includes/generated-seed-table.tex').write_text(
    '\\begin{tabular}{lrrrrr}\\toprule\n'
    'Configuration & Recovered & \\% & Absent & Unknown & CPU s/reaction\\\\\\midrule\n'
    +'\n'.join(table)+'\n\\bottomrule\\end{tabular}\n')

# Figure 1: a conceptual diagram, explicitly separate from measured trace data.
fig=plt.figure(figsize=(9,4.5)); ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,10),ylim=(0,5));ax.axis('off')
def box(x,y,w,h,title,body,color=GREEN):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.09,rounding_size=.08',fc='#F2F7F6',ec=color,lw=1.1))
    ax.text(x+.13,y+h-.25,title,color=color,weight='bold',fontsize=10)
    ax.text(x+.13,y+h-.57,body,color=INK,fontsize=8.6,va='top',linespacing=1.5)
def arrow(a,b,color=INK): ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,color=color,lw=1.2))
ax.text(.15,4.85,'a  Continuous fragment growth',weight='bold',fontsize=12,color=INK)
box(.2,2.9,2.75,1.5,'01  Extend the frontier','Choose a strong source edge.\nCheck all active fragment bonds.\nCarry compatible placements.')
box(3.65,2.9,2.75,1.5,'02  Preserve distinctions','Canonicalize candidate states.\nKeep locked roles and boundaries.\nRetain correlated symmetry.')
box(7.1,2.9,2.65,1.5,'03  Continue or branch','Grow until saturation.\nFork distinct fragment placements.\nContinue under each prefix.')
arrow((2.98,3.65),(3.57,3.65));arrow((6.45,3.65),(7.02,3.65))
ax.annotate('repeat after each accepted extension',xy=(1.55,2.88),xytext=(5.0,2.5),ha='center',fontsize=8.5,
            arrowprops=dict(arrowstyle='->',connectionstyle='angle,angleA=180,angleB=-90,rad=15',color=GREEN),color=GREEN)
ax.text(.15,2.05,'b  Shared search history',weight='bold',fontsize=12,color=INK)
pos=[(.6,.95),(2,1.55),(2,.45),(3.4,.95),(4.8,.95)]
for a,b in [(0,1),(0,2),(1,3),(2,3),(3,4)]: arrow(pos[a],pos[b],BLUE)
for i,(x,y) in enumerate(pos):ax.plot(x,y,'o',ms=15,color=BLUE);ax.text(x,y,str(i),ha='center',va='center',color='white',fontsize=8)
ax.text(1.25,.05,'Fork',ha='center',fontsize=9);ax.text(3.5,.05,'Exact state join',ha='center',fontsize=9)
box(5.8,.28,3.95,1.58,'04  Compile retained families','Merge equal or contained families.\nKeep incomparable alternatives.\nPreserve every discovering path.',BLUE)
ax.text(.2,-.3,'Schematic: joins require equal assignments, island partitions and deferred edges at a compatible continuation.',fontsize=8,color=INK)
save(fig,'fig1_algorithm')

# Figure 2: full-denominator recovery and paired successful-work compute.
fig,axes=plt.subplots(1,2,figsize=(10,4.6),gridspec_kw={'width_ratios':[1.15,1]},layout='constrained')
labels=['SLAP: uncut control','SLAP: earlier expanded','SLAP: single-edge sweep','AAM: 1 seed + sweep','AAM: 2 seeds + sweep','AAM: 3 seeds + sweep','AAM: 10 seeds + sweep']
counts=[slap['uncut_recovered'],slap['expanded_baseline_recovered'],slap['sweep_union_recovered']]+[methods[k]['golden_outcomes']['recovered'] for k in ['seeds1','seeds2','seeds3','seeds10']]
colors=[ORANGE]*3+[GREEN]*4;y=np.arange(len(labels))
axes[0].barh(y,np.asarray(counts)/N*100,color=colors,height=.65)
for i,v in enumerate(counts):axes[0].text(v/N*100+1,i,f'{v:,}  ({v/N*100:.2f}%)',va='center',fontsize=8)
axes[0].set(yticks=y,yticklabels=labels,xlim=(0,125),xticks=[0,25,50,75,100],xlabel='Verified reference recovery (%)')
axes[0].invert_yaxis();axes[0].set_title('a  All 1,851 Golden records',loc='left');axes[0].grid(axis='x',alpha=.16);axes[0].set_axisbelow(True)
ax=axes[1]
for k,label in [('seeds1','1 seed'),('seeds2','2 seeds'),('seeds3','3 seeds'),('seeds10','10 seeds'),('slap_sweep','SLAP + sweep')]:
    d=methods[k]; x=d['common_mean_cpu_seconds']; y=d['golden_recovery_percent']; c=ORANGE if k=='slap_sweep' else GREEN
    ax.scatter(x,y,s=65,color=c,zorder=3)
    ax.annotate(label,(x,y),xytext=(0,10 if k!='seeds2' else -17),textcoords='offset points',ha='center',fontsize=8.5,color=c)
ax.set(xscale='log',xlim=(8,130),ylim=(96.5,99.65),xlabel='Search CPU seconds / reaction (log scale)',ylabel='Verified recovery on all 1,851 (%)')
ax.set_xticks([10,20,50,100],labels=['10','20','50','100']);ax.grid(alpha=.2);ax.set_title('b  Compute on 1,807 common completed cases',loc='left')
fig.get_layout_engine().set(rect=(0,.07,1,1))
fig.text(.02,.012,'Both directions. CPU excludes measured persistence/loading; instrumentation and interrupted-work accounting differ. No latency claim.',fontsize=8,color=INK)
save(fig,'fig2_golden')

# Figure 3: reference-free holdout, strict cap 100 only.
rows=read('holdout_scores.json'); ov=read('holdout_overlap.json')
assert len(rows)==140
x=np.array([r['slap_events'] for r in rows]); y=np.array([r['new_events'] for r in rows])
assert ((y==x).sum(),(y<x).sum(),(y>x).sum())==(133,6,1)
fig,axes=plt.subplots(1,3,figsize=(10,3.6),gridspec_kw={'width_ratios':[1.15,1,1]},layout='constrained')
from collections import Counter
ct=Counter(zip(x,y))
axes[0].plot([0,22],[0,22],color=GRAY,lw=1,ls='--')
for (a,b),n in ct.items():axes[0].scatter(a,b,s=24+18*np.sqrt(n),color=GREEN if b<=a else ORANGE,alpha=.8,edgecolors='white',linewidths=.5)
axes[0].annotate('case 123\nstrict cap 100',(3,19),xytext=(8,19),fontsize=8,arrowprops=dict(arrowstyle='-',color=GRAY))
axes[0].set(xlim=(-1,23),ylim=(-1,23),xlabel='SLAP: best saved event count',ylabel='AAM: best saved event count')
axes[0].set_title('a  133 ties; 6 AAM-lower; 1 SLAP-lower',loc='left',fontsize=9)
st=ov['target_class_statuses']; vals=[st['represented'],st['excluded_from_saved_families'],st['unresolved']]
axes[1].barh([2,1,0],vals,color=[GREEN,ORANGE,GRAY],height=.6)
for yy,n in zip([2,1,0],vals):axes[1].text(n+3,yy,str(n),va='center',fontsize=9)
axes[1].set(yticks=[2,1,0],yticklabels=['Represented','Excluded','Unresolved'],xlim=(0,185),xlabel='SLAP heavy mapping classes')
axes[1].set_title('b  Membership in saved AAM families',loc='left',fontsize=9)
wins=[0,1,2,5]; nums=[ov['aam_extra_cases_by_event_window'][str(k)] for k in wins]
axes[2].bar(range(4),nums,color=BLUE,width=.6)
for j,n in enumerate(nums):axes[2].text(j,n+3,str(n),ha='center',fontsize=9)
axes[2].set(xticks=range(4),xticklabels=['0','+1','+2','+5'],ylim=(0,155),xlabel='Events above best saved AAM score',ylabel='Cases with an extra AAM heavy pattern')
axes[2].set_title('c  Alternatives absent from saved SLAP',loc='left',fontsize=9)
fig.get_layout_engine().set(rect=(0,.08,1,1));fig.text(.015,.01,'140 XYZ/WBO cases, 10 seed orders, tolerance 1.0, cap 100. No reference mappings. Marker area reflects repeated score pairs.',fontsize=8)
save(fig,'fig3_holdout')

# Figure 4: concrete recovery under the separate publication engine.
coll=read('collection.json');fig,ax=plt.subplots(figsize=(6.7,3.6),layout='constrained')
for k,c,label in [('single',BLUE,'Smaller to larger'),('bidirectional',GREEN,'Both directions')]:
    yy=[coll['modes'][k]['event_windows'][str(w)] for w in range(7)]
    ax.plot(range(7),np.array(yy)/N*100,'o-',color=c,label=label,lw=1.8,ms=5)
    ax.annotate(f'{yy[-1]:,}/{N:,} ({yy[-1]/N*100:.2f}%)',(6,yy[-1]/N*100),xytext=(-5,9 if k=='bidirectional' else -15),textcoords='offset points',ha='right',fontsize=8.5,color=c)
ax.set(xlabel='Additional bond events above best collected candidate',ylabel='Concrete reference recovery (%)',xticks=range(7),ylim=(89.5,100.2))
ax.legend(loc='lower right',frameon=False);ax.grid(alpha=.2);ax.set_title('Earlier frozen publication engine: concrete candidate collection',loc='left',fontsize=10)
save(fig,'fig4_event_windows')

def normalize_trace(raw):
    frames=[];locked={};deferred=[];island=0;count=0;local={}
    for e in raw['events']:
        typ=e['type']
        if typ=='seed_start':
            island+=1; local={};count=e['init_cands'];patterns=e.get('cand_patterns',[])
            local=patterns[0]['witness'] if patterns else {}; title=f"Start fragment {island} at atom {e['seed']}"
        elif typ=='commit':
            patterns=e.get('cands_pattern_after',[]);local=e['cands_sample_after'][0];count=e['cands_after']
            title=f"Extend to {e['element']}{e['added']} · WBO {e['edge']['wbo']:.3f}"
        elif typ=='consumed':
            edge=[e['edge']['frag_atom'],e['edge']['ext_atom']];deferred.append(edge)
            title=f"Defer boundary {edge[0]}–{edge[1]} · continue growth"
            patterns=[]
        elif typ=='seed_end':
            if not e.get('iso'): continue
            local=e['iso'];count=e['final_cands'];patterns=e.get('cand_patterns',[])
            title=f"Fragment {island} saturated · {count} retained placement"+('s' if count!=1 else '')
        elif typ=='island_locked':
            locked.update({str(k):v for k,v in e['pairs']});continue
        else:continue
        mapping={**locked,**local};active=list(map(int,local))
        frames.append(dict(title=title,event=typ,fragment=island,candidates=count,mapping=mapping,
            active=active,locked=list(map(int,locked)),deferred=[list(x) for x in deferred],
            samples=[{**locked,**p['witness']} for p in patterns],
            highlight=e.get('added',e.get('seed')),edge=e.get('edge')))
    # Final alternatives come directly from terminal states, not trace samples.
    for n,state in enumerate(raw['graph']['states']):
        if state['id'] not in {s['state'] for s in raw['graph']['stops'] if s['reason']=='objective_met'}:continue
        frames.append(dict(title=f"Retained terminal {state['id']} · complete explicit-atom mapping",event='terminal',
            fragment=island,candidates=count,mapping=dict(state['mapping']),active=[r for r,p in state['mapping']],
            locked=[],deferred=state['deferred_edges'],samples=[],highlight=None,edge=None))
    return dict(case=raw['case'],name=raw['name'],input=raw['input'],frames=frames,graph=raw['graph'],scope=raw['scope'])

traces=[normalize_trace(read(n)) for n in ['growth_trace.json','growth_trace_1.json']]
for t in traces:
    for endpoint in ['reactant','product']:
        t['input'][endpoint]['view_basis']=view_basis(t['input'][endpoint]).tolist()
(E/'animation_data.json').write_text(json.dumps(traces,separators=(',',':'))+'\n')

def molecule(ax,ep,active=(),mapping=None,deferred=(),label=True):
    xy=projection(ep);w=np.asarray(ep['wbo']);active=set(active);deferred={tuple(sorted(e)) for e in deferred}
    for i,j in zip(*np.where(np.triu(w,1)>=.2)):
        color=ORANGE if (i,j) in deferred else GREEN if i in active and j in active else '#CFD7D9'
        ax.plot(xy[[i,j],0],xy[[i,j],1],color=color,lw=1.6 if color==GREEN else 1.0,ls='--' if color==ORANGE else '-',zorder=1)
    for i,(a,b) in enumerate(xy):
        is_h=ep['elements'][i]=='H';on=i in active
        ax.scatter(a,b,s=63 if is_h else 82,color=GREEN if on else '#E1E6E7',edgecolors='white',linewidths=.5,zorder=3)
        if label and (on or mapping is None):
            ax.text(a,b,str(mapping.get(i,i) if mapping is not None else i),fontsize=5.1,ha='center',va='center',color='white' if on else '#65777C',zorder=4)
    ax.set_aspect('equal');ax.axis('off');ax.margins(.18)

# Figure 5: genuine molecular snapshots with the measured candidate trajectory.
t=traces[0];selected=[0,6,14,17]
fig=plt.figure(figsize=(9,5.1));gs=fig.add_gridspec(3,4,height_ratios=[1,1,.7],hspace=.25,wspace=.12)
for j,i in enumerate(selected):
    f=t['frames'][i];ax=fig.add_subplot(gs[0,j]);molecule(ax,t['input']['reactant'],f['active']);ax.set_title(f"{len(f['active'])} atom"+('s' if len(f['active'])!=1 else '')+f" · {f['candidates']} candidates",fontsize=8.5)
    ax=fig.add_subplot(gs[1,j]);m=dict((int(a),b) for a,b in f['mapping'].items());molecule(ax,t['input']['product'],m.values(),{p:r for r,p in m.items()})
fig.text(.02,.75,'R',color=INK,weight='bold');fig.text(.02,.43,'P',color=INK,weight='bold')
ax=fig.add_subplot(gs[2,:]);ff=[f for f in t['frames'] if f['event'] in ('seed_start','commit')]
ax.plot([len(f['active']) for f in ff],[f['candidates'] for f in ff],'-o',color=GREEN,ms=3)
ax.set(xlabel='Atoms in the growing fragment',ylabel='Live candidates',xticks=[1,4,8,12,16,18]);ax.grid(alpha=.15)
fig.suptitle('Continuous growth on an 18-atom carbocation endpoint pair',fontsize=11,weight='bold',y=.99)
fig.subplots_adjust(left=.09,right=.98,bottom=.12,top=.9)
save(fig,'fig5_growth')
print('Generated five figures (PDF/SVG/PNG), seed table, and animation data.')

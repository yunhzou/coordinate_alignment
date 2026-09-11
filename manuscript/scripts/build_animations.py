"""Build the offline viewer and MP4/GIF movies from recorded events."""

from rxn_core.viewers import style_document
import json
import os
from pathlib import Path
import numpy as np
import imageio.v2 as imageio
from PIL import Image

MAN=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(MAN/'.mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from molecular_layout import projection
DATA=json.loads((MAN/'evidence/animation_data.json').read_text())
OUT=MAN/'animations';OUT.mkdir(exist_ok=True)
html=(MAN/'scripts/animation_template.html').read_text().replace('__DATA__',json.dumps(DATA).replace('</','<\\/'))
(OUT/'index.html').write_text(style_document(html, layout='growth_animation'))

def project(ep):
    return projection(ep)

def molecule(ax,ep,f,side):
    xy=project(ep);m={int(k):int(v) for k,v in f['mapping'].items()};inv={v:k for k,v in m.items()}
    active=set(f['active']);locked=set(f['locked']);deferred={tuple(sorted(e)) for e in f['deferred']}
    source=lambda i:i if side=='R' else inv.get(i)
    w=np.asarray(ep['wbo'])
    for i,j in zip(*np.where(np.triu(w,1)>=.2)):
        isdefer=side=='R' and (i,j) in deferred
        color='#c96b2c' if isdefer else '#007f68' if source(i) in active and source(j) in active else '#a9c0b9'
        ax.plot(xy[[i,j],0],xy[[i,j],1],color=color,lw=2 if color=='#007f68' else 1.25,ls='--' if isdefer else '-',zorder=1)
    for i,(x,y) in enumerate(xy):
        a=source(i);col='#dfb54e' if a==f['highlight'] else '#007f68' if a in active else '#3269a8' if a in locked else '#d7dfda'
        size=36 if ep['elements'][i]=='H' else 83
        ax.scatter(x,y,s=size,c=col,ec='white',lw=.7,zorder=3)
        if (len(xy)<=25 or ep['elements'][i]!='H') and (side=='R' or a is not None):
            ax.annotate(ep['elements'][i]+str(a if side=='P' else i),(x,y),xytext=(0,7),textcoords='offset points',ha='center',fontsize=7,color='#203b37')
    radius=np.max(np.linalg.norm(xy,axis=1))*1.2
    ax.set(xlim=(-radius,radius),ylim=(-radius,radius));ax.set_aspect('equal');ax.axis('off')
    ax.set_title('REACTANT R' if side=='R' else 'PRODUCT P · labels follow source identities',loc='left',fontsize=10,color='#617771',pad=12)

metadata=[]
for ci,t in enumerate(DATA):
    name='movie1_continuous_growth' if ci==0 else 'movie2_deferred_boundaries'
    frames=t['frames'];gif=[];movie=OUT/(name+'.mp4')
    with imageio.get_writer(movie,fps=10,codec='libx264',quality=8,macro_block_size=1,
                           ffmpeg_params=['-movflags','+faststart']) as writer:
        for i,f in enumerate(frames):
            fig=plt.figure(figsize=(12.8,7.2),dpi=100,facecolor='#f7faf7')
            fig.text(.04,.94,'GRAFT · CONTINUOUS FRAGMENT GROWTH',fontsize=11,color='#007f68',weight='bold')
            fig.text(.04,.88,f['title'],fontsize=18,color='#203b37',weight='bold')
            fig.text(.04,.835,t['name']+f"  ·  case {t['case']}  ·  event {i+1}/{len(frames)}",fontsize=10,color='#617771')
            ax=fig.add_axes([.035,.22,.445,.55]);molecule(ax,t['input']['reactant'],f,'R')
            ax=fig.add_axes([.53,.22,.445,.55]);molecule(ax,t['input']['product'],f,'P')
            fig.text(.045,.15,f"{len(f['mapping'])}/{len(t['input']['reactant']['elements'])} assigned atoms    |    fragment {f['fragment']}    |    {f['candidates']} live compressed candidates",fontsize=12,color='#203b37')
            for j,(c,label) in enumerate([('#007f68','Growing fragment'),('#3269a8','Earlier atoms'),('#dfb54e','New atom'),('#c96b2c','Deferred boundary')]):
                fig.text(.045+j*.23,.1,'● '+label,color=c,fontsize=10)
            fig.text(.045,.045,'Recorded algorithm events, not physical reaction time. One displayed witness; other candidates remain represented.',fontsize=9,color='#617771')
            # A measured algorithm-event progress bar, not elapsed runtime.
            fig.add_artist(Line2D([.04,.04+.92*(i+1)/len(frames)],[.018,.018],transform=fig.transFigure,color='#007f68',lw=4))
            fig.canvas.draw();rgb=np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy();plt.close(fig)
            hold=18 if f['event'] in ('seed_end','terminal','seed_start') else 7
            for _ in range(hold):writer.append_data(rgb)
            gif.append((Image.fromarray(rgb).resize((768,432)),hold*100))
        gif[0][0].save(OUT/(name+'.gif'),save_all=True,append_images=[im for im,d in gif[1:]],duration=[d for im,d in gif],loop=0,optimize=False)
        gif[-1][0].save(OUT/(name+'_poster.png'))
    metadata.append(dict(file=movie.name,case=t['case'],events=len(frames),duration_seconds=sum(d for im,d in gif)/1000,
                         scope=t['scope']))
    print(name,metadata[-1],flush=True)
(OUT/'movies.json').write_text(json.dumps(metadata,indent=2)+'\n')

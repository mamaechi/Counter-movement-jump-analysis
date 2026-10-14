import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .core import PERCENT, channel_label


def comparison(subject,signal,axis,data,unit,events=None,show_trials=True):
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True,gridspec_kw={'height_ratios':[3,1]})
    for limb,color in [('Left','#0072B2'),('Right','#D55E00')]:
        y=data.get((limb,signal,axis))
        if y is None: continue
        mean=y.mean(0); sd=y.std(0,ddof=1)
        axes[0].plot(PERCENT,mean,color=color,label=f'{limb} mean ± SD (n={len(y)})')
        axes[0].fill_between(PERCENT,mean-sd,mean+sd,color=color,alpha=.16)
        if show_trials:
            axes[0].plot(PERCENT,y.T,color=color,alpha=.25,lw=.65)
    l=data.get(('Left',signal,axis)); r=data.get(('Right',signal,axis))
    if l is not None and r is not None:
        axes[1].plot(PERCENT,r.mean(0)-l.mean(0),color='#333333',label='Right − Left')
        axes[1].legend(fontsize=8)
    else:
        axes[1].text(.05,.5,'No matched bilateral channel',transform=axes[1].transAxes)
    axes[1].axhline(0,color='gray',lw=.6)
    if events:
        for key,value in events.items():
            for ax in axes: ax.axvline(value,color='gray',ls=':',lw=.7)
            axes[0].text(value,1.01,key,rotation=45,fontsize=7,transform=axes[0].get_xaxis_transform())
    axes[0].set_title(f'{subject} · {channel_label(signal,axis)}',pad=38 if events else 8)
    axes[0].set_ylabel(unit); axes[1].set_ylabel('Difference')
    axes[1].set_xlabel('Normalized movement (%)'); axes[1].set_xlim(0,100)
    axes[0].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.18)
    fig.tight_layout()
    plt.close(fig)
    return fig


def spm_figure(title,inference):
    fig,ax=plt.subplots(figsize=(9,3.5))
    inference.plot(ax=ax)
    ax.set_title(title); ax.set_xlabel('Normalized movement (%)'); ax.set_ylabel('SPM{t}: Right − Left')
    fig.tight_layout()
    plt.close(fig)
    return fig


def phase_review_figure(subject,data,trial_index,trial_id,threshold_fraction,events=None):
    fig,axes=plt.subplots(2,1,figsize=(9,5),sharex=True)
    force=[]
    for limb,color in [('Left','#0072B2'),('Right','#D55E00')]:
        knee=data.get((limb,'Knee Angles','X')); grf=data.get((limb,'GRF','Z'))
        if knee is not None: axes[0].plot(PERCENT,knee[trial_index],label=f'{limb} knee X',color=color)
        if grf is not None:
            force.append(grf[trial_index]); axes[1].plot(PERCENT,grf[trial_index],label=f'{limb} GRF Z',color=color)
    if len(force)==2:
        total=np.sum(force,axis=0)
        axes[1].plot(PERCENT,total,label='Total GRF Z',color='black',lw=1)
        axes[1].axhline(threshold_fraction*total.max(),color='gray',ls='--',label='Candidate threshold')
    for key,value in (events or {}).items():
        for ax in axes: ax.axvline(value,color='gray',ls=':',lw=.7)
        axes[0].text(value,1.01,key,rotation=45,fontsize=7,transform=axes[0].get_xaxis_transform())
    axes[0].set_title(f'{subject} · trial {trial_id} · event review',pad=35)
    axes[0].set_ylabel('Knee X (source units)'); axes[1].set_ylabel('GRF (source units)')
    axes[1].set_xlabel('Normalized movement (%)')
    for ax in axes:
        ax.set_xlim(0,100); ax.grid(alpha=.2)
        if ax.lines: ax.legend(fontsize=7)
    fig.tight_layout(); plt.close(fig)
    plt.close(fig)
    return fig

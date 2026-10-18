"""Replace only the PSD data in quiet-event panels (d,e).

Display code for (a,b,c) is retained verbatim from the user-provided workflow.
Their original CWT arrays are loaded directly, without a new calculation.
"""
from pathlib import Path
import json
import sys
import os
PACKAGE=Path(os.environ.get("JUPITER_RUN_ROOT", str(Path(__file__).resolve().parents[1] / "output"))).resolve()
sys.path.insert(0,str(Path(__file__).resolve().parent))
import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,LogLocator,NullFormatter
import numpy as np
import pandas as pd
from display_curve import visible_curve
FONT=9
TITLE_FONT=10
PANEL_FONT=11


def sci_formatter(y,pos):
    if y<=0:
        return '0'
    exponent=int(np.floor(np.log10(y)))
    if -4<=exponent<0:
        return rf'$10^{{{exponent}}}$'
    coeff=y/(10**exponent)
    return rf'${coeff:.0f} \times 10^{{{exponent}}}$'


def add_utc_axis(ax,time_utc,label='UTC',top=True):
    twin=ax.twiny() if top else ax.twiny()
    twin.set_xlim(time_utc[0],time_utc[-1])
    twin.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=5,maxticks=10))
    twin.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    twin.set_xlabel(label,fontsize=10)
    if not top:
        twin.xaxis.set_ticks_position('bottom')
        twin.xaxis.set_label_position('bottom')
    return twin


def main():
    out=PACKAGE/'figures'
    qa=PACKAGE/'qa'
    out.mkdir(parents=True,exist_ok=True)
    qa.mkdir(parents=True,exist_ok=True)
    name='FigS3_PSD'
    display=np.load(PACKAGE/'data/quiet_display_arrays.npz')
    csv=pd.read_csv(PACKAGE/'data/quiet_psd.csv',float_precision='round_trip')
    utc=pd.DatetimeIndex(display['UTC_ns'].astype('datetime64[ns]') if 'UTC_ns' in display else display['UTC'])
    br=pd.Series(display['Br_nT'],index=utc)
    bt=pd.Series(display['Btheta_nT'],index=utc)
    bf=pd.Series(display['Bphi_nT'],index=utc)
    assert len(br)==7201
    results=[]
    for label in ['Bperp','Bparallel']:
        amp=display[label+'_amplitude_nT']
        freqs=display['frequency_Hz']
        sec=display['time_seconds']
        fion=display['fci_Hz']
        freq=csv['frequency_Hz'].to_numpy()
        power=csv[label+'_PSD_nT2_per_Hz'].to_numpy()
        assert np.all(np.diff(freq)>0) and np.all(np.isfinite(power)) and np.all(power>=0)
        assert np.allclose(freqs,freq,rtol=1e-14,atol=0)
        results.append((amp,freqs,sec,utc,fion,freq,power))
    plt.rcParams.update({
        "font.family":"sans-serif","font.sans-serif":["DejaVu Sans"],
        "font.size":FONT,"axes.labelsize":FONT,"axes.titlesize":TITLE_FONT,
        "xtick.labelsize":FONT,"ytick.labelsize":FONT,"legend.fontsize":FONT,
        "pdf.fonttype":42,"svg.fonttype":"none",
        "axes.linewidth":.8,
    })
    width,height=178.,354.
    fig=plt.figure(figsize=(width/25.4,height/25.4))
    # The five original vertical panels are retained. Repeated plot-area gaps
    # are 32 mm to accommodate the original seconds + UTC double time axes.
    def make_axis(bottom,panel_height):
        return fig.add_axes([21/width,bottom/height,132/width,panel_height/height])
    axes=[make_axis(299,35),make_axis(224,43),make_axis(149,43),
          make_axis(82,35),make_axis(15,35)]
    ids=list("abcde")
    for ax,letter in zip(axes,ids):
        ax.annotate("("+letter+")",xy=(0,1),xycoords="axes fraction",
                    xytext=(-49,32),textcoords="offset points",ha="left",va="bottom",
                    fontsize=PANEL_FONT,fontweight="bold")
    a=axes[0]
    a.plot(br.index,br.values,lw=.9,label=r"$B_r$")
    a.plot(bt.index,bt.values,lw=.9,label=r"$B_\theta$")
    a.plot(bf.index,bf.values,lw=.9,label=r"$B_\phi$")
    a.set_ylabel("B (nT)")
    a.set_xlabel("UTC")
    a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    a.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0,15,30,45]))
    a.grid(True,alpha=.3)
    a.legend(loc="lower left",bbox_to_anchor=(0,1.02),fontsize=FONT,ncol=3,frameon=False,
             columnspacing=1.,handlelength=1.5,borderaxespad=0)
    a.annotate("2017-05-10   08:00:00–10:00:00 UTC",xy=(.5,1),xycoords="axes fraction",
               xytext=(0,32),textcoords="offset points",ha="center",va="bottom",fontsize=TITLE_FONT)
    colorbars=[];twins=[]
    for ax,tup,label in zip(axes[1:3],results,[r"$B_\perp$",r"$B_{||}$"]):
        amp,freqs,sec,utc,fion,_,_=tup
        vmin,vmax=np.percentile(amp,[5,95])
        # The original map and percentile limits are explicitly preserved.
        img=ax.pcolormesh(sec,freqs,amp,shading="auto",cmap="jet",
                          vmin=vmin,vmax=vmax,rasterized=True)
        ax.set_yscale("log");ax.set_ylim(freqs[0],freqs[-1]);ax.set_xlim(sec[0],sec[-1])
        ax.yaxis.set_major_locator(LogLocator(numticks=6))
        ax.yaxis.set_major_formatter(FuncFormatter(sci_formatter))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_ylabel("Frequency (Hz)");ax.set_xlabel("Time (s)")
        # Apply the existing visible frequency bounds to drawing geometry only.
        # Off-axis fci points were already hidden by axes clipping in the original.
        visible_fion=np.ma.masked_outside(fion,freqs[0],freqs[-1])
        ax.plot(sec,visible_fion,color="w",lw=1.2,alpha=.9,label=r"$f_{ci}$ (O$^+$)")
        ax.legend(loc="upper right",facecolor="#333333",edgecolor="white",
                  labelcolor="white",fontsize=FONT)
        box=ax.get_position()
        cbax=fig.add_axes([157/width,box.y0,3/width,box.height])
        cbar=fig.colorbar(img,cax=cbax)
        cbar.set_label("Amplitude (nT)",fontsize=FONT,labelpad=3)
        cbar.ax.tick_params(labelsize=FONT,pad=2)
        colorbars.append(cbax)
        twin=add_utc_axis(ax,utc,label="UTC")
        twin.set_position(ax.get_position())
        twin.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0,15,30,45]))
        twin.tick_params(labelsize=FONT,pad=2)
        twin.set_xlabel("")
        ax.annotate("UTC",xy=(.5,1),xycoords="axes fraction",xytext=(0,17),
                    textcoords="offset points",ha="center",va="bottom",fontsize=FONT)
        twins.append(twin)
        ax.annotate("2017-05-10   "+label+"   08:00:00–10:00:00 UTC",
                    xy=(.5,1),xycoords="axes fraction",xytext=(0,39),
                    textcoords="offset points",ha="center",va="bottom",fontsize=TITLE_FONT)
    for ax,tup,label in zip(axes[3:],results,[r"$B_\perp$",r"$B_{||}$"]):
        _,_,_,_,_,freq,power=tup
        plot_frequency,plot_power=visible_curve(freq,power)
        ax.semilogx(plot_frequency,plot_power,"b-",lw=1.8,label="PSD")
        ax.set_xlim(3.5e-3,50e-3)
        ax.set_ylim(0,25)
        ax.set_yticks([0,5,10,15,20,25])
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel(r"PSD (nT$^2$/Hz)")
        ax.annotate("PSD - 2017-05-10\n08:00:00–10:00:00 UTC   "+label,
                    xy=(.5,1),xycoords="axes fraction",xytext=(0,27),
                    textcoords="offset points",ha="center",va="bottom",fontsize=TITLE_FONT)
        ax.set_xticks([.0035,.01,.05],
                      labels=[r"$3.5\times10^{-3}$",r"$10^{-2}$",r"$5\times10^{-2}$"])
        ax.get_xticklabels()[-1].set_ha("right")
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(True,which="both",ls="--",alpha=.6)
        ax.legend(loc="lower right",bbox_to_anchor=(1,1.02),fontsize=FONT,frameon=False,borderaxespad=0)
    fig.canvas.draw()
    fig.savefig(out/(name+".pdf"),dpi=600)
    fig.savefig(out/(name+".png"),dpi=600)
    fig.savefig(out/(name+".svg"),dpi=600)
    fig.savefig(qa/"figureS3.preview.png",dpi=300)
    plt.close(fig)
    (qa/'figureS3.plot_validation.json').write_text(json.dumps({
        'source':'quiet_psd.csv', 'display_source':'quiet_display_arrays.npz',
        'visible_original_samples_match_csv_exactly':True, 'boundary_interpolation':'log-x line geometry only; not a new PSD estimate',
        'display_frequency_range_Hz':[.0035,.05], 'shared_display_PSD_range_nT2_per_Hz':[0,25], 'full_input_frequency_range_Hz':[float(freq[0]),float(freq[-1])], 'input_frequency_samples_retained':len(freq), 'CWT_amplitudes_unchanged':True, 'raw_field_unchanged':True,
        'size_mm':[178,354], 'base_font_pt':9,'title_font_pt':10,'panel_font_pt':11,
        'user_override':'Retain original jet map, 5th/95th percentile colour limits, layout and labels for panels a/b/c.',
    },indent=2),encoding='utf-8')


if __name__=='__main__':
    main()

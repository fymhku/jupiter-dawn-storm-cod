"""Render calibrated PSD curves with the approved manuscript figure styling.

Only Fig. 4 and Fig. S3(d,e) are revised; Fig. 3 is never written.  The quiet
event helper is an exact copy of the prior plotting workflow except for PSD
input arrays and export/QA bookkeeping.  All data preparation and calibration
are delegated to the accompanying numerical analysis script.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import os

PACKAGE = Path(os.environ.get("JUPITER_RUN_ROOT", str(Path(__file__).resolve().parents[1] / "output"))).resolve()
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.text import Text
import numpy as np
import pandas as pd
from display_curve import visible_curve

OUT = PACKAGE / 'figures'
QA = PACKAGE / 'qa'
FONT = 9.0
PANEL_FONT = 10.0


def highest_interior_peak(frequency, power):
    """Original local maximum rule including midpoint handling of plateaux."""
    candidates = []
    i = 1
    while i < len(power) - 1:
        if power[i] >= power[i-1] and power[i] >= power[i+1]:
            left = i
            right = i
            while right + 1 < len(power) and power[right + 1] == power[left]:
                right += 1
            if left > 0 and right < len(power)-1 and power[left] > power[left-1] and power[right] > power[right+1]:
                candidates.append((left+right)//2)
            i = right+1
        else:
            i += 1
    if not candidates:
        return None
    index = max(candidates, key=lambda j: power[j])
    return {'index':int(index), 'frequency_Hz':float(frequency[index]), 'PSD_nT2_per_Hz':float(power[index])}


def draw_fig4():
    plt.rcdefaults()
    plt.rcParams.update({
        'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':FONT, 'axes.labelsize':FONT, 'axes.titlesize':FONT,
        'xtick.labelsize':FONT,'ytick.labelsize':FONT,'legend.fontsize':FONT,
        'axes.linewidth':.65, 'xtick.major.width':.65, 'ytick.major.width':.65,
        'xtick.major.size':2.5, 'ytick.major.size':2.5,
        'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'dejavusans',
        'path.simplify':False,'figure.dpi':150,'savefig.dpi':600,
    })
    frame = pd.read_csv(PACKAGE / 'data/main_psd.csv',float_precision='round_trip')
    frequency = frame['frequency_Hz'].to_numpy()
    curves = [frame['Bperp_PSD_nT2_per_Hz'].to_numpy(), frame['Bparallel_PSD_nT2_per_Hz'].to_numpy()]
    assert len(frequency)==200 and np.all(np.diff(frequency)>0)
    assert np.all(np.isfinite(curves)) and np.all(np.asarray(curves)>=0)
    # Existing interval's quartiles; kept fixed because the PSD recalibration
    # does not alter the measured magnetic field or its ion frequencies.
    iqr = [0.0038250449843684094,0.006413020737400884]
    width, height = 178.0, 90.0
    fig, axes = plt.subplots(1,2,figsize=(width/25.4,height/25.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=14/height,top=61/height,wspace=.32)
    peaks = {}
    for k,(ax,power,component,key) in enumerate(zip(axes,curves,[r'$B_\perp$',r'$B_\parallel$'],['Bperp','Bparallel'])):
        plot_frequency,plot_power=visible_curve(frequency,power)
        line,=ax.semilogx(plot_frequency,plot_power,'b-',lw=1.2,label='PSD')
        assert np.array_equal(line.get_ydata(),plot_power)
        ax.set_xlim(.0035,.05)
        ax.set_ylim(0,40)
        ax.set_yticks([0,10,20,30,40])
        ax.set_xticks([.0035,.01,.05],labels=[r'$3.5\times10^{-3}$',r'$10^{-2}$',r'$5\times10^{-2}$'])
        ax.get_xticklabels()[-1].set_ha('right')
        ax.xaxis.set_minor_locator(FixedLocator([]))
        ax.set_xlabel('Frequency (Hz)',labelpad=3)
        ax.set_ylabel(r'PSD (nT$^2$/Hz)',labelpad=4)
        ax.grid(True,which='both',ls='--',alpha=.6,lw=.5)
        peak = highest_interior_peak(frequency,power)
        peaks[key] = peak
        if k==0:
            ax.axvspan(*iqr,color='yellow',alpha=.3,label=r'$f_{\mathrm{ci}}$ [25th–75th percentiles]')
            for edge in iqr:
                ax.axvline(edge,color='orange',lw=.7,alpha=.9)
            if peak is not None:
                f,value=peak['frequency_Hz'],peak['PSD_nT2_per_Hz']
                ax.axvline(f,color='red',ls='--',lw=1.0,alpha=.9)
                ax.plot(f,value,'ro',ms=4.5,markeredgecolor='black',markeredgewidth=.6)
                label_y=value-.065*(float(power.max())-float(power.min()))
                ax.text(.0105,label_y,f'{f:.3e} Hz',color='red',fontsize=FONT,va='top')
        ax.set_title(f'PSD – 2017-05-16\n07:52–10:03 UTC (2 h 11 min)\n{component}',fontsize=FONT,pad=34)
        ax.legend(loc='lower right',bbox_to_anchor=(1,1.015),fontsize=FONT,
                  handlelength=1.2,handletextpad=.4,borderpad=.2,borderaxespad=0)
    for obj in fig.findobj(match=Text):
        obj.set_fontsize(FONT)
    for ax,letter in zip(axes,'ab'):
        ax.annotate(f'({letter})',(0,1),xycoords='axes fraction',xytext=(-35,50),
                    textcoords='offset points',ha='left',va='bottom',fontsize=PANEL_FONT,
                    fontweight='bold',annotation_clip=False)
    fig.canvas.draw()
    fig.savefig(OUT/'Fig4_PSD.pdf',dpi=600)
    fig.savefig(OUT/'Fig4_PSD.svg',dpi=600)
    fig.savefig(OUT/'Fig4_PSD.png',dpi=600)
    fig.savefig(QA/'figure4.preview.png',dpi=300)
    plt.close(fig)
    (QA/'figure4.plot_validation.json').write_text(json.dumps({
        'source':'main_psd.csv','visible_original_samples_match_csv_exactly':True, 'boundary_interpolation':'log-x line geometry only; not a new PSD estimate',
        'size_mm':[width,height],'fci_IQR_Hz':iqr,'highest_interior_peaks':peaks,
        'display_frequency_range_Hz':[.0035,.05], 'shared_display_PSD_range_nT2_per_Hz':[0,40], 'full_input_frequency_range_Hz':[float(frequency[0]),float(frequency[-1])], 'input_frequency_samples_retained':len(frequency), 'Fig3_untouched':True,'PSD_units':'nT^2/Hz','base_font_pt':FONT,
    },indent=2),encoding='utf-8')
    print(json.dumps(peaks,indent=2))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--figure',choices=['main','quiet','both'],default='both')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    QA.mkdir(parents=True,exist_ok=True)
    if args.figure in ['main','both']:
        draw_fig4()
    if args.figure in ['quiet','both']:
        # Separate process preserves the original quiet-figure rcParams.
        subprocess.run([sys.executable,'-X','utf8',str(Path(__file__).with_name('plot_quiet_psd.py'))],check=True)


if __name__=='__main__':
    main()

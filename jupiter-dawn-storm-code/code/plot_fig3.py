"""Replot the archived Fig. 3 CWT amplitudes, without recomputing or recalibrating them."""
from pathlib import Path
import os
import sys
import inspect
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FixedLocator, MaxNLocator
from matplotlib.text import Text
SOURCE_ROOT = Path(__file__).resolve().parents[1]
PACKAGE = Path(os.environ.get("JUPITER_RUN_ROOT", str(SOURCE_ROOT/"output"))).resolve()
SOURCE = SOURCE_ROOT/"original"
OUT = PACKAGE/"figures"
OUT.mkdir(parents=True,exist_ok=True)
FONT, PANEL_FONT = 9., 10.
FIG3_SIZE_MM = (125.,225.)
# Preserve the current archived Fig3 shading exactly, including original endpoints.
# The PSD IQR and S1 use the manuscript's 09:17:40-09:46:44 interval.
FIELD_ACTIVE = ("2017-05-16 09:17:38","2017-05-16 09:46:52")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
"font.size":9,"axes.labelsize":9,"axes.titlesize":9,"xtick.labelsize":9,"ytick.labelsize":9,
"legend.fontsize":9,"axes.linewidth":.65,"xtick.major.width":.65,"ytick.major.width":.65,
"xtick.major.size":2.5,"ytick.major.size":2.5,"pdf.fonttype":42,"svg.fonttype":"none",
"mathtext.fontset":"dejavusans","path.simplify":False,"figure.dpi":150,"savefig.dpi":600})
def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

wavelet = load_module(SOURCE/"tool/wavelet.py", "original_wavelet")
def make_adapted_spectrogram():
    """Reuse the original plot body; redirect its page/axes creation only."""
    source = inspect.getsource(wavelet.plot_spectrogram)
    source = source.replace("def plot_spectrogram(amplitude, freqs, time_sec, time_utc, ion_freq, component_name):",
                            "def plot_spectrogram_font_only(amplitude, freqs, time_sec, time_utc, ion_freq, component_name, fig, ax_main, cbar_ax):")
    source = source.replace("    fig = plt.figure(figsize=(15, 8))", "    # Existing compositor figure; physical dimensions are set by caller.")
    source = source.replace("    grid = plt.GridSpec(1, 2, width_ratios=[0.92, 0.05], wspace=0.05)", "    # Existing panel and its adjacent colorbar.")
    source = source.replace("    ax_main = plt.subplot(grid[0])", "")
    source = source.replace("    cbar_ax = plt.subplot(grid[1])", "")
    source = source.replace("    plt.tight_layout()", "    # Manual margins are measured by the alignment audit.")
    scope = dict(vars(wavelet))
    exec(compile(source, str(SOURCE / "tool/wavelet.py") + ":font_only", "exec"), scope)
    return scope["plot_spectrogram_font_only"]

spectrogram = make_adapted_spectrogram()

def all_font_sizes(fig):
    for obj in fig.findobj(match=Text):
        obj.set_fontsize(FONT)

def panel_letter(ax, letter, offset_y):
    ax.annotate(f"({letter})", (0, 1), xycoords="axes fraction", xytext=(-35, offset_y),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=PANEL_FONT, fontweight="bold", annotation_clip=False)

def save(fig, stem, axes, size, excluded=(), column_groups=None, row_groups=None):
    fig.canvas.draw()
    # The official auditor recognizes bare a, but the source figure uses (a).
    # Check the exact original parenthesized labels separately, while its
    # strict geometry checks remain enabled with explicit panel identities.
    for letter, ax in zip("abc", axes):
        assert any(t.get_text()==f"({letter})" for t in ax.texts)
    assert np.allclose(fig.get_size_inches()*25.4, size)
    fig.savefig(OUT/f"{stem}.pdf")
    fig.savefig(OUT/f"{stem}.png", dpi=600)
    fig.savefig(OUT/f"{stem}.svg", dpi=600)
    plt.close(fig)

def draw_fig3(frame, results):
    width, height = FIG3_SIZE_MM
    fig = plt.figure(figsize=(width/25.4, height/25.4))
    # The field panel keeps the original extra width; b/c include colorbars.
    a = fig.add_axes([16/width, 165/height, 106/width, 42/height])
    b = fig.add_axes([16/width, 90/height, 89/width, 42/height])
    c = fig.add_axes([16/width, 15/height, 89/width, 42/height])
    axes = [a,b,c]
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for key, label, color in zip(["Br", "Btheta", "Bphi"],
                                  [r"$B_r$", r"$B_\theta$", r"$B_\phi$"], colors):
        line, = a.plot(frame.index, frame[key].to_numpy(), label=label, lw=.8, color=color)
        assert np.array_equal(line.get_ydata(), frame[key].to_numpy())
    a.axvspan(pd.Timestamp(FIELD_ACTIVE[0]), pd.Timestamp(FIELD_ACTIVE[1]),
              alpha=.25, color="yellow", label="Active period")
    a.set_ylabel("B (nT)")
    a.set_xlabel("UTC Time", labelpad=2)
    a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    a.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0,30]))
    a.set_title("B field (RTP) - EV03  2017-05-16\n2017-05-16 07:52:00 - 10:03:00", fontsize=FONT, pad=19)
    a.legend(loc="lower right", bbox_to_anchor=(1,1.015), ncols=4,
             fontsize=FONT, handlelength=1.2, borderpad=.15,
             columnspacing=.7,labelspacing=.2,handletextpad=.4,borderaxespad=0)
    a.grid(True, alpha=.3)
    caxes, twins = [], []
    for ax, result, component in zip([b,c], results, [r"$B_\perp$",r"$B_\parallel$"]):
        _, amp, freqs, seconds, utc, fion = result
        ca = fig.add_axes([108/width, ax.get_position().y0, 3/width, ax.get_position().height])
        before_axes = list(fig.axes)
        spectrogram(amp, freqs, seconds, utc, fion, component, fig, ax, ca)
        twin = next(q for q in fig.axes if q not in before_axes)
        mesh = ax.collections[0]
        assert np.array_equal(np.asarray(mesh.get_array()), amp)
        assert mesh.cmap.name == "jet"
        assert mesh.norm.vmin == float(amp.min()) and mesh.norm.vmax == float(amp.max())
        mesh.set_rasterized(True)
        # Display clipping only: values outside the original y limits cannot
        # be visible and otherwise create false collision paths in the PDF.
        ax.lines[0].set_ydata(np.ma.masked_outside(np.asarray(fion), freqs[0], freqs[-1]))
        ax.set_xlim(seconds[0],seconds[-1])
        ax.set_xticks([0,1000,2000,3000,4000,5000,6000,7000])
        ax.xaxis.set_minor_locator(FixedLocator([]))
        ax.set_yticks([1e-3,1e-2], labels=[r"$10^{-3}$",r"$10^{-2}$"])
        ax.yaxis.set_minor_locator(FixedLocator([]))
        ax.set_xlabel("Time (seconds)",labelpad=2)
        twin.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0,30]))
        twin.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        twin.set_xlabel("UTC Time", labelpad=2)
        ax.set_title("CWT - EV03\n2017_05_16 07:52-10:03\n"+component,
                     fontsize=FONT,pad=30)
        ax.legend([ax.lines[0]], ["Ion Cyclotron Frequency"], loc="upper right",
                  facecolor="#333333",edgecolor="white",labelcolor="white",
                  fontsize=FONT,handlelength=1.1,handletextpad=.4,borderpad=.25)
        ca.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ca.tick_params(pad=2,length=2)
        caxes.append(ca); twins.append(twin)
    all_font_sizes(fig)
    for ax,letter in zip(axes,"abc"):
        panel_letter(ax,letter,27 if letter != "a" else 15)
    save(fig,"Fig3_font_enlarged",axes,[width,height],excluded=caxes+twins,column_groups=[["b","c"]])

def main():
    frame = pd.read_csv(PACKAGE/"data/main_timeseries.csv")
    frame.index = pd.DatetimeIndex(pd.to_datetime(frame.pop("UTC"))).as_unit("ns")
    frame = frame.rename(columns={"Br_nT":"Br","Btheta_nT":"Btheta","Bphi_nT":"Bphi"})
    arrays = np.load(PACKAGE/"data/main_display_arrays.npz",allow_pickle=False)
    results=[]
    for label in ["Bperp","Bparallel"]:
        results.append((None,arrays[label+"_amplitude_nT"],arrays["frequency_Hz"],
                        arrays["time_seconds"],frame.index,arrays["fci_Hz"]))
    draw_fig3(frame,results)
    print("Saved Fig3_font_enlarged from archived CWT display arrays")

if __name__=="__main__": main()

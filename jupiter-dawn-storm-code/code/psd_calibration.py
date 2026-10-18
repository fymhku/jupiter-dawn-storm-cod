"""Reproducible wavelet-based one-sided PSD for the 2017 Juno events.

This authorized revision changes only the 1-D spectral analysis: the original
20-min field decomposition is followed by subtraction of each complete-window
component mean. For each cmor1.5-1.0 scale, exclude ceil(sqrt(1.5)*scale) samples
at each end and calibrate the remaining mean raw coefficient power by the
actual discrete impulse-response energy. The existing CWT display figures are
not read as input or modified.

Run with a Python environment containing numpy, pandas, scipy, PyWavelets:
    python psd_calibration.py --mag-dir /path/to/JSS_1s_data
PyWavelets must support cwt(..., precision=...), e.g. distribution >=1.9.0.
The exact installed distribution metadata, function signature and source hash
are recorded because the local distribution and module version strings differ.

Writes to the selected output directory; the bundled reference inputs are not changed.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import inspect
import json
import platform
import os
import sys
import time

HERE = Path(__file__).resolve().parent
INPUT_ROOT = HERE.parent
PACKAGE = Path(os.environ.get("JUPITER_RUN_ROOT", str(INPUT_ROOT / "output"))).resolve()
SAVE_COEFFICIENTS = False

import numpy as np
import pandas as pd
import scipy
from scipy.signal import welch, fftconvolve
import pywt

WAVELET = "cmor1.5-1.0"
BETA = 1.5
DT = 1.0
FINAL_PRECISION = 20
PRECISIONS = (12, 16, 18, 20)
EVENTS = {
    "main": {"date": "2017_05_16", "start": "2017-05-16 07:52:00",
             "stop": "2017-05-16 10:03:00", "fmin": .003, "N": 7861},
    "quiet": {"date": "2017_05_10", "start": "2017-05-10 08:00:00",
              "stop": "2017-05-10 10:00:00", "fmin": .001, "N": 7201},
}
COMPONENTS = ("Bperp", "Bparallel")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")


def log(value):
    print(value, flush=True)


def decompose(frame):
    """Original component-wise 1200-s block projection; no scalar mean removal."""
    br, bt, bf = [frame[k].to_numpy() for k in ("Br_nT", "Btheta_nT", "Bphi_nT")]
    par, perp = np.empty(len(frame)), np.empty(len(frame))
    for first in range(0, len(frame), 1200):
        sl = slice(first, min(first+1200, len(frame)))
        mean = np.array([br[sl].mean(), bt[sl].mean(), bf[sl].mean()])
        norm = np.linalg.norm(mean)
        direction = mean/norm if norm > 1e-12 else np.array([0.,0.,1.])
        bp = br[sl]*direction[0]+bt[sl]*direction[1]+bf[sl]*direction[2]
        total = np.sqrt(br[sl]**2+bt[sl]**2+bf[sl]**2)
        par[sl] = bp
        perp[sl] = np.sqrt(np.maximum(total**2-bp**2, 0.))
    return par, perp


def load_event(mag_dir, event):
    spec = EVENTS[event]
    path = Path(mag_dir)/f"JUNO_MAG_JSS_1s_{spec['date']}.txt"
    if path.exists():
        frame = pd.read_csv(path, sep=r"\s+", skiprows=1,
                            names=["UTC","Br_nT","Btheta_nT","Bphi_nT"])
        frame["UTC"] = pd.to_datetime(frame.UTC, format="%Y-%m-%d/%H:%M:%S")
    else:
        # A distributed revision package can be reproduced from its identical
        # exported field window without access to the authors' drive layout.
        path = INPUT_ROOT/"data"/f"{event}_timeseries.csv"
        frame = pd.read_csv(path,usecols=["UTC","Br_nT","Btheta_nT","Bphi_nT"])
        frame["UTC"] = pd.to_datetime(frame.UTC)
    frame = frame[(frame.UTC >= spec["start"]) & (frame.UTC <= spec["stop"])].reset_index(drop=True)
    assert len(frame) == spec["N"]
    assert np.all(np.diff(frame.UTC.to_numpy(dtype="datetime64[s]").astype("int64")) == 1)
    assert np.isfinite(frame.iloc[:,1:].to_numpy()).all()
    par, perp = decompose(frame)
    frame["Bparallel_nT"], frame["Bperp_nT"] = par, perp
    frame["Btotal_nT"] = np.sqrt(frame.Br_nT**2+frame.Btheta_nT**2+frame.Bphi_nT**2)
    frame["fci_Hz"] = .0152*frame.Btotal_nT/16.
    for c in COMPONENTS:
        frame[c+"_demeaned_nT"] = frame[c+"_nT"]-frame[c+"_nT"].mean()
    return frame, path


def coefficients(x, scales, precision, method="conv"):
    return pywt.cwt(np.asarray(x, dtype=float), scales, WAVELET,
                    sampling_period=DT, method=method, precision=precision)[0]


def actual_kernels(scales, precision):
    """Measure the CWT impulse response, including its exact even/odd alignment."""
    half = int(np.ceil(8*np.max(scales)))+8
    impulse = np.zeros(2*half+1)
    impulse[half] = 1.
    h = coefficients(impulse, scales, precision, method="fft")
    lags = np.arange(len(impulse))-half
    return h, lags


def coi_margins(scales, n):
    margin = np.ceil(np.sqrt(BETA)*np.asarray(scales)).astype(int)
    valid = n-2*margin
    if np.any(valid <= 0):
        raise ValueError("One or more scales have no samples outside the COI")
    return margin, valid


def interior_power(c, margin):
    return np.array([np.mean(np.abs(row[...,m:-m])**2, axis=-1)
                     for row,m in zip(c, margin)])


def spectrum(x, scales, precision):
    """Return PSD, uncalibrated interior power, actual filter energy and CWT."""
    c = coefficients(x, scales, precision)
    h, lags = actual_kernels(scales, precision)
    q = np.sum(np.abs(h)**2, axis=1)
    margin, valid = coi_margins(scales, len(x))
    g = interior_power(c, margin)
    return 2*DT*g/q, g, q, c, h, lags, margin, valid


def peak_report(f, p):
    """Original plateau-aware interior local-maximum rule; no interpolation."""
    n = len(p)
    candidate = np.zeros(n, dtype=bool)
    candidate[1:-1] = (p[1:-1] >= p[:-2]) & (p[1:-1] >= p[2:])
    found = []
    i = 1
    while i < n-1:
        if not candidate[i]:
            i += 1
            continue
        first = i
        while i < n-1 and candidate[i]:
            i += 1
        last = i-1
        if p[first] > p[first-1] and p[last] > p[last+1]:
            k = (first+last)//2
            found.append({"index":int(k),"frequency_Hz":float(f[k]),"value":float(p[k])})
    ranked = sorted(found,key=lambda row:row["value"],reverse=True)
    k = int(np.argmax(p))
    return {"all_local_maxima_in_frequency_order":found,
            "highest_interior_local_maximum":ranked[0] if ranked else None,
            "all_local_maxima_above_3mHz":[row for row in found if row["frequency_Hz"]>.003],
            "global_maximum":{"index":k,"frequency_Hz":float(f[k]),"value":float(p[k]),
                              "is_endpoint":k in (0,n-1)}}


def difference(a,b):
    a,b = np.asarray(a),np.asarray(b)
    return {"max_absolute":float(np.max(np.abs(a-b))),
            "relative_L2":float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-300)),
            "max_pointwise_relative":float(np.max(np.abs(b-a)/np.maximum(np.abs(a),1e-300)))}


def old_cache(event, component, frequencies):
    """Optional historical comparison; never a prerequisite for the new PSD.

    Prefer the spectra shipped with the package. The former workspace cache is
    only an optional fallback while authoring; an absent comparison returns
    (None, None), and the new calculation and its independent tests continue.
    """
    archived = INPUT_ROOT/"data"/f"{event}_{component}_original_raw_spectrum.csv"
    packaged_table = INPUT_ROOT/"data"/f"{event}_psd.csv"
    if archived.exists():
        path=archived
        d=pd.read_csv(path,float_precision="round_trip")
        f,p=d.frequency_Hz.to_numpy(),d.mean_wavelet_power_nT2.to_numpy()
    elif packaged_table.exists():
        d=pd.read_csv(packaged_table,float_precision="round_trip")
        column=component+"_old_raw_mean_power_nT2"
        if column not in d or not np.isfinite(d[column]).all():
            return None,None
        path=packaged_table
        f,p=d.frequency_Hz.to_numpy(),d[column].to_numpy()
    else:
        return None,None
    if not np.array_equal(f,frequencies) or not np.isfinite(p).all():
        return None,None
    if not archived.exists():
        pd.DataFrame({"frequency_Hz":f,"mean_wavelet_power_nT2":p}).to_csv(archived,index=False)
    return p.copy(),path


def filter_checks(h,lags,f,n,margin,q):
    nfft = 2**18
    ff = np.fft.fftfreq(nfft,d=DT)
    rows=[]
    for kernel,fc,m,qj in zip(h,f,margin,q):
        hp=np.abs(np.fft.fft(kernel,n=nfft))**2
        total=hp.sum()
        central=(ff>.5*fc)&(ff<1.5*fc)
        center=len(kernel)//2
        # Exact diagonal of the finite CWT covariance for iid input after
        # global demeaning: sum_k |A_nk|^2 - |sum_k A_nk|^2/N.
        row_energy=fftconvolve(np.ones(n),np.abs(kernel)**2)[center:center+n]
        row_sum=fftconvolve(np.ones(n),kernel)[center:center+n]
        effective=np.mean((row_energy-np.abs(row_sum)**2/n)[m:-m])
        rows.append({"frequency_Hz":float(fc),"Q":float(qj),
                     "frequency_response_peak_Hz":float(ff[np.argmax(hp)]),
                     "mainband_energy_fraction":float(hp[central].sum()/total),
                     "negative_frequency_energy_fraction":float(hp[ff<0].sum()/total),
                     "effective_white_noise_gain_after_demean_COI":float(effective),
                     "white_noise_expected_PSD_over_ideal":float(effective/qj)})
    return rows


def verify_kernel_alignment():
    rng=np.random.default_rng(20261016)
    scales=np.array([2.,2.25,3.,19.25,20.,20.0625,237.421534])
    h,lags=actual_kernels(scales,FINAL_PRECISION)
    center=int(np.flatnonzero(lags==0)[0])
    x=rng.normal(size=257)
    reference=coefficients(x,scales,FINAL_PRECISION)
    manual=np.array([np.convolve(x,k)[center:center+len(x)] for k in h])
    fft=coefficients(x,scales,FINAL_PRECISION,method="fft")
    out={"scales":scales.tolist(),"pulse_reconstruction_vs_pywt_conv":difference(reference,manual),
         "pywt_fft_vs_conv":difference(reference,fft)}
    out["pass"]=out["pulse_reconstruction_vs_pywt_conv"]["max_absolute"]<1e-10 and out["pywt_fft_vs_conv"]["max_absolute"]<1e-10
    assert out["pass"]
    return out


def synthetic_tests():
    rng=np.random.default_rng(2026101602)
    f=np.array([.001,.003,.004211917906742496,.01,.05])
    s=1/f
    n,replicates=8192,128
    h,lags=actual_kernels(s,FINAL_PRECISION)
    q=np.sum(np.abs(h)**2,axis=1)
    margin,valid=coi_margins(s,n)
    x=rng.normal(size=(replicates,n))
    x-=x.mean(axis=-1,keepdims=True)
    c=coefficients(x,s,FINAL_PRECISION,method="fft")
    estimates=2*DT*interior_power(c,margin)/q[:,None]
    means=estimates.mean(axis=1)
    se=estimates.std(axis=1,ddof=1)/np.sqrt(replicates)
    diag=filter_checks(h,lags,f,n,margin,q)
    expected=2*np.array([row["white_noise_expected_PSD_over_ideal"] for row in diag])
    z=(means-expected)/se
    noise={"seed":2026101602,"N":n,"replicates":replicates,"input_population_variance":1.,
           "frequency_Hz":f.tolist(),"ideal_one_sided_PSD":2.,
           "expected_after_exact_finite_demean_COI":expected.tolist(),
           "mean_estimate":means.tolist(),"standard_error":se.tolist(),"z_scores":z.tolist(),
           "pass":bool(np.all(np.abs(z)<4.))}
    assert noise["pass"]
    # A pure tone is broadened by the Morlet response. A small peak offset from
    # its input frequency is expected for this constant-Q spectral smoother.
    f_tone=.004211917906742496
    fgrid=np.logspace(np.log10(.003),np.log10(.05),200)
    phases=np.linspace(0,2*np.pi,16,endpoint=False)
    tone=np.sqrt(2.)*np.cos(2*np.pi*f_tone*np.arange(n)[None,:]+phases[:,None])
    tone-=tone.mean(axis=-1,keepdims=True)
    st=1/fgrid
    ht,lt=actual_kernels(st,FINAL_PRECISION)
    qt=np.sum(abs(ht)**2,axis=1)
    mt,vt=coi_margins(st,n)
    ct=coefficients(tone,st,FINAL_PRECISION,method="fft")
    measured=2*DT*interior_power(ct,mt).mean(axis=1)/qt
    # Prediction from the actual finite kernels for phase-averaged demeaned
    # cosines, independently of pywt's transform of the test signal.
    positive=np.exp(2j*np.pi*f_tone*np.arange(n))
    positive-=positive.mean()
    predicted=[]
    for kernel,m in zip(ht,mt):
        center=len(kernel)//2
        a=fftconvolve(positive,kernel)[center:center+n]
        b=fftconvolve(positive.conj(),kernel)[center:center+n]
        predicted.append(.5*np.mean((abs(a)**2+abs(b)**2)[m:-m]))
    predicted=2*DT*np.array(predicted)/qt
    delta=difference(predicted,measured)
    peak=peak_report(fgrid,measured)["highest_interior_local_maximum"]
    expected_peak=peak_report(fgrid,predicted)["highest_interior_local_maximum"]
    sine={"input_frequency_Hz":f_tone,"input_amplitude_nT":float(np.sqrt(2.)),
          "phase_count":len(phases),"N":n,"measured_peak":peak,"predicted_peak":expected_peak,
          "peak_relative_offset_from_input":float(peak["frequency_Hz"]/f_tone-1),
          "PSD_vs_independent_kernel_prediction":delta,
          "pass":bool(delta["relative_L2"]<1e-10 and peak["index"]==expected_peak["index"]),
          "interpretation":"Frequency-dependent Morlet smoothing can shift a sinusoid's PSD peak slightly; test against filter prediction, not an invented requirement of exact input-frequency recovery."}
    assert sine["pass"]
    return {"white_noise":noise,"sinusoid":sine}


def run_event(event,mag_dir,data_dir,qa_dir,results,validation):
    frame,source=load_event(mag_dir,event)
    f=np.logspace(np.log10(EVENTS[event]["fmin"]),np.log10(.05),200)
    scales=pywt.central_frequency(WAVELET)/(f*DT)
    n=len(frame)
    margin,valid=coi_margins(scales,n)
    table={"frequency_Hz":f,"scale_samples":scales,"COI_margin_samples":margin,
           "N_valid":valid,"valid_start_UTC":frame.UTC.iloc[margin].astype(str).to_numpy(),
           "valid_stop_UTC":frame.UTC.iloc[n-margin-1].astype(str).to_numpy()}
    npz={"frequency_Hz":f,"scale_samples":scales,"COI_margin_samples":margin,"N_valid":valid}
    event_result={"date":EVENTS[event]["date"],"analysis_window_UTC":[EVENTS[event]["start"],EVENTS[event]["stop"]],
                  "N":n,"input_path":str(source),"input_sha256":digest(source),
                  "PSD_frequency_range_Hz":[float(f[0]),float(f[-1])],"n_frequencies":len(f),
                  "N_valid_min_max":[int(valid.min()),int(valid.max())],
                  "COI_margin_samples_min_max":[int(margin.min()),int(margin.max())],"components":{}}
    if event == "main":
        event_result["fci_IQR_Hz"] = [.0038250449843684094,.006413020737400884]
        event_result["fci_IQR_provenance"] = "Preserved manuscript IQR for 09:17:40-09:46:44 UT; no PSD-dependent change"
    event_validation={"components":{}}
    frame.to_csv(data_dir/f"{event}_timeseries.csv",index=False)
    final_h=None
    final_q=None
    # Kernels are common to components, computed once at each precision.
    kernel_by_precision={}
    for precision in PRECISIONS:
        h,lags=actual_kernels(scales,precision)
        kernel_by_precision[precision]=(h,lags,np.sum(np.abs(h)**2,axis=1))
    for component in COMPONENTS:
        x_original=frame[component+"_nT"].to_numpy()
        x=frame[component+"_demeaned_nT"].to_numpy()
        old,path=old_cache(event,component,f)
        old_available=old is not None
        old_values=old if old_available else np.full(len(f),np.nan)
        table[component+"_old_raw_mean_power_nT2"]=old_values
        npz[component+"_old_raw_mean_power_nT2"]=old_values
        # Exact reproduction of the source algorithm is checked separately.
        if old_available:
            original_coeff=coefficients(x_original,scales,12)
            original_recomputed=np.mean(abs(original_coeff)**2,axis=1)
            old_comparison=difference(old,original_recomputed)
            old_comparison["status"]="checked"
            old_comparison["pass"]=bool(np.allclose(old,original_recomputed,rtol=1e-10,atol=1e-10))
            old_comparison["note"]="Historical reproduction is reported independently; a different library's legacy behavior does not prevent the calibrated calculation or its tests."
        else:
            old_comparison={"status":"skipped","reason":"No usable historical spectrum supplied; it is not required for the calibrated PSD."}
        precision_powers={}
        for precision in PRECISIONS:
            c=coefficients(x,scales,precision)
            h,lags,q=kernel_by_precision[precision]
            raw=interior_power(c,margin)
            psd=2*DT*raw/q
            precision_powers[precision]=psd
            npz[component+f"_PSD_precision{precision}"]=psd
            table[component+f"_PSD_precision{precision}_nT2_per_Hz"]=psd
            if precision==FINAL_PRECISION:
                final_h,final_q=h,q
                table[component+"_PSD_nT2_per_Hz"]=psd
                table[component+"_valid_raw_mean_power_nT2"]=raw
                npz[component+"_PSD_nT2_per_Hz"]=psd
                npz[component+"_valid_raw_mean_power_nT2"]=raw
                # Preserve complete coefficients for exact reproducibility;
                # the mask is defined by the per-scale margin arrays.
                if SAVE_COEFFICIENTS:
                    np.savez_compressed(data_dir/f"{event}_{component}_final_coefficients.npz",
                                        frequency_Hz=f,scale_samples=scales,coefficients=c,
                                        COI_margin_samples=margin,N_valid=valid)
        psd=precision_powers[FINAL_PRECISION]
        new_peaks=peak_report(f,psd)
        old_peaks=peak_report(f,old) if old_available else None
        old_top=old_peaks["highest_interior_local_maximum"] if old_available else None
        new_top=new_peaks["highest_interior_local_maximum"]
        event_result["components"][component]={
            "subtracted_full_window_mean_nT":float(x_original.mean()),
            "old_raw_spectrum_peaks":old_peaks,"new_PSD_peaks":new_peaks,
            "old_cache_path":str(path) if old_available else None,
            "old_cache_sha256":digest(path) if old_available else None,
            "highest_local_peak_frequency_change_Hz":None if old_top is None or new_top is None else new_top["frequency_Hz"]-old_top["frequency_Hz"],
            "new_PSD_value_at_old_highest_local_peak":None if old_top is None else float(psd[old_top["index"]]),
        }
        fw,pw=welch(x,fs=1.,window="hann",nperseg=2048,noverlap=1024,nfft=2048,
                    detrend=False,return_onesided=True,scaling="density",average="mean")
        selected=(fw>=f.min())&(fw<=f.max())
        pd.DataFrame({"frequency_Hz":fw,"PSD_nT2_per_Hz":pw}).to_csv(data_dir/f"{event}_{component}_Welch.csv",index=False)
        event_validation["components"][component]={
            "original_cache_reproduction":old_comparison,
            "precision12_vs_final":difference(precision_powers[12],psd),
            "precision16_vs_final":difference(precision_powers[16],psd),
            "precision18_vs_final":difference(precision_powers[18],psd),
            "precision_peaks":{str(k):peak_report(f,v) for k,v in precision_powers.items()},
            "Welch_peaks_in_same_frequency_band":peak_report(fw[selected],pw[selected]),
            "Welch_band_integral_nT2":float(np.trapezoid(pw[selected],fw[selected])),
            "wavelet_band_integral_nT2":float(np.trapezoid(psd,f)),
            "Welch_preprocessing":"Same full-record mean-subtracted component. Hann segments 2048 samples; 50% overlap; no per-segment detrending. Wavelet uses frequency-dependent COI masks, so effective time/frequency smoothing differs.",
        }
    table["Q_filter_energy"]=final_q
    table["calibration_factor_seconds"]=2*DT/final_q
    npz["Q_filter_energy"]=final_q
    npz["calibration_factor_seconds"]=2*DT/final_q
    pd.DataFrame(table).to_csv(data_dir/f"{event}_psd.csv",index=False)
    np.savez_compressed(data_dir/f"{event}_spectra.npz",**npz)
    diag=filter_checks(final_h,lags,f,n,margin,final_q)
    pd.DataFrame(diag).to_csv(qa_dir/f"numerical_{event}_filter_response.csv",index=False)
    event_validation["filter_response"]={
        "minimum_mainband_energy_fraction":min(r["mainband_energy_fraction"] for r in diag),
        "maximum_negative_frequency_energy_fraction":max(r["negative_frequency_energy_fraction"] for r in diag),
        "Q_min_max":[float(final_q.min()),float(final_q.max())],
        "white_noise_gain_after_COI_min_max":[min(r["white_noise_expected_PSD_over_ideal"] for r in diag),max(r["white_noise_expected_PSD_over_ideal"] for r in diag)],
        "remaining_edge_effect":"COI removes the conventional power e^-2 boundary region, not all nonzero wavelet tails. Residual finite-record white-noise gain is reported, not renormalized away."}
    event_validation["precision18_vs_final_max_relative_difference"]=max(event_validation["components"][c]["precision18_vs_final"]["max_pointwise_relative"] for c in COMPONENTS)
    event_validation["precision18_vs_final_below_0_1percent"]=event_validation["precision18_vs_final_max_relative_difference"]<.001
    results["events"][event]=event_result
    validation["events"][event]=event_validation
    dump(data_dir/"numerical_results.json",results)
    dump(qa_dir/"numerical_validation.json",validation)
    log(f"DATA READY {event}: "+json.dumps({c:event_result["components"][c]["new_PSD_peaks"] for c in COMPONENTS}))


def main():
    global SAVE_COEFFICIENTS
    parser=argparse.ArgumentParser()
    parser.add_argument("--mag-dir",type=Path,default=INPUT_ROOT/"data/mag")
    parser.add_argument("--save-coefficients",action="store_true",help="also save the large complex CWT coefficient arrays")
    args=parser.parse_args()
    SAVE_COEFFICIENTS = args.save_coefficients
    if "precision" not in inspect.signature(pywt.cwt).parameters:
        raise RuntimeError("This script requires PyWavelets cwt with an explicit precision argument; install a supported distribution (>=1.9).")
    data_dir,qa_dir=PACKAGE/"data",PACKAGE/"qa"
    data_dir.mkdir(parents=True,exist_ok=True)
    qa_dir.mkdir(parents=True,exist_ok=True)
    started=time.time()
    cwt_file=Path(inspect.getfile(pywt.cwt))
    results={"status":"running","final_precision":FINAL_PRECISION,"wavelet":WAVELET,
             "sampling_interval_seconds":DT,"PSD_units":"nT^2/Hz",
             "PSD_formula":"2*dt*mean_valid(abs(W_s[n])**2)/sum_k(abs(h_s[k])**2)",
             "preprocessing":"Original nonoverlapping 1200-s vector-mean decomposition, then subtract the scalar mean of each complete-window component before CWT. No band-pass filter or taper.",
             "COI_formula":"margin_samples=ceil(sqrt(1.5)*scale); retain indices margin through N-margin-1",
             "COI_convention":"Complex-Morlet power envelope exp(-2*u^2/1.5); excluded boundary extends to e^-2 of its peak power.",
             "scope":"PSD curves only; existing Fig3 and FigS3 CWT amplitude maps are unchanged and use their original preprocessing.",
             "peak_rule":"Highest plateau-aware interior local maximum; endpoints excluded; no interpolation, fitting, or significance threshold.",
             "units_and_limit":"Frequency-smoothed one-sided PSD estimator calibrated by actual filter energy. Its numerical integration over center frequencies is not asserted to reproduce variance exactly.",
             "versions":{"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,
                         "PyWavelets_distribution":importlib.metadata.version("PyWavelets"),"pywt_module_string":pywt.__version__},
             "runtime":{"cwt_signature":str(inspect.signature(pywt.cwt)),"cwt_source_file":str(cwt_file),"cwt_source_SHA256":digest(cwt_file),
                        "script_SHA256":digest(__file__),"version_note":"The local installed PyWavelets distribution metadata and module string differ; both plus source hash are recorded. The script requires the public precision keyword and does not monkeypatch PyWavelets."},
             "events":{}}
    validation={"status":"running","kernel_alignment":verify_kernel_alignment(),"events":{},
                "Welch_parameters":{"fs":1.,"window":"hann","nperseg":2048,"noverlap":1024,"nfft":2048,"detrend":False,"scaling":"density","return_onesided":True,"average":"mean"}}
    dump(data_dir/"numerical_results.json",results)
    for event in EVENTS:
        log(f"Computing {event}, precisions={PRECISIONS}")
        run_event(event,args.mag_dir,data_dir,qa_dir,results,validation)
    log("Validating white noise and sinusoidal filter response")
    validation["synthetic_tests"]=synthetic_tests()
    validation["elapsed_seconds"]=time.time()-started
    results["status"]="complete"
    validation["status"]="complete"
    validation["convergence_pass"]=all(v["precision18_vs_final_below_0_1percent"] for v in validation["events"].values())
    dump(data_dir/"numerical_results.json",results)
    dump(qa_dir/"numerical_validation.json",validation)
    files=[Path(__file__),*sorted(data_dir.glob("*")),*sorted(qa_dir.glob("numerical*"))]
    dump(qa_dir/"numerical_SHA256.json",{("code/"+p.name if p==Path(__file__) else str(p.relative_to(PACKAGE))):digest(p) for p in files if p.is_file() and p.name!="numerical_SHA256.json"})
    log(f"COMPLETE in {time.time()-started:.1f}s; convergence_pass={validation['convergence_pass']}")


if __name__=="__main__":
    main()

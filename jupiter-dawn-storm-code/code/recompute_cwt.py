"""Optional CWT amplitude recomputation; does not overwrite the archived display arrays."""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import pywt
from psd_calibration import decompose
SOURCE_ROOT=Path(__file__).resolve().parents[1]
PACKAGE=Path(os.environ.get("JUPITER_RUN_ROOT",str(SOURCE_ROOT/"output"))).resolve()
out=PACKAGE/"data"
out.mkdir(parents=True,exist_ok=True)
for event in ["main","quiet"]:
    frame=pd.read_csv(SOURCE_ROOT/"data"/f"{event}_timeseries.csv")
    par,perp=decompose(frame)
    frequency=np.logspace(np.log10(.001),np.log10(.05),200)
    scales=pywt.central_frequency("cmor1.5-1.0")/frequency
    result={"frequency_Hz":frequency,"time_seconds":np.arange(len(frame),dtype=float),
            "fci_Hz":.0152*np.sqrt(sum(frame[c].to_numpy()**2 for c in ["Br_nT","Btheta_nT","Bphi_nT"]))/16}
    for label,series in [("Bperp",perp),("Bparallel",par)]:
        coefficients,_=pywt.cwt(series,scales,"cmor1.5-1.0",sampling_period=1.,method="conv",precision=12)
        result[label+"_amplitude_nT"]=np.abs(coefficients)/np.sqrt(scales)[:,None]
    np.savez_compressed(out/f"{event}_display_recomputed.npz",**result)
    print(f"Saved {event} CWT amplitudes (precision=12, no PSD demeaning or COI mask)")

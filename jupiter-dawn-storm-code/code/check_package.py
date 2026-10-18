"""Check release-file integrity and the observation/PSD input contracts."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
results={"checksums":0,"events":{},"license":"not specified by author"}
for line in (ROOT/"SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
    digest,name=line.split("  ",1)
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    results["checksums"]+=1
for event,n in [("main",7861),("quiet",7201)]:
    frame=pd.read_csv(ROOT/"data"/f"{event}_timeseries.csv")
    utc=pd.DatetimeIndex(pd.to_datetime(frame.UTC))
    assert len(frame)==n and np.all(np.diff(utc.to_numpy(dtype="datetime64[s]").astype("int64"))==1)
    assert set(utc.year)=={2017}
    assert np.isfinite(frame[["Br_nT","Btheta_nT","Bphi_nT"]].to_numpy()).all()
    psd=pd.read_csv(ROOT/"data"/f"{event}_psd.csv",float_precision="round_trip")
    assert len(psd)==200 and np.all(np.diff(psd.frequency_Hz)>0)
    assert np.isfinite(psd[["Bperp_PSD_nT2_per_Hz","Bparallel_PSD_nT2_per_Hz"]]).all().all()
    assert (psd[["Bperp_PSD_nT2_per_Hz","Bparallel_PSD_nT2_per_Hz"]]>=0).all().all()
    results["events"][event]={"samples":n,"start_UTC":str(utc[0]),"stop_UTC":str(utc[-1]),"PSD_samples":200}
    arrays=np.load(ROOT/"data"/f"{event}_display_arrays.npz",allow_pickle=False)
    assert arrays["Bperp_amplitude_nT"].shape==(200,n)
    assert arrays["Bparallel_amplitude_nT"].shape==(200,n)
results["status"]="pass"
out=Path(os.environ.get("JUPITER_RUN_ROOT",str(ROOT/"output")))/"qa"
out.mkdir(parents=True,exist_ok=True)
(out/"input_validation.json").write_text(json.dumps(results,indent=2)+"\n",encoding="utf-8")
print(json.dumps(results,indent=2))

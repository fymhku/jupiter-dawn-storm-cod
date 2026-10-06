# Juno ULF-wave analysis

Reproducible analysis code accompanying the paper on **ion-cyclotron-frequency
ultra-low-frequency (ULF) waves in the Jovian magnetosphere** observed by the
Juno magnetometer (MAG).

The package covers the four analyses used in the paper:

| Script | Analysis |
|--------|----------|
| `run_cwt.py`      | Continuous wavelet transform (CWT) spectrograms |
| `run_psd.py`      | Time-averaged wavelet PSD + ion-cyclotron frequency band |
| `run_bandpass.py` | Band-passed perpendicular field and the active period |
| `run_field_plot.py` | Magnetic-field components (JSS-RTP) with the active period |

---

## 1. Data

**Source.** Juno MAG JSS 1-s spherical (RTP) products, publicly available from
the NASA Planetary Data System.  One ASCII file per day.

**Placement.** Put the daily files in one folder and point the code at it:

```
JUNO_MAG_JSS_1s_2017_03_17.txt
JUNO_MAG_JSS_1s_2017_05_16.txt
...
```

Either edit `MAG_DIR` in `config.py` or set the environment variable

```bash
set JUNO_MAG_DIR=D:\path\to\juno_B1s_JSSRTP     # Windows
export JUNO_MAG_DIR=/path/to/juno_B1s_JSSRTP    # Linux / macOS
```

**Format** (header line followed by 1-s records):

```
UT_Time                        Br(nT)          Btheta(nT)          Bphi(nT)
2017-03-17/00:00:00       0.40889433       -1.6055669       -2.5356380
2017-03-17/00:00:01       0.37830267       -1.6117249       -2.4714024
```

* `UT_Time` — `YYYY-MM-DD/HH:MM:SS`
* `Br`, `Btheta`, `Bphi` — spherical components in nT
* cadence — 1 s

---

## 2. Installation

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows
# source .venv/bin/activate       # Linux / macOS
pip install -r requirements.txt
```

Tested with Python 3.12, numpy 1.26, scipy 1.13, pandas 2.x,
matplotlib 3.9, PyWavelets 1.4.

---

## 3. Configuration

Everything used in the paper lives in `config.py`:

* `EVENTS` — the two events, their 2-h analysis windows, their active periods
  and the band-pass bands plotted for each;
* `CWT` — wavelet, number of scales, frequency range, ion mass;
* `PSD` — frequency band, normalisation, ion-cyclotron quantiles;
* `BANDPASS` — filter order, background block length, envelope settings.

### Events used in the paper

| Event | Window (UTC) | Active period (UTC) |
|-------|--------------|---------------------|
| 2017-03-17 | 08:56:30 – 10:56:30 | 10:13:29 – 10:28:13 |
| 2017-05-16 | 07:52:00 – 10:03:00 | 09:17:38 – 09:46:52 |

For 2017-05-16 the wave packet detected in the band-passed perpendicular field
is 09:27:38–09:36:52; it was **widened by ±10 min** so that the standard
[0.25, 0.75] ion-cyclotron interval contains the B_perp spectral peak.

---

## 4. Usage

Run from the project root; figures are written to `output/` (override with the
`JUNO_ULF_OUT` environment variable).

```bash
python run_cwt.py         # -> output/cwt/
python run_psd.py         # -> output/psd/
python run_bandpass.py    # -> output/bandpass/
python run_field_plot.py  # -> output/field/
```

Each script loops over all events in `config.EVENTS` and prints the numbers it
derives (peak frequencies, ion-cyclotron band edges, detected packets).

---

## 5. Methods

### 5.1 Parallel / perpendicular decomposition

`juno_ulf/components.py`.  The background (parallel) direction is the mean field
vector of a **running 20-min block**.  For each sample in the block

```
b0        = <B> / |<B>|
B_par(t)  = B(t) . b0
B_perp(t) = sqrt(|B(t)|^2 - B_par(t)^2)
```

A block mean (rather than one average over the whole window) keeps the
background definition local, which matters when the field rotates.

### 5.2 Continuous wavelet transform

`juno_ulf/cwt.py`.  Complex Morlet wavelet `cmor1.5-1.0`, 200 logarithmically
spaced scales between 3 and 50 mHz, 1-s sampling.  Coefficients are normalised
by `1/sqrt(scale)` before the amplitude is taken.

### 5.3 Ion cyclotron frequency

```
f_ci = e B / (2 pi m_i)  =  0.0152 * B[nT] / (m_i [amu])   Hz
```

The default ion mass is **16 amu (O+)**, i.e. `f_ci = 9.5e-4 * B[nT]` Hz.
Change `CWT["ion_mass_amu"]` for a different ion (1 amu = proton).

### 5.4 PSD

`juno_ulf/psd.py`.  The "PSD" is the time average of the squared wavelet
coefficients, `E(f) = <|W(f,t)|^2>_t`, evaluated in 3–50 mHz.  Local maxima are
found with an explicit plateau rule; the highest-energy one is the reported
spectral peak.

### 5.5 Active period and band-passing

`juno_ulf/bandpass.py`.  The perpendicular magnitude `|B_perp|` is filtered with
a zero-phase Butterworth band-pass (`sosfiltfilt`).  Its Hilbert envelope,
smoothed over 60 s, defines wave packets as the contiguous runs above
`0.30 x max(envelope)`; the last run is the reported active period.

---

## 6. Outputs

| Folder | Content |
|--------|---------|
| `output/cwt/` | CWT spectrograms of B_parallel and B_perp |
| `output/psd/` | PSD of both components; ion-cyclotron band on B_perp only |
| `output/bandpass/` | Band-passed \|B_perp\| with the active period shaded |
| `output/field/` | B_r, B_theta, B_phi with the active period shaded |

---

## 7. Notes

* B_parallel has no local maximum in 3–50 mHz for either event, so its PSD is
  plotted without the ion-cyclotron band.
* Band-passing a non-negative quantity such as `|B_perp|` removes its mean, so
  the filtered curve oscillates about zero; the envelope of that curve, not its
  value, carries the wave amplitude.
* The JSS-RTP product already provides spherical components, so `Br/Btheta/Bphi`
  are used directly without further rotation.

## 8. License

Add the license required by the journal / your institution before public release.

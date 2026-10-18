# Implemented methods and interpretation limits

## Magnetic-field decomposition

Use nonoverlapping, consecutive 1200-s vector-mean blocks, including the last
partial block. Bparallel is the signed field projection on the block-mean
direction. Bperp is the nonnegative magnitude of the perpendicular residual.
It is a scalar magnitude, not two signed transverse components. PSD then
subtracts the full analysis-window scalar mean of each component.

## Calibrated one-sided wavelet PSD

The wavelet is `cmor1.5-1.0`; dt=1 s. Final PSD uses explicit PyWavelets CWT
precision=20 with convolution. For scale s_j, exclude
`ceil(sqrt(1.5)*s_j)` samples at each end. For the retained N_j samples,

```text
Q_j = sum_k |h_j[k]|^2
S_B(f_j) = 2*dt/(N_j*Q_j) * sum_n |W_j[n]|^2
```

h_j is the actual discrete implementation's impulse response, not an ideal
continuous Morlet formula. Units are nT²/Hz. Do not additionally divide by
frequency, scale, or frequency-grid spacing. The output is frequency-smoothed
by the wavelet response; integration over centre frequencies is not asserted
to recover time-series variance exactly. The conventional e^-2 COI exclusion
does not make all finite-record boundary effects vanish.

Both event grids have 200 logarithmic samples. Main: 3–50 mHz; quiet: 1–50 mHz.
Figures display 3.5–50 mHz. PSD computation uses the complete analysis window,
not only the active period. The preserved yellow IQR band uses the main
active period 09:17:40–09:46:44 UT; its endpoints are 0.0038250449843684094 and
0.006413020737400884 Hz. The source formula is fci=0.0152*|B|/16 in Hz for B in nT.

## Peaks and the display range

Interior local maxima are discrete, plateau-aware and exclude endpoints.
No interpolation/fitting is used for peak frequency, and the local-maximum
rule is not a statistical significance test. The main perpendicular local
peak is 4.211917906742496 mHz, about 25.1397217457 nT²/Hz. The parallel local
peak near 3.359238216634793 mHz remains in the complete CSV but falls outside
the displayed range. At the perpendicular peak's frequency, the parallel
PSD is approximately 24.6503936099 nT²/Hz. These numbers alone do not establish
strong perpendicular power dominance or uniquely identify the driving mechanism.

`display_curve.py` interpolates intersections of the plotted log-frequency
line segments with the axis edges. This is clipping geometry only; the
intersections are never new PSD estimates or inputs to peak identification.

## CWT amplitude displays

Fig3 and S3(a–c) retain `abs(W)/sqrt(scale)` in nT from the original non-demeaned
components, with no PSD COI masking or filter-energy calibration. Optional
`cwt` recomputation uses explicit precision 12. Archived arrays are retained
for stable replotting and are not silently replaced if a library version
changes its results. Fig3 retains the original displayed shading endpoints
09:17:38–09:46:52 UT; this small display difference is recorded rather than
silently changing the current manuscript figure. The S1 and PSD IQR interval
is the manuscript's 09:17:40–09:46:44 UT interval.

## Band-pass and theory

S1 uses the original block projection and a fourth-order 2–5 mHz Butterworth
SOS filter, applied forward and backward (`sosfiltfilt`). Its highlighted
interval is supplied from the manuscript's visual identification of the
largest perturbations, not newly estimated by a packet-widening algorithm.

The supplied three- and eighth-degree root scripts are retained without
rederiving their polynomial coefficients or validating the selected branch
against final Fig1. Their original comments include generic/example wording;
execution alone does not certify physical-model or final-figure equivalence.

## Validation

The PSD implementation checks discrete impulse alignment, precision 18 vs20,
128 white-noise replicates, sinusoidal filter response, and a Welch reference.
Welch uses 1 Hz sampling, Hann, nperseg=nfft=2048, noverlap=1024, detrend=False,
and one-sided density scaling after full-window mean removal. Wavelet and
Welch smoothing differ; identical peak grids or integrals are not required.
Actual execution results for this release candidate are in VALIDATION.md.

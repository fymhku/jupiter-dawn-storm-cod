# Implemented analysis definitions

This note describes what the packaged code calculates. It does not change
the original numerical definitions or resolve methodological questions.

## Field decomposition

Consecutive, **non-overlapping** 20-minute blocks begin at the analysis-window
start. Within each block, `b0 = mean(B)/|mean(B)|`, `B_parallel = B dot b0`,
and `B_perp = sqrt(max(|B|^2 - B_parallel^2, 0))`.
`B_perp` is a non-negative magnitude; it is not a signed transverse vector
component. It is not a continuously moving mean. Block changes can produce
jumps in the projected series and require consideration when interpreting
spectra and filtered signals.

## CWT and spectra

Complex Morlet `cmor1.5-1.0`, 200 logarithmically spaced frequencies and
assumed sampling rate 1 Hz. Spectrogram amplitudes use
`abs(coef)/sqrt(scale)`. The routine called PSD instead computes
`mean(abs(raw_coef)**2, axis=time)` over the **full input window**.
No additional per-Hz density normalization is implemented. Existing plot
labels are retained; this code quantity should not be described as a
calibrated per-Hz PSD without establishing the required normalization.

The main grid remains 3-50 mHz, as in the supplied modular configuration;
the existing quiet-period scripts use 1-50 mHz. Selecting another grid
changes the discrete frequency bins and can change reported peak values.
Local maxima are found by the supplied plateau rule. The highest local
maximum and the global largest value are different reported quantities.

## Ion frequency and percentile band

`f_ci = 0.0152 * |B[nT]| / ion_mass_amu` Hz, with the singly charged
equivalent mass defaulting to 16 amu. This gives the same charge-to-mass
ratio for O+ and S2+. The yellow band is the interval between the 25th and
75th percentiles of the local `f_ci` samples in the configured active period,
not 0.25-0.75 times a representative cyclotron frequency.

## Band-pass and selected interval

Fourth-order Butterworth SOS filtering at 2-5 mHz, applied forward and
backward (`sosfiltfilt`). The May 16 shading follows the visually selected
interval 09:17:40-09:46:44 UTC. A Hilbert envelope, smoothed over 60 seconds,
and the last run above 30% of its maximum are retained as diagnostics;
they do not set the configured shading automatically.

## Theory scripts

`theory/solve.py` tracks a selected root of
`x^3 - (1-q^2)x^2 - 2q^2 x + 1.375q^2 = 0`, for q in [0,2].
`theory/solve_8.py` uses coefficients
`[1,0,-(1+2q^2+q^4),0,q^4+2q^2,0,q^2-q^4,0,-q^4]`.
These are the actual coefficient arrays in the supplied code. Their
relation to the manuscript equations and model parameters must be checked
before claiming an exact reconstruction of the final theory figure.

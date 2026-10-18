# Data dictionary

All timestamps are UTC. Magnetic fields are in nT, frequencies in Hz,
sampling/scale quantities in samples at dt=1 s. No missing-value code is
used in the included event windows.

## Time series

`Br_nT`, `Btheta_nT`, `Bphi_nT`: prepared spherical field components.
`Bparallel_nT`: signed block projection; `Bperp_nT`: nonnegative perpendicular
residual magnitude; `Btotal_nT`: vector magnitude; `fci_Hz`: original ion
cyclotron frequency conversion. `*_demeaned_nT` subtract the full-window mean.

## PSD tables

`frequency_Hz`: original centre-frequency grid, unchanged by display limits.
`scale_samples`: CWT scale; `COI_margin_samples`: per-end exclusion count;
`N_valid`: retained sample count; `valid_start_UTC` and `valid_stop_UTC`:
frequency-dependent retained interval. `*_PSD_nT2_per_Hz`: final precision20
calibrated one-sided density. `*_PSD_precision*_nT2_per_Hz`: convergence grids.
`*_old_raw_mean_power_nT2` and `*_valid_raw_mean_power_nT2`: uncalibrated
squared-coefficient powers, not densities. `Q_filter_energy`: actual discrete
response energy; `calibration_factor_seconds`: 2*dt/Q.

## Display NPZ files

`frequency_Hz` and `time_seconds` give axes; `Bperp_amplitude_nT` and
`Bparallel_amplitude_nT` have shape (200, event sample count) and contain
abs(W)/sqrt(scale), not PSD. `fci_Hz` is the overlay. The quiet NPZ also
contains UTC_ns and field columns. Load with NumPy; do not enable pickle.

## Other tables

Welch CSVs contain frequency_Hz and PSD_nT2_per_Hz. Historical raw-spectrum
CSVs contain frequency_Hz and mean_wavelet_power_nT2. S1 reference contains
UTC, Bperp_nT and bandpassed_Bperp_nT.

## Exact distributed CSV column lists

- `FigS1_reference.csv`: `UTC`, `Bperp_nT`, `bandpassed_Bperp_nT`
- `main_Bparallel_Welch.csv`: `frequency_Hz`, `PSD_nT2_per_Hz`
- `main_Bparallel_original_raw_spectrum.csv`: `frequency_Hz`, `mean_wavelet_power_nT2`
- `main_Bperp_Welch.csv`: `frequency_Hz`, `PSD_nT2_per_Hz`
- `main_Bperp_original_raw_spectrum.csv`: `frequency_Hz`, `mean_wavelet_power_nT2`
- `main_psd.csv`: `frequency_Hz`, `scale_samples`, `COI_margin_samples`, `N_valid`, `valid_start_UTC`, `valid_stop_UTC`, `Bperp_old_raw_mean_power_nT2`, `Bperp_PSD_precision12_nT2_per_Hz`, `Bperp_PSD_precision16_nT2_per_Hz`, `Bperp_PSD_precision18_nT2_per_Hz`, `Bperp_PSD_precision20_nT2_per_Hz`, `Bperp_PSD_nT2_per_Hz`, `Bperp_valid_raw_mean_power_nT2`, `Bparallel_old_raw_mean_power_nT2`, `Bparallel_PSD_precision12_nT2_per_Hz`, `Bparallel_PSD_precision16_nT2_per_Hz`, `Bparallel_PSD_precision18_nT2_per_Hz`, `Bparallel_PSD_precision20_nT2_per_Hz`, `Bparallel_PSD_nT2_per_Hz`, `Bparallel_valid_raw_mean_power_nT2`, `Q_filter_energy`, `calibration_factor_seconds`
- `main_timeseries.csv`: `UTC`, `Br_nT`, `Btheta_nT`, `Bphi_nT`, `Bparallel_nT`, `Bperp_nT`, `Btotal_nT`, `fci_Hz`, `Bperp_demeaned_nT`, `Bparallel_demeaned_nT`
- `quiet_Bparallel_Welch.csv`: `frequency_Hz`, `PSD_nT2_per_Hz`
- `quiet_Bparallel_original_raw_spectrum.csv`: `frequency_Hz`, `mean_wavelet_power_nT2`
- `quiet_Bperp_Welch.csv`: `frequency_Hz`, `PSD_nT2_per_Hz`
- `quiet_Bperp_original_raw_spectrum.csv`: `frequency_Hz`, `mean_wavelet_power_nT2`
- `quiet_psd.csv`: `frequency_Hz`, `scale_samples`, `COI_margin_samples`, `N_valid`, `valid_start_UTC`, `valid_stop_UTC`, `Bperp_old_raw_mean_power_nT2`, `Bperp_PSD_precision12_nT2_per_Hz`, `Bperp_PSD_precision16_nT2_per_Hz`, `Bperp_PSD_precision18_nT2_per_Hz`, `Bperp_PSD_precision20_nT2_per_Hz`, `Bperp_PSD_nT2_per_Hz`, `Bperp_valid_raw_mean_power_nT2`, `Bparallel_old_raw_mean_power_nT2`, `Bparallel_PSD_precision12_nT2_per_Hz`, `Bparallel_PSD_precision16_nT2_per_Hz`, `Bparallel_PSD_precision18_nT2_per_Hz`, `Bparallel_PSD_precision20_nT2_per_Hz`, `Bparallel_PSD_nT2_per_Hz`, `Bparallel_valid_raw_mean_power_nT2`, `Q_filter_energy`, `calibration_factor_seconds`
- `quiet_timeseries.csv`: `UTC`, `Br_nT`, `Btheta_nT`, `Bphi_nT`, `Bparallel_nT`, `Bperp_nT`, `Btotal_nT`, `fci_Hz`, `Bperp_demeaned_nT`, `Bparallel_demeaned_nT`

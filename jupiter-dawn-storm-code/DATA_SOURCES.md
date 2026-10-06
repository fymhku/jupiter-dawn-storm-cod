# Data sources and input format

## Public observation archives

- Juno magnetometer Jupiter archive, NASA Planetary Data System, Planetary
  Plasma Interactions Node: https://doi.org/10.17189/1519711
  (mission index: https://pds-ppi.igpp.ucla.edu/mission/juno).
- Hubble Space Telescope ultraviolet observations, program GO-14634,
  Mikulski Archive for Space Telescopes: https://doi.org/10.17909/t9-69wa-4q90.

These links identify observation sources. HST reduction and auroral-map
construction are not implemented by the scripts supplied in this package.

## MAG input contract

Supported scripts consume locally prepared 1-second MAG time series in
spherical JSS/RTP components. They do not parse arbitrary PDS-native files.
No raw-to-JSS/RTP conversion script was found among the supplied files.
Before claiming reproduction directly from raw archive products, document
the selected native products and the preprocessing/coordinate conversion.

Expected daily names:

```text
JUNO_MAG_JSS_1s_2017_05_10.txt
JUNO_MAG_JSS_1s_2017_05_16.txt
```

Files must contain one header line followed by four whitespace-separated
columns, in this order:

```text
UT_Time                  Br(nT)  Btheta(nT)  Bphi(nT)
YYYY-MM-DD/HH:MM:SS       <Br>    <Btheta>    <Bphi>
```

The second line describes the format and is not a numerical sample. Times
are UTC; fields are in nT. The CWT and filter assume regular 1 Hz sampling.
Window boundaries are inclusive. Validate gaps, duplicate timestamps,
missing samples and non-finite values before running the analysis.
Supply genuine measurements, not an interpolated synthetic example, when
comparing with the manuscript. Place inputs in `data/mag/` or set
`JUNO_MAG_DIR`. No observation data are included in the release ZIP.

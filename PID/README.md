# PID: Particle Identification

Particle ID separates conversion electrons from muons and other particles that leave tracks in the signal momentum window. The analysis uses the track PID score from EventNtuple (`trkpid.result`), which combines tracker and calorimeter information. The selection also requires a calorimeter cluster with nonzero energy.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Cut: `good_trkpid` in `analyze.py`. It requires PID score > `trkpid_threshold` and an event calorimeter cluster with `energyDep_` > 0
- Threshold: `trkpid_threshold` in `VERSION_CUTS` (`config.py`, `config_1BB.py`). Current values range from 0.525 to 0.61, depending on the version
- Systematic: `PID_Efficiency`, ±5%, planned (`documentation/UncertaintyAnalysis.md`)

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- PID efficiency for CE and its uncertainty (`PID_Efficiency` is still planned)
- Threshold choice and its stability across EventNtuple versions
- Muon and pion rejection in the signal window

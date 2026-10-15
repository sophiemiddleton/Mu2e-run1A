# Tracker: Tracking, Momentum Scale and Resolution

The tracker measures the electron's momentum and time. Its resolution and absolute momentum scale set the width and position of the CE peak, and how far the DIO tail reaches into the signal region. This directory covers the track-level selection and the detector response that every momentum PDF depends on. The track quality classifier has its own directory, [`trkqual/`](../trkqual/).

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Cuts in `analyze.py`, with thresholds in `VERSION_CUTS` (`config.py`, `config_1BB.py`):
  - `has_hits`: at least `nhits_threshold` active hits (17 to 21)
  - `within_t0err`: `t0err` < `t0err_threshold` (0.78 to 0.85)
  - `pz_over_pt`: `tandip_lower` < tanDip < `tandip_upper` (lower 0.46 to 0.575, upper 0.85 to 0.99)
  - `within_d0`: `d0_cut_lower` < $d_0$ < `d0_cut_upper` (currently open)
- Detector response: the energy loss ⊗ resolution kernel used to build the CE and DIO momentum PDFs (`helper.py`, `momentum_pdf_builder.py`)
- Systematics (`documentation/UncertaintyAnalysis.md`): `Abs_Mom_Scale` (±3%, implemented), `Mom_Res_Scale` (±10%, planned), `Track_Efficiency` (±2%, planned), `Timing_Cal` (±50 ps, planned)

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- Momentum resolution systematic (`Mom_Res_Scale` is still planned)
- Tracking efficiency and its uncertainty
- Timing calibration systematic
- Optimizing the hit-count, `t0err` and tanDip cuts

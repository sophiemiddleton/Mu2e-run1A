# CRV: Cosmic Ray Veto

The cosmic ray veto surrounds the detector solenoid. It vetoes tracks that come with a CRV coincidence close in time. Its inefficiency, and any gaps in coverage, set the residual cosmic background. The veto window also costs signal efficiency through accidental coincidences from beam-induced activity.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Cut: `no_crv_veto` in `analyze.py`. A track is vetoed if any CRV coincidence (`crvcoincs.time`) has $0 < t_{\text{trk}} - t_{\text{CRV}} < 150$ ns, using the track time at the tracker front
- Window: `crv_veto_dt_lower` and `crv_veto_dt_upper` in `VERSION_CUTS` (`config.py`, `config_1BB.py`), currently 0 and 150 ns for all versions
- Systematic: `CRV_Efficiency`, ±10%, implemented (`documentation/UncertaintyAnalysis.md`)
- Related: the cosmic yield and shape (`Cosmic_mix`, OffSpill control region). See `documentation/ControlRegions.md`

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- Veto efficiency measurement and its systematic
- Optimizing the veto time window
- Signal dead time from accidental coincidences

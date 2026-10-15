# RMC: Radiative Muon Capture Background

In radiative muon capture, $\mu^- \text{Al} \rightarrow \nu_\mu \gamma\, \text{Mg}^*$, a photon converts to $e^+e^-$, externally or internally. The photon endpoint is below the CE energy, so RMC only enters the signal region through resolution tails. Its time follows the muonic-aluminum lifetime. Conversions give $e^+$ and $e^-$ in equal numbers, and DIO gives no $e^+$, so late-time $e^+$ can constrain RMC.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Single-component fit: `component_fits/rmc.py`
- Normalization: `RMCExt_mix` and `RMCInt_mix` (0N phase-space samples). The 1N samples give nothing above the 80 MeV generator threshold.
- Control region: late-time $e^+$ constraint in `run_chain.py --rmc-control` (`documentation/ControlRegions.md`)

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- Photon spectrum endpoint and shape (Plestid 0N model)
- $e^+$/$e^-$ charge-efficiency systematic (`charge_syst`) for the control region
- Size of the RMC contribution in the signal window

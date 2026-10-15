# RPC: Radiative Pion Capture Background

Pions in the beam stop in the target and capture: $\pi^- N \rightarrow \gamma N^*$. The photon either converts in material (**external**) or produces an $e^+e^-$ pair directly (**internal** conversion). RPC is prompt and follows the pion lifetime, so the delayed signal time window suppresses it. Pion survival to late times sets the residual rate.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Single-component fit: `component_fits/rpc.py`
- Normalization: `RPCExt_mix` and `RPCInt_mix`. Expected yield = POT × stopped π/POT × survival weights × RPC branching fraction × the fraction of the Bistirlich spectrum above 50 MeV, times 0.0069 for internal
- Control region: $e^+$ fit in `control_region.py` (`documentation/ControlRegions.md`)

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- Momentum and time shape: replace the Gaussian with a data-driven model
- Pion survival and extinction systematics
- Use of the $e^+$ control region

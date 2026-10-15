# CE: Conversion Electron Signal

The signal is a mono-energetic electron at about 104.97 MeV/c from $\mu^- \text{Al} \rightarrow e^- \text{Al}$. After energy loss in the target and detector, and after tracker resolution, the reconstructed momentum peak is shifted and has a low-side tail. Its time follows the muonic-aluminum lifetime (about 864 ns).

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Single-component fit: `component_fits/ce.py`
- Model definition: `model/physics_components*.py` (look for CE)
- Normalization: `CeM_mix` in the dataset registry. Expected yield = $N_{\text{stop}\,\mu}$ × captures per stopped muon × $R_{\mu e}$
- Background: `documentation/LikelihoodDefinition.md`

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- Momentum response function (energy loss ⊗ resolution) and its parameterization
- Momentum scale and resolution systematics
- Acceptance × efficiency as a function of momentum and time

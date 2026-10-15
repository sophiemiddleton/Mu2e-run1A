# DIO: Decay-in-Orbit Background

A muon bound in aluminum can decay in orbit. Nuclear recoil gives the electron spectrum a tail that reaches the CE endpoint, falling roughly as $(E_{\text{end}} - E)^5$. This makes DIO the irreducible background just below the signal peak. Its time distribution matches CE.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Single-component fit: `component_fits/dio.py`
- Custom momentum PDF: `DIO_custom_model_2025` in `momentum_pdf_builder.py`
- Normalization: `DIO_mix` (`DIOtail95Mix1BB`). Expected yield = $N_{\text{stop}\,\mu}$ × DIO fraction × the fraction of the Czarnecki/Szafron spectrum above 95 MeV
- Background: `documentation/LikelihoodDefinition.md`

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- Shape model near the endpoint and its theory uncertainty
- Correlation with the CE momentum response (shared energy loss and resolution)
- Yield systematics

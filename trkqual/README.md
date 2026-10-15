# trkqual: Track Quality

The track quality score is a classifier output stored in EventNtuple (`trkqual.result`). It separates well-reconstructed tracks from tracks with poor fits, wrong hit assignments or badly measured momentum. Badly measured tracks fill the high-momentum tail, so this cut controls how much DIO leaks into the signal region.

## Inputs to the main fit

Results from this directory feed the pyfitter fit
([github.com/sophiemiddleton/pyfitter](https://github.com/sophiemiddleton/pyfitter)).
The current implementation there:

- Cut: `good_trkqual` in `analyze.py`. It requires quality > `trkqual_threshold`
- Threshold: `trkqual_threshold` in `VERSION_CUTS` (`config.py`, `config_1BB.py`). Current values range from 0.145 to 0.228, depending on the version
- Systematic: there is no dedicated term yet. The effect is partly covered by `Track_Efficiency` (±2%, planned)

## Studies

Add each study as a subdirectory with its own `README.md`. See
[Adding your analysis](../README.md#adding-your-analysis) in the top-level README.

| Study | Description | Contact |
| :--- | :--- | :--- |
| | | |

## Open items

- Threshold optimization: CE efficiency against DIO tail suppression
- Validating the classifier on data, and the efficiency systematic
- Threshold stability across EventNtuple versions

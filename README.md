# Mu2e Run-1A Muon Conversion Analysis

This is the organizational repository for the Mu2e Run-1A search for neutrinoless
muon-to-electron conversion in the field of an aluminum nucleus,
$\mu^- N \rightarrow e^- N$.

This repository holds the studies that feed the analysis: per-process background
and signal studies, control-region work, and related checks. The main fit and
limit-setting framework is **not** here. It lives in
[**pyfitter**](https://github.com/sophiemiddleton/pyfitter). For anything about
running the fit, see that repository and its `documentation/` directory.

---

## The analysis in brief

### Signal

A muon stopped in the aluminum target forms muonic aluminum. If it converts, it
emits a mono-energetic electron (CE) with momentum about 104.97 MeV/c. The
electron time follows the muonic-aluminum lifetime, $\tau_{\mu\text{Al}} \approx 864$ ns.

### Backgrounds

| Process | Origin | Momentum shape | Time shape |
| :--- | :--- | :--- | :--- |
| **DIO**, decay in orbit | Bound muon decays $\mu^- \rightarrow e^- \bar\nu_e \nu_\mu$; the recoil tail reaches the CE endpoint | Falling spectrum near the endpoint (Czarnecki/Szafron) | $e^{-t/\tau_{\mu\text{Al}}}$, same as CE |
| **RPC**, radiative pion capture | $\pi^- N \rightarrow \gamma N^*$; the photon converts (external) or produces an $e^+e^-$ pair directly (internal) | Broad, up to about 140 MeV/c | Prompt, follows the pion lifetime |
| **RMC**, radiative muon capture | $\mu^- N \rightarrow \nu_\mu \gamma N^*$, with $\gamma \rightarrow e^+e^-$ (external or internal) | Ends below the CE endpoint | $e^{-t/\tau_{\mu\text{Al}}}$ |
| **Cosmics** | Cosmic-ray-induced electrons that get past the cosmic ray veto | Roughly flat across the window | Uniform |

### Strategy

The analysis is an unbinned, extended maximum-likelihood **shape fit** in
reconstructed track momentum and time (1D in momentum, or 2D in momentum × time):

$$
\mathcal{L} = \mathcal{P}_{\text{Poisson}}(N_{\text{obs}}; N_{\text{exp}})
\prod_{i}^{N_{\text{obs}}} \sum_j w_j\, f_j(p_i, t_i)
\cdot \prod_k \mathcal{C}_k(\theta_k)
$$

Each process $j$ has a PDF $f_j$ that combines its theory spectrum, its
acceptance and efficiency, and the detector response: energy loss plus resolution.
Systematic uncertainties enter as constrained nuisance parameters $\theta_k$. The
signal time window starts after the prompt beam flash, at 640 ns, and runs to 1650 ns.

The analysis chain in pyfitter has three stages:

1. **Selection.** Apply the analysis cuts to EventNtuple samples and record cut
   flows and efficiencies.
2. **Fit.** Build the per-process PDFs and yields, then fit with
   [`zfit`](https://zfit.readthedocs.io).
3. **Limit.** Run a frequentist CLs scan that gives a 90% CL upper limit on
   $R_{\mu e}$.

Expected yields are normalized to the run plan (beam mode, live time, POT,
stopped muons). Control regions give data-driven inputs: OffSpill cosmics for the
cosmic shape, and late-time $e^+$ for RMC.

Development currently uses the MDC2025 mock datasets (1BB beam mode).

### Where to read more (in pyfitter)

| Topic | pyfitter document |
| :--- | :--- |
| Environment setup and a first job | `documentation/GettingStarted.md` |
| Running selection → fit → limit | `documentation/RunChain.md` |
| Datasets, acceptances, expected yields | `documentation/DatasetRegistry.md` |
| Likelihood and per-process models | `documentation/LikelihoodDefinition.md` |
| Control regions | `documentation/ControlRegions.md` |
| Systematics | `documentation/UncertaintyAnalysis.md`, `documentation/SubsidiaryMeasurements.md` |
| Limits and interpretation | `documentation/CombineLikeSensitivityLimits.md`, `documentation/ResultsInterpretation.md` |

---

## Repository layout

Each physics process, and each detector or selection area, gets its own
top-level directory. Studies go in subdirectories underneath it.

### Physics processes

| Directory | Scope | pyfitter component | Dataset key(s) |
| :--- | :--- | :--- | :--- |
| [`CE/`](CE/) | Conversion electron signal: momentum response, acceptance, signal shape | `component_fits/ce.py` | `CeM_mix` |
| [`DIO/`](DIO/) | Decay-in-orbit background: endpoint spectrum, shape model, yield | `component_fits/dio.py` | `DIO_mix` (`DIOtail95Mix1BB`) |
| [`RPC/`](RPC/) | Radiative pion capture, internal and external conversion, $e^+$ control region | `component_fits/rpc.py` | `RPCExt_mix`, `RPCInt_mix` |
| [`RMC/`](RMC/) | Radiative muon capture, internal and external conversion, late-time $e^+$ control region | `component_fits/rmc.py` | `RMCExt_mix`, `RMCInt_mix` |

### Detector and selection

These areas cover the selection cuts and detector response that every process
shares. Cuts are applied in pyfitter's `analyze.py`. Their thresholds are set
per EventNtuple version in `VERSION_CUTS` (`config.py`, `config_1BB.py`).

| Directory | Scope | pyfitter cut(s) | Systematic(s) |
| :--- | :--- | :--- | :--- |
| [`Tracker/`](Tracker/) | Track selection, momentum scale and resolution, energy loss, timing | `has_hits`, `within_t0err`, `pz_over_pt`, `within_d0` | `Abs_Mom_Scale`, `Mom_Res_Scale`, `Track_Efficiency`, `Timing_Cal` |
| [`trkqual/`](trkqual/) | Track quality classifier threshold and efficiency | `good_trkqual` | covered by `Track_Efficiency` for now |
| [`PID/`](PID/) | Particle ID: electron selection, muon and pion rejection | `good_trkpid` | `PID_Efficiency` |
| [`CRV/`](CRV/) | Cosmic ray veto: efficiency, veto window, dead time | `no_crv_veto` | `CRV_Efficiency` |

### Directory READMEs

Every directory has a `README.md` with the same sections:

- **Summary.** The process or detector area, and what it does to the
  momentum and time shapes.
- **Inputs to the main fit.** Where it's implemented in pyfitter. For a
  process: the component fit, the PDF or model code, how the yield is
  normalized, and any control region. For a detector area: the cuts, their
  thresholds, and the systematics.
- **Studies.** A table listing each study subdirectory with a description and a
  contact. Add a row when you add a study.
- **Open items.** Work that still needs an owner. Check this list before you
  start a new study.

### Adding your analysis

1. **Create a new top-level directory** for your study, named after the process
   or topic (for example `Cosmics/`, `Momentum-Scale/`, `Trigger-Eff/`). If your
   study belongs to one of the existing areas, add a subdirectory there instead,
   for example `DIO/<your-study>/`.
2. **Add a `README.md`** to it that covers:
   - what the study measures and how it feeds the main fit (a shape, a yield, a
     constraint or a systematic)
   - the input datasets, by full dataset name and version
   - how to run it, with the exact commands
   - where the outputs go and what format they're in
   - the people responsible
3. **Keep large files out of git.** Don't commit ROOT, parquet or ntuple files, or
   large `.npz` outputs. Put them on `/exp/mu2e/data` or `/pnfs` and record the
   path in your README.
4. **Use the standard tools.** Read EventNtuple through
   [`pyutils`](https://github.com/Mu2e/pyutils) in the `pyenv ana` environment,
   so your selections and units match pyfitter.
5. **Open a pull request.** All changes go through review.

### Environment

```bash
mu2einit
pyenv ana
```

---

## Contact

Mu2e Analysis Group. To be added as a participant, contact sophie@fnal.gov.

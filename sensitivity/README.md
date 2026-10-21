# Sensitivity: Running pyfitter

The sensitivity (90% CL upper limit on $R_{\mu e}$) comes from
[pyfitter](https://github.com/Mu2e/pyfitter). One script, `run_chain.py`, runs
the whole chain:

1. **Selection.** Apply the analysis cuts to each dataset.
2. **Fit.** Build the per-process PDFs and fit them.
3. **Limit.** Run the CLs scan and compute the upper limit.

For full details see `documentation/RunChain.md` in pyfitter.

## 1. Set up

```bash
mu2einit
pyenv ana
git clone https://github.com/Mu2e/pyfitter.git
cd pyfitter
```

## 2. Run the chain

```bash
python run_chain.py --tag run1a --version run1a --group mix --loc tape --jobs 16
```

| Option | Meaning |
| :--- | :--- |
| `--tag` | Name of the output directory |
| `--version` | Cut set from `VERSION_CUTS` in `config.py`. Use `run1a` for the Run-1A analysis cuts |
| `--group mix` | Run all the mixed signal and background samples in `datasets.yaml` |
| `--loc` | Where the input files are: `tape`, `disk` or `local` |
| `--jobs` | Number of parallel workers |

The datasets, acceptances, expected yields and fit type (`1d` or `2d`) come from
`datasets.yaml`. Don't edit it by hand: regenerate it with
`python make_registry.py`.

Always pass `--version run1a`, including when you rerun part of the chain, so
every stage uses the same cuts.

Add `--dry-run` to print the commands without running them.

## 3. Find the results

Everything is written under `run1a/` (the `--tag` name):

| Path | Contents |
| :--- | :--- |
| `run1a/cutflows/` | Cut flow per dataset |
| `run1a/manifest.yaml` | Acceptances, cut efficiencies and expected yields |
| `run1a/fit/plots/` | Fit plots |
| `run1a/fit/output_card.yaml` | Fitted values passed to the limit |
| `run1a/limits/` | The upper limit |
| `run1a/logs/` | Logs for each stage, check here if something fails |

## Rerunning part of the chain

Selection is the slow step. Once it has run, you can redo later stages on their own:

```bash
# Redo the fit and limit, reusing the selection
python run_chain.py --tag run1a --version run1a --skip-selection

# Redo only the limit
python run_chain.py --tag run1a --version run1a --limit-only

# Selection and fit only, no limit
python run_chain.py --tag run1a --version run1a --group mix --skip-limit
```

Other useful options:

| Option | Meaning |
| :--- | :--- |
| `--fit-type 1d` / `2d` | Override the fit type in the registry |
| `--toys N` | Number of toys for the limit |
| `--asimov` | Fast asymptotic limit (Asimov) instead of toys |
| `--control-region` | Take the cosmic momentum shape from the OffSpill sample |
| `--rmc-control` | Constrain the RMC yield with late-time $e^+$ |

## Contact

Sophie Middleton (sophie@fnal.gov)

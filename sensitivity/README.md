# Sensitivity: Running pyfitter

The sensitivity (90% CL upper limit on $R_{\mu e}$) comes from
[pyfitter](https://github.com/Mu2e/pyfitter). One script, `run_chain.py`, runs
the whole chain:

1. **Selection.** Apply the analysis cuts to each dataset.
2. **Fit.** Build the per-process PDFs and fit them.
3. **Limit.** Run the CLs scan and compute the upper limit.

For full details see `documentation/RunChain.md` in pyfitter.

## 1. Set up

From the top of this repository, in every new shell:

```bash
source sensitivity/setup.sh
getToken            # only needed to read files from tape or disk
cd $PYFITTER_DIR    # for run_chain.py
```

`setup.sh` does the following:

1. Runs `mu2einit` if you haven't already.
2. Activates `pyenv ana 2.8.0`.
3. Clones pyfitter into `sensitivity/pyfitter` the first time. Git ignores that
   directory.
4. Exports `PYFITTER_DIR` and checks that every package and pyfitter file the scripts
   need is present.

Mu2e/pyfitter is private, so the clone uses ssh and needs a GitHub ssh key: test it
with `ssh -T git@github.com`. To change the defaults, set these before sourcing:

| Variable | Default | Use |
| :--- | :--- | :--- |
| `PYFITTER_DIR` | `sensitivity/pyfitter` | Use an existing pyfitter checkout instead of cloning |
| `PYFITTER_REPO` | `git@github.com:Mu2e/pyfitter.git` | Clone a fork instead |
| `PYFITTER_REF` | repository default | Branch, tag or commit to check out after cloning |
| `MU2E_ANA_VERSION` | `2.8.0` | Another `pyenv ana` version |

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

## Cut-and-count cross-check

`cutcount/` gives the limit from a simple count in a fixed momentum and time window,
instead of the full shape fit. It uses the same selection and the same statistics:

1. **Selection.** Calls pyfitter's `run_chain.py --skip-fit`, so `process.py` and
   `analyze.py` apply the cuts of `--version`.
2. **Count.** Counts the post-cut tracks of each dataset in the window. It scales them
   to the run as `run_chain.py` does: `expected × acceptance × N_window / N_generated`
   for MC, and `N_window × (on-spill time / cosmic livetime)` for cosmics.
3. **Limit.** `counting_stats.py` is a one-bin version of pyfitter's
   `simple_limit_combine_card.py`. It uses the same conventions: POI μ = N_CE / card
   CE yield, Combine lnN nuisances, and a one-sided q_μ. It reports Asimov expected
   limits and bands, a toy cross-check, or an empirical toy CLs scan
   (`--method empirical`).

By default the window is the `signal_region` of the cut set: 103.34–104.74 MeV/c and
640–1650 ns for `run1a`. Use `--p-window` and `--t-window` to change it.

```bash
source sensitivity/setup.sh   # sets PYFITTER_DIR

# Full chain: selection, count, limit
python sensitivity/cutcount/run_cutcount.py --tag run1a_count --version run1a --group mix --loc tape --jobs 16

# Reuse an existing run_chain tag (no selection), try another window
python sensitivity/cutcount/run_cutcount.py --tag $PYFITTER_DIR/run1a --version run1a --skip-selection \
    --p-window 103.6 104.9 --t-window 700 1650

# Limit only, from an existing card
python sensitivity/cutcount/counting_stats.py --card run1a_count/cutcount/count_card.yaml --method empirical
```

The window has to lie inside the preselection: 100–110 MeV/c and the
`within_t0_475` window (475–1650 ns for `run1a`). The script warns if it does not.

Everything is written to `<tag>/cutcount/`:

| File | Contents |
| :--- | :--- |
| `count_table.csv` | Counts, efficiencies and yields in the window, per dataset |
| `count_card.yaml` | One-bin card: yields, lnN systematics, R_μe normalisation |
| `count_limits_<method>.yaml` | Limits in events and R_μe |

Other options:

| Option | Meaning |
| :--- | :--- |
| `--mc-stat` | Add a 1/√N lnN for each dataset's MC statistics in the window |
| `--n-obs N` | Also compute an observed limit for N events in the window |
| `--freeze-nuisances` | No systematics |
| `--skip-limit` | Stop after writing the card |

## Summary web page

`make_page.py` writes `index.html`, a single page with no external dependencies. It shows
the momentum and time distributions per process, scaled to the run and with the counting
window shaded. It also shows the expected R_μe limits from the cut and count and from the
shape fit. Run it after `run_chain.py` and `cutcount/run_cutcount.py` on the same tag:

```bash
python sensitivity/make_page.py --tag $PYFITTER_DIR/run1a            # writes sensitivity/index.html
python sensitivity/make_page.py --tag $PYFITTER_DIR/run1a --out run1a.html --cutcount-dir run1a_count/cutcount
```

## Contact

Sophie Middleton (sophie@fnal.gov)

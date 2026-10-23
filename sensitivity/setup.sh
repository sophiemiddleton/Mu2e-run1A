# Set up the environment for the sensitivity scripts. Source it, from any directory:
#
#     source sensitivity/setup.sh
#
# It
#   1. sets up the Mu2e environment (mu2einit) if that has not been done,
#   2. activates the Mu2e python environment (pyenv ana <version>),
#   3. clones pyfitter if there is no checkout yet, and exports PYFITTER_DIR,
#   4. checks that the python packages and pyfitter pieces the scripts use are there.
#
# Override any of these before sourcing:
#   MU2E_ANA_VERSION  pyenv ana version (default 2.8.0, the version the scripts are tested with)
#   PYFITTER_DIR      existing pyfitter checkout to use (default: sensitivity/pyfitter, cloned if missing)
#   PYFITTER_REPO     repository to clone (default git@github.com:Mu2e/pyfitter.git, needs a GitHub ssh key)
#   PYFITTER_REF      branch, tag or commit to check out after cloning (default: the repository default)

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "setup.sh must be sourced: source ${0}" >&2
    exit 1
fi

_sens_setup() {
    local sens_dir
    sens_dir="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
    export SENSITIVITY_DIR="$sens_dir"

    # 1. Mu2e environment
    if ! declare -F pyenv > /dev/null; then
        echo "[setup] Setting up the Mu2e environment (mu2einit)"
        source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh || return 1
    fi

    # 2. Python environment
    local version="${MU2E_ANA_VERSION:-2.8.0}"
    echo "[setup] Activating pyenv ana ${version}"
    pyenv ana "$version" || return 1

    # 3. pyfitter checkout
    export PYFITTER_DIR="${PYFITTER_DIR:-$sens_dir/pyfitter}"
    if [ ! -f "$PYFITTER_DIR/run_chain.py" ]; then
        local repo="${PYFITTER_REPO:-git@github.com:Mu2e/pyfitter.git}"
        echo "[setup] Cloning $repo into $PYFITTER_DIR"
        git clone "$repo" "$PYFITTER_DIR" || {
            echo "[setup] Clone failed. Mu2e/pyfitter is private: check your GitHub ssh key (ssh -T git@github.com)" >&2
            return 1
        }
        if [ -n "$PYFITTER_REF" ]; then
            git -C "$PYFITTER_DIR" checkout "$PYFITTER_REF" || return 1
        fi
    fi
    echo "[setup] pyfitter: $PYFITTER_DIR ($(git -C "$PYFITTER_DIR" log -1 --format='%h %cs' 2>/dev/null))"

    # 4. Checks
    python - "$PYFITTER_DIR" << 'EOF' || return 1
import importlib, sys
from pathlib import Path
missing = []
for module in ('pyutils', 'zfit', 'hepstats', 'uproot', 'awkward', 'numpy', 'scipy', 'pandas', 'pyarrow', 'yaml', 'matplotlib'):
    try:
        importlib.import_module(module)
    except ImportError:
        missing.append(module)
pyfitter = Path(sys.argv[1])
for name in ('run_chain.py', 'process.py', 'analyze.py', 'config.py', 'normalization.py',
             'uncertainties_config.py', 'datasets.yaml', 'sensitivity_scan/simple_limit_combine_card.py'):
    if not (pyfitter / name).exists():
        missing.append(f'pyfitter/{name}')
if missing:
    print('[setup] MISSING: ' + ', '.join(missing), file=sys.stderr)
    sys.exit(1)
print(f'[setup] python {sys.version.split()[0]}: all packages and pyfitter files found')
EOF

    echo "[setup] Ready. For example:"
    echo "    python $sens_dir/cutcount/run_cutcount.py --tag run1a --version run1a --group mix --loc tape --jobs 16"
    echo "    python $sens_dir/make_page.py --tag run1a"
    echo "[setup] Reading files from tape or disk (--loc tape/disk) also needs a grid token: getToken"
}

_sens_setup
_sens_status=$?
unset -f _sens_setup
return $_sens_status

"""Cut-and-count sensitivity in a fixed momentum and time window.

Stage 1: selection. Runs pyfitter's run_chain.py --skip-fit, i.e. process.py --nofit with the
         analyze.py cuts of --version, so the selection is identical to the full fit. Skip it
         with --skip-selection to reuse any existing tag (including one made by run_chain.py).
Stage 2: count. From <tag>/manifest.yaml, counts the post-cut tracks of every dataset inside
         the window and scales them exactly as run_chain.py scales its fit yields:
             MC:     yield = expected * acceptance * N_window / N_generated_in_sample
             Cosmic: yield = N_window * acceptance (on-spill time / cosmic sample livetime)
         Writes a one-bin card with the pyfitter standard lnN systematics.
Stage 3: limit. Combine-style CLs on the card (counting_stats.py).

Outputs (in <tag>/cutcount/):
    count_table.csv          per-dataset counts, efficiencies and yields
    count_card.yaml          one-bin datacard
    count_limits_<method>.yaml, plots and npz

Example:
    python run_cutcount.py --pyfitter ~/pyfitter --tag run1a --version run1a --group mix --loc tape --jobs 16
    python run_cutcount.py --pyfitter ~/pyfitter --tag ~/pyfitter/run1a --version run1a --skip-selection
"""

import argparse
import csv
import logging
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

from counting_stats import add_stats_arguments, run as run_limits

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('run_cutcount')


def import_pyfitter(pyfitter_dir):
    """Make pyfitter's config and uncertainties_config importable."""
    if not (pyfitter_dir / 'run_chain.py').exists():
        raise FileNotFoundError(f'No run_chain.py in {pyfitter_dir}: pass --pyfitter or set PYFITTER_DIR')
    sys.path.insert(0, str(pyfitter_dir))


def default_window(version):
    """Signal region of the cut set: the same window as analyze.py's signal_region cut."""
    from config import get_version_cuts
    cuts = get_version_cuts(version)
    return ((cuts['signal_region_p_lower'], cuts['signal_region_p_upper']),
            (cuts['signal_region_t_lower'], cuts['signal_region_t_upper']), cuts)


def run_selection(pyfitter_dir, tag_dir, args, p_window, t_window):
    """Selection only, through run_chain.py so the cuts and manifest match the fit chain."""
    cmd = [sys.executable, str(pyfitter_dir / 'run_chain.py'),
           '--tag', str(tag_dir), '--version', args.version, '--registry', args.registry,
           '--loc', args.loc, '--jobs', str(args.jobs), '--skip-fit',
           # The preselection momentum and time range has to contain the counting window
           '--fitrange_low', str(min(args.fitrange_low[0], p_window[0])), str(min(args.fitrange_low[1], t_window[0])),
           '--fitrange_hi', str(max(args.fitrange_hi[0], p_window[1])), str(max(args.fitrange_hi[1], t_window[1]))]
    if args.group:
        cmd += ['--group', args.group]
    if args.datasets:
        cmd += ['--datasets', *args.datasets]
    if args.dry_run:
        cmd.append('--dry-run')
    logger.info('Selection: %s', ' '.join(cmd))
    return subprocess.run(cmd, cwd=str(pyfitter_dir)).returncode == 0


def generated_in_sample(cutflow_path):
    """Events before any cut ('No cuts' row): the denominator of the run_chain cut efficiency."""
    with open(cutflow_path, newline='') as f:
        rows = list(csv.DictReader(f))
    return float(rows[0]['events_passing'])


def canonical_process_name(process):
    """Same mapping as run_chain.py: RPCInt/RPCExt -> RPC, etc."""
    name = str(process).strip().lower()
    if name.startswith('rpc'):
        return 'RPC'
    if name.startswith('rmc'):
        return 'RMC'
    if name in ('cosmic', 'cosmics'):
        return 'Cosmic'
    if name in ('ce', 'cem'):
        return 'CE'
    return str(process)


def count_window(manifest, p_window, t_window):
    """Per-dataset counts and yields in the window, scaled as in run_chain.write_input_card."""
    rows = []
    for entry in manifest['datasets']:
        if entry.get('status') != 'ok' or not entry.get('process'):
            logger.warning('[%s] skipped (status %s, process %s)', entry['name'], entry.get('status'), entry.get('process'))
            continue
        # momentum is stored as float32: compare in float64 so the window edges are exact
        df = pd.read_parquet(entry['parquet'], columns=['momentum', 'time']).astype('float64')
        # Strict inequalities, as the signal_region cut in analyze.py
        in_window = ((df['momentum'] > p_window[0]) & (df['momentum'] < p_window[1]) &
                     (df['time'] > t_window[0]) & (df['time'] < t_window[1]))
        n_window = int(in_window.sum())
        process = canonical_process_name(entry['process'])
        acceptance = float(entry['acceptance'])
        if process == 'Cosmic':
            # Cosmic acceptance is on-spill time / sample livetime: scale the raw count
            efficiency = None
            net_acceptance = None
            expected_yield = n_window * acceptance
        else:
            efficiency = n_window / generated_in_sample(entry['cutflow'])
            net_acceptance = acceptance * efficiency
            expected_yield = float(entry['expected']) * net_acceptance if entry.get('expected') is not None else None
        if expected_yield is None:
            logger.warning('[%s] no expected yield in the manifest, omitted', entry['name'])
            continue
        rows.append({
            'name': entry['name'], 'process': process, 'n_postcut': len(df), 'n_window': n_window,
            'window_efficiency': efficiency, 'acceptance': acceptance, 'net_acceptance': net_acceptance,
            'expected': entry.get('expected'), 'yield': expected_yield,
        })
        if 0 < n_window < 10:
            logger.warning('[%s] only %d MC tracks in the window: the yield has a large MC statistical error', entry['name'], n_window)
    return rows


def write_card(tag_dir, rows, args, p_window, t_window, run_cfg):
    """One-bin card in the input_card.yaml layout, with the pyfitter standard lnN systematics."""
    from uncertainties_config import STANDARD_SYSTEMATICS

    processes = {}
    for row in rows:
        processes[row['process']] = processes.get(row['process'], 0.0) + row['yield']

    systematics = {}
    for name, syst in STANDARD_SYSTEMATICS.items():
        effects = {p: v for p, v in syst['processes'].items() if p in processes}
        if effects:
            systematics[name] = {'type': syst['type'], 'description': syst['description'], 'effects': effects}
    if args.mc_stat:
        # Optional: MC statistical error of each dataset's count, 1/sqrt(N), weighted into its process
        for row in rows:
            if row['n_window'] > 0 and processes[row['process']] > 0:
                frac = row['yield'] / processes[row['process']] / row['n_window'] ** 0.5
                systematics[f'mcstat_{row["name"]}'] = {
                    'type': 'lnN', 'description': f'MC statistics of {row["name"]} ({row["n_window"]} tracks)',
                    'effects': {row['process']: round(1.0 + frac, 4)}}

    signal_efficiency = sum(r['net_acceptance'] for r in rows if r['process'] == args.signal_process)
    run_cfg = run_cfg or {}
    card = {
        'name': f'{tag_dir.name}_count_card',
        'metadata': {'experiment': 'Mu2e', 'signal_process': args.signal_process, 'version': args.version,
                     'description': f'Cut-and-count card from run_cutcount.py for tag {tag_dir.name}'},
        'window': {'mom': [float(p_window[0]), float(p_window[1])], 'time': [float(t_window[0]), float(t_window[1])]},
        'processes': {name: {'yield': float(value)} for name, value in processes.items()},
        'systematics': systematics,
        'poi': {'name': args.signal_process, 'initial': processes.get(args.signal_process, 1.0)},
        'normalization': {
            'signal_efficiency': float(signal_efficiency),
            'stopped_muons': float(run_cfg['stopped_muons']) if run_cfg.get('stopped_muons') is not None else 5.588e15,
            'captures_per_stopped_muon': float(run_cfg.get('captures_per_stopped_muon', 0.609)),
        },
    }
    if run_cfg.get('stopped_muons') is None:
        logger.warning('No "run" section in the registry: R_mue uses the default N_muons = 5.588e15')
    if args.n_obs is not None:
        card['observed'] = float(args.n_obs)

    card_path = tag_dir / 'cutcount' / 'count_card.yaml'
    with open(card_path, 'w') as f:
        yaml.safe_dump(card, f, sort_keys=False)
    logger.info('Wrote %s', card_path)
    return card_path


def main(args):
    pyfitter_dir = Path(args.pyfitter).expanduser().resolve()
    import_pyfitter(pyfitter_dir)
    args.registry_given = args.registry is not None
    if args.registry is None:
        args.registry = str(pyfitter_dir / 'datasets.yaml')

    p_default, t_default, cuts = default_window(args.version)
    p_window = tuple(args.p_window) if args.p_window else p_default
    t_window = tuple(args.t_window) if args.t_window else t_default
    logger.info('Counting window: %.2f < p < %.2f MeV/c, %.0f < t < %.0f ns', *p_window, *t_window)

    tag_dir = Path(args.tag).expanduser().resolve()
    if not args.skip_selection:
        if not run_selection(pyfitter_dir, tag_dir, args, p_window, t_window):
            logger.error('Selection failed, see %s/logs', tag_dir)
            return 1
        if args.dry_run:
            return 0

    manifest_path = tag_dir / 'manifest.yaml'
    if not manifest_path.exists():
        logger.error('No manifest in %s: run the selection first (drop --skip-selection)', tag_dir)
        return 1
    with open(manifest_path) as f:
        manifest = yaml.safe_load(f)
    if manifest.get('version') not in (None, args.version):
        logger.warning('Tag was selected with --version %s, counting with --version %s', manifest['version'], args.version)

    # The window has to sit inside the preselection, or the counts are silently truncated
    mom_lo, mom_hi = manifest['fitrange_low'][0], manifest['fitrange_hi'][0]
    t0_lo, t0_hi = cuts['mom_range_t0_lower'], cuts['mom_range_t0_upper']
    if p_window[0] < mom_lo or p_window[1] > mom_hi:
        logger.warning('Momentum window %s extends outside the selected range [%s, %s]', p_window, mom_lo, mom_hi)
    if t_window[0] < t0_lo or t_window[1] > t0_hi:
        logger.warning('Time window %s extends outside the within_t0 cut [%s, %s]', t_window, t0_lo, t0_hi)

    (tag_dir / 'cutcount').mkdir(parents=True, exist_ok=True)
    rows = count_window(manifest, p_window, t_window)
    if not any(r['process'] == args.signal_process and r['yield'] > 0 for r in rows):
        logger.error('No %s signal in the window', args.signal_process)
        return 1
    table = pd.DataFrame(rows)
    table_path = tag_dir / 'cutcount' / 'count_table.csv'
    table.to_csv(table_path, index=False)
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        print('\n================== Cut-and-count window yields ==================')
        print(table.to_string(index=False))
        print(table.groupby('process')['yield'].sum().to_string())
        print('=' * 66 + '\n')
    logger.info('Wrote %s', table_path)

    # Stopped muons for R_mue: from the registry the tag was selected with, unless overridden
    registry = args.registry if args.registry_given else manifest.get('registry', args.registry)
    with open(registry) as f:
        run_cfg = (yaml.safe_load(f) or {}).get('run')
    card_path = write_card(tag_dir, rows, args, p_window, t_window, run_cfg)

    if args.skip_limit:
        return 0
    args.output_dir = str(tag_dir / 'cutcount')
    run_limits(card_path, args)
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Cut-and-count sensitivity in a fixed momentum and time window')
    parser.add_argument('--pyfitter', default=os.environ.get('PYFITTER_DIR'), required=os.environ.get('PYFITTER_DIR') is None,
                        help='pyfitter checkout (default: $PYFITTER_DIR)')
    parser.add_argument('--tag', required=True, help='Output tag directory (an existing run_chain tag works with --skip-selection)')
    parser.add_argument('--version', default='run1a', help='Cut set from VERSION_CUTS in pyfitter config.py')
    parser.add_argument('--registry', default=None, help='Dataset registry (default: <pyfitter>/datasets.yaml)')
    parser.add_argument('--group', default=None, help='Registry group, e.g. mix')
    parser.add_argument('--datasets', nargs='+', default=None, help='Registry dataset names')
    parser.add_argument('--loc', default='tape', help='File location: disk, tape or local')
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--fitrange_low', type=float, nargs=2, default=[100, 475], help='Preselection minimum mom, time (widened to contain the window)')
    parser.add_argument('--fitrange_hi', type=float, nargs=2, default=[110, 1650], help='Preselection maximum mom, time (widened to contain the window)')
    parser.add_argument('--p-window', type=float, nargs=2, default=None, metavar=('LO', 'HI'), help='Momentum window in MeV/c (default: signal_region_p of the cut set)')
    parser.add_argument('--t-window', type=float, nargs=2, default=None, metavar=('LO', 'HI'), help='Time window in ns (default: signal_region_t of the cut set)')
    parser.add_argument('--signal-process', default='CE')
    parser.add_argument('--mc-stat', action='store_true', help='Add a 1/sqrt(N) lnN per dataset for the MC statistics in the window')
    parser.add_argument('--skip-selection', action='store_true', help='Reuse <tag>/manifest.yaml, cut flows and parquet')
    parser.add_argument('--skip-limit', action='store_true', help='Stop after writing the card')
    parser.add_argument('--dry-run', action='store_true', help='Print the selection command without running it')
    add_stats_arguments(parser)
    sys.exit(main(parser.parse_args()))

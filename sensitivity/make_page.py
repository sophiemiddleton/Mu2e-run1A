"""Build a self-contained web page with the momentum and time distributions and the
cut-and-count and shape-fit sensitivities of one tag.

Inputs (all from the tag made by run_chain.py and cutcount/run_cutcount.py):
    <tag>/manifest.yaml                     datasets, acceptances, expected yields, parquet paths
    <tag>/cutcount/count_card.yaml          counting window
    <tag>/cutcount/count_limits_asymptotic.yaml   cut-and-count limit
    <tag>/limits/card_{2d,1d}_cls_limits.npz      shape-fit limit

Each track is weighted to the run exactly as in the counting (run_cutcount.py):
    MC:     expected * acceptance / N_generated_in_sample
    Cosmic: acceptance (on-spill time / cosmic sample livetime)

The page has no external dependencies, so it opens from disk or any web server.

Usage:
    python make_page.py --tag ~/pyfitter/run1a
    python make_page.py --tag ~/pyfitter/run1a --out run1a.html --cutcount-dir run1a_count/cutcount
"""

import argparse
import csv
import datetime
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# Fixed order: each process keeps its color whatever is present
PROCESS_ORDER = ['CE', 'DIO', 'Cosmic', 'RPC', 'RMC']
MOM_BIN = 0.1    # MeV/c
TIME_BIN = 25.0  # ns


def canonical_process_name(process):
    """Same mapping as run_chain.py and run_cutcount.py."""
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


def generated_in_sample(cutflow_path):
    with open(cutflow_path, newline='') as f:
        return float(next(csv.DictReader(f))['events_passing'])


def histograms(manifest, p_window, t_window, mom_edges, time_edges):
    """Weighted histograms per process: all tracks, and tracks inside the other variable's window."""
    out = {}
    for entry in manifest['datasets']:
        if entry.get('status') != 'ok' or not entry.get('process'):
            continue
        process = canonical_process_name(entry['process'])
        df = pd.read_parquet(entry['parquet'], columns=['momentum', 'time'])
        if process == 'Cosmic':
            weight = float(entry['acceptance'])
        elif entry.get('expected') is None:
            continue
        else:
            weight = float(entry['expected']) * float(entry['acceptance']) / generated_in_sample(entry['cutflow'])
        # float64, as in run_cutcount.py: a float32 compare rounds the window edges and moves tracks
        p, t = df['momentum'].to_numpy(np.float64), df['time'].to_numpy(np.float64)
        in_p = (p > p_window[0]) & (p < p_window[1])
        in_t = (t > t_window[0]) & (t < t_window[1])
        h = out.setdefault(process, {'mom_all': 0, 'mom_win': 0, 'time_all': 0, 'time_win': 0,
                                     'fit_yield': 0.0, 'window_yield': 0.0, 'window_tracks': 0, 'tracks': 0})
        h['mom_all'] = h['mom_all'] + np.histogram(p, mom_edges)[0] * weight
        h['mom_win'] = h['mom_win'] + np.histogram(p[in_t], mom_edges)[0] * weight
        h['time_all'] = h['time_all'] + np.histogram(t, time_edges)[0] * weight
        h['time_win'] = h['time_win'] + np.histogram(t[in_p], time_edges)[0] * weight
        h['fit_yield'] += len(df) * weight
        h['window_yield'] += int(np.sum(in_p & in_t)) * weight
        h['window_tracks'] += int(np.sum(in_p & in_t))
        h['tracks'] += len(df)
    for h in out.values():
        for key in ('mom_all', 'mom_win', 'time_all', 'time_win'):
            h[key] = [float(v) for v in np.asarray(h[key], dtype=float)]
    return {p: out[p] for p in PROCESS_ORDER if p in out} | {p: v for p, v in out.items() if p not in PROCESS_ORDER}


def count_limit(cutcount_dir):
    path = cutcount_dir / 'count_limits_asymptotic.yaml'
    if not path.exists():
        return None
    with open(path) as f:
        res = yaml.safe_load(f)
    asimov = res['asimov']
    return {
        'method': 'Asymptotic CLs (Asimov), one bin',
        'events': asimov['band_+0']['events'], 'rmue': asimov['band_+0']['rmue'],
        'rmue_lo': asimov['band_-1']['rmue'], 'rmue_hi': asimov['band_+1']['rmue'],
        'efficiency': res['signal_efficiency'], 'background': res['background_yield'],
        'source': str(path),
    }


def shape_limit(tag_dir):
    for dim in ('2d', '1d'):
        path = tag_dir / 'limits' / f'card_{dim}_cls_limits.npz'
        if path.exists():
            d = np.load(path)
            events, rmue = d['asimov_events'], d['asimov_rmue']
            return {
                'method': f'Asymptotic CLs (Asimov), {dim.upper()} shape fit',
                'events': float(events[2]), 'rmue': float(rmue[2]),
                'rmue_lo': float(rmue[1]), 'rmue_hi': float(rmue[3]),
                'efficiency': float(d['signal_efficiency']), 'background': None,
                'source': str(path),
            }
    return None


def main():
    parser = argparse.ArgumentParser(description='Web page with the distributions and sensitivities of a tag')
    parser.add_argument('--tag', required=True, help='Tag directory (run_chain.py output)')
    parser.add_argument('--cutcount-dir', default=None, help='Default: <tag>/cutcount')
    parser.add_argument('--out', default=str(Path(__file__).resolve().parent / 'index.html'))
    parser.add_argument('--title', default='Run-1A sensitivity')
    args = parser.parse_args()

    tag_dir = Path(args.tag).expanduser().resolve()
    cutcount_dir = Path(args.cutcount_dir).expanduser().resolve() if args.cutcount_dir else tag_dir / 'cutcount'
    if not (tag_dir / 'manifest.yaml').exists():
        raise SystemExit(f'No manifest.yaml in {tag_dir}: run the selection for this tag first '
                         '(is PYFITTER_DIR set? source sensitivity/setup.sh)')
    with open(tag_dir / 'manifest.yaml') as f:
        manifest = yaml.safe_load(f)
    card_path = cutcount_dir / 'count_card.yaml'
    if not card_path.exists():
        raise SystemExit(f'No counting card at {card_path}: run cutcount/run_cutcount.py first')
    with open(card_path) as f:
        card = yaml.safe_load(f)
    p_window, t_window = card['window']['mom'], card['window']['time']

    mom_range = (float(manifest['fitrange_low'][0]), float(manifest['fitrange_hi'][0]))
    time_range = (float(manifest['fitrange_low'][1]), float(manifest['fitrange_hi'][1]))
    mom_edges = np.arange(mom_range[0], mom_range[1] + MOM_BIN / 2, MOM_BIN)
    time_edges = np.arange(time_range[0], time_range[1] + TIME_BIN / 2, TIME_BIN)

    data = {
        'title': args.title,
        'tag': manifest.get('tag', tag_dir.name),
        'version': card.get('metadata', {}).get('version'),
        'generated': datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
        'window': {'mom': p_window, 'time': t_window},
        'range': {'mom': mom_range, 'time': time_range},
        'edges': {'mom': [float(v) for v in mom_edges], 'time': [float(v) for v in time_edges]},
        'processes': histograms(manifest, p_window, t_window, mom_edges, time_edges),
        'count': count_limit(cutcount_dir),
        'shape': shape_limit(tag_dir),
        'sources': {'manifest': str(tag_dir / 'manifest.yaml'), 'card': str(card_path)},
    }
    payload = json.dumps(data).replace('</', '<\\/')
    html = TEMPLATE.replace('__TITLE__', args.title).replace('__DATA__', payload)
    Path(args.out).write_text(html)
    print(f'Wrote {args.out}')


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7;
  --surface-1: #fcfcfb;
  --text-primary: #0b0b0b;
  --text-secondary: #52514e;
  --text-muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --border: rgba(11,11,11,0.10);
  --window: rgba(11,11,11,0.05);
  --s-CE: #2a78d6; --s-DIO: #eb6834; --s-Cosmic: #1baf7a; --s-RPC: #eda100; --s-RMC: #e87ba4;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface-1: #1a1a19;
    --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --window: rgba(255,255,255,0.06);
    --s-CE: #3987e5; --s-DIO: #d95926; --s-Cosmic: #199e70; --s-RPC: #c98500; --s-RMC: #d55181;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface-1: #1a1a19;
  --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --window: rgba(255,255,255,0.06);
  --s-CE: #3987e5; --s-DIO: #d95926; --s-Cosmic: #199e70; --s-RPC: #c98500; --s-RMC: #d55181;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--text-primary);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1080px; margin: 0 auto; padding: 32px 16px 48px; }
h1 { font-size: 24px; font-weight: 600; margin: 0 0 4px; }
h2 { font-size: 16px; font-weight: 600; margin: 0 0 2px; }
.sub { color: var(--text-secondary); margin: 0 0 24px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; margin-bottom: 12px; }
.tile, .panel { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px; padding: 20px; }
.tile .label { color: var(--text-secondary); font-size: 14px; }
.tile .value { font-size: 36px; font-weight: 600; margin: 4px 0; }
.tile .detail { color: var(--text-secondary); font-size: 14px; }
.tile .method { color: var(--text-muted); font-size: 13px; margin-top: 8px; }
.compare { color: var(--text-secondary); margin: 0 0 24px; }
.controls { display: flex; flex-wrap: wrap; gap: 12px 24px; align-items: center; margin-bottom: 16px; }
.seg { display: inline-flex; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
.seg button { font: inherit; font-size: 14px; background: var(--surface-1); color: var(--text-secondary);
  border: 0; padding: 6px 12px; cursor: pointer; }
.seg button[aria-pressed="true"] { background: var(--text-primary); color: var(--surface-1); }
.legend { display: flex; flex-wrap: wrap; gap: 4px 16px; font-size: 14px; color: var(--text-secondary); }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
.key { width: 16px; height: 2px; border-radius: 1px; display: inline-block; }
.windowkey { width: 14px; height: 12px; background: var(--window); border: 1px solid var(--axis); display: inline-block; }
.charts { display: grid; gap: 16px; margin-bottom: 24px; }
.panel .note { color: var(--text-secondary); font-size: 14px; margin: 0 0 8px; }
.chart { position: relative; }
.chart svg { display: block; width: 100%; height: 320px; }
.tooltip { position: absolute; pointer-events: none; background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 8px; padding: 8px 10px; font-size: 13px; box-shadow: 0 2px 8px rgba(0,0,0,0.12);
  display: none; min-width: 160px; z-index: 2; }
.tooltip .row { display: flex; justify-content: space-between; gap: 12px; font-variant-numeric: tabular-nums; }
.tooltip .row span:first-child { display: inline-flex; align-items: center; gap: 6px; color: var(--text-secondary); }
.tooltip .head { color: var(--text-primary); font-weight: 600; margin-bottom: 4px; }
table { width: 100%; border-collapse: collapse; font-size: 14px; }
th, td { text-align: right; padding: 8px 10px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
th { color: var(--text-secondary); font-weight: 600; }
tfoot td { font-weight: 600; border-bottom: 0; }
.tablewrap { overflow-x: auto; }
footer { color: var(--text-muted); font-size: 13px; margin-top: 24px; word-break: break-all; }
</style>
</head>
<body>
<main>
  <h1 id="title"></h1>
  <p class="sub" id="subtitle"></p>

  <div class="tiles" id="tiles"></div>
  <p class="compare" id="compare"></p>

  <div class="controls">
    <div class="seg" role="group" aria-label="Y scale">
      <button data-scale="log" aria-pressed="true">Log</button>
      <button data-scale="linear" aria-pressed="false">Linear</button>
    </div>
    <div class="seg" role="group" aria-label="Tracks shown">
      <button data-sel="all" aria-pressed="true">All selected tracks</button>
      <button data-sel="win" aria-pressed="false">Inside the other window</button>
    </div>
    <div class="legend" id="legend"></div>
  </div>

  <div class="charts">
    <section class="panel">
      <h2>Reconstructed momentum</h2>
      <p class="note" id="mom-note"></p>
      <div class="chart" id="chart-mom"></div>
    </section>
    <section class="panel">
      <h2>Reconstructed time</h2>
      <p class="note" id="time-note"></p>
      <div class="chart" id="chart-time"></div>
    </section>
  </div>

  <section class="panel">
    <h2>Expected yields</h2>
    <p class="note">Events expected in the run. MC tracks is the raw number of simulated tracks in the counting window.</p>
    <div class="tablewrap"><table id="yields"></table></div>
  </section>

  <footer id="footer"></footer>
</main>

<script>
const DATA = __DATA__;
const NS = 'http://www.w3.org/2000/svg';
const state = { scale: 'log', sel: 'all' };
const procs = Object.keys(DATA.processes);

function sci(x, digits = 2) {
  if (x == null || !isFinite(x)) return '–';
  const e = Math.floor(Math.log10(Math.abs(x)));
  return `${(x / 10 ** e).toFixed(digits)} × 10<sup>${e}</sup>`;
}
function num(x) {
  if (x === 0) return '0';
  if (Math.abs(x) >= 100) return x.toFixed(0);
  if (Math.abs(x) >= 1) return x.toFixed(2);
  if (Math.abs(x) >= 0.001) return x.toPrecision(2);
  return x.toExponential(1);
}
function decade(x) {
  const e = Math.round(Math.log10(x));
  if (e >= 0 && e <= 3) return String(10 ** e);
  const sup = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  return '10' + String(e).split('').map(c => sup[c]).join('');
}
function el(tag, attrs = {}, parent) {
  const node = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (parent) parent.appendChild(node);
  return node;
}
const color = p => `var(--s-${p}, var(--text-muted))`;

// ---------- header, tiles, legend, table ----------
document.getElementById('title').textContent = DATA.title;
const w = DATA.window;
document.getElementById('subtitle').textContent =
  `Tag ${DATA.tag}${DATA.version ? ', cut set ' + DATA.version : ''}. ` +
  `Counting window ${w.mom[0]}–${w.mom[1]} MeV/c and ${w.time[0]}–${w.time[1]} ns. Expected 90% CL upper limits.`;

function tile(name, res) {
  const div = document.createElement('div');
  div.className = 'tile';
  if (!res) {
    div.innerHTML = `<div class="label">${name}</div><div class="value">Not run</div>`;
    return div;
  }
  div.innerHTML = `<div class="label">${name}: expected R<sub>μe</sub> limit</div>
    <div class="value">${sci(res.rmue)}</div>
    <div class="detail">${res.events.toFixed(2)} signal events · 68% band ${sci(res.rmue_lo)} to ${sci(res.rmue_hi)}</div>
    <div class="detail">CE efficiency ${(100 * res.efficiency).toFixed(1)}%${res.background != null ? ` · background ${res.background.toFixed(2)} events in the window` : ''}</div>
    <div class="method">${res.method}</div>`;
  return div;
}
const tiles = document.getElementById('tiles');
tiles.appendChild(tile('Cut and count', DATA.count));
tiles.appendChild(tile('Shape fit', DATA.shape));
if (DATA.count && DATA.shape) {
  const r = DATA.count.rmue / DATA.shape.rmue;
  document.getElementById('compare').textContent = r >= 1
    ? `The shape fit limit is ${r.toFixed(2)}× stronger than the cut and count.`
    : `The cut and count limit is ${(1 / r).toFixed(2)}× stronger than the shape fit.`;
}

const legend = document.getElementById('legend');
for (const p of procs) legend.insertAdjacentHTML('beforeend', `<span><i class="key" style="background:${color(p)}"></i>${p}</span>`);
legend.insertAdjacentHTML('beforeend', `<span><i class="windowkey"></i>Counting window</span>`);

const table = document.getElementById('yields');
let rows = '<thead><tr><th>Process</th><th>Counting window</th><th>Fit range</th><th>MC tracks in window</th></tr></thead><tbody>';
let bW = 0, bF = 0;
for (const p of procs) {
  const h = DATA.processes[p];
  if (p !== 'CE') { bW += h.window_yield; bF += h.fit_yield; }
  rows += `<tr><td><span class="legend"><span><i class="key" style="background:${color(p)}"></i>${p}</span></span></td>
    <td>${num(h.window_yield)}</td><td>${num(h.fit_yield)}</td><td>${h.window_tracks.toLocaleString()}</td></tr>`;
}
rows += `</tbody><tfoot><tr><td>Total background</td><td>${num(bW)}</td><td>${num(bF)}</td><td></td></tr></tfoot>`;
table.innerHTML = rows;

const r = DATA.range;
document.getElementById('footer').innerHTML =
  `Generated ${DATA.generated} by make_page.py. Selection range ${r.mom[0]}–${r.mom[1]} MeV/c, ${r.time[0]}–${r.time[1]} ns.<br>` +
  `Manifest: ${DATA.sources.manifest}<br>Counting card: ${DATA.sources.card}` +
  (DATA.count ? `<br>Cut-and-count limit: ${DATA.count.source}` : '') +
  (DATA.shape ? `<br>Shape-fit limit: ${DATA.shape.source}` : '');

// ---------- charts ----------
function niceStep(range, target) {
  const raw = range / target, mag = 10 ** Math.floor(Math.log10(raw)), f = raw / mag;
  return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag;
}
function linearTicks(lo, hi, target) {
  const step = niceStep(hi - lo, target), out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}

function draw(key) {
  const host = document.getElementById(`chart-${key}`);
  host.innerHTML = '';
  const edges = DATA.edges[key], win = DATA.window[key];
  const field = `${key}_${state.sel}`;
  const series = procs.map(p => ({ p, v: DATA.processes[p][field] }));
  const W = Math.max(host.clientWidth, 280), H = 320;
  const m = { l: 64, r: 12, t: 12, b: 40 };
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img',
    'aria-label': `${key === 'mom' ? 'Momentum' : 'Time'} distribution per process` }, host);

  const x0 = edges[0], x1 = edges[edges.length - 1];
  const X = v => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
  const all = series.flatMap(s => s.v), vmax = Math.max(...all, 1e-9);
  let Y, yTicks, yLo;
  if (state.scale === 'log') {
    const pos = all.filter(v => v > 0);
    const lo = Math.max(Math.min(...pos), vmax * 1e-5);
    const dLo = Math.floor(Math.log10(lo)), dHi = Math.ceil(Math.log10(vmax * 1.5));
    yLo = 10 ** dLo;
    Y = v => H - m.b - (Math.log10(v) - dLo) / (dHi - dLo) * (H - m.t - m.b);
    yTicks = [];
    const every = Math.ceil((dHi - dLo + 1) / 7);
    for (let d = dLo; d <= dHi; d += every) yTicks.push(10 ** d);
  } else {
    yTicks = linearTicks(0, vmax * 1.1, 5);
    const top = Math.max(yTicks[yTicks.length - 1], vmax * 1.1);
    yLo = 0;
    Y = v => H - m.b - v / top * (H - m.t - m.b);
  }

  // window band, behind everything
  el('rect', { x: X(win[0]), y: m.t, width: X(win[1]) - X(win[0]), height: H - m.t - m.b, fill: 'var(--window)' }, svg);

  // grid and y axis
  for (const t of yTicks) {
    const y = Y(t);
    el('line', { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: 'var(--grid)', 'stroke-width': 1 }, svg);
    const lab = el('text', { x: m.l - 8, y: y + 4, 'text-anchor': 'end', fill: 'var(--text-muted)', 'font-size': 12 }, svg);
    lab.textContent = state.scale === 'log' ? decade(t) : num(t);
  }
  el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: 'var(--axis)', 'stroke-width': 1 }, svg);
  for (const t of linearTicks(x0, x1, W < 500 ? 5 : 10)) {
    const lab = el('text', { x: X(t), y: H - m.b + 18, 'text-anchor': 'middle', fill: 'var(--text-muted)', 'font-size': 12 }, svg);
    lab.textContent = t;
  }
  const xl = el('text', { x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', fill: 'var(--text-secondary)', 'font-size': 13 }, svg);
  xl.textContent = key === 'mom' ? 'Momentum at tracker front [MeV/c]' : 'Time at tracker front [ns]';
  const binw = edges[1] - edges[0];
  const yl = el('text', { transform: `translate(14 ${(m.t + H - m.b) / 2}) rotate(-90)`, 'text-anchor': 'middle',
    fill: 'var(--text-secondary)', 'font-size': 13 }, svg);
  yl.textContent = `Events / ${+binw.toFixed(3)} ${key === 'mom' ? 'MeV/c' : 'ns'}`;

  // step lines; on log scale bins below the axis floor break the line
  for (const s of series) {
    let d = '', open = false;
    s.v.forEach((v, i) => {
      const ok = state.scale === 'linear' || v >= yLo;
      if (!ok) { open = false; return; }
      const y = Y(v);
      d += open ? `V${y}H${X(edges[i + 1])}` : `M${X(edges[i])} ${y}H${X(edges[i + 1])}`;
      open = true;
    });
    el('path', { d, fill: 'none', stroke: color(s.p), 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }, svg);
  }

  // hover: bin highlight and tooltip
  const hl = el('rect', { y: m.t, height: H - m.t - m.b, fill: 'var(--text-primary)', opacity: 0.06, visibility: 'hidden' }, svg);
  const tip = document.createElement('div');
  tip.className = 'tooltip';
  host.appendChild(tip);
  const hit = el('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent' }, svg);
  const unit = key === 'mom' ? 'MeV/c' : 'ns';
  function show(evt) {
    const box = svg.getBoundingClientRect();
    const px = (evt.clientX - box.left) * W / box.width;
    const i = Math.min(edges.length - 2, Math.max(0, Math.floor((px - m.l) / (W - m.l - m.r) * (edges.length - 1))));
    hl.setAttribute('x', X(edges[i])); hl.setAttribute('width', Math.max(1, X(edges[i + 1]) - X(edges[i])));
    hl.setAttribute('visibility', 'visible');
    let html = `<div class="head">${+edges[i].toFixed(2)}–${+edges[i + 1].toFixed(2)} ${unit}</div>`;
    let bkg = 0;
    for (const s of series) {
      if (s.p !== 'CE') bkg += s.v[i];
      html += `<div class="row"><span><i class="key" style="background:${color(s.p)}"></i>${s.p}</span><span>${num(s.v[i])}</span></div>`;
    }
    html += `<div class="row"><span>Background</span><span>${num(bkg)}</span></div>`;
    tip.innerHTML = html;
    tip.style.display = 'block';
    const left = (evt.clientX - box.left) + 16;
    tip.style.left = `${left + tip.offsetWidth > box.width ? left - tip.offsetWidth - 32 : left}px`;
    tip.style.top = '16px';
  }
  hit.addEventListener('pointermove', show);
  hit.addEventListener('pointerleave', () => { tip.style.display = 'none'; hl.setAttribute('visibility', 'hidden'); });
}

function notes() {
  const inner = state.sel === 'win';
  document.getElementById('mom-note').textContent = inner
    ? `Tracks with ${w.time[0]} < t < ${w.time[1]} ns. Shaded: momentum window.`
    : 'All tracks passing the selection. Shaded: momentum window.';
  document.getElementById('time-note').textContent = inner
    ? `Tracks with ${w.mom[0]} < p < ${w.mom[1]} MeV/c. Shaded: time window.`
    : 'All tracks passing the selection. Shaded: time window.';
}
function render() { notes(); draw('mom'); draw('time'); }

for (const b of document.querySelectorAll('[data-scale]')) b.addEventListener('click', () => {
  state.scale = b.dataset.scale;
  document.querySelectorAll('[data-scale]').forEach(o => o.setAttribute('aria-pressed', o === b));
  render();
});
for (const b of document.querySelectorAll('[data-sel]')) b.addEventListener('click', () => {
  state.sel = b.dataset.sel;
  document.querySelectorAll('[data-sel]').forEach(o => o.setAttribute('aria-pressed', o === b));
  render();
});
let resizeTimer;
window.addEventListener('resize', () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(render, 100); });
render();
</script>
</body>
</html>
"""


if __name__ == '__main__':
    main()

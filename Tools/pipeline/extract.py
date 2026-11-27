#!/usr/bin/env python3
"""Captures -> models: a folder of tape captures in, one .ttm model (and a QC report) per tape and speed out.

  python3 Tools/pipeline/extract.py <captures folder or files...> [-o Resources/Models] [--jobs N]
                                    [--channel L|R] [--dry-run] [--only ATR_15]

Every capture is the Technical Tape test file (Tools/capture/make_capture.py, 333 s) recorded off the
repro head. A tape at a speed needs a normal pass and, ideally, a +6 dB pass (see docs/CAPTURES.md).
File names are read loosely: the tape (ATR, Maxell, or any word you use consistently), the speed
(15, 7.5 / 7,5 / 7½ / 7_5, 3.75 / 3¾ / 3_75, optionally followed by "ips") and the pass ("+6", "6dB",
"hot" or "plus6" for the hot pass; anything else is the normal pass). Examples:
  ATR_15ips_0dB.wav   Maxell 7.5 +6dB.wav   atr-3.75-hot.wav   MAXELL_7½ips_normal.wav
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'capture'))

SPEEDS = (('3.75', r'3\s*[.,_]\s*75|3\s*¾|3\s+3\s*/\s*4|375(?=\s*ips)'),
          ('7.5', r'7\s*[.,_]\s*5|7\s*½|7\s+1\s*/\s*2|75(?=\s*ips)'),
          ('15', r'15'))
HOT = r'(\+\s*6|plus\s*6|\b6\s*db|hot|_6db|\+6db)'
KNOWN_TAPES = {'atr': 'ATR', 'maxell': 'Maxell'}


def parse_name(path):
    """(tape, speed, pass) from a file name; pass is 'normal' or 'hot'. None for what can't be read."""
    stem = os.path.splitext(os.path.basename(path))[0]
    low = stem.lower().replace('ips', ' ips')
    speed = None
    for name, pat in SPEEDS:
        if re.search(r'(?<![\d.])(' + pat + r')(?![\d])', low):
            speed = name
            break
    hot = re.search(HOT, low) is not None
    tape = None
    for k, v in KNOWN_TAPES.items():
        if k in low:
            tape = v
            break
    if tape is None:
        words = [w for w in re.split(r'[\s_\-.,+]+', stem) if re.match(r'^[A-Za-z][A-Za-z0-9]*$', w)
                 and w.lower() not in ('ips', 'db', 'hot', 'normal', 'regular', 'nominal', 'plus', 'pass', 'take')]
        tape = words[0] if words else None
    return tape, speed, ('hot' if hot else 'normal')


def collect(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            for root, _, names in os.walk(p):
                files += [os.path.join(root, n) for n in sorted(names) if n.lower().endswith(('.wav', '.wave'))]
        else:
            files.append(p)
    groups, skipped = {}, []
    for f in files:
        tape, speed, pas = parse_name(f)
        if tape is None or speed is None:
            skipped.append(f)
            continue
        g = groups.setdefault((tape, speed), {})
        if pas in g:
            older, newer = sorted([g[pas], f], key=os.path.getmtime)
            print(f'  note: two {pas} passes for {tape} {speed} ips; using the newer one, {os.path.basename(newer)}')
            g[pas] = newer
        else:
            g[pas] = f
    return groups, skipped


def default_jobs():
    """One build per ~5 GB of memory and per two cores (a build peaks at ~3.5 GB)."""
    try:
        mem_gb = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / 2 ** 30
    except (ValueError, OSError, AttributeError):
        mem_gb = 8.0
    return max(1, min((os.cpu_count() or 2) // 2, int(mem_gb // 5)))


def run_one(tape, speed, passes, out_dir, channel, simulated=False):
    """Analyse both passes and build the model (runs in a worker process). Returns a summary dict."""
    import analyse
    import build
    import report
    lines = []

    def log(s=''):
        lines.append(s)

    t0 = time.time()
    stem = f'{tape}_{speed}'
    try:
        log(f'{tape} {speed} ips')
        r0 = analyse.analyse(passes['normal'], channel, log=log)
        r6 = analyse.analyse(passes['hot'], channel, log=log) if 'hot' in passes else None
        if r6 is None:
            log('  no +6 dB pass: the model stops at the top of the normal pass (+12 VU), and is continued above by '
                'the fitted saturation curve')
        out = os.path.join(out_dir, f'{stem}.ttm')
        model, dbg = build.build(tape, speed, r0, r6, out, log=log,
                                 extra_meta={'simulated': True} if simulated else None)
        rep = os.path.join(out_dir, f'{stem}_report.html')
        try:
            report.write_report(rep, model, dbg, r0, r6, lines)
        except Exception:  # noqa: BLE001 - the model is good even if a plot isn't
            lines.append('  QC report failed:\n' + traceback.format_exc())
            rep = None
        return dict(ok=True, name=stem, model=out, report=rep, seconds=time.time() - t0, log=lines,
                    warnings=model.meta.get('warnings', []))
    except Exception as e:  # noqa: BLE001 - one bad capture must not stop the others
        lines.append(traceback.format_exc())
        return dict(ok=False, name=stem, error=str(e), log=lines, seconds=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs', nargs='+', help='capture WAVs or folders of them')
    ap.add_argument('-o', '--out', default=os.path.join(HERE, '..', '..', 'Resources', 'Models'))
    ap.add_argument('--jobs', type=int, default=default_jobs(),
                    help='models built at once (each needs ~4 GB of memory; default from this machine)')
    ap.add_argument('--channel', choices=['L', 'R'], help='use this channel of stereo files (default: the louder)')
    ap.add_argument('--only', help='build just these, e.g. ATR_15,Maxell_7.5')
    ap.add_argument('--dry-run', action='store_true', help='show how the files were read, build nothing')
    ap.add_argument('--simulated', action='store_true', help=argparse.SUPPRESS)   # captures of Tools/sim's virtual machine
    a = ap.parse_args()

    groups, skipped = collect(a.inputs)
    for f in skipped:
        print(f'  skipped (no tape or speed in the name): {f}')
    if a.only:
        want = set(a.only.split(','))
        groups = {k: v for k, v in groups.items() if f'{k[0]}_{k[1]}' in want}
    if not groups:
        print('no captures found')
        return 1
    print('captures:')
    for (tape, speed), g in sorted(groups.items()):
        n = os.path.basename(g.get('normal', '-- missing --'))
        h = os.path.basename(g.get('hot', '-- missing --'))
        print(f'  {tape:8s} {speed:>4s} ips   normal: {n:40s} +6 dB: {h}')
    if a.dry_run:
        return 0
    os.makedirs(a.out, exist_ok=True)
    channel = {'L': 0, 'R': 1}.get(a.channel)
    todo = [(t, s, g) for (t, s), g in sorted(groups.items()) if 'normal' in g]
    for (t, s), g in groups.items():
        if 'normal' not in g:
            print(f'  {t} {s} ips: only a +6 dB pass; a normal pass is needed. Skipped.')
    results = []

    def report_one(r):
        results.append(r)
        print('\n'.join(r['log']))
        print(f'  -> {"built " + r["model"] if r["ok"] else "FAILED: " + r["error"]} ({r["seconds"]:.0f} s)\n', flush=True)

    if a.jobs <= 1:
        for t, s_, g in todo:
            report_one(run_one(t, s_, g, a.out, channel, a.simulated))
    else:
        done = set()
        try:
            with ProcessPoolExecutor(max_workers=a.jobs) as ex:
                futs = {ex.submit(run_one, t, s_, g, a.out, channel, a.simulated): (t, s_) for t, s_, g in todo}
                for fu in as_completed(futs):
                    report_one(fu.result())
                    done.add(futs[fu])
        except BrokenProcessPool:
            # a worker died (usually out of memory): build what's left one at a time
            print('  a build process stopped (out of memory?); building the rest one at a time', flush=True)
            for t, s_, g in todo:
                if (t, s_) not in done:
                    report_one(run_one(t, s_, g, a.out, channel, a.simulated))
    summary = {r['name']: {k: v for k, v in r.items() if k != 'log'} for r in results}
    with open(os.path.join(a.out, 'extract_summary.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    bad = [r['name'] for r in results if not r['ok']]
    print(f'{len(results) - len(bad)} model(s) built in {a.out}' + (f'; failed: {", ".join(bad)}' if bad else ''))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())

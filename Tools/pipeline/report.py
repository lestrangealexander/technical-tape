"""QC report for a built model: one self-contained HTML page (plots embedded) next to the .ttm.

What to look at: the calibration and pass offsets (did the capture go as planned), the warnings, the
frequency response and compression curves (do they look like tape), the harmonic curves (smooth, no
noise-floor plateaus), the extrapolation above the hottest capture, the wow & flutter spectrum.
"""
from __future__ import annotations

import base64
import html
import io
import math

import numpy as np

import analyse as an

FK = an.FK
INK = '#1f2328'
MUTED = '#6e7781'
SERIES = ['#2f6fb3', '#c4553a', '#3d8c5a', '#8a5cc2', '#b8860b', '#1b998b', '#c2185b']


def _fig(w=7.2, h=3.4):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#c9d1d9', 'axes.labelcolor': INK, 'xtick.color': MUTED,
                         'ytick.color': MUTED, 'axes.grid': True, 'grid.color': '#eaeef2', 'grid.linewidth': 0.8,
                         'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
    fig, ax = plt.subplots(figsize=(w, h), dpi=110)
    return fig, ax


def _png(fig):
    import matplotlib.pyplot as plt
    b = io.BytesIO()
    fig.tight_layout()
    fig.savefig(b, format='png')
    plt.close(fig)
    return base64.b64encode(b.getvalue()).decode('ascii')


def _band_db(H, f, frac=1 / 6):
    sel = (FK >= f * 2 ** (-frac / 2)) & (FK <= f * 2 ** (frac / 2))
    return 10 * np.log10(np.mean(np.abs(H[..., sel]) ** 2, axis=-1) + 1e-30)


def plot_response(dbg):
    levels, Hg = dbg['levels'], dbg['Hg']
    fig, ax = _fig()
    f = np.geomspace(20, 21000, 400)
    ref = _band_db(Hg[np.argmin(np.abs(levels + 20)), 0], 1000.0)
    for i, vu in enumerate((-20, -10, 0, 6, 12, 18)):
        j = int(np.argmin(np.abs(levels - vu)))
        if abs(levels[j] - vu) > 1.6:
            continue
        y = np.array([_band_db(Hg[j, 0], x, 1 / 12) for x in f]) - ref
        ax.semilogx(f, y, color=SERIES[i % len(SERIES)], lw=1.4, label=f'{levels[j]:+.0f} VU')
    ax.set_xlim(20, 21000)
    ax.set_ylim(-24, 6)
    ax.set_xlabel('Hz')
    ax.set_ylabel('dB re 1 kHz at -20 VU')
    ax.set_title('Frequency response (fundamental) by input level', loc='left', fontsize=10, color=INK)
    ax.legend(ncol=6, fontsize=8, loc='lower left')
    return _png(fig)


def plot_compression(dbg):
    levels, Hg = dbg['levels'], dbg['Hg']
    lx, Hx = dbg.get('levels_x', levels), dbg.get('Hg_x', Hg)
    fig, ax = _fig()
    for i, fr in enumerate((100.0, 1000.0, 5000.0, 10000.0, 15000.0)):
        g = np.array([_band_db(Hx[j, 0], fr) for j in range(len(lx))])
        g -= g[0]
        n = len(levels)
        c = SERIES[i % len(SERIES)]
        ax.plot(lx[:n], g[:n], color=c, lw=1.6, marker='o', ms=2.5, label=f'{fr / 1000:g} kHz' if fr >= 1000 else f'{fr:.0f} Hz')
        ax.plot(lx[n - 1:], g[n - 1:], color=c, lw=1.2, ls='--')
    ax.axvline(levels[-1], color=MUTED, lw=0.8, ls=':')
    ax.text(levels[-1] + 0.5, 0.5, 'hottest capture', color=MUTED, fontsize=8)
    ax.set_xlabel('input level (VU)')
    ax.set_ylabel('gain change (dB)')
    ax.set_ylim(-30, 2)
    ax.set_title('Compression: fundamental gain vs level (dashed: fitted continuation)', loc='left', fontsize=10, color=INK)
    ax.legend(ncol=5, fontsize=8, loc='lower left')
    return _png(fig)


def plot_harmonics(dbg):
    levels, Hg, H = dbg['levels'], dbg['Hg'], dbg['H']
    fig, ax = _fig()
    f0 = 1000.0
    h1 = np.array([_band_db(Hg[j, 0], f0) for j in range(len(levels))])
    for i, m in enumerate((2, 3, 4, 5, 7, 9)):
        if m > Hg.shape[1]:
            continue
        hm = np.array([_band_db(Hg[j, m - 1], m * f0) for j in range(len(levels))])
        raw = np.array([_band_db(H[j, m - 1], m * f0) for j in range(len(levels))])
        c = SERIES[i % len(SERIES)]
        ax.plot(levels, hm - h1, color=c, lw=1.5, label=f'H{m}')
        ax.plot(levels, raw - h1, color=c, lw=0, marker='.', ms=3, alpha=0.6)
    ax.set_xlabel('input level (VU)')
    ax.set_ylabel('dB re fundamental')
    ax.set_ylim(-110, 0)
    ax.set_title('Harmonics of 1 kHz (dots: measured, line: as used after noise gating)', loc='left', fontsize=10,
                 color=INK)
    ax.legend(ncol=6, fontsize=8, loc='lower right')
    return _png(fig)


def plot_preemphasis(model):
    fig, ax = _fig(7.2, 2.6)
    ax.semilogx(model.pemp_freqs, model.pemp_db, color=SERIES[0], lw=1.6)
    ax.set_xlim(30, 17000)
    ax.set_xlabel('Hz')
    ax.set_ylabel('dB')
    ax.set_title('Record pre-emphasis estimate |P(f)| (from where compression starts, re 1 kHz)', loc='left',
                 fontsize=10, color=INK)
    return _png(fig)


def plot_flutter(model):
    fig, ax = _fig(7.2, 2.8)
    c = np.asarray(model.flutter_cont, float).reshape(-1, 2)
    if len(c):
        dev = (2 * np.pi * c[:, 0]) ** 2 * c[:, 1]                # speed-deviation PSD (1/Hz)
        ax.loglog(c[:, 0], np.sqrt(dev * c[:, 0] * 0.23) * 100, color=SERIES[0], lw=1.4, label='continuum (1/3 oct)')
    for f, a in np.asarray(model.flutter_peaks, float).reshape(-1, 2):
        ax.loglog([f, f], [1e-5, 2 * np.pi * f * a * 100], color=SERIES[1], lw=1.4)
    ax.set_xlim(0.2, 200)
    ax.set_ylim(1e-5, 1)
    ax.set_xlabel('Hz')
    ax.set_ylabel('speed deviation, % rms')
    ax.set_title('Wow & flutter as modelled (lines: discrete components)', loc='left', fontsize=10, color=INK)
    ax.legend(fontsize=8, loc='upper right')
    return _png(fig)


def plot_autogain(model):
    fig, ax = _fig(7.2, 2.4)
    g = np.asarray(model.ag_gain_db, float)
    d = model.ag_drive0 + model.ag_step * np.arange(len(g))
    ax.plot(d, g, color=SERIES[2], lw=1.6)
    ax.set_xlabel('INPUT (dB)')
    ax.set_ylabel('AG make-up (dB)')
    ax.set_title('Auto gain: loudness-matched make-up for each INPUT setting', loc='left', fontsize=10, color=INK)
    return _png(fig)


def plot_music(model):
    """The music fit: model minus machine per 1/3 octave on each music pass (total power, the session's drift
    out; bands lost in the machine's noise left out)."""
    mu = (model.meta or {}).get('music')
    if not mu or not mu.get('band_errors_db'):
        return None
    cent = 1000.0 * 2 ** (np.arange(-13, 13) / 3.0)
    fig, ax = _fig(7.2, 2.8)
    for k, (lv, e) in enumerate(zip(mu['passes_db'], mu['band_errors_db'])):
        ev = np.array([np.nan if v is None else v for v in e], float)
        ax.semilogx(cent, ev, color=SERIES[k % len(SERIES)], lw=1.6, marker='o', ms=2.5, label=f'{lv:+.1f} dB')
    ax.axhline(0.0, color=MUTED, lw=0.8)
    ax.set_xlim(50, 20000)
    ax.set_ylim(-6, 6)
    ax.set_xlabel('Hz')
    ax.set_ylabel('model - machine (dB)')
    ax.legend(loc='lower left', fontsize=8, ncol=3, title='music onto the tape', title_fontsize=8)
    ax.set_title('Music through the model against the machine, per 1/3 octave', loc='left', fontsize=10, color=INK)
    return _png(fig)


def _music_summary(mu):
    if not mu:
        return 'no music passes: the treble stage was fitted on the pink noise'
    rows = [f'{lv:+.1f} dB: {e:.2f} dB rms, above 1.8 kHz {t:+.2f} dB, top end {tm[0]:+.1f} / {tm[1]:+.1f}, presence '
            f'{pr[0]:+.1f} / {pr[1]:+.1f} dB (machine / model)'
            for lv, e, t, tm, pr in zip(mu['passes_db'], mu['error_per_pass_db'], mu['treble_error_per_pass_db'],
                                        mu['top_machine_model_db'], mu['presence_machine_model_db'])]
    return html.escape(f'{mu["error_db"]:.2f} dB rms over 1/3 octaves 63 Hz-16 kHz and passes, total power'
                       + (f' ({mu["outside_machine_db"]:.2f} dB outside the machine\'s: what the fit minimises)'
                          if mu.get('outside_machine_db') is not None else '') + '; ') + '<br>' + \
        '<br>'.join(html.escape(r) for r in rows)


def _hot_summary(p):
    if 'hot_in_db' not in p:
        return html.escape(str(p))
    fit = p.get('hot_fit', {})
    return html.escape(f'onto the tape {p["hot_in_db"]:+.2f} dB hotter than the normal pass, played back '
                       f'{p["hot_out_db"]:+.2f} dB, {p.get("hot_tau_us", 0.0):+.1f} us; odd harmonics line up to '
                       f'{fit.get("odd_rms_db")} dB rms and compression to {fit.get("compression_rms_db")} dB over '
                       f'{fit.get("points")} points (2nd harmonic: {fit.get("h2_rms_db")} dB, it follows the playback '
                       f'level). Used above the normal pass\'s top level.')


def write_report(path, model, dbg, r0, r6, log_lines=()):
    m = model.meta
    esc = html.escape
    rows = [('Tape / speed', f'{m.get("tape")} {m.get("speed_ips")} ips'),
            ('Built', m.get('built', '')),
            ('Captures', '<br>'.join(esc(str(s)) for s in m.get('sources', []))),
            ('Calibration tone', f'{r0.cal_dbfs:+.2f} dBFS RMS at {r0.cal_freq:.2f} Hz, THD {r0.cal_thd_db:.1f} dB '
                                 f'(0 VU is expected at -18.00 dBFS)'),
            ('Noise floor', f'{r0.noise_floor_dbfs:.1f} dBFS'),
            ('Levels captured', f'{model.captured_low_vu:+.1f} to {model.captured_high_vu:+.1f} VU'),
            ('Harmonic orders', str(m.get('orders'))),
            ('Branches', esc(', '.join(f'{i + 1}: {s}' for i, s in enumerate(m.get('branch_modes', []))))),
            ('Kernel truncation (worst)', f'{max(m.get("truncation_db", [-999])):.1f} dB re fundamental'),
            ('Pre-emphasis fit', esc(str(m.get('preemphasis')))),
            ('Treble stage', ('; '.join((f'{tb.fc_lo / 1000:g}-{tb.fc_hi / 1000:g} kHz' if tb.fc_hi > 0 else
                                         f'from {tb.fc_lo / 1000:g} kHz') for tb in model.treble_bands)
                              + ': each band saturates on its own level, read by '
                              + (f'its energy over {model.treble_bands[0].rms_ms:g} ms'
                                 + (f' with {model.treble_bands[0].peak_mix:g} of its peaks mixed in'
                                    if getattr(model.treble_bands[0], 'peak_mix', 0.0) > 0 else '')
                                 if getattr(model.treble_bands[0], 'rms_ms', 0.0) > 0 else 'its peaks')
                              + (f', its detector also hearing {model.treble_bands[0].leak:g} of the whole signal'
                                 if getattr(model.treble_bands[0], 'leak', 0.0) > 0 else '')
                              if getattr(model, 'treble_bands', None)
                              else 'none')),
            ('Music fit', _music_summary(m.get('music'))),
            ('Pink noise' + (' (a check: the music decided)' if m.get('music') else ' (broadband fit)'),
             esc(str(m.get('broadband')))),
            ('Above the hottest capture', esc(str(m.get('extrapolation')))),
            ('Wow & flutter (capture)', esc(str(m.get('wow_flutter')))),
            ('Hot pass', _hot_summary(m.get('passes', {})) if r6 is not None else 'none')]
    warn = m.get('warnings', [])
    imgs = [plot_response(dbg), plot_compression(dbg), plot_harmonics(dbg), plot_music(model), plot_preemphasis(model),
            plot_flutter(model), plot_autogain(model)]
    imgs = [b for b in imgs if b]
    css = """body{font:14px/1.45 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;color:#1f2328;max-width:860px;
margin:32px auto;padding:0 20px}h1{font-size:22px;margin:0 0 4px}p.sub{color:#6e7781;margin:0 0 20px}
table{border-collapse:collapse;width:100%;margin:8px 0 20px}td{padding:5px 8px;border-bottom:1px solid #eaeef2;
vertical-align:top}td:first-child{color:#6e7781;width:210px}.warn{background:#fff8e6;border:1px solid #f0d58c;
border-radius:6px;padding:10px 14px;margin:0 0 18px}img{width:100%;margin:6px 0 14px}pre{background:#f6f8fa;
padding:10px;border-radius:6px;font-size:12px;overflow:auto}"""
    out = [f'<!DOCTYPE html><html><head><meta charset="utf-8"><title>Technical Tape model QC: {esc(str(m.get("tape")))} '
           f'{esc(str(m.get("speed_ips")))} ips</title><style>{css}</style></head><body>',
           f'<h1>{esc(str(m.get("tape")))} at {esc(str(m.get("speed_ips")))} ips</h1>',
           '<p class="sub">Technical Tape model QC report</p>']
    if warn:
        out.append('<div class="warn"><b>Warnings</b><ul>' + ''.join(f'<li>{esc(w)}</li>' for w in warn) + '</ul></div>')
    out.append('<table>' + ''.join(f'<tr><td>{k}</td><td>{v}</td></tr>' for k, v in rows) + '</table>')
    out += [f'<img alt="" src="data:image/png;base64,{b}">' for b in imgs]
    if log_lines:
        out.append('<h3>Build log</h3><pre>' + esc('\n'.join(log_lines)) + '</pre>')
    out.append('</body></html>')
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(out))
    return path

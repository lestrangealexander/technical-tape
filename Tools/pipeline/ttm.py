"""Technical Tape model files (.ttm): little-endian tagged chunks. Read by Source/engine/TapeModel.cpp.

  'TTM1' u32 version
  HEAD  f64 sampleRate, f64 refPeak, i32 numLevels, f32 level0Db, f32 levelStepDb, i32 numBranches,
        i32 preDelay, f32 lookAheadMs, f32 releaseMs, f32 capturedLowVu, f32 capturedHighVu, i32 flags,
        f32 margin
  PEMP  i32 n, f32 freqHz[n], f32 magDb[n]                      (record pre-emphasis estimate |P|)
  BRCH  i32 order, i32 mode, i32 count, i32 length, f16 kernels[count*length], then for mode 0 only
        f32 weights[numLevels*count]
          mode 0: a basis of `count` kernels, mixed per level by the weights table
          mode 1: one kernel per detector level (count == numLevels), interpolated directly
  TRBL  f32 fcLoHz [, f32 fcHiHz [, f32 leak [, f32 rmsMs, f32 peakMix]]]
                                       (version 3+, optional, one per treble band, lowest first: the band's
                                        detector hears through a high-pass at fcLo and, if fcHi > 0, a
                                        low-pass at fcHi, plus `leak` times the rest of the shared
                                        detector's pre-emphasised signal (that signal less the same
                                        band of it); version 3 files have fcLo only, version 4 no leak,
                                        version 5 no rmsMs. rmsMs > 0: the band's detector reads the
                                        energy over a centred window of that length, scaled so a sine
                                        reads its peak, instead of the peak envelope (0), mixed with
                                        peakMix (0..1) of the peak envelope)
  BRHF  as BRCH, order 1: that band's share of the fundamental, indexed by its detector's level (the
        k-th BRHF belongs to the k-th TRBL)
  FLUT  i32 nPeaks, (f32 hz, f32 rmsSec)[nPeaks], i32 nCont, (f32 hz, f32 psd s^2/Hz)[nCont], f32 maxDelaySec
  AGTB  i32 n, f32 drive0Db, f32 stepDb, f32 gainDb[n]
  META  UTF-8 JSON (provenance and QC; the engine ignores it)
"""
from __future__ import annotations

import json
import struct
from dataclasses import dataclass, field

import numpy as np

VERSION = 6


@dataclass
class Branch:
    order: int
    kernels: np.ndarray        # (count, length) float32 (stored as float16)
    weights: np.ndarray | None # (numLevels, count) float32 for a basis; None: one kernel per level


@dataclass
class TrebleBand:
    fc_lo: float               # Hz: the band's detector high-pass, where its share starts
    fc_hi: float               # Hz: its low-pass, where the next band takes over (0 for the top band)
    branch: Branch
    leak: float = 0.0          # how much of the rest of the shared detector's signal its detector also hears
    rms_ms: float = 0.0        # > 0: its detector reads the band's energy over this window (ms), not its peaks
    peak_mix: float = 0.0      # ...mixed with this much (0..1) of its peak envelope


@dataclass
class Model:
    sample_rate: float
    ref_peak: float
    level0_db: float
    level_step_db: float
    num_levels: int
    pre_delay: int
    look_ahead_ms: float
    release_ms: float
    captured_low_vu: float
    captured_high_vu: float
    pemp_freqs: np.ndarray
    pemp_db: np.ndarray
    branches: list
    flutter_peaks: np.ndarray  # (n, 2): hz, rms seconds
    flutter_cont: np.ndarray   # (n, 2): hz, psd s^2/Hz
    flutter_max_delay: float
    ag_drive0: float
    ag_step: float
    ag_gain_db: np.ndarray
    meta: dict = field(default_factory=dict)
    flags: int = 0
    margin: float = 1.02
    treble_bands: list = field(default_factory=list)   # TrebleBand, lowest first; empty: no treble stage

    def levels_db(self):
        return self.level0_db + self.level_step_db * np.arange(self.num_levels)


def _chunk(tag, payload):
    return tag.encode('ascii') + struct.pack('<I', len(payload)) + payload


def _f32(a):
    return np.asarray(a, '<f4').tobytes()


def write(path, m: Model):
    out = [b'TTM1', struct.pack('<I', VERSION)]
    out.append(_chunk('HEAD', struct.pack('<ddiffiiffffif', m.sample_rate, m.ref_peak, m.num_levels, m.level0_db,
                                          m.level_step_db, len(m.branches), m.pre_delay, m.look_ahead_ms,
                                          m.release_ms, m.captured_low_vu, m.captured_high_vu, m.flags, m.margin)))
    out.append(_chunk('PEMP', struct.pack('<i', len(m.pemp_freqs)) + _f32(m.pemp_freqs) + _f32(m.pemp_db)))
    def branch_payload(b):
        r, n = b.kernels.shape
        mode = 1 if b.weights is None else 0
        if mode == 1:
            assert r == m.num_levels
        else:
            assert b.weights.shape == (m.num_levels, r)
        payload = struct.pack('<iiii', b.order, mode, r, n) + np.asarray(b.kernels, '<f2').tobytes()
        if mode == 0:
            payload += _f32(b.weights)
        return payload

    for b in m.branches:
        out.append(_chunk('BRCH', branch_payload(b)))
    for tb in m.treble_bands:
        assert tb.branch.order == 1 and tb.fc_lo > 0
        out.append(_chunk('TRBL', struct.pack('<fffff', tb.fc_lo, tb.fc_hi, tb.leak, tb.rms_ms, tb.peak_mix)))
        out.append(_chunk('BRHF', branch_payload(tb.branch)))
    fp, fc = np.asarray(m.flutter_peaks, float).reshape(-1, 2), np.asarray(m.flutter_cont, float).reshape(-1, 2)
    out.append(_chunk('FLUT', struct.pack('<i', len(fp)) + _f32(fp) + struct.pack('<i', len(fc)) + _f32(fc)
                      + struct.pack('<f', m.flutter_max_delay)))
    out.append(_chunk('AGTB', struct.pack('<iff', len(m.ag_gain_db), m.ag_drive0, m.ag_step) + _f32(m.ag_gain_db)))
    out.append(_chunk('META', json.dumps(m.meta, separators=(',', ':')).encode('utf-8')))
    with open(path, 'wb') as fh:
        fh.write(b''.join(out))


def read(path) -> Model:
    data = open(path, 'rb').read()
    assert data[:4] == b'TTM1', 'not a Technical Tape model'
    assert struct.unpack_from('<I', data, 4)[0] in (2, 3, 4, 5, 6), 'unsupported model version'
    treble_fc, treble_br = [], []
    pos = 8
    branches, kw = [], {}
    head = None
    while pos < len(data):
        tag = data[pos:pos + 4].decode('ascii')
        n = struct.unpack_from('<I', data, pos + 4)[0]
        p = data[pos + 8:pos + 8 + n]
        pos += 8 + n
        if tag == 'HEAD':
            head = struct.unpack('<ddiffiiffffif', p)
        elif tag == 'PEMP':
            k = struct.unpack_from('<i', p)[0]
            a = np.frombuffer(p, '<f4', 2 * k, 4)
            kw['pemp_freqs'], kw['pemp_db'] = a[:k].copy(), a[k:].copy()
        elif tag in ('BRCH', 'BRHF'):
            order, mode, r, ln = struct.unpack_from('<iiii', p)
            kern = np.frombuffer(p, '<f2', r * ln, 16).reshape(r, ln).astype(np.float32)
            w = None
            if mode == 0:
                nl = head[2]
                w = np.frombuffer(p, '<f4', nl * r, 16 + 2 * r * ln).reshape(nl, r).copy()
            if tag == 'BRCH':
                branches.append(Branch(order, kern, w))
            else:
                treble_br.append(Branch(order, kern, w))
        elif tag == 'TRBL':
            lo = struct.unpack_from('<f', p)[0]
            hi = struct.unpack_from('<f', p, 4)[0] if len(p) >= 8 else 0.0
            lk = struct.unpack_from('<f', p, 8)[0] if len(p) >= 12 else 0.0
            rm = struct.unpack_from('<f', p, 12)[0] if len(p) >= 16 else 0.0
            pk = struct.unpack_from('<f', p, 16)[0] if len(p) >= 20 else 0.0
            treble_fc.append((lo, hi, lk, rm, pk))
        elif tag == 'FLUT':
            k = struct.unpack_from('<i', p)[0]
            fp = np.frombuffer(p, '<f4', 2 * k, 4).reshape(k, 2).copy()
            o = 4 + 8 * k
            c = struct.unpack_from('<i', p, o)[0]
            fc = np.frombuffer(p, '<f4', 2 * c, o + 4).reshape(c, 2).copy()
            kw['flutter_peaks'], kw['flutter_cont'] = fp, fc
            kw['flutter_max_delay'] = struct.unpack_from('<f', p, o + 4 + 8 * c)[0]
        elif tag == 'AGTB':
            k, d0, st = struct.unpack_from('<iff', p)
            kw['ag_drive0'], kw['ag_step'] = d0, st
            kw['ag_gain_db'] = np.frombuffer(p, '<f4', k, 12).copy()
        elif tag == 'META':
            kw['meta'] = json.loads(p.decode('utf-8'))
    sr, rp, nl, l0, ls, nb, pd, la, rel, lo, hi, flags, margin = head
    assert len(treble_fc) == len(treble_br), 'treble bands and their branches do not pair up'
    kw['treble_bands'] = [TrebleBand(lo_, hi_, br, lk_, rm_, pk_) for (lo_, hi_, lk_, rm_, pk_), br in zip(treble_fc, treble_br)]
    return Model(sr, rp, l0, ls, nl, pd, la, rel, lo, hi, branches=branches, flags=flags, margin=margin, **kw)

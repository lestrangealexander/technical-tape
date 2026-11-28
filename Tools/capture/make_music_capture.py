#!/usr/bin/env python3
"""Music capture template: tones, a pink-noise check, sync markers and a 60 s gap for a mix of your own.

Recorded through the machine at three levels (-10, 0 and +10 dB into it), it shows how the machine treats
real music as the level changes: what pink noise and sweeps can't pin down (how the highs move with the
level on dense material with a mix's own balance, and how loud it gets). 96 kHz, 24-bit, mono.

  0:00.000  silence 2 s
  0:02.000  1 kHz, 0 VU (-18 dBFS RMS), 8 s          calibration check
  0:11.000  1 kHz, -20 VU, 6 s                       the pass's level, where the tape is linear
  0:18.000  pink noise, 0 VU, 8 s                    the same noise as the main capture's 0 VU section
  0:27.000  sync markers (3 chirps)
  0:30.000  MUSIC: paste a mix here (up to 60 s; mono), loud parts about 0 VU
  1:32.000  sync markers (3 chirps)
  1:36.270  end

usage: make_music_capture.py out.wav [out.json]
"""
import json
import sys

import numpy as np

import make_capture as mc

MUSIC_START_S = 30.0
MUSIC_LEN_S = 60.0


def main(out_wav, out_json=None):
    parts, cues = [], []
    pos = 0.0

    def add(x, label=None, level=None, kind=None):
        nonlocal pos
        if label:
            cues.append({'start_s': round(pos, 3), 'length_s': round(len(x) / mc.FS, 3), 'section': label,
                         'level': level, 'kind': kind})
        parts.append(x)
        pos += len(x) / mc.FS

    add(mc.silence(2.0), 'Silence (noise floor)', None, 'silence')
    add(mc.tone(1000, 8.0, 0.0), '1 kHz calibration tone', '0 VU (-18 dBFS RMS)', 'cal')
    add(mc.silence(1.0))
    add(mc.tone(1000, 6.0, -20.0), '1 kHz level tone', '-20 VU', 'level')
    add(mc.silence(1.0))
    add(mc.pink(8.0, 0, 102), 'Pink noise (as the main capture)', '0 VU (RMS)', 'noise')
    add(mc.silence(1.0))
    add(mc.marker_burst(), 'Sync markers (3 chirps)', '-6 VU', 'sync')
    add(mc.silence(MUSIC_START_S - pos))
    add(mc.silence(MUSIC_LEN_S), 'MUSIC: paste a mono mix here, loud parts about 0 VU', None, 'music')
    add(mc.silence(2.0))
    add(mc.marker_burst(), 'Sync markers (3 chirps)', '-6 VU', 'sync')
    add(mc.silence(3.0), 'Silence (end)', None, 'silence')
    x = mc.upsample_from_48k(np.concatenate(parts), 96000)
    assert float(np.max(np.abs(x))) < 0.99
    mc.write_wav24(out_wav, x, 96000)
    if out_json:
        json.dump({'sample_rate_hz': 96000, 'cues': cues}, open(out_json, 'w'), indent=1)
    for c in cues:
        m, s = divmod(c['start_s'], 60)
        print(f'{int(m)}:{s:06.3f}  {c["section"]}' + (f' ({c["level"]})' if c['level'] else ''))


if __name__ == '__main__':
    main(*sys.argv[1:3])

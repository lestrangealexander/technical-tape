# Recording the captures

Twelve files make the six models: two tapes (ATR, Maxell) × three speeds (3¾, 7½, 15 ips) ×
two passes (normal and hot). The six ATR files are in (November 2026); the Maxell is to come.

## The test file

Play **`Technical Tape Capture 96k 24-bit.wav`** (96 kHz, 24-bit, mono, 5 min 33 s) and record the
return at 96 kHz. It is the same signal as the Unique Vintage Processor test file
(`UVP Tape Capture 48k.wav`), only at 96 kHz: identical to 22.7 kHz, nothing above 24 kHz. The
48 kHz file works as well; record at the rate of the file you play. The pipeline regenerates the
signal internally and checks every capture against it sample by sample, so it must be one of these
two files, unedited. To make fresh copies:

```
python3 Tools/capture/make_capture.py "Technical Tape Capture 96k 24-bit.wav" cues.json 96000
python3 Tools/capture/make_capture.py "UVP Tape Capture 48k.wav" cues.json
```

It contains: silence, a 1 kHz calibration tone at 0 VU (−18 dBFS RMS), two sync bursts,
15 sine sweeps from −30 to +12 VU in 3 dB steps (10 Hz–22 kHz), pink noise at −20, −10, 0 and
+3 VU, a 3150 Hz tone for wow & flutter, two more sync bursts and silence. Its peak is −3 dBFS.

## Setting up the machine

1. Clean and demagnetise the heads; let the machine warm up for 20 minutes.
2. Bias and record EQ as you would normally use that tape at that speed. No noise reduction.
3. Calibrate once per tape and speed: play the test file's 1 kHz tone and set the record level
   so the machine reads **0 VU**, then set the repro (or interface input) so the tone comes back
   at **−18 dBFS RMS** in the DAW. The pipeline measures the calibration from the tone, so a
   few tenths of a dB off is fine; it is reported.
4. Monitor **off the repro head** and record the return while the machine records, in one pass.

## The two passes

- **Normal pass**: the test file as calibrated.
- **Hot pass**: the same file hotter *into the machine*, made in the analogue domain (the
  machine's record level, or the interface's output level if it has the headroom) — not by
  turning up the file, which would clip it (the file peaks at −3 dBFS). Turn the return down so
  nothing clips on the way back in. The exact amounts don't matter: the pipeline measures both
  offsets from where the passes overlap (it searches −3 to +36 dB), so +6 dB is a minimum, not a
  target. Nothing may clip in the converters (watch the interface's input meters on the loudest
  sweep).

The hot pass carries the model above the normal pass's +12 VU, where tape really saturates.
Without it the model still works (the saturation curve is continued from +12 VU), but less of it
is measured.

**What the ATR hot passes were.** All three went onto the tape 22.4–22.7 dB hotter than the
normal pass, with the return about 13 dB lower, the same at every speed: the pipeline measured
it from the odd harmonics and the compression, which line up between the passes to 0.5–1.5 dB
and 0.1–0.2 dB. That is far more than +6 dB, and better for it: the ATR models are measured all
the way to +35 VU, where a 1 kHz tone comes back nearly square. **For the Maxell, repeat exactly
what you did for the ATR hot pass.**

## Files

- 96 kHz (or 48 kHz with the 48 kHz file), 24-bit or 32-bit float, mono or stereo (the louder
  channel is used; `--channel L` or `R` picks one). A 96 kHz capture is converted to 48 kHz with a
  filter that is flat to 22 kHz and lets nothing the tape puts above 24 kHz fold back.
- Start recording a few seconds before playing the test file and stop a few seconds after it
  ends. Don't trim, normalise, fade or process the files in any way.
- Name each file so the tape, speed and pass can be read from it, for example:

```
ATR_3.75ips_0dB.wav     ATR_3.75ips_+6dB.wav
ATR_7.5ips_0dB.wav      ATR_7.5ips_+6dB.wav
ATR_15ips_0dB.wav       ATR_15ips_+6dB.wav
Maxell_3.75ips_0dB.wav  Maxell_3.75ips_+6dB.wav
Maxell_7.5ips_0dB.wav   Maxell_7.5ips_+6dB.wav
Maxell_15ips_0dB.wav    Maxell_15ips_+6dB.wav
```

The names are read loosely: `7½`, `7,5` or `7_5` work for 7.5, `3¾` for 3.75, and `+6`, `6dB`,
`hot` or `plus6` mark the hot pass (an upload that turns `+` into `_`, as in `ATR_15ips__6dB.wav`,
still reads right).

## From captures to the plugin

```
python3 Tools/pipeline/extract.py <folder with the 12 files> -o Resources/Models --dry-run   # check the names were read right
python3 Tools/pipeline/extract.py <folder with the 12 files> -o Resources/Models
```

For each tape and speed this writes `ATR_15.ttm` (etc.) and `ATR_15_report.html`, a QC page
with the calibration, warnings, responses, compression curves, harmonics, the pre-emphasis
estimate and wow & flutter. Then either rebuild the plugin (the models in `Resources/Models` are
built into it), or copy the `.ttm` files into the user models folder to try them at once:

- Mac: `~/Library/Application Support/Sound Aquarium/Technical Tape/Models/`
- Windows: `%APPDATA%\Sound Aquarium\Technical Tape\Models\`

A file there replaces the built-in model of the same name the next time the plugin loads.

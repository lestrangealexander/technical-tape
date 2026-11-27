# Technical Tape — Special 1500

Sound Aquarium Recording. A tape machine plugin made from captures of one real machine: ATR and
Maxell tape at 3¾, 7½ and 15 ips. It works the way Acustica's Nebula and Taupe do — sampled,
level-dependent convolution of the fundamental and its harmonics — with the known weak points of
that approach redone: continuous blending across levels, a look-ahead detector on the
pre-emphasised signal (so high frequencies saturate first, as on tape), alias-free harmonics,
measured wow and flutter, and a saturation curve that keeps going past the hottest capture.
`docs/DESIGN.md` has the research, the design and the validation figures.

**Status:** the three ATR models (3¾, 7½ and 15 ips) are built from the real captures, measured
from −30 VU to about +35 VU: the hot passes went onto the tape 22½ dB above the normal ones, so the
whole drive range of the INPUT knob is measured rather than extrapolated. A treble stage lets the
highs saturate on their own level, as on tape (in two bands at 3¾ ips, where the top octave
saturates apart from the octave below it), so busy material at the slower speeds behaves like the
machine. MAXELL is greyed out until its captures are in. MENU > Tape models shows what each
model is.

The Mac build is AU and VST3, universal for Apple Silicon and Intel, and needs macOS 11 or later.
The Windows build is a 64-bit VST3 for Windows 10 or 11. Builds are for personal use, so they are
not notarised or code-signed.

## Installing on a Mac

1. Download `TechnicalTape-mac.zip` from the latest build on the repository's **Releases** page
   and double-click it to unzip it.
2. Copy the plugin files:
   - `Technical Tape.component` goes in `~/Library/Audio/Plug-Ins/Components/`.
   - `Technical Tape.vst3` goes in `~/Library/Audio/Plug-Ins/VST3/`.

   In Finder, use **Go › Go to Folder…** and paste the path. `~` is your home folder.
3. The plugins aren't notarised, so clear the download quarantine once. Paste this into Terminal:

   ```
   xattr -dr com.apple.quarantine ~/Library/Audio/Plug-Ins/Components/"Technical Tape.component" ~/Library/Audio/Plug-Ins/VST3/"Technical Tape.vst3"
   ```

4. Restart your DAW and rescan plugins. In Logic, AU plugins are validated automatically at launch.

When you update to a newer build, quit the DAW, replace both files, run the same `xattr` line again
and restart the DAW.

## Installing on Windows

1. Download `TechnicalTape-win64.zip`. Before unzipping, right-click it, choose **Properties**, tick
   **Unblock** if it's shown, then click **OK**.
2. Unzip it and copy the whole `Technical Tape.vst3` folder into `C:\Program Files\Common Files\VST3\`.
   Windows will ask for administrator permission; that's normal.
3. Restart your DAW and rescan plugins.

Presets are stored in `~/Library/Audio/Presets/Sound Aquarium/Technical Tape` (Mac) and
`%APPDATA%\Sound Aquarium\Technical Tape\Presets` (Windows).

## The faceplate

| Control | What it does |
|---|---|
| ATR / MAXELL | the tape (the supply reel changes with it); MAXELL is greyed out until its captures are in |
| 3¾ / 7½ / 15 | tape speed |
| INPUT | how hard the tape is driven (−18 … +18 dB); 0 VU = −18 dBFS RMS. At 15 ips the ATR stays under 1 % distortion to about +6 VU, saturates from +9 VU and is close to square-wave clipping above +25 VU; all of it is measured |
| OUTPUT | output level (−18 … +18 dB) |
| AG | auto gain: holds the loudness as INPUT changes (on by default) |
| BYPASS | latency-matched, click-free bypass |
| OS | oversampling 1X … 16X (default 8X: alias-free for every harmonic) |
| MENU | presets, wow & flutter amount (0–200 %, default as measured), interface size, model info |
| INPUT VU / OUTPUT VU | level onto the tape and level out, VU ballistics |

Double-click a knob to type a value; Alt/Option-click returns it to its default. Changing tape,
speed or OS crossfades between two engine cores, so it never clicks. The reels turn while the
host's transport runs.

## Building

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --target TechnicalTape_VST3      # also TechnicalTape_AU on a Mac
```

CMake fetches JUCE 9.0.6 (or use `-DTT_JUCE_PATH=...`). The models in `Resources/Models/*.ttm` are
built into the plugin. On Linux, JUCE needs the usual X11/ALSA/freetype development packages.
GitHub Actions (`.github/workflows/build.yml`) builds a universal Mac AU + VST3 and a Windows VST3
on every push, runs pluginval/auval, renders the editor and publishes the zips.

Developer tools (`-DTT_BUILD_TOOLS=ON`): `TTSnapshot` renders the editor to PNG and clicks through
it; `tt_render` (always built) runs the engine offline on a WAV:

```
tt_render Resources/Models/ATR_15.ttm in.wav out.wav --os 8 --drive 6
```

## Building models from captures

```
python3 Tools/pipeline/extract.py <folder of WAVs> -o Resources/Models --dry-run   # check the names were read right
python3 Tools/pipeline/extract.py <folder of WAVs> -o Resources/Models
```

Each tape and speed gets a `.ttm` model and an HTML QC report. Rebuild the plugin, or drop the
`.ttm` files into the user models folder (`~/Library/Application Support/Sound Aquarium/Technical
Tape/Models` on a Mac, `%APPDATA%\Sound Aquarium\Technical Tape\Models` on Windows) to use them
straight away. A tape lights up in TAPE SELECT as soon as one of its models is there. Python 3.10+
with numpy, scipy and matplotlib. `docs/CAPTURES.md` says how to record the captures.

## Checks

| Script | What it checks |
|---|---|
| `Tests/validate_capture.py <model> <normal.wav> <hot.wav>` | a model against its own captures: the 1 kHz and 3.15 kHz tones of both passes (level, harmonics, waveform match) |
| `Tests/validate_model.py <model> <tape> <speed>` | a model against the virtual machine: sines by level, pink-noise bands, two-tone IMD, transient null |
| `Tests/compare_engine.py <model> <tt_render>` | the C++ engine against the Python reference, sample by sample (−103 dB) |
| `Tests/engine_checks.py <tt_render> <model> [<model2>]` | aliasing by OS factor, clicks when switching, wow & flutter spectrum |

## Layout

```
Source/engine/    the DSP (plain C++, no JUCE): model reader, convolution, oversampling, detector, flutter
Source/           plugin, editor, model bank, presets;  Source/ui/  faceplate drawing and controls
Tools/pipeline/   captures -> models (analyse, build, extract, report, reference engine)
Tools/capture/    the test-signal generator;  Tools/sim/  the virtual tape machine used for testing
Tools/gui/        the approved faceplate mockup generator
Tests/            renderer, comparisons, validation, editor snapshot
Resources/        fonts (DM Sans, Yellowtail), the Pegasus badge, Models/
```

Fonts: DM Sans (SIL Open Font License, `Resources/DMSans-OFL.txt`) and Yellowtail (Apache License 2.0,
`Resources/Yellowtail-LICENSE.txt`).
FFT: PFFFT (BSD-like licence, `Source/engine/pffft/LICENSE.txt`).

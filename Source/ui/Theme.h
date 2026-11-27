// Technical Tape — faceplate theme: palette, type and drawing helpers (Sound Aquarium faceplate style)
#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

namespace tt::ui
{
namespace col
{
    inline const juce::Colour salmon   { 0xFFFEA993 };
    inline const juce::Colour salmonDk { 0xFFEE917B };
    inline const juce::Colour yellow   { 0xFFFDE6B0 };
    inline const juce::Colour yellowLn { 0xFFEFCB86 };
    inline const juce::Colour ink      { 0xFF3A302C };
    inline const juce::Colour white    { 0xFFFFFFFF };
    inline const juce::Colour sky      { 0xFFD1EEFE };
    inline const juce::Colour pack     { 0xFF4A3C35 };
    inline const juce::Colour packLine { 0xFF57483F };
} // namespace col

enum class Weight { Medium, SemiBold, Bold };

/** DM Sans at a CSS-style size (px) with CSS-style letter-spacing (px). */
juce::Font font (Weight w, float sizePx, float letterSpacingPx = 0.0f);
/** Yellowtail (the "Special 1500" script) at a CSS-style size. */
juce::Font scriptFont (float sizePx);

/** DM Sans at about `sizePx`, letter-spaced so `text` is exactly `width` wide (justified lines). */
juce::Font justifiedFont (Weight w, float sizePx, const juce::String& text, float width);
/** Yellowtail sized so `text` is exactly `width` wide. */
juce::Font justifiedScript (const juce::String& text, float width);

/** One line of a justified column: its visible ink runs exactly from x to x + width (the first glyph's
    left edge and the last glyph's right edge, not their advances, so lines in different faces line up).
    DM Sans is letter-spaced to fit; the script (script = true) is sized to fit. */
void drawJustifiedLine (juce::Graphics& g, const juce::String& text, float x, float baselineY, float width,
                        bool script, Weight w, float sizePx, juce::Colour c);

enum class Anchor { Start, Middle, End };

/** Width of a single line as drawText lays it out (trailing letter-spacing left out). */
float textWidth (const juce::Font& f, const juce::String& s);

/** Single line of text with its baseline at y (like SVG text). */
void drawText (juce::Graphics& g, const juce::String& s, float x, float baselineY, const juce::Font& f,
               juce::Colour c, Anchor a = Anchor::Middle);

/** Angles in degrees, 0 = 12 o'clock, clockwise (the faceplate's convention). */
juce::Point<float> polar (float cx, float cy, float r, float deg);
juce::Path arc (float cx, float cy, float r, float a0Deg, float a1Deg);
void strokeRound (juce::Graphics& g, const juce::Path& p, float width);

/** The Pegasus badge, tinted. */
std::unique_ptr<juce::Drawable> pegasus (juce::Colour c);

} // namespace tt::ui

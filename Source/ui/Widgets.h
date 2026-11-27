// Technical Tape — faceplate controls
#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include "Theme.h"

namespace tt::ui
{
//==============================================================================
/** Knobs (white disc, salmon track, white value arc, charcoal pointer), menus, dialogs. */
class TapeLookAndFeel final : public juce::LookAndFeel_V4
{
public:
    TapeLookAndFeel();

    void drawRotarySlider (juce::Graphics&, int x, int y, int w, int h, float pos,
                           float startAngle, float endAngle, juce::Slider&) override;
    juce::Font getPopupMenuFont() override;
    void drawPopupMenuBackground (juce::Graphics&, int width, int height) override;
    juce::Font getAlertWindowTitleFont() override;
    juce::Font getAlertWindowMessageFont() override;
    juce::Font getAlertWindowFont() override;
    juce::Font getTextButtonFont (juce::TextButton&, int buttonHeight) override;
    void drawButtonBackground (juce::Graphics&, juce::Button&, const juce::Colour&, bool over, bool down) override;
    void fillTextEditorBackground (juce::Graphics&, int width, int height, juce::TextEditor&) override;
    void drawTextEditorOutline (juce::Graphics&, int width, int height, juce::TextEditor&) override;
    juce::Rectangle<int> getTooltipBounds (const juce::String& tip, juce::Point<int> screenPos, juce::Rectangle<int> parentArea) override;
    void drawTooltip (juce::Graphics&, const juce::String& text, int width, int height) override;
};

//==============================================================================
/** A slider whose double-click opens type-in instead of resetting (Alt/Option-click still resets). */
class ValueSlider final : public juce::Slider
{
public:
    std::function<void()> onDoubleClick;
    void mouseDoubleClick (const juce::MouseEvent& e) override
    {
        if (onDoubleClick != nullptr) onDoubleClick();
        else juce::Slider::mouseDoubleClick (e);
    }
};

/** Small type-in box that opens over a control. Enter or clicking away applies it, Esc cancels. */
class ValueEntry final : public juce::TextEditor
{
public:
    ValueEntry();
    void open (juce::Rectangle<int> area, const juce::String& text, std::function<void (const juce::String&)> onCommit);
    void dismiss (bool apply);
private:
    std::function<void (const juce::String&)> commit;
};

//==============================================================================
/** A rotary knob with its label underneath; the label shows the value while you touch it. */
class ParamKnob final : public juce::Component
{
public:
    ParamKnob (const juce::String& name, float bodyRadius, float ringWidth, bool bipolar,
               float labelBaselineFromCentre, float labelSize, float labelSpacing, Weight labelWeight = Weight::Bold);
    ~ParamKnob() override;

    ValueSlider slider;
    /** Places the knob so its centre lands on (cx, cy) in the parent. */
    void placeAt (float cx, float cy);
    juce::Point<float> getKnobCentre() const { return { (float) getX() + (float) getWidth() * 0.5f, (float) getY() + ringOuter() }; }
    float getKnobRadius() const { return bodyR; }

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    void mouseEnter (const juce::MouseEvent&) override { repaint(); }
    void mouseExit (const juce::MouseEvent&) override { repaint(); }
    void mouseDrag (const juce::MouseEvent&) override { repaint(); }
    void mouseUp (const juce::MouseEvent&) override { repaint(); }

    juce::String name;
    float bodyR, ringW, labelOffset, labelSize, labelSpacing;
    Weight labelWeight;
    float ringOuter() const { return bodyR + ringW * 0.5f + std::max (4.0f, bodyR * 0.06f) + ringW * 0.5f; }
};

//==============================================================================
/** TAPE SELECT / SPEED SELECT: a white-bordered box like the badge's, a header over one row per choice.
    The chosen row is solid white with salmon letters. Drawn 0.5 px in, so its edges land where the
    mockup's do on half-pixel positions (x = 183.5, 666.5). */
class SelectBox final : public juce::Component, public juce::TooltipClient
{
public:
    SelectBox (const juce::String& header, juce::StringArray labels, float textSize, float spacing);
    std::function<void (int)> onPick;
    void setSelected (int index);
    int getSelected() const noexcept { return selected; }
    int getNumChoices() const noexcept { return labels.size(); }
    juce::String getHeader() const { return header; }
    /** Picks a row as a click would; a greyed-out row does nothing. */
    void pick (int index) { if (isRowEnabled (index) && onPick != nullptr) onPick (index); }

    /** A row that can't be chosen yet is drawn faint and ignores the mouse; its tooltip says why. */
    void setRowEnabled (int index, bool enabled, const juce::String& tooltip = {});
    bool isRowEnabled (int index) const noexcept;
    juce::String getTooltip() override;

    void paint (juce::Graphics&) override;
    void mouseMove (const juce::MouseEvent&) override;
    void mouseExit (const juce::MouseEvent&) override;
    void mouseDown (const juce::MouseEvent&) override;
    void mouseUp (const juce::MouseEvent&) override;

    static constexpr float kSize = 150.0f, kHead = 38.0f, kBorder = 3.0f, kRadius = 12.5f;
    static constexpr float kDisabledAlpha = 0.4f;   // a row that can't be chosen yet

private:
    int rowAt (juce::Point<float> p) const;
    juce::Rectangle<float> rowArea (int index) const;
    juce::String header;
    juce::StringArray labels;
    float textSize, spacing;
    int selected = 0, hover = -1, pressed = -1;
    std::vector<bool> disabled;
    juce::StringArray tips;
};

/** Rounded pill (BYPASS, OS): white, charcoal when on. */
class PillButton final : public juce::Button
{
public:
    PillButton (const juce::String& label, float textSize = 12.0f, float spacing = 1.4f, Weight w = Weight::SemiBold);
    void paintButton (juce::Graphics&, bool over, bool down) override;
    void setLabel (const juce::String& s) { setButtonText (s); repaint(); }
private:
    float size, spacing;
    Weight weight;
};

/** AG and the speed buttons: a 56 px roller with a white ring. On: solid white, salmon letters.
    Off: salmon, white letters. */
class RoundButton final : public juce::Button
{
public:
    RoundButton (const juce::String& label, bool toggles);
    void paintButton (juce::Graphics&, bool over, bool down) override;
    bool hitTest (int x, int y) override;
};

/** MENU / SAVE halves of the white pill in the top strip. */
class ZoneButton final : public juce::Button
{
public:
    ZoneButton (const juce::String& label, bool leftHalf, bool burger);
    void paintButton (juce::Graphics&, bool over, bool down) override;
private:
    bool left, burger;
};

//==============================================================================
/** Round VU meter (UVP's meter face, a VU scale): -20 .. +3 VU, 0 VU = -18 dBFS RMS. */
class TapeVu final : public juce::Component
{
public:
    explicit TapeVu (const juce::String& label) : name (label) { setInterceptsMouseClicks (false, false); }
    /** Needle position in VU (already ballistic). */
    void setVu (float vu);
    static float angleFor (float vu) noexcept;
    void paint (juce::Graphics&) override;
private:
    juce::String name;
    float value = -60.0f;
};

//==============================================================================
/** One tape reel (NAB hub, three windows or two slots), turned to an angle. */
class Reel final : public juce::Component
{
public:
    enum class Style { ThreeWindows, TwoSlots };
    struct Colours { juce::Colour flange, back, rim, screw; };
    static Colours skyBlue();
    static Colours white();

    Reel (Style s, Colours c, float packRadius, float turnDeg);
    void setLook (Style s, Colours c);
    void setTurn (float deg);
    float getTurn() const noexcept { return turn; }
    void paint (juce::Graphics&) override;

    static constexpr float kRadius = 110.0f;

private:
    void build();
    Style style;
    Colours colours;
    float packR, turn;
    juce::Path flange, keys, back;
    std::vector<juce::Point<float>> screws;
};

} // namespace tt::ui

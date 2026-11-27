// Technical Tape — faceplate controls
#include "Widgets.h"

#include <map>
#include <set>

namespace tt::ui
{
//==============================================================================
TapeLookAndFeel::TapeLookAndFeel()
{
    setColour (juce::PopupMenu::backgroundColourId, col::white);
    setColour (juce::PopupMenu::textColourId, col::ink);
    setColour (juce::PopupMenu::headerTextColourId, col::ink.withAlpha (0.6f));
    setColour (juce::PopupMenu::highlightedBackgroundColourId, col::yellow);
    setColour (juce::PopupMenu::highlightedTextColourId, col::ink);
    setColour (juce::AlertWindow::backgroundColourId, col::white);
    setColour (juce::AlertWindow::textColourId, col::ink);
    setColour (juce::AlertWindow::outlineColourId, col::salmon);
    setColour (juce::TextButton::buttonColourId, col::yellow);
    setColour (juce::TextButton::buttonOnColourId, col::ink);
    setColour (juce::TextButton::textColourOffId, col::ink);
    setColour (juce::TextButton::textColourOnId, col::white);
    setColour (juce::TextEditor::backgroundColourId, col::white);
    setColour (juce::TextEditor::textColourId, col::ink);
    setColour (juce::TextEditor::outlineColourId, col::salmonDk);
    setColour (juce::TextEditor::focusedOutlineColourId, col::ink);
    setColour (juce::TextEditor::highlightColourId, col::yellow);
    setColour (juce::CaretComponent::caretColourId, col::ink);
}

void TapeLookAndFeel::drawRotarySlider (juce::Graphics& g, int x, int y, int w, int h, float pos, float, float, juce::Slider& s)
{
    const auto& pr = s.getProperties();
    const float bodyR = (float) pr.getWithDefault ("bodyR", 20.0);
    const float ringW = (float) pr.getWithDefault ("ringW", 4.0);
    const bool bipolar = (bool) pr.getWithDefault ("bipolar", false);
    const float cx = (float) x + (float) w * 0.5f, cy = (float) y + (float) h * 0.5f;
    const float a0 = -135.0f, a1 = 135.0f, av = a0 + (a1 - a0) * pos;
    const float rr = bodyR + ringW * 0.5f + std::max (4.0f, bodyR * 0.06f);

    g.setColour (col::salmonDk);
    strokeRound (g, arc (cx, cy, rr, a0, a1), ringW);
    g.setColour (col::white);
    if (bipolar)
    {
        const float s0 = std::min (0.0f, av), s1 = std::max (0.0f, av);
        if (s1 - s0 > 0.5f)
            strokeRound (g, arc (cx, cy, rr, s0, s1), ringW);
        else
        {
            const auto p = polar (cx, cy, rr, 0.0f);
            g.fillEllipse (p.x - ringW * 0.5f, p.y - ringW * 0.5f, ringW, ringW);
        }
    }
    else if (av - a0 > 0.5f)
        strokeRound (g, arc (cx, cy, rr, a0, av), ringW);

    g.fillEllipse (cx - bodyR, cy - bodyR, 2.0f * bodyR, 2.0f * bodyR);
    const auto p1 = polar (cx, cy, bodyR * 0.42f, av), p2 = polar (cx, cy, bodyR * 0.8f, av);
    g.setColour (col::ink);
    juce::Path ptr;
    ptr.startNewSubPath (p1);
    ptr.lineTo (p2);
    strokeRound (g, ptr, std::max (2.6f, bodyR * 0.06f));
}

juce::Font TapeLookAndFeel::getPopupMenuFont() { return font (Weight::Medium, 14.0f, 0.4f); }

void TapeLookAndFeel::drawPopupMenuBackground (juce::Graphics& g, int width, int height)
{
    g.fillAll (col::white);
    g.setColour (col::salmon.withAlpha (0.6f));
    g.drawRect (0, 0, width, height, 1);
}

// tooltips: a small white pill with charcoal text, like the top strip's controls
juce::Rectangle<int> TapeLookAndFeel::getTooltipBounds (const juce::String& tip, juce::Point<int> screenPos, juce::Rectangle<int> parentArea)
{
    const int w = (int) std::ceil (textWidth (font (Weight::SemiBold, 13.0f, 0.4f), tip)) + 24, h = 26;
    return juce::Rectangle<int> (screenPos.x > parentArea.getCentreX() ? screenPos.x - (w + 12) : screenPos.x + 16,
                                 screenPos.y > parentArea.getCentreY() ? screenPos.y - (h + 8) : screenPos.y + 18, w, h)
               .constrainedWithin (parentArea);
}

void TapeLookAndFeel::drawTooltip (juce::Graphics& g, const juce::String& text, int width, int height)
{
    const auto r = juce::Rectangle<float> ((float) width, (float) height);
    g.setColour (col::white);
    g.fillRoundedRectangle (r, (float) height * 0.5f);
    g.setColour (col::ink.withAlpha (0.18f));
    g.drawRoundedRectangle (r.reduced (0.5f), (float) height * 0.5f - 0.5f, 1.0f);
    drawText (g, text, (float) width * 0.5f, (float) height * 0.5f + 13.0f * 0.36f, font (Weight::SemiBold, 13.0f, 0.4f), col::ink);
}

juce::Font TapeLookAndFeel::getAlertWindowTitleFont() { return font (Weight::Bold, 16.0f, 1.2f); }
juce::Font TapeLookAndFeel::getAlertWindowMessageFont() { return font (Weight::Medium, 14.0f); }
juce::Font TapeLookAndFeel::getAlertWindowFont() { return font (Weight::Medium, 14.0f); }
juce::Font TapeLookAndFeel::getTextButtonFont (juce::TextButton&, int) { return font (Weight::SemiBold, 13.0f, 1.2f); }

void TapeLookAndFeel::drawButtonBackground (juce::Graphics& g, juce::Button& b, const juce::Colour& bg, bool over, bool down)
{
    auto r = b.getLocalBounds().toFloat().reduced (1.0f);
    auto c = bg;
    if (down) c = c.darker (0.08f);
    else if (over) c = c.brighter (0.04f);
    g.setColour (c);
    g.fillRoundedRectangle (r, r.getHeight() * 0.5f);
}

void TapeLookAndFeel::fillTextEditorBackground (juce::Graphics& g, int w, int h, juce::TextEditor& te)
{
    g.setColour (te.findColour (juce::TextEditor::backgroundColourId).withAlpha (1.0f));
    g.fillRoundedRectangle (0.0f, 0.0f, (float) w, (float) h, std::min (10.0f, (float) h * 0.5f));
}

void TapeLookAndFeel::drawTextEditorOutline (juce::Graphics& g, int w, int h, juce::TextEditor& te)
{
    const auto c = te.findColour (te.hasKeyboardFocus (true) ? juce::TextEditor::focusedOutlineColourId
                                                              : juce::TextEditor::outlineColourId);
    if (c.isTransparent()) return;
    g.setColour (c);
    g.drawRoundedRectangle (0.75f, 0.75f, (float) w - 1.5f, (float) h - 1.5f, std::min (10.0f, (float) h * 0.5f) - 0.75f, 1.5f);
}

//==============================================================================
ValueEntry::ValueEntry()
{
    setJustification (juce::Justification::centred);
    setFont (font (Weight::SemiBold, 14.0f, 0.3f));
    setIndents (6, 0);
    setBorder ({});
    setSelectAllWhenFocused (true);
    setScrollbarsShown (false);
    setPopupMenuEnabled (false);
    setInputRestrictions (16);
    setColour (juce::TextEditor::backgroundColourId, col::white.withAlpha (0.995f));
    setColour (juce::TextEditor::outlineColourId, col::salmonDk);
    setColour (juce::TextEditor::focusedOutlineColourId, col::ink.withAlpha (0.45f));
    setColour (juce::TextEditor::textColourId, col::ink);
    setColour (juce::TextEditor::highlightColourId, col::yellow);
    setColour (juce::TextEditor::highlightedTextColourId, col::ink);
    onReturnKey = [this] { dismiss (true); };
    onEscapeKey = [this] { dismiss (false); };
    onFocusLost = [this] { dismiss (true); };
    setVisible (false);
}

void ValueEntry::open (juce::Rectangle<int> area, const juce::String& text, std::function<void (const juce::String&)> onCommit)
{
    dismiss (true);
    commit = std::move (onCommit);
    setBounds (area);
    setText (text, false);
    setVisible (true);
    toFront (false);
    grabKeyboardFocus();
    selectAll();
}

void ValueEntry::dismiss (bool apply)
{
    if (! isVisible()) return;
    auto cb = std::move (commit);
    commit = nullptr;
    const auto t = getText();
    setVisible (false);
    if (apply && cb != nullptr) cb (t);
}

//==============================================================================
ParamKnob::ParamKnob (const juce::String& n, float bodyRadius, float ringWidth, bool bipolar, float labelBaselineFromCentre,
                      float lSize, float lSpacing, Weight lWeight)
    : name (n), bodyR (bodyRadius), ringW (ringWidth), labelOffset (labelBaselineFromCentre),
      labelSize (lSize), labelSpacing (lSpacing), labelWeight (lWeight)
{
    slider.setSliderStyle (juce::Slider::RotaryHorizontalVerticalDrag);
    slider.setTextBoxStyle (juce::Slider::NoTextBox, false, 0, 0);
    slider.setRotaryParameters (juce::MathConstants<float>::pi * 1.25f, juce::MathConstants<float>::pi * 2.75f, true);
    slider.setMouseDragSensitivity (240);
    slider.getProperties().set ("bodyR", bodyRadius);
    slider.getProperties().set ("ringW", ringWidth);
    slider.getProperties().set ("bipolar", bipolar);
    slider.setScrollWheelEnabled (true);
    addAndMakeVisible (slider);
    slider.addMouseListener (this, false);
}

ParamKnob::~ParamKnob() { slider.removeMouseListener (this); }

void ParamKnob::placeAt (float cx, float cy)
{
    const float ro = ringOuter();
    const float halfW = std::max (ro, 52.0f);
    const float hgt = ro + labelOffset + 8.0f;
    setBounds (juce::roundToInt (cx - halfW), juce::roundToInt (cy - ro), juce::roundToInt (2.0f * halfW), juce::roundToInt (hgt));
}

void ParamKnob::resized()
{
    const float ro = ringOuter();
    const float cx = (float) getWidth() * 0.5f;
    slider.setBounds (juce::Rectangle<float> (cx - ro, 0.0f, 2.0f * ro, 2.0f * ro).getSmallestIntegerContainer());
}

void ParamKnob::paint (juce::Graphics& g)
{
    const bool touching = slider.isMouseOverOrDragging() || slider.isMouseButtonDown();
    const auto s = touching ? slider.getTextFromValue (slider.getValue()).toUpperCase() : name;
    drawText (g, s, (float) getWidth() * 0.5f, ringOuter() + labelOffset, font (labelWeight, labelSize, labelSpacing), col::ink);
}

//==============================================================================
SelectBox::SelectBox (const juce::String& h, juce::StringArray l, float size, float sp)
    : header (h), labels (std::move (l)), textSize (size), spacing (sp) {}

void SelectBox::setSelected (int index)
{
    if (index != selected)
    {
        selected = index;
        repaint();
    }
}

void SelectBox::setRowEnabled (int index, bool enabled, const juce::String& tooltip)
{
    if (index < 0 || index >= labels.size()) return;
    disabled.resize ((size_t) labels.size(), false);
    tips.ensureStorageAllocated (labels.size());
    while (tips.size() < labels.size()) tips.add ({});
    if (disabled[(size_t) index] == ! enabled && tips[index] == tooltip) return;
    disabled[(size_t) index] = ! enabled;
    tips.set (index, tooltip);
    if (! enabled && hover == index) hover = -1;
    repaint();
}

bool SelectBox::isRowEnabled (int index) const noexcept
{
    if (index < 0 || index >= labels.size()) return false;
    return (size_t) index >= disabled.size() || ! disabled[(size_t) index];
}

juce::String SelectBox::getTooltip()
{
    const int r = rowAt (getMouseXYRelative().toFloat());
    return r >= 0 && r < tips.size() ? tips[r] : juce::String();
}

juce::Rectangle<float> SelectBox::rowArea (int index) const
{
    // inside the border: the header, then the rows with a 3 px rule between each
    const int n = std::max (1, labels.size());
    const float x0 = 0.5f + kBorder, w = kSize - 2.0f * kBorder;
    const float top = kHead, bottom = kSize - kBorder;
    const float rowH = (bottom - top - kBorder * (float) (n - 1)) / (float) n;
    return { x0, top + (float) index * (rowH + kBorder), w, rowH };
}

int SelectBox::rowAt (juce::Point<float> p) const
{
    for (int i = 0; i < labels.size(); ++i)
        if (rowArea (i).expanded (0.0f, kBorder * 0.5f).contains (p)) return i;
    return -1;
}

void SelectBox::mouseMove (const juce::MouseEvent& e)
{
    int h = rowAt (e.position);
    if (! isRowEnabled (h)) h = -1;
    if (h != hover) { hover = h; repaint(); }
    setMouseCursor (h >= 0 ? juce::MouseCursor::PointingHandCursor : juce::MouseCursor::NormalCursor);
}

void SelectBox::mouseExit (const juce::MouseEvent&)
{
    if (hover != -1) { hover = -1; repaint(); }
}

void SelectBox::mouseDown (const juce::MouseEvent& e)
{
    pressed = rowAt (e.position);
    if (! isRowEnabled (pressed)) pressed = -1;
    repaint();
}

void SelectBox::mouseUp (const juce::MouseEvent& e)
{
    const int r = rowAt (e.position);
    const bool hit = r >= 0 && r == pressed;
    pressed = -1;
    repaint();
    if (hit) pick (r);
}

void SelectBox::paint (juce::Graphics& g)
{
    const juce::Rectangle<float> outer (0.5f, 0.0f, kSize, kSize);
    const auto inner = outer.reduced (kBorder);
    juce::Path clip;
    clip.addRoundedRectangle (inner, kRadius - kBorder);

    g.setColour (col::salmon);
    g.fillRoundedRectangle (outer, kRadius);
    {
        juce::Graphics::ScopedSaveState save (g);
        g.reduceClipRegion (clip);
        for (int i = 0; i < labels.size(); ++i)
        {
            if (i == selected) g.setColour (col::white);
            else if (i == pressed) g.setColour (col::salmonDk);
            else if (i == hover) g.setColour (col::salmon.interpolatedWith (col::salmonDk, 0.5f));
            else continue;
            g.fillRect (rowArea (i));
        }
    }
    // rules under the header and between the rows, then the border over everything
    g.setColour (col::white);
    g.fillRect (juce::Rectangle<float> (outer.getX(), kHead - kBorder, kSize, kBorder));
    for (int i = 0; i + 1 < labels.size(); ++i)
        g.fillRect (juce::Rectangle<float> (outer.getX(), rowArea (i).getBottom(), kSize, kBorder));
    g.drawRoundedRectangle (outer.reduced (kBorder * 0.5f), kRadius - kBorder * 0.5f, kBorder);

    const float cx = outer.getCentreX();
    drawText (g, header, cx + 1.1f, kHead * 0.5f + 12.5f * 0.36f, font (Weight::Bold, 12.5f, 2.2f), col::white);
    for (int i = 0; i < labels.size(); ++i)
        drawText (g, labels[i], cx + spacing * 0.5f, rowArea (i).getCentreY() + textSize * 0.36f, font (Weight::Bold, textSize, spacing),
                  i == selected ? col::salmon : (isRowEnabled (i) ? col::white : col::white.withAlpha (kDisabledAlpha)));
}

PillButton::PillButton (const juce::String& label, float textSize, float sp, Weight w)
    : juce::Button (label), size (textSize), spacing (sp), weight (w) {}

void PillButton::paintButton (juce::Graphics& g, bool over, bool down)
{
    const bool on = getToggleState();
    auto b = getLocalBounds().toFloat();
    auto body = on ? col::ink : col::white;
    if (! on && (over || down)) body = col::white.interpolatedWith (col::yellow, down ? 0.7f : 0.45f);
    if (on && over) body = col::ink.brighter (0.12f);
    g.setColour (body);
    g.fillRoundedRectangle (b, b.getHeight() * 0.5f);
    drawText (g, getButtonText(), b.getCentreX(), b.getCentreY() + size * 0.36f, font (weight, size, spacing), on ? col::white : col::ink);
}

RoundButton::RoundButton (const juce::String& label, bool toggles) : juce::Button (label)
{
    setClickingTogglesState (toggles);
}

bool RoundButton::hitTest (int x, int y)
{
    const auto c = getLocalBounds().toFloat().getCentre();
    return c.getDistanceFrom ({ (float) x, (float) y }) <= (float) getWidth() * 0.5f;
}

void RoundButton::paintButton (juce::Graphics& g, bool over, bool down)
{
    const bool on = getToggleState();
    const auto b = getLocalBounds().toFloat();
    auto body = on ? col::white : col::salmon;
    if (over || down) body = on ? col::white.interpolatedWith (col::yellow, down ? 0.6f : 0.35f)
                                : col::salmon.interpolatedWith (col::salmonDk, down ? 0.8f : 0.5f);
    g.setColour (body);
    g.fillEllipse (b);
    g.setColour (col::white);
    g.drawEllipse (b.reduced (1.5f), 3.0f);
    drawText (g, getButtonText(), b.getCentreX() + 0.6f, b.getCentreY() + 18.0f * 0.36f, font (Weight::Bold, 18.0f, 1.2f),
              on ? col::salmon : col::white);
}

ZoneButton::ZoneButton (const juce::String& label, bool leftHalf, bool withBurger)
    : juce::Button (label), left (leftHalf), burger (withBurger) {}

void ZoneButton::paintButton (juce::Graphics& g, bool over, bool down)
{
    auto b = getLocalBounds().toFloat();
    if (over || down)
    {
        juce::Path p;
        p.addRoundedRectangle (b.getX(), b.getY(), b.getWidth(), b.getHeight(), 12.0f, 12.0f, left, ! left, left, ! left);
        g.setColour (col::yellow.withAlpha (down ? 0.9f : 0.6f));
        g.fillPath (p);
    }
    const auto f = font (Weight::SemiBold, 12.0f, 1.6f);
    if (burger)
    {
        g.setColour (col::ink);
        for (int i = 0; i < 3; ++i)
        {
            juce::Path l;
            l.startNewSubPath (15.0f, 7.0f + 5.0f * (float) i);
            l.lineTo (29.0f, 7.0f + 5.0f * (float) i);
            strokeRound (g, l, 2.0f);
        }
        // long preset names first lose some letter-spacing, then get cut short with an ellipsis
        auto label = getButtonText();
        const float room = b.getWidth() - 40.0f - 8.0f;
        auto nf = f;
        if (textWidth (nf, label) > room)
            nf = font (Weight::SemiBold, 12.0f, 0.5f);
        if (textWidth (nf, label) > room)
        {
            const auto dots = juce::String::fromUTF8 ("\xe2\x80\xa6");
            while (label.isNotEmpty() && textWidth (nf, label.trimEnd() + dots) > room)
                label = label.dropLastCharacters (1);
            label = label.trimEnd() + dots;
        }
        drawText (g, label, 40.0f, 16.6f, nf, col::ink, Anchor::Start);
    }
    else
        drawText (g, getButtonText(), b.getCentreX(), 16.6f, f, col::ink);
}

//==============================================================================
float TapeVu::angleFor (float vu) noexcept
{
    return std::min (58.0f, -50.0f + 100.0f * std::pow (10.0f, (vu - 3.0f) / 20.0f));
}

void TapeVu::setVu (float vu)
{
    if (std::abs (angleFor (vu) - angleFor (value)) > 0.05f)
    {
        value = vu;
        repaint();
    }
}

void TapeVu::paint (juce::Graphics& g)
{
    const float R = (float) getWidth() * 0.5f;
    const float cx = R, cy = R;
    g.setColour (col::white);
    g.fillEllipse (0.0f, 0.0f, 2.0f * R, 2.0f * R);
    const float px = cx, py = cy + R * 0.6f;
    const float band = R * 1.04f - 6.0f;

    g.setColour (col::yellow);
    strokeRound (g, arc (px, py, band, angleFor (-20.0f), angleFor (0.0f)), 13.0f);
    g.setColour (col::salmon);
    strokeRound (g, arc (px, py, band, angleFor (0.0f), angleFor (3.0f)), 13.0f);
    for (int v : { -20, -10, -7, -5, -3, -2, -1, 0, 1, 2, 3 })
    {
        const float a = angleFor ((float) v);
        const bool big = v == -20 || v == -10 || v == -5 || v == 0 || v == 3;
        juce::Path t;
        t.startNewSubPath (polar (px, py, band - 13.0f, a));
        t.lineTo (polar (px, py, band - (big ? 22.0f : 18.0f), a));
        g.setColour (col::ink.withAlpha (big ? 0.4f : 0.22f));
        strokeRound (g, t, big ? 1.8f : 1.2f);
        if (big)
        {
            const auto pt = polar (px, py, band - 36.0f, a);
            const juce::String s = v == 3 ? "+3" : juce::String (std::abs (v));
            drawText (g, s, pt.x, pt.y + 4.5f, font (Weight::SemiBold, 13.0f), col::ink.withAlpha (0.6f));
        }
    }
    drawText (g, name, cx, py - 24.0f, font (Weight::Bold, 12.0f, 2.4f), col::ink.withAlpha (0.55f));
    const float an = angleFor (value);
    juce::Path needle;
    needle.startNewSubPath (polar (px, py, 13.0f, an));
    needle.lineTo (polar (px, py, band + 4.0f, an));
    g.setColour (col::ink.withAlpha (0.85f));
    strokeRound (g, needle, 2.4f);
    g.setColour (col::salmon);
    g.fillEllipse (px - 9.0f, py - 9.0f, 18.0f, 18.0f);
    g.setColour (col::white);
    g.fillEllipse (px - 3.4f, py - 3.4f, 6.8f, 6.8f);
}

//==============================================================================
namespace
{
    juce::Point<float> rot (juce::Point<float> p, float deg)
    {
        const float a = juce::degreesToRadians (deg);
        return { p.x * std::cos (a) - p.y * std::sin (a), p.x * std::sin (a) + p.y * std::cos (a) };
    }

    /** Closed polygon with some corners rounded (quadratic, `radius` along each side): as the mockup. */
    juce::Path roundedPoly (const std::vector<juce::Point<float>>& P, const std::map<int, float>& radii)
    {
        const int n = (int) P.size();
        auto walk = [&] (int i, float dist, int d, std::set<int>& passed)
        {
            int j = i;
            float rem = dist;
            for (int guard = 0; guard < n; ++guard)
            {
                const int k = ((j + d) % n + n) % n;
                const float seg = P[(size_t) k].getDistanceFrom (P[(size_t) j]);
                if (seg >= rem - 1e-6f)
                    return P[(size_t) j] + (P[(size_t) k] - P[(size_t) j]) * (rem / seg);
                rem -= seg;
                passed.insert (k);
                j = k;
            }
            return P[(size_t) j];
        };
        std::set<int> skip;
        std::map<int, std::pair<juce::Point<float>, juce::Point<float>>> special;
        for (const auto& [c, r] : radii)
        {
            const auto a = walk (c, r, -1, skip);
            const auto b = walk (c, r, 1, skip);
            special[c] = { a, b };
        }
        juce::Path p;
        bool first = true;
        for (int i = 0; i < n; ++i)
        {
            if (skip.count (i) != 0) continue;
            if (auto it = special.find (i); it != special.end())
            {
                if (first) p.startNewSubPath (it->second.first); else p.lineTo (it->second.first);
                p.quadraticTo (P[(size_t) i], it->second.second);
            }
            else
            {
                if (first) p.startNewSubPath (P[(size_t) i]); else p.lineTo (P[(size_t) i]);
            }
            first = false;
        }
        p.closeSubPath();
        return p;
    }
} // namespace

Reel::Colours Reel::skyBlue() { return { col::sky, juce::Colour (0xFFA6CBE8), juce::Colour (0xFFB5D8F1), juce::Colour (0xFF9DC3E2) }; }
Reel::Colours Reel::white() { return { col::white, juce::Colour (0xFFE6D8D1), juce::Colour (0xFFEDE3DD), juce::Colour (0xFFDCCBC3) }; }

Reel::Reel (Style s, Colours c, float packRadius, float turnDeg) : style (s), colours (c), packR (packRadius), turn (turnDeg)
{
    setInterceptsMouseClicks (false, false);
    build();
}

void Reel::setLook (Style s, Colours c)
{
    style = s;
    colours = c;
    build();
    repaint();
}

void Reel::setTurn (float deg)
{
    deg = std::fmod (deg, 360.0f);
    if (std::abs (deg - turn) > 0.01f)
    {
        turn = deg;
        repaint();
    }
}

void Reel::build()
{
    // everything at turn 0, centred on the origin; paint() turns and places it
    constexpr float R = kRadius;
    const float hole = 0.31f * R;
    flange.clear();
    flange.setUsingNonZeroWinding (false);
    flange.addEllipse (-R, -R, 2.0f * R, 2.0f * R);
    if (style == Style::ThreeWindows)
    {
        const float ri = 0.38f * R, ro = 0.86f * R, hi = 0.13f * R;
        const float tho = juce::radiansToDegrees (std::asin (0.33f * R / ro));
        std::vector<juce::Point<float>> base { { -hi, -ri }, { hi, -ri } };
        for (int k = 0; k < 24; ++k)
            base.push_back (polar (0.0f, 0.0f, ro, tho - 2.0f * tho * (float) k / 23.0f));
        const int n = (int) base.size();
        for (float a : { 0.0f, 120.0f, 240.0f })
        {
            std::vector<juce::Point<float>> pts;
            for (const auto& q : base) pts.push_back (rot (q, a));
            flange.addPath (roundedPoly (pts, { { 0, 5.0f }, { 1, 5.0f }, { 2, 7.0f }, { n - 1, 7.0f } }));
        }
    }
    else
    {
        const float x0 = 0.38f * R, x1 = 0.84f * R, h = 0.12f * R, r = 6.0f;
        for (float side : { 0.0f, 180.0f })
        {
            juce::Path slot;
            slot.addRoundedRectangle (x0, -h, x1 - x0, 2.0f * h, r);
            slot.applyTransform (juce::AffineTransform::rotation (juce::degreesToRadians (side)));
            flange.addPath (slot);
        }
    }
    flange.addEllipse (-hole, -hole, 2.0f * hole, 2.0f * hole);

    keys.clear();
    for (float a : { 0.0f, 120.0f, 240.0f })
    {
        const juce::Point<float> q[] = { rot ({ -4.0f, -hole }, a), rot ({ 4.0f, -hole }, a), rot ({ 4.0f, -0.255f * R }, a),
                                         rot ({ -4.0f, -0.255f * R }, a) };
        keys.startNewSubPath (q[0]);
        for (int i = 1; i < 4; ++i) keys.lineTo (q[i]);
        keys.closeSubPath();
    }
    screws.clear();
    if (style == Style::ThreeWindows)
        for (float a : { 60.0f, 180.0f, 300.0f }) screws.push_back (rot ({ 0.0f, -0.36f * R }, a));
    else
        for (float a : { 45.0f, 135.0f, 225.0f, 315.0f }) screws.push_back (rot ({ 0.0f, -0.36f * R }, a));
}

void Reel::paint (juce::Graphics& g)
{
    constexpr float R = kRadius;
    const auto place = juce::AffineTransform::rotation (juce::degreesToRadians (turn))
                           .translated ((float) getWidth() * 0.5f, (float) getHeight() * 0.5f);
    juce::Graphics::ScopedSaveState save (g);
    g.addTransform (place);

    g.setColour (colours.back);
    g.fillEllipse (-(R - 2.0f), -(R - 2.0f), 2.0f * (R - 2.0f), 2.0f * (R - 2.0f));
    g.setColour (col::pack);
    g.fillEllipse (-packR, -packR, 2.0f * packR, 2.0f * packR);
    g.setColour (col::packLine);
    for (float k : { 8.0f, 18.0f, 30.0f })
        if (packR - k > 0.4f * R)
            g.drawEllipse (-(packR - k), -(packR - k), 2.0f * (packR - k), 2.0f * (packR - k), 1.0f);
    g.setColour (colours.flange);
    g.fillPath (flange);
    g.setColour (colours.rim);
    g.drawEllipse (-(R - 1.25f), -(R - 1.25f), 2.0f * (R - 1.25f), 2.0f * (R - 1.25f), 2.5f);

    // NAB centre: the reel's drive keys, then the deck's white adapter
    const float hole = 0.31f * R;
    g.setColour (col::ink);
    g.fillEllipse (-hole, -hole, 2.0f * hole, 2.0f * hole);
    g.setColour (colours.flange);
    g.fillPath (keys);
    const float ad = 0.245f * R;
    g.setColour (col::white);
    g.fillEllipse (-ad, -ad, 2.0f * ad, 2.0f * ad);
    g.setColour (col::ink.withAlpha (0.75f));
    for (float a : { 60.0f, 180.0f, 300.0f })
    {
        const auto p = polar (0.0f, 0.0f, 0.155f * R, a);
        g.fillEllipse (p.x - 2.6f, p.y - 2.6f, 5.2f, 5.2f);
    }
    g.setColour (col::ink);
    g.fillEllipse (-7.0f, -7.0f, 14.0f, 14.0f);
    g.setColour (col::white);
    g.fillEllipse (-2.6f, -2.6f, 5.2f, 5.2f);
    g.setColour (colours.screw);
    for (const auto& s : screws)
        g.fillEllipse (s.x - 2.4f, s.y - 2.4f, 4.8f, 4.8f);
}

} // namespace tt::ui

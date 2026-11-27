// Technical Tape — editor (the Special 1500 faceplate)
//
// Everything is laid out on a fixed 1000 x 660 faceplate in the same coordinates as the approved
// mockup (Tools/gui/gen_tt.py), then scaled as a whole for the interface size option.
#include "PluginEditor.h"
#include "ParamText.h"

using namespace tt::ui;

namespace
{
    const char* const kOsNames[] = { "1X", "2X", "4X", "8X", "16X" };
    constexpr float kUiScales[] = { 0.75f, 0.9f, 1.0f, 1.25f, 1.5f, 2.0f };

    // ---- the mockup's geometry
    constexpr float R = 110.0f;
    const juce::Point<float> kReelL { 138.0f, 180.0f }, kReelR { 379.0f, 180.0f }, kVuIn { 621.0f, 180.0f }, kVuOut { 862.0f, 180.0f };
    constexpr float kPackL = 0.80f * R, kPackR = 0.55f * R;
    constexpr float kBoxW = 150.0f, kBoxH = 150.0f, kBoxY = 310.0f, kBoxCY = kBoxY + kBoxH / 2.0f;
    // middle row, three boxes of one size, each centred between two of the top row: tape select between the
    // reels, the badge between the take-up reel and the input meter (the faceplate's centre), speed select
    // between the meters
    constexpr float kSelX = (138.0f + 379.0f) / 2.0f - kBoxW / 2.0f;     // 183.5
    constexpr float kLogoX = (379.0f + 621.0f) / 2.0f - kBoxW / 2.0f;    // 425
    constexpr float kSpdX = (621.0f + 862.0f) / 2.0f - kBoxW / 2.0f;     // 666.5
    // the wordmark under the badge, set justified: the first two lines exactly as wide as the box above
    // them, "Special 1500" larger, 25 px past the box on each side (as in Tools/gui/gen_tt.py)
    constexpr float kScriptW = kBoxW + 50.0f;                            // 200
    constexpr float kScriptX = kLogoX + (kBoxW - kScriptW) / 2.0f;       // 400
    constexpr float kTextBaselines[] = { 484.0f, 504.0f, 546.0f };
    constexpr float kRowY = 548.0f;
    // bottom row: INPUT, AG beside it, OUTPUT
    constexpr float kSmallR = 28.0f, kSmallY = 572.0f;
    constexpr float kAgX = 238.0f;
    constexpr float kIps[] = { 3.75f, 7.5f, 15.0f };

    juce::RangedAudioParameter& param (TechnicalTapeProcessor& p, const char* id) { return *p.apvts.getParameter (id); }

    void setFromTyped (juce::RangedAudioParameter& p, float value)
    {
        value = p.getNormalisableRange().snapToLegalValue (value);
        p.beginChangeGesture();
        p.setValueNotifyingHost (p.convertTo0to1 (value));
        p.endChangeGesture();
    }

    /** Tangent point on circle (c, r) for a line from an outside point; side +1 picks the right-hand one. */
    juce::Point<float> tangentFrom (juce::Point<float> from, juce::Point<float> c, float r, float side)
    {
        const float dx = from.x - c.x, dy = from.y - c.y;
        const float d = std::hypot (dx, dy), phi = std::atan2 (dy, dx), a = std::acos (r / d);
        const juce::Point<float> p1 { c.x + r * std::cos (phi + a), c.y + r * std::sin (phi + a) };
        const juce::Point<float> p2 { c.x + r * std::cos (phi - a), c.y + r * std::sin (phi - a) };
        return side * p1.x >= side * p2.x ? p1 : p2;
    }
} // namespace

//==============================================================================
class TechnicalTapeEditor::Faceplate final : public juce::Component
{
public:
    Faceplate (TechnicalTapeProcessor& p, TechnicalTapeEditor& e) : proc (p), editor (e)
    {
        logo = pegasus (col::white);

        // ---- top strip
        bypass.setClickingTogglesState (true);
        addAndMakeVisible (bypass);
        bypass.setBounds (18, 9, 84, 24);
        bypassAtt = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment> (proc.apvts, ParamIDs::bypass, bypass);

        addAndMakeVisible (osPill);
        osPill.setBounds (680, 9, 56, 24);
        osPill.onClick = [this] { editor.showOsMenu(); };
        osAtt = std::make_unique<juce::ParameterAttachment> (param (proc, ParamIDs::oversampling), [this] (float v)
        {
            osPill.setLabel (kOsNames[juce::jlimit (0, 4, juce::roundToInt (v))]);
        });

        addAndMakeVisible (menuZone);
        addAndMakeVisible (saveZone);
        menuZone.setBounds (752, 9, 148, 24);
        saveZone.setBounds (900, 9, 82, 24);
        menuZone.onClick = [this] { editor.showPresetMenu(); };
        saveZone.onClick = [this] { editor.chooseAndSavePreset(); };

        // ---- reels (the supply reel follows the tape), meters
        addAndMakeVisible (reelL);
        addAndMakeVisible (reelR);
        reelL.setBounds (juce::Rectangle<float> (2.0f * R, 2.0f * R).withCentre (kReelL).toNearestInt());
        reelR.setBounds (juce::Rectangle<float> (2.0f * R, 2.0f * R).withCentre (kReelR).toNearestInt());
        addAndMakeVisible (vuInMeter);
        addAndMakeVisible (vuOutMeter);
        vuInMeter.setBounds (juce::Rectangle<float> (2.0f * R, 2.0f * R).withCentre (kVuIn).toNearestInt());
        vuOutMeter.setBounds (juce::Rectangle<float> (2.0f * R, 2.0f * R).withCentre (kVuOut).toNearestInt());

        // ---- tape select between the reels (the tape runs into it), speed select between the meters
        addAndMakeVisible (tapeSelect);
        addAndMakeVisible (speedSelect);
        tapeSelect.setBounds ((int) std::floor (kSelX), (int) kBoxY, (int) SelectBox::kSize + 1, (int) SelectBox::kSize);
        speedSelect.setBounds ((int) std::floor (kSpdX), (int) kBoxY, (int) SelectBox::kSize + 1, (int) SelectBox::kSize);
        tapeSelect.onPick = [this] (int i) { tapeAtt->setValueAsCompleteGesture ((float) i); };
        speedSelect.onPick = [this] (int i) { speedAtt->setValueAsCompleteGesture ((float) i); };
        // a tape without any model yet (Maxell until it is captured) is greyed out, and plays the other tape
        const auto& bank = proc.getModels();
        for (int t = 0; t < ModelBank::kTapes; ++t)
            tapeSelect.setRowEnabled (t, bank.tapeAvailable (t),
                                      bank.tapeAvailable (t) ? juce::String() : juce::String (ModelBank::kTapeNames[t]) + ": not captured yet");
        tapeAtt = std::make_unique<juce::ParameterAttachment> (param (proc, ParamIDs::tape), [this] (float v)
        {
            const int t = proc.getModels().effectiveTape (juce::roundToInt (v));
            tapeSelect.setSelected (t);
            if (t == 0) reelL.setLook (Reel::Style::ThreeWindows, Reel::skyBlue());
            else reelL.setLook (Reel::Style::TwoSlots, Reel::white());
        });
        speedAtt = std::make_unique<juce::ParameterAttachment> (param (proc, ParamIDs::speed), [this] (float v)
        {
            speedSelect.setSelected (juce::roundToInt (v));
        });

        // ---- AG beside INPUT
        addAndMakeVisible (ag);
        ag.setBounds (juce::Rectangle<float> (2.0f * kSmallR, 2.0f * kSmallR).withCentre ({ kAgX, kSmallY }).toNearestInt());
        agAtt = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment> (proc.apvts, ParamIDs::autoGain, ag);
        ag.setTooltip ("Auto gain: keeps the loudness level as INPUT drives the tape harder");

        // ---- INPUT under the supply reel, OUTPUT under the output meter
        for (auto* k : { &input, &output })
        {
            k->setInterceptsMouseClicks (false, true);
            addAndMakeVisible (k);
        }
        input.placeAt (kReelL.x, kRowY);
        output.placeAt (kVuOut.x, kRowY);
        attach (input, ParamIDs::input);
        attach (output, ParamIDs::output);

        addChildComponent (entry);
        addMouseListener (&clickAway, true);   // a click anywhere else applies an open type-in box

        for (auto* a : { osAtt.get(), tapeAtt.get(), speedAtt.get() })
            a->sendInitialUpdate();
    }

    ~Faceplate() override { removeMouseListener (&clickAway); }

    void paint (juce::Graphics& g) override
    {
        g.fillAll (col::salmon);
        g.setColour (col::yellow);
        g.fillRect (0, 0, kWidth, 42);
        drawText (g, "OS", 662.0f, 26.0f, font (Weight::SemiBold, 13.0f, 1.6f), col::ink);
        g.setColour (col::white);
        g.fillRoundedRectangle (752.0f, 9.0f, 230.0f, 24.0f, 12.0f);
        juce::Path divider;
        divider.startNewSubPath (900.0f, 14.0f);
        divider.lineTo (900.0f, 28.0f);
        g.setColour (col::ink.withAlpha (0.3f));
        strokeRound (g, divider, 1.5f);

        // tape: off the supply pack's right side and the take-up pack's left side, down into the tape select box
        const juce::Point<float> e1 { kSelX + kBoxW / 2.0f - 25.0f, kBoxY + 2.0f }, e2 { kSelX + kBoxW / 2.0f + 25.0f, kBoxY + 2.0f };
        const auto t1 = tangentFrom (e1, kReelL, kPackL, 1.0f), t2 = tangentFrom (e2, kReelR, kPackR, -1.0f);
        juce::Path tape;
        tape.startNewSubPath (t1);
        tape.lineTo (e1);
        tape.startNewSubPath (e2);
        tape.lineTo (t2);
        g.setColour (col::ink);
        strokeRound (g, tape, 3.0f);

        // the badge box in the middle of the faceplate
        const juce::Rectangle<float> box (kLogoX + 1.5f, kBoxY + 1.5f, kBoxW - 3.0f, kBoxH - 3.0f);
        g.setColour (col::salmon);
        g.fillRoundedRectangle (box, 11.0f);
        g.setColour (col::white);
        g.drawRoundedRectangle (box, 11.0f, 3.0f);
        if (logo != nullptr)
        {
            const float lw = 118.0f * 1218.0f / 1249.0f;
            logo->drawWithin (g, { kLogoX + kBoxW / 2.0f - lw / 2.0f, kBoxCY - 59.0f, lw, 118.0f },
                              juce::RectanglePlacement::centred, 1.0f);
        }

        // wordmark under the badge, set justified: each line's ends on the column's edges (a magazine column)
        const juce::String l1 ("Sound Aquarium Recording"), l2 ("TECHNICAL TAPE"), l3 ("Special 1500");
        drawJustifiedLine (g, l1, kLogoX, kTextBaselines[0], kBoxW, false, Weight::Medium, 11.13f, col::white);
        drawJustifiedLine (g, l2, kLogoX, kTextBaselines[1], kBoxW, false, Weight::Bold, 16.26f, col::white);
        drawJustifiedLine (g, l3, kScriptX, kTextBaselines[2], kScriptW, true, Weight::Medium, 0.0f, col::white);
    }

    void attach (ParamKnob& k, const char* id)
    {
        sliderAtts.add (new juce::AudioProcessorValueTreeState::SliderAttachment (proc.apvts, id, k.slider));
        auto& p = param (proc, id);
        k.slider.setDoubleClickReturnValue (true, (double) p.convertFrom0to1 (p.getDefaultValue()));
        k.slider.onDoubleClick = [this, &k, &p] { openEntry (p, k.getKnobCentre(), 90); };
    }

    /** Type-in box centred on `at` (faceplate coordinates); the typed value is clamped to the parameter. */
    void openEntry (juce::RangedAudioParameter& p, juce::Point<float> at, int width)
    {
        constexpr int h = 30;
        auto area = juce::Rectangle<int> (width, h).withCentre (at.roundToInt());
        area = area.constrainedWithin ({ 4, 4, kWidth - 8, kHeight - 8 });
        entry.open (area, p.getCurrentValueAsText(), [&p] (const juce::String& text)
        {
            if (auto v = tt::parseTypedValue (text, false))
                setFromTyped (p, *v);
        });
    }

    /** Turns the reels: `spin` is 0 (stopped) .. 1 (running at speed). The take-up pack is smaller, so
        that reel turns faster. */
    void advanceReels (float dtSeconds, float spin)
    {
        if (spin <= 1e-4f) return;
        const int s = juce::jlimit (0, 2, juce::roundToInt (proc.apvts.getRawParameterValue (ParamIDs::speed)->load()));
        const float surface = 845.0f * kIps[s] * spin;       // 15 ips: supply reel at 0.4 rev/s
        reelL.setTurn (reelL.getTurn() + dtSeconds * surface / kPackL);
        reelR.setTurn (reelR.getTurn() + dtSeconds * surface / kPackR);
    }

    TechnicalTapeProcessor& proc;
    TechnicalTapeEditor& editor;
    std::unique_ptr<juce::Drawable> logo;

    PillButton bypass { "BYPASS" }, osPill { "8X", 13.0f, 0.6f, Weight::Bold };
    ZoneButton menuZone { "MENU", true, true }, saveZone { "SAVE", false, false };
    Reel reelL { Reel::Style::ThreeWindows, Reel::skyBlue(), kPackL, 12.0f };
    Reel reelR { Reel::Style::ThreeWindows, Reel::white(), kPackR, 47.0f };
    TapeVu vuInMeter { "INPUT VU" }, vuOutMeter { "OUTPUT VU" };
    SelectBox tapeSelect { "TAPE SELECT", { "ATR", "MAXELL" }, 20.0f, 2.5f };
    SelectBox speedSelect { "SPEED SELECT", { juce::String::fromUTF8 ("3\xc2\xbe IPS"), juce::String::fromUTF8 ("7\xc2\xbd IPS"), "15 IPS" },
                            18.0f, 2.0f };
    RoundButton ag { "AG", true };
    ParamKnob input { "INPUT", 46.0f, 6.0f, true, 86.0f, 13.0f, 2.6f };
    ParamKnob output { "OUTPUT", 46.0f, 6.0f, true, 86.0f, 13.0f, 2.6f };   // the same size as INPUT
    ValueEntry entry;

    struct ClickAway final : juce::MouseListener
    {
        explicit ClickAway (ValueEntry& e) : entry (e) {}
        void mouseDown (const juce::MouseEvent& e) override
        {
            if (entry.isVisible() && e.eventComponent != &entry && ! entry.isParentOf (e.eventComponent))
                entry.dismiss (true);
        }
        ValueEntry& entry;
    } clickAway { entry };

    juce::OwnedArray<juce::AudioProcessorValueTreeState::SliderAttachment> sliderAtts;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> bypassAtt, agAtt;
    std::unique_ptr<juce::ParameterAttachment> osAtt, tapeAtt, speedAtt;
};

//==============================================================================
TechnicalTapeEditor::TechnicalTapeEditor (TechnicalTapeProcessor& p) : juce::AudioProcessorEditor (p), proc (p)
{
    setLookAndFeel (&lnf);
    face = std::make_unique<Faceplate> (p, *this);
    addAndMakeVisible (*face);
    tooltips.toFront (false);     // over the faceplate
    face->setBounds (0, 0, kWidth, kHeight);
    applySize();
    lastTick = juce::Time::getMillisecondCounterHiRes();
    startTimerHz (30);
}

TechnicalTapeEditor::~TechnicalTapeEditor()
{
    stopTimer();
    fileChooser.reset();
    face.reset();
    setLookAndFeel (nullptr);
}

void TechnicalTapeEditor::applySize()
{
    const float s = juce::jlimit (0.5f, 3.0f, proc.uiScale.load());
    face->setTransform (juce::AffineTransform::scale (s));
    setSize (juce::roundToInt ((float) kWidth * s), juce::roundToInt ((float) kHeight * s));
}

void TechnicalTapeEditor::timerCallback()
{
    const double now = juce::Time::getMillisecondCounterHiRes();
    const float dt = (float) juce::jlimit (0.0, 0.2, (now - lastTick) * 0.001);
    lastTick = now;

    // VU needles: the engine averages with VU timing; the needle adds a little mechanical inertia
    auto toVu = [] (float ms) { return 10.0f * std::log10 (std::max (ms, 1e-12f) / TechnicalTapeProcessor::kZeroVuMeanSquare); };
    const float k = 1.0f - std::exp (-dt / 0.045f);
    vuIn += (std::max (-40.0f, toVu (proc.getEngine().inputMeanSquare())) - vuIn) * k;
    vuOut += (std::max (-40.0f, toVu (proc.getEngine().outputMeanSquare())) - vuOut) * k;
    face->vuInMeter.setVu (vuIn);
    face->vuOutMeter.setVu (vuOut);

    // reels: spin up and run down like a transport (~0.4 s)
    const float target = proc.isRolling() ? 1.0f : 0.0f;
    spin += (target - spin) * (1.0f - std::exp (-dt / 0.4f));
    face->advanceReels (dt, spin);

    if (savedFlashTicks > 0 && --savedFlashTicks == 0)
        face->saveZone.setButtonText ("SAVE");
    refreshPresetLabel();
}

void TechnicalTapeEditor::showOsMenu()
{
    juce::PopupMenu m;
    m.addSectionHeader ("OVERSAMPLING");
    const int current = juce::roundToInt (proc.apvts.getRawParameterValue (ParamIDs::oversampling)->load());
    for (int i = 0; i < 5; ++i)
        m.addItem (i + 1, kOsNames[i], true, i == current);
    m.addSeparator();
    m.addItem (-1, "8X and 16X are alias-free for every captured harmonic", false, false);
    m.showMenuAsync (juce::PopupMenu::Options().withTargetComponent (&face->osPill),
                     [this] (int r)
                     {
                         if (r <= 0) return;
                         setFromTyped (param (proc, ParamIDs::oversampling), (float) (r - 1));
                     });
}

void TechnicalTapeEditor::showPresetMenu()
{
    const auto files = proc.presets.list();
    const auto current = proc.presets.currentFile;
    juce::PopupMenu m;
    m.addSectionHeader ("PRESETS");
    if (files.isEmpty())
        m.addItem (-1, "No presets saved yet", false, false);
    for (int i = 0; i < files.size(); ++i)
        m.addItem (1000 + i, PresetManager::displayName (files[i]), true, files[i] == current);
    m.addSeparator();
    m.addItem (1, "Previous preset", ! files.isEmpty());
    m.addItem (2, "Next preset", ! files.isEmpty());
    m.addItem (3, "Save preset...");
    m.addItem (7, "Load preset file...");
    m.addItem (4, "Delete \"" + proc.presets.currentName + "\"", current.existsAsFile());
    m.addItem (5, "Reset to default");
    m.addItem (6, "Show presets folder");
    m.addSeparator();

    juce::PopupMenu wf;
    const float wfNow = proc.apvts.getRawParameterValue (ParamIDs::wowFlutter)->load();
    for (int pc : { 0, 50, 100, 150, 200 })
        wf.addItem (300 + pc, pc == 0 ? juce::String ("Off") : juce::String (pc) + " %" + (pc == 100 ? " (as measured)" : ""),
                    true, std::abs (wfNow - (float) pc) < 0.5f);
    m.addSubMenu ("Wow & flutter", wf);

    juce::PopupMenu sizes;
    for (int i = 0; i < (int) std::size (kUiScales); ++i)
        sizes.addItem (200 + i, juce::String (juce::roundToInt (kUiScales[i] * 100.0f)) + " %", true,
                       std::abs (proc.uiScale.load() - kUiScales[i]) < 0.01f);
    m.addSubMenu ("Interface size", sizes);

    juce::PopupMenu info;
    const auto& bank = proc.getModels();
    for (int i = 0; i < ModelBank::kCount; ++i)
    {
        const auto& s = bank.slot (i);
        auto line = juce::String (ModelBank::kTapeNames[i / ModelBank::kSpeeds]) + " " + ModelBank::kSpeedNames[i % ModelBank::kSpeeds]
                    + " ips: ";
        if (s.present) line += (s.simulated ? "simulated machine" : "captured") + juce::String (s.origin == "built in" ? "" : " (your file)");
        else if (! bank.tapeAvailable (i / ModelBank::kSpeeds)) line += "not captured yet";
        else if (bank.resolved (i) >= 0) line += "not captured yet, plays " + ModelBank::fileNameFor (bank.resolved (i)).upToLastOccurrenceOf (".", false, false);
        else line += "missing";
        info.addItem (-1, line, false, false);
    }
    info.addSeparator();
    info.addItem (8, "Show models folder");
    m.addSubMenu ("Tape models", info);

    m.showMenuAsync (juce::PopupMenu::Options().withTargetComponent (&face->menuZone),
                     [this, files, current] (int r)
                     {
                         if (r >= 1000 && r - 1000 < files.size()) proc.presets.load (files[r - 1000]);
                         else if (r >= 300 && r <= 500) setFromTyped (param (proc, ParamIDs::wowFlutter), (float) (r - 300));
                         else if (r >= 200 && r < 200 + (int) std::size (kUiScales)) { proc.uiScale.store (kUiScales[r - 200]); applySize(); }
                         else if (r == 1) proc.presets.step (-1);
                         else if (r == 2) proc.presets.step (1);
                         else if (r == 3) chooseAndSavePreset();
                         else if (r == 4) proc.presets.remove (current);
                         else if (r == 5) proc.presets.resetToDefault();
                         else if (r == 6) PresetManager::getFolder().startAsProcess();
                         else if (r == 7) chooseAndLoadPreset();
                         else if (r == 8)
                         {
                             const auto dir = ModelBank::userFolder();
                             dir.createDirectory();
                             dir.startAsProcess();
                         }
                         refreshPresetLabel();
                     });
}

void TechnicalTapeEditor::chooseAndSavePreset()
{
    const auto current = proc.presets.currentFile;
    const auto start = current.existsAsFile() ? current
                                              : PresetManager::getFolder().getChildFile (PresetManager::suggestName() + PresetManager::kExtension);
    fileChooser = std::make_unique<juce::FileChooser> ("Save preset", start, juce::String ("*") + PresetManager::kExtension);
    fileChooser->launchAsync (juce::FileBrowserComponent::saveMode | juce::FileBrowserComponent::canSelectFiles
                                  | juce::FileBrowserComponent::warnAboutOverwriting,
                              [safe = juce::Component::SafePointer<TechnicalTapeEditor> (this)] (const juce::FileChooser& fc)
                              {
                                  const auto file = fc.getResult();
                                  if (safe != nullptr && file != juce::File())
                                      safe->savePresetTo (file);
                              });
}

void TechnicalTapeEditor::chooseAndLoadPreset()
{
    fileChooser = std::make_unique<juce::FileChooser> ("Load preset", PresetManager::getFolder(), juce::String ("*") + PresetManager::kExtension);
    fileChooser->launchAsync (juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
                              [safe = juce::Component::SafePointer<TechnicalTapeEditor> (this)] (const juce::FileChooser& fc)
                              {
                                  const auto file = fc.getResult();
                                  if (safe == nullptr || file == juce::File()) return;
                                  if (! safe->proc.presets.load (file))
                                      safe->showProblem ("Preset not loaded", file.getFileName() + " isn't a Technical Tape preset.");
                                  safe->refreshPresetLabel();
                              });
}

bool TechnicalTapeEditor::savePresetTo (const juce::File& file)
{
    const auto result = proc.presets.saveTo (file);
    if (result.failed())
    {
        showProblem ("Preset not saved", result.getErrorMessage());
        return false;
    }
    face->saveZone.setButtonText ("SAVED");
    savedFlashTicks = 45; // 1.5 s
    refreshPresetLabel();
    return true;
}

void TechnicalTapeEditor::showProblem (const juce::String& title, const juce::String& message)
{
    juce::NativeMessageBox::showForOptionsAsync (juce::MessageBoxOptions()
                                                     .withIconType (juce::MessageBoxIconType::WarningIcon)
                                                     .withTitle (title)
                                                     .withMessage (message)
                                                     .withButton ("OK", 1)
                                                     .withAssociatedComponent (this),
                                                 std::function<void (int)> {});
}

void TechnicalTapeEditor::refreshPresetLabel()
{
    const auto& f = proc.presets.currentFile;
    const auto label = f != juce::File() ? f.getFileNameWithoutExtension().toUpperCase() : juce::String ("MENU");
    if (face->menuZone.getButtonText() != label)
        face->menuZone.setButtonText (label);
}

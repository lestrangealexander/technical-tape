// Developer tool: render the Technical Tape editor to PNG with live audio running through it, and
// click its controls the way a user would.
//   TTSnapshot <output folder>
#include "PluginEditor.h"
#include <cmath>

namespace
{
template <typename T>
T* find (juce::Component& root, std::function<bool (T&)> pred = {})
{
    for (auto* c : root.getChildren())
    {
        if (auto* t = dynamic_cast<T*> (c))
            if (pred == nullptr || pred (*t))
                return t;
        if (auto* deeper = find<T> (*c, pred))
            return deeper;
    }
    return nullptr;
}

int failures = 0;
void check (bool ok, const juce::String& what)
{
    std::printf ("  %s  %s\n", ok ? "ok  " : "FAIL", what.toRawUTF8());
    if (! ok) ++failures;
}

void click (juce::Component& c)
{
    // what a mouse click does to a button (onClick, and the toggle for toggling buttons); delivered
    // through the message loop like a real one
    if (auto* b = dynamic_cast<juce::Button*> (&c))
        b->triggerClick();
}
} // namespace

int main (int argc, char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;
    const juce::File outDir (argc > 1 ? juce::String (argv[1]) : juce::File::getCurrentWorkingDirectory().getFullPathName());
    outDir.createDirectory();

    TechnicalTapeProcessor proc;
    proc.setPlayConfigDetails (2, 2, 48000.0, 512);
    const auto t0 = juce::Time::getMillisecondCounterHiRes();
    proc.prepareToPlay (48000.0, 512);
    std::printf ("prepareToPlay: %.0f ms, latency %d samples\n", juce::Time::getMillisecondCounterHiRes() - t0, proc.getLatencySamples());
    auto get = [&proc] (const char* id) { return proc.apvts.getRawParameterValue (id)->load(); };

    // programme: a drum-ish pulse over a chord, around 0 VU
    juce::AudioBuffer<float> buf (2, 512);
    juce::MidiBuffer midi;
    double ph[3] = { 0, 0, 0 };
    const double fr[3] = { 110.0, 277.2, 659.3 };
    long n = 0;
    auto runAudio = [&] (int blocks)
    {
        for (int b = 0; b < blocks; ++b)
        {
            for (int i = 0; i < 512; ++i, ++n)
            {
                float x = 0.0f;
                for (int k = 0; k < 3; ++k)
                {
                    ph[k] += 2.0 * juce::MathConstants<double>::pi * fr[k] / 48000.0;
                    x += 0.07f * (float) std::sin (ph[k]);
                }
                const double beat = std::fmod ((double) n / 48000.0, 0.5);
                x *= (float) (0.6 + 0.4 * std::exp (-beat * 8.0));
                buf.setSample (0, i, x);
                buf.setSample (1, i, x);
            }
            proc.processBlock (buf, midi);
        }
    };
    auto pump = [] (int ms) { juce::MessageManager::getInstance()->runDispatchLoopUntil (ms); };

    std::unique_ptr<juce::AudioProcessorEditor> ed (proc.createEditor());
    ed->setVisible (true);
    auto shoot = [&] (const juce::String& name, int rounds = 40)
    {
        for (int i = 0; i < rounds; ++i) { runAudio (3); pump (15); }
        auto img = ed->createComponentSnapshot (ed->getLocalBounds(), true, 1.0f);
        juce::PNGImageFormat png;
        auto f = outDir.getChildFile (name);
        f.deleteFile();
        juce::FileOutputStream os (f);
        png.writeImageToStream (img, os);
        std::printf ("wrote %s\n", f.getFullPathName().toRawUTF8());
    };

    std::printf ("Technical Tape editor checks\n");
    check (ed->getWidth() == 1000 && ed->getHeight() == 660, "editor is 1000 x 660 at 100 %");
    shoot ("tt_atr_15.png");

    auto round = [&] (const juce::String& label)
    {
        return find<tt::ui::RoundButton> (*ed, [&label] (tt::ui::RoundButton& b) { return b.getButtonText() == label; });
    };
    auto box = [&] (const juce::String& header)
    {
        return find<tt::ui::SelectBox> (*ed, [&header] (tt::ui::SelectBox& b) { return b.getHeader() == header; });
    };
    auto* tapes = box ("TAPE SELECT");
    auto* speeds = box ("SPEED SELECT");
    check (tapes != nullptr && tapes->getSelected() == 0 && tapes->getNumChoices() == 2, "TAPE SELECT shows ATR / MAXELL, ATR chosen");
    check (speeds != nullptr && speeds->getSelected() == 2 && speeds->getNumChoices() == 3, "SPEED SELECT shows three speeds, 15 chosen");
    check (find<tt::ui::RoundButton> (*ed, [] (tt::ui::RoundButton& b) { return b.getButtonText() != "AG"; }) == nullptr,
           "the bottom row has only INPUT, AG and OUTPUT");
    if (proc.getModels().tapeAvailable (1))
    {
        if (tapes != nullptr) { tapes->pick (1); pump (30); }
        check (juce::roundToInt (get (ParamIDs::tape)) == 1 && tapes != nullptr && tapes->getSelected() == 1, "MAXELL selects the Maxell tape");
    }
    else
    {
        // not captured yet: the row is greyed out, a click does nothing, and a host or preset asking for
        // Maxell keeps ATR lit (and playing)
        check (tapes != nullptr && ! tapes->isRowEnabled (1) && tapes->isRowEnabled (0), "MAXELL is greyed out until it is captured");
        if (tapes != nullptr) { tapes->pick (1); pump (30); }
        check (juce::roundToInt (get (ParamIDs::tape)) == 0 && tapes != nullptr && tapes->getSelected() == 0, "a click on MAXELL is ignored");
        auto* tp = proc.apvts.getParameter (ParamIDs::tape);
        tp->setValueNotifyingHost (tp->convertTo0to1 (1.0f));
        pump (30);
        check (tapes != nullptr && tapes->getSelected() == 0, "Maxell from the host shows (and plays) ATR");
        tp->setValueNotifyingHost (tp->convertTo0to1 (0.0f));
        pump (30);
    }
    if (speeds != nullptr) { speeds->pick (1); pump (30); }
    check (juce::roundToInt (get (ParamIDs::speed)) == 1 && speeds != nullptr && speeds->getSelected() == 1,
           juce::String::fromUTF8 ("7\xc2\xbd IPS selects 7.5 ips"));
    shoot (proc.getModels().tapeAvailable (1) ? "tt_maxell_7_5.png" : "tt_atr_7_5.png");

    if (auto* ag = round ("AG")) { click (*ag); pump (30); }
    check (get (ParamIDs::autoGain) < 0.5f, "AG turns auto gain off");
    if (speeds != nullptr) { speeds->pick (0); pump (30); }
    check (juce::roundToInt (get (ParamIDs::speed)) == 0, juce::String::fromUTF8 ("3\xc2\xbe IPS selects 3.75 ips"));
    if (tapes != nullptr) { tapes->pick (0); pump (30); }
    check (juce::roundToInt (get (ParamIDs::tape)) == 0, "ATR selects the ATR tape");
    auto* in = proc.apvts.getParameter (ParamIDs::input);
    in->setValueNotifyingHost (in->convertTo0to1 (9.0f));
    shoot ("tt_atr_375_ag_off_drive.png");

    auto* bypass = find<tt::ui::PillButton> (*ed, [] (tt::ui::PillButton& p) { return p.getButtonText() == "BYPASS"; });
    if (bypass != nullptr) { click (*bypass); pump (30); }
    check (get (ParamIDs::bypass) > 0.5f, "BYPASS engages");
    shoot ("tt_bypass.png", 10);

    proc.uiScale.store (1.5f);
    ed.reset (proc.createEditor());
    check (ed->getWidth() == 1500 && ed->getHeight() == 990, "interface size 150 %");
    shoot ("tt_150.png", 5);
    ed.reset();

    std::printf ("%s (%d failure%s)\n", failures == 0 ? "all good" : "PROBLEMS", failures, failures == 1 ? "" : "s");
    return failures == 0 ? 0 : 1;
}

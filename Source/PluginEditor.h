// Technical Tape — editor (the Special 1500 faceplate)
#pragma once

#include "PluginProcessor.h"
#include "ui/Widgets.h"

class TechnicalTapeEditor final : public juce::AudioProcessorEditor, private juce::Timer
{
public:
    explicit TechnicalTapeEditor (TechnicalTapeProcessor&);
    ~TechnicalTapeEditor() override;

    void paint (juce::Graphics&) override {}
    void resized() override {}

    static constexpr int kWidth = 1000, kHeight = 660;

    /** Saves the current settings as a preset file (what SAVE does once a file is chosen). */
    bool savePresetTo (const juce::File& file);

private:
    class Faceplate;

    void timerCallback() override;
    void applySize();
    void showPresetMenu();
    void showOsMenu();
    void chooseAndSavePreset();
    void chooseAndLoadPreset();
    void showProblem (const juce::String& title, const juce::String& message);
    void refreshPresetLabel();

    TechnicalTapeProcessor& proc;
    tt::ui::TapeLookAndFeel lnf;
    std::unique_ptr<Faceplate> face;
    juce::TooltipWindow tooltips { this, 600 };
    int savedFlashTicks = 0;
    std::unique_ptr<juce::FileChooser> fileChooser;
    float vuIn = -60.0f, vuOut = -60.0f, spin = 0.0f;
    double lastTick = 0.0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TechnicalTapeEditor)
};

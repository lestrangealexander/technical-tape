// Technical Tape — Sound Aquarium
#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include "ModelBank.h"
#include "PresetManager.h"
#include "engine/TapeEngine.h"

namespace ParamIDs
{
inline constexpr const char* tape = "tape";                // ATR, Maxell
inline constexpr const char* speed = "speed";              // 3.75, 7.5, 15 ips
inline constexpr const char* input = "input";
inline constexpr const char* output = "output";
inline constexpr const char* autoGain = "autoGain";
inline constexpr const char* bypass = "bypass";
inline constexpr const char* oversampling = "oversampling";
inline constexpr const char* wowFlutter = "wowFlutter";    // % of the measured amount
} // namespace ParamIDs

class TechnicalTapeProcessor final : public juce::AudioProcessor
{
public:
    TechnicalTapeProcessor();
    ~TechnicalTapeProcessor() override = default;

    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override {}
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    using AudioProcessor::processBlock;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

    const juce::String getName() const override { return "Technical Tape"; }
    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }
    bool isMidiEffect() const override { return false; }
    double getTailLengthSeconds() const override { return 0.25; }

    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram (int) override {}
    const juce::String getProgramName (int) override { return "Default"; }
    void changeProgramName (int, const juce::String&) override {}

    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    juce::AudioProcessorParameter* getBypassParameter() const override { return bypassParam; }

    static juce::AudioProcessorValueTreeState::ParameterLayout createLayout();

    juce::AudioProcessorValueTreeState apvts;
    PresetManager presets;
    std::atomic<float> uiScale { 1.0f };

    const ModelBank& getModels() const noexcept { return *models; }
    tt::TapeEngine& getEngine() noexcept { return engine; }
    /** True while the host's transport runs (or, without one, while audio is coming in). */
    bool isRolling() const noexcept { return rolling.load(); }

    /** 0 VU = -18 dBFS RMS (a sine at -15 dBFS peak), as in the captures. */
    static constexpr float kZeroVuMeanSquare = 0.015848932f;

private:
    std::shared_ptr<const ModelBank> models;
    tt::TapeEngine engine;      // after the models: it refers to them
    juce::AudioProcessorParameter* bypassParam = nullptr;
    std::atomic<bool> rolling { false };
    int silentBlocks = 0;

    std::atomic<float>* pTape = nullptr;   std::atomic<float>* pSpeed = nullptr;
    std::atomic<float>* pInput = nullptr;  std::atomic<float>* pOutput = nullptr;
    std::atomic<float>* pAutoGain = nullptr; std::atomic<float>* pBypass = nullptr;
    std::atomic<float>* pOs = nullptr;     std::atomic<float>* pWf = nullptr;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TechnicalTapeProcessor)
};

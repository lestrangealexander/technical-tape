// Technical Tape — Sound Aquarium
#include "PluginProcessor.h"
#include "PluginEditor.h"

namespace
{
juce::String dbText (float v, int)
{
    return (v > 0.04f ? "+" : "") + juce::String (std::abs (v) < 0.05f ? 0.0f : v, 1) + " dB";
}

std::optional<float> parseNumber (juce::String t)
{
    t = t.trim().toLowerCase().replace (juce::String::fromUTF8 ("\xe2\x88\x92"), "-").replace (",", ".");
    const auto num = t.retainCharacters ("0123456789.-+");
    if (num.isEmpty() || ! num.containsAnyOf ("0123456789")) return std::nullopt;
    return num.getFloatValue();
}
} // namespace

juce::AudioProcessorValueTreeState::ParameterLayout TechnicalTapeProcessor::createLayout()
{
    using namespace juce;
    AudioProcessorValueTreeState::ParameterLayout layout;
    auto db = [] (const char* id, const char* name)
    {
        return std::make_unique<AudioParameterFloat> (ParameterID { id, 1 }, name, NormalisableRange<float> { -18.0f, 18.0f, 0.1f }, 0.0f,
                                                      AudioParameterFloatAttributes().withStringFromValueFunction (dbText)
                                                          .withValueFromStringFunction ([] (const String& s) { return parseNumber (s).value_or (0.0f); }));
    };
    layout.add (std::make_unique<AudioParameterChoice> (ParameterID { ParamIDs::tape, 1 }, "Tape", StringArray { "ATR", "Maxell" }, 0));
    layout.add (std::make_unique<AudioParameterChoice> (ParameterID { ParamIDs::speed, 1 }, "Speed",
                                                        StringArray { String::fromUTF8 ("3\xc2\xbe ips"), String::fromUTF8 ("7\xc2\xbd ips"), "15 ips" }, 2));
    layout.add (db (ParamIDs::input, "Input"));
    layout.add (db (ParamIDs::output, "Output"));
    layout.add (std::make_unique<AudioParameterBool> (ParameterID { ParamIDs::autoGain, 1 }, "Auto Gain", true));
    layout.add (std::make_unique<AudioParameterBool> (ParameterID { ParamIDs::bypass, 1 }, "Bypass", false));
    layout.add (std::make_unique<AudioParameterChoice> (ParameterID { ParamIDs::oversampling, 1 }, "Oversampling",
                                                        StringArray { "1x", "2x", "4x", "8x", "16x" }, 3));
    layout.add (std::make_unique<AudioParameterFloat> (ParameterID { ParamIDs::wowFlutter, 1 }, "Wow & Flutter",
                                                       NormalisableRange<float> { 0.0f, 200.0f, 1.0f }, 100.0f,
                                                       AudioParameterFloatAttributes()
                                                           .withStringFromValueFunction ([] (float v, int) { return String (roundToInt (v)) + " %"; })
                                                           .withValueFromStringFunction ([] (const String& s) { return parseNumber (s).value_or (100.0f); })));
    return layout;
}

TechnicalTapeProcessor::TechnicalTapeProcessor()
    : AudioProcessor (BusesProperties()
                          .withInput ("Input", juce::AudioChannelSet::stereo(), true)
                          .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      apvts (*this, nullptr, "TechnicalTapeState", createLayout()),
      presets (apvts)
{
    pTape = apvts.getRawParameterValue (ParamIDs::tape);
    pSpeed = apvts.getRawParameterValue (ParamIDs::speed);
    pInput = apvts.getRawParameterValue (ParamIDs::input);
    pOutput = apvts.getRawParameterValue (ParamIDs::output);
    pAutoGain = apvts.getRawParameterValue (ParamIDs::autoGain);
    pBypass = apvts.getRawParameterValue (ParamIDs::bypass);
    pOs = apvts.getRawParameterValue (ParamIDs::oversampling);
    pWf = apvts.getRawParameterValue (ParamIDs::wowFlutter);
    bypassParam = apvts.getParameter (ParamIDs::bypass);

    models = ModelBank::shared();
    engine.setModels (models->enginePointers());
}

bool TechnicalTapeProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    const auto in = layouts.getMainInputChannelSet(), out = layouts.getMainOutputChannelSet();
    if (in != out) return false;
    return out == juce::AudioChannelSet::mono() || out == juce::AudioChannelSet::stereo();
}

void TechnicalTapeProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    const int ch = std::max (1, std::min (2, getTotalNumOutputChannels()));
    tt::EngineParams p;
    p.model = ModelBank::indexOf ((int) std::lround (pTape->load()), (int) std::lround (pSpeed->load()));
    p.osIndex = (int) std::lround (pOs->load());
    p.inputDb = pInput->load();
    p.outputDb = pOutput->load();
    p.autoGain = pAutoGain->load() > 0.5f;
    p.bypass = pBypass->load() > 0.5f;
    p.wowFlutter = pWf->load() * 0.01f;
    engine.setParams (p);
    engine.prepare (sampleRate, std::max (32, samplesPerBlock), ch);
    setLatencySamples (engine.latencySamples());
}

void TechnicalTapeProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    const int numIn = getTotalNumInputChannels(), numOut = getTotalNumOutputChannels();
    for (int c = numIn; c < numOut; ++c) buffer.clear (c, 0, buffer.getNumSamples());

    tt::EngineParams p;
    p.model = ModelBank::indexOf ((int) std::lround (pTape->load()), (int) std::lround (pSpeed->load()));
    p.osIndex = (int) std::lround (pOs->load());
    p.inputDb = pInput->load();
    p.outputDb = pOutput->load();
    p.autoGain = pAutoGain->load() > 0.5f;
    p.bypass = pBypass->load() > 0.5f;
    p.wowFlutter = pWf->load() * 0.01f;
    engine.setParams (p);

    const int nCh = std::min (2, buffer.getNumChannels());
    engine.process (buffer.getArrayOfWritePointers(), nCh, buffer.getNumSamples());

    // the reels turn while the transport runs; without transport information, while there is signal
    bool run = false;
    bool known = false;
    if (auto* ph = getPlayHead())
        if (auto pos = ph->getPosition())
        {
            known = true;
            run = pos->getIsPlaying() || pos->getIsRecording();
        }
    if (! known)
    {
        silentBlocks = engine.inputMeanSquare() > 1e-9f ? 0 : silentBlocks + 1;
        run = silentBlocks < 50;
    }
    rolling.store (run && ! p.bypass);
}

void TechnicalTapeProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    auto st = apvts.copyState();
    st.setProperty ("presetName", presets.currentName, nullptr);
    st.setProperty ("presetFile", presets.currentFile.getFullPathName(), nullptr);
    st.setProperty ("uiScale", uiScale.load(), nullptr);
    if (auto xml = st.createXml()) copyXmlToBinary (*xml, destData);
}

void TechnicalTapeProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    if (auto xml = getXmlFromBinary (data, sizeInBytes))
        if (xml->hasTagName (apvts.state.getType()))
        {
            auto tree = juce::ValueTree::fromXml (*xml);
            presets.currentName = tree.getProperty ("presetName", "default").toString();
            const auto path = tree.getProperty ("presetFile", {}).toString();
            presets.currentFile = juce::File::isAbsolutePath (path) ? juce::File (path) : juce::File();
            uiScale.store ((float) (double) tree.getProperty ("uiScale", 1.0));
            apvts.replaceState (tree);
        }
}

juce::AudioProcessorEditor* TechnicalTapeProcessor::createEditor()
{
    return new TechnicalTapeEditor (*this);
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new TechnicalTapeProcessor();
}

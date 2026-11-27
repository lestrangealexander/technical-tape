// Technical Tape — the six tape models (2 tapes x 3 speeds)
//
// Each model ships inside the plugin. A model file placed in the user models folder replaces the shipped
// one with the same name, so new captures can be tried without rebuilding:
//   Mac      ~/Library/Application Support/Sound Aquarium/Technical Tape/Models/ATR_15.ttm
//   Windows  %APPDATA%\Sound Aquarium\Technical Tape\Models\ATR_15.ttm
#pragma once

#include <juce_core/juce_core.h>
#include "engine/TapeModel.h"

#include <array>
#include <memory>
#include <vector>

class ModelBank
{
public:
    static constexpr int kTapes = 2, kSpeeds = 3, kCount = kTapes * kSpeeds;
    static constexpr const char* kTapeNames[kTapes] = { "ATR", "Maxell" };
    static constexpr const char* kSpeedNames[kSpeeds] = { "3.75", "7.5", "15" };

    static int indexOf (int tape, int speed) noexcept { return tape * kSpeeds + speed; }
    static juce::String fileNameFor (int index) { return juce::String (kTapeNames[index / kSpeeds]) + "_" + kSpeedNames[index % kSpeeds] + ".ttm"; }
    static juce::File userFolder();

    /** Loads every model: the user folder's file if there is one, otherwise the shipped one. */
    void loadAll();

    /** One bank for the whole process: every instance of the plugin shares the models (and so the
        prepared kernels). Loaded when the first instance needs it, freed with the last. */
    static std::shared_ptr<const ModelBank> shared();

    struct Slot
    {
        tt::ModelData data;
        bool present = false;
        bool simulated = false;            // built from the virtual machine, not a real capture
        juce::String origin;               // "built in" or the file's path
        juce::String problem;              // why a user file was not used
    };
    const Slot& slot (int i) const { return slots[(size_t) i]; }

    /** Model pointers for the engine, one per slot: a missing model is stood in for by the nearest one
        (same tape, nearest speed; then the other tape), so a choice always plays something. */
    std::vector<const tt::ModelData*> enginePointers() const;
    bool anyPresent() const noexcept;
    /** The slot that actually plays for slot i. */
    int resolved (int i) const noexcept;
    /** A tape can be chosen once at least one of its speeds has a model (shipped or in the user folder).
        Until then its TAPE SELECT row is greyed out, and a host or preset that asks for it plays the
        other tape. */
    bool tapeAvailable (int tape) const noexcept;
    /** The tape that plays when `tape` is asked for. */
    int effectiveTape (int tape) const noexcept;

private:
    std::array<Slot, kCount> slots;
};

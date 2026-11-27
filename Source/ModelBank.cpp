// Technical Tape — the six tape models
#include "ModelBank.h"
#include "BinaryData.h"

#include <mutex>

juce::File ModelBank::userFolder()
{
   #if JUCE_MAC
    return juce::File::getSpecialLocation (juce::File::userHomeDirectory)
               .getChildFile ("Library/Application Support/Sound Aquarium/Technical Tape/Models");
   #else
    return juce::File::getSpecialLocation (juce::File::userApplicationDataDirectory)
               .getChildFile ("Sound Aquarium/Technical Tape/Models");
   #endif
}

namespace
{
    bool isSimulated (const tt::ModelData& m)
    {
        // the pipeline records its sources; models built from the virtual machine say so
        const auto meta = juce::String (m.metaJson).removeCharacters (" \t\r\n");
        return meta.contains ("\"simulated\":true");
    }
} // namespace

void ModelBank::loadAll()
{
    for (int i = 0; i < kCount; ++i)
    {
        auto& s = slots[(size_t) i];
        s = Slot();
        const auto name = fileNameFor (i);

        const auto user = userFolder().getChildFile (name);
        if (user.existsAsFile())
        {
            juce::MemoryBlock mb;
            if (user.loadFileAsData (mb))
            {
                const auto err = tt::parseModel (mb.getData(), mb.getSize(), s.data);
                if (err.empty())
                {
                    s.present = true;
                    s.origin = user.getFullPathName();
                    s.simulated = isSimulated (s.data);
                    continue;
                }
                s.problem = user.getFileName() + ": " + juce::String (err);
                s.data = tt::ModelData();
            }
        }

        for (int r = 0; r < BinaryData::namedResourceListSize; ++r)
        {
            if (juce::String (BinaryData::originalFilenames[r]) != name) continue;
            int size = 0;
            if (const char* d = BinaryData::getNamedResource (BinaryData::namedResourceList[r], size))
            {
                const auto err = tt::parseModel (d, (size_t) size, s.data);
                if (err.empty())
                {
                    s.present = true;
                    s.origin = "built in";
                    s.simulated = isSimulated (s.data);
                }
                else
                    s.problem = name + " (built in): " + juce::String (err);
            }
        }
    }
}

std::shared_ptr<const ModelBank> ModelBank::shared()
{
    static std::mutex mtx;
    static std::weak_ptr<const ModelBank> current;
    const std::lock_guard<std::mutex> lock (mtx);
    if (auto sp = current.lock())
        return sp;
    auto bank = std::make_shared<ModelBank>();
    bank->loadAll();
    current = bank;
    return bank;
}

bool ModelBank::anyPresent() const noexcept
{
    for (const auto& s : slots)
        if (s.present) return true;
    return false;
}

int ModelBank::resolved (int i) const noexcept
{
    if (slots[(size_t) i].present) return i;
    const int tape = i / kSpeeds, speed = i % kSpeeds;
    for (int t : { tape, 1 - tape })
        for (int d = 0; d < kSpeeds; ++d)
            for (int sp : { speed - d, speed + d })
                if (sp >= 0 && sp < kSpeeds && slots[(size_t) indexOf (t, sp)].present)
                    return indexOf (t, sp);
    return -1;
}

bool ModelBank::tapeAvailable (int tape) const noexcept
{
    if (tape < 0 || tape >= kTapes) return false;
    for (int sp = 0; sp < kSpeeds; ++sp)
        if (slots[(size_t) indexOf (tape, sp)].present) return true;
    return false;
}

int ModelBank::effectiveTape (int tape) const noexcept
{
    tape = juce::jlimit (0, kTapes - 1, tape);
    if (tapeAvailable (tape)) return tape;
    for (int t = 0; t < kTapes; ++t)
        if (tapeAvailable (t)) return t;
    return tape;
}

std::vector<const tt::ModelData*> ModelBank::enginePointers() const
{
    std::vector<const tt::ModelData*> out;
    if (! anyPresent()) return out;
    for (int i = 0; i < kCount; ++i)
        out.push_back (&slots[(size_t) resolved (i)].data);
    return out;
}

cmake_minimum_required(VERSION 3.22)

if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")
    set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64" CACHE STRING "Universal binary")
endif()

project(TechnicalTape VERSION 0.9.2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Windows: link the C++ runtime statically so the plugin needs no Visual C++ redistributable.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

option(TT_BUILD_TOOLS "Build the developer tools (snapshot, engine renderer)" OFF)

# ---------------------------------------------------------------- the engine (plain C++, no JUCE)
set(TT_ENGINE_SOURCES
    Source/engine/Fft.cpp
    Source/engine/Design.cpp
    Source/engine/Convolver.cpp
    Source/engine/TapeModel.cpp
    Source/engine/TapeEngine.cpp
    Source/engine/pffft/pffft.c
    Source/engine/pffft/pffft_common.c)

add_library(tt_engine STATIC ${TT_ENGINE_SOURCES})
target_include_directories(tt_engine PUBLIC Source/engine Source/engine/pffft)
# PFFFT is compiled into the engine: no dllimport/dllexport on Windows
target_compile_definitions(tt_engine PUBLIC PFFFT_STATIC_DEFINE)
set_target_properties(tt_engine PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
    target_compile_options(tt_engine PRIVATE /O2 /fp:precise)
    target_compile_definitions(tt_engine PRIVATE _USE_MATH_DEFINES)
else()
    target_compile_options(tt_engine PRIVATE -O3)
endif()

# offline renderer: tests, validation and benchmarks without a DAW
add_executable(tt_render Tests/tt_render.cpp)
target_link_libraries(tt_render PRIVATE tt_engine)

# ---------------------------------------------------------------- the plugin
if(DEFINED TT_JUCE_PATH)
    add_subdirectory(${TT_JUCE_PATH} JUCE EXCLUDE_FROM_ALL)
else()
    include(FetchContent)
    FetchContent_Declare(JUCE
        GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
        GIT_TAG 9.0.6
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(JUCE)
endif()

if(APPLE)
    set(TT_FORMATS AU VST3)
else()
    set(TT_FORMATS VST3)
endif()

# the six tape models (Tools/pipeline/extract.py writes them here)
file(GLOB TT_MODELS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/Resources/Models/*.ttm")

juce_add_binary_data(TTAssets SOURCES
    Resources/DMSans-Medium.ttf
    Resources/DMSans-SemiBold.ttf
    Resources/DMSans-Bold.ttf
    Resources/Yellowtail-Regular.ttf
    Resources/Pegasus.svg
    ${TT_MODELS})

juce_add_plugin(TechnicalTape
    COMPANY_NAME "Sound Aquarium"
    COMPANY_WEBSITE ""
    BUNDLE_ID com.soundaquarium.technicaltape
    PLUGIN_MANUFACTURER_CODE SqAu
    PLUGIN_CODE TtSp
    FORMATS ${TT_FORMATS}
    PRODUCT_NAME "Technical Tape"
    DESCRIPTION "Sampled tape machine: ATR and Maxell tape at 3.75, 7.5 and 15 ips"
    IS_SYNTH FALSE
    NEEDS_MIDI_INPUT FALSE
    NEEDS_MIDI_OUTPUT FALSE
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    AU_MAIN_TYPE kAudioUnitType_Effect
    VST3_CATEGORIES Fx Distortion
    COPY_PLUGIN_AFTER_BUILD FALSE)

set(TT_PLUGIN_SOURCES
    Source/PluginProcessor.cpp
    Source/PluginEditor.cpp
    Source/ModelBank.cpp
    Source/ui/Theme.cpp
    Source/ui/Widgets.cpp)

target_sources(TechnicalTape PRIVATE ${TT_PLUGIN_SOURCES})
target_include_directories(TechnicalTape PRIVATE Source)

if(MSVC)
    target_compile_options(TechnicalTape PRIVATE /utf-8 /bigobj)
endif()

target_compile_definitions(TechnicalTape PUBLIC
    JUCE_WEB_BROWSER=0
    JUCE_USE_CURL=0
    JUCE_VST3_CAN_REPLACE_VST2=0
    JUCE_DISPLAY_SPLASH_SCREEN=0)

target_link_libraries(TechnicalTape
    PRIVATE
        TTAssets
        tt_engine
        juce::juce_audio_utils
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags)

if(TT_BUILD_TOOLS)
    # renders the editor to PNG and drives its controls (run under xvfb-run on a headless box)
    juce_add_console_app(TTSnapshot PRODUCT_NAME "TTSnapshot")
    target_sources(TTSnapshot PRIVATE Tests/Snapshot.cpp ${TT_PLUGIN_SOURCES})
    target_include_directories(TTSnapshot PRIVATE Source)
    target_compile_definitions(TTSnapshot PRIVATE JUCE_WEB_BROWSER=0 JUCE_USE_CURL=0 JUCE_MODAL_LOOPS_PERMITTED=1)
    target_link_libraries(TTSnapshot PRIVATE TTAssets tt_engine juce::juce_audio_utils juce::juce_recommended_config_flags)
    if(MSVC)
        target_compile_options(TTSnapshot PRIVATE /utf-8 /bigobj)
    endif()
endif()

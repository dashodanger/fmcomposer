# MUDTracker

MUDTracker is a music and sound creation tool, forked from [FM Composer](https://web.archive.org/web/20180727172405/http://fmcomposer.org/en/),
featuring a custom FM synthesizer engine and a tracker-like interface.
Released with 150+ FM instruments and drums, ranging from synth to acoustic sounds, covering the whole MIDI instrument set.

# Features
- 6 operator, 24 voice polyphony FM sound engine with intuitive drag&drop interface
- Tracker-style sequencer
- Lots of effects available : vibrato, tremolo, arpeggio, pitch slides, real time FM parameter modification, loop points...
- MIDI integration : MIDI file import with partial XG/GS support, MIDI keyboard support
- MUS integration: DMX MUS file import
- WAV export
- Piano roll view

# Technical overview
GUI is written in C++ and using the SFML library to get low level access to the display, keyboard and mouse input. Each GUI element has its class. They are used by Views which are the actual program screens (Pattern screen, Instrument screen, Settings etc.). Views, for most of them, contain everything they need for providing the features to the user. This results in quite big classes, which are split in several files for readability.

The audio engine is written in pure C and was optimized for lowest CPU usage possible (cache-friendly structures, loops on small data sets, no function calls in critical parts, taking advantage of unsigned int wrapping for phase accumulators and so on). A part of the optimization relies on the use of FISTP instruction instead of the standard C ftol() for int/float conversions.

MUDTracker uses its own binary format, .MDTI for storing individual instruments, .MDTB for a single-file General MIDI compatible instrument set, and .MDTS for songs that contain embedded instruments.

# Compiling
For Unix builds, you'll need the following additional development libraries :
- SFML 2 (package names may vary; for Debian/Ubuntu this is likely `libsfml-dev`)
- ALSA (package names may vary; for Debian/Ubuntu this is likely `libasound2-dev`)

The program and necessary resources will be located in the `/bin` subdirectory of your build folder when compilation is complete.

On Unix, the program can also be installed via `cmake --install <build folder>` to an appropriate location. This usually requires root permissions.

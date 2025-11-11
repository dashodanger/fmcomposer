# User Manual

## Tutorial

Let's launch MUDTracker.exe. The window below appears, showing the song editor. Hover the screenshot with mouse to see help.

Click on the first column of a channel then press an alphabetical key to add a note. You can edit the Key/Note mapping in the config page if they don't suit to your keyboard.

Adding a note will also add the instrument number (it depends on the instrument you select in the instrument list), and a default volume. You can select/move all those data with the mouse, edit them with right click or double click and more, please see keyboard shortcuts.

Notes can be slurred by removing their instrument number.

To stop a note, press '='

The last column is reserved for effects. They can change the panning, delay a note, change the tempo and many other things, see effect list.

Try composing a simple tune. To play several notes at the same time, put them on different channels :

A C major chord.
Bored with the default piano sound ? Click on the Instrument tab . The Add button above the instrument list is now available :

Click on Add to add a second instrument into the list, then Load to open one of the instruments shipped with the program. You'll found them in the instruments folder.

See the next chapter to learn how to create your own sounds.

## Creating an FM Instrument

### How it Works

FM (Frequency Modulation) synthesis works by modifying a wave's frequency depending on the amplitude of another wave :

The Carrier get its frequency changed as the Modulator evolves.
To get interesting sounds, more than a single modulation is needed. FM Composer provides 6 of them, which can be configured in any imaginable way. They are called operators, each one can generate a basic shape : sine wave, square, triangle...

When they are at the bottommost position, you listen to their direct output : changing their volume impacts as you expect, the output volume. When they are above, their sound isn't directly heared : their output modulates the frequency of the operator below . The volume of those operators will impact the tone of the below ones. When several operators are next to each other, there is no frequency modulation thing, their output is simply summed like if they were separate sounds. The way the operators are arranged is called an algorithm.

### Editing the Algorithm

Drag&drop the blocks to arrange them in the way you like
When an operator is moved, its bottom links are destroyed
Right-click to create additional links between operators

### Editing the Parameters

Click on the Instruments tab. Don't be impressed by the amount of parameters. They are organised per operator. When your mouse hovers a parameter block, the corresponding operator is highlighted so you can see what you are actually editing :

| Parameter | Description |
|-----------|-------------|
| Volume envelope | Volume-over-time parameters. It's a classic ADSR envelope but with added delay, initial volume and hold (held time before decay) additional parameters. |
| Frequency | Frequency parameters. Can be a ratio calculated from the playing note, or a fixed frequency. Integer ratios are best for harmonious sounds. 1/4 tone and Fine parameters allow to break those integer values and create inharmonic sounds. |
| Waveform | Changes the waveform type and its offset. |
| Pitch envelope | Changes the way the frequency evolves over time. Initial freq/decay time defines a starting frequency and a time to reach the normal frequency. Release freq/release time is the same but triggered on note release. |
| Keyboard scalings | Allow to modify operator parameters depending on the note position on the keyboard. The higher the values, the higher effect strength will be applied. |
| LFO | Choose how the LFO affects the operator's frequency and amplitude |

### Making Your First Sound

Start by experimenting how FM works. Put the operators 1 and 6 on a bottom position, Drag&drop the 2, 3, 4 and 5 on the 6 then mute the 6th. This will mute all the operators above it so you can listen to Op1 :

Play a note on your keyboard, you should hear a simple sine wave sound which is operator 1's output. Play with the envelope parameters and the waveform choice. Now put an operator above it. Listen how the sound is modified as you change its volume and frequency multiplier. This is very basic 2-operators FM and you can already make a lot of sounds with it. But they won't be much refined, the 4 others operators are here so you can build more precise harmonics to get the sound you want.

Splitting a sound make it easier to program
Let's say you want to create a pan flute sound. Try to split it in several parts. You have the 'T' noisy attack sound, and the main sustained sound of the pipe. You can use two blocks of 3 operators to make this sound : the first 3 will do the attack sound, the 3 others the pipe sound. It's a simple way of creating sounds when you don't know where to start : split the complex final sound into separate simple ones.

**Hints:**
- To create noise, air-type sounds in FM you have to use the Feedback parameter on the operator #1 pushed to the max.
- To create a closed pipe sound (all closed pipes have Odd harmonics in their sound) use a carrier with ratio=1 and modulator with ratio=2

**Improving the pan flute sound**
By listening to a real pan flute, you'll probably notice that it's not just a noisy attack + a resonating sound. There is still an airy sound in the sustained part, and the attack has some harmonics that are not pure noise. Try to recreate those subtleties by tweaking your current sound.

A great way of understanding how to make sounds is to look at some of the built-in FM Composer sounds (/instruments/ folder). You can see how they are made and use them as a starting point if creating your sound from scratch seems too annoying.

### Making the sound you're thinking about

Creating a particular sound can be achieved easily after some practice. Try creating a rough approximation with few operators, then refine the sound by adding others operators and fiddling with multipliers and waveforms.

#### General sound making tips

- Don't use too high modulation values except if deliberate. This will lead to the typical unrefined sounds we find in some old games. Try playing with more subtle modulations.
- If a sound has aliasing artifacts as you play in the upper range, you need to turn down some modulators as the note is getting higher. This is done with two parameters : Center note and Volume scale. Setting a negative value for Volume will decrease the operator volume as the key played is higher than the Center note.
- Mute the operators you aren't working on. This is easily done by putting all unused ops above one that is muted. Then you can pick from this stack as you build your sound
- If you're stuck, try some algorithm presets. You may find new sounds and ideas this way
- Don't overuse the special waveforms, they sound harsh. Sine and Smooth Saw are the most useful ones. The others may be interesting for synths, percussions, or as way to add brightness to dull sounds.
- For ensemble type sounds you may need a chorus effect. It's done by detuning the operators (Detune parameter). Different types of choruses can be achieved depending if the detuned operators are linked together or not (classic chorus or a timbre-changing chorus)

## Create a Song

As seen in the tutorial tutoriel, adding notes is done by selecting an instrument in the list, clicking in the pattern then pressing a key. 24 channels are available, each one having a volume, panning and reverb. Shift+Mouse wheel allows to scroll horizontally so you can see all the channels. To move/swap channels, click on their number in red then drag&drop :

Right clicking on the pattern list shows several options for organizing, creating and removing them :

Moving a pattern is done, like the channels, by drag&drop.

Pattern editing functions are also available with the right click. You can transpose the selected notes, create fade in/outs, insert/delete rows, add effects...

## Shortcuts

### General
- CTRL+O : Open file
- CTRL+S : Save file
- Esc : Close popup
- Enter/Space bar : Play/Stop the song

### Song editor
- CTRL+Z : Undo
- CTRL+Y : Redo
- CTRL+F : Search
- F3 : Search next item
- Shift+E/D : Transpose selection +1 / -1 semitone
- Shift+R/F : Transpose selection +1 / -1 octave
- CTRL+C: Copy
- CTRL+V: Paste
- CTRL+X: Cut
- CTRL+A : Select all
- Arrows : Move cursor
- Shift+Arrows : Select things
- Tab/Shift+Tab : Move cursor to the next/previous channel
- PageUp/PageDown : Move the cursor to the first/last row
- +/- : Modify a value

### Instrument editor
- CTRL+Z : Undo
- CTRL+Y : Redo
- CTRL+C : Copy operator parameters (if an operator is selected)
- CTRL+V : Paste operator parameters (if an operator is selected)
- Right click : Copy/Paste parameters, or Link operators (algorithm edition)

## Effects

The effect column is divided in two parts, a letter representing the effect type, and a value. If you are familiar with trackers, you'll notice some similarities. A * suffix denotes global effects, they can be put on any channel

| Letter | Effect Name | Description | Range |
|--------|-------------|-------------|-------|
| A | Arpeggio | Create a fast, chiptune-like 3 notes arpeggio. The first note is the note on the pattern, the two others are defined by the value of this effect : the units represent the second note, in semitones from the first (0-9), the tens represent the third node (0-25 although the common usage range ) | 0-255 |
| B* | Pattern jump | Jump to another pattern. Use in combination with C (row jump) to create a loop point | 0-255 |
| C* | Row jump | Jump to another row (on the same pattern if no pattern jump is specified) | 0-255 |
| D | Note delay | Allows notes to be played with precise timing (1/8th row steps) | 0-7 |
| E | Portamento up | Increase the pitch of the current note. Value affects the slide speed | 0-255 |
| F | Portamento down | Decrease the pitch of the current note. Value affects the slide speed | 0-255 |
| G | Glissando | Create a smooth pitch transition between notes. Value controls the glissando speed | 0-255 |
| H | Vibrato | Add vibrato. It is persistent, so you can use other effects while keeping vibrato active. Put a H0 to stop it. | 0-255 |
| I | Pitch bend | A MIDI-like pitch bend. Unlike portamento up/down which continuously change pitch without any limit, pitch bend sets the pitch relative to the current note. 0=-2 semitones, 64 = -1 semitone, 127=nothing, 192=+1 semitone, 255=+2semitones. | 0-255 |
| J | Tremolo | Create a tremolo effect. Like vibrato it's a persistent effect, use J0 to stop it. | 0-255 |
| K | Instrument control | Modify a parameter of the current instrument. Use the instrument editor to select which parameter this effect controls | 0-255 |
| M | Channel volume | Change the channel volume. 0=muted, 99=full volume | 0-99 |
| N | Channel volume slide | Smooth channel volume slide. 0=fast decrease, 127=no change, 255=fast increase. | 0-255 |
| P | Panning slide | Create a smooth panning transition. 127 does nothing, lower values pans to the left faster as it comes close to 0, higher values pans to the right faster as it comes close to 255 | 0-255 |
| Q | Note retrigger | Fastly repeats the note on the same row, up to 8 times | 0-7 |
| R | Reverb send | Set the reverb amount. 0=no reverb, 99=maximum | 0-99 |
| T* | Set tempo | Change the tempo (BPM). | 0-255 |
| V* | Global volume | Set the global volume. 0=muted, 99=full volume | 0-99 |
| W* | Global volume slide | Smooth global volume slide. 0=fast decrease, 127=no change, 255=fast increase. | 0-255 |
| X | Channel panning | Set panning to the desired value. 0=left, 127=middle, 255=right. | 0-255 |
| Y* | Sync marker | Special markers for synchronizing events to the song. See the FM Lib API | 0-255 |

## Import/Export

### Import
- **MIDI** : good support of midi format 0/1 files. If the tempo is buggy, it's because of the quantization (default 1 beat = 8 rows). You need to find the right quantization value to do a proper MIDI import. 1 beat for 8 rows works for most songs. For ternary, try 1 for 12. You can also check the option "Preserve unquantized notes" in the configuration page. This will add a 'D' effect to the notes that doesn't fit perfectly in a row. D has 8 possible delay values, so it makes the import precision 8 times higher.

### Export
- **MIDI format 0** : export is here but quite primitive, don't use FM Composer as a serious MIDI editor.
- **WAVE** : renders the song as 16 bit stereo audio, with the sample rate of your choice. Please note that, due to the way FM works, some percussions or highly distorted sounds may sound a bit different depending on the sample rate. Best is to export at the same sample rate you created the song then resample with an audio editor.

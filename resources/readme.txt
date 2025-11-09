|	___  ____   _______ _____              _             
|	|  \/  | | | |  _  \_   _|            | |            
|	| .  . | | | | | | | | |_ __ __ _  ___| | _____ _ __ 
|	| |\/| | | | | | | | | | '__/ _` |/ __| |/ / _ \ '__|
|	| |  | | |_| | |/ /  | | | | (_| | (__|   <  __/ |   
|	\_|  |_/\___/|___/   \_/_|  \__,_|\___|_|\_\___|_|                                                    
|                                                     
|                                 v1.0 (????-??-??)
|
| Developed by Dashodanger -- https://github.com/dashodanger/mudtracker
| 
| Forked from FM Composer by Stephane "Phanoo" Damo -- 
|	https://web.archive.org/web/20180727172405/http://fmcomposer.org/en/ (Original Site)
|
|___________________________________________________________________________________________________________________


Thank you for downloading MUDTracker !

*** Changelog ***

v1.0 (????-??-??)
    - [Feature] Working CMake for GCC/Clang + Unix (with optional install target)
	- [Feature] Added support for --appdir parameter for separate resource directory
	- [Feature] Added support for DMX MUS (Doom format) import
	- [Feature] Will load entire General MIDI instrument set when importing MIDI/MUS
	- [Feature] Added "Load GM" button to Instrument Editor view to load General MIDI defaults
	- [Fix] Added missing <sstream> include causing Clang builds to fail
	- [Fix] Addressed GCC/Clang compilation warnings
	- [Fix] Update link changed to Github repository
	- [Maintenance] Updated portaudio/portmidi/simpleini/tinyfiledialogs to latest releases
	- [Maintenance] Deprecated SFML2 functions migrated to their replacements
	- [Removal] FLAC and MP3 dependencies/export options removed
	- [Removal] Instruments and songs are no longer LZ4 compressed

-- Below changelogs are for the original FM Composer --

v1.7 (2018-05-11)
	- [Fix] Crash after undoing a pasted selection that was out of the pattern bounds
	- [Fix] Under some circumstances, the song playback position was changed when navigating to the Pattern screen
	- [FM engine fix] Portamento effect (G) was wrong when Global Transposition was used

v1.6 (2018-04-14)
	- [Feature] Multi-track export !
	- [Fix] Glitchy 24-bit FLAC export / Wrong duration on VBR MP3 export
	- [Fix] Rendering is now anti-aliased for non-round DPI screens (eg. 120 DPI)
	- [Fix] The textbox had some edge cases that could lead to a crash

v1.5 (2018-04-04)
	- [Feature] Export bit depth choice for WAVE and FLAC exports
	- [Feature] Actions done through the right click menu in patterns are now undoable
	- [Feature] DirectSound devices are now shown in the available sound devices
	- [Fix] Can't add instruments to a new song if the default one wasn't found
	- [Fix] Crash when 'Remove unused' was clicked if the song had 1 instrument
	- [Fix] Messing with sound devices could crash under some circumstances

v1.4 (2018-03-03)
	- [Feature] FLAC export
	- [Feature] Editing step : entering a note automatically skips 0...16 rows for fast beat making (see 'Editing step' slider on the right side of the Pattern view)
	- [Feature] Support for high-DPI monitors
	- [Optimization] Reduced CPU usage
	- [Fix] Export bugs : 'stop to pattern xxx' was ignored for MP3 exports / song kept playing after exporting
	- [Fix] LAME MP3 is now used as dynamic library (needed by its license)

v1.3 (2018-02-18)
	- [Fix] Selection cursor was misplaced/missized after a Paste action 
	- [Fix] Player engine : notes without volume used the previous note volume, but ignored volumes changes inbetween
	- [Fix] Minor graphical glitch : channel header on the right side sometimes needs the user to scroll a bit more to show up

v1.2 (2018-02-04)
	- [Fix] 'Remove rows' function did nothing
	- [Fix] Small click noise occured on the next note after a channel was unmuted
	- [Fix] Recent songs menu width not immediately updated
	- [Fix] Blurry pattern top (channel params)

v1.1 (2018-01-31)
	- [Feature] Added some new demo songs, mostly from imported MIDIs
	- [Fix] Crash on Undo action under certain circumstances 
	- [Fix] Pattern list was selected on mouse release even if it wasn't focused
	- [Fix] small 1-frame graphical glitch when scrolling in pattern list

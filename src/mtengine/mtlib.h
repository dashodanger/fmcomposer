#ifdef __cplusplus
extern "C"{
#endif

#ifndef MTLIB_H
#define MTLIB_H

	/* Number of channels (polyphony) */
#define FM_ch 24
	/* Number of operators */
#define FM_op 6


	enum{ FM_NOTE, FM_INSTR, FM_VOL, FM_FXTYPE, FM_FXVALUE };
	enum { MT_ERR_FILEIO = -1, MT_ERR_FILECORRUPTED = -2, MT_ERR_FILEVERSION = -3 };
	enum fmInstrumentFlags{FM_INSTR_LFORESET=1, FM_INSTR_SMOOTH=2, FM_INSTR_TRANSPOSABLE=4};
	enum mtRenderTypes{MT_RENDER_8, MT_RENDER_16, MT_RENDER_24, MT_RENDER_32, MT_RENDER_FLOAT, MT_RENDER_PAD32=64};
	typedef struct fm_instrument_operator
	{
		unsigned char mult;
		unsigned char finetune;
		char detune;
		char vol;
		char connectOut;
		char connect;
		char connect2;

		char delay, i, a, h, d, s, r;

		char kbdVolScaling;

		char kbdAScaling;
		char kbdDScaling;
		char velSensitivity;
		char lfoFM;
		char lfoAM;
		char waveform;
		char fixedFreq;

		char kbdCenterNote;
		char kbdPitchScaling;

		char pitchInitialRatio;
		char pitchDecay;
		char pitchFinalRatio;
		char pitchRelease;

		char envLoop;
		char muted;

		unsigned char offset;

	}fm_instrument_operator;


	/* MDTI file format */
	typedef struct fm_instrument{
		char magic[4]; /* "MDTI" identifier */
		unsigned char dummy;
		unsigned char version;
		char name[24];
		char toMix[4];
		char feedback;
		char lfoDelay;
		char lfoSpeed;
		char lfoA;
		char lfoWaveform;
		char lfoOffset;
		char volume;
		char feedbackSource;
		char envReset;
		char tuning;
		char transpose;
		char phaseReset;
		char flags;
		char temperament[12];
		unsigned char kfx;
		fm_instrument_operator op[FM_op];
	}fm_instrument;



	typedef struct Cell{
		unsigned char note;
		unsigned char instr;
		unsigned char vol;
		unsigned char fx;
		unsigned char fxdata;
	}Cell;


	typedef struct ChannelState{
		float time;
		unsigned char tempo;
		unsigned char vol[FM_ch];
		unsigned char pan[FM_ch];
	}ChannelState;

	typedef struct fm_operator{
		// dynamic operator data
		float *connect, *connect2, *connectOut, *toMix;
		float out;
		float *waveform;
		unsigned int phase; // 10.10 bit phase accumulator (10 MSB used for sine lookup table)
		unsigned int pitch;
		float amp, ampDelta;



		float			env;
		unsigned		state;
		float				incr;
		float				pitchMod, pitchTime, pitchDestRatio;
		// static operator data
		float		kbdVolScaling;
		unsigned		envCount;
		float	vol;

		unsigned	delay, offset;
		float	i, a, h, d, s, r;
		float prevAmp, realAmp, ampTarget;
		float	lfoFM, lfoAM;
		float portaDestIncr;
		unsigned char id, mult, baseVol, kbdCenterNote;
		char baseA, baseD;
		char finetune, detune, fixedFreq, envLoop, pitchFinalRatio;
		float velSensitivity, volScaling;
	}fm_operator;


	typedef struct fm_channel{
		int newNote;

		// real time panning/volume
		float pan, vol;

		// used for smoothing panning/volume changes (avoid clicks)
		float destPan, destVol;

		// initial song channel pannings/volumes/reverb amounts
		unsigned char initial_pan, initial_vol, initial_reverb;


		float instrVol;
		float reverbSend;
		int muted;

		float *feedbackSource;
		float mixer;


		// DAHDSR envelope
		float currentEnvLevel;
		float	feedbackLevel;
		int op0;
		int realChannelNumber;
		float lfo;
		// channels
		unsigned int	lfoPhase, lfoMask, lfoIncr, lfoWaveform;
		unsigned algo;
		unsigned char note, baseArpeggioNote;
		float arpTimer;
		int arpIter;
		fm_instrument* instr;
		unsigned char instrNumber;
		float ramping;
		float rampingPicture;
		float lfoEnv, lfoA, lfoFMCurrentValue, lfoAMCurrentValue;
		unsigned lfoDelayCpt, lfoDelayCptMax, lfoOffset;
		unsigned char fxActive, fxData, noteVol, untransposedNote;
		char transpose;
		unsigned active;


		float fadeFrom, fadeFrom2, fadeIncr, fade, tuning;
		float lastRender, lastRender2, delta;
		float pitchBend;
		fm_instrument* cInstr;
		fm_operator op[FM_op];

	}fm_channel;


	typedef struct mtsynth{
		char songName[64], author[64], comments[256];
		float globalVolume;
		float playbackVolume;
		unsigned char _globalVolume;
		fm_instrument *instrument; // here are stored your instruments
		unsigned char instrumentCount;
		char transpose, looping;
		unsigned char loopCount;
		int channelStatesDone;

		// reverb
		unsigned reverbPhaseL, reverbPhaseL2, reverbPhaseR, reverbPhaseR2;
		float* revBuf;
		unsigned revBufSize;

		unsigned allpassPhaseL, allpassPhaseR, allpassPhaseL2, allpassPhaseR2;

		unsigned allpassMod, allpassMod2, reverbMod1, reverbMod2, reverbMod3, reverbMod4;
		unsigned revOffset2, revOffset3, revOffset4, revOffset5, revOffset6, revOffset7;

		float reverbRoomSize, initialReverbRoomSize;
		float reverbLength, initialReverbLength;

		unsigned char tempo, initial_tempo;
		unsigned row, order, playing, saturated;

		float noConnect;
		ChannelState **channelStates;
		Cell(**pattern)[FM_ch];
		unsigned patternCount;
		unsigned *patternSize;
		unsigned frameTimer;
		float frameTimerFx;

		float outL, outR;
		unsigned char diviseur;
		float sampleRateRatio;


		unsigned sampleRate;
		float noteIncr[128];


		fm_channel ch[FM_ch];

		float transitionSpeed;
		int tempRow, tempOrder;
		unsigned readSeek, totalFileSize;
	}mtsynth;


	/** Creates the synth
		@param samplerate : sample rate in Hz
		@return pointer to the allocated fm synth
		*/
	mtsynth* mt_create(int samplerate);

	void mt_destroy(mtsynth* mt);

	/** Load a song
		@param filename
		@return 1 if ok, 0 if failed
		*/
	int mt_loadSong(mtsynth* mt, const char* filename);



	/** Load a song
		@param filename
		@return 1 if ok, 0 if failed
		*/
	int mt_loadSongFromMemory(mtsynth* mt, char* data, unsigned len);

	/** Get total song length
		@return length in seconds
		*/
	float mt_getSongLength(mtsynth* mt);

	/** Play the song */
	void mt_play(mtsynth* mt);

	/** Stop the song
		@param mode : 0 = force note off, 1 = hard cut
		*/
	void mt_stop(mtsynth* mt, int mode);


	void mt_setPlaybackVolume(mtsynth *mt, int volume);

	/** Set the global volume
		@param volume : volume, 0-99
		*/
	void mt_setVolume(mtsynth *mt, int volume);

	/** Render the sound
		@param buffer : audio buffer, left and right channels are interleaved
		@param length : number of samples to render
		*/
	void mt_render(mtsynth* mt, void* buffer, unsigned length, unsigned type);

	/** Play a note
		@param instrument : instrument number, 0-255
		@param note : midi note number, 0-127 (C0 - G10)
		@param channel : channel number, 0-23
		@param volume : volume, 0-99
		*/
	void mt_playNote(mtsynth* mt, unsigned instrument, unsigned note, unsigned channel, unsigned volume);

	/** Stops a note on a channel
		@param channel : channel number, 0-23
		*/
	void mt_stopNote(mtsynth* mt, unsigned channel);





	/** Set the playing position
		@param pattern : pattern number
		@param row : row number
		@param mode : 0 = keep playing notes, 1 = force note off, 2 = hard cut
		*/
	void mt_setPosition(mtsynth* mt, int pattern, int row, int mode);

	/** Set the playing position in seconds
		@param time : position in seconds
		@param mode : 0 = keep playing notes, 1 = force note off, 2 = hard cut
		*/
	void mt_setTime(mtsynth* mt, int time, int mode);

	/** Set the tempo
		@param tempo : tempo in BPM, 0-255
		*/
	void mt_setTempo(mtsynth* mt, int tempo);


	/** Get the current playing position
		@param *order : current pattern number
		@param *row : current row number
		*/
	void mt_getPosition(mtsynth* mt, int *pattern, int *row);

	/** Get the playing time in seconds
		@return current playing time in seconds
		*/
	float mt_getTime(mtsynth* mt);

	/** Set the sample rate
		@param samplerate : sample rate in Hz
		@return 1 if ok, 0 if failed (keeps the previous sample rate)
		*/
	int mt_setSampleRate(mtsynth* mt, int samplerate);

	/** Set channel volume
		@param channel : channel number, 0-23
		@param volume : volume, 0-99
		*/
	void mt_setChannelVolume(mtsynth *mt, int channel, int volume);

	/** Set channel panning
		@param channel : channel number, 0-23
		@param panning : panning, 0-255
		*/
	void mt_setChannelPanning(mtsynth *mt, int channel, int panning);

	/** Set channel reverb
		@param channel : channel number, 0-23
		@param panning : reverb amount, 0-99
		*/
	void mt_setChannelReverb(mtsynth *mt, int channel, int reverb);

	/** Write data to a pattern
		@param pattern : pattern number
		@param row : row number
		@param channel : channel number
		@param type : the column to write into (FM_NOTE, FM_INSTR, FM_VOL, FM_FXTYPE, FM_FXVALUE)
		@param value : the value to write. 255 is considered empty
		@return 1 if success, 0 if failed (pattern/row/channel/type out of bounds)
		*/
	int mt_write(mtsynth *mt, unsigned pattern, unsigned row, unsigned channel, Cell data);

	/** Create a new pattern at the desired position
		@param rows : number of rows
		@param position : the position where to insert the pattern
		@return 1 if success, 0 if failed
		*/
	int mt_insertPattern(mtsynth* mt, unsigned rows, unsigned position);

	/** Remove a pattern
		@param pattern : the pattern number
		@return 1 if success, 0 if failed
		*/
	int mt_removePattern(mtsynth* mt, unsigned pattern);

	/** Resize a pattern. Contents are not stretched/scaled.
		@param pattern : the pattern number
		@param size : the new size
		@return 1 if success, 0 if failed
		*/
	int mt_resizePattern(mtsynth* mt, unsigned pattern, unsigned size, unsigned scaleContent);







	/* ########################################################### */

	/* Stops a note by its number */
	void mt_stopNoteID(mtsynth* mt, unsigned note);



	void mt_clearSong(mtsynth* mt);

	int mt_getPatternSize(mtsynth *mt, int pattern);

	int mt_resizeInstrumentList(mtsynth* mt, unsigned size);
	void mt_portamento(mtsynth* mt, unsigned channel, float value);


	void mt_patternClear(mtsynth* mt);
	int mt_resizePatterns(mtsynth* mt, unsigned count);
	int mt_loadInstrument(mtsynth* mt, const char *filename, unsigned slot);
	int mt_loadInstrumentFromMemory(mtsynth* mt, char *data, unsigned slot);
	int mt_loadInstrumentBank(mtsynth* mt, const char *filename);
	int mt_loadInstrumentBankFromMemory(mtsynth* mt, char *data);
	int mt_saveInstrument(mtsynth* mt, const char* filename, unsigned slot);
	int mt_saveInstrumentBank(mtsynth* mt, const char* filename);
	void mt_removeInstrument(mtsynth* mt, unsigned slot, int removeOccurences);
	void mt_movePattern(mtsynth* mt, int from, int to);
	/* Saves the song to file */
	int mt_saveSong(mtsynth* mt, const char* filename);


	void mt_buildStateTable(mtsynth* mt, unsigned orderStart, unsigned orderEnd, unsigned channelStart, unsigned channelEnd);
	int mt_initReverb(mtsynth *mt, float roomSize);

	/* Forces all sound to stop. Cut notes and reverb. */
	void mt_stopSound(mtsynth* mt);
	void mt_moveChannels(mtsynth* mt, int from, int to);
	int mt_clearPattern(mtsynth* mt, unsigned pattern, unsigned rowStart, unsigned count);
	int mt_insertRows(mtsynth *mt, unsigned pattern, unsigned row, unsigned count);
	int mt_removeRows(mtsynth *mt, unsigned pattern, unsigned row, unsigned count);


	float mt_volumeToExp(int volume);

	int mt_isInstrumentUsed(mtsynth *mt, unsigned id);
	void mt_createDefaultInstrument(mtsynth* mt, unsigned slot);
#endif

#ifdef __cplusplus
}
#endif
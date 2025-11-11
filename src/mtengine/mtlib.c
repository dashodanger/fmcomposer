#include "mtlib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <xmmintrin.h>
#include <float.h>
#include <math.h>

/* Current version of instrument/song formats */
#define MUDTRACKER_VERSION 1

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define _99TO1 1/99
#define _24TO1 1/24
#define _2400TO1 1/2400
#define SEMITONE_RATIO 0.059463 * 0.01 /* 0.059463 = ratio between two semitones https://en.wikipedia.org/wiki/Twelfth_root_of_two */

#ifdef _MSC_VER
#undef min
#undef max
#endif

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define clamp(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

/* Sine wave lookup table size */

#define LUTsize 2048
#define LUTratio (LUTsize / 1024)

/* Reverb delays, in number of samples at 48000Hz. Automatically scaled for other samples rates. */

#define REVERB_DELAY_L1 1.6*4096 // 85ms
#define REVERB_DELAY_L2 1.5*2485 // 72 
#define REVERB_DELAY_R1 1.6*3801 // 79
#define REVERB_DELAY_R2 1.5*2333 // 69
#define REVERB_ALLPASS2 1.5*1170 // 5.5
#define REVERB_ALLPASS1 1.5*2508 // 7.7ms


static float wavetable[8][LUTsize];


/* Exponential tables for envelopes and volumes scales */

static float expEnv[100], expVol[100], expVolOp[100];

// waveforms
static const unsigned int lfoMasks[28] = {
	0xffc00 * LUTratio, // sine
	0xffc00 * LUTratio, // sine
	0xffc00 * LUTratio, // sine
	0xffc00 * LUTratio, // sine
	0xffc00 * LUTratio, // sine
	0xffc00 * LUTratio, // sine
	0xffc00 * LUTratio, // sine
	0xffc00 * LUTratio, // sine
	0xf0000 * LUTratio,  // less reso
	0xefc00 * LUTratio,  // less reso
	0xdfc00 * LUTratio,  // less
	0xbfc00 * LUTratio,  // squarelike
	0x88000 * LUTratio, // high freq
	0x40000 * LUTratio, // square
	0x60000 * LUTratio, // pulse/near square
	0x7fc00 * LUTratio, // abs sine
	0x78000 * LUTratio, // abs sine less reso
	0x70000 * LUTratio, // less reso
	0x3fc00 * LUTratio, // saw
	0xa0000 * LUTratio, // ? 10
	0xfffc00 * LUTratio, // sine
	0x2ffc00 * LUTratio, // sine
};

static const unsigned int lfoWaveforms[28] = {
	0,
	1,
	2,
	3,
	4,
	5,
	6,
	7,
	0,  // less reso
	0,  // less reso
	0,  // less
	0,  // squarelike
	0, // high freq
	0, // square
	0, // pulse/near square
	0, // abs sine
	0, // abs sine less reso
	0, // less reso
	0, // saw
	0, // ? 10
	0, // sine
	0, // sine
};


/* Band-limited triangle, square and sawtooth generators */

float trg(float x, float theta) { return 1 - 2 * acos((1 - theta) * sin(2 * M_PI*x)) / M_PI; }

float sqr(float x, float theta) { return 2 * atan(sin(2 * M_PI* x) / theta) / M_PI; }

float swt(float x, float theta) { return (1 + trg((2 * x - 1) / 4, theta) * sqr(x / 2, theta)) / 2; }

void mt_setDefaults(mtsynth* mt)
{

	for (unsigned ch = 0; ch < FM_ch; ++ch)
	{
		mt->ch[ch].note = 255;
		mt->ch[ch].instrNumber = 255;
		mt->ch[ch].vol = expVol[99];
		mt->ch[ch].initial_vol = 99;
		mt->ch[ch].reverbSend = 0;
		mt->ch[ch].destPan = mt->ch[ch].pan = mt->ch[ch].initial_pan = 127;
		mt->ch[ch].noteVol = 99;
	}

	mt_setVolume(mt, 60);
	mt->initial_tempo = 120;
	mt->diviseur = 4;
	mt->initialReverbLength = mt->reverbLength = 0.875;
	mt->initialReverbRoomSize = 0.55;
	mt->looping = -1;
	mt->channelStatesDone = 0;
	mt->playbackVolume = 1;
}



static unsigned int g_seed = 0;
int fast_rand(void)
{
	g_seed = (214013 * g_seed + 2531011);
	return (g_seed >> 16) & 0x7FFF;
}

void mt_destroy(mtsynth* mt)
{
	free(mt->revBuf);
	free(mt->instrument);
	for (unsigned i = 0; i < mt->patternCount; i++)
	{
		free(mt->pattern[i]);
		free(mt->channelStates[i]);
	}
	free(mt->patternSize);
	free(mt->pattern);
	free(mt->channelStates);
	free(mt);
}

mtsynth* mt_create(int _sampleRate)
{
	mtsynth *mt = calloc(1, sizeof(mtsynth));

	if (mt)
	{

		/* avoid slow float denormals + round floats down (needed if program compiled with QIfirst option (FISTP) fast int/float conversion) */
		#ifdef _WIN32
		_control87(_RC_DOWN, _MCW_RC);
		#endif
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		/* Build waveform tables */

		for (unsigned i = 0; i < LUTsize; i++)
			wavetable[0][i] = sin(i * 2 * M_PI / LUTsize);					// 0 sine

		for (unsigned i = 0; i < LUTsize; i++)
			wavetable[1][i] = (swt((float)(i + LUTsize / 2) / LUTsize, 0.2) - 0.5)*2.5*(1 / 0.464670);	// 3 soft saw

		for (unsigned i = 0; i < LUTsize; i++)
			wavetable[2][i] = (swt((float)(i + LUTsize / 2) / LUTsize, 0.05) - 0.5) * 2 * (1 / 0.649969);	// 4 saw

		for (unsigned i = 0; i < LUTsize; i++)
			wavetable[3][i] = trg((float)i / LUTsize, 0.01)*(1 / 0.909893);			// 1 triangle

		for (unsigned i = 0; i < LUTsize; i++)
			wavetable[4][i] = sqr((float)i / LUTsize, 0.1)*0.7*(1 / 0.655584);			// 2 square

		for (unsigned i = 0; i < LUTsize / 2; i++)
			wavetable[5][i] = sin(i * 2 * M_PI / (LUTsize / 2));				// 5 double sine

		for (unsigned i = 0; i < LUTsize / 2; i++)
			wavetable[6][i] = sin(i * 2 * M_PI / LUTsize);					// 6 half period sine

		for (unsigned i = 0; i < LUTsize; i++)
			wavetable[7][i] = fast_rand() / 16383.5 - 0.5;

		/*for (int i = 0; i < 7; i++){
			float max=0;
			for (int j = 0; j< LUTsize; j++){
			if (wavetable[i][j]>max){
			max = wavetable[i][j];

			}
			}
			printf("max %d is %f\n", i, max);
			}*/



		/* Build exponential tables for volume/envelopes */

		float ini = 0.00001;
		for (unsigned i = 1; i < 99; i++)
		{
			expVol[i] = pow(10, (log(100.0 / (i + 1)) * (-10)) / 20.0);
			expEnv[i] = ini;
			ini *= 1.1;
			expVolOp[i] = expVol[i] * (i*0.01);
		}

		expEnv[96] = 0.1;
		expEnv[97] = 0.2;
		expEnv[98] = 0.5;
		expEnv[99] = expVol[99] = expVolOp[99] = 1;

		mt_setDefaults(mt);

		if (!mt_setSampleRate(mt, _sampleRate))
		{
			free(mt);
			return 0;
		}
	}

	return mt;
}

int mt_initReverb(mtsynth *mt, float roomSize)
{
	/* Initialize reverb parameters and buffers */

	mt->reverbPhaseL = mt->reverbPhaseL2 = mt->reverbPhaseR = mt->reverbPhaseR2 = mt->allpassPhaseL = mt->allpassPhaseR = mt->allpassPhaseL2 = mt->allpassPhaseR2 = 0;

	unsigned mod1 = roomSize*REVERB_DELAY_L1 / mt->sampleRateRatio; // 85ms
	unsigned mod2 = roomSize* REVERB_DELAY_L2 / mt->sampleRateRatio; // 72 
	unsigned mod3 = roomSize*REVERB_DELAY_R1 / mt->sampleRateRatio; // 79
	unsigned mod4 = roomSize*REVERB_DELAY_R2 / mt->sampleRateRatio; // 69
	unsigned mod5 = (roomSize*REVERB_ALLPASS1) / mt->sampleRateRatio; // 5.5
	unsigned mod6 = (roomSize*REVERB_ALLPASS2) / mt->sampleRateRatio; // 7.7ms

	unsigned revBufSize = mod1 + mod2 + mod3 + mod4 + 2 * (mod5 + mod6);

	float* newR = realloc(mt->revBuf, sizeof(float)*revBufSize);

	if (!newR)
	{
		return 0;
	}

	mt->reverbRoomSize = roomSize;
	mt->revBufSize = revBufSize;
	mt->revBuf = newR;

	memset(mt->revBuf, 0, sizeof(float)*revBufSize);



	mt->reverbMod1 = mod1;
	mt->reverbMod2 = mod2;
	mt->reverbMod3 = mod3;
	mt->reverbMod4 = mod4;
	mt->allpassMod = mod5;
	mt->allpassMod2 = mod6;

	mt->revOffset2 = mt->reverbMod1 + mt->reverbMod2;
	mt->revOffset3 = mt->revOffset2 + mt->reverbMod3;
	mt->revOffset4 = mt->revOffset3 + mt->reverbMod4;
	mt->revOffset5 = mt->revOffset4 + mt->allpassMod;
	mt->revOffset6 = mt->revOffset5 + mt->allpassMod;
	mt->revOffset7 = mt->revOffset6 + mt->allpassMod2;
	return 1;
}

int mt_setSampleRate(mtsynth* mt, int sampleRate)
{

	mt->sampleRate = sampleRate;
	mt->sampleRateRatio = 48000.0 / sampleRate;
	mt->transitionSpeed = 20 * (1 / mt->sampleRateRatio);

	/* initialize the MIDI note frequencies table (converted into phase accumulator increments) */

	for (unsigned x = 0; x < 128; ++x)
		mt->noteIncr[x] = pow(2, (x - 9.0) / 12.0) / sampleRate * 32840 * 440 * LUTratio;

	/* Reset instruments so their values can be regenerated */
	for (unsigned ch = 0; ch < FM_ch; ++ch)
		mt->ch[ch].cInstr = 0;

	return mt_initReverb(mt, mt->initialReverbRoomSize);
}

/* Calculates the volume of each operator */
void mt_calcOpVol(fm_operator *o, int note, int volume)
{

	float noteScaling = 1 + (note - o->kbdCenterNote)*o->volScaling;
	float opVol = (expVol[volume] * o->velSensitivity + (1 - o->velSensitivity))*expVolOp[o->baseVol];
	o->vol = clamp(opVol * noteScaling, 0, 1) * 5000 * LUTratio;

}

/* Calculates the pitch of each operator */
void mt_calcPitch(mtsynth* mt, int ch, int note)
{

	note = clamp(note + mt->ch[ch].transpose + mt->transpose * ((mt->ch[ch].instr->flags & FM_INSTR_TRANSPOSABLE) >> 2), 0, 127);

	mt->ch[ch].note = note;

	float frequency = mt->noteIncr[note] + mt->noteIncr[note] * SEMITONE_RATIO* mt->ch[ch].instr->temperament[note % 12];

	for (unsigned op = 0; op < FM_op; ++op)
	{
		fm_operator* o = &mt->ch[ch].op[op];

		if (o->fixedFreq == 0)
		{
			o->incr = frequency *(o->mult + (float)o->finetune*_24TO1 + (float)o->detune*_2400TO1);
		}
		/* Fixed frequency */
		else
			o->incr = (o->mult * o->mult + (float)o->mult *(float)o->finetune*_24TO1) * LUTratio*mt->sampleRateRatio;

		o->incr += o->incr*mt->ch[ch].tuning;
	}
}


void mt_initChannels(mtsynth* mt)
{
	if (mt->order >= mt->patternCount || mt->row >= mt->patternSize[mt->order])
		return;

	mt->tempo = mt->channelStates[mt->order][mt->row].tempo;
	mt->globalVolume = expVol[mt->_globalVolume] * 4096 / LUTsize;
	mt->reverbLength = mt->initialReverbLength;
	if (mt->initialReverbRoomSize != mt->reverbRoomSize)
	{
		mt_initReverb(mt, mt->initialReverbRoomSize);
	}

	for (unsigned ch = 0; ch < FM_ch; ++ch)
	{
		mt->ch[ch].cInstr = 0;
		mt->ch[ch].pan = mt->ch[ch].destPan = mt->channelStates[mt->order][mt->row].pan[ch];
		mt->ch[ch].vol = expVol[mt->channelStates[mt->order][mt->row].vol[ch]];
		mt->ch[ch].reverbSend = expVol[mt->ch[ch].initial_reverb];
		mt->ch[ch].pitchBend = 1;
		mt->ch[ch].fadeFrom=0;
		mt->ch[ch].fadeFrom2=0;
		mt->ch[ch].currentEnvLevel=0;
	}
}




void _mt_render(mtsynth* mt, float* buffer, unsigned length)
{

	unsigned b = 0;
	while (b < length)
	{
		// player
		if (mt->playing)
		{
			/* Song frame tick */
			if (mt->frameTimer == 0)
			{

				for (unsigned ch = 0; ch < FM_ch; ++ch)
				{
					Cell* row = &mt->pattern[mt->order][mt->row][ch];
					mt->ch[ch].fxActive = 0;

					/* Note stop ? */
					if (row->note == 128 && row->fx != 'D')
					{
						mt_stopNote(mt, ch);
					}
					/* Note play ? */
					else if (row->note != 255)
					{
						// portamento or note delay : don't play the note now !
						if (row->fx != 'D' && (row->fx != 'G' || row->fx == 'G' && mt->ch[ch].note == 255))
						{
							mt_playNote(mt, row->instr, row->note, ch, row->vol);
							mt->ch[ch].baseArpeggioNote = row->note;
						}
					}
					// only volume change
					else if (row->vol != 255 && mt->ch[ch].instr)
					{
						mt->ch[ch].noteVol = row->vol;
						for (unsigned op = 0; op < FM_op; ++op)
						{
							mt_calcOpVol(&mt->ch[ch].op[op], mt->ch[ch].note, row->vol);
						}
					}

					/* Handle effects (after note actions) */

					mt->ch[ch].fxData = row->fxdata;
					switch (row->fx)
					{

						case 'B': // jump pattern
							mt->tempOrder = mt->ch[ch].fxData;
							break;
						case 'C': // jump row
							mt->tempRow = mt->ch[ch].fxData;
							break;
						case 'G': // portamento
							// update portamento dest frequency if note set
							if (row->note != 255)
							{
								for (unsigned op = 0; op < FM_op; ++op)
								{
									float pitchScaling = 1 + ((int)row->note - mt->ch[ch].instr->op[op].kbdCenterNote)*mt->ch[ch].instr->op[op].kbdPitchScaling*0.001;

									if (mt->ch[ch].instr->op[op].fixedFreq == 0)
									{
										mt->ch[ch].op[op].portaDestIncr = mt->noteIncr[clamp(row->note+mt->transpose,0,127)] * pitchScaling*(mt->ch[ch].op[op].mult + (double)mt->ch[ch].op[op].finetune*0.041666666667 + (double)mt->ch[ch].op[op].detune*0.00041666666667);

									}
									else // fixed frequency
										mt->ch[ch].op[op].portaDestIncr = (mt->ch[ch].op[op].mult * (mt->ch[ch].op[op].mult) + (double)mt->ch[ch].op[op].mult*(double)mt->ch[ch].op[op].finetune*0.041666666667) * LUTratio;
								}
							}
							// repeated effects
						case 'A': // arpeggio
						case 'D': // delay
						case 'E': // portamento up
						case 'F': // portamento down
						case 'W': // global volume slide
						case 'P': // panning slide

							mt->ch[ch].arpTimer = 0;
							mt->ch[ch].arpIter = 0;
							mt->ch[ch].fxActive = row->fx;
							break;
						case 'Q': // retrigger note
							if (mt->ch[ch].fxData > 0)
							{
								mt->ch[ch].arpTimer = 24 / mt->ch[ch].fxData;
								mt->ch[ch].arpIter = 0;
								mt->ch[ch].fxActive = row->fx;
							}
							break;
						case 'H': // vibrato
							mt->ch[ch].lfoEnv = 1;
							mt->ch[ch].lfoIncr = (mt->ch[ch].fxData / 16 * 128)*LUTratio;
							for (unsigned op = 0; op < FM_op; ++op)
							{
								mt->ch[ch].op[op].lfoFM = (mt->ch[ch].fxData % 16)*0.003;
							}
							break;
						case 'I': // pitch bend
							mt->ch[ch].fxActive = 'I';
							mt->ch[ch].pitchBend = 1 - (float)(128 - mt->ch[ch].fxData) * 0.00092852373168154813872606848242328;
							break;
						case 'J': // tremolo
							mt->ch[ch].lfoEnv = 1;
							mt->ch[ch].lfoIncr = (mt->ch[ch].fxData / 16 * 128)*LUTratio;
							for (unsigned op = 0; op < FM_op; ++op)
							{
								mt->ch[ch].op[op].lfoAM = (mt->ch[ch].fxData % 16)*(1.0 / 16);
							}
							break;
						case 'K':
							if (!mt->ch[ch].instr)
								break;
							/* Global instrument edit*/
							if (mt->ch[ch].instr->kfx / 32 == 0)
							{
								switch (mt->ch[ch].instr->kfx)
								{
									case 0:
										mt->ch[ch].instrVol = expVol[min(99, mt->ch[ch].fxData)];
										break;
									case 1:
										mt->ch[ch].transpose = clamp((char)mt->ch[ch].fxData, -12, 12);
										mt_calcPitch(mt, ch, mt->ch[ch].untransposedNote);
										break;
									case 2:
										mt->ch[ch].tuning = 0.0006*clamp((char)mt->ch[ch].fxData, -100, 100);
										mt_calcPitch(mt, ch, mt->ch[ch].untransposedNote);
										break;
									case 3:
										mt->ch[ch].lfoIncr = 1 + expVol[clamp(mt->ch[ch].fxData, 0, 99)] * expVol[clamp(mt->ch[ch].fxData, 0, 99)] * 5000 * mt->sampleRateRatio*LUTratio;
										break;
									case 4:
										mt->ch[ch].lfoDelayCptMax = expVol[clamp(mt->ch[ch].fxData, 0, 99)] * expVol[clamp(mt->ch[ch].fxData, 0, 99)] * 200000 * mt->sampleRateRatio;
										break;
									case 5:
										mt->ch[ch].lfoA = expEnv[clamp(mt->ch[ch].fxData, 0, 99)] * mt->sampleRateRatio;
										break;
									case 6:
										mt->ch[ch].lfoMask = lfoMasks[clamp(mt->ch[ch].fxData, 0, 19)];
										mt->ch[ch].lfoWaveform = lfoWaveforms[clamp(mt->ch[ch].fxData, 0, 19)];
										break;
									case 7:
										mt->ch[ch].lfoOffset = clamp(mt->ch[ch].fxData, 0, 31) * LUTsize / 32;
										break;

								}
							}
							/* Operator edit */
							else
							{
								fm_operator *o = &mt->ch[ch].op[mt->ch[ch].instr->kfx / 32 - 1];

								switch (mt->ch[ch].instr->kfx % 32)
								{
									case 0:
										o->baseVol = min(99, mt->ch[ch].fxData);
										mt_calcOpVol(o, mt->ch[ch].note, mt->ch[ch].noteVol);
										break;
									case 1:
										o->baseVol = o->vol * !mt->ch[ch].instr->op[mt->ch[ch].instr->kfx / 32 - 1].muted;
										mt_calcOpVol(o, mt->ch[ch].note, mt->ch[ch].noteVol);
										break;
									case 2:
										o->waveform = wavetable[clamp(mt->ch[ch].fxData, 0, 7)];
										break;
									case 3:{
											   o->mult = clamp(mt->ch[ch].fxData, 0, 40);
											   float frequency = mt->noteIncr[mt->ch[ch].note] + mt->noteIncr[mt->ch[ch].note] * SEMITONE_RATIO * mt->ch[ch].instr->temperament[mt->ch[ch].note % 12];
											   o->incr = frequency *(o->mult + (float)o->finetune*_24TO1 + (float)o->detune*_2400TO1) * (1 + mt->ch[ch].tuning);
											   break;
									}
									case 4:
										o->mult = clamp(mt->ch[ch].fxData, 0, 255);
										o->incr = (o->mult * o->mult + (float)o->mult *(float)o->finetune*_2400TO1) * LUTratio*mt->sampleRateRatio * (1 + mt->ch[ch].tuning);
										break;
									case 5:{
											   o->finetune = clamp(mt->ch[ch].fxData, 0, 24);
											   float frequency = mt->noteIncr[mt->ch[ch].note] + mt->noteIncr[mt->ch[ch].note] * SEMITONE_RATIO * mt->ch[ch].instr->temperament[mt->ch[ch].note % 12];
											   o->incr = frequency *(o->mult + (float)o->finetune*_24TO1 + (float)o->detune*_2400TO1) * (1 + mt->ch[ch].tuning);
											   break;
									}
									case 6:{
											   o->detune = clamp((char)mt->ch[ch].fxData, -100, 100);
											   float frequency = mt->noteIncr[mt->ch[ch].note] + mt->noteIncr[mt->ch[ch].note] * SEMITONE_RATIO * mt->ch[ch].instr->temperament[mt->ch[ch].note % 12];
											   o->incr = frequency *(o->mult + (float)o->finetune*_24TO1 + (float)o->detune*_2400TO1) * (1 + mt->ch[ch].tuning);
											   break;
									}
									case 7:
										o->delay = expEnv[mt->ch[ch].fxData] * 3000000 / mt->sampleRateRatio;
										break;
									case 8:
										o->i = expVol[mt->ch[ch].fxData];
										break;
									case 9:
										o->baseA = clamp(mt->ch[ch].fxData, 0, 99);
										break;
									case 10:
										o->h = expEnv[clamp(mt->ch[ch].fxData, 0, 80)] * 700000 / mt->sampleRateRatio;
										break;
									case 11:
										o->baseD = clamp(mt->ch[ch].fxData, 0, 99);
										break;
									case 12:
										o->s = expVol[clamp(mt->ch[ch].fxData, 0, 99)];
										break;
									case 13:{
												char value = clamp((char)mt->ch[ch].fxData, -99, 99);
												o->r = (value >= 0) ? exp(-(expEnv[value])*mt->sampleRateRatio) : 2 - exp(-(expEnv[abs(value)])*mt->sampleRateRatio);
												break;
									}
									case 14:
										o->envLoop = clamp(mt->ch[ch].fxData, 0, 1);
										break;
									case 15:
										o->lfoFM = expVol[clamp(mt->ch[ch].fxData, 0, 99)] * expVol[clamp(mt->ch[ch].fxData, 0, 99)];
										break;

									case 16:
										o->lfoAM = expVol[clamp(mt->ch[ch].fxData, 0, 99)];
										break;




								}
							}


							break;
						case 'M': // channel volume
							mt->ch[ch].vol = expVol[mt->ch[ch].fxData];
							break;
						case 'R': // reverb send
							mt->ch[ch].reverbSend = expVol[mt->ch[ch].fxData];
							break;
						case 'S': // global reverb params
							if (mt->ch[ch].fxData <= 40)
							{
								mt->reverbLength = 0.5 + mt->ch[ch].fxData*0.0125;
							}
							else
							{
								mt_initReverb(mt, clamp(mt->ch[ch].fxData - 40, 1, 40)*0.025);
							}
							break;
						case 'T': // tempo
							mt->tempo = max(1, mt->ch[ch].fxData);
							break;


						case 'X': // panning
							mt->ch[ch].destPan = mt->ch[ch].fxData;
							break;
					}

				}
			}
			mt->frameTimer += 8;
			if (mt->frameTimer >= (60.0 / mt->diviseur) * mt->sampleRate / mt->tempo)
			{

				mt->frameTimer = 0;

				if (++mt->row >= mt->patternSize[mt->order])
				{ // jump to next pattern
					mt->row = 0;
					mt->order++;
				}

				if (mt->tempOrder != -1 || mt->tempRow != -1)
				{
					mt->loopCount++;

					if (mt->tempOrder != -1)
						mt->order = min(mt->tempOrder, mt->patternCount - 1);

					if (mt->tempRow != -1)
						mt->row = min(mt->tempRow, mt->patternSize[mt->order] - 1);

					mt->tempOrder = mt->tempRow = -1;
				}

				if (mt->order >= mt->patternCount)
				{
					mt->loopCount++;
					mt->order = 0;
				}

				if (mt->looping != -1 && mt->loopCount > mt->looping)
				{
					mt->playing = 0;
				}

			}

			if (mt->frameTimerFx >= 0.005*(60.0 / mt->diviseur) * mt->sampleRate / mt->tempo)
			{

				for (unsigned ch = 0; ch < FM_ch; ++ch)
				{
					switch (mt->ch[ch].fxActive)
					{
						case 'A': // arpeggio
						{
									  mt->ch[ch].arpTimer++;
									  if (mt->ch[ch].arpTimer >= 8)
									  {
										  mt->ch[ch].arpTimer -= 8;
										  mt->ch[ch].arpIter = (mt->ch[ch].arpIter + 1) % 3;
										  mt_playNote(mt, 255, mt->ch[ch].arpIter == 0 ? mt->ch[ch].baseArpeggioNote : mt->ch[ch].arpIter == 1 ? (mt->ch[ch].baseArpeggioNote + mt->ch[ch].fxData % 16) : (mt->ch[ch].baseArpeggioNote + mt->ch[ch].fxData / 16), ch, 255);
									  }

						}
							break;
						case 'Q': // retrigger note
						{
									  mt->ch[ch].arpTimer++;
									  if (mt->ch[ch].arpTimer >= 24 / mt->ch[ch].fxData && mt->ch[ch].arpIter < mt->ch[ch].fxData)
									  {
										  mt->ch[ch].arpTimer -= 24 / mt->ch[ch].fxData;
										  mt_playNote(mt, mt->ch[ch].instrNumber, mt->ch[ch].untransposedNote, ch, 255);
										  mt->ch[ch].arpIter++;
									  }
						}
							break;
						case 'D': // delay
						{
									  int delay = mt->frameTimer / ((60.0 / mt->diviseur) * mt->sampleRate / mt->tempo / 8);
									  if (delay >= mt->ch[ch].fxData)
									  {

										  if (mt->pattern[mt->order][mt->row][ch].note < 127)
											  mt_playNote(mt, mt->pattern[mt->order][mt->row][ch].instr, mt->pattern[mt->order][mt->row][ch].note, ch, mt->pattern[mt->order][mt->row][ch].vol);
										  else if (mt->pattern[mt->order][mt->row][ch].note == 128)
											  mt_stopNote(mt, ch);

										  mt->ch[ch].fxActive = 0;
									  }
						}
							break;

						case 'E': // portamento up
							for (unsigned op = 0; op < FM_op; ++op)
							{
								mt->ch[ch].op[op].incr += mt->ch[ch].fxData*mt->ch[ch].op[op].incr*0.0001;
							}
							break;
						case 'F': // portamento down
							for (unsigned op = 0; op < FM_op; ++op)
							{
								mt->ch[ch].op[op].incr += -mt->ch[ch].fxData*mt->ch[ch].op[op].incr*0.0001;
							}
							break;
						case 'G': // portamento
							for (unsigned op = 0; op < FM_op; ++op)
							{
								mt->ch[ch].op[op].incr += (mt->ch[ch].op[op].portaDestIncr - mt->ch[ch].op[op].incr)*mt->ch[ch].fxData*0.001;
							}
							break;
						case 'I':{ // pitch bend

									 //int pos = mt->order*mt->patternSize[mt->order]+mt->row+1;

									 /*if (mt->pattern[pos / mt->patternSize[mt->order]][pos%mt->patternSize[mt->order]].fx == 'I'){
										 float nextPB = 1-(float)(64-mt->pattern[pos / mt->patternSize[mt->order]][pos%mt->patternSize[mt->order]].fxdata[ch]) / 538.1489198433845617116833784366;
										 mt->pitchBend[ch]=(mt->pitchBend[ch]*10+nextPB)/11;
										 }*/
									 break;
						}
						case 'N': // channel volume slide
							mt->ch[ch].vol = clamp(mt->ch[ch].vol + ((int)mt->ch[ch].fxData - 127)*0.0001, 0, 1);
							break;
						case 'P': // panning slide
							mt->ch[ch].pan = clamp(mt->ch[ch].pan + (127 - (int)mt->ch[ch].fxData)*-0.05, 0, 255);
							break;
						case 'W': // global volume slide
							mt->globalVolume = clamp(mt->globalVolume + ((int)mt->ch[ch].fxData - 127)*0.0001, 0, 1);
							break;

					}
				}
				mt->frameTimerFx -= 0.005*(60.0 / mt->diviseur) * mt->sampleRate / mt->tempo;
			}
			mt->frameTimerFx++;
		}


		for (unsigned ch = 0; ch < FM_ch; ++ch)
		{
			if (!mt->ch[ch].active)
				continue;


			mt->ch[ch].pan = (mt->ch[ch].pan*(mt->transitionSpeed - 1) + mt->ch[ch].destPan) / mt->transitionSpeed;
			//mt->ch[ch].vol = (mt->ch[ch].vol*(speed-1)+mt->ch[ch].destVol)/speed;


			// Update lfo
			if (mt->ch[ch].lfoDelayCpt++ >= mt->ch[ch].lfoDelayCptMax)
			{
				mt->ch[ch].lfoPhase += mt->ch[ch].lfoIncr;
				mt->ch[ch].lfoEnv += (1.f - mt->ch[ch].lfoEnv)*mt->ch[ch].lfoA;
				mt->ch[ch].lfo = wavetable[mt->ch[ch].lfoWaveform][((mt->ch[ch].lfoPhase & mt->ch[ch].lfoMask) >> 10) % LUTsize] * mt->ch[ch].lfoEnv;
			}
			int opOutUsed = 0;
			mt->ch[ch].currentEnvLevel = 0;
			for (unsigned op = 0; op < FM_op; ++op)
			{
				fm_operator* o = &mt->ch[ch].op[op];
				if (o->connectOut != &mt->noConnect)
				{
					opOutUsed += mt->ch[ch].op[o->id].state;
					mt->ch[ch].currentEnvLevel += mt->ch[ch].op[o->id].env;
				}

				/* Handle envelope */

				switch (o->state)
				{
					/* Delay */
					case 1:
						if (o->envCount++ >= o->delay)
						{

							if (mt->ch[ch].instr->op[op].pitchInitialRatio > 0)
								o->pitchMod = 1 + expVol[mt->ch[ch].instr->op[op].pitchInitialRatio] * expVol[mt->ch[ch].instr->op[op].pitchInitialRatio] * 12;
							else if (mt->ch[ch].instr->op[op].pitchInitialRatio < 0)
								o->pitchMod = 1 + (float)mt->ch[ch].instr->op[op].pitchInitialRatio*_99TO1;
							else
								o->pitchMod = 1;

							o->pitchTime = expEnv[mt->ch[ch].instr->op[op].pitchDecay] * mt->sampleRateRatio;
							o->pitchDestRatio = 1;

							if (mt->ch[ch].instr->phaseReset || o->env < 0.1)
							{
								o->phase = o->offset;
							}

							if (o->envCount >= 99999999)
							{
								o->env = o->s;
							}
							else if (mt->ch[ch].instr->envReset)
								o->env = o->i;



							o->env += (1.4f - o->env) * o->a;
							if (o->env >= 1.f)
							{
								o->env = 1.f;
								o->state = o->h > 0 ? 3 : 4;
							}
							else
								o->state = 2;
						}
						break;
						/* Attack */
					case 2:
						o->env += (1.4f - o->env) * o->a;
						if (o->env >= 1.f)
						{
							o->env = 1.f;
							o->state = o->h > 0 ? 3 : 4;
						}
						break;
						/* Hold */
					case 3:
						if (o->envCount++ >= o->h)
							o->state++;
						break;
						/* Decay - Sustain */
					case 4:
						o->env -= (o->env - o->s) * o->d;
						if (o->env - o->s < 0.001f)
						{
							o->env = o->s;


							if (o->s < 0.001f)
							{
								if (o->envLoop)
								{
									o->envCount = 99999999;
									o->state = 1;
								}
								else
								{
									o->state = o->env = o->amp = 0;
								}
							}
							else
							{
								o->envCount = 99999999;
								if (o->envLoop)
								{

									o->state = 1;
								}
								else
								{
									o->state = 5;
								}
							}
						}

						break;
						/* Release */
					case 6:
						o->env *= o->r;

						if (o->r <= 1)
						{
							if (o->env < 0.001f)
								o->state = o->env = o->amp = 0;
						}
						else
						{
							if (o->env >= 1.f)
							{
								o->env = 1.f;
								o->state = 5;
							}
						}
						break;
				}

				o->pitchMod -= (o->pitchMod - o->pitchDestRatio)*o->pitchTime;
				o->ampDelta = (o->env * o->vol *(1.f - mt->ch[ch].lfo * o->lfoAM) - o->amp) / 8;
				//o->amp = o->env * o->vol *(1.f - mt->ch[ch].lfo * o->lfoAM );
				o->pitch = o->incr * o->pitchMod * mt->ch[ch].pitchBend *(1 + mt->ch[ch].lfo * o->lfoFM);

			}
			mt->ch[ch].active = opOutUsed;
		}

		/* Previous stuff didnt need to be updated for every sample, we do 8 rendering steps for 1 update step to save CPU */

		for (unsigned iter = 0; iter < 8; iter++)
		{
			float rendu = 0, renduL = 0, renduR = 0, fxL = 0, fxR = 0;

			for (unsigned ch = 0; ch < FM_ch; ++ch)
			{
				if (!mt->ch[ch].active || mt->ch[ch].muted)
					continue;

				/* FM calculations, unrolled to be sure the compiler doesn't generate a loop */

				mt->ch[ch].op[0].phase += mt->ch[ch].op[0].pitch;
				mt->ch[ch].op[0].amp += mt->ch[ch].op[0].ampDelta;
				mt->ch[ch].op[0].out = mt->ch[ch].op[0].waveform[((mt->ch[ch].op[0].phase >> 10) + (unsigned)*mt->ch[ch].op[0].connect + (unsigned)*mt->ch[ch].op[0].connect2 + (unsigned)(*mt->ch[ch].feedbackSource*mt->ch[ch].feedbackLevel)) % LUTsize] * mt->ch[ch].op[0].amp;


				mt->ch[ch].op[1].phase += mt->ch[ch].op[1].pitch;
				mt->ch[ch].op[1].amp += mt->ch[ch].op[1].ampDelta;
				mt->ch[ch].op[1].out = mt->ch[ch].op[1].waveform[((mt->ch[ch].op[1].phase >> 10) + (unsigned)*mt->ch[ch].op[1].connect + (unsigned)*mt->ch[ch].op[1].connect2) % LUTsize] * mt->ch[ch].op[1].amp;


				mt->ch[ch].op[2].phase += mt->ch[ch].op[2].pitch;
				mt->ch[ch].op[2].amp += mt->ch[ch].op[2].ampDelta;
				mt->ch[ch].op[2].out = mt->ch[ch].op[2].waveform[((mt->ch[ch].op[2].phase >> 10) + (unsigned)*mt->ch[ch].op[2].connect + (unsigned)*mt->ch[ch].op[2].connect2) % LUTsize] * mt->ch[ch].op[2].amp;


				mt->ch[ch].op[3].phase += mt->ch[ch].op[3].pitch;
				mt->ch[ch].op[3].amp += mt->ch[ch].op[3].ampDelta;
				mt->ch[ch].op[3].out = mt->ch[ch].op[3].waveform[((mt->ch[ch].op[3].phase >> 10) + (unsigned)*mt->ch[ch].op[3].connect + (unsigned)*mt->ch[ch].op[3].connect2) % LUTsize] * mt->ch[ch].op[3].amp;


				mt->ch[ch].op[4].phase += mt->ch[ch].op[4].pitch;
				mt->ch[ch].op[4].amp += mt->ch[ch].op[4].ampDelta;
				mt->ch[ch].op[4].out = mt->ch[ch].op[4].waveform[((mt->ch[ch].op[4].phase >> 10) + (unsigned)*mt->ch[ch].op[4].connect + (unsigned)*mt->ch[ch].op[4].connect2) % LUTsize] * mt->ch[ch].op[4].amp;


				mt->ch[ch].op[5].phase += mt->ch[ch].op[5].pitch;
				mt->ch[ch].op[5].amp += mt->ch[ch].op[5].ampDelta;
				mt->ch[ch].op[5].out = mt->ch[ch].op[5].waveform[((mt->ch[ch].op[5].phase >> 10) + (unsigned)*mt->ch[ch].op[5].connect + (unsigned)*mt->ch[ch].op[5].connect2) % LUTsize] * mt->ch[ch].op[5].amp;


				mt->ch[ch].mixer = *mt->ch[ch].op[0].toMix + *mt->ch[ch].op[1].toMix + *mt->ch[ch].op[2].toMix + *mt->ch[ch].op[3].toMix;

				rendu = (*mt->ch[ch].op[0].connectOut + *mt->ch[ch].op[1].connectOut + *mt->ch[ch].op[2].connectOut + *mt->ch[ch].op[3].connectOut + *mt->ch[ch].op[4].connectOut + *mt->ch[ch].op[5].connectOut)*mt->ch[ch].vol*mt->ch[ch].instrVol;

				mt->ch[ch].lastRender2 = mt->ch[ch].lastRender;
				mt->ch[ch].lastRender = rendu;

				/* Is a smooth transition needed between two notes ? */

				if (mt->ch[ch].fade > 0.00001)
				{
					rendu = rendu*(1 - mt->ch[ch].fade) + mt->ch[ch].fadeFrom*mt->ch[ch].fade;
					mt->ch[ch].fadeFrom += mt->ch[ch].delta*mt->ch[ch].fade;
					mt->ch[ch].fade *= mt->ch[ch].fadeIncr;
				}

				float trenduL = rendu*wavetable[0][LUTsize / 4 + (unsigned)mt->ch[ch].pan*LUTratio];
				float trenduR = rendu*wavetable[0][(unsigned)mt->ch[ch].pan*LUTratio];


				renduL += trenduL;
				renduR += trenduR;
				fxL += trenduL*mt->ch[ch].reverbSend;
				fxR += trenduR*mt->ch[ch].reverbSend;
			}

			/* Reverb phases */

			unsigned prevPhaseL = mt->reverbPhaseL;
			mt->reverbPhaseL = (mt->reverbPhaseL + 1) % mt->reverbMod1;
			unsigned prevPhaseL2 = mt->reverbPhaseL2;
			mt->reverbPhaseL2 = (mt->reverbPhaseL2 + 1) % mt->reverbMod2;
			unsigned prevPhaseR = mt->reverbPhaseR;
			mt->reverbPhaseR = (mt->reverbPhaseR + 1) % mt->reverbMod3;
			unsigned prevPhaseR2 = mt->reverbPhaseR2;
			mt->reverbPhaseR2 = (mt->reverbPhaseR2 + 1) % mt->reverbMod4;

			/* Two comb filters, left */

			mt->outL = ((mt->revBuf[mt->reverbPhaseL] + mt->revBuf[mt->reverbMod1 + mt->reverbPhaseL2]))*0.5;
			mt->revBuf[mt->reverbPhaseL] =  fxR + (mt->revBuf[mt->reverbPhaseL] + mt->revBuf[prevPhaseL])*0.5*mt->reverbLength;
			mt->revBuf[mt->reverbMod1 + mt->reverbPhaseL2] =fxL + (mt->revBuf[mt->reverbMod1 + mt->reverbPhaseL2] + mt->revBuf[mt->reverbMod1 + prevPhaseL2])*0.5*mt->reverbLength;

			/* Two comb filters, right */

			mt->outR = ((mt->revBuf[mt->revOffset2 + mt->reverbPhaseR] + mt->revBuf[mt->revOffset3 + mt->reverbPhaseR2]))*0.5;
			mt->revBuf[mt->revOffset2 + mt->reverbPhaseR] = fxL + (mt->revBuf[mt->revOffset2 + mt->reverbPhaseR] + mt->revBuf[mt->revOffset2 + prevPhaseR])*0.5*mt->reverbLength;
			mt->revBuf[mt->revOffset3 + mt->reverbPhaseR2] =  fxR + (mt->revBuf[mt->revOffset3 + mt->reverbPhaseR2] + mt->revBuf[mt->revOffset3 + prevPhaseR2])*0.5*mt->reverbLength;

			/* First allpass */

			float outL2 = 0.5*mt->outL + mt->revBuf[mt->revOffset4 + mt->allpassPhaseL];
			mt->revBuf[mt->revOffset4 + mt->allpassPhaseL] = mt->outL - 0.5 * outL2;
			mt->allpassPhaseL = (mt->allpassPhaseL + 1) % mt->allpassMod;

			float outR2 = 0.5*mt->outR + mt->revBuf[mt->revOffset5 + mt->allpassPhaseR];
			mt->revBuf[mt->revOffset5 + mt->allpassPhaseR] = mt->outR - 0.5 * outR2;
			mt->allpassPhaseR = (mt->allpassPhaseR + 1) % mt->allpassMod;

			/* Second allpass */

			float outL22 = 0.5*outL2 + mt->revBuf[mt->revOffset6 + mt->allpassPhaseL2];
			mt->revBuf[mt->revOffset6 + mt->allpassPhaseL2] = outL2 - 0.5 * outL22;
			mt->allpassPhaseL2 = (mt->allpassPhaseL2 + 1) % mt->allpassMod2;

			float outR22 = 0.5*outR2 + mt->revBuf[mt->revOffset7 + mt->allpassPhaseR2];
			mt->revBuf[mt->revOffset7 + mt->allpassPhaseR2] = outR2 - 0.5 * outR22;
			mt->allpassPhaseR2 = (mt->allpassPhaseR2 + 1) % mt->allpassMod2;


			/* Final mix */
			buffer[b] = (renduL + outL22) * mt->globalVolume * mt->playbackVolume;
			buffer[b + 1] =(renduR + outR22) * mt->globalVolume * mt->playbackVolume;
			b += 2;
			if (b>=length)
				return;
		}
	}
}

void mt_playNote(mtsynth* mt, unsigned _instrument, unsigned note, unsigned ch, unsigned volume)
{
	if (ch >= FM_ch || _instrument == 255 && !mt->ch[ch].instr || _instrument != 255 && _instrument >= mt->instrumentCount)
		return;

	/* Instrument changed, update parameters */
	if (_instrument != 255 && _instrument < mt->instrumentCount && mt->ch[ch].cInstr != &mt->instrument[_instrument])
	{

		mt->ch[ch].cInstr = &mt->instrument[_instrument];
		mt->ch[ch].instrNumber = _instrument;

		mt->ch[ch].instr = &mt->instrument[_instrument];
		mt->ch[ch].instrVol = expVol[mt->ch[ch].instr->volume];
		mt->ch[ch].lfoMask = lfoMasks[mt->ch[ch].instr->lfoWaveform];
		mt->ch[ch].lfoWaveform = lfoWaveforms[mt->ch[ch].instr->lfoWaveform];
		mt->ch[ch].feedbackLevel = expVol[mt->ch[ch].instr->feedback];
		mt->ch[ch].lfoA = expEnv[mt->ch[ch].instr->lfoA] * mt->sampleRateRatio;
		mt->ch[ch].lfoIncr = 1 + expVol[mt->ch[ch].instr->lfoSpeed] * expVol[mt->ch[ch].instr->lfoSpeed] * 5000 * mt->sampleRateRatio*LUTratio;
		mt->ch[ch].lfoDelayCptMax = expVol[mt->ch[ch].instr->lfoDelay] * expVol[mt->ch[ch].instr->lfoDelay] * 200000 * mt->sampleRateRatio;
		mt->ch[ch].lfoEnv = mt->ch[ch].lfoDelayCpt = mt->ch[ch].lfo = mt->ch[ch].lfoPhase = 0;
		mt->ch[ch].pitchBend = 1;
		mt->ch[ch].transpose = mt->ch[ch].instr->transpose;
		mt->ch[ch].tuning = 0.0006 * mt->ch[ch].instr->tuning;
		mt->ch[ch].lfoOffset = mt->ch[ch].instr->lfoOffset * LUTsize / 32;
		for (unsigned op = 0; op < FM_op; ++op)
		{
			fm_operator* o = &mt->ch[ch].op[op];
			o->env = 0;
			o->connectOut = (mt->ch[ch].instr->op[op].connectOut >= 0) ? &mt->ch[ch].op[mt->ch[ch].instr->op[op].connectOut].out : &mt->noConnect;
			o->id = mt->ch[ch].instr->op[op].connectOut;
			o->connect = (mt->ch[ch].instr->op[op].connect >= 0) ? &mt->ch[ch].op[mt->ch[ch].instr->op[op].connect].out : &mt->noConnect;
			o->connect2 = (mt->ch[ch].instr->op[op].connect2>5) ? &mt->ch[ch].mixer :
				(mt->ch[ch].instr->op[op].connect2 >= 0 ? &mt->ch[ch].op[mt->ch[ch].instr->op[op].connect2].out : &mt->noConnect);

			o->waveform = wavetable[mt->ch[ch].instr->op[op].waveform];
			o->lfoFM = expVol[mt->ch[ch].instr->op[op].lfoFM] * expVol[mt->ch[ch].instr->op[op].lfoFM];
			o->lfoAM = expVol[mt->ch[ch].instr->op[op].lfoAM];

			o->delay = expEnv[mt->ch[ch].instr->op[op].delay] * 3000000 / mt->sampleRateRatio;

			o->i = expVol[mt->ch[ch].instr->op[op].i];
			o->h = expEnv[mt->ch[ch].instr->op[op].h] * 700000 / mt->sampleRateRatio;
			o->s = expVol[mt->ch[ch].instr->op[op].s];
			o->r = (mt->ch[ch].instr->op[op].r >= 0) ? exp(-(expEnv[mt->ch[ch].instr->op[op].r])*mt->sampleRateRatio) : 2 - exp(-(expEnv[abs(mt->ch[ch].instr->op[op].r)])*mt->sampleRateRatio);

			o->finetune = mt->ch[ch].instr->op[op].finetune;
			o->detune = mt->ch[ch].instr->op[op].detune;
			o->mult = mt->ch[ch].instr->op[op].mult;
			o->baseVol = mt->ch[ch].instr->op[op].vol * !mt->ch[ch].instr->op[op].muted;
			o->baseA = mt->ch[ch].instr->op[op].a;
			o->baseD = mt->ch[ch].instr->op[op].d;
			o->fixedFreq = mt->ch[ch].instr->op[op].fixedFreq;
			o->offset = ((unsigned int)mt->ch[ch].instr->op[op].offset)* LUTsize * 32;
			o->envLoop = mt->ch[ch].instr->op[op].envLoop;
			o->pitchFinalRatio = mt->ch[ch].instr->op[op].pitchFinalRatio;
			o->velSensitivity = (float)mt->ch[ch].instr->op[op].velSensitivity*_99TO1;
			o->volScaling = mt->ch[ch].instr->op[op].kbdVolScaling*0.001;
			o->kbdCenterNote = mt->ch[ch].instr->op[op].kbdCenterNote;
		}
		for (unsigned op = 0; op < FM_op - 2; ++op)
		{
			mt->ch[ch].op[op].toMix = (mt->ch[ch].instr->toMix[op] >= 0) ? &mt->ch[ch].op[mt->ch[ch].instr->toMix[op]].out : &mt->noConnect;
		}
		mt->ch[ch].feedbackSource = &mt->ch[ch].op[mt->ch[ch].instr->feedbackSource].out;

	}

	/* Note changed */
	if (note < 128 && mt->ch[ch].instr)
	{
		mt->ch[ch].untransposedNote = note;


		mt_calcPitch(mt, ch, note);

		if (volume < 100)
			mt->ch[ch].noteVol = volume;

		if (mt->ch[ch].instr->flags & FM_INSTR_LFORESET)
		{
			mt->ch[ch].lfoEnv = mt->ch[ch].lfoDelayCpt = mt->ch[ch].lfo = 0;
			mt->ch[ch].lfoPhase = mt->ch[ch].lfoOffset * LUTsize / 2;
		}


		/* Trigger note transition smoothing algorithm to avoid clicks/pops */
		if (mt->ch[ch].instr->flags & FM_INSTR_SMOOTH && mt->ch[ch].currentEnvLevel > 0.1 && (mt->ch[ch].instr->envReset || mt->ch[ch].instr->phaseReset))
		{

			mt->ch[ch].fade = 1;
			mt->ch[ch].fadeFrom = mt->ch[ch].lastRender;
			mt->ch[ch].delta = clamp((mt->ch[ch].lastRender - mt->ch[ch].lastRender2), -2000, 2000)*mt->sampleRateRatio;

			mt->ch[ch].fadeIncr = 0.95 - mt->ch[ch].note*0.001;
		}

		for (unsigned op = 0; op < FM_op; ++op)
		{
			fm_operator* o = &mt->ch[ch].op[op];

			mt_calcOpVol(o, mt->ch[ch].note, volume == 255 ? mt->ch[ch].noteVol : volume);
			o->amp = 0;
			o->a = expEnv[(int)max(0, min(99, (o->baseA + mt->ch[ch].instr->op[op].kbdAScaling*((int)mt->ch[ch].note - mt->ch[ch].instr->op[op].kbdCenterNote)*0.07f)))] * mt->sampleRateRatio;
			o->d = 1 - exp(-expEnv[(int)max(0, min(99, (o->baseD + mt->ch[ch].instr->op[op].kbdDScaling*((int)mt->ch[ch].note - mt->ch[ch].instr->op[op].kbdCenterNote)*0.07f)))] * mt->sampleRateRatio);

			if (_instrument != 255)
			{
				if (mt->ch[ch].instr->envReset)
				{
					o->env = 0;
					o->out = 0;
				}

				mt->ch[ch].op0 = o->envCount = o->pitchTime = 0;
				o->pitchMod = o->pitchDestRatio = 1;
				o->state = 1;
			}
		}



	}



	mt->ch[ch].active = 1;
}

/* Creates a table containing all current pannings/volumes/tempo/time info for each row, for fast seeking */

void mt_buildStateTable(mtsynth* mt, unsigned orderStart, unsigned orderEnd, unsigned channelStart, unsigned channelEnd)
{

	orderStart = clamp(orderStart, 0, mt->patternCount);
	orderEnd = clamp(orderEnd, 0, mt->patternCount);
	channelStart = clamp(channelStart, 0, FM_ch);
	channelEnd = clamp(channelEnd, 0, FM_ch);


	for (int order = orderStart; order < orderEnd; order++)
	{

		if (order == 0)
		{
			for (unsigned ch = 0; ch < FM_ch; ch++)
			{
				mt->channelStates[order][0].pan[ch] = mt->ch[ch].initial_pan;
				mt->channelStates[order][0].vol[ch] = mt->ch[ch].initial_vol;
			}
			mt->channelStates[order][0].tempo = mt->initial_tempo;
			mt->channelStates[order][0].time = 0;
		}
		for (int j = 0; j < mt->patternSize[order]; j++)
		{

			/* Replicate previous row data (tempo/time) */
			if (j>0)
			{
				mt->channelStates[order][j].tempo = mt->channelStates[order][j - 1].tempo;
				mt->channelStates[order][j].time = mt->channelStates[order][j - 1].time + 60.f / (mt->channelStates[order][j].tempo*mt->diviseur);
			}
			else if (order > 0)
			{
				mt->channelStates[order][j].tempo = mt->channelStates[order - 1][mt->patternSize[order - 1] - 1].tempo;
				mt->channelStates[order][j].time = mt->channelStates[order - 1][mt->patternSize[order - 1] - 1].time + 60.f / (mt->channelStates[order][j].tempo*mt->diviseur);
			}
			for (unsigned ch = channelStart; ch< channelEnd; ch++)
			{
				/* Replicate previous row data (pan/vol for each channel) */
				if (j>0)
				{
					mt->channelStates[order][j].vol[ch] = mt->channelStates[order][j - 1].vol[ch];
					mt->channelStates[order][j].pan[ch] = mt->channelStates[order][j - 1].pan[ch];
				}
				else if (order > 0)
				{
					mt->channelStates[order][j].vol[ch] = mt->channelStates[order - 1][mt->patternSize[order - 1] - 1].vol[ch];
					mt->channelStates[order][j].pan[ch] = mt->channelStates[order - 1][mt->patternSize[order - 1] - 1].pan[ch];
				}

				switch (mt->pattern[order][j][ch].fx)
				{
					case 'T':
						mt->channelStates[order][j].tempo = mt->pattern[order][j][ch].fxdata == 0 ? 1 : mt->pattern[order][j][ch].fxdata;
						break;
					case 'X':
						mt->channelStates[order][j].pan[ch] = mt->pattern[order][j][ch].fxdata;
						break;
					case 'M':
						mt->channelStates[order][j].vol[ch] = mt->pattern[order][j][ch].fxdata;
						break;
				}
			}
		}
	}
	mt->channelStatesDone = 1;
}



void mt_render(mtsynth* mt, void* buffer, unsigned length, unsigned type)
{

	float *rendered = malloc(4*length); // float = 4bytes

	if (!rendered)
		return;

	_mt_render(mt, (float*)rendered, length);

	switch (type%64)
	{
		case MT_RENDER_FLOAT:
		{
			float *buf_f = buffer;
			for (unsigned i = 0; i < length; i++)
			{
				buf_f[i] = clamp(rendered[i]/32768,-1.0,1.0);
			}
			break;
		}
		case MT_RENDER_8:
		{
			if(type & MT_RENDER_PAD32)
			{
				int *buf_32 = buffer;
				for (unsigned i = 0; i < length; i++)
				{
					buf_32[i] = (clamp((signed char)(rendered[i]/256), -128,127));
				}
			}
			else
			{
				unsigned char *buf_8 = buffer;
				for (unsigned i = 0; i < length; i++)
				{
					buf_8[i] = clamp(128+rendered[i]/256, 0,255);
				}
			}
			
			break;
		}
		case MT_RENDER_16:{
			if(type & MT_RENDER_PAD32)
			{
				int *buf_32 = buffer;
				for (unsigned i = 0; i < length; i++)
				{
					buf_32[i] = clamp(rendered[i], -32768,32767);
				}
			}
			else
			{
				signed short *buf_16 = buffer;
				for (unsigned i = 0; i < length; i++)
				{
					buf_16[i] = clamp(rendered[i], -32768,32767);
				}
			}
			
			break;
		}
		case MT_RENDER_24:
		{
			unsigned char *buf_24 = buffer;

			if(type & MT_RENDER_PAD32)
			{
				for (unsigned i = 0; i < length; i++)
				{
					int val = clamp(rendered[i]*256, -8388608,8388607);

					// negative 24bit values should stay negative 32 bit values ! (negative has the top byte to 255)
					buf_24[i*4+3] = (val  < 0) ? 255 : 0;
					
					buf_24[i*4+2] = (unsigned char)((val&0x00ff0000) >> 16);
					buf_24[i*4+1] = (unsigned char)((val&0x00ff00)>>8);
					buf_24[i*4] = (unsigned char)(val & 0xff);
					
					
				}
			}
			else
			{
				for (unsigned i = 0; i < length; i++)
				{
					int val = clamp(rendered[i]*256, -8388608,8388607);

					buf_24[i*3+2] = (unsigned char)((val&0x00ff0000) >> 16);
					buf_24[i*3+1] = (unsigned char)((val&0x00ff00)>>8);
					buf_24[i*3] = (unsigned char)(val & 0xff);

				}
			}
			break;
		}
		case MT_RENDER_32:
		{
			int *buf_32 = buffer;
			for (unsigned i = 0; i < length; i++)
			{
				buf_32[i] = (signed int)clamp(((double)rendered[i]*256*256),-2147483648.f,2147483647.f);
			}
			break;
		}
	}

	free(rendered);
}

void mt_stopNote(mtsynth* mt, unsigned ch)
{
	if (ch >= FM_ch || !mt->ch[ch].active || mt->ch[ch].note > 127)
		return;


	for (unsigned op = 0; op < FM_op; ++op)
	{
		fm_operator *o = &mt->ch[ch].op[op];

		o->state = 6;
		o->pitchTime = expEnv[mt->ch[ch].instr->op[op].pitchRelease] * mt->sampleRateRatio;

		if (o->pitchFinalRatio>0)
			o->pitchDestRatio = 1 + expVol[o->pitchFinalRatio] * expVol[o->pitchFinalRatio] * 12;
		else if (mt->ch[ch].instr->op[op].pitchFinalRatio < 0)
			o->pitchDestRatio = 1 + (float)o->pitchFinalRatio*_99TO1;
		else
			o->pitchDestRatio = 1;
	}
	mt->ch[ch].note = 255;

}

void mt_stopSound(mtsynth* mt)
{
	for (unsigned ch = 0; ch < FM_ch; ++ch)
	{
		mt->ch[ch].active = 0;
		mt->ch[ch].lastRender = mt->ch[ch].lastRender2 = 0;
		mt->ch[ch].note = 255;
		mt->ch[ch].cInstr = 0;
		mt->ch[ch].instrNumber = 255;
		mt->ch[ch].currentEnvLevel = 0;
		for (unsigned op = 0; op < FM_op; ++op)
		{
			mt->ch[ch].op[op].state = mt->ch[ch].op[op].env = mt->ch[ch].op[op].amp = 0;
		}
	}
	memset(mt->revBuf, 0, mt->revBufSize*sizeof(float));
}

void mt_play(mtsynth* mt)
{
	if (mt->playing)
	{
		mt_stop(mt,1);
		mt_setPosition(mt,0,0,2);
	}
	if (!mt->channelStatesDone)
		mt_buildStateTable(mt, 0, mt->patternCount, 0, FM_ch);
	mt->playing = mt->patternCount > 0;
	mt->frameTimer = mt->frameTimerFx = 0;
	mt->tempRow = mt->tempOrder = -1;
	mt->looping = -1;
	mt->loopCount = 0;
	mt_initChannels(mt);
}


void mt_stop(mtsynth* mt, int cut)
{
	if (cut)
	{
		mt_stopSound(mt);
	}
	for (unsigned ch = 0; ch < FM_ch; ++ch)
	{
		mt_stopNote(mt, ch);
		mt->ch[ch].cInstr = 0;
	}
	mt->playing = 0;
}

void mt_setPosition(mtsynth* mt, int order, int row, int cutNotes)
{
	if (cutNotes == 1)
	for (unsigned ch = 0; ch < FM_ch; ++ch) mt_stopNote(mt, ch);
	else if (cutNotes == 2)
		mt_stopSound(mt);

	mt->order = clamp(order, 0, mt->patternCount - 1);
	mt->row = clamp(row, 0, (int)mt->patternSize[order] - 1);
	mt->frameTimer = mt->frameTimerFx = 0;
	if (mt->playing)
		mt_initChannels(mt);
}

#include <stdint.h>

static uint32_t adler32(const void *buf, size_t buflength)
{
	const uint8_t *buffer = (const uint8_t*)buf;

	uint32_t s1 = 1;
	uint32_t s2 = 0;

	for (size_t n = 0; n < buflength; n++)
	{
		s1 = (s1 + buffer[n]) % 65521;
		s2 = (s2 + s1) % 65521;
	}
	return (s2 << 16) | s1;
}

int mt_saveSong(mtsynth* mt, const char* filename)
{
	FILE *fp = fopen(filename, "wb+");
	if (!fp)
	{
		return 0;
	}
	fputc('M', fp);
	fputc('D', fp);
	fputc('T', fp);
	fputc('S', fp);
	fputc(0x00, fp); // unused byte
	fputc(MUDTRACKER_VERSION, fp); // version
	unsigned char temp = strlen(&mt->songName[0]);
	fwrite(&temp, sizeof(temp), 1, fp);
	fwrite(&mt->songName[0], temp, 1, fp);

	temp = strlen(&mt->author[0]);
	fwrite(&temp, sizeof(temp), 1, fp);
	fwrite(&mt->author[0], temp, 1, fp);

	temp = strlen(&mt->comments[0]);
	fwrite(&temp, sizeof(temp), 1, fp);
	fwrite(&mt->comments[0], temp, 1, fp);

	fwrite(&mt->initial_tempo, sizeof(mt->initial_tempo), 1, fp); // tempo
	fwrite(&mt->diviseur, sizeof(mt->diviseur), 1, fp); // quarter note
	fwrite(&mt->_globalVolume, sizeof(mt->_globalVolume), 1, fp);
	fwrite(&mt->transpose, sizeof(mt->transpose), 1, fp);

	temp = round(mt->initialReverbLength * 160);
	fwrite(&temp, sizeof(temp), 1, fp);

	temp = round(mt->initialReverbRoomSize * 160);
	fwrite(&temp, sizeof(temp), 1, fp);

	for (unsigned ch = 0; ch < FM_ch; ++ch)
	{
		fwrite(&mt->ch[ch].initial_pan, sizeof(mt->ch[ch].initial_pan), 1, fp); // ch panning
		fwrite(&mt->ch[ch].initial_vol, sizeof(mt->ch[ch].initial_vol), 1, fp); // ch volume
		fwrite(&mt->ch[ch].initial_reverb, sizeof(mt->ch[ch].initial_reverb), 1, fp); // ch volume
	}

	temp = mt->patternCount;

	fwrite(&temp, sizeof(temp), 1, fp);
	for (unsigned i = 0; i < mt->patternCount; i++)
	{
		temp = mt->patternSize[i];
		fwrite(&temp, sizeof(temp), 1, fp);
		fwrite((char*)&mt->pattern[i][0], sizeof(Cell)*mt->patternSize[i] * FM_ch, 1, fp);
	}
	fwrite(&mt->instrumentCount, 1, 1, fp);


	for (int slot = 0; slot < mt->instrumentCount; slot++)
	{
		mt->instrument[slot].version = MUDTRACKER_VERSION;
	}

	fwrite((char*)&mt->instrument[0], sizeof(fm_instrument)*mt->instrumentCount, 1, fp);

	int totalSize = ftell(fp);
	char *all = malloc(totalSize);
	fseek(fp, 0, SEEK_SET);
	if (fread(all, totalSize, 1, fp) != 1)
		return 0;

	unsigned checksum = adler32(all, totalSize);

	fseek(fp, 0, SEEK_END);
	fwrite((char*)&checksum, 4, 1, fp);

	fclose(fp);

	return 1;
}

void mt_patternClear(mtsynth* mt)
{
	if (mt->pattern && mt->patternCount > 0)
	{

		mt->row = mt->order = 0;

		for (unsigned i = 0; i < mt->patternCount; i++)
		{
			free(mt->pattern[i]);
			free(mt->channelStates[i]);
		}

		mt->patternCount = 0;
	}
}

void mt_instrumentRecovery(fm_instrument * i)
{
	i->magic[0] = 'F';
	i->magic[1] = 'M';
	i->magic[2] = 'C';
	i->magic[3] = 'I';

	i->lfoWaveform = clamp(i->lfoWaveform, 0, 19);
	i->volume = clamp(i->volume, 0, 99);
	i->feedbackSource = clamp(i->feedbackSource, 0, 5);
	i->transpose = clamp(i->transpose, -12, 12);
	i->tuning = clamp(i->tuning, -100, 100);

	int nbOuts = 0;

	for (int j = 0; j < 4; j++)
		i->toMix[j] = clamp(i->toMix[j], -1, 5);

	for (int op = 0; op < 6; op++)
	{
		i->op[op].vol = clamp(i->op[op].vol, 0, 99);
		i->op[op].delay = clamp(i->op[op].delay, 0, 70);
		i->op[op].a = clamp(i->op[op].a, 0, 99);
		i->op[op].h = clamp(i->op[op].h, 0, 80);
		i->op[op].d = clamp(i->op[op].d, 0, 99);
		i->op[op].s = clamp(i->op[op].s, 0, 99);
		i->op[op].r = clamp(i->op[op].r, -99, 99);
		if (i->op[op].fixedFreq)
			i->op[op].mult = clamp(i->op[op].mult, 0, 255);
		else
			i->op[op].mult = clamp(i->op[op].mult, 0, 40);
		i->op[op].finetune = clamp(i->op[op].finetune, 0, 24);
		i->op[op].detune = clamp(i->op[op].detune, -100, 100);
		i->op[op].waveform = clamp(i->op[op].waveform, 0, 7);
		i->op[op].offset = clamp(i->op[op].offset, 0, 31);
		i->op[op].pitchDecay = clamp(i->op[op].pitchDecay, 0, 99);
		i->op[op].pitchRelease = clamp(i->op[op].pitchRelease, 0, 99);
		i->op[op].pitchInitialRatio = clamp(i->op[op].pitchInitialRatio, -99, 99);
		i->op[op].pitchFinalRatio = clamp(i->op[op].pitchFinalRatio, -99, 99);


		i->op[op].connect = clamp(i->op[op].connect, -1, 5);
		i->op[op].connect2 = clamp(i->op[op].connect2, -1, 6);
		i->op[op].connectOut = clamp(i->op[op].connectOut, -1, 5);

		if (i->op[op].connectOut >= 0)
		{
			nbOuts++;
		}

		if (i->op[op].connect == op)
			i->op[op].connect = -1;

		if (i->op[op].connect2 == op)
			i->op[op].connect2 = -1;

		for (int op2 = 0; op2 < 6; op2++)
		{
			if (op != op2)
			{
				if (i->op[op].connect == op2 && i->op[op2].connect == op)
				{
					i->op[op2].connect = -1;
				}
				if (i->op[op].connect2 == op2 && i->op[op2].connect2 == op)
				{
					i->op[op2].connect2 = -1;
				}
			}
		}
	}
	if (nbOuts == 0)
	{
		i->op[0].connectOut = 0;
		i->op[1].connectOut = 1;
		i->op[2].connectOut = 2;
		i->op[3].connectOut = 3;
		i->op[4].connectOut = 4;
		i->op[5].connectOut = 5;
	}
}




static int readFromMemory(mtsynth *mt, char *dst, int len, char *from)
{
	if (mt->readSeek >= mt->totalFileSize)
		return 0;

	memcpy(dst, from + mt->readSeek, len);
	mt->readSeek += len;
	return 1;
}

static char* readFromMemoryPtr(mtsynth *mt, int len, char *from)
{
	if (mt->readSeek >= mt->totalFileSize)
		return 0;

	mt->readSeek += len;
	return &from[mt->readSeek - len];
}

int mt_loadSongFromMemory(mtsynth* mt, char* data, unsigned len)
{
	mt->totalFileSize = len;
	unsigned char nbOrd, nbRow, temp;
	unsigned int error = 0;

	if (mt->totalFileSize < 3 * FM_ch + 6)
	{
		return 2;
	}


	mt->readSeek = 5;
	readFromMemory(mt, (char *)&temp, 1, data);

	if (temp != MUDTRACKER_VERSION)
	{
		return MT_ERR_FILEVERSION;
	}


	mt->order = mt->row = 0;
	mt_patternClear(mt);

	readFromMemory(mt, (char *)&temp, 1, data);

	readFromMemory(mt, &mt->songName[0], temp, data);
	mt->songName[temp] = 0;

	readFromMemory(mt, (char *)&temp, 1, data);
	readFromMemory(mt, &mt->author[0], temp, data);
	mt->author[temp] = 0;

	readFromMemory(mt, (char *)&temp, 1, data);
	readFromMemory(mt, &mt->comments[0], temp, data);
	mt->comments[temp] = 0;




	readFromMemory(mt, (char *)&mt->initial_tempo, sizeof(mt->initial_tempo), data);
	mt->initial_tempo = max(1, mt->initial_tempo);

	readFromMemory(mt, (char *)&mt->diviseur, sizeof(mt->diviseur), data);
	mt->diviseur = clamp(mt->diviseur, 1, 32);

	readFromMemory(mt, (char *)&mt->_globalVolume, sizeof(mt->_globalVolume), data);

	mt_setVolume(mt, mt->_globalVolume);

	readFromMemory(mt, &mt->transpose, sizeof(mt->transpose), data);

	readFromMemory(mt, (char *)&temp, sizeof(temp), data);
	mt->initialReverbLength = (float)temp / 160;

	readFromMemory(mt, (char *)&temp, sizeof(temp), data);
	mt->initialReverbRoomSize = (float)temp / 160;

	mt_initReverb(mt, mt->initialReverbRoomSize);

	for (unsigned ch = 0; ch < FM_ch; ++ch)
	{
		mt->ch[ch].cInstr = 0;
		readFromMemory(mt, (char *)&mt->ch[ch].initial_pan, sizeof(mt->ch[ch].initial_pan), data); // ch panning

		readFromMemory(mt, (char *)&mt->ch[ch].initial_vol, sizeof(mt->ch[ch].initial_vol), data); // ch volume
		mt->ch[ch].initial_vol = min(mt->ch[ch].initial_vol, 99);

		readFromMemory(mt, (char *)&mt->ch[ch].initial_reverb, sizeof(mt->ch[ch].initial_reverb), data); // ch volume
		mt->ch[ch].initial_reverb = min(mt->ch[ch].initial_reverb, 99);
	}

	readFromMemory(mt, (char *)&nbOrd, sizeof(nbOrd), data);
	mt_resizePatterns(mt, nbOrd);

	for (unsigned i = 0; i < nbOrd; i++)
	{
		readFromMemory(mt, (char *)&nbRow, sizeof(nbRow), data);

		mt_resizePattern(mt, i, max(1, nbRow), 0);

		readFromMemory(mt, (char*)&mt->pattern[i][0], sizeof(Cell) * mt->patternSize[i] * FM_ch, data);
	}

	readFromMemory(mt, (char *)&mt->instrumentCount, 1, data);
	mt_resizeInstrumentList(mt, mt->instrumentCount);

	if (mt->instrumentCount <= 0 || mt->instrumentCount > 255)
	{
		mt_resizeInstrumentList(mt, 1);
	}
	if (mt->patternCount == 0)
	{
		mt_resizePatterns(mt, 1);
	}

	readFromMemory(mt, (char*)&mt->instrument[0], sizeof(fm_instrument) * mt->instrumentCount, data);

	unsigned checksum;

	if (!readFromMemory(mt, (char*)&checksum, 4, data))
	{
		error++;
	}

	if (checksum != adler32(data, max(0, (int)mt->totalFileSize - 4)))
	{
		error++;
	}


	if (error)
	{
		for (int i = 0; i < mt->instrumentCount; i++)
		{
			mt_instrumentRecovery(&mt->instrument[i]);
		}
		return MT_ERR_FILECORRUPTED;
	}

	mt_buildStateTable(mt, 0, mt->patternCount, 0, FM_ch);

	return 0;
}

char* mt_fileToMemory(mtsynth *mt, const char* filename)
{
	FILE *fp = fopen(filename, "rb");
	if (!fp)
	{
		return 0;
	}

	fseek(fp, 0, SEEK_END);
	mt->totalFileSize = ftell(fp);
	char *all = malloc(mt->totalFileSize);
	if (!all)
	{
		return 0;
	}

	fseek(fp, 0, SEEK_SET);
	if (fread(all, mt->totalFileSize, 1, fp) != 1)
	{
		return 0;
	}
	fclose(fp);
	mt->readSeek = 0;
	return all;
}

int mt_loadSong(mtsynth* mt, const char* filename)
{
	char *data = mt_fileToMemory(mt, filename);

	if (!data)
		return MT_ERR_FILEIO;

	int result = mt_loadSongFromMemory(mt, data, mt->totalFileSize);
	free(data);


	return result;
}

void mt_clearSong(mtsynth* mt)
{
	mt_resizePatterns(mt, 0);
	mt->order = mt->row = mt->transpose = 0;
	mt_setDefaults(mt);
	memset(mt->songName, 0, 64);
	memset(mt->author, 0, 64);
	memset(mt->comments, 0, 256);
}

void mt_createDefaultInstrument(mtsynth* mt, unsigned slot)
{
	strncpy((char*)&mt->instrument[slot].name[0], "Default", 7);
	strncpy((char*)&mt->instrument[slot].magic[0], "MDTI", 4);
	mt->instrument[slot].dummy = 0;
	mt->instrument[slot].version = MUDTRACKER_VERSION;
	for (unsigned op = 0; op < FM_op; ++op)
	{
		mt->instrument[slot].op[op].connectOut = op;
		mt->instrument[slot].op[op].connect = -1;
		mt->instrument[slot].op[op].connect2 = -1;
	}
	mt->instrument[slot].volume = 99;
	mt->instrument[slot].op[0].a = 99;
	mt->instrument[slot].op[0].mult = 1;
	mt->instrument[slot].op[0].vol = 99;
	mt->instrument[slot].op[0].r = 99;
}

int mt_resizeInstrumentList(mtsynth* mt, unsigned size)
{
	if (size > 255)
	{
		return 0;
	}
	if (mt->instrumentCount > 0 && size == 0)
	{
		mt->instrumentCount = 0;
		return 1;
	}

	fm_instrument* newI = realloc(mt->instrument, sizeof(fm_instrument)*size);

	if (!newI)
	{
		return 0;
	}

	mt->instrument = newI;

	if (size > mt->instrumentCount)
	{
		memset((char*)&mt->instrument[mt->instrumentCount], 0, (size - mt->instrumentCount)*sizeof(fm_instrument));
		for (unsigned i = mt->instrumentCount; i < size; i++)
		{
			mt_createDefaultInstrument(mt,i);
		}
	}
	mt->instrumentCount = size;
	return 1;
}

int mt_resizePatterns(mtsynth* mt, unsigned count)
{
	if (count > 256)
		return 0;

	if (mt->patternCount > 0 && count == 0)
	{
		mt->patternCount = 0;
		return 1;
	}

	if (count < mt->patternCount && mt->pattern)
	{
		for (unsigned i = count; i < mt->patternCount; i++)
		{
			free(mt->pattern[i]);
			free(mt->channelStates[i]);
		}
		if (mt->order >= count)
			mt->order = max(0, count - 1);
	}

	unsigned int* newPs = realloc(mt->patternSize, sizeof(unsigned)*count);
	Cell(** newPa)[24] = realloc(mt->pattern, sizeof(Cell*)*count);
	ChannelState** newC = realloc(mt->channelStates, sizeof(ChannelState*)*count);

	if (!newPs || !newPa || !newC)
	{
		return 0;
	}

	mt->patternSize = newPs;
	mt->pattern = newPa;
	mt->channelStates = newC;

	unsigned oldPatternCount = mt->patternCount;
	mt->patternCount = count;

	if (count > oldPatternCount)
	{
		for (unsigned i = oldPatternCount; i < count; i++)
		{
			mt->pattern[i] = 0;
			mt->channelStates[i] = 0;
			mt->patternSize[i] = 0;
			if (!mt_resizePattern(mt, i, 1, 0))
			{
				return 0;
			}
		}
	}


	mt->channelStatesDone = 0;
	return 1;
}


int mt_clearPattern(mtsynth* mt, unsigned pattern, unsigned rowStart, unsigned count)
{
	if (pattern >= mt->patternCount || rowStart > 255 || count > 256)
		return 0;
	memset(&mt->pattern[pattern][rowStart], 255, count*sizeof(Cell)*FM_ch);
	memset(&mt->channelStates[pattern][rowStart], 255, count*sizeof(ChannelState));
	mt->channelStatesDone = 0;
	return 1;
}


int mt_insertPattern(mtsynth* mt, unsigned rows, unsigned pos)
{
	{
		if (pos > mt->patternCount || !mt_resizePatterns(mt, mt->patternCount + 1))
		return 0;
	}
	Cell(*newPtr)[24] = mt->pattern[mt->patternCount - 1];
	ChannelState* newPtri = mt->channelStates[mt->patternCount - 1];

	for (unsigned i = mt->patternCount - 1; i > pos; i--)
	{
		mt->pattern[i] = mt->pattern[i - 1];
		mt->channelStates[i] = mt->channelStates[i - 1];
		mt->patternSize[i] = mt->patternSize[i - 1];
	}
	mt->pattern[pos] = newPtr;
	mt->channelStates[pos] = newPtri;

	if (!mt_resizePattern(mt, pos, rows, 0))
		return 0;
	mt_clearPattern(mt, pos, 0, rows);
	mt->channelStatesDone = 0;
	return 1;
}

int mt_removePattern(mtsynth* mt, unsigned order)
{
	if (order >= mt->patternCount)
		return 0;

	if (mt->patternCount == 1)
	{
		mt_clearPattern(mt, mt->patternCount - 1, 0, mt->patternSize[0]);
	}
	else
	{
		free(mt->pattern[order]);
		free(mt->channelStates[order]);
		for (unsigned i = order; i < mt->patternCount - 1; i++)
		{
			mt->pattern[i] = mt->pattern[i + 1];
			mt->channelStates[i] = mt->channelStates[i + 1];
			mt->patternSize[i] = mt->patternSize[i + 1];
		}
		mt->patternCount--;

	}
	mt->order = min(mt->order, mt->patternCount - 1);
	mt->row = min(mt->row, mt->patternSize[mt->order] - 1);
	mt->channelStatesDone = 0;
	return 1;
}

int mt_resizePattern(mtsynth* mt, unsigned order, unsigned size, unsigned scaleContent)
{
	if (order >= mt->patternCount || size == 0)
	{
		return 0;
	}

	int oldPatternSize = mt->patternSize[order];

	size = clamp(size, 1, 256);

	float scaleRatio = 0;
	if (scaleContent)
	{
		scaleRatio = (float)size / oldPatternSize;
	}

	/* Shrink content */
	if (scaleContent && scaleRatio < 1)
	{
		for (int i = 0; i < mt->patternSize[order]; i++)
		{
			for (int ch = 0; ch < FM_ch; ch++)
				mt->pattern[order][(unsigned)round(i*0.5)][ch] = mt->pattern[order][i][ch];
		}
	}

	Cell(*newP)[24] = realloc(mt->pattern[order], sizeof(Cell)*size*FM_ch);
	ChannelState* newC = realloc(mt->channelStates[order], sizeof(ChannelState)*size);

	if (newP && newC)
	{
		mt->pattern[order] = newP;
		mt->channelStates[order] = newC;
	}
	else
	{
		return 0;
	}

	if (size > mt->patternSize[order])
		mt_clearPattern(mt, order, mt->patternSize[order], size - mt->patternSize[order]);


	mt->patternSize[order] = size;
	mt->row = min(mt->row, mt->patternSize[mt->order] - 1);


	/* Expand content */
	if (scaleContent && scaleRatio > 1)
	{
		for (int i = oldPatternSize - 1; i >= 0; i--)
		{
			for (int ch = 0; ch < FM_ch; ch++)
				mt->pattern[order][(unsigned)(i*scaleRatio)][ch] = mt->pattern[order][i][ch];

			memset(&mt->pattern[order][(unsigned)round(i*scaleRatio) + 1], 255, sizeof(Cell)*FM_ch);
		}
	}
	mt->channelStatesDone = 0;
	return 1;
}

float mt_getTime(mtsynth* mt)
{
	if (mt->order >= mt->patternCount || mt->row > mt->patternSize[mt->order])
		return 0;

	return mt->channelStates[mt->order][mt->row].time;
}

int mt_saveInstrument(mtsynth* mt, const char* filename, unsigned slot)
{
	if (slot < mt->instrumentCount)
	{
		FILE *fp = fopen(filename, "wb");
		if (!fp)
			return 0;
		fputc('M', fp);
		fputc('D', fp);
		fputc('T', fp);
		fputc('I', fp);
		fputc(0x00, fp); // unused byte
		fputc(MUDTRACKER_VERSION, fp); // version
		fwrite((char*)&mt->instrument[slot].name[0], sizeof(fm_instrument)-6, 1, fp);
		fclose(fp);
		return 1;
	}
	return 0;
}

int mt_saveInstrumentBank(mtsynth* mt, const char *filename)
{
	FILE *fp = fopen(filename, "wb");

	if (!fp)
		return 0;

	fputc('M', fp);
	fputc('D', fp);
	fputc('T', fp);
	fputc('B', fp);
	fputc(0x00, fp); 				//unused byte
	fputc(MUDTRACKER_VERSION, fp);
	fputc(0x00, fp); 				//unused byte
	fputc(mt->instrumentCount, fp);
	for (int i = 0; i < mt->instrumentCount; ++i)
	{
		fputc('S', fp);
		fputc('L', fp);
		fputc('O', fp);
		fputc('T', fp);
		fputc(i, fp);
		fputc('M', fp);
		fputc('D', fp);
		fputc('T', fp);
		fputc('I', fp);
		fputc(0x00, fp); // unused byte
		fputc(MUDTRACKER_VERSION, fp);
		fwrite((char*)&mt->instrument[i].name[0], sizeof(fm_instrument)-6, 1, fp);
	}
	fclose(fp);
}

int mt_loadInstrumentFromMemory(mtsynth* mt, char *data, unsigned slot)
{
	if (slot >= mt->instrumentCount)
	{
		mt_resizeInstrumentList(mt, slot + 1);
	}

	readFromMemory(mt, (char*)&mt->instrument[slot].magic, 4, data);
	readFromMemory(mt, (char*)&mt->instrument[slot].dummy, 1, data);
	readFromMemory(mt, (char*)&mt->instrument[slot].version, 1, data);

	if (mt->instrument[slot].version != MUDTRACKER_VERSION)
	{
		return MT_ERR_FILEVERSION;
	}

	readFromMemory(mt, (char*)&mt->instrument[slot].name[0], sizeof(fm_instrument)-6, data);

	return 0;
}


int mt_loadInstrument(mtsynth* mt, const char* filename, unsigned slot)
{

	char *data = mt_fileToMemory(mt, filename);

	if (!data)
		return MT_ERR_FILEIO;

	int result = mt_loadInstrumentFromMemory(mt, data, slot);
	free(data);

	return result;
}

int mt_loadInstrumentBankFromMemory(mtsynth* mt, char *data)
{
	char magic_check[4] = {0, 0, 0, 0};
	uint8_t version = 0;
	uint8_t instruments = 0;
	uint8_t slot = 0;

	// check bank magic
	readFromMemory(mt, magic_check, 4, data);
	if (memcmp(magic_check, "MDTB", 4) != 0)
	{
		return MT_ERR_FILECORRUPTED;
	}
	readFromMemory(mt, (char *)&version, 1, data); //padding
	readFromMemory(mt, (char *)&version, 1, data); //version
	if (version != MUDTRACKER_VERSION)
	{
		return MT_ERR_FILEVERSION;
	}
	readFromMemory(mt, (char *)&instruments, 1, data); //padding
	readFromMemory(mt, (char *)&instruments, 1, data); //instrument count
	if (8 /* header */ + (instruments * (5 /* 'SLOT' + # */ + sizeof(fm_instrument))) > mt->totalFileSize)
	{
		return MT_ERR_FILECORRUPTED;
	}

	mt_resizeInstrumentList(mt, 0);

	for (int i = 0; i < instruments; ++i)
	{
		readFromMemory(mt, magic_check, 4, data);
		if (memcmp(magic_check, "SLOT", 4) != 0)
		{
			return MT_ERR_FILECORRUPTED;
		}
		readFromMemory(mt, (char *)&slot, 1, data);
		if (slot >= instruments)
		{
			return MT_ERR_FILECORRUPTED;
		}
		if (mt_loadInstrumentFromMemory(mt, data, slot) != 0)
		{
			return MT_ERR_FILECORRUPTED;
		}
	}

	return 0;

}

int mt_loadInstrumentBank(mtsynth* mt, const char* filename)
{
	char *data = mt_fileToMemory(mt, filename);

	if (!data)
		return MT_ERR_FILEIO;

	int result = mt_loadInstrumentBankFromMemory(mt, data);
	free(data);

	return result;
}

void mt_removeInstrument(mtsynth* mt, unsigned slot, int removeOccurences)
{
	if (slot >= mt->instrumentCount)
		return;

	if (removeOccurences)
	{
		for (unsigned i = 0; i < mt->patternCount; i++)
		{
			for (unsigned j = 0; j < mt->patternSize[i]; j++)
			{
				for (unsigned ch = 0; ch < FM_ch; ch++)
				{
					if (mt->pattern[i][j][ch].instr == slot)
					{
						mt->pattern[i][j][ch].instr = mt->pattern[i][j][ch].vol = mt->pattern[i][j][ch].note = mt->pattern[i][j][ch].fx = mt->pattern[i][j][ch].fxdata = 255;
					}
					else if (mt->pattern[i][j][ch].instr < 255 && mt->pattern[i][j][ch].instr > slot)
					{
						mt->pattern[i][j][ch].instr--;
					}
				}
			}
		}
	}

	if (mt->instrumentCount == 1)
		return;

	for (unsigned i = slot; i < mt->instrumentCount - 1; i++)
	{
		mt->instrument[i] = mt->instrument[i + 1];
	}
	mt_resizeInstrumentList(mt, mt->instrumentCount - 1);
}

void mt_setVolume(mtsynth *mt, int volume)
{
	volume = clamp(volume, 0, 99);
	mt->_globalVolume = volume;
	mt->globalVolume = expVol[volume] * 4096 / LUTsize;
}


void mt_setPlaybackVolume(mtsynth *mt, int volume)
{
	mt->playbackVolume = expVol[volume] * 4096 / LUTsize;
}

void mt_getPosition(mtsynth* mt, int *order, int *row)
{
	*order = mt->order;
	*row = mt->row;
}

void mt_setTime(mtsynth* mt, int time, int cutNotes)
{
	for (int i = 0; i < mt->patternCount; i++)
	{
		for (int j = 0; j < mt->patternSize[i]; j++)
		{
			if (mt->channelStates[i][j].time >= time)
			{
				mt_setPosition(mt, i, j, cutNotes);
				return;
			}
		}
	}
	mt_setPosition(mt, mt->patternCount - 1, mt->patternSize[mt->patternCount - 1] - 1, cutNotes);
}

void mt_movePattern(mtsynth* mt, int from, int to)
{
	if (from < 0 || from >= mt->patternCount || to<0 || to >= mt->patternCount)
		return;

	/* Move patterns by pairs until the elements are in the right position */
	for (int j = from; j >to; j--)
	{

		Cell(*ptr)[24] = mt->pattern[j];
		mt->pattern[j] = mt->pattern[j - 1];
		mt->pattern[j - 1] = ptr;

		ChannelState* ptr2 = mt->channelStates[j];
		mt->channelStates[j] = mt->channelStates[j - 1];
		mt->channelStates[j - 1] = ptr2;

		unsigned int size = mt->patternSize[j];
		mt->patternSize[j] = mt->patternSize[j - 1];
		mt->patternSize[j - 1] = size;
	}

	for (int j = from; j < to; j++)
	{

		Cell(*ptr)[24] = mt->pattern[j];
		mt->pattern[j] = mt->pattern[j + 1];
		mt->pattern[j + 1] = ptr;

		ChannelState* ptr2 = mt->channelStates[j];
		mt->channelStates[j] = mt->channelStates[j + 1];
		mt->channelStates[j + 1] = ptr2;

		unsigned int size = mt->patternSize[j];
		mt->patternSize[j] = mt->patternSize[j + 1];
		mt->patternSize[j + 1] = size;
	}
	mt->channelStatesDone = 0;
}

void mt_moveChannels(mtsynth* mt, int from, int to)
{
	if (from < 0 || from >= FM_ch || to < 0 || to >= FM_ch)
		return;

	/* Move pattern contents */
	for (int i = 0; i < mt->patternCount; i++)
	{
		for (int j = 0; j < mt->patternSize[i]; j++)
		{

			for (int ch = from; ch >to; ch--)
			{

				Cell ptr = mt->pattern[i][j][ch];
				mt->pattern[i][j][ch] = mt->pattern[i][j][ch - 1];
				mt->pattern[i][j][ch - 1] = ptr;

				unsigned char ptr2 = mt->channelStates[i][j].pan[ch];
				mt->channelStates[i][j].pan[ch] = mt->channelStates[i][j].pan[ch - 1];
				mt->channelStates[i][j].pan[ch - 1] = ptr2;

				unsigned char ptr3 = mt->channelStates[i][j].vol[ch];
				mt->channelStates[i][j].vol[ch] = mt->channelStates[i][j].vol[ch - 1];
				mt->channelStates[i][j].vol[ch - 1] = ptr3;
			}

			for (int ch = from; ch < to; ch++)
			{

				Cell ptr = mt->pattern[i][j][ch];
				mt->pattern[i][j][ch] = mt->pattern[i][j][ch + 1];
				mt->pattern[i][j][ch + 1] = ptr;

				unsigned char ptr2 = mt->channelStates[i][j].pan[ch];
				mt->channelStates[i][j].pan[ch] = mt->channelStates[i][j].pan[ch + 1];
				mt->channelStates[i][j].pan[ch + 1] = ptr2;

				unsigned char ptr3 = mt->channelStates[i][j].vol[ch];
				mt->channelStates[i][j].vol[ch] = mt->channelStates[i][j].vol[ch + 1];
				mt->channelStates[i][j].vol[ch + 1] = ptr3;

			}

		}
	}

	mt_stopSound(mt);

	/* Move channels */

	for (int ch = from; ch > to; ch--)
	{

		fm_channel channel = mt->ch[ch];
		mt->ch[ch] = mt->ch[ch - 1];
		mt->ch[ch - 1] = channel;
	}

	for (int ch = from; ch < to; ch++)
	{

		fm_channel channel = mt->ch[ch];
		mt->ch[ch] = mt->ch[ch + 1];
		mt->ch[ch + 1] = channel;
	}


	mt->channelStatesDone = 0;
}

void mt_setChannelVolume(mtsynth *mt, int channel, int volume)
{
	volume = clamp(volume, 0, 99);
	mt->ch[channel].initial_vol = volume;
	mt->ch[channel].vol = expVol[volume];
	mt->channelStatesDone = 0;
}
void mt_setChannelPanning(mtsynth *mt, int channel, int panning)
{
	panning = clamp(panning, 0, 255);
	mt->ch[channel].initial_pan = panning;
	mt->ch[channel].destPan = panning;
	mt->channelStatesDone = 0;
}

void mt_setChannelReverb(mtsynth *mt, int channel, int reverb)
{
	reverb = clamp(reverb, 0, 99);
	mt->ch[channel].initial_reverb = reverb;
	mt->ch[channel].reverbSend = expVol[reverb];
}

void mt_setTempo(mtsynth* mt, int tempo)
{
	tempo = clamp(tempo, 1, 255);
	mt->tempo = mt->initial_tempo = tempo;
	mt->channelStatesDone = 0;
}

float mt_getSongLength(mtsynth* mt)
{
	if (mt->patternCount == 0)
		return 0;

	if (!mt->channelStatesDone)
		mt_buildStateTable(mt, 0, mt->patternCount, 0, FM_ch);

	return mt->channelStates[mt->patternCount - 1][mt->patternSize[mt->patternCount - 1] - 1].time + 1.0 / mt->channelStates[mt->patternCount - 1][mt->patternSize[mt->patternCount - 1] - 1].tempo*(60.0 / mt->diviseur);
}

float mt_volumeToExp(int volume)
{
	return expVol[volume];
}

int mt_write(mtsynth *mt, unsigned pattern, unsigned row, unsigned channel, Cell data)
{
	if (pattern >= mt->patternCount || row >= mt->patternSize[pattern] || channel >= FM_ch)
		return 0;

	struct Cell *current = &mt->pattern[pattern][row][channel];

	if (data.note != 255)
		current->note = data.note;

	if (data.instr != 255)
		current->instr = data.instr;

	if (data.vol != 255)
		current->vol = data.vol;

	if (data.fx != 255)
	{
		current->fx = data.fx;
		mt->channelStatesDone = 0;
	}

	if (data.fxdata != 255)
	{
		current->fxdata = data.fxdata;
		mt->channelStatesDone = 0;
	}

	return 1;
}

int mt_getPatternSize(mtsynth *mt, int pattern)
{
	if (pattern >= mt->patternCount)
		return 0;
	return mt->patternSize[pattern];
}

int mt_insertRows(mtsynth *mt, unsigned pattern, unsigned row, unsigned count)
{
	if (pattern >= mt->patternCount || row >= mt->patternSize[pattern] || mt->patternSize[pattern] + count > 256)
		return 0;

	if (!mt_resizePattern(mt, pattern, mt->patternSize[pattern] + count, 0))
		return 0;

	mt->channelStatesDone = 0;

	for (int i = mt->patternSize[pattern] - 1; i >= row; i--)
	{
		for (unsigned ch = 0; ch < FM_ch; ch++)
		{
			mt->pattern[mt->order][i][ch] = mt->pattern[mt->order][i - count][ch];
		}
	}

	if (!mt_clearPattern(mt, pattern, row, count))
		return 0;
	return 1;
}

int mt_removeRows(mtsynth *mt, unsigned pattern, unsigned row, unsigned count)
{
	if (pattern >= mt->patternCount || row + count > mt->patternSize[pattern])
		return 0;

	for (int i = row; i < mt->patternSize[pattern] - count; i++)
	{
		for (unsigned ch = 0; ch < FM_ch; ch++)
		{
			mt->pattern[pattern][i][ch] = mt->pattern[pattern][i + count][ch];
		}
	}

	mt->channelStatesDone = 0;

	if (!mt_resizePattern(mt, pattern, mt->patternSize[pattern] - count, 0))
		return 0;

	return 1;
}

int mt_isInstrumentUsed(mtsynth *mt, unsigned id)
{
	unsigned instrCount=0;

	for (int j = 0; j < mt->patternCount; j++)
	{
		for (int k = 0; k < mt->patternSize[j]; k++)
		{
			for (int l = 0; l < FM_ch; l++)
			{
				if (mt->pattern[j][k][l].instr == id)
					instrCount++;
			}
		}
	}
	return instrCount >0;
}

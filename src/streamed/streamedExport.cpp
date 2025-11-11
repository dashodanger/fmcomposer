#include "../mtengine/mtlib.h"
#include "streamedExport.h"
#include <stdio.h>
#include "../gui/popup/popup.hpp"
#include "../views/settings/configEditor.hpp"
#include "portaudio.h"
#include "tinyfiledialogs.h"

extern mtsynth *phanoo; 
extern Popup *popup;
extern RenderWindow* window;
extern PaStream *stream;

extern Thread thread;
extern ConfigEditor *config;

struct StreamedExport streamedExport = {};

sf::Thread waveExportThread(&waveExportFunc);

void promptStreamedExport()
{
	static const char *filterName = "*.wav";
	const char *descName = "Export as WAVE";
	const char *typeName = "WAVE file";
					
	const char *fileName = tinyfd_saveFileDialog(descName, NULL, 1, &filterName, typeName);
	if (fileName)
	{
		
		string fileNameOk = forceExtension(fileName, "wav");
		song_stop();
		Pa_StopStream(stream);
		Pa_CloseStream(stream);
		popup->show(POPUP_WORKING);
		streamedExport.fileName = streamedExport.originalFileName = fileNameOk;
		waveExportThread.launch();						
	}
}


void stopExport(){
	streamedExport.running=0;
	song_stop();
}

std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot); 
}

void exportFinished(){
	
	song_stop();
	

	/* multi-track export */
	if (streamedExport.multitrackAssoc.size() > 0)
	{
		streamedExport.multiTrackIter++;
		streamedExport.fileName = remove_extension(streamedExport.originalFileName)+"-"+std::to_string(streamedExport.multiTrackIter+1)+"."+string("wav");

		if (streamedExport.multiTrackIter < streamedExport.multitrackAssoc.size())
		{
			waveExportFunc();
		}
	}

	if (streamedExport.multitrackAssoc.size() == 0 || streamedExport.multiTrackIter >= streamedExport.multitrackAssoc.size()) {
		popup->close();
		if (!windowFocus){
			tinyfd_notifyPopup("MUDTracker", "Export finished !","info");
		}
		streamedExport.running=0;
		config->selectSoundDevice(config->approvedDeviceId,config->approvedSampleRate, config->currentLatency, true);
		Pa_StartStream( stream );
		fm->looping=-1;
		for (unsigned i = 0; i < FM_ch; i++)
		{
			fm->ch[i].muted = streamedExport.mutedChannels[i];
		}
	}
}

void exportStart(){
	streamedExport.running=1;
	mt_setPosition(fm, streamedExport.fromPattern,0,2);
	song_play();
	fm->looping=streamedExport.nbLoops; // disable loop points so we aren't stuck forever

	/* multi-track export */
	if (streamedExport.multitrackAssoc.size() > 0)
	{
		for (unsigned i = 0; i < FM_ch; i++)
		{
			fm->ch[i].muted = 1;
		}

		for (unsigned i = 0; i < streamedExport.multitrackAssoc[streamedExport.multiTrackIter].size(); i++)
		{
			if (streamedExport.mutedChannels[streamedExport.multitrackAssoc[streamedExport.multiTrackIter][i]] != 1)
				fm->ch[streamedExport.multitrackAssoc[streamedExport.multiTrackIter][i]].muted = 0;
		}
	}
}

// size in bytes for 1 sample, for each fmRenderType (8, 16, 24 bits, 32 bits, float)
int bitDepths_bytes[5] = {1,2,3,4,4};

int waveExportFunc(){

	unsigned int size=0;
	int out[16384];

	int format;

	if (streamedExport.bitDepth == MT_RENDER_FLOAT)
	{
		format =  3; // IEEE float
	}
	else
	{
		format =  1; // integer
	}

	int bits=16,channels=2,bytes_per_sample=bitDepths_bytes[streamedExport.bitDepth];
	int block_align=channels*bytes_per_sample;
	int bitrate=fm->sampleRate*channels*bytes_per_sample;
	int bits_sample=8*bytes_per_sample;
	FILE *fp = fopen(streamedExport.fileName.c_str(), "wb");
	if (!fp){
		popup->show(POPUP_SAVEFAILED);
		return 0;
	}

	fwrite("RIFF    WAVEfmt ",16,1,fp);
	fwrite((char*)&bits,4,1,fp); // SubChunk1Size 
	fwrite((char*)&format,2,1,fp); // pcm format
	fwrite((char*)&channels,2,1,fp); // nb channels
	fwrite((char*)&fm->sampleRate,4,1,fp); // sample rate
	fwrite((char*)&bitrate,4,1,fp); // byte rate =sample_rate*num_channels*bytes_per_sample
	fwrite((char*)&block_align,2,1,fp); // block align
	fwrite((char*)&bits_sample,2,1,fp); // bits/sample
	fwrite("data    ",8,1,fp);
	exportStart();

	while(fm->playing && streamedExport.running && fm->order<=streamedExport.toPattern){
		
		mt_render(fm, &out[0],16384,streamedExport.bitDepth);
		fwrite(&out[0],bitDepths_bytes[streamedExport.bitDepth]*16384,1,fp);
		size+=bitDepths_bytes[streamedExport.bitDepth]*16384;// bits per sample * num samples
		popup->sliders[0].setValue(((float)fm->order/fm->patternCount)*100);
	}
	song_stop();
	fseek(fp,40,0);
	size-=4;
	fwrite(&size,sizeof(int),1,fp);
	fseek(fp,4,0);
	size+=36;
	fwrite(&size,sizeof(int),1,fp);
	fclose(fp);
	exportFinished();
	return 1;
}
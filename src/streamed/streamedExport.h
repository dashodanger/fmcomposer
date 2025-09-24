#ifndef STREAMEDEXPORT_H
#define STREAMEDEXPORT_H

#include <string>
#include <vector>

typedef struct StreamedExport{
	int param;
	int fromPattern;
	int toPattern;
	int nbLoops;
	int bitDepth;
	bool mutedChannels[FM_ch];
	std::vector< std::vector< int> > multitrackAssoc;
	int multiTrackIter;
	int running;
	std::string fileName;
	std::string originalFileName;
}StreamedExport;

extern struct StreamedExport streamedExport;

void promptStreamedExport();

int waveExportFunc();
void stopExport();

#endif

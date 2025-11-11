#include "instrEditor.hpp"
#include "../../midi/midi.h"
#include "tinyfiledialogs.h"


void InstrEditor::instrument_load(const char *filename)
{
	if (loadInstrument(filename, instrList->value))
	{
		updateFromFM();
		updateInstrListFromFM();
		updateToFM();
		valueChanged = 1;
		addToUndoHistory();
	}
}
void InstrEditor::instrument_open()
{
	mouse.clickLock2 = 5;
	static const char * filters[4] = { "*.mdti" };
	const char *fileName = tinyfd_openFileDialog("Load an instrument", instrDir.c_str(), 1, filters, "MUDTracker instrument", false);
	if (fileName)
	{
		instrDir = dirnameOf(fileName);
		instrument_load(fileName);
	}
}
void InstrEditor::instrument_load_default_gm()
{
	mt_resizeInstrumentList(fm, 0);
	for (int i = 0; i < 128; ++i)
	{
		addInstrument(i, 0);
	}
	for (int i = 24; i < 88; ++i)
	{
		addInstrument(i, 1);
	}
	updateFromFM();
	updateInstrListFromFM();
	updateToFM();
	valueChanged = 1;
	addToUndoHistory();
}
void InstrEditor::instrument_save()
{
	mouse.clickLock2 = 5;
	static const char * filters[4] = { "*.mdti" };
	const char *fileName = tinyfd_saveFileDialog("Save an instrument", instrDir.c_str(), 1, filters, "MUDTracker instrument");
	if (fileName)
	{
		instrDir = dirnameOf(fileName);
		string fileNameOk = forceExtension(fileName, "mdti");
		mt_saveInstrument(fm, fileNameOk.c_str(), instrList->value);
	}
}

void InstrEditor::instrument_bank_load()
{
	mouse.clickLock2 = 5;
	static const char * filters[4] = { "*.mdtb" };
	const char *fileName = tinyfd_openFileDialog("Load instrument bank", instrDir.c_str(), 1, filters, "MUDTracker bank", false);
	if (fileName)
	{
		instrDir = dirnameOf(fileName);
		mt_loadInstrumentBank(fm, fileName);
		updateFromFM();
		updateInstrListFromFM();
		updateToFM();
	}
}

void InstrEditor::instrument_bank_save()
{
	mouse.clickLock2 = 5;
	static const char * filters[4] = { "*.mdtb" };
	const char *fileName = tinyfd_saveFileDialog("Save instrument bank", instrDir.c_str(), 1, filters, "MUDTracker bank");
	if (fileName)
	{
		instrDir = dirnameOf(fileName);
		string fileNameOk = forceExtension(fileName, "mdtb");
		mt_saveInstrumentBank(fm, fileNameOk.c_str());
	}
}

int InstrEditor::loadInstrument(string filename, int slot)
{

	int result = mt_loadInstrument(fm, filename.c_str(), slot);

	if (result == 0)
		return 1;

	if (result == MT_ERR_FILEIO)
	{
		popup->show(POPUP_OPENFAILED);
	}
	else if (result == MT_ERR_FILEVERSION)
	{
		popup->show(POPUP_WRONGVERSION);
	}
	else if (result == MT_ERR_FILECORRUPTED)
	{
		popup->show(POPUP_FILECORRUPTED);
	}
	return 0;
}
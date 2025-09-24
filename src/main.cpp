#include <SFML/Audio.hpp>
#include <SFML/Graphics.hpp>
#include "fmengine/fmlib.h"
#include "globalFunctions.hpp"
#include "views/pianoroll/pianoRoll.hpp"
#include "views/instrument/instrEditor.hpp"
#include "views/pattern/songEditor.hpp"
#include "views/general/generalEditor.hpp"
#include "views/settings/configEditor.hpp"
#include "midi/midi.h"
#include "portaudio.h"
#include <fstream>

#define SI_CONVERT_GENERIC
#include "SimpleIni.h"
#include <fenv.h>
#include "tinyfiledialogs.h"
#include "input/input.hpp"
#include "state.hpp"
#include "gui/mainmenu.hpp"
#include "views/pattern/songFileActions.hpp"
#include "gui/sidebar.hpp"
#include "ProgramOptions.hxx"
#include "whereami.h"


Uint32 textEntered[32];

Texture *tileset;

View patternView, patternTopView, borderView, globalView, patNumView, instrView;
void *focusedElement;

State *state;
int windowFocus = 1;


int main(int argc, char *argv[])
{
	po::parser parser;
	po::option &app_dir = parser["appdir"];
	app_dir.bind(appdir);
	std::string song_to_play;
	po::option &song = parser["song"];
	song.bind(song_to_play);
	po::option &unknown = parser[""];

	if (!parser(argc, argv) || !app_dir.was_set())
	{
		int length = wai_getExecutablePath(NULL, 0, NULL);
		appdir.resize(length);
		int dirname_length = 0;
		wai_getExecutablePath(&appdir[0], length, &dirname_length);
		appdir.resize(dirname_length);
		if (appdir.back() != '/')
			appdir.push_back('/');
	}
	else
	{
		if (appdir.back() != '/')
			appdir.push_back('/');
	}

	global_initialize();

	gui_initialize();

	/* Create main UI */
	sidebar = new Sidebar();
	menu = new Menu(14);

	/* Create pages */
	config = new ConfigEditor();
	songEditor = new SongEditor();
	generalEditor = new GeneralEditor();
	instrEditor = new InstrEditor(46);
	popup = new Popup();
	pianoRoll = new Pianoroll(0);

	updateViews(WindowWidth, WindowHeight);

	/* Start the app on the song editor */
	menu->goToPage(PAGE_SONG);

	/* Select a default sound device or the previous selected one */
	config->selectBestSoundDevice();

	/* Load recent songs list */
	config->loadRecentSongs();


	if (song.was_set())
	{
		song_load(song_to_play.c_str(), false);
	}
	else
	{
		config->loadLastSong();
	}

	while (window->isOpen())
	{

		
		midi_getEvents();

		input_update();




		/* Don't respond to any input when window lose focus */

		if (!windowFocus || mouse.clickLock2 > 0)
		{

			if (mouse.clickLock2 > 0)
				mouse.clickLock2--;
			goto drawing;
		}


		if (contextMenu)
			contextMenu->update();

		if (contextMenu == recentSongs)
		{
			int value = contextMenu->clicked();
			if (value >= 0)
			{
				lastSongOpened = recentSongs->text[value].getString().toAnsiString();
				song_load(recentSongs->text[value].getString().toAnsiString().c_str());
			}
		}
		if (contextMenu == copyPasta)
		{
			int value = contextMenu->clicked();

			switch(value)
			{
				case 0:
					copiedSliderValueOk = copiedSliderValue;
					break;
				case 1:
					copiedSlider->setValue(copiedSliderValueOk);
					break;
			}
		}
		
		if (!popup->visible)
		{

			static float dropdowntimer=0, dropdownclosetimer=0;
			/* Recent songs menu */
			if (menu->hovered() == 1)
			{
				dropdowntimer+=frameTime60;
				if (contextMenu!=recentSongs  && !mouse.clickg && dropdowntimer>10)
					recentSongs->show(menu->items[2*4].position.x - 4, 28);
			}
			/* Edit tools menu */
			else if (menu->hovered() == 17)
			{
				if (state == songEditor)
					songEditor->editTools.show(menu->items[18*4].position.x - 4, 28);
				else if (state == instrEditor)
					instrEditor->editTools.show(menu->items[18*4].position.x - 4, 28);
			}
			/* Auto close those menus when mouse out */
			else if ( (contextMenu == recentSongs && !recentSongs->hover()
				|| contextMenu == &songEditor->editTools && !songEditor->editTools.hover()
				|| contextMenu == &instrEditor->editTools && !instrEditor->editTools.hover()
				))
			{
				dropdownclosetimer+=frameTime60;
				if (dropdownclosetimer > 10)
				{
					contextMenu = NULL;
					dropdownclosetimer=0;
				}
			}
			else
			{
				dropdowntimer=0;
				dropdownclosetimer=0;
			}
			menu->update();
			// re check clicklock in case it was set by a menu action
			if (!mouse.clickLock2)
			{
				state->update();
				instrEditor->handleGlobalEvents();
				mouse.pos = input_getmouse(borderView);
				sidebar->update();
			}
		}

		popup->handleEvents();


	drawing:
		window->clear(colors[BACKGROUND]);

		state->draw();

		sidebar->draw();
		window->setView(globalView);

		showInstrumentLights();

		if (popup->visible)
		{
			popup->draw();
		}

		if (contextMenu)
		{
			contextMenu->draw();
		}

	
		window->display();

		updateMouseCursor();

	}

	global_exit();
	return(0);
}

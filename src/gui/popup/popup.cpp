#include "popup.hpp"
#include "../../streamed/streamedExport.h"
#include "../drawBatcher.hpp"

Popup *popup;

Popup::Popup() :title("", font, charSize), view(FloatRect(0.f, 0.f, windowWidth, windowHeight)), visible(0), moving(0)
{
	bg.setFillColor(colors[BACKGROUND]);
	bg.setOutlineThickness(1.f);
	bg.setOutlineColor(colors[FOCUSOUTLINE]);
	titlebar.setFillColor(colors[FOCUSOUTLINE]);
	title.setFillColor(colors[BLOCKTEXT]);
	shadow.setPosition(10, 10 - 32);
	shadow.setFillColor(colors[POPUPSHADOW]);
}



void Popup::setSize(int _w, int _h)
{
	w = _w;
	h = _h;
	bg.setSize(Vector2f(w, h));

}


void Popup::setMelodicList()
{
	lists[1].clear();
	for (int i = 0; i < 128; i++)
	{
		lists[1].add(int2str[i] + " " + string(midiProgramNames[i]), false);
	}
	lists[1].updateSize();
}

void Popup::setPercussionList()
{
	lists[1].clear();
	for (int i = 0; i < 62; i++)
	{
		lists[1].add(noteName(i + 23) + "	" + midiPercussionNames[i], false);
	}
	lists[1].updateSize();
}

void Popup::handleEvents()
{
	if (!visible)
		return;

	mouse.pos = input_getmouse(view);
	if (mouse.clickg)
	{


		// enable window move
		if (mouse.pos.x > titlebar.getPosition().x && mouse.pos.x < titlebar.getPosition().x + titlebar.getSize().x && mouse.pos.y >= titlebar.getPosition().y &&
			mouse.pos.y <= titlebar.getPosition().y + titlebar.getSize().y)
		{
			moving = 1;
			deltaX = mouse.pos.x - titlebar.getPosition().x;
			deltaY = mouse.pos.y - titlebar.getPosition().y;
		}

		// click outside popup : hide it
		if (delay>0)
		{
			close();
		}
	}
	if (!Mouse::isButtonPressed(Mouse::Left) && moving)
	{
		moving = 0;
	}

	// move the window
	if (Mouse::isButtonPressed(Mouse::Left) && moving)
	{
		mouse.cursor = CURSOR_SWAP;
		mouse.pos = input_getmouse(globalView);
		view.setCenter((float)windowWidth / 2 - (mouse.pos.x - deltaX + 1), (float)windowHeight / 2 - (mouse.pos.y + 32 - deltaY));
	}

	// actions when button clicked
	if (!moving)
	{

		for (int i = 0; i < checkboxes.size(); i++)
		{
			if (checkboxes[i].clicked())
			{
				switch (type)
				{
					case POPUP_FADE:
						if (i == 0 && checkboxes[i].checked)
						{
							checkboxes[1].checked = 0;
							checkboxes[2].checked = 0;

						}
						else if (i == 1 && checkboxes[i].checked)
						{
							checkboxes[0].checked = 0;
							checkboxes[2].checked = 0;
						}
						else
						{
							checkboxes[0].checked = 0;
							checkboxes[1].checked = 0;
						}
						if (!checkboxes[i].checked)
							checkboxes[i].checked = 1;
						break;
					case POPUP_SEARCH:
						/* Search in the whole song : uncheck other checkboxes */
						if (i == 5)
						{
							checkboxes[6].checked = checkboxes[7].checked = 0;
							if (!checkboxes[i].checked)
								checkboxes[i].checked = 1;
						}
						/* Search in current pattern : uncheck other checkboxes */
						else if (i == 6)
						{
							checkboxes[5].checked = checkboxes[7].checked = 0;
							if (!checkboxes[i].checked)
								checkboxes[i].checked = 1;
						}
						/* Search in selection : uncheck other checkboxes */
						else if (i == 7)
						{
							checkboxes[5].checked = checkboxes[6].checked = 0;
							if (!checkboxes[i].checked)
								checkboxes[i].checked = 1;
						}
						/* Auto check 'Replace' if replace types are checked */
						else if (i > 8)
						{
							if (checkboxes[i].checked)
								checkboxes[8].checked = 1;
						}

						/* Make replace note/transpose mutually exclusive */
						if (i == 14 && checkboxes[i].checked)
							checkboxes[9].checked = 0;
						if (i == 9 && checkboxes[i].checked)
							checkboxes[14].checked = 0;


						break;
					case POPUP_STREAMEDEXPORT:

						break;

					case POPUP_TRANSPOSE:
						if (i == 0)
						{
							checkboxes[0].checked = 1;
							checkboxes[1].checked = 0;
						}
						if (i == 1)
						{
							checkboxes[1].checked = 1;
							checkboxes[0].checked = 0;
						}
						break;
				}
			}
		}
		// button click actions
		for (int i = 0; i < buttons.size(); i++)
		{
			if (buttons[i].clicked())
			{

				buttonActions(i);
				break;
			}
		}
		for (int i = 0; i < lists.size(); i++)
		{
			if (lists[i].clicked())
			{

				switch (type)
				{
					case POPUP_MULTITRACKEXPORT:

						

						if (i == 0)
						{
							updateMultitrackExportList();
						}
						else if (i == 1)
						{
							streamedExport.multitrackAssoc[lists[0].value].clear();
							for (unsigned i = 0; i < lists[1].selecteds.size(); i++)
							{
								if (lists[1].selecteds[i])
								{
									streamedExport.multitrackAssoc[lists[0].value].push_back(i);
								}
							}

							

						}

						break;
					case POPUP_MIDIEXPORT:
						if (i == 0)
						{
							// update midi instr export list for melodic/percussion types
							midiExportAssocChannels[lists[0].value] == 9 ? setPercussionList() : setMelodicList();

							lists[1].select(midiExportAssoc[lists[0].value]);
							lists[2].select(midiExportAssocChannels[lists[0].value]);



						}

						if (i == 1)
						{
							midiExportAssoc[lists[0].value] = lists[1].value;
						}

						if (i == 2)
						{

							midiExportAssocChannels[lists[0].value] = lists[2].value;

							if (lists[2].value == 9)
							{
								setPercussionList();

							}
							else
							{
								setMelodicList();
							}
							lists[1].select(midiExportAssoc[lists[0].value]);

						}
						break;
					case POPUP_EFFECTS:
						updateEffectDescription();

						break;
					case POPUP_SEARCH:
						switch (i)
						{
							case 0: checkboxes[2].checked = 1; break;
							case 1: checkboxes[3].checked = 1; break;
							case 2: checkboxes[11].checked = 1; checkboxes[8].checked = 1; break;
							case 3: checkboxes[12].checked = 1; checkboxes[8].checked = 1; break;
						}

						break;


				}

			}
		}
		for (int i = 0; i < sliders.size(); i++)
		{
			if (sliders[i].update())
			{
				switch (type)
				{
					case POPUP_TEMPERAMENT:
						songModified(1);
						fm->instrument[instrList->value].temperament[i] = sliders[i].value;
						break;
					case POPUP_TRANSPOSE:


						if (i == 0)
						{
							checkboxes[0].checked = 1;
							checkboxes[1].checked = 0;
						}
						if (i == 1)
						{
							checkboxes[1].checked = 1;
							checkboxes[0].checked = 0;
						}

						updateIntervalDescription();

						break;
					case POPUP_SEARCH:
						switch (i)
						{
							case 0: checkboxes[0].checked = 1; break;
							case 1: checkboxes[1].checked = 1; break;
							case 2: checkboxes[4].checked = 1; break;
							case 3:
								checkboxes[9].checked = 1; checkboxes[8].checked = 1;

								/* Make replace note/transpose mutually exclusive */
								checkboxes[14].checked = 0;

								break;
							case 4: checkboxes[10].checked = 1; checkboxes[8].checked = 1; break;
							case 5: checkboxes[13].checked = 1; checkboxes[8].checked = 1; break;
							case 6:
								checkboxes[14].checked = 1; checkboxes[8].checked = 1;

								/* Make replace note/transpose mutually exclusive */
								checkboxes[9].checked = 0;

								break;
						}


						if (i < 3 && sliders[i].value == -1)
						{
							sliders[i].setDisplayedValueOnly("Any");
						}
						else if (i >= 3 && i < 6 && sliders[i].value == -1)
						{
							sliders[i].setDisplayedValueOnly("None");
						}

						break;
					case POPUP_FADE:
						/* Click on % volume scaling */
						if (i == 2)
						{
							checkboxes[0].checked = 1;
							checkboxes[1].checked = 0;
							checkboxes[2].checked = 0;
						}
						/* Click on fade sliders */
						if (i == 0 || i == 1)
						{
							checkboxes[0].checked = 0;
							checkboxes[1].checked = 1;
							checkboxes[2].checked = 0;
						}
						/* Click on set volume */
						if (i == 3)
						{
							checkboxes[0].checked = 0;
							checkboxes[1].checked = 0;
							checkboxes[2].checked = 1;
						}
						break;
					case POPUP_EFFECTS:
						updateEffectDescription();
						break;
					case POPUP_STREAMEDEXPORT:
						if (i == 0)
						{
							if (sliders[0].value > sliders[1].value)
							{
								sliders[1].setValue(sliders[0].value);
							}
						}
						/* Pattern : to */
						else if (i==1)
						{

							if (sliders[0].value > sliders[1].value)
							{
								sliders[0].setValue(sliders[1].value);
							}
						}
						else if (i == 3)
						{
							updateBitDepthDescription();
						}
						break;

				}


			}
		}
	}
}


void Popup::draw()
{
	
	window->setView(view);

	titlebar.setSize(Vector2f(bg.getSize().x + 2, 32));
	titlebar.setPosition(bg.getPosition().x - 1, bg.getPosition().y - 32);
	title.setPosition(bg.getPosition().x + 8, bg.getPosition().y - 32 + 8);
	if (delay > 0)
	{
		delay -= frameTime60;
		if (delay <= 0)
			visible = 0;
	}

	drawBatcher.initialize();
	drawBatcher.addItem(&shadow);
	drawBatcher.addItem(&bg);


	

	
	for (int i = 0; i < texts.size(); i++)
		drawBatcher.addItem(&texts[i]);
	for (int i = 0; i < buttons.size(); i++)
	{
		drawBatcher.addItem(&buttons[i]);
	}
	for (int i = 0; i < lists.size(); i++)
	{
		drawBatcher.addItem(&lists[i]);
	}


	for (int i = 0; i < sliders.size(); i++)
	{
		drawBatcher.addItem(&sliders[i]);
	}

	
	for (int i = 0; i < shapes.size(); i++)
	{
		drawBatcher.addItem(&shapes[i]);
	}
	drawBatcher.addItem(&titlebar);
	drawBatcher.addItem(&title);
	drawBatcher.draw();

	for (int i = 0; i < sprites.size(); i++)
		window->draw(sprites[i]);

	for (int i = 0; i < checkboxes.size(); i++)
	{
		checkboxes[i].draw();
	}
}

void Popup::close(bool pressOK)
{
	if (!visible)
	{
		return;
	}

	if (pressOK)
	{
		buttonActions(0);
	}

	visible = 0;

	/* Avoid buffered keypresses to stay */
	textEnteredCount = 0;

	// save widget states to restore them the next time the popup shows up

	if (type >= savedState.size())
	{
		savedState.resize(type + 1);
	}
	savedState[type].resize(lists.size() + sliders.size() + checkboxes.size());
	int item = 0;
	for (int i = 0; i < lists.size(); i++)
	{
		savedState[type][item] = lists[i].value;
		item++;
	}
	for (int i = 0; i < sliders.size(); i++)
	{
		savedState[type][item] = sliders[i].value;
		item++;
	}
	for (int i = 0; i < checkboxes.size(); i++)
	{
		savedState[type][item] = checkboxes[i].checked;
		item++;
	}

}



void Popup::updateIntervalDescription()
{
	texts[0].setString("Interval : " + intervals[abs(sliders[1].value)]);
}

void Popup::updateBitDepthDescription()
{
	string bitDepths[5] = {"8 bit", "16 bit","24 bit", "32 bit", "Float"};
	sliders[3].setMinMax(0,4);
	sliders[3].setDisplayedValueOnly(bitDepths[sliders[3].value]);
}

void Popup::updateMultitrackExportList()
{
	lists[1].unselectAll();
	for (unsigned i = 0; i < streamedExport.multitrackAssoc[lists[0].value].size(); i++)
	{
		lists[1].select(streamedExport.multitrackAssoc[lists[0].value][i],true,true);
	}
}

# Adding MUDTracker song playback to your software

- Copy the files from this folder to your project's folder
- Include "mtlib.h" to your C/C++ code
- Create the MDT engine with the desired playback frequency
```
mtsynth* mt = mt_create(44100);
```

- Load some song 
```
// from a file
mt_loadSong(mt,"mysong.mdts");

// from memory
mt_loadSongFromMemory(mt, char* data, unsigned length)
```
- Play !
```
mt_play(mt);
```

To get the sound output, you need to request it using the mt_render function. Usually, the framework/sound library you use will give you some callback function for that.

```
void myAudioCallback(short *out, int nbFrames){

  // nbFrames*2 because the output is stereo
  // MT_RENDER_16 will output 16 bit signed short
  // There are also MT_RENDER_24, MT_RENDER_FLOAT and other other formats available
  
  mt_render(mt,out,nbFrames*2,MT_RENDER_16);

}
```

- Change the output volume
```
mt_setPlaybackVolume(mt,int volume); // 0 to 99
```

- Get the song length in seconds
```
float songLength = mt_getSongLength(mt);
```

- Set the playback position
```
// pattern is the index of the pattern to seek at
// row is the row number
// cutMode tells how the playing notes are affected : 0 = keep playing notes, 1 = force note off, 2 = hard cut
mt_setPosition(mt, int pattern, int row, int cutMode);
```

- Get the playback position
```
int currentPattern, currentRow;
mt_getPosition(mt, &currentPattern, &currentRow);
```

- Change the song tempo
```
// tempo is in BPM, from 1 to 255. Setting may be overrided by the song if some pattern use Tempo commands
mt_setTempo(mt, int tempo);
```

- Once you are tired of this
```
mt_destroy(mt); // free resources allocated with mt_create

```

There are a lot more functions in mtlib.h, take a look at it for more informations

#ifndef WAVIO_H
#define WAVIO_H
// by Evan X. Merz
// www.thisisnotalabel.com

// Example Wav file input and output
// this was written for educational purposes, but feel free to use it for anything you like 
// as long as you credit me appropriately ("wav IO based on code by Evan Merz")

// if you catch any bugs in this, or improve upon it significantly, send me the changes
// at evan at thisisnotalabel dot com, so we can share your changes with the world

#include "global.h"
using namespace std;

/**
*  WAV File Specification
*  FROM http://ccrma.stanford.edu/courses/422/projects/WaveFormat/
*  The canonical WAVE format starts with the RIFF header:
*  0         4   ChunkID          Contains the letters "RIFF" in ASCII form
*  (0x52494646 big-endian form).
*  4         4   ChunkSize        36 + SubChunk2Size, or more precisely:
*  4 + (8 + SubChunk1Size) + (8 + SubChunk2Size)
*  This is the size of the rest of the chunk 
*  following this number.  This is the size of the 
*  entire file in bytes minus 8 bytes for the
*  two fields not included in this count:
*  ChunkID and ChunkSize.
*  8         4   Format           Contains the letters "WAVE"
*  (0x57415645 big-endian form).
*
*  The "WAVE" format consists of two subchunks: "fmt " and "data":
*  The "fmt " subchunk describes the sound data's format:
*  12        4   Subchunk1ID      Contains the letters "fmt "
*  (0x666d7420 big-endian form).
*  16        4   Subchunk1Size    16 for PCM.  This is the size of the
*  rest of the Subchunk which follows this number.
*  20        2   AudioFormat      PCM = 1 (i.e. Linear quantization)
*  Values other than 1 indicate some 
*  form of compression.
*  22        2   NumChannels      Mono = 1, Stereo = 2, etc.
*  24        4   SampleRate       8000, 44100, etc.
*  28        4   ByteRate         == SampleRate * NumChannels * BitsPerSample/8
*  32        2   BlockAlign       == NumChannels * BitsPerSample/8
*  The number of bytes for one sample including
*  all channels. I wonder what happens when
*  this number isn't an integer?
*  34        2   BitsPerSample    8 bits = 8, 16 bits = 16, etc.
*
*  The "data" subchunk contains the size of the data and the actual sound:
*  36        4   Subchunk2ID      Contains the letters "data"
*  (0x64617461 big-endian form).
*  40        4   Subchunk2Size    == NumSamples * NumChannels * BitsPerSample/8
*  This is the number of bytes in the data.
*  You can also think of this as the size
*  of the read of the subchunk following this 
*  number.
*  44        *   Data             The actual sound data.
*/
class WavFileIO
{
    public:
        WavFileIO(){}

        WavFileIO(string &tmpPath)
        {
            myPath = tmpPath;
            read();
        }

        ~WavFileIO()
        {
            //if (myData) delete [] myData;
            if (myData_short) delete [] myData_short;
        }

        string getPath() {return myPath;}
        void setPath(string &newPath) {myPath = newPath;}

        // read a wav file into this class
        bool read();

        // write out the wav file
        bool save();

        // return a printable summary of the wav file
        string getSummary();

        // I made this public so that you can toss whatever you want in here
        // maybe a recorded buffer, maybe just whatever you want
        //char* myData;
        short* myData_short;
        int myDataSize; // because the read() function has been "reinterpreted"
                        // from "char" into "short", the original size is for "char"!!!
                        // i.e. must be divided by 2!!!
        int mySampleRate;

    private:
        string	myPath;
        int 	myChunkSize;
        int	    mySubChunk1Size;
        short 	myFormat;
        short 	myChannels;
        int   	myByteRate;
        short 	myBlockAlign;
        short 	myBitsPerSample;
};

#endif

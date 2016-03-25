#include "wav_io.h"

using namespace std;

// read a wav file into this class
bool WavFileIO::read()
{
	ifstream inFile( myPath.c_str(), ios::in | ios::binary);

	//printf("Reading wav file...\n"); // for debugging only

	inFile.seekg(4, ios::beg);
	inFile.read( (char*) &myChunkSize, 4 ); // read the ChunkSize

	inFile.seekg(16, ios::beg);
	inFile.read( (char*) &mySubChunk1Size, 4 ); // read the SubChunk1Size

	//inFile.seekg(20, ios::beg);
	inFile.read( (char*) &myFormat, sizeof(short) ); // read the file format.  This should be 1 for PCM

	//inFile.seekg(22, ios::beg);
	inFile.read( (char*) &myChannels, sizeof(short) ); // read the # of channels (1 or 2)

	//inFile.seekg(24, ios::beg);
	inFile.read( (char*) &mySampleRate, sizeof(int) ); // read the samplerate

	//inFile.seekg(28, ios::beg);
	inFile.read( (char*) &myByteRate, sizeof(int) ); // read the byterate

	//inFile.seekg(32, ios::beg);
	inFile.read( (char*) &myBlockAlign, sizeof(short) ); // read the blockalign

	//inFile.seekg(34, ios::beg);
	inFile.read( (char*) &myBitsPerSample, sizeof(short) ); // read the bitspersample

	inFile.seekg(40, ios::beg);
	inFile.read( (char*) &myDataSize, sizeof(int) ); // read the size of the data

	// read the data chunk
//	myData = new char[myDataSize];
	myData_short = new short[(int)((double)myDataSize/2)];
	inFile.seekg(44, ios::beg);
//	inFile.read(myData, myDataSize);
	inFile.read(reinterpret_cast<char*>(myData_short), myDataSize);

	inFile.close(); // close the input file

	return true; // this should probably be something more descriptive
}

// write out the wav file
bool WavFileIO::save()
{
	fstream myFile (myPath.c_str(), ios::out | ios::binary);

	// write the wav file per the wav file format
	myFile.seekp (0, ios::beg); 
	myFile.write ("RIFF", 4);
	myFile.write ((char*) &myChunkSize, 4);
	myFile.write ("WAVE", 4);
	myFile.write ("fmt ", 4);
	myFile.write ((char*) &mySubChunk1Size, 4);
	myFile.write ((char*) &myFormat, 2);
	myFile.write ((char*) &myChannels, 2);
	myFile.write ((char*) &mySampleRate, 4);
	myFile.write ((char*) &myByteRate, 4);
	myFile.write ((char*) &myBlockAlign, 2);
	myFile.write ((char*) &myBitsPerSample, 2);
	myFile.write ("data", 4);
	myFile.write ((char*) &myDataSize, 4);
	myFile.write (reinterpret_cast<char*>(myData_short), myDataSize);

	return true;
}

// return a printable summary of the wav file
string WavFileIO::getSummary()
{
	stringstream ss;
	ss<<" Format: "<<myFormat
      <<"\n Channels: "<<myChannels
      <<"\n SampleRate: "<<mySampleRate
      <<"\n ByteRate: "<<myByteRate
      <<"\n BlockAlign: "<<myBlockAlign
      <<"\n BitsPerSample: "<<myBitsPerSample
      <<"\n DataSize: "<<myDataSize;
	return ss.str();
}

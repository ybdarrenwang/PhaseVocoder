#include "phasevocoder.h"

#define FRAME_LENGTH 4096
#define FRAME_SHIFT 1024

using namespace std;

void usage( const char *prog ) {
    cout << endl
        << "Usage: "<<prog<<" [options] ..."<<endl<<endl
        << "options:"<<endl
        << " -t TIME_SCALE_FACTOR  (>1 for slowdown, <1 for speedup) [default:1]"<<endl
        << " -p PITCH_SCALE_FACTOR (>1 for higher, <1 for lower)     [default:1]"<<endl
        << " -i INPUT_WAVE_FILE_PATH"<<endl
        << " -o OUTPUT_WAVE_FILE_PATH"<<endl<<endl
        << " --phaseLock           (enable phase-locking)            [default:off]"<<endl
        << " --fdTimeStretch       (enable time-stretching by"<<endl
        << "                        frequency domain manipulation)   [default:off]"<<endl
        << " --fdPitchShift        (enable pitch-shifting by"<<endl
        << "                        frequency domain manipulation)   [default:off]"<<endl<<endl;
    exit(1);
}

void readConfig(vector<string> &Args, double &ts_factor, double &ps_factor, string &input_file, string &output_file, bool &phase_lock, bool &fd_time_stretch, bool &fd_pitch_shift)
{
    for(unsigned int i=0; i<Args.size(); ++i) {
        if( Args[i] == "-t" && Args.size() > i ) {
            ++i;
            stringstream ss(Args[i]);
            ss >> ts_factor;
        }
        else if( Args[i] == "-p" && Args.size() > i )
        {
            ++i;
            stringstream ss(Args[i]);
            ss >> ps_factor;
        }
        else if( Args[i] == "-i" && Args.size() > i ) {
            ++i;
            input_file = Args[i];
        }
        else if( Args[i] == "-o" && Args.size() > i ) {
            ++i;
            output_file = Args[i];
        }
        else if( Args[i] == "--phaseLock" && Args.size() > i ) {
            phase_lock = true;
        }
        else if( Args[i] == "--fdTimeStretch" && Args.size() > i ) {
            fd_time_stretch = true;
        }
        else if( Args[i] == "--fdPitchShift" && Args.size() > i ) {
            fd_pitch_shift = true;
        }
    }
}

int main(int argc, char **argv) {
    if( argc < 2 ){
        usage( argv[0] );
        exit(0);
    }

    // Read arguments
    string input_file="", output_file="";
    double ts_factor=1, ps_factor=1;
    bool phase_lock=false, fd_time_stretch=false, fd_pitch_shift=false;

    vector<string> Args;
    for(int i=1; i<argc; ++i)
    {
        string tmpstr(argv[i]);
        Args.push_back( tmpstr );
    }
    readConfig(Args, ts_factor, ps_factor, input_file, output_file, phase_lock, fd_time_stretch, fd_pitch_shift);
    if (input_file=="" || output_file=="") {
        cerr<<"ERROR: input or output file not specified"<<endl;
        exit(1);
    }
    cout<<"[ "<<input_file<<" -> "<<output_file<<" ]"<<endl;

    // Execute
    PhaseVocoder *pv = new PhaseVocoder(FRAME_LENGTH, FRAME_SHIFT, phase_lock, fd_time_stretch, fd_pitch_shift, ts_factor, ps_factor);
    pv->ReadWave(input_file);
    pv->Analysis();
    pv->PitchShifting();
    pv->TimeStretching();
    pv->Synthesis();
    pv->WriteWave(output_file);
    delete pv;

    #ifdef DEBUG
    cout<<"Complete!"<<endl<<endl;
    #endif

    return 0;
}

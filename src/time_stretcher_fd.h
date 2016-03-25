#ifndef TIMESTRETCHER_FD_H
#define TIMESTRETCHER_FD_H

#include "time_stretcher.h"

using namespace std;

class TimeStretcherFD : public TimeStretcher
{
    public:
        TimeStretcherFD(int n, int s) : TimeStretcher(n, s) {}
        virtual ~TimeStretcherFD(){}

        virtual void UpdatePhase(vector<double> mag, vector<double> prev_phase, vector<double> phase, vector<double>& synth_ph);
        virtual void SynthesizeFrame(vector<double>&, vector<double>&, Frame*);
        virtual void Stretch(double rate, vector<Frame*>& input_spec, vector<Frame*>& output_spec, bool reset_phase);
};

#endif

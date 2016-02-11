#include "pitch_shifter.h"

// Find Max. Frequency Bin
int find_max(float* a, int n, float threshold)
{
    int i, max_bin;
    float max=0;
    for (i=0; i<n; i++)
    {
        if (a[i]>max)
        {
            max = a[i];
            max_bin = i;
        }
    }

    if (max > threshold)
        return max_bin;
    else
        return -1;
}

// ============================================================================ //
// < Fixed KTH approach >
// ============================================================================ //

float *MUL2(float r1, float i1, float r2, float i2) //multiplation of complex numbers
{
    float *ans = new float[2];

    ans[0] = r1*r2-i1*i2;
    ans[1] = r1*i2+r2*i1;

    return ans;
}

void PitchShifter::PS_KTH_new(float **Spectrum, float PitchQuotient, int frame_shift, bool reset_ph, 
        float *window_real, float *window_imag)
{
    int i, j;
    int DelFreqBin;
    float PQResidual; // for interpolation
    float energy=0, new_energy=0; // for energy preservation
    int max_bin;
    bool *mark_phasor = new bool[FFT_SIZE/2+1];
    float max_amp_real, max_amp_imag, threshold;
    //	float *ROI, *ROI_shift; // region of influence

    float **SpecTemp = new float *[2];
    SpecTemp[0] = new float[FFT_SIZE/2+1];
    SpecTemp[1] = new float[FFT_SIZE/2+1];
    float *MagTemp = new float[FFT_SIZE/2+1];
    //AfxMessageBox("start");
    // initialization
    for (i=0; i<FFT_SIZE/2+1; i++)
    {
        SpecTemp[0][i] = 0;
        SpecTemp[1][i] = 0;
        mark_phasor[i] = false;
    }

    if (reset_ph == true)
    {
        for (i=0; i<FFT_SIZE/2+1; i++)
        {
            phasor[i] = 0;
            prev_subband[i] = i;
        }
    }

    for (i=0; i<FFT_SIZE/2+1; i++)
    {
        StudentMag[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
        StudentPha[i] = atan2(Spectrum[1][i], Spectrum[0][i]);
        energy += StudentMag[i]*StudentMag[i];
    }
    //FILE* myfile = fopen("ROI.txt","w");
    // iteratively subtract window spectrum from speech spectrum
    max_bin = find_max(StudentMag, FFT_SIZE/2+1, 0);
    threshold = 0.007*StudentMag[max_bin];
    while (max_bin != -1)
    {
        // calculate phasor
        DelFreqBin = floor((float)max_bin*PitchQuotient) - max_bin;
        //		if ((float)max_bin*PitchQuotient-floor((float)max_bin*PitchQuotient) > 0.5)
        //			DelFreqBin += 1;
        PQResidual = (float)max_bin*PitchQuotient - floor((float)max_bin*PitchQuotient);

        if (mark_phasor[max_bin] == false)
        {
            mark_phasor[max_bin] = true;
            phasor[max_bin] = phasor[prev_subband[max_bin]] + (max_bin*(PitchQuotient-1))*(2.0*PI*frame_shift)/FFT_SIZE;

            while(phasor[max_bin] >= PI)
                phasor[max_bin] -= 2.0 * PI;
            while(phasor[max_bin] < -1.0*PI)
                phasor[max_bin] += 2.0 * PI;
        }

        // shift the entire window spectrum
        max_amp_real = Spectrum[0][max_bin];
        max_amp_imag = Spectrum[1][max_bin];
        for(j=0; j<FFT_SIZE/2+1; j++)
        {
            if (j-max_bin+FFT_SIZE/2 >= 0 && j-max_bin+FFT_SIZE/2 < FFT_SIZE) // subtract window spectrum from original spectrum
            {
                float *ROI = MUL2(max_amp_real, max_amp_imag, window_real[j-max_bin+FFT_SIZE/2], window_imag[j-max_bin+FFT_SIZE/2]);
                //fprintf(myfile,"%f + %f*i\n",ROI[0],ROI[1]);
                if (j+DelFreqBin >=0 && j+DelFreqBin < FFT_SIZE/2) // shift the entire window spectrum
                {
                    float *ROI_shift = MUL2(ROI[0], ROI[1], cos(phasor[max_bin]), sin(phasor[max_bin]));

                    /*					SpecTemp[0][j+DelFreqBin] += ROI_real_shift;
                                                        SpecTemp[1][j+DelFreqBin] += ROI_imag_shift;*/

                    // interpolation
                    SpecTemp[0][j+DelFreqBin] += (1-PQResidual)*ROI_shift[0];
                    SpecTemp[0][j+DelFreqBin+1] += PQResidual*ROI_shift[0];
                    SpecTemp[1][j+DelFreqBin] += (1-PQResidual)*ROI_shift[1];
                    SpecTemp[1][j+DelFreqBin+1] += PQResidual*ROI_shift[1];

                    delete [] ROI_shift;
                }

                Spectrum[0][j] -= ROI[0];
                Spectrum[1][j] -= ROI[1];
                StudentMag[j] -= ABS2(ROI[0], ROI[1]);

                delete [] ROI;
            }
        }
        max_bin = find_max(StudentMag, FFT_SIZE/2+1, threshold);
    }

    //FILE* myfile = fopen("output_spec.txt","w");
    // generate output spectrum
    for(i=0;i<FFT_SIZE/2+1;i++)
    {
        Spectrum[0][i] = SpecTemp[0][i];
        Spectrum[1][i] = SpecTemp[1][i];
        //fprintf(myfile,"%f\n",ABS2(Spectrum[0][i],Spectrum[1][i]));
    }
    //fclose(myfile);
    //AfxMessageBox("end");
    // allocate region-of-influence for sinusoid freq. tracking
    for (i=0; i<FFT_SIZE/2+1; i++)
        MagTemp[i] = ABS2(Spectrum[0][i],Spectrum[1][i]);
    ChannelGrouping(MagTemp, prev_subband);

    // energy preservation; cumulate the phasors
    for (i=0; i<FFT_SIZE/2+1; i++)
        new_energy += Spectrum[0][i]*Spectrum[0][i] + Spectrum[1][i]*Spectrum[1][i];
    float AmpQuotient = sqrt(energy/new_energy);

    for (i=0; i<FFT_SIZE/2+1; i++)
    {
        Spectrum[0][i] *= AmpQuotient;
        Spectrum[1][i] *= AmpQuotient;

        phasor[i] = phasor[prev_subband[i]]; // cumulate the phasors
    }

    delete [] SpecTemp[0];
    delete [] SpecTemp[1];
    delete [] SpecTemp;
    delete [] mark_phasor;
}

#include "vocoder_functions.h"

float VocoderFunctions::ABS2(float a,float b) {return sqrt(a*a+b*b);}

vector<float> VocoderFunctions::vectorWeightedSum(vector<float> v1, vector<float> v2, float w1, float w2) {
    if (w2==0) return v1;
    if (w1==0) return v2;
    vector<float> ans;
    if (v1.size()!=v2.size()) return ans;
    for (int i=0; i<v1.size(); ++i)
        ans.push_back(v1[i]*w1+v2[i]*w2);
    return ans;
}

/**
 * Input:  a float vector
 *         e.g. [0,1,2,1,0,3,6,5,0]
 * Output: an int vector indication the corresponding peak index of each bin;
 *         e.g. [2,2,2,2,7,7,7,7,7]
 * Definition:
 * - peak: larger then preceding 2 and following 2 elements
 * - valley: the minimum between 2 peaks
 */
vector<int> VocoderFunctions::groupChannel(vector<float>& spec) {
    vector<int> peak_idx;
    int peak_ptr = 0;
    int valley_ptr = 0;
    peak_idx.push_back(0);
    peak_idx.push_back(0);
    for (int i=2; i<spec.size()-2; ++i) {
        if (spec[i] > spec[i-1] && spec[i] > spec[i-2] && spec[i] > spec[i+1] && spec[i] > spec[i+2]) {
            // peak: overwrite from valley to current position
            for (int j=valley_ptr; j<i; ++j)
                peak_idx[j] = i;
            peak_ptr = i;
        }
        else // otherwise: update valley if necessary
            if (spec[i]<spec[valley_ptr])
                valley_ptr = i;
        peak_idx.push_back(peak_ptr);
    }
    peak_idx.push_back(peak_ptr);
    peak_idx.push_back(peak_ptr);
    return peak_idx;
}
/*
void VocoderFunctions::new_ChannelGrouping(float *Spec, int *ChannelGroupFlag, float pitch) {
    int i,j;
    int GroupBoundary, PrevGroupBoundary=0, PrevPeak=0;
    float MinTemp;

    // determine how much bins to look forward/backward
    float bin_period = FFT_SIZE/(SamplingRate/pitch); // f0 interval width counted by bin number
    int find_max_range = (int)floor(bin_period*0.5);
    bool find_max;

    // determine region of influence
    for(i=find_max_range; i<FFT_SIZE/2+1-find_max_range; i++)
    {
        // find local mixima
        find_max = true;
        for (j=i-find_max_range; j<i+find_max_range; j++)
        {
            if (Spec[i]<Spec[j])
                find_max = false;
        }

        // determine region of influence
        if (find_max == true)
        {
            if (PrevPeak > 0)
            {
                // find minima between present and previous peak
                MinTemp = Spec[PrevPeak];
                for (j=PrevPeak+1; j<i; j++)
                {
                    if (Spec[j] < MinTemp)
                    {
                        MinTemp = Spec[j];
                        GroupBoundary = j;
                    }
                }

                // annotate the previous group around the previous peak
                for (j = PrevGroupBoundary; j<GroupBoundary; j++)
                    ChannelGroupFlag[j] = PrevPeak;

                PrevGroupBoundary = GroupBoundary;
            }
            else // annotate the first group
            {
                for (j=0; j<i; j++)
                    ChannelGroupFlag[j] = i;
            }
            PrevPeak = i;
        }
    }

    // annotate the last group
    for (j = PrevGroupBoundary; j<FFT_SIZE/2+1; j++)
        ChannelGroupFlag[j] = PrevPeak;
}

// ============================================================================ //
// < Voice cut-off frequency estimation >
// Reference: Applying the Harmonic Plus Noise Model in Concatenative Speech Synthesis, page3
// ============================================================================ //
int VocoderFunctions::VCO_estimation(float *Mag, int *harmonic_peak, float f0)
{
    // parameter setting
    float VCO_lb = 1000, VCO_ub = 8000; // the upper/lower bound of VCO
    int n = 3; // if n consecutive harmonics are unvoiced, then set VCO

    int i, j, harmonic_no=0, bin_index, VCO;
    int ChannelGroup[FFT_SIZE/2+1];
    float max_minor_mag, central_freq, residual;

    // record the bin index of harmonic bins
    for (i=0; i<FFT_SIZE/2+1; i++)
    {
        if (i == harmonic_peak[i])
            harmonic_no++;
    }

    int *harmonic_bin = new int[harmonic_no];
    bool *harmonic_voiced = new bool[harmonic_no];
    for (i=0; i<harmonic_no; i++)
        harmonic_voiced[i] = 0;

    j = 0;
    for (i=0; i<FFT_SIZE/2+1; i++)
    {
        if (i == harmonic_peak[i])
        {
            harmonic_bin[j] = i;
            j++;
        }
    }

    // determine how much bins to look forward/backward
    float bin_period = FFT_SIZE/(SamplingRate/f0); // f0 interval width counted by bin number
    int find_max_range = (int)floor(bin_period*0.5);

    // determine voiced/unvoiced harmonic peaks using peakiness
    ChannelGrouping(Mag, ChannelGroup); // locate minor peaks
    for (i=0; i<harmonic_no; i++)
    {
        // find max{A(omega-i)}
        max_minor_mag = 0;
        for (j=-1*find_max_range; j<=find_max_range; j++)
        {
            bin_index = harmonic_bin[i]+j;
            if (bin_index == ChannelGroup[bin_index] && j!=0)
            {
                if (Mag[bin_index] > max_minor_mag)
                    max_minor_mag = Mag[bin_index];
            }
        }

        if (Mag[harmonic_bin[i]]/max_minor_mag > 4 // A(omega-c)-max(A(omega-i)) > 6dB
                || max_minor_mag == 0) // there's no other peaks
            harmonic_voiced[i] = 1;
    }

    // determine voiced/unvoiced harmonic peaks using the location of harmonics
    //	myfile=fopen("harmonic-2.txt","w");
    for (i=0; i<harmonic_no; i++)
    {
        central_freq = SamplingRate*((float)harmonic_bin[i]/FFT_SIZE);
        residual = central_freq/f0 - floor(central_freq/f0);
        if (residual > 0.5)
            residual = 1-residual;

        if (residual < 0.1)
        {
            harmonic_voiced[i] = 1;
            //			fprintf(myfile,"%d, %d\n",harmonic_bin[i], 1);
        }
        //		else
        //			fprintf(myfile,"%d, %d\n",harmonic_bin[i], 0);
    }
    //	fclose(myfile);

    // set the harmonics beneath lower bound be voiced
    // i = 0;
    //    while(harmonic_bin[i] < (VCO_lb/SamplingRate)*FFT_SIZE)
    //    {
    //    harmonic_voiced[i] = 1;
    //    i++;
    //    }
    // set the harmonics exceed upper bound be unvoiced
    //i = harmonic_no-1;
    //while(harmonic_bin[i] > (VCO_ub/SamplingRate)*FFT_SIZE)
    //{
    //harmonic_voiced[i] = 0;
    //i--;
    //}

    // if n consecutive harmonics are unvoiced, then set VCO
    i = 0;
    while(harmonic_voiced[i] == 1 && i+n-1<harmonic_no)
    {
        i++;
        for (j=0; j<n; j++)
        {
            if (harmonic_voiced[i+j] == 1)
                harmonic_voiced[i] = 1;
        }
    }
    //	VCO = SamplingRate*(float)harmonic_bin[i-1]/FFT_SIZE;
    VCO = harmonic_bin[i-1];
    if (VCO < (VCO_lb/SamplingRate)*FFT_SIZE || // lower bound
            VCO > (VCO_ub/SamplingRate)*FFT_SIZE) // upper bound
        VCO = 0;

    //AfxMessageBox("VCO");

    delete [] harmonic_bin;
    delete [] harmonic_voiced;

    return VCO;
}

void VocoderFunctions::VCO_contour(float ***Spectrum, int n, float *pitch)
{
    int i, j, counter;
    float a,average_bin_th = 0;
    int *subband;

    voiced_bin_th = new int[n];

    // calculate VCO
    for (i=0; i<n; i++)
    {
        for (j=0; j<FFT_SIZE/2+1; j++)
            StudentMag[j] = ABS2(Spectrum[i][0][j],Spectrum[i][1][j]);

        subband = new int[FFT_SIZE/2+1];
        new_ChannelGrouping(StudentMag, subband, pitch[i]);
        voiced_bin_th[i] = VCO_estimation(StudentMag, subband, pitch[i]);
        delete [] subband;
    }

    // VCO contour smoothing
    for (i=0; i<n; i++)
    {
        if (voiced_bin_th[i] > 0)
        {
            average_bin_th += voiced_bin_th[i];
            counter++;
        }
    }
    average_VCO/=counter;

    // interpolation
    for (i=0; i<n; i++)
    {
        if (voiced_bin_th[i] == 0)
        {
            if (i == 0) // if the first voiced_bin_th has no value
            {
                while(voiced_bin_th[i] == 0)
                    i++;
                while(i > 0 && i < n)
                {
                    voiced_bin_th[i-1] = voiced_bin_th[i];
                    i--;
                }	
            }	
            else // interpolation
            {
                counter = 0;
                while(i+counter < n && voiced_bin_th[i+counter] == 0)
                    counter++;

                if (i+counter < n)
                    a = (float)(voiced_bin_th[i+counter]-voiced_bin_th[i-1])/(counter+1);
                else // if the last frame has no pitch
                    a = 0;

                while(i < n && voiced_bin_th[i] == 0)
                {
                    voiced_bin_th[i] = (int)floor(voiced_bin_th[i-1] + a);
                    i++;
                }
            }
        }
    }
}
*/

#include "lpc.h"

// calculate how many fundamental frequency vowel_marks needed
// return: the number of vowel_marks
int LPC::CalPath(int frame_no, int frame_shift, double *Pitch)
{
	int frame_ptr = 0, n = 0;
	double vowel_marker = 0, T;

	while(vowel_marker < frame_no*frame_shift)
	{
		if (Pitch[frame_ptr] == 0)
		{
			if (frame_ptr == 0)// if the first frame has no pitch
			{
				while(Pitch[frame_ptr] == 0)
					frame_ptr++;
				while(frame_ptr > 0)
				{
					Pitch[frame_ptr-1] = Pitch[frame_ptr];
					frame_ptr--;
				}
			}
			else
				Pitch[frame_ptr] = Pitch[frame_ptr-1];
		}

		vowel_marker += (double)SamplingRate/Pitch[frame_ptr];
		frame_ptr = floor(vowel_marker/frame_shift);
		++n;
	}
	
	return n;
}

// record the fundamental frequency vowel_marks
void LPC::vowel_markPath(int *vowel_mark, int n, int frame_shift, double *Pitch)
{
	int frame_ptr = 0;
	double vowel_marker = 0, T;

	for (int i=0; i<n; ++i)
	{
		vowel_mark[i] = floor(vowel_marker);
		vowel_marker += (double)SamplingRate/Pitch[frame_ptr];
		frame_ptr = floor(vowel_marker/frame_shift);
	}
}

void LPC::SpecEnv(double *coeff, double *SpecEnv_mag, double *SpecEnv_pha)
{
    double A_real, A_imag;
    int i,j;

    for (i=0; i<FFT_SIZE/2+1; ++i)
    {
        A_real = 1;
        A_imag = 0;

        for (j=1; j<=order; j++)
        {
            A_real += coeff[j-1]*cos(i*j*2.0*PI/FFT_SIZE);
            A_imag -= coeff[j-1]*sin(i*j*2.0*PI/FFT_SIZE);
        }

        SpecEnv_mag[i] = 1.0/ABS2(A_real, A_imag);
        SpecEnv_pha[i] = -1*atan2(A_imag, A_real);
    }
}

/* Input : n elements of time domain data
Output: m lpc coefficients, excitation energy */
double LPC::From_Data(double *data,double *lpc,int n,int m)
{
    double *aut= new double[m+1];
    double error;
    int i,j;

    /* autocorrelation, p+1 lag coefficients */

    j=m+1;
    while(j--){
        double d=0;
        for(i=j;i<n;i++)d+=data[i]*data[i-j];
        aut[j]=d;
    }

    /* Generate lpc coefficients from autocorr values */

    error=aut[0];
    if(error==0){
        memset(lpc,0,m*sizeof(double));
        return 0;
    }

    for(i=0;i<m;i++){
        double r=-aut[i+1];

        /* Sum up this iteration's reflection coefficient; note that in
           Vorbis we don't save it.  If anyone wants to recycle this code
           and needs reflection coefficients, save the results of 'r' from
           each iteration. */

        for(j=0;j<i;j++)r-=lpc[j]*aut[i-j];
        r/=error; 

        /* Update LPC coefficients and total error */

        lpc[i]=r;
        for(j=0;j<i/2;j++){
            double tmp=lpc[j];
            lpc[j]+=r*lpc[i-1-j];
            lpc[i-1-j]+=r*tmp;
        }
        if(i%2)lpc[j]+=lpc[j]*r;

        error*=1.0-r*r;
    }

    /* we need the error value to know how big an impulse to hit the
       filter with later */
    delete [] aut;
    return error;
}

void LPC::Predict(double *coeff,double *prime,int m,double *data,long n)
{
    /* in: coeff[0...m-1] LPC coefficients 
       prime[0...m-1] initial values (allocated size of n+m-1)
out: data[0...n-1] data samples */

    long i,j,o,p;
    double y;
    double *work = new double[m+n];

    if(!prime)
        for(i=0;i<m;i++)
            work[i]=0.;
    else
        for(i=0;i<m;i++)
            work[i]=prime[i];

    for(i=0;i<n;i++){
        y=0;
        o=i;
        p=m;
        for(j=0;j<m;j++)
            y-=work[o++]*coeff[--p];
        data[i]=work[o]=y;
    }

    delete [] work;
}

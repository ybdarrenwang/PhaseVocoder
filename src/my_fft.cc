#include "my_fft.h"

void MyFFT::fft_float(unsigned  NumSamples, bool InverseTransform, float *RealIn, float *ImagIn, float *RealOut, float *ImagOut)
{
	unsigned NumBits;    /* Number of bits needed to store indices */
	unsigned i, j, k, n;
	unsigned BlockSize, BlockEnd;

	float angle_numerator = 2.0 * DDC_PI;
	float tr, ti;     /* temp real, temp imaginary */

	if(InverseTransform)
		angle_numerator = -angle_numerator;

	NumBits = NumberOfBitsNeeded ( NumSamples );

	/* Do simultaneous data copy and bit-reversal ordering into outputs... */
	for ( i=0; i < NumSamples; i++ )
	{
		j = ReverseBits ( i, NumBits );
		RealOut[j] = RealIn[i];
		ImagOut[j] = (ImagIn == NULL) ? 0.0 : ImagIn[i];
	}

	/*  Do the FFT itself... */
	BlockEnd = 1;
	for ( BlockSize = 2; BlockSize <= NumSamples; BlockSize <<= 1 )
	{
		float delta_angle = angle_numerator / (float)BlockSize;
		float sm2 = sin ( -2 * delta_angle );
		float sm1 = sin ( -delta_angle );
		float cm2 = cos ( -2 * delta_angle );
		float cm1 = cos ( -delta_angle );
		float w = 2 * cm1;
		float ar[3], ai[3];

		for ( i=0; i < NumSamples; i += BlockSize )
		{
			ar[2] = cm2;
			ar[1] = cm1;

			ai[2] = sm2;
			ai[1] = sm1;

			for ( j=i, n=0; n < BlockEnd; j++, n++ )
			{
				ar[0] = w*ar[1] - ar[2];
				ar[2] = ar[1];
				ar[1] = ar[0];

				ai[0] = w*ai[1] - ai[2];
				ai[2] = ai[1];
				ai[1] = ai[0];

				k = j + BlockEnd;
				tr = ar[0]*RealOut[k] - ai[0]*ImagOut[k];
				ti = ar[0]*ImagOut[k] + ai[0]*RealOut[k];

				RealOut[k] = RealOut[j] - tr;
				ImagOut[k] = ImagOut[j] - ti;

				RealOut[j] += tr;
				ImagOut[j] += ti;
			}
		}
		BlockEnd = BlockSize;
	}
	if(InverseTransform) {
		float denom = (float)NumSamples;
		for(i=0;i<NumSamples;++i)
        {
			RealOut[i] /= denom;
			ImagOut[i] /= denom;
		}
	}
}

unsigned MyFFT::NumberOfBitsNeeded ( unsigned PowerOfTwo )
{
	for (unsigned int i=0; ; ++i)
		if ( PowerOfTwo & (1 << i) )
			return i;
}

unsigned MyFFT::ReverseBits (unsigned index, unsigned NumBits)
{
	unsigned i, rev;
	for (i=rev=0; i < NumBits; ++i)
    {
		rev = (rev << 1) | (index & 1);
		index >>= 1;
	}
	return rev;
}

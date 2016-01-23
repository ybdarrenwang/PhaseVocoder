#ifndef HAMMING_WINDOW_H
#define HAMMING_WINDOW_H

#include "global.h"

class HammingWindow : public Window
{
    public:
        HammingWindow(int l) : Window(l) {
            int i, half_length = floor((double)length/2);
            if (length%2)
                for (i=0; i<=half_length; i++)
                    _window[half_length+i] = _window[half_length-i] = 0.53836 - 0.46164 * cos(2*PI*(half_length-i)/(length-1));
            else
                for (i=0; i<half_length; i++)
                    _window[half_length+i] = _window[half_length-i-1] = 0.53836 - 0.46164 * cos(2*PI*(half_length-i)/length);
        }
};

#endif

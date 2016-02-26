#ifndef __WINDOW_H__
#define __WINDOW_H__

#include "global.h"

class Window
{
    public:
        Window(int l) : length(l) { _window = new double[length]; }
        virtual ~Window() { delete _window; }
        virtual void applyWindow(double* input){
            for (int i=0; i<length; ++i)
                input[i]*=_window[i];
        }

    protected:
        int length;
        double* _window;
};

#endif

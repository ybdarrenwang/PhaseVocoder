#ifndef __WINDOW_H__
#define __WINDOW_H__

#include "global.h"

class Window
{
    public:
        Window(int l) : length(l) {
            _window = new float[length];
        }
        virtual ~Window() {
            delete _window;
        }
        virtual void applyWindow(float* input);

    protected:
        int length;
        float* _window;
};

#endif

#ifndef __WINDOW_H__
#define __WINDOW_H__

#include "global.h"

class Window
{
    public:
        Window(int l) : length(l) {_window = new float[length];}
        virtual ~Window() {delete _window;}
        float* getHamming();

    private:
        int length;
        float* _window;
};

#endif

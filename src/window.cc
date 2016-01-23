#include "window.h"

void Window::applyWindow(float* input) {
    for (int i=0; i<length; ++i)
        input[i]*=_window[i];
}

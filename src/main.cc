#include "time_stretcher.h"
#include "pitch_shifter.h"
#include "overlap_adder.h"
#include "complex.h"
#include "window.h"
#include "frame.h"
#include <iostream>

using namespace std;

int main()
{
    Complex* c = new Complex(1.0, 2.0);
    Frame* f = new Frame(256);
    Window* w = new Window(256);
    w->getHamming();
    delete c;
    delete w;
    delete f;
    cout<<"end"<<endl;
    return 0;
}

#!/bin/sh

# Half speed
../PhaseVocoder.exe -t 2.0 -i HungarianDanceNo5.wav -o half_speed.wav

# Double speed
../PhaseVocoder.exe -t 0.5 -i HungarianDanceNo5.wav -o double_speed.wav

# 1 Whole step higher
../PhaseVocoder.exe -p 1.125 -i HungarianDanceNo5.wav -o higher.wav

# 1 Whole step lower
../PhaseVocoder.exe -p 0.889 -i HungarianDanceNo5.wav -o lower.wav

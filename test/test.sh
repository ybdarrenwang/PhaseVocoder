#!/bin/sh

# Half speed
../PhaseVocoder.exe -t 2.0 -i HungarianDanceNo5.wav -o half_speed.wav

# Double speed
../PhaseVocoder.exe -t 0.5 -i HungarianDanceNo5.wav -o double_speed.wav

# Half speed, phase locked
../PhaseVocoder.exe -t 2.0 -i HungarianDanceNo5.wav -o half_speed_pl.wav --phaselock

# Double speed, phase locked
../PhaseVocoder.exe -t 0.5 -i HungarianDanceNo5.wav -o double_speed_pl.wav --phaselock

# 1 Whole step higher
../PhaseVocoder.exe -p 1.125 -i HungarianDanceNo5.wav -o higher_pitch.wav

# 1 Whole step lower
../PhaseVocoder.exe -p 0.889 -i HungarianDanceNo5.wav -o lower_pitch.wav

../PhaseVocoder.exe -p 0.889 -t 1.5 -i HungarianDanceNo5.wav -o tmp1.wav

../PhaseVocoder.exe -p 1.125 -t 0.7 -i HungarianDanceNo5.wav -o tmp2.wav --phaselock

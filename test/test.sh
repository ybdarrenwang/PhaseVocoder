#!/bin/sh

mkdir -p output/

# Half speed
../bin/PhaseVocoder -t 2.0 -i HungarianDanceNo5.wav -o output/half_speed.wav
# Double speed
../bin/PhaseVocoder -t 0.5 -i HungarianDanceNo5.wav -o output/double_speed.wav
# Half speed, phase locked
../bin/PhaseVocoder -t 2.0 -i HungarianDanceNo5.wav -o output/half_speed_pl.wav --phaselock
# Double speed, phase locked
../bin/PhaseVocoder -t 0.5 -i HungarianDanceNo5.wav -o output/double_speed_pl.wav --phaselock
# Half speed, spec interpolate
../bin/PhaseVocoder -t 2.0 -i HungarianDanceNo5.wav -o output/half_speed_fd.wav --specInterpolate
# Double speed, spec interpolate
../bin/PhaseVocoder -t 0.5 -i HungarianDanceNo5.wav -o output/double_speed_fd.wav --specInterpolate
# Half speed, spec interpolate
../bin/PhaseVocoder -t 2.0 -i HungarianDanceNo5.wav -o output/half_speed_fd_pl.wav --specInterpolate --phaselock
# Double speed, spec interpolate
../bin/PhaseVocoder -t 0.5 -i HungarianDanceNo5.wav -o output/double_speed_fd_pl.wav --specInterpolate --phaselock
# 1 Whole step higher
../bin/PhaseVocoder -p 1.125 -i HungarianDanceNo5.wav -o output/higher_pitch.wav
# 1 Whole step lower
../bin/PhaseVocoder -p 0.889 -i HungarianDanceNo5.wav -o output/lower_pitch.wav

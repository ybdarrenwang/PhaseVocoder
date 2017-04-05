Phase Vocoder
=============

This is a phase vocoder written in C++. The main purpose of a phase vocoder is to change the speed and pitch of a given voice recording file.

Current version supports the following time stretching algorithms
- Default: number of frames remain the same, while the output synthesis-rate/frame-shift/hop-size is modified correspondingly.
- --phaseLock: perform phase locking as described in [Laroche99]
- --fdTimeStretch: perform interpolate/extrapolate on spectrogram, thus the number of frames changed; while the output synthesis-rate/frame-shift/hop-size remain the same as analysis.
- Both --phaseLocking and --fdTimeStretch

Heuristically, use --phaseLock for human speech and --fdTimeStretch for music yield the best result.

Current version supports the following pitch shifting algorithms
- Default: perform time stretching for compensation, then modify output synthesis-rate/frame-shift/hop-size so the total length remains but pitch changed (due to up/down sampling)
- --fdPitchShift: directly modify spectrogram; both the number of frames and output synthesis-rate/frame-shift/hop-size remain the same as analysis.
Default approach is suggested. --fdPitchShift hasn't been optimized yet.

Platforms
---------
Linux, OS X

Requirements
------------
g++ 4.2 or higher

Installation
------------
`make` to create PhaseVocoder.exe.

Test
----
After PhaseVocoder.exe is complied, go to test/ directory and run `sh test.sh`.

Reference
---------
- Moulines, Eric, and Jean Laroche. "Non-parametric techniques for pitch-scale and time-scale modification of speech." Speech communication 16.2 (1995): 175-205.
- Laroche, Jean, and Mark Dolson. "Improved phase vocoder time-scale modification of audio." Speech and Audio Processing, IEEE Transactions on 7.3 (1999): 323-332.

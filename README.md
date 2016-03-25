Phase Vocoder
=============

This is a phase vocoder written in C++. The main purpose of a phase vocoder is to change the speed and pitch of a given voice recording file.

Current version supports the following time stretching algorithms
- Default: number of frames remain the same, while the output synthesis-rate/frame-shift/hop-size is modified correspondingly.
- --phaseLock: perform phase locking as described in [Laroche1999]
- --specInterpolate: perform interpolate/extrapolate on spectrogram, thus the number of frames changed; while the output synthesis-rate/frame-shift/hop-size remain the same as analysis.
- Both --phaseLocking and --specInterpolate

Heuristically, use --phaseLock for human speech and --specInterpolate for music yield the best result.

As to pitch shifting, note the implementation of current version is different from usual.
* Popular approach: perform time stretching for compensation, then modify output synthesis-rate/frame-shift/hop-size so the total length remains but pitch changed (due to up/down sampling)
* My version: directly modify spectrogram; both the number of frames and output synthesis-rate/frame-shift/hop-size remain the same as analysis.
Shall implement the common approach in the future.

Also note this version does not include source-filter decomposition, nor other advanced pitch shifting algorithms.

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

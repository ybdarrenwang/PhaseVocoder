CXX=g++
CXXFLAGS=-O3 -Wall

SRC=src
BUILD=build

SOURCES=main.cc \
    vocoder_functions.cc \
    time_stretcher.cc \
    time_stretcher_pl.cc \
    time_stretcher_fd.cc \
    time_stretcher_fd_pl.cc \
    my_fft.cc \
    frame.cc \
    wav_io.cc \
    pitch_shifter.cc \
    phasevocoder.cc

OBJ:=$(addprefix $(BUILD)/, $(addsuffix .o,$(basename $(SOURCES))))

.PHONY: dir clean

all: dir program

debug: CXXFLAGS += -DDEBUG -g
debug: all

dir:
	mkdir -p $(BUILD)

program: $(OBJ)
	$(CXX) $(CXXFLAGS) -o PhaseVocoder.exe $^

$(BUILD)/%.o: $(SRC)/%.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(BUILD)

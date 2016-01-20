GPP = g++
GPPFlag = -O2 -Wall

SRC = src
BUILD = build

CDIR = cd
MKDIR = mkdir -p
RM = rm -f

objects = $(BUILD)/main.o $(BUILD)/lpc.o $(BUILD)/VocoderFunctions.o $(BUILD)/SourceFilter.o

.PHONY : dir

all: dir program

dir:
	$(MKDIR) $(BUILD)

program: $(objects)
	$(GPP) $(GPPFlag) -o PhaseVocoder.exe $^

$(BUILD)/%.o: $(SRC)/%.cpp
	$(GPP) -c -o $@ $<

.PHONY : clean
clean:
	$(RM) $(BUILD)/*

GPP = g++
GPPFlag = -O2 -Wall

SRC = src
BUILD = build
DOC = doc
DOXYGEN = doxygen
DOXYGENCFG = $(DOC)/doxygen.cfg

CDIR = cd
MKDIR = mkdir -p
RM = rm -f

objects = $(BUILD)/main.o $(BUILD)/SocialNet.o $(BUILD)/User.o $(BUILD)/Message.o $(BUILD)/Photo.o $(BUILD)/FacebookLikeDecorator.o $(BUILD)/FacebookReplyDecorator.o $(BUILD)/TwitterReplyDecorator.o $(BUILD)/TwitterTweetDecorator.o $(BUILD)/json.o $(BUILD)/Facebook.o $(BUILD)/Twitter.o $(BUILD)/AutoCompleteBasicQ.o $(BUILD)/AutoCompleteAdvQ.o

.PHONY : dir

all: dir program

dir:
	$(MKDIR) $(BUILD)

program: $(objects)
	$(GPP) $(GPPFlag) -o SocialNetworkSearchEngine.exe $^

$(BUILD)/%.o: $(SRC)/%.cpp
	$(GPP) -c -o $@ $<

.PHONY : doc
doc:
	$(CDIR) $(DOC);
	$(DOXYGEN) $(DOXYGENCFG)

.PHONY : clean
clean:
	$(RM) $(BUILD)/*

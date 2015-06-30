ROSS:=/usr/local/
SRC_BASE:=/Users/Mark/Desktop/RPI-Research
CODESBASE:=$(SRC_BASE)/codes-base/INSTALL
CODESNET:=$(SRC_BASE)/codes-net/INSTALL
ASPEN:=$(SRC_BASE)/aspen-sim


ifndef CODESBASE
$(error CODESBASE is undefined, see README.txt)
endif
ifndef CODESNET
$(error CODESNET is undefined, see README.txt)
endif
ifndef ROSS
$(error ROSS is undefined, see README.txt)
endif

# ross conflates CFLAGS with CPPFLAGS, so use just this one
override CPPFLAGS += $(shell $(ROSS)/bin/ross-config --cflags) -I$(CODESBASE)/include -I$(CODESNET)/include -I$(ASPEN)/aspen
CC = $(shell $(ROSS)/bin/ross-config --cc)
LDFLAGS = $(shell $(ROSS)/bin/ross-config --ldflags) -L$(CODESBASE)/lib -L$(CODESNET)/lib -L$(ASPEN)/lib
LDLIBS = $(shell $(ROSS)/bin/ross-config --libs) -lcodes-net -lcodes-base -laspen -lc++

AspenNet:
	$(CC) -g $(CPPFLAGS) -v -O2 test.c AspenNet_AspenUtils.cpp $(LDFLAGS) $(LDLIBS) -o AspenNet

clean:   
	rm -f AspenNet *.o

SRC_BASE:=/home/blancm3/RPI-Research
ROSS:=$(SRC_BASE)/ROSS/INSTALL
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
LDLIBS = -laspen -lcodes-net -lcodes-base $(shell $(ROSS)/bin/ross-config --libs) -lcodes-net#-lc++

# Actual compilation directions:

AspenNet_AspenUtils.o: AspenNet_AspenUtils.cpp
	mpic++ -g $(CPPFLAGS) -c -O2 AspenNet_AspenUtils.cpp -o AspenNet_AspenUtils.o

AspenNet.o: AspenNet.c
	$(CC) -g $(CPPFLAGS) -c -O2 AspenNet.c -o AspenNet.o

test.o: test.c
	$(CC) -g $(CPPFLAGS) -c -O2 test.c -o test.o

AspenNet: AspenNet.o AspenNet_AspenUtils.o
	mpic++ -g $(CPPFLAGS)  AspenNet.o AspenNet_AspenUtils.o  $(LDFLAGS) $(LDLIBS) -o AspenNet

testfile: test.o AspenNet_AspenUtils.o
	mpic++ -g $(CPPFLAGS) -O2 test.o AspenNet_AspenUtils.o $(LDFLAGS) $(LDLIBS) -o test

clean:   
	rm -f *.o
	rm -f AspenNet
	rm -f test

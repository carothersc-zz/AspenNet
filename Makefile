SRC_BASE:=/home/blancm3/RPI-Research
ROSS:=$(SRC_BASE)/ROSS/INSTALL
CODESBASE:=$(SRC_BASE)/codes/INSTALL
#CODESNET:=$(SRC_BASE)/codes-net/INSTALL
ASPEN:=$(SRC_BASE)/aspen-2016

#if test -d $CODESBASE ;
# then $(error CODESBASE is undefined, see README.txt);
#fi
#ifndef ROSS
#$(error ROSS is undefined, see README.txt)
#endif

# ross conflates CFLAGS with CPPFLAGS, so use just this one
override CPPFLAGS += $(shell $(ROSS)/bin/ross-config --cflags) -I$(CODESBASE)/include -I$(ASPEN)/aspen -I$(ASPEN)/c
#override CPPFLAGS += -O3
CC = $(shell $(ROSS)/bin/ross-config --cc)
LDFLAGS = $(shell $(ROSS)/bin/ross-config --ldflags) -L$(CODESBASE)/lib -L$(ASPEN)/c -L$(ASPEN)/lib 
LDLIBS = -laspenc -laspen -lcodes $(shell $(ROSS)/bin/ross-config --libs) -rdynamic


# Actual compilation directions:

AspenNet.o: AspenNet.c
	$(CC) $(CPPFLAGS) -g -c AspenNet.c -o AspenNet.o

test: test.o
	mpic++  -o test test.o $(LDFLAGS) $(LDLIBS)
	
test.o:
	$(CC) -fPIC -g $(CPPFLAGS) -c test.c -o test.o

AspenNet: AspenNet.o
	mpic++ $(CPPFLAGS) -g  AspenNet.o AspenNet_AspenUtils.o  $(LDFLAGS) $(LDLIBS) -o AspenNet

testfile: test.o
	$(CC) $(CPPFLAGS) -g test.o $(LDFLAGS) $(LDLIBS) -o test

clean:   
	rm -f *.o
	rm -f AspenNet
	rm -f test

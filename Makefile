#
# Makefile for AO transition densities
# ctc 8-5-14
#

# Compiler Choice and target executable name
CC   = g++
EXEC = aoproject.out

# Options for linux 
LDFLAGS = -L/usr/lib -L/usr/local/lib -lm -lgsl -lgslcblas
CPPFLAGS = -fopenmp -I/usr/include -O3 -Wno-unused-variable -Wno-deprecated -fpermissive
#CPPFLAGS = -I/usr/include -g3 -ggdb -Wall -Wno-unused-variable -Wno-deprecated
#CPPFLAGS = -O2 -Wno-deprecated
#CPPFLAGS = -O3 -march=pentium4 -malign-double -Wno-deprecated
#CPPFLAGS = -g3 -ggdb -Wall -Wno-deprecated

SHELL = /bin/sh

HEADERS1 = main.h atom.h molecule.h util.h
HEADERS  = $(HEADERS1)

OBJS = main.o atom.o molecule.o

#fix : fixPulses.o wpiGridParameters.o wpiInterpolation.o
#	$(CC) fixPulses.o wpiGridParameters.o wpiInterpolation.o -o fixer $(LDFLAGS) $(CLIBS) $(CPPFLAGS)

new : $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) $(LDFLAGS) $(CLIBS) $(CPPFLAGS)


$(OBJS)	: $(HEADERS)

$(HEADERS) : 

tarfile :
	tar cf aodens.tar Makefile *.h *.cc

.PHONY: .clean
clean	:
	rm -f $(OBJS) $(EXEC)

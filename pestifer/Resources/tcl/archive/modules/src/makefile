# BONDSTRUCT -- a simple bond topology manager for Tcl
#
# Cameron F Abrams 2016-2018
# Drexel University, Philadelphia, Pennsylvania
#
# cfa22@drexel.edu
#
# makefile for bondstruct.so TcL module

CC=gcc

LIBS += -lm

CFLAGS += -O3 -I/usr/include/tcl

USER_FLAGS =

OBJ += bondstruct.o
OBJ += linkcell.o

all: bondstruct.so 

bondstruct.so:  $(OBJ) bondstruct_wrap.o
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LIBS)
#	cp $@ ../lib/

%.o: %.c
	$(CC) -c -fpic $(CFLAGS) $< $(LIBS)

%_wrap.c: %.i %.c %.h
	swig -tcl8 $<

bondstruct.o : bondstruct.h

linkcell.o : linkcell.h

.PHONY: clean force

clean:
	rm -f *.o *.so *_wrap.c


# Makefile for Rote sequences generator.

CC = g++
CFLAGS = -Wall -O3 -ffast-math -fast
LDFLAGS = -lcgal -lgmp
INCLUDES = 

objects = rote_generation.o 

all: generator

generator:
	$(CC) $(CFLAGS) $(INCLUDES) -o generator rote_generation.cpp $(LDFLAGS) -frounding-math

doc:
	doxygen -u
	doxygen 

pdf:
	make -C doc/latex/

clean:
	rm -f generator *.o Doxyfile.bak
	rm -rf doc/

.PHONY: clean

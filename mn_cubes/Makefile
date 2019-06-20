# Makefile for m,n-cubes counting. Compiles Bentley-Ottmann first.

CC = g++
CFLAGS = -O2 -Wall
LDFLAGS = -lgmp -lgmpxx
INCLUDES =

objects = bentley_ottmann/geometry.o bentley_ottmann/bentley_ottmann.o mncubes.o 

all: mncubes

mncubes: $(objects)
	$(CC) $(CFLAGS) $(INCLUDES) -o mncubes bentley_ottmann/*.o mncubes.o $(LDFLAGS)

$(objects): %.o: %.cpp
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

doc:
	make -C bentley_ottmann/ doc

doc_tex:
	make -C bentley_ottmann/ tex

clean:
	rm -f bentley_ottmann/*.o *.o mncubes
	make -C bentley_ottmann/ clean

.PHONY: clean

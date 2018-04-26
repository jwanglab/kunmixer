CC=gcc
CFLAGS=-lhts -std=c99

OBJECTS = fastdiff fasttype

all: $(OBJECTS)

fastdiff: fastdiff.c
	$(CC) $(CFLAGS) fastdiff.c -o fastdiff

fasttype: fasttype.c
	$(CC) $(CFLAGS) bed.c fasttype.c -o fasttype

.PHONY: clean
clean:
	-rm $(OBJECTS)

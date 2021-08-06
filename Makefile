CC=gcc
CFLAGS=-std=c99

OBJECTS = fastdiff fasttype ktype

all: $(OBJECTS)

fastdiff: fastdiff.c
	$(CC) $(CFLAGS) fastdiff.c -o fastdiff -lhts

fasttype: fasttype.c
	$(CC) $(CFLAGS) bed.c fasttype.c -o fasttype -lhts

ktype: ktype.c
	$(CC) $(CFLAGS) bed.c ktype.c -o ktype -lz

.PHONY: clean
clean:
	-rm $(OBJECTS)

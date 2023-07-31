CC=gcc
CFLAGS=-std=c99

OBJECTS = alntype ktype fastdiff

all: $(OBJECTS)

fastdiff: fastdiff.c
	$(CC) $(CFLAGS) fastdiff.c -o fastdiff -lhts

alntype: alntype.c
	$(CC) $(CFLAGS) bed.c alntype.c -o alntype -lhts

ktype: ktype.c
	$(CC) $(CFLAGS) bed.c ktype.c -o ktype -lz -lbz2

.PHONY: clean
clean:
	-rm $(OBJECTS)

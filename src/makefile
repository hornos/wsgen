CC=gcc
CFLAGS=-Wall
LDFLAGS=-lm
SOURCES=wsgen.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=wsgen

.PHONY = all clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)

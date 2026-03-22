CC = gcc
LD = gcc
CFLAGS = -Ofast -Wall -fopenmp -march=native -lm
LDFLAGS = -fopenmp -lm
RM = /bin/rm -f
OBJS = main.o quicksort.o
EXECUTABLE = main

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
		$(LD) $(OBJS) $(LDFLAGS) -o $(EXECUTABLE)

quicksort.o: quicksort.h quicksort.c
		$(CC) $(CFLAGS) -c quicksort.c

main.o: main.c quicksort.h
		$(CC) $(CFLAGS) -c main.c

clean:
		$(RM) $(EXECUTABLE) $(OBJS)

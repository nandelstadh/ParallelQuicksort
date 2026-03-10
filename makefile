CC = gcc-15
LD = gcc-15
CFLAGS = -O3 -Wall -Werror -fopenmp -fsanitize=address
LDFLAGS = -fopenmp
RM = /bin/rm -f
OBJS = main.o quicksort.o
EXECUTABLE = main

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE)

quicksort.o: quicksort.h quicksort.c
	$(CC) $(CFLAGS) -c quicksort.c

main.o: main.c quicksort.h
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)

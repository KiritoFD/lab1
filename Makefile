CC = gcc
CFLAGS = -Wall -Wextra -g

all: repeat_finder

repeat_finder: main.c dag_algorithms.c
	$(CC) $(CFLAGS) -o repeat_finder main.c dag_algorithms.c

clean:
	rm -f repeat_finder
	rm -f *.o

.PHONY: all clean

CC = gcc
CFLAGS = -Wall -Wextra -g

all: dag_demo

dag_demo: main.c dag_algorithms.c dag_algorithms.h
	$(CC) $(CFLAGS) -o dag_demo main.c dag_algorithms.c

clean:
	rm -f dag_demo
	rm -f *.o

.PHONY: all clean

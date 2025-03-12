CC = gcc
CFLAGS = -Wall -Wextra -O2

all: dna_repeat_finder

dna_repeat_finder: dna_repeat_finder.c
	$(CC) $(CFLAGS) -o dna_repeat_finder dna_repeat_finder.c

clean:
	rm -f dna_repeat_finder

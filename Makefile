CC ?= gcc
CFLAGS ?= -Wall -g -O3 -pthread

cosmotype: src/cosmotype.c
		$(CC) $(CFLAGS) src/cosmotype.c -lhts -o cosmotype 

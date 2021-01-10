CC ?= gcc
CFLAGS ?= -Wall -g -O3 -pthread

somaticfreq: src/somaticfreq.c
		$(CC) $(CFLAGS) src/somaticfreq.c -lhts -o somaticfreq


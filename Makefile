CC ?= gcc
CFLAGS ?= -Wall -g -O3 -pthread

basecounts: src/basecounts.c
		$(CC) $(CFLAGS) -lhts -o basecounts src/basecounts.c

test: basecounts
		./basecounts test/loci.bed /Volumes/datadrive/hg19/hg19.fa test/lmo2_loci.bam > test/res.tsv
		cat test/res.tsv

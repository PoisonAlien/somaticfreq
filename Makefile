CC ?= gcc
CFLAGS ?= -Wall -g -O3 -pthread

cosmotype: src/cosmotype.c
		$(CC) $(CFLAGS) -lhts -o cosmotype src/cosmotype.c

test: basecounts
		./cosmotype -o myreport test/loci.bed test/lmo2_loci.bam > test/res.tsv
		cat test/res.tsv

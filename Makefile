CC ?= gcc
CFLAGS ?= -Wall -g -O3 -pthread

cosmotype: src/cosmotype.c
		$(CC) $(CFLAGS) src/cosmotype.c -lhts -o cosmotype 

test: basecounts
		./cosmotype -o myreport test/loci.bed test/lmo2_loci.bam > test/res.tsv
		cat test/res.tsv

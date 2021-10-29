## somaticfreq: knowledge-based genotyping of targetted somatic variants from the tumor BAM file

`somaticfreq` allows rapid (kind of) genotyping of known somatic variants from the tumor BAM files. This facilitates to get a quick overlook of known somatic hot-spots in a matter of minutes, without spending hours on variant calling and annotation. In simple words, it fetches nucleotide frequencies of known somatic hotspots and prioritizes them based on allele frequency. Output includes a browsable/sharable HTML report of candidate variants (an [example](https://poisonalien.github.io/STAR_wrapper_script/)). See below for data sources. 

### Installation: 

Clone the repository and run `make`. Requires [htslib](https://github.com/samtools/htslib). 

```
git clone https://github.com/PoisonAlien/somaticfreq
cd somaticfreq
make
```

### Usage:

`somaticfreq` takes a tsv file of known somatic variants and a BAM file as input.

```
$ somaticreq
somaticfreq: A tool to extract nucleotide counts/variant allele frequencies
             of targeted (somatic) variants from the BAM file
USAGE:
    somaticfreq [OPTIONS] <loci> <bam>
    e.g; somaticfreq cancerhotspots_v2_GRCh37.tsv Tumor.bam
OPTIONS:
    -f  Indexed fasta file. If provided, extracts and adds reference base to the output tsv
    -q  Mapping quality threshold. Default 10 [Read filter]
    -F  Exclude reads with FLAGS >= INT. Default 1024 (i.e, read is PCR or optical duplicate) [Read filter]
    -v  VAF threshold. Default 0.05 [Variant filter]
    -d  Depth of coverage threshold. Default 30 [Variant filter]
    -t  Min. number of reads supporting tumor allele. Default 8 [Variant filter]
    -o  Output file basename. Default parses from basename of BAM file
ARGS:
    <loci>   variant file
    <bam>    Indexed BAM file
OUPUT:
    <output>.html   A browsable html report of variants post filtering
    <output>.tsv    TSV file with nucletide counts and VAF for all variants
```

Example:
```
$ somaticfreq data/cancerhotspots_v2_GRCh37.tsv Tumor.bam 
```

### Output

Output includes a browsable HTML file with variants passing the VAF/read depth filters and, a TSV file including nucleotide counts of all variants analyzed.

* Example report from a [metastatic breast cancer](https://poisonalien.github.io/STAR_wrapper_script/) sample

### Data sources: 

1. Known [cancerhotspots](https://www.cancerhotspots.org/#/download) are [included](https://github.com/PoisonAlien/somaticfreq/tree/main/data) as a part of this repository for both GRCh37 and GRCh38 assemblies (3180 variants). This should be sufficient and cover most of the known driver genes/events.
2. COSMIC database provides a comprehensive list of PATHOGENIC, non-SNP, and confirmed SOMATIC variants (approx. 90K variants). You can download the tsv file from [here](https://cancer.sanger.ac.uk/cosmic/download) and use this [R script](https://gist.github.com/PoisonAlien/fa4199e34a089a873820fd46eba028df) to filter and convert it to the input format for `somaticfreq`. (Requires free registration and free for academic usage).
3. You can also create your own custom database of variants. Look at one of the [example files](https://github.com/PoisonAlien/somaticfreq/tree/main/data) for an example. Positions are in the 1-based coordinate system.

### References

1. Chang MT, Asthana S, Gao SP, et al. Identifying recurrent mutations in cancer reveals widespread lineage diversity and mutational specificity. Nat Biotechnol. 2016;34(2):155-163. doi:10.1038/nbt.3391

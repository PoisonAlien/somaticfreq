## somaticfreq: knowledge based genotyping of targetted somatic variants from tumor BAM file

`somaticfreq` allows rapid genotyping (kind of) known somatic variants from BAM files. This facilitates to get a quick overlook of known somatic hot-spots in a matter of minutes, without spending hours on variant calling and annotation. In simple words, it fetches nucleotide frequencies of known somatic hotspots and prioritizes them based on allele frequency. Output includes a browsable/sharable HTML report of candidate variants. See below for data sources. 

### Installation: 

Clone the repository and run `make`. Requires [htslib](https://github.com/samtools/htslib). 

```
git clone https://github.com/PoisonAlien/somaticfreq
cd somaticfreq
make
```

If you have [htslib](https://github.com/samtools/htslib) installed in non standard directory (for ex. from Conda installtion) add the `lib` and `include` directories to make command.

```
make CFLAGS="-I/home/user/miniconda2/include -L/home/user/miniconda2/lib"
```

### Usage:

`somaticfreq` takes a tsv file of known somatic variants and a BAM file as input.

```
$ somaticreq
somaticfreq: A tool to extract nucleotide counts/variant allele frequencies
             of targeted (somatic) variants from BAM file
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
$ somaticfreq cancerhotspots_v2_GRCh37.tsv Tumor.bam 
```

### Output

Output includes a browsable HTML file with variants passing the VAF/read depth filters and, a TSV file including nucleotide counts of all variants analyzed.

    Example report from a [metastatic brest cancer](https://poisonalien.github.io/STAR_wrapper_script/) sample

### Data sources: 

1. Known [cancerhotspots](https://www.cancerhotspots.org/#/download) are included as a part of this repository for both GRCh37 and GRCh38 assemblies (3180 variants).
2. COSMIC database provides a comprehensive list of PATHOGENIC, non-SNP, and confirmed SOMATIC variants (approx. 90K variants). You can download the tsv file from [here](https://cancer.sanger.ac.uk/cosmic/download) and use this [R script](https://gist.github.com/PoisonAlien/fa4199e34a089a873820fd46eba028df) to filter and convert it to the input format for `somaticfreq` (Requires free registration and free for academic usage).
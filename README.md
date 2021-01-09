## somaticfreq: A tool to extract nucleotide counts/variant allele frequencies of targeted (somatic) variants from BAM file

`somaticfreq` allows rapid genotyping (sort of) known somatic variants from BAM files. This facilitates to get a quick overlook of known somatic hot-spots in a matter of minutes, without spending hours on variant calling and annotation. Output includes a browsable/sharable HTML report of candidate variants. 

Installation: 

Clone the repository and run `make`

```
git clone https://github.com/PoisonAlien/somaticfreq
cd somaticfreq
make
```

Usage:

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

Output includes a browsable HTML file with variants passing the VAF/read depth filters and TSV file including nucleotide counts of all variants analyzed.

Data sources: 

1. Known [cancerhotspots](https://www.cancerhotspots.org/#/download) are included as a part of this repository for both GRCh37 and GRCh38 assemblies (3180 variants).
2. COSMIC database provides a comprehensive list of PATHOGENIC, non-SNP, and confirmed SOMATIC variants (approx. 90K variants). You can download the tsv file from [here](https://cancer.sanger.ac.uk/cosmic/download) and use this [R script](https://gist.github.com/PoisonAlien/fa4199e34a089a873820fd46eba028df) to filter and convert it to the input format for `somaticfreq`  
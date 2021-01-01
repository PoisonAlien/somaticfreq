#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

void extbases(char *region, int32_t target_pos, char *bam, char *norm);
void is_bam(char *fname);

void usage(){
    fprintf(stderr, "Program: basecounts\n");
    fprintf(stderr, "A tool to generate nucleotide counts and VAF table from user specific loci or gentotypes\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "    gtftools basecounts <loci> <fasta> <bam>...\n");
    fprintf(stderr, "ARGS:\n");
    fprintf(stderr, "    <loci>   A 2 (or more) column tsv file with chr\\tpos\n");
    fprintf(stderr, "    <fasta>  Indexed fasta file. Extracts and adds reference base to the ouput\n");
    fprintf(stderr, "    <bam>... One or more BAM files\n");
}

int main(int argc, char *argv[]){

    if(argc < 3){
        usage();
        return 0;
    }
    
    //Check if they are BAM files
    for(int b = 3; b < argc; b++){
        is_bam(argv[b]);
    }

    FILE *bed_fp;
    bed_fp = fopen(argv[1], "r");
    char buff[1000];
    char *fract = "true";

    faidx_t *fa = fai_load(argv[2]);

    while(fgets(buff,1000,bed_fp) != NULL){
        //Remove trailing new line chars
        int len = strlen(buff);
        if(buff[len-1] == '\n' ){
            buff[len-1] = 0;
        }
            
        char *chrom = strtok(buff,"\t");
        char *start = strtok(NULL,"\t");
        //char *end = strtok(NULL,"\t");
        char loci[250] = "";
        
        int32_t target_pos = 0;

        strcat(loci, chrom); strcat(loci, ":"); strcat(loci, start); strcat(loci, "-"); strcat(loci, start);
        //printf("%s\n", loci);
        int templen = 100;
        char *seq = fai_fetch(fa, loci, &templen);
        
        target_pos = atoi(start);

        printf("%s:%s\t%s", chrom, start, seq);
        for(int b = 3; b < argc; b++){
            char *bamfile = argv[b];
            extbases(loci, target_pos, bamfile, fract);
        }
        printf("\n");
        free(seq);
    } 
    fclose(bed_fp);

    return 0;
}

void is_bam(char *fname){
    int l = strlen(fname);
    if (l>=4 && strcasecmp(fname+l-4, ".bam") != 0){
        fprintf(stderr, "%s not a BAM file!\n", fname);
        exit(0);
    }
}

void extbases(char *region, int32_t target_pos, char *bam, char *norm){
    
    samFile *fp_in = hts_open(bam,"r"); //open bam file
    hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header

    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_itr_t *samitr = sam_itr_querys(fp_idx, bamHdr, region);    

    int32_t tot_reads = 0;
    float nt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    while(sam_itr_next(fp_in, samitr, aln) > 0){
            
            int32_t pos = aln->core.pos ; //left most position of alignment in zero based coordianate (+1)
            //char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
            uint32_t len = aln->core.l_qseq; //length of the read.
            uint32_t* cig = bam_get_cigar(aln);
                
            uint8_t *q = bam_get_seq(aln); //quality string
            uint32_t q2 = aln->core.qual ; //mapping quality


            if(q2 <= 10){
                break;
            }

            char *qseq = (char *)malloc(len);
            
            for(int i=0; i< len ; i++){
                qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
            }
            int32_t pos_onread = 0;

            for( int k=0;k< aln->core.n_cigar ;++k){
                int cop =cig[k] & BAM_CIGAR_MASK; // CIGAR string
                int cl = cig[k] >> BAM_CIGAR_SHIFT; // CIGAR length

                if(BAM_CIGAR_STR[cop] == 'M'){
                    pos_onread = pos_onread + cl;
                    pos = pos + cl;
                }else if(BAM_CIGAR_STR[cop] == 'S'){
                    pos_onread = pos_onread + cl;
                }else if(BAM_CIGAR_STR[cop] == 'I'){
                    pos_onread = pos_onread + cl;
                }
                

                if(pos > target_pos){
                    //printf("%c\n", BAM_CIGAR_STR[cop]);
                    if(BAM_CIGAR_STR[cop] == 'M'){
                        pos_onread = pos_onread - (pos - target_pos);
                        //printf("%d\t%lld\t%d\t%c\n", pos, aln->core.pos, tar_pos_onread, qseq[pos - target_pos]);
                        if(qseq[pos_onread] == 'A'){
                            nt[0] = nt[0] + 1;
                        }else if(qseq[pos_onread] == 'T'){
                            nt[1] = nt[1] + 1;
                        }else if(qseq[pos_onread] == 'G'){
                            nt[2] = nt[2] + 1;
                        }else if(qseq[pos_onread] == 'C'){
                            nt[3] = nt[3] + 1;
                        }
                        break;
                    }
                }else if(pos == target_pos){
                    if(BAM_CIGAR_STR[cop] == 'I'){
                        nt[4] = nt[4] + 1;
                        // insertion sequence
                        // for(int i = 0; i < cl; i++){
                        //     printf("%c", qseq[pos_onread+i]);
                        // }
                        break;
                    }
                }else if(pos == target_pos){
                    if(BAM_CIGAR_STR[cop] == 'D'){
                        nt[5] = nt[5] + 1;
                        // deletion sequence
                        // for(int i = 0; i < cl; i++){
                        //     printf("%c", qseq[pos_onread+i]);
                        // }
                        break;
                    }
                }                
            }
        }    
        if(strcmp(norm, "true")){
            printf("\t%.2f|%.2f|%.2f|%.2f|%.2f|%.2f", nt[0]/tot_reads, nt[1]/tot_reads, nt[2]/tot_reads, nt[3]/tot_reads, nt[4]/tot_reads, nt[5]/tot_reads);               
        }else{
            printf("\t%.f|%.f|%.f|%.f|%.f|%.f", nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);               
        }
        
        hts_itr_destroy(samitr);
        bam_destroy1(aln);
        sam_hdr_destroy(bamHdr);
        sam_close(fp_in);
}
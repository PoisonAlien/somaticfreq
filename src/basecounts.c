#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

void extbases(char *region, int32_t target_pos, char *bam, char *norm, FILE *html);
void is_bam(char *fname);
void printhead(FILE *fn);
void printrow(FILE *fn,   char *loci, float counts[]);
void printtail(FILE *fn);

void usage(){
    fprintf(stderr, "Program: basecounts\n");
    fprintf(stderr, "A tool to generate nucleotide counts and VAF table from user specific loci or gentotypes\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "    basecounts <loci> <fasta> <bam>...\n");
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
    hts_verbose = 0; //suppresses htslib warnings
    
    //Check if they are BAM files
    for(int b = 3; b < argc; b++){
        is_bam(argv[b]);
    }

    FILE *bed_fp;
    bed_fp = fopen(argv[1], "r");
    char buff[1000];
    char *fract = "true";

    faidx_t *fa = fai_load(argv[2]);

    FILE *html_fp;
    html_fp = fopen("report.html", "w" );
    printhead(html_fp);

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
        
        target_pos = atoi(start)-1;

        printf("%s:%s\t%s", chrom, start, seq);
        for(int b = 3; b < argc; b++){
            char *bamfile = argv[b];
            extbases(loci, target_pos, bamfile, fract, html_fp);
        }
        printf("\n");
        free(seq);
    } 
    fclose(bed_fp);
    printtail(html_fp);

    fclose(html_fp);

    return 0;
}

void is_bam(char *fname){
    int l = strlen(fname);
    if (l>=4 && strcasecmp(fname+l-4, ".bam") != 0){
        fprintf(stderr, "%s not a BAM file!\n", fname);
        exit(0);
    }
}

void extbases(char *region, int32_t target_pos, char *bam, char *norm, FILE *html){
    
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
                }else if(BAM_CIGAR_STR[cop] == 'D'){
                    pos_onread = pos_onread + cl;
                }
                //printf("%c\t", BAM_CIGAR_STR[cop]);
                

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
                    }else if(BAM_CIGAR_STR[cop] == 'D'){
                        nt[5] = nt[5] + 1;
                        // deletion sequence
                        // for(int i = 0; i < cl; i++){
                        //     printf("%c", qseq[pos_onread+i]);
                        // }
                        break;
                    }
                }               
            }
            //printf("\n");
        }    
        if(strcmp(norm, "true")){
            printf("\t%.2f|%.2f|%.2f|%.2f|%.2f|%.2f", nt[0]/tot_reads, nt[1]/tot_reads, nt[2]/tot_reads, nt[3]/tot_reads, nt[4]/tot_reads, nt[5]/tot_reads);               
        }else{
            printf("\t%.f|%.f|%.f|%.f|%.f|%.f", nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);               
        }
        printrow(html, region, nt);
        
        hts_itr_destroy(samitr);
        bam_destroy1(aln);
        sam_hdr_destroy(bamHdr);
        sam_close(fp_in);
}



void printhead(FILE *fn){
    fprintf(fn , "<!DOCTYPE html>\n");
    fprintf(fn, "<html lang=\"en\">");

    fprintf(fn, "<head>\n");
        fprintf(fn, "<meta charset=\"UTF-8\">\n");
        fprintf(fn, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");
        fprintf(fn , "<title>COSMIC genotypes</title>\n");
        fprintf(fn, "<link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css\">\n");
        fprintf(fn, "<script src=\"https://code.jquery.com/jquery-3.5.1.js\"></script>\n");
        fprintf(fn, "<script src=\"https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js\"></script>\n");
        fprintf(fn, "<style>\n");
        fprintf(fn, "body {padding: 0px; margin: 0; font-family: Verdana, Geneva, Tahoma, sans-serif;}\n");
        //fprintf(fn, "table {position: absolute; left: 50%%; top: 50%%; transform: translate(-50%%, -50%%); border-collapse: collapse; width: 800px; height: 200px; border: 1px solid #bdc3c7; box-shadow: 2px 2px 12px rgba(0, 0, 0, 0.2), -1px -1px 8px rgba(0, 0, 0, 0.2);}\n");
        fprintf(fn, "tr {transition: all .2s ease-in; cursor: pointer;}\n");
        fprintf(fn, "th, td {padding: 12px; text-align: left; border-bottom: 1px solid #ddd;}\n");
        fprintf(fn, "#header {background-color: #16a085; color: #fff;}\n");
        fprintf(fn, "h1 {font-weight: 600; text-align: center; background-color: #16a085; color: #fff; padding: 10px 0px;}    \n");
        fprintf(fn, "tr:hover {background-color: #f5f5f5; transform: scale(1.02); box-shadow: 2px 2px 12px rgba(0, 0, 0, 0.2), -1px -1px 8px rgba(0, 0, 0, 0.2);}\n");
        fprintf(fn, "@media only screen and (max-width: 768px) {table {width: 90%%;}}\n");
        fprintf(fn, "</style>\n");
        fprintf(fn, "<script>\n");
        fprintf(fn, "$(document).ready(function() {\n");
        fprintf(fn, "$('#cosmic').DataTable();});\n");
        fprintf(fn, "</script>\n");
    fprintf(fn, "</head>\n");

    fprintf(fn, "<h1>COSMIC genotypes </h1>\n");
    fprintf(fn, "<table id=\"cosmic\" class=\"display\">\n");
    fprintf(fn, "<thead><tr id=\"header\"><th>Loci</th><th>A</th><th>T</th><th>G</th><th>C</th><th>Ins</th><th>Del</th></tr></thead>\n");
    fprintf(fn, "<tbody>\n");
}

void printrow(FILE *fn,   char *loci, float counts[]){
    fprintf(fn, "<tr><td class='col1'>%s</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td></tr>\n", loci, counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]);
}

void printtail(FILE *fn){
    fprintf(fn, "</tbody>\n");
    fprintf(fn, "</table>\n");
    fprintf(fn, "</body>\n");
    fprintf(fn, "</html>\n");
}


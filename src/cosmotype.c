#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <math.h>

//void extbases(char *region, int32_t target_pos, char *bam, char *norm, FILE *html);
void is_bam(char *fname);
void printhead(FILE *fn);
void printrow(FILE *fn,  char *chrom, char *start, char *ref, char *alt, char *gene, char *vc, char *pc, char *cosid, float vaf, float counts[]);
void printtail(FILE *fn);

void usage(){
    fprintf(stderr, "-------------------------------------------------------------------------\n");
    fprintf(stderr, "cosmotype: A tool to genotype known COSMIC variants from BAM file\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "    cosmotype [OPTIONS] <loci> <bam>\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "    -f  Indexed fasta file. Extracts and adds reference base to the ouput\n");
    fprintf(stderr, "    -m  Map quality. Default 10\n");
    fprintf(stderr, "    -o  Output file. Default \"report\"\n");
    fprintf(stderr, "ARGS:\n");
    fprintf(stderr, "    <loci>   A 2 (or more) column tsv file with chr\\tpos\n");
    fprintf(stderr, "    <bam>    Indexed BAM file\n");
    fprintf(stderr, "-------------------------------------------------------------------------\n");
}

int main (int argc, char **argv){
  int m = 10;
  char *fafile = NULL;
  char *op = "report";
  int c;
  char *bedfile = NULL;
  char *bam = NULL;

  opterr = 0;

  while ((c = getopt (argc, argv, "m:f:o:")) != -1)
    switch (c)
      {
      case 'm':
        m = atoi(optarg);
        break;
      case 'f':
        fafile = optarg;
        break;
      case 'o':
        op = optarg;
        break;
      case '?':
        if (optopt == 'm'){
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            usage();
        }
          
        if (optopt == 'f'){
          usage();
          fprintf (stderr, "Input fasta file must provided when -%c is specified.\n", optopt); 
        } else if (isprint (optopt)){
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          usage();
        }else{
            fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            usage();
        }
          
        return 1;
      default:
        abort ();
      }

  
    bedfile = argv[optind];
    if(bedfile == NULL){
        usage();
        fprintf(stderr, "Missing loci file!\n");
        return 0;
    }
    bam = argv[optind+1];
    if(bam == NULL){
        usage();
        fprintf(stderr, "Missing BAM file!\n");
        return 0;
    }

    is_bam(bam);
    printf ("m = %d, fa = %s, ouput = %s, loci = %s, bam = %s\n", m, fafile, op, bedfile, bam);


    //Open bed file
    FILE *bed_fp;
    bed_fp = fopen(bedfile, "r");
    char buff[1000];
    //char *fract = "true";    

    //Open HTML report file and print HTML header
    FILE *html_fp;
    html_fp = fopen(strcat(op, ".html"), "w" );
    printhead(html_fp);

    //fasta file
    char *seq;
    faidx_t *fa = fai_load(fafile);

    //BAM file
    samFile *fp_in = hts_open(bam,"r"); //open bam file
    hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    

    //For every loci in the BED file
    while(fgets(buff,1000,bed_fp) != NULL){
        //Remove trailing new line chars
        int len = strlen(buff);
        if(buff[len-1] == '\n' ){
            buff[len-1] = 0;
        }
            
        char *chrom = strtok(buff,"\t");
        char *start = strtok(NULL,"\t");
        char *ref = strtok(NULL,"\t");
        char *alt = strtok(NULL,"\t");
        char *gene = strtok(NULL,"\t");
        char *vc = strtok(NULL,"\t");
        char *pc = strtok(NULL,"\t");
        char *cosid = strtok(NULL,"\t");

        //char *end = strtok(NULL,"\t");
        char loci[250] = "";

        strcat(loci, chrom); strcat(loci, ":"); strcat(loci, start); strcat(loci, "-"); strcat(loci, start);
        //printf("%s\n", loci);
        
        //Fetch base at target loci from fasta file 
        if(fa != NULL){
            int templen = 100;
            seq = fai_fetch(fa, loci, &templen);
            printf("%s:%s\t%s\t%s>%s\t%s\t%s\t%s\t%s", chrom, start, seq, ref, alt, gene, vc, pc, cosid);
            free(seq);
        }else{
            printf("%s:%s\tNA\t%s>%s\t%s\t%s\t%s\t%s", chrom, start, ref, alt, gene, vc, pc, cosid);
        }
        
        int32_t target_pos = atoi(start) -1; //input position are 1 based
        //extbases(loci, target_pos, bam, fract, html_fp);
        
        //load reads in target loci
        hts_itr_t *samitr = sam_itr_querys(fp_idx, bamHdr, loci);
        
        //Keep track of total reads and nt counts per loci
        int32_t tot_reads = 0;
        float nt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        char *gt = "NA";
        float vaf = 0.0;

        //For every read in the BAM file of target region
        while(sam_itr_next(fp_in, samitr, aln) > 0){
            
            int32_t pos = aln->core.pos ; //left most position of alignment in zero based coordianate (0-based)
            //char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
            uint32_t len = aln->core.l_qseq; //length of the read.
            uint32_t* cig = bam_get_cigar(aln);
                
            uint8_t *q = bam_get_seq(aln); //quality string
            uint32_t q2 = aln->core.qual ; //mapping quality


            //MAPQ filter
            if(q2 <= m){
                break;
            }
            
            tot_reads = tot_reads +1;

            //get nucleotide id and converts them into IUPAC id.
            char *qseq = (char *)malloc(len);
            for(int i=0; i< len ; i++){
                qseq[i] = seq_nt16_str[bam_seqi(q,i)]; 
            }
            
            //target position on the read
            int32_t pos_onread = 0;

            //For every CIGAR string
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

        if(strcmp(alt, "A")){
            vaf = nt[0]/tot_reads;
        }else if(strcmp(alt, "T")){
            vaf = nt[1]/tot_reads;
        }else if(strcmp(alt, "G")){
            vaf = nt[2]/tot_reads;
        }else if(strcmp(alt, "C")){
            vaf = nt[3]/tot_reads;
        }else{
            vaf = 0.0;
        }

        hts_itr_destroy(samitr);
        printf("\t%.f|%.f|%.f|%.f|%.f|%.f|%.f\n", vaf, nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);
        //wite nt counts as table row
        if(!isnan(vaf)){
            if(vaf > 0.05){
                printrow(html_fp, chrom, start, ref, alt, gene, vc, pc, cosid, vaf, nt);
            }
        }
    } 
    
    //Close all open connections and destroy destroy objects    
    printtail(html_fp);
    bam_destroy1(aln);
    sam_hdr_destroy(bamHdr);
    sam_close(fp_in);
    fclose(bed_fp);
    fai_destroy(fa);
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
    fprintf(fn, "<thead><tr id=\"header\"><th>Chr</th><th>Pos</th><th>Ref</th><th>Alt</th><th>Gene</th><th>Type</th><th>HGVSp</th><th>COSMIC ID</th><th>VAF</th><th>A</th><th>T</th><th>G</th><th>C</th><th>Ins</th><th>Del</th></tr></thead>\n");
    fprintf(fn, "<tbody>\n");
}

void printrow(FILE *fn,  char *chrom, char *start, char *ref, char *alt, char *gene, char *vc, char *pc, char *cosid, float vaf, float counts[]){
    fprintf(fn, "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td></tr>\n", 
        chrom, start, ref, alt, gene, vc, pc, cosid, vaf, counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]);
}

void printtail(FILE *fn){
    fprintf(fn, "</tbody>\n");
    fprintf(fn, "</table>\n");
    fprintf(fn, "</body>\n");
    fprintf(fn, "</html>\n");
}


// void extbases(char *region, int32_t target_pos, char *bam, char *norm, FILE *html){
    
//     samFile *fp_in = hts_open(bam,"r"); //open bam file
//     hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
//     bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header

//     bam1_t *aln = bam_init1(); //initialize an alignment
//     hts_itr_t *samitr = sam_itr_querys(fp_idx, bamHdr, region);    

//     int32_t tot_reads = 0;
//     float nt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//     while(sam_itr_next(fp_in, samitr, aln) > 0){
            
//             int32_t pos = aln->core.pos ; //left most position of alignment in zero based coordianate (+1)
//             //char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
//             uint32_t len = aln->core.l_qseq; //length of the read.
//             uint32_t* cig = bam_get_cigar(aln);
                
//             uint8_t *q = bam_get_seq(aln); //quality string
//             uint32_t q2 = aln->core.qual ; //mapping quality


//             if(q2 <= 10){
//                 break;
//             }

//             char *qseq = (char *)malloc(len);
            
//             for(int i=0; i< len ; i++){
//                 qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
//             }
//             int32_t pos_onread = 0;

//             for( int k=0;k< aln->core.n_cigar ;++k){
//                 int cop =cig[k] & BAM_CIGAR_MASK; // CIGAR string
//                 int cl = cig[k] >> BAM_CIGAR_SHIFT; // CIGAR length

//                 if(BAM_CIGAR_STR[cop] == 'M'){
//                     pos_onread = pos_onread + cl;
//                     pos = pos + cl;
//                 }else if(BAM_CIGAR_STR[cop] == 'S'){
//                     pos_onread = pos_onread + cl;
//                 }else if(BAM_CIGAR_STR[cop] == 'I'){
//                     pos_onread = pos_onread + cl;
//                 }
                

//                 if(pos > target_pos){
//                     //printf("%c\n", BAM_CIGAR_STR[cop]);
//                     if(BAM_CIGAR_STR[cop] == 'M'){
//                         pos_onread = pos_onread - (pos - target_pos);
//                         //printf("%d\t%lld\t%d\t%c\n", pos, aln->core.pos, tar_pos_onread, qseq[pos - target_pos]);
//                         if(qseq[pos_onread] == 'A'){
//                             nt[0] = nt[0] + 1;
//                         }else if(qseq[pos_onread] == 'T'){
//                             nt[1] = nt[1] + 1;
//                         }else if(qseq[pos_onread] == 'G'){
//                             nt[2] = nt[2] + 1;
//                         }else if(qseq[pos_onread] == 'C'){
//                             nt[3] = nt[3] + 1;
//                         }
//                         break;
//                     }
//                 }else if(pos == target_pos){
//                     if(BAM_CIGAR_STR[cop] == 'I'){
//                         nt[4] = nt[4] + 1;
//                         // insertion sequence
//                         // for(int i = 0; i < cl; i++){
//                         //     printf("%c", qseq[pos_onread+i]);
//                         // }
//                         break;
//                     }
//                 }else if(pos == target_pos){
//                     if(BAM_CIGAR_STR[cop] == 'D'){
//                         nt[5] = nt[5] + 1;
//                         // deletion sequence
//                         // for(int i = 0; i < cl; i++){
//                         //     printf("%c", qseq[pos_onread+i]);
//                         // }
//                         break;
//                     }
//                 }                
//             }
//         }    
//         if(strcmp(norm, "true")){
//             printf("\t%.2f|%.2f|%.2f|%.2f|%.2f|%.2f", nt[0]/tot_reads, nt[1]/tot_reads, nt[2]/tot_reads, nt[3]/tot_reads, nt[4]/tot_reads, nt[5]/tot_reads);               
//         }else{
//             printf("\t%.f|%.f|%.f|%.f|%.f|%.f", nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);               
//         }
//         printrow(html, region, nt);
        
//         hts_itr_destroy(samitr);
//         bam_destroy1(aln);
//         sam_hdr_destroy(bamHdr);
//         sam_close(fp_in);
// }




#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <math.h>
#include <argp.h>
#include <time.h>

void usage(){
    fprintf(stderr, "-------------------------------------------------------------------------\n");
    fprintf(stderr, "ntfreq: A tool to extract nucleotide counts of targeted loci from BAM file(s)\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "    ntfreq [OPTIONS] <loci> <bam>...\n");
    fprintf(stderr, "    e.g; ntfreq cancer_hotspots_hg19.tsv Normal.bam Primary_tumor.bam Relapse_tumor.bam\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "    -f  Indexed fasta file. If provided, extracts and adds reference base to the ouput\n");
    fprintf(stderr, "    -q  Minimum mapping quality threshold. Default 10 [Read filter]\n");
    fprintf(stderr, "    -F  Exclude reads with FLAGS >= INT. Default 1024 (i.e, read is PCR or optical duplicate)\n");
    fprintf(stderr, "ARGS:\n");
    fprintf(stderr, "    <loci>   A 2 (or more) column tsv file with chr\\tpos\n");
    fprintf(stderr, "    <bam>... One or more BAM files\n");
    fprintf(stderr, "-------------------------------------------------------------------------\n");
}

void is_bam(char *fname){
    int l = strlen(fname);
    if (l>=4 && strcasecmp(fname+l-4, ".bam") != 0){
        usage();
        fprintf(stderr, "%s not a BAM file!\n", fname);
        exit(0);
    }
}

char *basename(char const *path){
    char *s = strrchr(path, '/');
    if (!s){
        return strdup(path);
    }else{
        return strdup(s + 1);
    }   
}

void extbases(char *region, int32_t target_pos, char *bam, int q, int flag){
    
    samFile *fp_in = hts_open(bam,"r"); //open bam file
    hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header

    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_itr_t *samitr = sam_itr_querys(fp_idx, bamHdr, region);    

    float nt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    while(sam_itr_next(fp_in, samitr, aln) > 0){
            
            int32_t pos = aln->core.pos ; //left most position of alignment in zero based coordianate (+1)
            uint32_t len = aln->core.l_qseq; //length of the read.
            uint32_t* cig = bam_get_cigar(aln);
            uint8_t *qs = bam_get_seq(aln); //quality string


            //MAPQ and FLAG filter
            if(aln->core.qual <= q){
                break;
            }

            if(aln->core.qual >= flag){
                break;
            }
            

            char *qseq = (char *)malloc(len);
            
            for(int i=0; i< len ; i++){
                qseq[i] = seq_nt16_str[bam_seqi(qs,i)]; //gets nucleotide id and converts them into IUPAC id.
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
                    pos = pos + cl;
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
        printf("\t%.f|%.f|%.f|%.f|%.f|%.f", nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);
        //printrow(html, region, nt);
        
        hts_itr_destroy(samitr);
        bam_destroy1(aln);
        sam_hdr_destroy(bamHdr);
        sam_close(fp_in);
}


int main(int argc, char *argv[]){

  char *fafile = NULL;  
  int c;
  char *bedfile = NULL;
  char *bam = NULL;

  uint32_t q = 10; //MAPQ
  uint16_t F = 1024; //BAM default FLAG

  opterr = 0;
  hts_verbose = 0; //suppresses htslib warnings

  //Parse cl args
  while ((c = getopt (argc, argv, "m:f:F:t:d:v:o:")) != -1){
      switch (c)
      {
      case 'm':
        q = atoi(optarg);
        break;
      case 'f':
        fafile = optarg;
        break;
      case 'F':
        F = atoi(optarg);
      case '?':
        if (optopt == 'm'){
            fprintf (stderr, "MAPQ must be provided when -%c is specified.\n", optopt);
            usage();
        }
    
        if (optopt == 'F'){
            fprintf (stderr, "FLAG must be provided when -%c is specified.\n", optopt);
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
  }

    //fasta file
    char *seq;
    faidx_t *fa = fai_load(fafile);
  
    bedfile = argv[optind];
    if(bedfile == NULL){
        usage();
        fprintf(stderr, "Missing loci file!\n");
        return 0;
    }

    optind = optind+1;
    bam = argv[optind];
    if(bam == NULL){
        usage();
        fprintf(stderr, "Missing BAM file! At least one BAM is required.\n");
        return 0;
    }else{
        int i = 0;
        for(int i= optind ; i < argc ; i++){
            is_bam(argv[i]);
        }

        printf("chr\tpos\tref");
        int j = 0;
        for(int j= optind ; j < argc ; j++){
            printf("\t%s", basename(argv[j]));
        }
        printf("\n");
    }

    FILE *bed_fp;
    bed_fp = fopen(bedfile, "r");
    char buff[1000];

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
        printf("%s\t%s\t", chrom, start);
        //Fetch base at target loci from fasta file 
        if(fa != NULL){
            int templen = 100;
            seq = fai_fetch(fa, loci, &templen);
            printf("%s\t", seq);
        }else{
            printf("N\t");
        }

        int i = 0;
        for(int i= optind ; i < argc ; i++){
            extbases(loci, target_pos, argv[i], q, F);
        }
        printf("\n");
        free(seq);
    }

    fclose(bed_fp);
    fai_destroy(fa);

    return 0;
}



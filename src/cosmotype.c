#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <math.h>
#include <argp.h>
#include <time.h>

void is_bam(char *fname);
void printhead(FILE *fn, char *bam_fn);
void printrow(FILE *fn,  char *chrom, char *start, char *ref, char *alt, char *gene, char *vc, char *pc, char *cosid, float vaf, float counts[]);
void printargtbl(FILE *fn, int mapq, char *bam, int tot_var, int som_vars, float vaf, int cov, float doc, int treads);
char *basename(char const *path);
char *removeExt(char* myStr);
void usage();

int main (int argc, char **argv){
  uint32_t q = 10;
  int d = 30;
  int t = 8;
  char *fafile = NULL;
  char *op = NULL;
  int c;
  char *bedfile = NULL;
  char *bam = NULL;
  float v = 0.05;
  uint16_t F = 1024;

  int vars_gt = 0; //No. of BED entries
  int som_vars = 0; //vars with somatic evidance
  float mean_doc = 0.0; //mean depth of coverage
  int novaf = 0; //entries for which VAF could not be estimated

  //int nrows = 89496; //for progress bar

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
      case 'd':
        d = atoi(optarg);
        break;
      case 'o':
        op = optarg;
        break;
      case 'v':
        v = atof(optarg);
      case 'F':
        F = atoi(optarg);
      case 't':
        t = atoi(optarg);
        break;
      case '?':
        if (optopt == 'm'){
            fprintf (stderr, "MAPQ must be provided when -%c is specified.\n", optopt);
            usage();
        }
        if (optopt == 'v'){
            fprintf (stderr, "VAF must be provided when -%c is specified.\n", optopt);
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

    is_bam(bam); //Check the file extension to to figure out BAM file or not
    
    //Output files
    char html_file[1000];
    char tsv_file[1000];
    char *bn = basename(bam); //get basename of the file
    char *bnn = removeExt(bn); //remove extension from the basename
    if(op == NULL){
        strcpy(html_file, bnn);
        strcat(html_file, ".html");
        strcpy(tsv_file, bnn);
        strcat(tsv_file, ".tsv");
    }else{
        strcpy(html_file, op);
        strcat(html_file, ".html");
        strcpy(tsv_file, op);
        strcat(tsv_file, ".tsv");
    }

    //Open bed file
    FILE *bed_fp;
    bed_fp = fopen(bedfile, "r");
    char buff[1000];
    //char *fract = "true";

    //Open HTML report file and print HTML header
    FILE *html_fp;
    html_fp = fopen(html_file, "w" );
    printhead(html_fp, bnn);

    //Open TSV report file and print HTML header
    FILE *tsv_fp;
    tsv_fp = fopen(tsv_file, "w" );
    fprintf(tsv_fp, "loci\tfa_ref\tCOSMIC_nt_change\tHugo_Symbol\tVariant_Classification\tHGVSp\tCOSMIC_ID\tVAF\tA\tT\tG\tC\tIns\tDel\n");

    //fasta file
    char *seq;
    faidx_t *fa = fai_load(fafile);

    //BAM file
    samFile *fp_in = hts_open(bam,"r"); //open bam file
    hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment

    fprintf(stderr, "Input BAM file           :  %s\n", bam);
    fprintf(stderr, "COSMIC variants          :  %s\n", bedfile);
    fprintf(stderr, "VAF filter               :  %.2f\n", v);
    fprintf(stderr, "min reads for t_allele   :  %d\n", t);
    fprintf(stderr, "MAPQ filter              :  %d\n", q);
    fprintf(stderr, "FLAG filter              :  %d\n", F);
    fprintf(stderr, "Coverage filter          :  %d\n", d);
    fprintf(stderr, "HTSlib version           :  %s\n", hts_version());
    
    //For every loci in the BED file
    while(fgets(buff,1000,bed_fp) != NULL){
        //pb_len = pb_len +1;
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
            fprintf(tsv_fp, "%s:%s\t%s\t%s>%s\t%s\t%s\t%s\t%s", chrom, start, seq, ref, alt, gene, vc, pc, cosid);
            free(seq);
        }else{
            fprintf(tsv_fp, "%s:%s\tNA\t%s>%s\t%s\t%s\t%s\t%s", chrom, start, ref, alt, gene, vc, pc, cosid);
        }
        
        int32_t target_pos = atoi(start) -1; //input position are 1 based
        //extbases(loci, target_pos, bam, fract, html_fp);
        
        //load reads in target loci
        hts_itr_t *samitr = sam_itr_querys(fp_idx, bamHdr, loci);
        
        //Keep track of total reads and nt counts per loci
        int32_t tot_reads = 0;
        float nt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        //char *gt = "NA";
        float vaf = 0.0;

        vars_gt = vars_gt + 1;

        if(vars_gt % 1000 == 0){
            fprintf(stderr, "Processed %d entries..\n", vars_gt);
        }

        //For every read in the BAM file of target region
        while(sam_itr_next(fp_in, samitr, aln) > 0){
            
            int32_t pos = aln->core.pos ; //left most position of alignment in zero based coordianate (0-based)
            //char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
            uint32_t len = aln->core.l_qseq; //length of the read.
            uint32_t* cig = bam_get_cigar(aln);
            //uint16_t samflag = aln->core.flag; // flag
            //uint32_t q2 = aln->core.qual ; //mapping quality    
            uint8_t *qs = bam_get_seq(aln); //quality string
            

            //MAPQ and FLAG filter
            if(aln->core.qual <= q){
                break;
            }

            if(aln->core.qual >= F){
                break;
            }
            
            tot_reads = tot_reads +1;
            //char ins_seq[100];
            //char del_seq[100];

            //get nucleotide id and converts them into IUPAC id.
            char *qseq = (char *)malloc(len);
            int i = 0;
            for(i=0; i< len ; i++){
                qseq[i] = seq_nt16_str[bam_seqi(qs,i)]; 
            }
            
            //target position on the read
            int32_t pos_onread = 0;

            //For every CIGAR string
            int k = 0;
            for(k=0;k< aln->core.n_cigar ;++k){
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

                if(pos > target_pos){
                    if(BAM_CIGAR_STR[cop] == 'M'){
                        pos_onread = pos_onread - (pos - target_pos);
                        //printf("%d\t%lld\t%d\t%c\n", pos, aln->core.pos, pos_onread, qseq[pos - target_pos]);
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
                        //     //printf("%c", qseq[pos_onread+i]);
                        //     strcat(ins_seq, qseq[pos_onread+i]);
                        // }
                        break;
                    }else if(BAM_CIGAR_STR[cop] == 'D'){
                        nt[5] = nt[5] + 1;
                        // deletion sequence
                        // for(int i = 0; i < cl; i++){
                        //     //printf("%c", qseq[pos_onread+i]);
                        //     strcat(del_seq, qseq[pos_onread+i]);
                        // }
                        break;
                    }
                }                
            }

            free(qseq);
        }

        int pass_treads = 1; //Tumor supporting reads
        if(strcmp(vc, "Frame_Shift_Ins") == 0){
            vaf = nt[4]/tot_reads;
        }else if(strcmp(vc, "Frame_Shift_Del") == 0){
            vaf = nt[5]/tot_reads;
        }else{
            if(strcmp(alt, "A") == 0 && strcmp(ref, "-") != 0){
                if(nt[0] < t){
                    pass_treads = 0;
                }
                vaf = nt[0]/tot_reads;
            }else if(strcmp(alt, "T") == 0 && strcmp(ref, "-") != 0){
                if(nt[1] < t){
                    pass_treads = 0;
                }
                vaf = nt[1]/tot_reads;
            }else if(strcmp(alt, "G") == 0 && strcmp(ref, "-") != 0){
                if(nt[2] < t){
                    pass_treads = 0;
                }
                vaf = nt[2]/tot_reads;
            }else if(strcmp(alt, "C") == 0 && strcmp(ref, "-") != 0){
                if(nt[3] < t){
                    pass_treads = 0;
                }
                vaf = nt[3]/tot_reads;
            }else{
                vaf = -1; 
                novaf = novaf + 1;
            }

        }
        

        mean_doc = mean_doc + tot_reads;

        hts_itr_destroy(samitr);
        fprintf(tsv_fp, "\t%.3f\t%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", vaf, nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);
        //write nt counts as table row
        if(!isnan(vaf)){
            //fprintf(stderr, "%s\t%s\t%.3f\n", ref, alt, vaf);
            if(vaf > v && tot_reads > d && pass_treads == 1){
                som_vars = som_vars + 1;
                printrow(html_fp, chrom, start, ref, alt, gene, vc, pc, cosid, vaf, nt);
            }
        }
    } 
    
    mean_doc = mean_doc/vars_gt;
    //Close all open connections and destroy destroy objects    
    printargtbl(html_fp, q, bam, vars_gt, som_vars, v, d, mean_doc, t);
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);
    fclose(bed_fp);
    fai_destroy(fa);
    fclose(html_fp);
    fclose(tsv_fp);

    free(bn);
    free(bnn);

    fprintf(stderr, "Done!\n");
    fprintf(stderr, "Summary:\n");
    fprintf(stderr, "Total variants processed :  %d\n", vars_gt);
    fprintf(stderr, "Variants > %.2f threshold:  %d\n", v, som_vars);
    fprintf(stderr, "Avg. depth of coverage   :  %.2f\n", mean_doc);
    fprintf(stderr, "Output html report       :  %s\n", html_file);
    fprintf(stderr, "Output TSV file          :  %s\n", tsv_file);

    if(novaf > 0){
        fprintf(stderr, "Could not estimate VAF for %d variants. VAF has been set to -1 in %s\n.Manually inspect them in IGV.\n", novaf, tsv_file);
    }

    return 0;
}

void is_bam(char *fname){
    int l = strlen(fname);
    if (l>=4 && strcasecmp(fname+l-4, ".bam") != 0){
        usage();
        fprintf(stderr, "%s is not a BAM file!\n", fname);
        exit(0);
    }else if(l<4){
        fprintf(stderr, "%s is not a BAM file!\n", fname);
        usage();
        exit(0);
    }
}

//from: https://stackoverflow.com/questions/3288006/are-there-any-c-apis-to-extract-the-base-file-name-from-its-full-path-in-linux
char *basename(char const *path){
    char *s = strrchr(path, '/');
    if (!s){
        return strdup(path);
    }else{
        return strdup(s + 1);
    }   
}

//from: https://stackoverflow.com/questions/2736753/how-to-remove-extension-from-file-name
char *removeExt(char* myStr) {
    char *retStr;
    char *lastExt;
    if (myStr == NULL) return NULL;
    if ((retStr = malloc (strlen (myStr) + 1)) == NULL) return NULL;
    strcpy (retStr, myStr);
    lastExt = strrchr (retStr, '.');
    if (lastExt != NULL)
        *lastExt = '\0';
    return retStr;
}

void printhead(FILE *fn, char *bam_fn){
    fprintf(fn , "<!DOCTYPE html>\n");
    fprintf(fn, "<html lang=\"en\">");

    fprintf(fn, "<head>\n");
        fprintf(fn, "<meta charset=\"UTF-8\">\n");
        fprintf(fn, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");
        fprintf(fn , "<title>%s | COSMIC variants</title>\n", bam_fn);
        fprintf(fn, "<link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css\">\n");
        fprintf(fn, "<script src=\"https://code.jquery.com/jquery-3.5.1.js\"></script>\n");
        fprintf(fn, "<script src=\"https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js\"></script>\n");
        fprintf(fn, "<style>\n");
        fprintf(fn, "body {padding: 0px; margin: 0; font-family: Verdana, Geneva, Tahoma, sans-serif;}\n");
        //fprintf(fn, "table {position: absolute; left: 50%%; top: 50%%; transform: translate(-50%%, -50%%); border-collapse: collapse; width: 800px; height: 200px; border: 1px solid #bdc3c7; box-shadow: 2px 2px 12px rgba(0, 0, 0, 0.2), -1px -1px 8px rgba(0, 0, 0, 0.2);}\n");
        fprintf(fn, "tr {transition: all .2s ease-in; cursor: pointer; color: #2c3e50;}\n");
        fprintf(fn, "th, td {padding: 12px; text-align: left; border-bottom: 1px solid #ddd;}\n");
        fprintf(fn, "#header {background-color: #2c3e50; color: #fff;}\n");
        fprintf(fn, "h1 {font-weight: 600; text-align: left; color: #2c3e50; padding: 10px 0px;}    \n");
        fprintf(fn, "h3 {text-align: left; color: #c0392b}\n");
        fprintf(fn, "tr:hover {background-color: #f5f5f5; transform: scale(1.02); box-shadow: 2px 2px 12px rgba(0, 0, 0, 0.2), -1px -1px 8px rgba(0, 0, 0, 0.2);}\n");
        fprintf(fn, "@media only screen and (max-width: 768px) {table {width: 90%%;}}\n");
        fprintf(fn, ".details tr { line-height: 15px; }\n");
        fprintf(fn, ".details th, td {padding: 5px; text-align: left; border-bottom: 1px solid #ddd; font-family:'Courier New', Courier, monospace;}\n");
        fprintf(fn, "</style>\n");
        fprintf(fn, "<script>\n");
        fprintf(fn, "$(document).ready(function() {\n");
        fprintf(fn, "$('#cosmic').DataTable();});\n");
        fprintf(fn, "</script>\n");
    fprintf(fn, "</head>\n");

    fprintf(fn, "<h1>COSMIC variants </h1>\n");
    fprintf(fn, "<table id=\"cosmic\" class=\"display\">\n");
    fprintf(fn, "<thead><tr id=\"header\"><th>Chr</th><th>Pos</th><th>Ref</th><th>Alt</th><th>Gene</th><th>Type</th><th>HGVSp</th><th>COSMIC ID</th><th>VAF</th><th>A</th><th>T</th><th>G</th><th>C</th><th>Ins</th><th>Del</th></tr></thead>\n");
    fprintf(fn, "<tbody>\n");
}


void printargtbl(FILE *fn, int mapq, char *bam, int tot_var, int som_vars, float vaf, int cov, float doc, int treads){
    
    time_t t = time(NULL);
    struct tm cutime = *localtime(&t);

    fprintf(fn, "</table>\n"); //end previous table before printing summary
    fprintf(fn, "<h3 >Summary </h3>");
    fprintf(fn, "<table class=\"details\">\n");
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Date</td><td>%d-%02d-%02d %02d:%02d:%02d</td></tn>", cutime.tm_year + 1900, cutime.tm_mon + 1, cutime.tm_mday, cutime.tm_hour, cutime.tm_min, cutime.tm_sec);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">BAM</td><td>%s</td></tn>", bam);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">COSMIC variants queried</td><td>%d</td></tn>", tot_var);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">COSMIC variants with somatic evidance</td><td>%d [VAF > %.2f]</td></tn>", som_vars, vaf);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Avg. depth of coverage</td><td>%.2f</td></tn>", doc);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">MAPQ filter</td><td>%d</td></tn>", mapq);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Coverage filter</td><td>%d</td></tn>", cov);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">VAF filter</td><td>%.2f</td></tn>", vaf);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Min. number of reads supporting tumor allele </td><td>%d</td></tn>", treads);
    fprintf(fn, "</table>\n");
    fprintf(fn, "</tbody>\n");
    fprintf(fn, "</body>\n");
    fprintf(fn, "</html>\n");

}

void printrow(FILE *fn,  char *chrom, char *start, char *ref, char *alt, char *gene, char *vc, char *pc, char *cosid, float vaf, float counts[]){
    fprintf(fn, "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td><a target=\"_blank\" href=\"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=%s#variants\"> %s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td><td>%.f</td></tr>\n", 
        chrom, start, ref, alt, gene, gene, vc, pc, cosid, vaf, counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]);
}

void usage(){
    fprintf(stderr, "-------------------------------------------------------------------------\n");
    fprintf(stderr, "cosmotype: A tool to extract nucleotide/variant allele frequencies of\n");
    fprintf(stderr, "           known PATHOGENIC and non-SNP COSMIC variants from BAM file\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "    cosmotype [OPTIONS] <loci> <bam>\n");
    fprintf(stderr, "    e.g; cosmotype COSMIC_nsyn.tsv Tumor.bam\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "    -f  Indexed fasta file. If provided, extracts and adds reference base to the ouput tsv\n");
    fprintf(stderr, "    -q  Mapping quality threshold. Default 10 [Read filter]\n");
    fprintf(stderr, "    -F  Exclude reads with FLAGS >= INT. Default 1024 (i.e, read is PCR or optical duplicate) [Read filter]\n");
    fprintf(stderr, "    -v  VAF threshold. Default 0.05 [Variant filter]\n");
    fprintf(stderr, "    -d  Depth of coverage threshold. Default 30 [Variant filter]\n");
    fprintf(stderr, "    -t  Min. number of reads supporting tumor allele . Default 8 [Variant filter]\n");
    fprintf(stderr, "    -o  Output file basename. Default parses from basename of BAM file\n");
    fprintf(stderr, "ARGS:\n");
    fprintf(stderr, "    <loci>   COSMIC variant file\n");
    fprintf(stderr, "    <bam>    Indexed BAM file\n");
    fprintf(stderr, "OUPUT:\n");
    fprintf(stderr, "    <output>.html   A browsable html report with COSMIC variants post variant filter\n");
    fprintf(stderr, "    <output>.tsv    TSV file with nucletide counts and VAF for all COSMIC variants\n");
    fprintf(stderr, "\nDETAILS:\n");
    fprintf(stderr, "<loci> is a tsv file with 8 columns parsed from COSMIC database. Ex. row:\n");
    fprintf(stderr, "Positions are in one-based coordinate systems. Example rows:\n");
    fprintf(stderr, "1	115258747	C	T	NRAS	Missense	p.G12D	COSV54736383\n");
    fprintf(stderr, "1	115258747	C	A	NRAS	Missense	p.G12V	COSV54736974\n");
    fprintf(stderr, "-------------------------------------------------------------------------\n");
}
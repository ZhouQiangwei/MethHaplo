#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <map>
#include <math.h>
#include <limits.h>
#include "../submodules/libbm/binaMeth.h"
#define BATBUF 3000
FILE* File_Open(const char* File_Name,const char* Mode);
bool isunsigned (char* str)
{
	for(int i = 0; i < strlen(str); i++)
	{
		if(!isdigit(str[i]))
		{
			return false;
		}
	}
	return true;
}
char strand_str[] = {'+', '-', '.'};
char *context_str[] = {"C", "CG", "CHG", "CHH"};
int main(int argc,char *argv[])
{
	const char* Help_String="Command Format :  splitmr [options] -i <methratiofile>\n"
        "\nMethyHaplo::splitmr \nversion: v1.0\n\n"
		"\n======================== PROGRAM OPTIONS ========================\n\n"
		"\t-i                    input file.  This file is the result of methyhaplo.\n"
        "\t--bm                  Y/N, bm format of methratio\n"
        "\t-c [ALL/CG/CHG/CHH]   Meth context. [default: ALL]\n"
		"\t-m [int]              Meth coverage. [default: 2]\n"
        "\t-n [int]              Total coverage. [default: 6]\n"
		"\t-f|--Mf [double]      Methy level. [default: 0.2]\n"
		"\t-h|--help\n";
	
	char* methratioIN;
    double weight=0.2;
    int mCover = 2;
    int TotalCover=6;
    char Mcontext[100]; strcpy(Mcontext, "ALL");
    char bmformat[100];
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-i")   )
		{
			methratioIN=argv[++i];
		}else if(!strcmp(argv[i], "--bm")   )
        {
            strcpy(bmformat, argv[++i]);
        }else if(!strcmp(argv[i], "-f") || !strcmp(argv[i], "--Mf")   )
		{
			weight=atof(argv[++i]);
		}else if(!strcmp(argv[i], "-m") )
		{
			mCover=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-n") )
        {
            TotalCover=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "-c") )
        {
            strcpy(Mcontext, argv[++i]);
        }else
		{
			printf("%s \n",Help_String);
			exit(0);
		}
	}
	if(argc<=1) 
	{
		printf("%s \n",Help_String);
		exit(0);
	}
    char methOUTp[1024];
    char methOUTn[1024];
    sprintf(methOUTp, "%s.p", methratioIN);
    sprintf(methOUTn, "%s.n", methratioIN);
	FILE* OutfileP = File_Open(methOUTp,"w");
    FILE* OutfileN = File_Open(methOUTn,"w");
	printf("[MethyHaplo] Processing methratio file ...\n");
	
	
	char Dummy[BATBUF];
	char chrom[100];
    char strand;char context[10];
    int meth=0;int cover=0;float ml=0;
    int pos=0;
    int nline=0;
    float hweight = 1.0-weight;
    if(strcmp(bmformat, "Y")==0){
        binaMethFile_t *ifp = NULL;
        ifp = bmOpen(methratioIN, NULL, "r");
        ifp->type = ifp->hdr->version;
        int SEGlen = 1000000;
        int start = 0, end = SEGlen-1;
        bmOverlappingIntervals_t *o;
        int i=0,j=0;
        for(i=0;i<ifp->cl->nKeys;i++){
            strcpy(chrom, (char*)ifp->cl->chrom[i]);
            int len = (int)ifp->cl->len[i];
            start = 0, end = SEGlen-1;
            while(start<len){
                if(end>len){
                    end = len;
                }
                o = bmGetOverlappingIntervals(ifp, chrom, start, end+1);
                if(o->l) {
                    for(j=0; j<o->l; j++) {
                        if(ifp->hdr->version & BM_COVER){
                            if(!(o->coverage[j]>=TotalCover)){
                                continue;
                            }
                            meth = (int)(o->coverage[j]*o->value[j]+0.5);
                        }else{
                            fprintf(stderr, "Warning: your input BM file without covarage information, please check!");
                            meth = 3;
                        }
                        ml = o->value[j];
                        if(ifp->hdr->version & BM_CONTEXT){
                            strcpy(context, context_str[o->context[j]]);
                        }
                        nline++;
                        if(nline%1000000 == 0){
                            fprintf(stderr, "[MethProcess] Processed %d meth loci\n", nline);
                        }
                        if(meth < mCover) continue;
                        if(!( ml >= weight && ml <= hweight ) && !(meth>2 && mCover-meth>2) ) continue;

                        pos = o->start[j];
                        if(strcmp(Mcontext, "ALL") == 0 || strcmp(Mcontext, context) == 0){
                            strand = strand_str[o->strand[j]];
                            if(strand == '+'){
                                //$1,$2,".\tC\tT\t600\tPASS\t.\tGT\t0/1"
                                fprintf(OutfileP, "%s\t%d\t.\tC\tT\t600\tPASS\t.\tGT\t0/1\n", chrom, pos);
                            }else if(strand == '-'){
                                fprintf(OutfileN, "%s\t%d\t.\tG\tA\t600\tPASS\t.\tGT\t0/1\n", chrom, pos);
                            }else{
                                fprintf(stderr, "Unexpected strand value: %c", strand);
                                exit(1);
                            }
                        }
                    }
                }
                start += SEGlen;
                end += SEGlen;
            }
        }

        bmClose(ifp);
    }else{
        FILE* INfile = File_Open(methratioIN,"r");
	    while(fgets(Dummy,BATBUF,INfile)!=0)
	    {
	    	if(Dummy[0]=='#') 
	    	{
	    		continue;
	    	}

	    	//chrM    1       -       CHH     2       1795    0.001114
	    	sscanf(Dummy,"%s%d%s%s%d%d%f",&chrom,&pos,&strand,&context,&meth,&cover,&ml);
            nline++;
            if(nline%1000000 == 0){
                fprintf(stderr, "[MethProcess] Processed %d meth loci\n", nline);
            }
            //if(ml < weight || ml > hweight) continue;
            if(meth < mCover) continue;
            if(cover < TotalCover) continue;
            if(!( ml >= weight && ml <= hweight ) && !(meth>2 && mCover-meth>2) ) continue;

	    	if(strcmp(Mcontext, "ALL") == 0 || strcmp(Mcontext, context) == 0){
                if(strand == '+'){
                    //$1,$2,".\tC\tT\t600\tPASS\t.\tGT\t0/1"
                    fprintf(OutfileP, "%s\t%d\t.\tC\tT\t600\tPASS\t.\tGT\t0/1\n", chrom, pos);
                }else if(strand == '-'){
                    fprintf(OutfileN, "%s\t%d\t.\tG\tA\t600\tPASS\t.\tGT\t0/1\n", chrom, pos);
                }else{
                    fprintf(stderr, "Unexpected strand value: %c", strand);
                    exit(1);
                }
            }
	    }
	    //remove(methratioIN);
        fclose(INfile);
    }
    fclose(OutfileP);
    fclose(OutfileN);
	printf("Done!\n");
}

FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}

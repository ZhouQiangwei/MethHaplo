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

int main(int argc,char *argv[])
{
	const char* Help_String="Command Format :  splitmr [options] -i <methratiofile>\n"
        "\nMethyHaplo::splitmr \nversion: v1.0\n\n"
		"\n======================== PROGRAM OPTIONS ========================\n\n"
		"\t-i                    input file.  This file is the result of methyhaplo.\n"
        "\t-c [ALL/CG/CHG/CHH]   Meth context. [default: ALL]\n"
		"\t-m [int]              Meth coverage. [default: 2]\n"
        "\t-n [int]              Total coverage. [default: 6]\n"
		"\t-f|--Mf [double]      Methy level. [default: 0.2]\n"
		"\t-h|--help\n";
	
	char* methratioIN;
    double weight=0.2;
    int mCover = 2;
    int TotalCover=6;
    char* Mcontext;
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-i")   )
		{
			methratioIN=argv[++i];
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
            Mcontext=argv[++i];
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
	
	FILE* INfile = File_Open(methratioIN,"r");
	
	char Dummy[BATBUF];
	char chrom[100];
    char strand;char context[10];
    int meth=0;int cover=0;float ml=0;
    int pos=0;
    int nline=0;
    float hweight = 1.0-weight;
	while(fgets(Dummy,BATBUF,INfile)!=0)
	{
		if(Dummy[0]=='#') 
		{
			continue;
		}

		//chrM    1       -       CHH     2       1795    0.001114
		sscanf(Dummy,"%s%d%s%s%d%d%f",chrom,&pos,&strand,context,&meth,&cover,&ml);
        nline++;
        if(nline%1000000 == 0){
            fprintf(stderr, "[MethProcess] Processed %d meth loci", nline);
        }
        if(ml < weight || ml > hweight) continue;
        if(meth < mCover) continue;
        if(cover < TotalCover) continue;

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

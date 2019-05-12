#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <map>
#include <vector>
#include <math.h>
#include <algorithm>
#include <limits.h>
#define BATBUF 3000
#define CHROM_NAME_LEN 40
//#define cutoffNumber 4
FILE* File_Open(const char* File_Name,const char* Mode);
struct methyhaplo
{
	char chrom[CHROM_NAME_LEN];
	int coordinate1;
	int coordinate2;
	std::string s2t;
};
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
double caculate_weight(int a,int b,int c,int d,double& fc)
{
	int max=a;int second=b;
	
	if(a<b)
	{
		max=b;
		second=a;
	}
	
	if(c>max)
	{
		second=max;
		max=c;
	}
	else if(c>second && c<=max)
	{
		second=c;
	}
	
	if(d>max)
	{
		second=max;
		max=d;
	}
	else if(d>second && d<=max)
	{
		second=d;
	}
	if(second>0) fc = log((double)max/second);

	if(a+b+c+d == 0) return 0;
	return ((double)(max+second))/(a+b+c+d);

}

int main(int argc,char *argv[])
{
	printf("\nBatMeth2::ASM \nversion: v1.0\n\n");
	const char* Help_String="Command Format :  ASM [options] -i <haplotype result> -o <ASM out result>\n"
		"\n======================== PROGRAM OPTIONS ========================\n\n"
		"\t-o|--out              Output file prefix. Allele-specific methylation(ASM) output.\n"
		"\t-i                    input file.  This file is the result of methyhaplo.\n"
		"\t-b                    output file. Allelically methylated regions(AMRs) bed format file.\n"
		"\t-p|--pvale [double]   Pvalue cutoff. default: 0.01\n"
		"\t-c|--cutoff [int]     print the region that >= number of adjant Cs. default: 4\n"
		"\t-w|--weight [double]  (max+second_max)/all > double(example:0.9)\n"
		"\t-f [double]           log(max/second_max) < fc.(exaple:1)\n"
		"\t-ic [int]             min coverage, default 4.\n"
		"\t-xc [int]             max coverage, default 1000.\n"
		"\t-h|--help\n";
	
	char* methyhaploIN;
	char* methyhaploOUT;
	char* bedOUT;
	bool Print_Bed=false;
	std::string pvalue="0.01";
	int cutoffNumber=4;double weight=0.9;double FC=1.0; //1.0*INT_MAX;
        int minCover = 4;int maxCover=1000;
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-i")   )
		{
			methyhaploIN=argv[++i];
		}else if(!strcmp(argv[i], "-o") ||!strcmp(argv[i], "--out"))
		{
			methyhaploOUT=argv[++i];
		}else if(!strcmp(argv[i], "-b")   )
		{
			bedOUT=argv[++i];
			Print_Bed=true;
		}else if(!strcmp(argv[i], "-p") || !strcmp(argv[i], "--pvalue")   )
		{
			pvalue=argv[++i];
		}else if(!strcmp(argv[i], "-c") || !strcmp(argv[i], "--cutoff")   )
		{
			cutoffNumber=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-w") || !strcmp(argv[i], "--weight")   )
		{
			weight=atof(argv[++i]);
		}else if(!strcmp(argv[i], "-f") )
		{
			FC=atof(argv[++i]);
		}else if(!strcmp(argv[i], "-ic") )
                {
                        minCover=atoi(argv[++i]);
                }else if(!strcmp(argv[i], "-xc") )
                {
                        maxCover=atoi(argv[++i]);
                }
		else
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
	FILE* Outfile = File_Open(methyhaploOUT,"w");
	printf("Processing ...\n");
	FILE* Outfile_bed;
	if(Print_Bed)
		Outfile_bed= File_Open(bedOUT,"w");
	std::string awkcmd="awk '$(NF-1)<= ";
	awkcmd+=pvalue;
	awkcmd+=" ' ";
	awkcmd+= methyhaploIN ;
	strcat(methyhaploIN,".tmp");
	awkcmd += ">"; 
	awkcmd += methyhaploIN;
	system(awkcmd.c_str());
	
	FILE* INfile = File_Open(methyhaploIN,"r");
	
	char Dummy[BATBUF];char strTmp[BATBUF];bool isSNP=false;bool AllSNP=true;
	static std::vector<methyhaplo> METHYhaplo;
	int continueN=0;int read_first_line=0;
	while(fgets(Dummy,BATBUF,INfile)!=0)
	{
		if(Dummy[0]=='#') 
		{
			continue;
		}
		if(Dummy[0]=='*')
		{
			if(continueN>=cutoffNumber)
			{
				if(METHYhaplo.size() >= 1)
				{
					std::string meth_state;
					if(isSNP)
						if(AllSNP)
							meth_state="Hetero";
						else
							meth_state="hetero";
					else
						meth_state="Meth";
					fprintf(Outfile,"******span:\t%s\t%d\t%d\t%s:%d\t%d\n",METHYhaplo[0].chrom,METHYhaplo[0].coordinate1,METHYhaplo[METHYhaplo.size()-1].coordinate2,meth_state.c_str(),METHYhaplo[METHYhaplo.size()-1].coordinate2-METHYhaplo[0].coordinate1,continueN);
					if(Print_Bed) 
							fprintf(Outfile_bed,"%s\t%d\t%d\t%s:%d\t%d\n",METHYhaplo[0].chrom,METHYhaplo[0].coordinate1,METHYhaplo[METHYhaplo.size()-1].coordinate2,meth_state.c_str(),METHYhaplo[METHYhaplo.size()-1].coordinate2-METHYhaplo[0].coordinate1,continueN);
				}
				while(METHYhaplo.size()>0)
				{
					fprintf(Outfile,"%s",METHYhaplo[0].s2t.c_str());
					METHYhaplo.erase(METHYhaplo.begin());
				}
			}
			METHYhaplo.clear();
			isSNP=false;AllSNP=true;
			continueN=0;
			read_first_line=0;
			continue;
		}
		//chr10   71723   71765   3       3       2       3       1.000000        1.000000
		methyhaplo Meth,Meth_bak;
		Meth.s2t=Dummy;
		sscanf(Dummy,"%s%d%d%s",Meth.chrom,&Meth.coordinate1,&Meth.coordinate2,strTmp);
		
		isSNP= ! isunsigned(strTmp);
		int a,b,c,d;
		if(!isSNP)
		{
			AllSNP=false;
			sscanf(Dummy,"%*s%*d%*d%d%d%d%d",&a,&b,&c,&d);
			double fc=0;
			if(!( caculate_weight(a,b,c,d,fc) >= weight && fc <= FC ))
			{
				continue;
			}
			if(a+b+c+d < minCover || a+b+c+d > maxCover) continue;
		}
		
		if(read_first_line==0)
		{
			
			METHYhaplo.push_back(Meth);
			continueN=2;
			read_first_line=1;
		}else
		{
			if( !strcmp(Meth_bak.chrom,Meth.chrom) && Meth_bak.coordinate2 == Meth.coordinate1)
			{
					METHYhaplo.push_back(Meth);
					continueN++;
			}else
			{
				if(continueN>=cutoffNumber)
				{
					if(METHYhaplo.size() >= 1)
					{
						std::string meth_state;
						if(isSNP)
							if(AllSNP)
								meth_state="Hetero";
							else
								meth_state="hetero";
						else
							meth_state="Meth";
						fprintf(Outfile,"******span:\t%s\t%d\t%d\t%s:%d\t%d\n",METHYhaplo[0].chrom,METHYhaplo[0].coordinate1,METHYhaplo[METHYhaplo.size()-1].coordinate2,meth_state.c_str(),METHYhaplo[METHYhaplo.size()-1].coordinate2-METHYhaplo[0].coordinate1,continueN);
						if(Print_Bed) 
								fprintf(Outfile_bed,"%s\t%d\t%d\t%s:%d\t%d\n",METHYhaplo[0].chrom,METHYhaplo[0].coordinate1,METHYhaplo[METHYhaplo.size()-1].coordinate2,meth_state.c_str(),METHYhaplo[METHYhaplo.size()-1].coordinate2-METHYhaplo[0].coordinate1,continueN);
					}
					while(METHYhaplo.size()>0)
					{
						fprintf(Outfile,"%s",METHYhaplo[0].s2t.c_str());
						METHYhaplo.erase(METHYhaplo.begin());
					}
				}
				METHYhaplo.clear();
				METHYhaplo.push_back(Meth);
				isSNP=false;AllSNP=true;
				continueN=2;
			}
		}
		Meth_bak=Meth;
	}
	remove(methyhaploIN);
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

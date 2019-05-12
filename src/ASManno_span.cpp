//这个不对，整合基因以及上下游划分长度后，但是不知道到哪里是TSS、TES
#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include "limits.h"
#include <map>
#include <sstream>
#define CHROMSIZE 100
#define BATBUF 2000

struct Methy_Hash
{
	int *plusMethy;
	int *NegMethy;
	int *MethContext;
	int Index;
};
struct Methy_Gff
{
	long *C_Num;
	long AverC;
	int Index;
};
using namespace std;
bool Collision=false;
map <string,int> String_Hash;
map <string,int> Context_Hash;


//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_Progress(float Percentage);
void caculate(int start,int end,Methy_Hash MethyList,char Strand,Methy_Gff & methGff_List);
void caculateHeatmap(const char* type,int start,int end,Methy_Hash MethyList,char Strand,const char* id,FILE *methGffcg,char* chrom, int geneS,int geneE,FILE *methGffcg_matrix,bool printtitle=true);
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%% ++++++++]\n");//progress bar....
unsigned u=0;
unsigned Total_Reads;//total hits in the input file...
float binspan=0.025;
unsigned nLevel=ceil(1/(double)binspan)-1;//
void int2str(int &int_temp,string &string_temp)
{
        stringstream stream;
        stream<<int_temp;
        string_temp=stream.str();   
}
void unsigned2str(unsigned &int_temp,string &string_temp)
{
        stringstream stream;
        stream<<int_temp;
        string_temp=stream.str();   
}
void str2int(int &int_temp,string &string_temp)
{
	stringstream stream(string_temp);
	stream>>int_temp;
}
//}-----------------------------   GLOBAL VARIABLES  -------------------------------------------------
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..
long AverC=0;
int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	printf("\nMethyhaplo: ASManno v1.0\n");
	const char* Help_String="Command Format : ASManno [options] -o <OUT_PREFIX> -G GENOME -gff <GFF file>/-gtf <GTF file>/-b <bed file> -ap <asm plus file> -an <asm neg file> [-B] [-P]\n"
		"\nUsage:\n"
		"\t-o|--out         Output file prefix\n"
		"\t-G|--genome      Genome\n"
		"\t-ap|--asmplus    ASM plus file.\n"
		"\t-an|--asmneg     ASM neg file.\n"
		"\t-p|--pvale       Pvalue cutoff. default: 0.01\n"
		"\t-gtf|-gff        Gtf/gff file\n"
		"\t-b|--BED         Bed file, chrom start end (strand)\n"
		"\t-d|--distance    ASM distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000\n"
		"\t-s|--step        Gene body and their flanking sequences using an overlapping sliding window of 5% of the sequence length at a step of 2.5% of the sequence length. So default step: 0.025 (2.5%)\n"
		//"\t-S|--chromStep   Caculate the density of genes/TEs in chromsome using an overlapping sliding window of 100000bp at a step of 50000bp, must equal \"-s\" in Split.. default step: 50000(bp)\n"//
		"\t-h|--help";

	int Genome_CountX=0;
	char* Output_Name;
	//char* methInfileName;
	int InFileStart=0,InFileEnd=0;
	int NegInFileStart=0,NegInFileEnd=0;
	string Geno;
	char *InFile;
	bool InputGff=false;
	int GffInFileStart=0,GffInFileEnd=0;
	bool InputBed=false;
	int BedInFileStart=0,BedInFileEnd=0;
	int distance=2000;
	int distanceHeatmap=2000;
	int chromStep=50000;
	int binsStep=1000;
	bool Diff=false;
	bool PU=false;
	bool TSS=false;
	bool TTS=false;
	bool GENE=false;
	bool GTF=false;
	int beddistance=200;
	float pval_cut = 0.01;
	//
	if(argc<8)
	{
		printf("%s \n",Help_String);
		exit(0);
	}
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-o") ||!strcmp(argv[i], "--out")  )
		{
			Output_Name=argv[++i];
		}
		else if(!strcmp(argv[i], "-G") || !strcmp(argv[i], "--genome"))
		{
			Geno=argv[++i];
		}else if(!strcmp(argv[i], "-p") || !strcmp(argv[i], "--pvalue"))
                {
                        pval_cut=atof(argv[++i]);
                }
		else if(!strcmp(argv[i], "-s") || !strcmp(argv[i], "--step"))
		{
			binspan=atof(argv[++i]);
			nLevel=ceil(1/(double)binspan)-1;
			printf("binspn: %f, nLevel: %d\n", binspan, nLevel);
		}
		else if(!strcmp(argv[i], "--bins"))
		{
			binsStep=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-B") || !strcmp(argv[i], "--body"))
		{
			Diff=true;
		}
		else if(!strcmp(argv[i], "-P") || !strcmp(argv[i], "--promoter"))
		{
			PU=true;//promoter UP
		}
		else if( !strcmp(argv[i], "--TSS"))
		{
			TSS=true;//promoter UP
		}
		else if( !strcmp(argv[i], "--TTS"))
		{
			TTS=true;//promoter UP
		}
		else if( !strcmp(argv[i], "--GENE"))
                {
                        GENE=true;//promoter UP
                }
		else if(!strcmp(argv[i], "-S") || !strcmp(argv[i], "--chromStep"))
		{
			chromStep=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-d") || !strcmp(argv[i], "--distance"))
		{
			distance=atoi(argv[++i]);
			distanceHeatmap = distance;
		}
		else if(!strcmp(argv[i],"-ap") || !strcmp(argv[i],"--asmplus"))
		{
			InFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
		}else if(!strcmp(argv[i],"-an") || !strcmp(argv[i],"--asmneg"))
                {
                        NegInFileStart=++i;
                        while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
                        if(argv[i][0]=='-') {NegInFileEnd=--i;}else {NegInFileEnd=i ;}
		}else if(!strcmp(argv[i], "-gtf") || !strcmp(argv[i], "-gff"))
		{
			InputGff=true;
			if(!strcmp(argv[i], "-gtf")) GTF=true;
			GffInFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {GffInFileEnd=--i;}else {GffInFileEnd=i ;}
		}else if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--BED"))
		{
			InputBed=true;
			BedInFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {BedInFileEnd=--i;}else {BedInFileEnd=i ;}
		}
		else
		{
			printf("%s \n",Help_String);
			exit(0);
		}
	}
	if (argc >1) 
	{
		try
		{
			printf("Output prefix: %s\n", Output_Name);
			time(&Start_Time);

			string L=Geno;L+=".ann.location";

			FILE* Location_File=File_Open(L.c_str(),"r");
			
			printf("LocationFile: %s\n",L.c_str() );

			struct Offset_Record
			{
				char Genome[40];
				unsigned Offset;
			} Temp_OR; 

			while (fgets(Temp_OR.Genome,39,Location_File)!=0)//count genomes..
			{
				fgets(Temp_OR.Genome,39,Location_File);
				Genome_CountX++;	
			}
			rewind(Location_File);
			Offset_Record Genome_Offsets[Genome_CountX+1];

			int Genome_Count=0;
			while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<=Genome_CountX)
			{
				Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
				fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File);
				for(int i=0;i<40;i++) 
				{
					if (Genome_Offsets[Genome_Count].Genome[i] == '\n' ||Genome_Offsets[Genome_Count].Genome[i] == '\r')
					{ 
						Genome_Offsets[Genome_Count].Genome[i]=0;
						break;
					} 
				}
				Genome_Count++;	
			}
			Genome_Count--;

			Methy_Hash Methy_List[Genome_Count];
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				printf("%s, ",Genome_Offsets[i].Genome);
				String_Hash[Genome_Offsets[i].Genome]=i;
				//meth ini
				Methy_List[i].plusMethy = new int[Genome_Offsets[i+1].Offset];
				Methy_List[i].NegMethy = new int[Genome_Offsets[i+1].Offset];
				Methy_List[i].MethContext =new int[Genome_Offsets[i+1].Offset];
				Methy_List[i].Index=i;
			}
			printf("Loaded\n");
			string CG="CG",CHG="CHG",CHH="CHH";
			Context_Hash[CG.c_str()]=1;Context_Hash[CHG.c_str()]=2;Context_Hash[CHH.c_str()]=3;
			
			unsigned pos1,pos2;float pval=0,qval=0;int countMM=0,countMU=0,countUM=0,countUU=0;
			//int revG=0,revGA=0;
			//start to read batman hit file
			char Buf[BATBUF],Meth[BATBUF],Chrom[CHROMSIZE], noChrom[CHROMSIZE],effCT[BATBUF];
			printf("\nLoading plus ASM file...\n");
			FILE* INFILE;//=File_Open(methInfileName,"r");
		for(int f=InFileStart;f<=InFileEnd;f++)
		{
			INFILE=File_Open(argv[f],"r");
			int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
			fseek(INFILE, 0L, SEEK_END);off64_t File_Size=ftello64(INFILE);rewind(INFILE);
			Init_Progress();
			while (!feof(INFILE)) 
			{
				Total_Reads++;
				Progress++;
				if (Progress>Number_of_Tags) 
				{
					off64_t Current_Pos=ftello64(INFILE);
					off64_t Average_Length=Current_Pos/Total_Reads+1;//+1 avoids divide by zero..
					Number_of_Tags=(File_Size/Average_Length)/20;
					Progress=0;
					Show_Progress(Current_Pos*100/File_Size);
				}
				fgets(Meth,BATBUF,INFILE);
				if(Meth[0] == '*' || Meth[0] == '#') continue;
				sscanf(Meth,"%s%u%u%d%d%d%d%f%f",Chrom,&pos1,&pos2,&countMM,&countMU,&countUM,&countUU,&pval,&qval);
				if(pval > pval_cut) continue;
				int H=String_Hash[Chrom];
				pos1--;pos2--;//0-- for array start 0
				Methy_List[H].plusMethy[pos1]=1;
				Methy_List[H].plusMethy[pos2]=1;
				Methy_List[H].MethContext[pos1]=1;
				Methy_List[H].MethContext[pos2]=1;
			}//end read file while
			Done_Progress();
		} //end foreach meth file
			fclose(INFILE);
			//NEG ASM
			Total_Reads = 0;
			printf("\nLoading Neg ASM file...\n");
                        FILE* negINFILE;//=File_Open(methInfileName,"r");
        	        for(int f=InFileStart;f<=InFileEnd;f++)
                	{
	                        negINFILE=File_Open(argv[f],"r");
        	                int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
                	        fseek(negINFILE, 0L, SEEK_END);off64_t File_Size=ftello64(negINFILE);rewind(negINFILE);
                        	Init_Progress();
	                        while (!feof(negINFILE))
        	                {
                	                Total_Reads++;
                        	        Progress++;
                                	if (Progress>Number_of_Tags)
	                                {
        	                                off64_t Current_Pos=ftello64(negINFILE);
                	                        off64_t Average_Length=Current_Pos/Total_Reads+1;//+1 avoids divide by zero..
                        	                Number_of_Tags=(File_Size/Average_Length)/20;
                                	        Progress=0;
                                        	Show_Progress(Current_Pos*100/File_Size);
	                                }
        	                        fgets(Meth,BATBUF,negINFILE);
                	                if(Meth[0] == '*' || Meth[0] == '#') continue;
                        	        sscanf(Meth,"%s%u%u%d%d%d%d%f%f",Chrom,&pos1,&pos2,&countMM,&countMU,&countUM,&countUU,&pval,&qval);
                                	if(pval > pval_cut) continue;
	                                int H=String_Hash[Chrom];
        	                        pos1--;pos2--;//0-- for array start 0
                	                Methy_List[H].NegMethy[pos1]=1;
                        	        Methy_List[H].NegMethy[pos2]=1;
					Methy_List[H].MethContext[pos1]=1;
					Methy_List[H].MethContext[pos2]=1;
	                        }//end read file while
        	                Done_Progress();
	                } //end foreach meth file
                        fclose(negINFILE);

			printf("ASM file Loaded.\n");

			Methy_Gff methGff_List;
			{//0--UP 1--BODY 2--DOWN
				methGff_List.C_Num=new long[nLevel];
				for(int j=0;j<nLevel;j++)
					methGff_List.C_Num[j] = 0;
				methGff_List.AverC=0;
				methGff_List.Index=0;
			}
			string du[3];
			du[0]="UP";du[1]="BODY";du[2]="DOWN";
			int FileS=0,FileE=0;
			if(InputGff)
			{
				FileS=GffInFileStart;
				FileE=GffInFileEnd;
			}else if(InputBed)
			{
				FileS=BedInFileStart;
				FileE=BedInFileEnd;
			}else
			{
				fprintf(stderr,"\nNo input Gff/Bed file .\n");
				exit(0);
			}
			for(int f=FileS;f<=FileE;f++){
				AverC=0;
				printf("\nProcessing %d out of %d. InFile: %s\n", f-FileS+1,FileE-FileS+1, argv[f]);
					FILE* GFFINFILE=File_Open(argv[f],"r");
					bool bed4 = false; int showbed = 0;
					int filelen = strlen(argv[f]);
					if(filelen>3 && argv[f][filelen-4]=='b' && argv[f][filelen-3]=='e' && argv[f][filelen-2]=='d' && argv[f][filelen-1]=='4' ) 
						bed4 = true;
					char T[100];sprintf (T, "%s.AverMethy.%d.txt", Output_Name, f-FileS+1);
					//gene body && Promoter
					char G_cg[100];sprintf (G_cg, "%s.body.cg.%d.txt", Output_Name, f-FileS+1);
					char P_cg[100];sprintf (P_cg, "%s.Promoter.cg.%d.txt", Output_Name, f-FileS+1);
					char mode[10];
					strcpy(mode, "aw+");
					if(f == FileS) strcpy(mode, "w");
					//
					FILE* MethGffbodyOUTFILEcg;
					if(Diff)
					{
						MethGffbodyOUTFILEcg=File_Open(G_cg,mode);
					}
					FILE* MethGffpromoterOUTFILEcg;
					if(PU){
						MethGffpromoterOUTFILEcg=File_Open(P_cg,mode);
					}
					//heatmap
					//TSS
					char TSS_cg[100];sprintf (TSS_cg, "%s.TSS.cg.%d.txt", Output_Name, f-FileS+1);
					FILE* MethGff_TSS_OUTFILEcg;
					char TSS_cg_matrix[100];sprintf (TSS_cg_matrix, "%s.TSS.cg.%d.matrix", Output_Name, f-FileS+1);
                                        FILE* MethGff_TSS_OUTFILEcg_matrix;
					if(TSS)
					{
						MethGff_TSS_OUTFILEcg=File_Open(TSS_cg,mode);
                                                MethGff_TSS_OUTFILEcg_matrix=File_Open(TSS_cg_matrix,mode);
					}
					//TTS
					char TTS_cg[100];sprintf (TTS_cg, "%s.TTS.cg.%d.txt", Output_Name, f-FileS+1);
					FILE* MethGff_TTS_OUTFILEcg;
					if(TTS)
					{
						MethGff_TTS_OUTFILEcg=File_Open(TTS_cg,mode);
					}
                                        //GENE
                                        char GENE_cg[100];sprintf (GENE_cg, "%s.GENE.cg.%d.txt", Output_Name, f-FileS+1);
                                        FILE* MethGff_GENE_OUTFILEcg;
					if(GENE)
					{
                                                MethGff_GENE_OUTFILEcg=File_Open(GENE_cg,mode);
					}
					//---------------------------------------------------------------------------
					char N[100];sprintf (N, "%s.Methy.%d.txt", Output_Name, f-FileS+1);
					//-----------------------------------------------------------------------------------------------
					printf("\n------------------------------------------------------------------------------------------------\n");
					printf("\nASM Aver output file: %s\n",T);
					printf("ASM output file: %s\n",N);
					if(Diff) printf("Gff/Bed different analysis methylation file: %s.regions.context.%d.txt\n",Output_Name,f-FileS+1);
					printf("\n------------------------------------------------------------------------------------------------\n");
					//--------------------------------------------------------------------------------------------------

				int Progress=0;int Number_of_Tags=INITIAL_PROGRESS_READS;
				fseek(GFFINFILE, 0L, SEEK_END);off64_t File_Size=ftello64(GFFINFILE);rewind(GFFINFILE);//
				Init_Progress();Total_Reads=0;//
				char Gff[BATBUF],temp[BATBUF];
				string id="";
				char Symbol[BATBUF];
				unsigned start=0,end=0;
				char Strand = '.';
				int MLEN=1;
				while(!feof(GFFINFILE))
				{
					Total_Reads++;
					Progress++;
					if (Progress>Number_of_Tags) 
					{
						off64_t Current_Pos=ftello64(GFFINFILE);
						off64_t Average_Length=Current_Pos/Total_Reads+1;//+1 avoids divide by zero..
						Number_of_Tags=(File_Size/Average_Length)/20;
						Progress=0;
						Show_Progress(Current_Pos*100/File_Size);
					}
					fgets(Gff,BATBUF,GFFINFILE);
					//string id="";
					if(InputGff)
					{
						//sscanf(Gff,"%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%[^;]",Chrom,temp,temp,&start,&end,temp,&Strand,temp,Symbol);
						if(GTF) sscanf(Gff,"%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s\t%s",Chrom,temp,temp,&start,&end,temp,&Strand,temp,temp,Symbol);
			                        else sscanf(Gff,"%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%[^\n\t]",Chrom,temp,temp,&start,&end,temp,&Strand,temp,Symbol);
						id = Symbol;
					}else if(InputBed)
					{
						if(bed4){
							sscanf(Gff,"%s\t%u\t%u\t%c",Chrom,&start,&end,&Strand);
							if(showbed == 0) 
							{	
								fprintf(stderr, "bed format and defined strand %c\n", Strand);
								showbed = 1;
							}
						}else{
							sscanf(Gff,"%s\t%u\t%u",Chrom,&start,&end);
							Strand='.';
							if(showbed == 0){
								fprintf(stderr, "bed format without strand, if you want define strand, please use bed4 format.\n");
								showbed = 1;
							}
						}
					}
					if(!strcmp(noChrom, Chrom)) continue;
					if(end-start<beddistance) continue;
					if(!InputBed && end-start < beddistance) fprintf(stderr, "Waring: element %s is too shoter to caculate heatmap.\n", id.c_str());
					else if(end-start < beddistance) fprintf(stderr, "Waring: element %s:%d-%d is too shoter to caculate heatmap.\n", Chrom, start, end);
                        	        map<string, int>::iterator it= String_Hash.find(Chrom);
                	                if(it == String_Hash.end()) {
        	                                printf("%s not detected meth\n", Chrom);
	                                        strcpy(noChrom,Chrom);
                                        	continue;
                                	}
					int H=String_Hash[Chrom];
					float binspan_b = binspan;
					if(start>distance && end+distance <Genome_Offsets[H+1].Offset && end-start > nLevel ) //&& end-start > nLevel 
					{
						start--;end--;
						MLEN = end-start+2*distance;
						//body
						caculate(start-distance,end+distance,Methy_List[H],Strand,methGff_List);
						//
					}
					//if(Strand=='.') Strand='+';
					if(InputBed)
					{
						string temp="";
						unsigned2str(start,temp);
						id=Chrom;id+=("."+temp);
						unsigned2str(end,temp);
						id=id+"."+temp;
					}
				}
				
				
				FILE* MethGffOUTFILE=File_Open(N,"w");
				//CG
				fprintf(MethGffOUTFILE,"C");
				int Mstep=1;
				{
					Mstep = floor((double) MLEN*binspan);
					for(int j=0;j<nLevel;j++)
					{
						fprintf(MethGffOUTFILE,"\t%f",(double) methGff_List.C_Num[j]/Mstep);
					}
				}
				fprintf(MethGffOUTFILE,"\n");
				fclose(MethGffOUTFILE);
				Done_Progress();
				printf("\n%s\nC : %d\n",argv[f], AverC);
			}
			//delete 
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				delete[] Methy_List[i].plusMethy;
				delete[] Methy_List[i].NegMethy;
				
			}
			{//0--UP 1--BODY 2--DOWN
				delete[] methGff_List.C_Num;
			}
		}
		catch(char* Err)
		{
			printf(Err);
			exit(-1);
		}
		time(&End_Time);printf("\nTime Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}
	else
		printf("%s \n",Help_String);
}

//{----------------------------------- FILE HANDLING ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  File_Open
 *  Description:  Open a file:
 *  Mode - "w" create/Write text from top    - "wb" Write/Binary  -"w+" open text read/write            -"w+b" same as prev/Binary
 *         "r" Read text from top            - "rb" Read/Binary   -"r+" open text read/write/nocreate   -"r+b" same as prev/binary
 *       - "a" text append/write                                  -"a+" text open file read/append      -"a+b" open binary read/append
 *
 * =====================================================================================
 */
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

//}----------------------------------- FILE HANDLING ---------------------------------------------------------


void Show_Progress(float Percentage)
{
	if (Percentage >98) return;
	printf("+%.0f%\b\b\b",Percentage);
	fflush(stdout);
}

void caculate(int start,int end,Methy_Hash MethyList,char Strand,Methy_Gff & methGff_List){
        int step = floor((double(end - start))*binspan);// 5%  2.5% step
        int nbins=0;
        methGff_List.AverC=0;
        
            unsigned long countC=0;
            unsigned long countC_1=0;    
            for(int i=start;i<=end;i++)
            {
	                if(Strand=='+' || Strand=='.')
	                {
				if( MethyList.MethContext[i] ==1 ){
	                            countC+=MethyList.plusMethy[i];
	                            methGff_List.AverC +=MethyList.plusMethy[i];
	                            AverC +=MethyList.plusMethy[i];
				}
			}
			if(Strand=='-' || Strand=='.')
			{
				if( MethyList.MethContext[i] ==1 ){
	                            countC+=MethyList.NegMethy[i];
	                            methGff_List.AverC +=MethyList.NegMethy[i];
	                            AverC +=MethyList.NegMethy[i];
				}
			}
                if( nbins<=nLevel && (i-start) == ((nbins+1)*step)){
                	if(i==end) nbins=nLevel;
                    if(nbins<=nLevel &&nbins>0){
                    	if(Strand=='+' || Strand=='.')
                    	{
	                            methGff_List.C_Num[nbins-1] += (countC+countC_1);
                        }
                        else
                        {
	                            methGff_List.C_Num[nLevel-nbins] += (countC+countC_1);
                        }
                    }
                    countC_1=countC;
                    nbins++;
                    countC=0;
                    
                }
      }
}

//caculate gene heatmap //http://www.sciencedirect.com/science/article/pii/S0092867413002225  //Fig.6G
void caculateHeatmap(const char* type,int start,int end,Methy_Hash MethyList,char Strand,const char* id,FILE *methGffcg,char* chrom, int geneS, int geneE ,FILE *methGffcg_matrix, bool printtitle){
        int step = floor((double(end - start))*((double)binspan/2)); //5%  2.5% step
	unsigned nLevel=ceil(1/((double)binspan/2))-1;
//	int step = floor((double(end - start))*((double)binspan/1));
//	unsigned nLevel=ceil(1/((double)binspan/1))-1;
        int nbins=0;
        int* Smeth_c=new int[nLevel];
        
            int countC=0;
            int countC_1=0;

            for(int i=start;i<=end;i++)
            {
	                //context
	               if(Strand=='+' || Strand=='.')
	                {
				if( MethyList.MethContext[i] ==1 )
	                            countC+=MethyList.plusMethy[i];
			}
			else if(Strand=='-' )//|| Strand=='.'
			{
				if( MethyList.MethContext[i] ==1 )
	                            countC+=MethyList.NegMethy[i]; 
			}
                if( (nbins!=nLevel && (i-start) == ((nbins+1)*step)) ||  i==end){
                    if(nbins<=nLevel &&nbins>0){
                    	if(Strand=='+' || Strand=='.')
                    	{
					  Smeth_c[nbins-1]=countC+countC_1;
                        }
                        else
                        {
					  Smeth_c[nLevel-nbins]=countC+countC_1;
                        }
                    }
                    countC_1=countC;
                    nbins++;
                    countC=0;
                    if(i==end) nbins=nLevel;
                }
      }
        if(printtitle && (!strcmp(type,"GENE") || !strcmp(type,"TSS")))
	{
		if(!strcmp(type,"GENE")){
			fprintf(methGffcg,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
		}else{
        	        fprintf(methGffcg_matrix,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
        	        fprintf(methGffcg,"%s",id);
		}
	}else if(printtitle) {
		fprintf(methGffcg,"%s",id);
	}
        unsigned cut = ceil((double)nLevel/2);
        unsigned BeginNo_cg=0,BeginNo_chg=0,BeginNo_chh=0;
        //if(!strcmp(type,"TSS"))
        {
	        for(int i=0;i<nLevel;i++)
	        {
	        	if(BeginNo_cg==0 && !strcmp(type,"TSS"))
	        	{	
		        	if(i>=cut )
		       		{
		       		if(Smeth_c[i]>0) BeginNo_cg=i;
		       		else if(i==nLevel-1) BeginNo_cg=0;
		        	}
	        	}  
	        	//if(Smeth_c[i]>1) Smeth_c[i]=0;
			if(Smeth_c[i] < 0){
				if(!strcmp(type, "GENE")) fprintf(methGffcg,"\tnan");
				else if(!strcmp(type, "TSS")) {
					fprintf(methGffcg,"\t0");
					fprintf(methGffcg_matrix,"\tnan");
				}else fprintf(methGffcg,"\t0");
			}else fprintf(methGffcg,"\t%d",Smeth_c[i]);
	        }
        }
        
        if(!strcmp(type,"TTS"))
        {
	        for(int i=nLevel-1;i>=0;i--)
	        {
	        	if(BeginNo_cg==0)
	        	{	
		        	if(i< cut-1 )
		       	{
		       		if(Smeth_c[i]>0) BeginNo_cg=i;
		       		else if(i==0) BeginNo_cg=cut;
		        	}
	        	}  
	        }
        }
	if(strcmp(type,"GENE")!=0){
        	fprintf(methGffcg,"\t%d\n",BeginNo_cg);
		if(strcmp(type,"TSS")==0) {
			fprintf(methGffcg_matrix,"\n");
		}
	}
}

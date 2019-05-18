//g++ bam2md.cpp -o bam2md -m64 -I./samtools-0.1.18/ -L./samtools-0.1.18/ -lbam -lz
#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include "limits.h"
#include <map>
#include <algorithm>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>

#ifdef __cplusplus
extern "C"
{
#endif
        #include "./samtools-0.1.18/sam_header.h"
        #include "./samtools-0.1.18/sam.h"
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif
        bam_header_t *bam_header_dup(const bam_header_t *h0);
        samfile_t *samopen(const char *fn, const char *mode, const void *aux);
        void samclose(samfile_t *fp);
#ifdef __cplusplus
}
#endif


#define MAX_HITS_ALLOWED 1500
#define CHROMSIZE 100
#define BATBUF 50000
#define MAXTAG 500

#define QLIMIT_FLOAT 30.0f
//#define max(a,b) ( ((a)>(b)) ? (a):(b) )
//#define min(a,b) ( ((a)>(b)) ? (b):(a) )

const float QLIMIT=QLIMIT_FLOAT;

struct Mismatch
{
	char Base;
	int Pos;
};
struct Gene_Hash
{
	char* Genome;
	int Index;
};
struct Methy_Hash
{
	int *plusMethylated,*plusUnMethylated;//,*plusCover
	int *plusG,*plusA;
	int *NegMethylated,*NegUnMethylated;//,*NegCover
	int *NegG,*NegA;
	//int *MethContext;
	int Index;
};
struct Offset_Record
{
	char Genome[40];
	unsigned Offset;
} Temp_OR; 
typedef struct {
   Gene_Hash* Genome_List;
   Offset_Record* Genome_Offsets;
   char* Org_Genome;
   char* Marked_Genome;
   char* Marked_GenomeE;
   Methy_Hash Methy_List;
   FILE* OUTFILE;
   FILE* samINFILE;
   samfile_t* BamInFile;
   bam1_t *b;
   bam_header_t *header;
   int ThreadID;
   off64_t File_Size;
} ARGS;
bool RELESEM = false;
bool printheader = true;
using namespace std;
bool Collision=false;
map <string,int> String_Hash;
float ENTROPY_CUTOFF=0;
///g++ ./src/split.cpp -o ./src/split -lpthread
//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
int fprintf_time(FILE *stream, const char *format, ...);
int strSearch(char *str1,char *str2);
int Get_ED(std::string & S);
char* replace(char*src, char*sub, char*dst);
void MinandSec(unsigned a[MAX_HITS_ALLOWED][10], int left, int right, int&min, int&second,int & loci);
//void MaxandSec(int a[], int left, int right, int & max, int & second);
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_Progress(float Percentage);
void fetch(const char *str, char c1, char c2, char *buf);
string&   replace_all(string&   str,const   string&   old_value,const   string&   new_value);
std::string m_replace(std::string str,std::string pattern,std::string dstPattern,int count=-1);
float Pr(float Q);
void Build_Pow10();
void initMem(char* temp,int size,char* temp2);
void *Process_read(void *arg);
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%  ++++++++]\n");//progress bar....
#define random(x) (rand()%x)
inline unsigned char Hash(char* S);
int Get_String_Length(FILE* INFILE);
string remove_soft_split(string & cigar,int & Read_Len,int & pos);

void ReverseC_Context(char* Dest,const char* seq,int & stringlength);
//---------------------------------------------------------------------------------------------------------------
unsigned Total_Reads_all;
//----mapping count
unsigned Tot_Unique_Org=0;//total unique hits obtained
unsigned ALL_MAP_Org=0;
unsigned Tot_Unique_Remdup=0;//total unique hits obtained after removing dups...
unsigned ALL_Map_Remdup=0;
float UPPER_MAX_MISMATCH=0.06;
bool REMOVE_DUP=false; //true; //true to removeDup, false will not remove PCR-dup
unsigned Mismatch_Qual[255][255][255]; //[readLength][255][255]
int QualCut=10;
const int POWLIMIT=300;
float POW10[POWLIMIT];
int QUALITYCONVERSIONFACTOR=33;
//----------
bool Methratio=false;
//-------meth count
unsigned long non_met_CG=0;
unsigned long met_CG=0;
unsigned long non_met_CHG=0;
unsigned long met_CHG=0;
unsigned long non_met_CHH=0;
unsigned long met_CHH=0;
bool bamformat=false;

int Sam=1;//1 true 0 false
unsigned Number_of_Tags = 0;
//}-----------------------------   GLOBAL VARIABLES  -------------------------------------------------
char Char2Comp[255];
const int INITIAL_PROGRESS_READS =10000;//progress bar initial count..
float calc_Entropy(string readString, int L);
void Print_Mismatch_Quality(FILE* OUTFILE_MM, int L);
bool SamSeqBeforeBS=false;
long longestChr = 100000;
string processingchr = "NULL";
string newchr = "NULL";
//
int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	
	const char* Help_String="Command Format :  bam2md [options] -g GENOME  -i/-b <SamfileSorted/BamfileSorted> \n"
		"\nUsage:\n"
		"\t-g|--genome           Genome\n"
		"\t-i|--input            Sam format file, sorted by chrom.\n"
		"\t-b|--binput           Bam format file, sorted by chrom.\n"
		//"\t-p|--threads          the number of threads.\n"
		"\t-n|--Nmismatch [float]  Number of mismatches, default 0.06 percentage of read length. [0-1]\n"
		"\t-Q [int]              caculate the methratio while read QulityScore >= Q. default:10\n"
		"\t-r|--remove_dup       REMOVE_DUP, default:false\n"
		"\t-f|--sam [outfile]    f for sam format outfile contain methState. default: sam format.\n"
		"\t--sam-seq-beforeBS    Converting BS read to the genome sequences.\n"
		"\t-h|--help";

	Char2Comp['A']=Char2Comp['a']='T';
	Char2Comp['C']=Char2Comp['c']='G';
	Char2Comp['G']=Char2Comp['g']='C';
	Char2Comp['T']=Char2Comp['t']='A';
	Char2Comp['N']=Char2Comp['n']='N';
	int Genome_CountX=0;
	char Output_Name[100];
	strcpy(Output_Name, "None");
	string Prefix="None";
	string methOutfileName;
	string methWigOutfileName;
	string Geno;
//	int Current_Option=0;
	int InFileStart=0,InFileEnd=0;
	char *InFile;
	string binsOutfileName="";
	
//	int Par_Count=0;
	std::string CMD;

	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-f") ||!strcmp(argv[i], "--sam")  )
		{
			strcpy(Output_Name, argv[++i]);
			Sam=1;
		}
		else if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--genome"))
		{
			Geno=argv[++i];
		}else if(!strcmp(argv[i], "-Q"))
		{
			QualCut=atoi(argv[++i]);
        }
		else if(!strcmp(argv[i], "--sam-seq-beforeBS"))
		{
			SamSeqBeforeBS=true;
		}
		else if(!strcmp(argv[i], "-r") || !strcmp(argv[i], "--remove_dup"))
		{
			REMOVE_DUP=true;
		}else if(!strcmp(argv[i], "-n") || !strcmp(argv[i], "--Nmismatch"))
		{    
			UPPER_MAX_MISMATCH=atof(argv[++i]);
			if(UPPER_MAX_MISMATCH>1){
				fprintf(stderr, "\nError defined mismatch paramater, should be 0-1.\n");
			}
		}  
		else if(!strcmp(argv[i], "-i") || !strcmp(argv[i], "--input"))
		{
			InFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-')
			{
				i++;
				continue;
			}
			if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
		}else if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--binput"))
		{
			InFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-')
			{
				i++;
				continue;
			}
			if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
			bamformat=true;
		}
		else if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
			printf("\n%s\n",Help_String);
                        exit(0);
		}
		else
		{
			printf("\nError: %s\n\n%s \n",argv[i],Help_String);
			exit(0);
		}
		
	}
	for(int i = 0; i < argc; i++) {CMD.append(argv[i]); CMD.append(" ");}
	
	if(argc<=3) printf("%s \n",Help_String);
	if (argc >3  && InFileStart) 
	{
		printf("\nBatMeth2::bam2md v2.0\n");

      try
		{
			time(&Start_Time);
			Build_Pow10();
			
			string G=Geno;G+=".bin";
			string L=Geno;L+=".ann.location";

			FILE* BINFILE=File_Open(G.c_str(),"r");
			FILE* Location_File=File_Open(L.c_str(),"r");
			ARGS args;
			fseek(BINFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(BINFILE);rewind(BINFILE);//load original genome..
			args.Org_Genome=new char[Genome_Size];if(!args.Org_Genome) throw("Insufficient memory to load genome..\n"); 
			if(!fread(args.Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");
			if(REMOVE_DUP){
				args.Marked_Genome=new char[Genome_Size+1];if(!args.Marked_Genome) throw("Insufficient memory to Mark genome..\n"); 
				args.Marked_GenomeE=new char[Genome_Size+1];if(!args.Marked_GenomeE) throw("Insufficient memory to Mark genome..\n"); 
			}
			printf("Load Genome.. %s\n", Geno.c_str());
			while (fgets(Temp_OR.Genome,39,Location_File)!=0)//count genomes..
			{
				fgets(Temp_OR.Genome,39,Location_File);
				Genome_CountX++;	
			}
			rewind(Location_File);
		
			args.Genome_Offsets = new Offset_Record[Genome_CountX+1];
			int Genome_Count=0;
			while (fgets(args.Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<=Genome_CountX)
			{
				args.Genome_Offsets[Genome_Count].Offset=atoi(args.Genome_Offsets[Genome_Count].Genome);
				if(longestChr < args.Genome_Offsets[Genome_Count].Offset) longestChr = atoi(args.Genome_Offsets[Genome_Count].Genome);
				fgets(args.Genome_Offsets[Genome_Count].Genome,39,Location_File);
				for(int i=0;i<40;i++) 
				{
					if (args.Genome_Offsets[Genome_Count].Genome[i] == '\n' ||args.Genome_Offsets[Genome_Count].Genome[i] == '\r')
					{ 
						args.Genome_Offsets[Genome_Count].Genome[i]=0;
						break;
					} 
				}
				Genome_Count++;	
			}
			Genome_Count--;
			args.OUTFILE = NULL;
			assert(longestChr>0);
			printf("\nLongest chr: %d\n",longestChr);
			if(Sam && strcmp(Output_Name,"None") ) args.OUTFILE=File_Open(Output_Name,"w");

			char* Split_Point=args.Org_Genome;//split and write...
			args.Genome_List = new Gene_Hash[Genome_Count];
			//args.Methy_List = new Methy_Hash;
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				String_Hash[args.Genome_Offsets[i].Genome]=i;
				//unsigned char H=Hash(Genome_Offsets[i].Genome);
				args.Genome_List[i].Genome=(Split_Point+=args.Genome_Offsets[i].Offset);
				args.Genome_List[i].Index=i;
				//meth ini
			}
			
			fclose(BINFILE);
			fclose(Location_File);

			////read file
			
			for(int f=InFileStart;f<=InFileEnd;f++)
			{
				printf("\nProcessing %d out of %d. File: %s, %d\n\n", f-InFileStart+1,InFileEnd-InFileStart+1, argv[f], bamformat);
				//fseek(args.INFILE, 0L, SEEK_END);args.File_Size=ftello64(args.INFILE);rewind(args.INFILE);
				char s2t[BATBUF];
				//if(args.OUTFILE!=NULL)
 				{ // && printheader
					//samfile_t *bamin = 0;
					args.BamInFile = 0;
					args.header;
                    if(bamformat)
                    {
                            if ((args.BamInFile = samopen(argv[f], "rb", 0)) == 0) {
                                    fprintf(stderr, "fail to open \"%s\" for reading.\n", argv[f]);
                            }
                            args.header=bam_header_dup((const bam_header_t*)args.BamInFile->header);
                            if(InFileStart && args.OUTFILE!=NULL) fprintf(args.OUTFILE,"%s",args.header->text);
                    }
                    else if(f==InFileStart){
                    	args.samINFILE=File_Open(argv[f],"r");
			       if(args.OUTFILE!=NULL) {
						while (fgets(s2t,BATBUF,args.samINFILE)!=0 ){
				                	if(s2t[0]=='@') 
				                	{    
			                            	s2t[BATBUF]='\0';s2t[BATBUF-1]='\n';
                        			    	if(Sam){
                        			    	    fprintf(args.OUTFILE,"%s",s2t);
                        				    }    
			                    	        continue;
                					}else
								break;
						}
						rewind(args.samINFILE);
					}
				}
				}
				//nothreads
				Process_read(&args);
				Done_Progress();
				if(!bamformat) fclose(args.samINFILE);
            	if(bamformat) 
            	{
                	bam_header_destroy(args.header);
                	samclose(args.BamInFile);
            	}
			}

			printf("genome process done!\n");
			//if(args.OUTFILE!=NULL) fclose(args.OUTFILE);
			if(Sam && strcmp(Output_Name,"None") ) fclose(args.OUTFILE);
			
//delete
			printf("Done and release memory!\n");
			if(RELESEM){
	                        delete [] args.Genome_List;
	                        delete [] args.Genome_Offsets;
                        	delete [] args.Org_Genome;
			}
		}
		catch(char* Err)
		{
			printf(Err);
			fprintf(stderr, "\nError cigar\n");
			exit(-1);
		}

		time(&End_Time);printf("Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	
	}
	
}

int Read_Tag(FILE *INFILE,char s2t[],string hits[],int& cntHit)
{
	flockfile(INFILE);
	char Buf[BATBUF];
	if (!feof(INFILE))
	{
		Total_Reads_all++;
		fgets(s2t,BATBUF,INFILE);//read description..
		cntHit=0;
		fgets(Buf,BATBUF,INFILE);
		while(!feof(INFILE) && Buf[0]!='@')//while in the record.. 
		{
			hits[cntHit++]=Buf; //read hit info..
			fgets(Buf,BATBUF,INFILE);
			if(cntHit>=MAX_HITS_ALLOWED) 
			{
				cntHit=0;
		            while(!feof(INFILE) && Buf[0]!='@')
		            {
		                    fgets(Buf,BATBUF,INFILE);
		            }
				break;
			}
			assert(MAX_HITS_ALLOWED > cntHit);
		}
		funlockfile(INFILE);
		return 1;
	}
	else 
	{
		funlockfile(INFILE);
		return 0;
	}
}

char* process_cigar(const char* cig,int Read_Len)
{
	char temp[8];
	char* cigar_rm = new char[1000]();*cigar_rm=0;unsigned n=0;
	char* buffer_cigar=new char[100]();*buffer_cigar=0;
	while(*cig!='\0')
	{
		if(*cig>='0' && *cig<='9')
		{
			temp[n]=*cig;
			cig++;n++;
		}else if(*cig=='S')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				sprintf(buffer_cigar,"%dS",length);
				strcat(cigar_rm, buffer_cigar);
			}
            cig++;n=0;
		}else if(*cig=='M')
		{
			temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				sprintf(buffer_cigar,"%dM",length);
				strcat(cigar_rm, buffer_cigar);
			}
			cig++;n=0;
			char buf[1024]="\0";
			sprintf( buf , "%dM",Read_Len);
			char buf_tmp[1024]="\0";
			sprintf( buf_tmp, "%dM",length);
			if( !strcmp(buf,  buf_tmp ) )
			{
				if(!strcmp(buf, cig))
					break;
				else
				{
					strcpy(cigar_rm, buffer_cigar);
				}
				
			}
		}else if(*cig=='I')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1I");
				else if(length==2) sprintf(buffer_cigar,"2I");
				else sprintf(buffer_cigar,"%dI",length);
				strcat(cigar_rm, buffer_cigar);
			}
			cig++;n=0;
		}else if(*cig=='D')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1D");
				else if(length==2) sprintf(buffer_cigar,"2D");
				else sprintf(buffer_cigar,"%dD",length);
				strcat(cigar_rm, buffer_cigar);
			}
			cig++;n=0;
		}else
		{
			printf(" --%d%c-- ",atoi(temp),*cig);
    		continue;//break;
		}
	}
	delete []cigar_rm;
	delete []buffer_cigar;
	return cigar_rm;
}
string getstring(char* seq, int l, int len){
	char tmp[10];
	for(int i=0; i<len; i++){
//printf("\ns %d %d %c\n", l, i, seq[l+i]);
		tmp[i] = seq[l+i];
	}
	if(len>0) tmp[len]='\0';
	return tmp;
}

int conutMismatch(string CIGr, int chrLen, char* Genome_seq, string readString, int pos, int hitType)
{
	char temp[5];unsigned lens=0;int Glens=0;int RLens=0;
	unsigned n=0;
	int Nmismatch=0;  //chrLen; //((ARGS *)arg)->Genome_Offsets[Hash_Index+1].Offset
	const char* cigr=CIGr.c_str();
	while(*cigr!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
	{
		if(*cigr>='0' && *cigr<='9')
		{
			temp[n]=*cigr;
			cigr++;n++;
		}else if(*cigr=='S')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			lens+=length;
	        	RLens+=length;
        	    	cigr++;n=0;
		}else if(*cigr=='M')
		{
			temp[n]='\0';int length=atoi(temp);
			for(int k=lens,r=RLens,g=Glens;k<length+lens;r++,k++,g++)
			{
				if (pos+g-1 >= chrLen) break;

				char genome_Char = toupper(Genome_seq[pos+g-1]);
				if (hitType==1 || hitType==3) {
					if (readString[k] != genome_Char && !(readString[k]=='T' && genome_Char=='C')) 
					{
						Nmismatch++;
					}
				}
				else if (hitType==2 || hitType==4) {
					if (readString[k] != genome_Char && !(readString[k]=='A' && genome_Char=='G')) {
						
						Nmismatch++;
					}
				}
			}
			cigr++;n=0;lens+=length;Glens+=length;RLens+=length;
		}else if(*cigr=='I')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			lens+=length;
			RLens+=length;
			cigr++;n=0;
		}else if(*cigr=='D')
		{
			int i;temp[n]='\0';unsigned int length=atoi(temp);
			Glens+=length;RLens+=length;
			cigr++;n=0;
		}else
		{
			break;
		}
	}
	return Nmismatch;
}

int processbamread(const bam_header_t *header, const bam1_t *b, char* Dummy,int &Flag,char* Chrom,int &pos,int &mapQuality,char* CIG,char* Chrom_P,int &pos_P,int &Insert_Size,char* forReadString,char* forQuality, int &hitType)
{
        uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
        int i;
        const bam1_core_t *c = &b->core;
	strcpy(Dummy,  bam1_qname(b));
	Flag = c->flag;

	if( (Flag & 0x100) || (Flag & 0x200) || (Flag & 0x400) || (Flag & 0x800) || (Flag & 0x4))
        	return -1;
	if( !(Flag & 0x1) )
        {
        	if(Flag==0)
                	hitType=1;
                else if(Flag==16)
                	hitType=4;
        }else if( !(Flag & 0x10)  && (Flag & 0x40) )
        	hitType=1;
       	else if( !(Flag & 0x10)  && (Flag & 0x80) )
        	hitType=2;
        else if( (Flag & 0x10)  && (Flag & 0x80) )
        	hitType=3;
        else if( (Flag & 0x10)  && (Flag & 0x40) )
        	hitType=4;
        if(hitType==0) return -1;


        if (c->tid < 0) return -1;
        else {
                if (header) strcpy(Chrom, header->target_name[c->tid]); 
                else sprintf(Chrom, "%d", c->tid); 
        }
	pos = c->pos + 1;
	mapQuality = c->qual;

	//define
	if(mapQuality < QualCut || Flag==4 || (int)pos <= 0 ) return -1;
	char strtemp[256];int j=0;
        if (c->n_cigar == 0) return -1;
        else {
                for (i = 0; i < c->n_cigar; ++i) {
                    	sprintf(strtemp, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP=X"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
			for(int l= 0; l<strlen(strtemp); ++l,j++) CIG[j] = strtemp[l];
                }
        }
	CIG[j] = '\0';
        if (c->mtid < 0) Chrom_P[0] =  '*';
        else if (c->mtid == c->tid) sprintf(Chrom_P, "="); 
        else {
                if (header) strcpy(Chrom_P, header->target_name[c->mtid]); 
                else sprintf(Chrom_P, "%d", c->mtid); 
        }
	if(strcmp(Chrom_P, "*") == 0) {
        	pos_P = 0;
                Insert_Size = 0;
        }else{
		pos_P = c->mpos + 1;
		Insert_Size=c->isize;
	}
        if (c->l_qseq) {
                for (i = 0; i < c->l_qseq; ++i) forReadString[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
		forReadString[i] = '\0';
                if (t[0] == 0xff) forQuality[0] =  '*';
                else for (i = 0; i < c->l_qseq; ++i) forQuality[i] = (char)(t[i] + 33); 
		if(i!=0) forQuality[i] = '\0';
        } else {forReadString[0] = '*'; forQuality[0] =  '*';}
	
	return 1;

        s = bam1_aux(b);
	char read[100];
        while (s < b->data + b->data_len) {
                uint8_t type, key[2];
                key[0] = s[0]; key[1] = s[1];
                s += 2; type = *s; ++s;
                sprintf(read, "\t%s:", (char*)key); 
                if (type == 'A') { sprintf(read, "A:%c", *s); ++s; }
                else if (type == 'C') { sprintf(read, "i:%d", *s); ++s; }
                else if (type == 'c') { sprintf(read, "i:%d", *(int8_t*)s); ++s; }
                else if (type == 'S') { sprintf(read, "i:%d", *(uint16_t*)s); s += 2; }
                else if (type == 's') { sprintf(read, "i:%d", *(int16_t*)s); s += 2; }
                else if (type == 'I') { sprintf(read, "i:%d", *(uint32_t*)s); s += 4; }
                else if (type == 'i') { sprintf(read, "i:%d", *(int32_t*)s); s += 4; }
                else if (type == 'f') { sprintf(read, "f:%g", *(float*)s); s += 4; }
                else if (type == 'd') { sprintf(read, "d:%lg", *(double*)s); s += 8; }
                else if (type == 'Z' || type == 'H') { sprintf(read, "%c:", type); while (*s) sprintf(read, "%c", *s++); ++s; }
                else if (type == 'B') {
                        uint8_t sub_type = *(s++);
                        int32_t n;
                        memcpy(&n, s, 4);
                        s += 4; // no point to the start of the array
                    	sprintf(read, "%c:%c", type, sub_type);
                        for (i = 0; i < n; ++i) {
                                sprintf(read,",");
                                if ('c' == sub_type || 'c' == sub_type) { sprintf(read, "%d", *(int8_t*)s); ++s; }
                                else if ('C' == sub_type) { sprintf(read, "%d", *(uint8_t*)s); ++s; }
                                else if ('s' == sub_type) { sprintf(read, "%d", *(int16_t*)s); s += 2; }
                                else if ('S' == sub_type) { sprintf(read, "%d", *(uint16_t*)s); s += 2; }
                                else if ('i' == sub_type) { sprintf(read, "%d", *(int32_t*)s); s += 4; }
                                else if ('I' == sub_type) { sprintf(read, "%d", *(uint32_t*)s); s += 4; }
                                else if ('f' == sub_type) { sprintf(read, "%g", *(float*)s); s += 4; }
                        }
                }
        }
        
}

void *Process_read(void *arg)
{
	unsigned Total_Reads=0;
	int Progress=0;Number_of_Tags=INITIAL_PROGRESS_READS;
	Init_Progress();
	int mismatch=0;int pos=0;int Top_Penalty=0;int mapQuality=0;int Flag=-1;
	string readString="";
	int hitType=0;
        int fileprocess = 0;
	
	string hits[MAX_HITS_ALLOWED];
	char Comp_String[MAXTAG];for (int i=1;i<MAXTAG;Comp_String[i++]=0);
	//start to read batman hit file
	char *s2t = (char*) malloc(600);
	char read_Methyl_Info[600];char rawReadBeforeBS[600];char temp[5];
	char Dummy[BATBUF],forReadString[BATBUF],Chrom[CHROMSIZE];
	char Chrom_P[CHROMSIZE];int pos_P=0;int Insert_Size=0;int Qsingle=0; //Paired-end reads
	string CIGr;char CIG[BATBUF];
	char forQuality[BATBUF],rcQuality[BATBUF],Quality[BATBUF];
	int r=1;char* samaddress;
    struct timeval now;
    struct timespec outtime;
	bam1_t *b = bam_init1();
	while( (!bamformat && (samaddress = fgets(s2t,BATBUF,((ARGS *)arg)->samINFILE))!=NULL) || (bamformat && r>0 ))
	{
		if(r < -1) {
			fprintf(stderr, "\ntruncated file.\n");
		}
		hitType = 0;
		Total_Reads++;
                Progress++;
                fileprocess++;

		if(bamformat) 
		{
			(r = samread(( (ARGS *)arg)->BamInFile, b));
			//bam_tostring(((ARGS *)arg)->header , b, s2t);
			int ct = processbamread(((ARGS *)arg)->header, b, Dummy,Flag,Chrom,pos,mapQuality,CIG,Chrom_P,pos_P,Insert_Size,forReadString,forQuality, hitType);
			if(ct == -1) continue;
		}
		
		if ( fileprocess>=1000000  ) {
                        fprintf_time(stderr, "Processed %d reads.\n", Total_Reads);
                        fileprocess = 0;
        	}

		if(s2t[0]=='@') 
		{
			continue;
		}
                printheader = false;
		if(!bamformat)
			sscanf(s2t,"%s%d%s%d%d%s%s%d%d%s%s",Dummy,&Flag,Chrom,&pos,&mapQuality,CIG,Chrom_P,&pos_P,&Insert_Size,forReadString,forQuality);
		map<string, int>::iterator iter;
		int H = -1;
		if(strcmp(newchr.c_str(), Chrom) != 0 ) newchr = Chrom;
		
		if(strcmp(processingchr.c_str(), "NULL") == 0 ) processingchr = Chrom;
		if(strcmp(processingchr.c_str(), Chrom) != 0) {
            fprintf(stderr, "Print output of %s\n", processingchr.c_str());
            iter = String_Hash.find(processingchr.c_str());
			if(iter != String_Hash.end()){
				H = iter->second;
			}
	        processingchr = Chrom;
		}

		if(!bamformat){
 			if(mapQuality < QualCut || Flag==4 || (int)pos <= 0 ) continue;
			if(strcmp(Chrom_P, "*") == 0) {
				pos_P = 0;
				Insert_Size = 0;
			}
		}

		readString=forReadString;
		int Read_Len=readString.length();
		CIGr=CIG;
		//for(;forReadString[Read_Len]!=0 && forReadString[Read_Len]!='\n' && forReadString[Read_Len]!='\r';Read_Len++);
	    	iter = String_Hash.find(processingchr.c_str());
		H = -1;
		if(iter != String_Hash.end()){
			H = iter->second;
		}else continue;

		if(Flag==-1) printf("\n%s\n", Dummy);
		assert(Flag!=-1);
		if(!bamformat){
			if( (Flag & 0x100) || (Flag & 0x200) || (Flag & 0x400) || (Flag & 0x800) || (Flag & 0x4))
        	    		continue;
			if( !(Flag & 0x1) )
			{
				if(Flag==0)
					hitType=1;
				else if(Flag==16)
					hitType=4;
			}else if( !(Flag & 0x10)  && (Flag & 0x40) )
				hitType=1;
			else if( !(Flag & 0x10)  && (Flag & 0x80) )
				hitType=2;
			else if( (Flag & 0x10)  && (Flag & 0x80) )
				hitType=3;
			else if( (Flag & 0x10)  && (Flag & 0x40) )
				hitType=4;
			if(hitType==0) continue;
		}
		int Nmismatch=conutMismatch(CIGr, ((ARGS *)arg)->Genome_Offsets[H+1].Offset, ((ARGS *)arg)->Genome_List[H].Genome, readString, pos, hitType);
		if(Nmismatch > 0.5 + UPPER_MAX_MISMATCH * strlen(readString.c_str())) continue;
		//char* Genome;
		//
		//Genome=((ARGS *)arg)->Genome_List[H].Genome;//load current genome..
		
		unsigned G_Skip=((ARGS *)arg)->Genome_List[H].Genome-((ARGS *)arg)->Org_Genome;
            int Flag_rm=0;
		if(hitType == 1 || hitType == 3 ) Flag_rm=2; else Flag_rm=4;
		char Mark = '0';char MarkE='0';
		if(REMOVE_DUP){
            Mark=((ARGS *)arg)->Marked_Genome[pos+G_Skip];
            MarkE=((ARGS *)arg)->Marked_GenomeE[pos+G_Skip+readString.size()];
        }
            if( !REMOVE_DUP || (!Mark || !(Mark & Flag_rm)) || (!MarkE || !(MarkE & Flag_rm)) )
		{
			int Hash_Index=((ARGS *)arg)->Genome_List[H].Index;//load current genome..
			strcpy(rawReadBeforeBS,readString.c_str());
			unsigned lens=0;int Glens=0;int RLens=0;
			unsigned n=0;bool CONTINUE=false;
			const char* cigr=CIGr.c_str();
			while(*cigr!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
			{
				if(*cigr>='0' && *cigr<='9')
				{
					temp[n]=*cigr;
					cigr++;n++;
				}else if(*cigr=='S')
				{
					int i;temp[n]='\0';int length=atoi(temp);
					for(i=RLens;i<RLens+length;i++)
					{
						read_Methyl_Info[i] = 'S';
					}
					lens+=length;
	                RLens+=length;
	                cigr++;n=0;
				}else if(*cigr=='M')
				{
					temp[n]='\0';int length=atoi(temp);
					for(int k=lens,r=RLens,g=Glens;k<length+lens;r++,k++,g++)
					{
						read_Methyl_Info[r] = '=';
						if (pos+g-1 >= ((ARGS *)arg)->Genome_Offsets[Hash_Index+1].Offset) break;

						char genome_Char = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g-1]);//
						char genome_CharFor1 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g+1-1]);
						char genome_CharFor2 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g+2-1]);
						char genome_CharBac1,genome_CharBac2;
						if(pos+g-1 > 2)
						{
							genome_CharBac1 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g-1-1]);
							genome_CharBac2 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g-2-1]);
						}
						if (hitType==1 || hitType==3) {
							//pthread_mutex_lock(&meth_counter_mutex13);
							if (readString[k]=='C' && genome_Char=='C')
							{
								read_Methyl_Info[r] = 'M';
								if(Methratio ) ((ARGS *)arg)->Methy_List.plusMethylated[pos+g-1]++;
//if(pos+g-1 < 6) printf("\n%d %c\n", pos+g-1, Genome[pos+g-1]);
							//	if(hitType==1) {
									if(genome_CharFor1=='G') 
									{
										read_Methyl_Info[r] = 'Z';met_CG++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}//Z methylated C in CpG context
									else if(genome_CharFor1!='G' && genome_CharFor2!='G')
									{
										read_Methyl_Info[r] = 'H';met_CHH++;
										//if(Methratio && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}//H methylated C in CHH context
									else if(genome_CharFor1!='G' && genome_CharFor2=='G')
									{
										read_Methyl_Info[r] = 'X';met_CHG++;
										//if(Methratio && hitsToken[ind][4] > QualCut  )  ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}//X methylated C in CHG context
							}
							else if (readString[k]=='T' && genome_Char=='C')
							{
								read_Methyl_Info[r] = 'U';
								rawReadBeforeBS[k] = 'C';
								if(Methratio ) ((ARGS *)arg)->Methy_List.plusUnMethylated[pos+g-1]++;
							//	if(hitType==1) {
									if(genome_CharFor1=='G')
									{
										read_Methyl_Info[r] = 'z';non_met_CG++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}//z unmethylated C in CpG context
									else if(genome_CharFor1!='G' && genome_CharFor2!='G') 
									{
										read_Methyl_Info[r] = 'h';non_met_CHH++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}//h unmethylated C in CHH context
									else if(genome_CharFor1!='G' && genome_CharFor2=='G') 
									{
										read_Methyl_Info[r] = 'x';non_met_CHG++;
										//if(Methratio  && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}//x unmethylated C in CHG context
								
							}
							else if (readString[k] != genome_Char) 
							{
								read_Methyl_Info[r] = readString[k]; // genome_Char; for hypol
							}
							if(Methratio)
							{
								if(readString[k]=='G') ((ARGS *)arg)->Methy_List.plusG[pos+g-1]++;
								else if(readString[k]=='A') ((ARGS *)arg)->Methy_List.plusA[pos+g-1]++;
								//Methy_List[H].plusCover[pos+g-1]++;
							}
							//pthread_mutex_unlock(&meth_counter_mutex13);
						}
						else if (hitType==2 || hitType==4) {
							//pthread_mutex_lock(&meth_counter_mutex24);
							if (readString[k]=='G' && genome_Char=='G')
							{
								read_Methyl_Info[r] = 'M';
								if(Methratio ) ((ARGS *)arg)->Methy_List.NegMethylated[pos+g-1]++;
							//	if(hitType==2) 
							//	{
									if(genome_CharBac1=='C') 
									{
										read_Methyl_Info[r] = 'Z';met_CG++;
										//if(Methratio  && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2!='C') 
									{
										read_Methyl_Info[r] = 'H';met_CHH++;
										//if(Methratio  && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2=='C') 
									{
										read_Methyl_Info[r] = 'X';met_CHG++;
										//if(Methratio && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}
							}
							else if (readString[k]=='A' && genome_Char=='G')
							{
								read_Methyl_Info[r] = 'U';
								rawReadBeforeBS[k] = 'G';
								if(Methratio) ((ARGS *)arg)->Methy_List.NegUnMethylated[pos+g-1]++;
							//	if(hitType==2) 
							//	{
									if(genome_CharBac1=='C') 
									{
										read_Methyl_Info[r] = 'z';non_met_CG++;
										//if(Methratio && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2!='C') 
									{
										read_Methyl_Info[r] = 'h';non_met_CHH++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2=='C') 
									{
										read_Methyl_Info[r] = 'x';non_met_CHG++;
										//if(Methratio  && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}
							}
							else if (readString[k] != genome_Char) {
								
								read_Methyl_Info[r] = readString[k]; //genome_Char;
							}
							if(Methratio)
							{
								if(readString[k]=='C') ((ARGS *)arg)->Methy_List.NegG[pos+g-1]++;
								else if(readString[k]=='T') ((ARGS *)arg)->Methy_List.NegA[pos+g-1]++;
								//Methy_List[H].NegCover[pos+g-1]++;
							}
							//pthread_mutex_unlock(&meth_counter_mutex24);
						}
					}
					cigr++;n=0;lens+=length;Glens+=length;RLens+=length;
				}else if(*cigr=='I')
				{
					int i;temp[n]='\0';int length=atoi(temp);
					for(i=0;i<length;i++)
					{
						read_Methyl_Info[i+RLens] = 'I';
					}
					lens+=length;
					RLens+=length;
					cigr++;n=0;
				}else if(*cigr=='D')
				{
					int i;temp[n]='\0';unsigned int length=atoi(temp);
					for(i=0;i<length;i++)
					{
						read_Methyl_Info[i+RLens] = 'D';
					}
					Glens+=length;RLens+=length;
					cigr++;n=0;
				}else
				{
					//printf(" --%d%c-- ",atoi(temp),*cigr);
					CONTINUE=true;
					break;
					//continue;
					//exit(0);
				}
			}
			if(CONTINUE) continue;
							
			read_Methyl_Info[RLens]='\0';
			rawReadBeforeBS[lens]='\0';
			string mappingStrand="N";
			if(hitType==1) mappingStrand="YC:Z:CT\tYD:Z:f";
			else if(hitType==4) mappingStrand="YC:Z:CT\tYD:Z:r";
			else if(hitType==3) mappingStrand="YC:Z:GA\tYD:Z:r";
			else if(hitType==2) mappingStrand="YC:Z:GA\tYD:Z:f";
            if(Sam && ((ARGS *)arg)->OUTFILE != NULL) {
			    flockfile(((ARGS *)arg)->OUTFILE);
			    if(Sam ) //&& Nmismatch <= UPPER_MAX_MISMATCH )
			    {
				if(SamSeqBeforeBS) 
				{
					fprintf(((ARGS *)arg)->OUTFILE,"%s\t%u\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\t%s\tRA:Z:%s\n",Dummy,Flag,Chrom,pos,mapQuality,CIGr.c_str(),Chrom_P,pos_P,Insert_Size,rawReadBeforeBS,forQuality,Nmismatch,read_Methyl_Info,mappingStrand.c_str(), readString.c_str());
				}else 
				{
					fprintf(((ARGS *)arg)->OUTFILE,"%s\t%u\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\t%s\n",Dummy,Flag,Chrom,pos,mapQuality,CIGr.c_str(),Chrom_P,pos_P,Insert_Size,readString.c_str(),forQuality,Nmismatch,read_Methyl_Info,mappingStrand.c_str());
				}
			    }
			    funlockfile(((ARGS *)arg)->OUTFILE);
            }
		}
		if(REMOVE_DUP){
			((ARGS *)arg)->Marked_Genome[pos+G_Skip] |= Flag;
			((ARGS *)arg)->Marked_GenomeE[pos+G_Skip+readString.size()] |= Flag;
		}
	}
	free(s2t);
}

void Print_Mismatch_Quality(FILE* OUTFILE_MM, int L) {
	char bases[] = "ACGT";

	for (int i=0; i<L; i++) {
		for (int j=0; j<4; j++) {
			for (int k=0; k<4; k++) {
				fprintf(OUTFILE_MM,"%u\t", Mismatch_Qual[i][bases[j]][bases[k]]);
			}
		} fprintf(OUTFILE_MM,"\n");
	}

}

inline unsigned char Hash(char* S)
{
	assert(false);
	unsigned char C=0;
	for (int i=2;S[i];C+=i*S[i++]);
	return C;
}

inline float calc_Entropy (string readString, int L)
{ 
	short entropy_arr[255]={0};
	//int length=strlen(readString);
	for(int i=0; i<L; i++) entropy_arr[readString[i]]++;
	
	float entropy=0.0;
	for(int i=0; i<4; i++) {
		double p = 1.0*entropy_arr["ACGT"[i]]/L;
		if(p>0) entropy-=p*log(p);
	}
  //if( entropy_arr["ACGT"[1]] + entropy_arr["ACGT"[3]] ==L ||  entropy_arr["ACGT"[0]] + entropy_arr["ACGT"[2]] == L ) entropy = 0;
	return entropy;
}

int Get_String_Length(FILE* INFILE)
{
	char Buf[BATBUF],Dummy[BATBUF],Tag[BATBUF];int L;
	fgets(Buf,BATBUF,INFILE);
	fgets(Buf,BATBUF,INFILE);
	//space-trouble in read description
	for(int i=0;Buf[i];i++) if(Buf[i]==' ') Buf[i]='_';
	sscanf(Buf,"%s%s",Dummy,Tag);	
	for (L=0;Tag[L];L++);
	rewind (INFILE);
	return L;
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
/*
string Get_RealRead(char* cig,string & read)
{
	std::string ReadR="";
	char temp[5];unsigned lens=0;
	unsigned n=0;
	while(*cig!='\0')
	{
		if(*cig>='0' && *cig<='9')
		{
			temp[n]=*cig;
			cig++;n++;
		}else if(*cig=='S')
		{
			int i;temp[n]='\0';
			for(i=0;i<atoi(temp);i++)
			{
				ReadR.append("S");
			}
			lens+=atoi(temp);
			cig++;n=0;
		}else if(*cig=='M')
		{
			temp[n]='\0';std::cout<<read<<endl;
			ReadR.append(read.substr(lens,atoi(temp)));
			cig++;n=0;lens+=atoi(temp);
		}else if(*cig=='I')
		{
			int i;temp[n]='\0';
			for(i=0;i<atoi(temp);i++)
			{
				ReadR.append("I");
			}
			lens+=atoi(temp);
			cig++;n=0;
		}else if(*cig=='D')
		{
			int i;temp[n]='\0';
			for(i=0;i<atoi(temp);i++)
			{
				ReadR.append("D");
			}
			//lens+=atoi(temp);
			cig++;n=0;
		}else
		{
			printf(" --%c-- ",*cig);
			exit(0);
		}
	}
	return ReadR;
}
*/
float Pr(float Q)
{
        assert((int)Q>=0 && (int)Q<POWLIMIT-1);
        return(1-POW10[(int)Q]);
        //printf("table: %f\tlib: %f\n",POW10[(int)Q],1-pow(10,-Q/10));
        //return (1-pow(10,-Q/10));
}
void Build_Pow10()
{
	for(int Q=0;Q<POWLIMIT;Q++)
	{
		POW10[Q]=(pow(10,-float(Q)/10));
	}
}
void fetch(const char *str, char c1, char c2, char *buf){
    const char *pl = strchr(str, c1)+1,
        *ph = strchr(str, c2);
    strncpy(buf, pl, ph-pl-1);
    buf[ph-pl] = '\0';
}
void initMem(char* temp,int size,char* temp2)
{
	temp= ( char* )malloc(size);
	strcpy(temp,temp2);
}
/*
X   for methylated C in CHG context
x   for not methylated C CHG
H   for methylated C in CHH context
h   for not methylated C in CHH context
Z   for methylated C in CpG context
z   for not methylated C in CpG context
M   for methylated C in Unknown context (CN or CHN )
U   for not methylated C in Unknown context (CN or CHN)
=    for match bases
A/T/C/G   for mismatch bases
*/
string&   replace_all(string&   str,const   string&   old_value,const   string&   new_value)   
{   
    while(true)   {   
        string::size_type   pos(0);   
        if(   (pos=str.find(old_value))!=string::npos   )   
            str.replace(pos,old_value.length(),new_value);   
        else   break;   
    }   
    return   str;   
}  

std::string m_replace(std::string str,std::string pattern,std::string dstPattern,int count)
{
    std::string retStr="";
    string::size_type pos;

    int szStr=str.length();
    int szPattern=pattern.size();
    int i=0;
    int l_count=0;
    if(-1 == count) // replace all
        count = szStr;

    for(i=0; i<szStr; i++)
    {
        pos=str.find(pattern,i);

        if(std::string::npos == pos)
            break;
        if(pos < szStr)
        {
            std::string s=str.substr(i,pos-i);
            retStr += s;
            retStr += dstPattern;
            i=pos+pattern.length()-1;
            if(++l_count >= count)
            {
                i++;
                break;
            }
        }
    }
    retStr += str.substr(i);
    return retStr;
}
int strSearch(char *str1,char *str2)
{
       int at,flag=1;
       if (strlen(str2) > strlen(str1))
       {
           at = -1;
       }
       else if (!strcmp(str1,str2))
       {
           at = 0;
       }
       else
       {
            unsigned i=0,j=0;
            for (i=0;i < strlen(str1)&&flag;)
           {
                  for (j=0;j < strlen(str2)&&flag;)
                 {
                       if (str1[i]!=str2[j])
                       {
                              i++;
                              j=0;
                       }
                       else if (str1[i]==str2[j])
                       {
                              i++;
                              j++;
                       }
                       if (j==strlen(str2))
                       {
                             flag = 0;
                       }
                       if(i==strlen(str1)) break;
                 }
            }
            at = i-j;
            if(flag==1) at=-1;
       }
       return at;
}

char* replace(char*src, char*sub, char*dst)
{
    int pos =0;
    int offset =0;
    int srcLen, subLen, dstLen;
    char*pRet = NULL;

    srcLen = strlen(src);
    subLen = strlen(sub);
    dstLen = strlen(dst);
    pRet = (char*)malloc(srcLen + dstLen - subLen +1);//(memory)if (NULL != pRet)
    {
        pos = strstr(src, sub) - src;
        memcpy(pRet, src, pos);
        offset += pos;
        memcpy(pRet + offset, dst, dstLen);
        offset += dstLen;
        memcpy(pRet + offset, src + pos + subLen, srcLen - pos - subLen);
        offset += srcLen - pos - subLen;
        *(pRet + offset) ='\0';
    }
    return pRet;
}
int Get_ED(std::string & S)//edit distance
{
	int ed=0;
	int i;
	for(i=0;i<S.size();i++)//momomo
	{
		if(S[i]=='I' || S[i]=='D') ed++;
	}
	return ed;
}
// 
void MinandSec(unsigned a[MAX_HITS_ALLOWED][10], int left, int right, int & min, int & second,int & loci)
{
	if(left == right)
	{
		min = a[left][4] ;
		second =  INT_MIN;
		loci=left;
	}
	else if(left +1== right)
	{
		min = a[left][4] < a[right][4] ? a[left][4] : a[right][4] ;
		second = a[left][4] > a[right][4] ? a[left][4] : a[right][4] ;
		loci=a[left][4] < a[right][4] ? left : right ;
	}
	else
	{
		int mid =(right + left) /2 ;

		int leftmin ;
		int leftsecond ;
		int leftloci;
		MinandSec(a, left, mid, leftmin, leftsecond,leftloci) ;

		int rightmin ;
		int rightsecond ;
		int rightloci;
		MinandSec(a, mid +1, right, rightmin, rightsecond,rightloci) ;

		if (leftmin < rightmin)
		{
			min = leftmin ;
			second = leftsecond < rightmin ? leftsecond : rightmin ;
			loci=leftloci;
		}
		else
		{
			min = rightmin ;
			second = leftmin > rightsecond ? rightsecond : leftmin ;
			loci=rightloci;
		}
	}
	return;
}
/*
void MaxandSec(int a[], int left, int right, int&max, int&second)
{
	if(left == right)
	{
		max = a[left] ;
		second =  INT_MIN;
	}
	else if(left +1== right)
	{
		max = a[left] > a[right] ? a[left] : a[right] ;
		second = a[left] < a[right] ? a[left] : a[right] ;
	}
	else
	{
		int mid =(right + left) /2 ;

		int leftmax ;
		int leftsecond ;
		MaxandSec(a, left, mid, leftmax, leftsecond) ;

		int rightmax ;
		int rightsecond ;
		MaxandSec(a, mid +1, right, rightmax, rightsecond) ;

		if (leftmax > rightmax)
		{
			max = leftmax ;
			second = leftsecond > rightmax ? leftsecond : rightmax ;
		}
		else
		{
			max = rightmax ;
			second = leftmax < rightsecond ? rightsecond : leftmax ;
		}
	}
}
*/

string remove_soft_split(string & cigar,int & Read_Len,int & pos)
{

	const char* cig=cigar.c_str();
	char temp[20];
	int n=0,lens=0,length=0,RLens=0,Glens=0;
	int genome_move_size=0; int headClip=0;int moveSize=0;
	char cigar_rm[1000];*cigar_rm=0;
	char buffer_cigar[100];*buffer_cigar=0; bool cigar_remove=false;
	char upper='N';char upper_cigar[1024]="\0";int upper_length=0;
	while(*cig!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
	{
		if(*cig>='0' && *cig<='9')
		{
			temp[n]=*cig;
			cig++;n++;
		}else if(*cig=='S')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length==0) cigar_remove=true;
			if(length>0) 
			{
				if(upper=='M')
					sprintf(upper_cigar,"%dM",length+upper_length);
				else
					sprintf(upper_cigar,"%dM",length);
				
				if(upper=='N')
					pos-=length;
			}

			lens+=length;
            RLens+=length;
            cig++;n=0;
            upper='M';upper_length=length;
		}else if(*cig=='M')
		{
			temp[n]='\0';int length=atoi(temp);
			if(length==0) cigar_remove=true;
			if(length>0) 
			{
				if(upper=='M')
					sprintf(upper_cigar,"%dM",length+upper_length);
				else
					sprintf(upper_cigar,"%dM",length);
			}

			cig++;n=0;lens+=length;Glens+=length;RLens+=length;
			char buf[1024]="\0";
			sprintf( buf , "%dM",Read_Len);
			char buf_tmp[1024]="\0";
			sprintf( buf_tmp, "%dM",length);
			if( !strcmp(buf,  buf_tmp ) ) 
			{ 
				if(!strcmp(buf, cig))
					break;
				else
				{
					cigar_remove=true;
					strcpy(cigar_rm, buffer_cigar);
				}
				
			}
			upper='M';upper_length=length;
		}else if(*cig=='I')
		{
			if(upper=='M')
				strcat(cigar_rm, upper_cigar);
			int i;temp[n]='\0';int length=atoi(temp);

			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1I");
				else if(length==2) sprintf(buffer_cigar,"2I");
				else sprintf(buffer_cigar,"%dI",length);
				strcat(cigar_rm, buffer_cigar);
			}

			lens+=length;RLens+=length;
			cig++;n=0;genome_move_size-=length;moveSize++;
			upper='I';
		}else if(*cig=='D')
		{
			if(upper=='M')
				strcat(cigar_rm, upper_cigar);
			int i;temp[n]='\0';int length=atoi(temp);

			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1D");
				else if(length==2) sprintf(buffer_cigar,"2D");
				else sprintf(buffer_cigar,"%dD",length);
				strcat(cigar_rm, buffer_cigar);
			}

			Glens+=length;RLens+=length;
			cig++;n=0;genome_move_size+=length;moveSize++;
			upper='D';
		}else
		{
			printf(" --%d%c-- ",atoi(temp),*cig);
			break;
		}
	}
	if(upper=='M')
		strcat(cigar_rm, upper_cigar);
	string s(cigar_rm);
	return s;
}

void ReverseC_Context(char* Dest,const char* seq,int & stringlength)
{
	if(stringlength!=5 || strlen(seq)!=5) {
		fprintf(stderr, "\nError string %d %d %s\n", stringlength, strlen(seq), seq);
		exit(0);
	}
	
        for (int i=stringlength-1;i>=0;i--)
        {
                *Dest=Char2Comp[seq[i]];
                Dest++;
        }
        *Dest=0;
        //return Dest;
}

int fprintf_time(FILE *stream, const char *format, ...)
{
        time_t timer;
        char buffer[26];
        struct tm* tm_info;
        va_list arg;
        int done;

        time(&timer);
        tm_info = localtime(&timer);
        strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", tm_info);

        fprintf(stream, "[%s] ", buffer);

        va_start(arg, format);
        done = vfprintf(stream, format, arg);
        va_end(arg);

        return done;

}

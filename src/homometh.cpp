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

struct Gene_Hash
{
        char* Genome;
        int Index;
};
using namespace std;
bool Collision=false;
int samplecol = 10;
map <string,int> String_Hash;
map <string,int> Context_Hash;


//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_Progress(float Percentage);
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
typedef struct {
    char* id;
    char* chrom;
    int position;
    short altalleles;
    char* RA;
    char* AA; // alternate alleles
    double* GLL; // genotype likelihoods added 11/25/13
    char* genotype; // encoded as integers 0 1 2 3 4 5 6 7
    short type;
    char* allele1;
    char* allele2; // temporary for SNPs
    char heterozygous; // only heterozygous variants will be used for printing out HAIRS
    int depth;
    int A1, A2;
    int H1, H2;
} VARIANT;
int parse_variant(VARIANT* variant, char* buffer);

//}-----------------------------   GLOBAL VARIABLES  -------------------------------------------------
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..
long AverPerCG=0,AverPerCHG=0,AverPerCHH=0;
long AverCG=0,AverCHG=0,AverCHH=0;
int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	printf("\nBatMeth2:  MethyHomo v1.0\n");
	const char* Help_String="Command Format :   methyHomo [options] -G genome -o/--out <OUT> --vcf <VCF file> --strand +/- \n"
		"\nUsage:\n"
		"\t-G|--genome      Genome\n"
		"\t--vcf            VCF file.\n"
		"\t--out            homo out file\n"
		"\t--context        default all context\n"
		"\t--strand         strand\n"
		"\t-h|--help";

	int Genome_CountX=0;
	//char* methInfileName;
	int InFileStart=0,InFileEnd=0;
	string Geno;
	char *InFile;
	char* Outfile;
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
	int Coverage = 5;
	int nC=3;
	char STRAND;
	//
	if(argc<6)
	{
		printf("error args\n%s \n",Help_String);
		exit(0);
	}
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-o") ||!strcmp(argv[i], "--out")  )
		{
			Outfile=argv[++i];
                }else if( !strcmp(argv[i], "--strand") ){
                        STRAND=argv[++i][0];
                }
		else if(!strcmp(argv[i], "-G") || !strcmp(argv[i], "--genome"))
		{
			Geno=argv[++i];
		}else if(!strcmp(argv[i], "-C") || !strcmp(argv[i], "--coverage"))
		{
			Coverage=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-nC"))
		{
			nC=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-s") || !strcmp(argv[i], "--step"))
		{
			binspan=atof(argv[++i]);
			nLevel=ceil(1/(double)binspan)-1;
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
		else if(!strcmp(argv[i], "-S") || !strcmp(argv[i], "--chromStep"))
		{
			chromStep=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-d") || !strcmp(argv[i], "--distance"))
		{
			distance=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i],"--vcf") )
		{
			InFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
			//methInfileName=argv[++i];
		}
		else if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--GFF"))
		{
			InputGff=true;
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
			time(&Start_Time);

                        string G=Geno;G+=".bin";
                        string L=Geno;L+=".ann.location";

                        FILE* BINFILE=File_Open(G.c_str(),"r");
                        FILE* Location_File=File_Open(L.c_str(),"r");
                        fseek(BINFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(BINFILE);rewind(BINFILE);//load original genome..
			char* Org_Genome;
                        Org_Genome=new char[Genome_Size];if(!Org_Genome) throw("Insufficient memory to load genome..\n");
                        if(!fread(Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");

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

			char* Split_Point=Org_Genome;//split and write...

			Gene_Hash* Genome_List;
			Genome_List = new Gene_Hash[Genome_Count];
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				printf("%s, ",Genome_Offsets[i].Genome);
				String_Hash[Genome_Offsets[i].Genome]=i;
				Genome_List[i].Genome=(Split_Point+=Genome_Offsets[i].Offset);
                                Genome_List[i].Index=i;
			}
			printf("Loaded\n");

                        fclose(BINFILE);
                        fclose(Location_File);	

		
			char Buf[BATBUF],Meth[BATBUF],Chrom[CHROMSIZE],context[BATBUF],effCT[BATBUF];
			unsigned pos;float methratio=0;int countC=0,countCT=0;
			printf("\nLoading methyratio file...\n");
			FILE* INFILE;//=File_Open(methInfileName,"r");
			FILE* fpout = File_Open(Outfile, "w");
		for(int f=InFileStart;f<=InFileEnd;f++)
		{
			INFILE=File_Open(argv[f],"r");
			int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
			fseek(INFILE, 0L, SEEK_END);off64_t File_Size=ftello64(INFILE);rewind(INFILE);
			Init_Progress();
			fgets(Buf,BATBUF,INFILE);//read first header marker..
			int start = 0, end = 0;
			bool printCut = false, printH = false;
			int prevchrom = -1;
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
				if(Meth[0] == '#') continue;
				VARIANT *variant;
				variant = (VARIANT*) malloc(sizeof(VARIANT));
				int type = parse_variant(variant, Meth);
	//			sscanf(Meth,"%s\t%u\t%c\t%s\t%d\t%d",Chrom,&pos,&Strand,context,&countC,&countCT);
				pos = variant->position;
				end = pos-1;
				int H=String_Hash[Chrom];
				if(H != prevchrom || prevchrom == -1) {
					start=0;
					fprintf(fpout, "#=====\n");
				}
				prevchrom = H;

				printH = false;

				if( type != 0) {
					//fprintf(fpout, "#====\n");
					printH = true;
					printCut=true;
					start = pos;
					continue;
				}
				
				if(printCut) {
					fprintf(fpout, "#====\n");
					printCut = false;
				}
				fprintf(fpout, "1\t0\t1\t%s\t%d\t%s\t%s\t%s\n", variant->chrom, pos, variant->RA, variant->AA, variant->genotype);
				start = pos;
				
			}//end read file while
			Done_Progress();
		} //end foreach meth file
			printf("Methratio file Loaded.\n");
			fclose(INFILE);
			fclose(fpout);
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

int parse_variant(VARIANT* variant, char* buffer) {
    int i = 0, j = 0, k = 0, s = 0, e = 0;
    char* tempstring;
    int col = 10;
    int flag = 0;

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    variant->chrom = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) variant->chrom[j - s] = buffer[j];
    variant->chrom[j - s] = '\0';
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    variant->position = atoi(tempstring);
    free(tempstring);

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i; // varid
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    variant->RA = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) variant->RA[j - s] = buffer[j];
    variant->RA[j - s] = '\0';
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    variant->AA = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) variant->AA[j - s] = buffer[j];
    variant->AA[j - s] = '\0';

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;

    while (buffer[i] != '\n') {
        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t' && buffer[i] != '\n') i++;
        e = i;
        if (col == samplecol) {
            variant->genotype = (char*) malloc(e - s + 1);
            for (j = s; j < e; j++) variant->genotype[j - s] = buffer[j];
            variant->genotype[j - s] = '\0';
        } else col++;
    }

    if (variant->genotype[0] != '2' && variant->genotype[2] != '2') // both alleles are 0/1
    {
        variant->allele1 = (char*) malloc(strlen(variant->RA) + 1);
        strcpy(variant->allele1, variant->RA);
        j = 0;
        while (variant->AA[j] != ',' && j < strlen(variant->AA)) j++;
        variant->allele2 = (char*) malloc(j + 1);
        for (i = 0; i < j; i++) variant->allele2[i] = variant->AA[i];
        variant->allele2[i] = '\0';
        variant->type = strlen(variant->allele2) - strlen(variant->allele1);
    } else if (variant->genotype[0] == '0' || variant->genotype[2] == '0') // at least one allele is reference
    {
        variant->allele1 = (char*) malloc(strlen(variant->RA) + 1);
        strcpy(variant->allele1, variant->RA);
        j = 0;
        while (variant->AA[j] != ',' && j < strlen(variant->AA)) j++;
        k = j + 1;
        while (variant->AA[k] != ',' && k < strlen(variant->AA)) k++;
        variant->allele2 = (char*) malloc(k - j + 1);
        for (i = j + 1; i < k; i++) variant->allele2[i - j - 1] = variant->AA[i];
        variant->allele2[i - j - 1] = '\0';
        variant->type = strlen(variant->allele2) - strlen(variant->allele1);
    } else // reference allele is missing 1/2 case
    {
        j = 0;
        while (variant->AA[j] != ',' && j < strlen(variant->AA)) j++;
        variant->allele1 = (char*) malloc(j + 1);
        for (i = 0; i < j; i++) variant->allele1[i] = variant->AA[i];
        variant->allele1[i] = '\0';
        k = j + 1;
        while (variant->AA[k] != ',' && k < strlen(variant->AA)) k++;
        variant->allele2 = (char*) malloc(k - j + 1);
        for (i = j + 1; i < k; i++) variant->allele2[i - j - 1] = variant->AA[i];
        variant->allele2[i - j - 1] = '\0';
        variant->type = strlen(variant->allele2) - strlen(variant->allele1);
    }

    i = strlen(variant->allele1) - 1;
    j = strlen(variant->allele2) - 1;
    flag = 0;
    while (i > 0 && j > 0) {
        if (variant->allele1[i] != variant->allele2[j]) break;
        i--;
        j--;
        flag++;
    }
    variant->allele1[i + 1] = '\0';
    variant->allele2[j + 1] = '\0';

    if (variant->type != 0) variant->position++; // add one to position for indels

    if ((variant->genotype[0] == '0' && variant->genotype[2] == '1') || (variant->genotype[0] == '1' && variant->genotype[2] == '0')) {
        variant->heterozygous = '1'; // variant will be used for outputting hairs
        return 1;
    }
    if ((variant->genotype[0] == '0' && variant->genotype[2] == '2') || (variant->genotype[0] == '2' && variant->genotype[2] == '0')) {
        variant->heterozygous = '1'; // variant will be used for outputting hairs
        return 1;
    } else if (variant->genotype[0] != variant->genotype[2] && variant->genotype[0] != '.' && variant->genotype[2] != '.') {
        variant->heterozygous = '2';
        return 0;
    } else {
        variant->heterozygous = '0';
        return 0;
    }
    //free(variant->genotype); free(variant->AA); free(variant->RA); free(variant->chrom);
}

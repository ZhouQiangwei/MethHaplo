#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <map>
#include <deque>
#include <vector>
#include <time.h>
#include <algorithm>
#include "methyhaplo.hpp"
/**
 * Input format, sam format from batmeth2  MD:Z:=================h==========hh=hhh=====h=================h=====h================x==h=
*SRR179602.48    0       chr3    123850297       60      85M     *       0       0       AAGGAGTTTTTGAAAAATATTAAAAGGGTTATTTAAATATTAAGAAGTAAGATTGTATTTGAATTTGATGAAATTGATGGTAGTT   DEC@=CEBECEEAEE5DD<>CACAC@ADDBBDD?D5=CCCD?DD=ACAC@C,>;?DCDDA?C5CC?D-:DD-:=5DBAD?CBBBA   XT:A:U  NM:i:0  MD:Z:=================h==========hh=hhh=====h=================h=====h================x==h=      XB:Z:BSW
*SRR179602.62    16      chr2    212475086       57      85M     *       0       0       CTACAAAAGCAATATTTTATATCCTTAACTCTAGAAACAAATCATTAACAAAACCATTAAAATTTAACATATTAAAAATTACTAA   C4/?7?+:+B@?<D@/A?>>:5?D?C?@C-B,B<8<>-:>B-DDB3:?A?AC-AAD:DA:C:DD=EEAE:DEE:EDDDBD:CCCC   XT:A:U  NM:i:1  MD:Z:==x==xh=A==x======h=============xH=====x==================h=hh===h=======hh=====h==x=      XB:Z:BSC
 *
 *###############  flag_num ##############
 *0x0001	the read is paired in sequencing, no matter whether it is mapped in a pair
 *0x0002	the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)
 *0x0004	the query sequence itself is unmapped //if(flag_num & 0x4)  //no match   if(!(flag_num&0x4)) //match
 *0x0008	the mate is unmapped 
 *0x0010    strand of the query (0 for forward; 1 for reverse strand)  // if(flag_num & 0x10) ---reverse
 *0x0020    strand of the mate 1 // if(flag_num & 0x20) --reverse
 *0x0040	the read is the first read in a pair 
 *0x0080	the read is the second read in a pair 
 *
 *possible usage:
 *  if ( !(flag_num & 0x4) and !(flag_num & 0x8) ) //left and right all mapped
 *		if (flag_num & 0x2)   //proper_pairs
 *		else if(unique_discordant_pairs) //unique && !(flag_num&0x2)
 *			if(items[6] != '=')  //diffchr_pairs
 *			elsif((flag_num & 0x10) and (flag_num & 0x20))  //discordantRR_pairs  //first '-' second '-' 
 *			elsif(!(flag_num & 0x10) and !(flag_num & 0x20))  //discordantFF_pairs //first '+' second '+' 
 *			elsif( ((flag_num & 0x10) and (flag_num & 0x40)) or (!(flag_num & 0x10) and (flag_num & 0x80)) ) //discordantRF_pairs //first '-'  || second '+' 
 *			else //discordantFR_minus_pairs
 *  elsif (flag_num & 0x8) and !(flag_num & 0x4)  //left mapped and right unmapped  --------singletons
 *  elsif !(flag_num & 0x8) and (flag_num & 0x4)  //left unmapped and right mapped  --------singletons
 *  else //unmapped
 * 
 * if ( !(flag_num & 0x4) and (flag_num & 0x40) )// match and first left
 * if ( !(flag_num & 0x4) and (flag_num & 0x80) )// match and second right
 * 
 *In a string FLAG :
 *   p=0x1 (paired), P=0x2 (properly paired), u=0x4 (unmapped),
 *   U=0x8 (mate unmapped), r=0x10 (reverse), R=0x20 (mate reverse)
 *   1=0x40 (first), 2=0x80 (second), s=0x100 (not primary),
 *   f=0x200 (failure) and d=0x400 (duplicate). 
 *   
 *	if flag_num & 0x10:
 *		first = '-'
 *	else:
 *		first = '+'
 *	if flag_num & 0x20:
 *		second = '-'
 *	else:
 *		second = '+'
 *
 *###############################
 * Assume that 
 * 1) the sam file shouls sorted
 * 2) all the mapped reads are sorted by the coordinates in ascending order
 *
 *
 *heter snp format:
 *chrM    16259   rs376682258     C       A,T  || chrM    16259   rs376682258     C       AT || chrM    16259   rs376682258     C       A
 *
 * output
 *  methylation pair statistics
 *  chrom <tab> coordinate1 <tab> coordinate2 <tab> #MM <tab> #MU <tab> #UM <tab> #UU
 *  M: for methylated coverage
 *  U: for un-methylated coverage
 * if heter SNP file also definded, output
 * methylation and heter SNP pair statistics
 * chrom <tab> coordinate1 <tab> coordinate2 <tab> #MM <tab> #MU <tab> #UM <tab> #UU
 * chrom <tab> coordinate1 <tab> coordinate2 <tab> #MR <tab> #MV <tab> #UR <tab> #UV <tab> Ref <tab> var <tab> 1/2 (coordinate1/coordinate2 is SNP)
 * chrom <tab> coordinate1 <tab> coordinate2 <tab> #RM <tab> #VM <tab> #RU <tab> #VU <tab> Ref <tab> var <tab> 1/2
 * chrom <tab> coordinate1 <tab> coordinate2 <tab> #MV1 <tab> #MV2 <tab> #UV1 <tab> #UV2 <tab> Ref <tab> var12 <tab> 1/2
 * chrom <tab> coordinate1 <tab> coordinate2 <tab> #V1M <tab> #V2M <tab> #V1U <tab> #V2U <tab> Ref <tab> var12 <tab> 1/2
 *  M: for methylated coverage
 *  U: for un-methylated coverage
 *  R: Reference base
 *  V: Variant base
 *
 * pairing: the first valid C with methylation, and the next valid C with methyaltion
 * pairing: the first valid C with methylation or heter SNP, and the next valid C with methyaltion or heter SNP
 * procedure:
 * 1) read the first mapped result
 *    extract the first start coordinate (first_start_coordinate)
 *    extract the first end   coordinate (first_end_coordinate), from sequence length
 *    extract the potential Cs
 * 2) read another mapped result
 *    extract the latest start coordinate
 *    extract the latest end   coordinate
 *    extract the potential Cs
 * 3) compare the first_end_coordiante and the latest start coordinate
 * 4)   if first_end_coordinate > latest start coordinate, go back to 2): read another mapped result
 * 5)   otherwise,
 * 6)     process the first mapped result
 * 7)     remove the first mapped result, and assign a new first mapped result; go back to 2)
 * 8) if reach end of the input file
 * 9)   process the first mapped result
 * 10)  remove the first mapped result, and assign a new first mapped result; go back to 9)
 * 11) if all mapped results are processed, finish
 * read mapped results until the new coordinates larger than the top
 *
 * @author glli
 * @author qwzhou
 */
using namespace std;

bool hetero_SNP=false;
map <string,int> String_Hash;
map <string,int> pairStatMap;
map <int,string> intpairStatMap;
char Char2Trans[2000];
unsigned chrom_offset=1;
int Quality_cutoff=10;
unsigned int GenomeSize=400;
GenomeSNP chromSNPs[400];
int minIS=0;
int maxIS=1000;
//----------------  main function -------------
void extractPotentialCs_single(alignedread read,deque<int> & potentialCs,int & last_coordinate);
int process_meth(char* buffer,heterSNP* var_t, char* chrom, int NMETH, int NCOVER, float MFloat);
bool sort_methvar(const heterSNP & h1, const heterSNP & h2);
bool uniq_and_mod(std::deque< heterSNP > & v);

int main(int argc, char* argv[])
{
	printf("\nBatMeth2::methyhaplo: haplotype assembly for DNA methylation(BSseq) next-generation sequencing data. \nversion: v1.0\n\n");
	const char* Help_String="Command Format :  methyhaplo [options] -o <haplotype result> -m <methratio file> -s <sam format file>/-b <bam>\n"
		"\n======================== METHYHAPLO OPTIONS ========================\n\n"
		"\t-o|--out              Output file prefix\n"
		"\t-s|--sam              sam file from batmeth2-align.  This file should be coordinate sorted, using the <samtools sort> command, and must contain methylstatus[MD:Z:].\n"
		"\t-b|--bam              bam file, should be coordinate sorted. (use this option or -s option but not both)\n"
		"\t-q INT                only process reads with mapping quality >= INT [default >= 20].\n"
		"\t-v|--snp              SNP file with genotypes for a single individual in BS-SNPer format.\n" 
		"\t--vcf                 SNP file with genotypes for a single individual in VCF format, V4.1. (use this option or -v option but not both)\n"
		"\t--PE                  paired-end reads.[default:single-end]\n"
		"\t-m|--methf            methratio file from batmeth2.\n"
		"\t-M|--NMETH            Number of methylated reads cover cytosine site. default: 2 [m>=2]\n"
		"\t-N|--NCOVER           Number of coverage reads in cytosine site. default: 6 [n >= 6]\n"
		"\t-f|--MFloat           cutoff of methratio. default: 0.2 [ f =< meth <= 1-f]\n"
		"\t--minIS <INT>         minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0\n"
		"\t--maxIS <INT>         maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000\n"
		"\t-h|--help\n";
//-------------format------------------------
	Char2Trans['M']=Char2Trans['m']=Char2Trans['X']=Char2Trans['H']=Char2Trans['Z']='M';
	Char2Trans['U']=Char2Trans['u']=Char2Trans['x']=Char2Trans['h']=Char2Trans['z']='U';
	Char2Trans['G']=Char2Trans['g']='G';
	Char2Trans['T']=Char2Trans['t']='T';
	Char2Trans['A']=Char2Trans['a']='A';
	Char2Trans['C']=Char2Trans['c']='C';
	Char2Trans['N']=Char2Trans['n']=Char2Trans['=']='=';
	Char2Trans['I']=Char2Trans['D']='O';//open a gap
//-------------------------------------------------
	//char* SAM_fileName;
	char* Align_fileName;bool bamformat=false;
	bool Paired=false;
	char* Output_Name;
	char pos_Output_tmp[100];
	char neg_Output_tmp[100];
	char* SNP_fileName;bool VCF=false;
	char* Meth_fileName; bool hetero_meth = false;
	bool paired_end=false;
	//meth
    int NMETH = 2;
    int NCOVER = 6; 
    float MFloat = 0.2;
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-o") ||!strcmp(argv[i], "--out")  )
		{
			Output_Name=argv[++i];
		}
		else if(!strcmp(argv[i], "-s") || !strcmp(argv[i], "--samfile"))
		{
			Align_fileName=argv[++i];
		}else if(!strcmp(argv[i], "--PE"))
		{
			paired_end=true;
		}else if(!strcmp(argv[i], "-m") || !strcmp(argv[i], "--methf"))
		{
			Meth_fileName=argv[++i];
			hetero_meth = true;
		}else if(!strcmp(argv[i], "-M") || !strcmp(argv[i], "--NMETH"))
		{
			NMETH=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-N") || !strcmp(argv[i], "--NCOVER"))
		{
			NCOVER=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-f") || !strcmp(argv[i], "--MFloat"))
		{
			MFloat=atof(argv[++i]);
		}
		else if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bamfile"))
		{
			Align_fileName=argv[++i];
			bamformat=true;
		}
		else if(!strcmp(argv[i], "-q"))
		{
			Quality_cutoff=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--snp"))
		{
			SNP_fileName=argv[++i];
			hetero_SNP=true;
		}else if(!strcmp(argv[i], "--vcf") )
		{
			SNP_fileName=argv[++i];
			hetero_SNP=true;
			VCF=true;
		}
		else
		{
			printf("%s \n",Help_String);
			exit(0);
		}
	}
	//printf("\n%f %f\n",fishers_exact(9,0,1,3),gsl_cdf_hypergeometric_Q(8,9,4,10) );
	if(argc<=3 || !hetero_meth) printf("%s \n",Help_String);
	if (argc >1 && hetero_meth) 
	{
		try
		{
			long nReads=0;
			char s2t[BATBUF],Dummy[BATBUF]; //s2t read the file line, Dummy- tmp
			//process the hetero SNP file, and store the hetero SNP. //chrM    16259   rs376682258     C       AT
			FILE* SNP_FILE;string old_chr_snp;string now_chr_snp;
			int read_first_line=0; //indicator the first snp
			unsigned Genome_Count=0;
			
			if(hetero_SNP) 
			{
				//init_pairStatMap(pairStatMap,intpairStatMap);
				SNP_FILE=File_Open(SNP_fileName,"r");
				printf("Reading hetero SNP file \n\n");
				char Ref_seq[BATBUF];char var_seq[BATBUF];char chrom[CHROM_NAME_LEN];int Quality=0;char Filter[BATBUF];char Genotype[BATBUF];char Indel[BATBUF];

				while(fgets(s2t,BATBUF,SNP_FILE)!=0)
				{
					if(s2t[0]=='#') continue;
					struct heterSNP heterSNP_tmp;
					bool INDEL=false;
					if(!VCF) //---- SNPer
					{
						sscanf(s2t,"%s%d%s%s%s%d%s%s",chrom,&heterSNP_tmp.pos,Dummy,Ref_seq,var_seq,&Quality,Filter,Genotype);
						if( strcmp(Filter,"PASS") || (strlen(Genotype)==2 && Genotype[0]==Genotype[1]) ) continue; //only need hetero SNP
						if(strlen(Ref_seq)!=1 || strlen(var_seq) > 3) continue; // not a snp, may be indels. //this part need update, indels
					}else //-----------------------VCF 4.1 --------------------- 
					{
						int samplecol=10; // default if there is a single sample in the VCF file
						int i=0,j=0,k=0,s=0,e=0; int col=10; int flag =0;
						char* tempstring;

						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i;
						for (j=s;j<e;j++) chrom[j-s] = s2t[j]; chrom[j-s] = '\0';
						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i;
						tempstring = (char*)malloc(e-s+1); for (j=s;j<e;j++) tempstring[j-s] = s2t[j]; tempstring[j-s] = '\0';
						heterSNP_tmp.pos = atoi(tempstring); free(tempstring);

						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i; // varid
						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e = i;
						for (j=s;j<e;j++) Ref_seq[j-s] = s2t[j]; Ref_seq[j-s] = '\0';
						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i;
						for (j=s;j<e;j++) var_seq[j-s] = s2t[j];  var_seq[j-s] = '\0';
						
						//anotation this line, so now we not just deal with snp.
						//if( strlen(Ref_seq)!=1 || strlen(var_seq) > 3 ) continue;
						
						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i;
						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i;
						for (j=s;j<e;j++) Filter[j-s] = s2t[j];  Filter[j-s] = '\0';
						if( strcmp(Filter,"PASS") ) continue; //col 7
						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i;
						for (j=s;j<e && j<5;j++) Indel[j-s] = s2t[j];  Indel[j-s] = '\0'; //line 8
						if( !strcmp(Indel,"INDEL") ) INDEL=true;
						while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t') i++; e= i;
						
						while (s2t[i] != '\n')
						{
							while (s2t[i] == ' ' || s2t[i] == '\t') i++; s = i; while (s2t[i] != ' ' && s2t[i] != '\t' && s2t[i] != '\n') i++; e= i;
							if (col == samplecol)
							{
								for (j=s;j<e;j++) Genotype[j-s] = s2t[j]; Genotype[j-s] = '\0';
							}
							else col++;
						}

						if ( (Genotype[0] =='0' && Genotype[2] == '1') || (Genotype[0] =='1' && Genotype[2] == '0')); 
						else if ( (Genotype[0] =='0' && Genotype[2] == '2') || (Genotype[0] =='2' && Genotype[2] == '0')); 
						else if( (Genotype[0] =='1' && Genotype[2] == '2') );
						else continue;
					}//end readline of vcf/snp
					
					now_chr_snp=chrom;
					if(read_first_line==0)
					{
				    		old_chr_snp=now_chr_snp;
				    		read_first_line=-1;
						String_Hash[chrom]=Genome_Count; // the reflect of chrom and index
						chromSNPs[Genome_Count].index=Genome_Count;
						strcpy(chromSNPs[Genome_Count].chrom,chrom);
				   	}
					if( old_chr_snp!=now_chr_snp )
					{
						Genome_Count++;
						String_Hash[chrom]=Genome_Count; // the reflect of chrom and index
						chromSNPs[Genome_Count].index=Genome_Count;
						strcpy(chromSNPs[Genome_Count].chrom,chrom);
					}

					heterSNP_tmp.RefBase=Ref_seq[0];
					if(INDEL)
					{
						heterSNP_tmp.doubleVar=false;
						heterSNP_tmp.variantBase1='O';//Open a gap
					}
					else if( strlen(var_seq)==1 ) //snp
					{
						heterSNP_tmp.doubleVar=false;
						heterSNP_tmp.variantBase1=var_seq[0];
						/*
						if(Ref_seq[0]=='C' && var_seq[0]=='T'){
							nReads++;
							old_chr_snp=now_chr_snp;
							continue;
						} //because if this is  a SNP so that the meth can also get the validCs.
						*/ //but can't known the details
					}else if(strlen(var_seq) >1) 
					{
						heterSNP_tmp.doubleVar=true;
						if(var_seq[1]==',') // VCF V4.1
						{
							heterSNP_tmp.variantBase1=var_seq[0];
							heterSNP_tmp.variantBase2=var_seq[2];
						}else //program: SNPer result 
						{
							heterSNP_tmp.variantBase1=var_seq[0];
							heterSNP_tmp.variantBase2=var_seq[1];
						}
					}
					heterSNP_tmp.methv=0;
					heterSNP_tmp.strand='.';
					chromSNPs[Genome_Count].HeterSNPs.push_back(heterSNP_tmp);
					
			
					nReads++;
					old_chr_snp=now_chr_snp;
					if (nReads % 1000 == 0) {
						Show_Progress_SNP(nReads);
					}
				}//end of read snp file
				Show_Progress_SNP(nReads);
			}// end of snp
			//----------done the snp file
			
			//process meth file
			nReads = 0;
			FILE* METH_FILE;
			read_first_line=0; //indicator the first snp
			
			METH_FILE=File_Open(Meth_fileName,"r");
			printf("Reading methratio file \n");
			char chrom[CHROM_NAME_LEN];

			while(fgets(s2t,BATBUF,METH_FILE)!=0)
			{
				if(s2t[0]=='#') continue;
				struct heterSNP heterSNP_tmp;
				int methfilter = process_meth(s2t, &heterSNP_tmp, chrom, NMETH, NCOVER, MFloat);
				if(methfilter != 1) continue;
				if(read_first_line==0 || old_chr_snp != chrom) {
					read_first_line = -1;
					if(String_Hash.find(chrom) == String_Hash.end()) {
						Genome_Count++;
						String_Hash[chrom]=Genome_Count;
						chromSNPs[Genome_Count].index=Genome_Count;
						strcpy(chromSNPs[Genome_Count].chrom,chrom);
					}
					old_chr_snp = chrom;
				}
				int H = String_Hash[chrom];
				chromSNPs[H].HeterSNPs.push_back(heterSNP_tmp);
				nReads++;
				if (nReads % 1000 == 0) {
					Show_Progress_Meth(nReads);
				}
			}
			//end meth file
			//sort meth and variant by position
			map<string,int>::iterator itchr;
   			for(itchr=String_Hash.begin();itchr!=String_Hash.end();++itchr){
   				//string key=itchr->first; //chrom
   				int H = itchr->second;
   				sort(chromSNPs[H].HeterSNPs.begin(), chromSNPs[H].HeterSNPs.end(), sort_methvar);
   				if(itchr->first == "chr1") uniq_and_mod(chromSNPs[H].HeterSNPs);
   			}
   			/*
			//test meth function
			int T = String_Hash["chr1"];
			while(chromSNPs[T].HeterSNPs.size()>0){
				chromSNPs[T].HeterSNPs.pop_front();
			}
			return 0;
			*/
			if(paired_end)
			{
				printf("\nPaired-end ");
				paiedend_methyhaplo(bamformat,Align_fileName,Output_Name);
				return 1;
			}else
				printf("\nSingle-end ");
			
			//process SAM file
			FILE* SAM_File;
			samfile_t *bamin = 0;bam1_t *b;bam_header_t *header;
			if(bamformat)
			{
				if ((bamin = samopen(Align_fileName, "rb", 0)) == 0) {
					fprintf(stderr, "fail to open \"%s\" for reading.\n", Align_fileName);
				}
				b = bam_init1();
				header=bam_header_dup((const bam_header_t*)bamin->header);
				printf("\nReading bam file \n");
			}
			else
			{
				SAM_File=File_Open(Align_fileName,"r");
				printf("\nReading sam file \n");
			}
			deque<int> pos_validCs;
			deque<int> Neg_validCs;
			deque<alignedread> pos_vDnaMethylMap;
			deque<alignedread> Neg_vDnaMethylMap;
			//pvalues
			vector<Pvals> Pvalues;
			//read sam read mapping
		        int pos_first_start_coordinate = -1; //+ strand
		        int pos_first_end_coordinate   = -1;
		        int pos_processed_coordinate   = -1; // the last coordinate already processed
		        int pos_last_coordinate        = -1; // the last coordinate for all the results already read
		        int neg_first_start_coordinate = -1; // - strand
		        int neg_first_end_coordinate   = -1;
		        int neg_processed_coordinate   = -1; // the last coordinate already processed
		        int neg_last_coordinate        = -1; // the last coordinate for all the results already read
		        
		        int withFirstResult = 0; // indicator whether the first mapped resul is read
		        //processed the read number
		        nReads = 1;
		        int pos_print=0,neg_print=0;
		        //open output file and write the results
		        sprintf(pos_Output_tmp,"%s.tmp.plus.txt",Output_Name);
		        sprintf(neg_Output_tmp,"%s.tmp.neg.txt",Output_Name);
			FILE* pos_OUTFILE=File_Open(pos_Output_tmp,"w");
			FILE* neg_OUTFILE=File_Open(neg_Output_tmp,"w");
			fprintf(pos_OUTFILE,"#chrom\tcoordinate1\tcoordinate2\t{ #MM\t#MU\t#UM\t#UU }/{#MV\tRef|Var}\tpvalue\tadjust_pvalue\n");
			fprintf(neg_OUTFILE,"#chrom\tcoordinate1\tcoordinate2\t{ #MM\t#MU\t#UM\t#UU }/{#MV\tRef|Var}\tpvalue\tadjust_pvalue\n");
			//char s2t[BATBUF],Dummy[BATBUF];
			int position,Quality;string old_chr;string now_chr; 
			int H=0;//index to genome
			static deque<heterSNP> plus_HeterSNPs;
			static deque<heterSNP> neg_HeterSNPs;
			bool heterSNP_chrom=false;int r;
			printf("\nProcessing sam file \n");
			int coverjoint = 0;
			while( (!bamformat && fgets(s2t,BATBUF,SAM_File)!=0) || (bamformat && (r = samread(bamin, b)) >= 0 ))
			{
				if(bamformat) 
				{
					char *tmp = bam_format1_core( header , b, 0); //2 >>2&3
					strcpy(s2t, tmp);
					free(tmp);
				}
				
				if(s2t[0]=='@') continue;
				//process bar
				nReads++;
				if (nReads % 1000 == 0) {
					Show_Progress_float(nReads);
				}
				alignedread read;
				sscanf(s2t,"%s%d%s%d%d%s%s%s%s%s%s%*[^0-9]%i\tMD:Z:%s",read.readid,&(read.flag),read.chrom,&(read.position),&(read.Quality),read.cigar,Dummy,Dummy,Dummy,\
						read.sequence,read.quality,&(read.mismatches),read.methState);
				if(read.Quality < Quality_cutoff || (read.flag & 0x4) || (read.flag & 0x200)) continue;

				int read_valid_len=Get_read_valid_Length(read.cigar);
				now_chr=read.chrom;
				if(now_chr!=old_chr)
				{	
					if(String_Hash.find(now_chr)!=String_Hash.end())
					{
						H=String_Hash[now_chr];//chromsome index   --  //chromSNPs[H].HeterSNPs
					}else
						continue;
				}
				if(!(read.flag & 0x4)) //match
				{
					if(withFirstResult==0)
					{
				    		old_chr=read.chrom;
				    		withFirstResult=-1;
				   	}
				   	if( !(read.flag & 0x2) && !(read.flag & 0x10) \
				   			|| ( (read.flag & 0x2) && ( ((read.flag & 0x40) && (read.flag & 0x20)) || ((read.flag & 0x80) && (read.flag & 0x10 ))) ) ) //+ strand 
				   	{
				    	while(chromSNPs[H].HeterSNPs.size()>0 && chromSNPs[H].HeterSNPs[0].strand != '+'){
				    		chromSNPs[H].HeterSNPs.pop_front();
				    	}
				    	if(chromSNPs[H].HeterSNPs.size()<2) continue;
				    	
				   		if( pos_first_start_coordinate==-1 || pos_first_end_coordinate==-1)//diff chromsome
				   		{
				   			//if(hetero_SNP && heterSNP_chrom) plus_HeterSNPs=chromSNPs[H].HeterSNPs;
					       	pos_first_start_coordinate = read.position;
					    	pos_first_end_coordinate   = pos_first_start_coordinate + read_valid_len;
				    	}

				    	if(now_chr==old_chr)// the same chromsome
				    	{
				    		
				    		//if(chromSNPs[H].HeterSNPs[0].pos < read.position)
						}
						//deal reads
						while(pos_first_end_coordinate < read.position || now_chr!=old_chr || pos_vDnaMethylMap.size() > memSize) 
	            		{
	            			if(pos_vDnaMethylMap.size() < 4) {
								break;
							}
	            					if( (hetero_SNP && heterSNP_chrom || hetero_meth) && plus_HeterSNPs.size() >0 )
	            						processOneRead_heterSNPs('+',Pvalues, plus_HeterSNPs,pos_vDnaMethylMap, pos_processed_coordinate, pos_potentialCs,pos_validCs,pos_print,pos_OUTFILE);
	            					else
	            						processOneRead('+',Pvalues,pos_vDnaMethylMap, pos_processed_coordinate, pos_potentialCs,pos_validCs,pos_print,pos_OUTFILE);
	            					pos_processed_coordinate = pos_first_end_coordinate - 1;
	            					pos_vDnaMethylMap.pop_front();
	            					removePotentialCs(pos_potentialCs, pos_processed_coordinate);
	            					if(pos_vDnaMethylMap.size()<=0) break;
					              pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
					              pos_first_end_coordinate   = pos_first_start_coordinate + Get_read_valid_Length(pos_vDnaMethylMap[0].cigar);
					              // skip the reads already processed
				                	while(pos_first_end_coordinate - 1 <= pos_processed_coordinate) {
					                    pos_vDnaMethylMap.pop_front();
					                    if(pos_vDnaMethylMap.size()<=0) break;
					                    pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
					                    pos_first_end_coordinate   = pos_first_start_coordinate + Get_read_valid_Length(pos_vDnaMethylMap[0].cigar);
					               	}
					                if( (pos_vDnaMethylMap.size() <=1 || now_chr!=old_chr) && pos_print==1) 
					                {
					                	fprintf(pos_OUTFILE,"******\n");
					                	pos_print=0;
					                }
					                if(now_chr!=old_chr && neg_print==1)
					                {
					                	fprintf(neg_OUTFILE,"******\n");
					                	neg_print=0;
					                }
	            		}
				   	}else if( !(read.flag & 0x2) && (read.flag & 0x10) || \
				   			 ( (read.flag & 0x2) && ( ((read.flag & 0x40) && (read.flag & 0x10)) || ((read.flag & 0x80) && (read.flag & 0x20 ))) ) ) //- strand //- strand
				   	{
				   		if(neg_first_start_coordinate==-1 || neg_first_end_coordinate==-1)
				   		{
				   			if(hetero_SNP && heterSNP_chrom) neg_HeterSNPs=chromSNPs[H].HeterSNPs;
					       	neg_first_start_coordinate = read.position;
					    	neg_first_end_coordinate   = neg_first_start_coordinate + read_valid_len;
				    	}
				    	if(now_chr==old_chr)
				    	{
				    		extractPotentialCs_single(read,Neg_potentialCs,neg_last_coordinate);
					        if(neg_last_coordinate < read.position + read_valid_len) {
					        	neg_last_coordinate = read.position + read_valid_len;
					        }
		            		//add new reads
							//read.strand='-';
							Neg_vDnaMethylMap.push_back(read);
						}
						//deal reads
						while(neg_first_end_coordinate < read.position || now_chr!=old_chr || Neg_vDnaMethylMap.size() > memSize) 
	            		{
							if(Neg_vDnaMethylMap.size() < 4) break;
	            					if(hetero_SNP && heterSNP_chrom && neg_HeterSNPs.size()>0 )
	            					{
	            						processOneRead_heterSNPs('-',Pvalues,neg_HeterSNPs,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_potentialCs,Neg_validCs,neg_print,neg_OUTFILE);
	            					}else
		            					processOneRead('-',Pvalues,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_potentialCs,Neg_validCs,neg_print,neg_OUTFILE);
	            					neg_processed_coordinate = neg_first_end_coordinate - 1;
	            					Neg_vDnaMethylMap.pop_front();
	            					removePotentialCs(Neg_potentialCs, neg_processed_coordinate);
	            					if(Neg_vDnaMethylMap.size()<=0) break;
					              neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
					              neg_first_end_coordinate   = neg_first_start_coordinate + Get_read_valid_Length(Neg_vDnaMethylMap[0].cigar);
					              // skip the reads already processed
				                	while(neg_first_end_coordinate - 1 <= neg_processed_coordinate) {
					                    Neg_vDnaMethylMap.pop_front();
					                    if(Neg_vDnaMethylMap.size()<=0) break;
					                    neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
					                    neg_first_end_coordinate   = neg_first_start_coordinate + Get_read_valid_Length(Neg_vDnaMethylMap[0].cigar);
					                }
					                if( (Neg_vDnaMethylMap.size() <= 1 || now_chr!=old_chr )&& neg_print==1) 
					                {
					                	fprintf(neg_OUTFILE,"******\n");
					                	neg_print=0;
					                }
					                if( now_chr!=old_chr && pos_print==1 )
					                {
					                	fprintf(pos_OUTFILE,"******\n");
					                	pos_print=0;
					                }
	            		}
				   	}
					
            				//diff chromsome --- clear vector
					if( now_chr!=old_chr && withFirstResult==-1)
					{
						pos_vDnaMethylMap.clear();
						pos_potentialCs.clear();
						pos_validCs.clear();
						Neg_vDnaMethylMap.clear();
						Neg_potentialCs.clear();
						Neg_validCs.clear();
						//add new potentialCs
	            		//add new reads and init 
				        pos_first_start_coordinate = -1;
				        pos_first_end_coordinate   = -1;
				        pos_processed_coordinate   = -1;
				        pos_last_coordinate        = -1; 
				        neg_first_start_coordinate = -1; 
				        neg_first_end_coordinate   = -1;
				        neg_processed_coordinate   = -1; 
				        neg_last_coordinate        = -1; 

						if( !(read.flag & 0x2) && !(read.flag & 0x10) \
					   			|| ( (read.flag & 0x2) && ( ((read.flag & 0x40) && !(read.flag & 0x10)) || ((read.flag & 0x80) && (read.flag & 0x10 ))) ) ) //+ strand
						{
					        if(hetero_SNP && heterSNP_chrom)
					        {
					        	plus_HeterSNPs=chromSNPs[H].HeterSNPs;
					        }
					        extractPotentialCs_single(read,pos_potentialCs,pos_last_coordinate);
		            		pos_last_coordinate = read.position + read_valid_len;
					    	pos_first_start_coordinate = read.position;
					    	pos_first_end_coordinate   = pos_first_start_coordinate + read_valid_len;
							//read.strand='+';
							pos_vDnaMethylMap.push_back(read);
						}else if( !(read.flag & 0x2) && (read.flag & 0x10) || \
							   	  ( (read.flag & 0x2) && ( ((read.flag & 0x40) && (read.flag & 0x10)) || ((read.flag & 0x80) && !(read.flag & 0x10 ))) ) ) //- strand
						{
					        if(hetero_SNP && heterSNP_chrom)
					        {
					        	neg_HeterSNPs=chromSNPs[H].HeterSNPs;
					        }
					        extractPotentialCs_single(read,Neg_potentialCs,neg_last_coordinate);
	            			neg_last_coordinate = read.position + read_valid_len;
					       	neg_first_start_coordinate = read.position;
					    	neg_first_end_coordinate   = neg_first_start_coordinate + read_valid_len;
							//read.strand='-';
							Neg_vDnaMethylMap.push_back(read);
						}
					}
				}// end if
				old_chr=now_chr;
			}//end sam
			Show_Progress_float(nReads);
			//sam done!
			//foreach the last  +strand
            while(true) 
            {
            	if(pos_vDnaMethylMap.size() < 4) break;
	        	if(hetero_SNP && heterSNP_chrom && plus_HeterSNPs.size() >0 )
	        		processOneRead_heterSNPs('+',Pvalues,plus_HeterSNPs,pos_vDnaMethylMap, pos_processed_coordinate, pos_potentialCs,pos_validCs,pos_print,pos_OUTFILE);
	        	else
            		processOneRead('+',Pvalues,pos_vDnaMethylMap, pos_processed_coordinate, pos_potentialCs,pos_validCs,pos_print,pos_OUTFILE);
            	pos_processed_coordinate = pos_first_end_coordinate - 1;
            	pos_vDnaMethylMap.pop_front();
            	removePotentialCs(pos_potentialCs, pos_processed_coordinate);
            	if(pos_vDnaMethylMap.size()<=0) break;
				pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
				pos_first_end_coordinate   = pos_first_start_coordinate + Get_read_valid_Length(pos_vDnaMethylMap[0].cigar);
			    // skip the reads already processed
			   	while(pos_first_end_coordinate - 1 <= pos_processed_coordinate) {
					pos_vDnaMethylMap.pop_front();
					if(pos_vDnaMethylMap.size()<=0) break;
			    	pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
				    pos_first_end_coordinate   = pos_first_start_coordinate + Get_read_valid_Length(pos_vDnaMethylMap[0].cigar);
				}
				if(pos_vDnaMethylMap.size() <= 1 && pos_print==1) 
				{
					fprintf(pos_OUTFILE,"******\n");
					pos_print=0;
					break;
				}
            }
            //-strand
            		while(true) 
            		{
            			if(Neg_vDnaMethylMap.size() < 4) break;
	            		if(hetero_SNP && heterSNP_chrom && neg_HeterSNPs.size() >0  || true)
	            			processOneRead_heterSNPs('-',Pvalues,neg_HeterSNPs,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_potentialCs,Neg_validCs,neg_print,neg_OUTFILE);
      				else
            				processOneRead('-',Pvalues,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_potentialCs,Neg_validCs,neg_print,neg_OUTFILE);
            			neg_processed_coordinate = neg_first_end_coordinate - 1;
            			Neg_vDnaMethylMap.pop_front();
            			removePotentialCs(Neg_potentialCs, neg_processed_coordinate);
            			if(Neg_vDnaMethylMap.size()<=0) break;
				neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
				neg_first_end_coordinate   = neg_first_start_coordinate + Get_read_valid_Length(Neg_vDnaMethylMap[0].cigar);
				            // skip the reads already processed
			     	while(neg_first_end_coordinate - 1 <= neg_processed_coordinate) {
					Neg_vDnaMethylMap.pop_front();
					if(Neg_vDnaMethylMap.size()<=0) break;
				       neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
				       neg_first_end_coordinate   = neg_first_start_coordinate + Get_read_valid_Length(Neg_vDnaMethylMap[0].cigar);
				}
				if(Neg_vDnaMethylMap.size() <= 1 && neg_print==1) 
				{
					fprintf(neg_OUTFILE,"******\n");
					neg_print=0;
				}
            		}
            		neg_HeterSNPs.clear();
            		plus_HeterSNPs.clear();
            		
            		fclose(pos_OUTFILE);
            		fclose(neg_OUTFILE);
            		if(!bamformat) fclose(SAM_File);
            		if(bamformat) 
            		{
            			bam_header_destroy(header);
            			bam_destroy1(b);
            			samclose(bamin);
            		}
            		//adjust pvalues
				char pos_Output[100];
				char neg_Output[100];
		       	sprintf(pos_Output,"%s.plus.txt",Output_Name);
		      	FILE* pos_Final_Out=File_Open(pos_Output,"w");
		      	
		      	printf("\nAdjust pvalue ...\n");
		      	if(Pvalues.size()==0)
		      	{
		      		printf("Pvalue.size is 0. \nDone!\n");
		      		exit(0);
		      	}
            		adjust(Pvalues);
            		printf("Print result ...\n");
            		FILE* pos_FILE=File_Open(pos_Output_tmp,"r");
		      		vector<Pvals>::const_iterator cur_iter = Pvalues.begin();

            		while(fgets(s2t,BATBUF,pos_FILE)!=0)
            		{
            			s2t[strlen(s2t)-1]='\0';
            			if(s2t[0]=='*' || s2t[0]=='#') 
            			{
            				fprintf(pos_Final_Out,"%s\n",s2t);
            				continue;
            			}
            			assert(Pvalues.size()>0);
            			while(cur_iter->strand=='-') cur_iter++;
            			fprintf(pos_Final_Out,"%s\t%f\n",s2t,cur_iter->adjust_pval);
            			cur_iter++;
            		}
            		fclose(pos_Final_Out);
            		fclose(pos_FILE);
            		sprintf(neg_Output,"%s.neg.txt",Output_Name);
            		FILE* neg_Final_Out=File_Open(neg_Output,"w");
            		FILE* neg_FILE=File_Open(neg_Output_tmp,"r");
				cur_iter = Pvalues.begin();
            		while(fgets(s2t,BATBUF,neg_FILE)!=0)
            		{
            			s2t[strlen(s2t)-1]='\0';
            			if(s2t[0]=='*' || s2t[0]=='#') 
            			{
            				fprintf(neg_Final_Out,"%s\n",s2t);
            				continue;
            			}
            			assert(Pvalues.size()>0);
            			while(cur_iter->strand=='+') cur_iter++;
            			fprintf(neg_Final_Out,"%s\t%f\n",s2t,cur_iter->adjust_pval);
            			cur_iter++;
            		}
            		fclose(neg_FILE);
            		fclose(neg_Final_Out);
            		remove(pos_Output_tmp);remove(neg_Output_tmp);
            		Pvalues.clear();
			vector<Pvals>().swap(Pvalues);
            		printf("\nDone!\n");
		}//end try
		catch(char* Err)
		{
			printf(Err);
			exit(-1);
		}
	}//end if argc
}//end main

//---------------fisher exact test------------------------
static inline double log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

static double log_hyper_g(const size_t k, const size_t n1, const size_t n2, const size_t t) {
  return (gsl_sf_lnfact(n1) - gsl_sf_lnfact(k) - gsl_sf_lnfact(n1 - k) +
          gsl_sf_lnfact(n2) - gsl_sf_lnfact(t - k) - gsl_sf_lnfact(n2 - (t - k)) -
          (gsl_sf_lnfact(n1 + n2) - gsl_sf_lnfact(t) - gsl_sf_lnfact(n1 + n2 - t)));
}
static double fishers_exact(size_t a, size_t b, size_t c, size_t d) {
  const size_t m = a + c; // sum of first column
  const size_t n = b + d; // sum of second column
  const size_t k = a + b; // sum of first row
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return exp(p);
}

//processed one reads of sam file while with heter SNPs
void processOneRead_heterSNPs(char strand,vector<Pvals> &Pvalues,deque<heterSNP>& HeterSNPs,deque<alignedread> & vDnaMethylMap, int & processed_coordinate,deque<int> & potentialCs,deque<int> & validCs,int& isprint,FILE* OUTFILE) 
{
        // no potentail Cs, no need further processing
        if(potentialCs.size() <= 0) {
            return;
        }
	
        // coverage at potential Cs
        //validCs.clear(); // clear the valid Cs
        while(validCs.size()>0) {
            if(validCs.size()>1) { // keep one validC
                validCs.pop_front();
            }
            else {
                break;
            }
        }
	while(HeterSNPs.size() > 0 && HeterSNPs[0].pos < potentialCs[0] )
	{
		HeterSNPs.pop_front();
	}
        // region to process: from start to end
        int start = processed_coordinate + 1;
        if(start < vDnaMethylMap[0].position ) { // vDnaMethylMap.get(0).getCoordinate() is the smallest available coordinate
            start = vDnaMethylMap[0].position;
        }
        int end = vDnaMethylMap[0].position + Get_read_valid_Length(vDnaMethylMap[0].cigar);

        int Nstat = 7;
        // index 0: M, m - methylated
        // index 1: U, u - un-methylated
        // index 2: A, a - A
        // index 3: C, c - C
        // index 4: G, g - G
        // index 5: T, t - T
        // index 6: others

        char* chrom = vDnaMethylMap[0].chrom;
        //int* coverage = new int[Nstat];
        int coverage[7];
        int processed_snp=0;bool process_snp=false;
        for(int i=0; i < potentialCs.size(); i++) {
            int coordinate = potentialCs[i];
            if(coordinate < start) {
                continue; // skip already processed Cs
            }
            if(coordinate >= end) {
                break; // stop if over the end of the first read
            }
            
	     if(HeterSNPs.size() > processed_snp && HeterSNPs[processed_snp].pos == coordinate) process_snp=true;
	     
	     if(process_snp)
	     {
	     	 validCs.push_back(coordinate);
	     	 process_snp=false;
	     	 processed_snp++;
	     	 continue;
	     }
            // initialize coverage
            for(int j=0; j<Nstat; j++) {
                coverage[j] = 0;
            }
            
            for(int iMap=0; iMap < vDnaMethylMap.size(); iMap++) {
                alignedread read = vDnaMethylMap[iMap];
                int index = coordinate - read.position;
                if(index < 0) {
                    break;
                }
                int RealLen=0;//genome move
                int ReadStep=0;//read move
                if(index < Get_read_valid_Length(read.cigar) ) {
	                for(;ReadStep<strlen(read.methState);ReadStep++)
	                {
	                if(RealLen==index) break;
		        	if(read.methState[ReadStep] == 'S' || read.methState[ReadStep] == 'I')
		        	{
		        		continue;
		        	}else if(read.methState[ReadStep] == 'D')
		        	{
		        		RealLen++;
		        		continue;
		        	}else{
		        		RealLen++;
		        	}
	                }
	                //ReadStep--;
                    char methylationLabel = read.methState[ReadStep];
                    switch(methylationLabel) {
                        case 'M':
                        case 'm':
                        case 'X':
                        case 'H':
                        case 'Z':
                            coverage[0]++;
                            break;
                        case 'U':
                        case 'u':
                        case 'x':
                        case 'h':
                        case 'z':
                            coverage[1]++;
                            break;
                        case 'A':
                        case 'a':
                            coverage[2]++;
                            break;
                        case 'C':
                        case 'c':
                            coverage[3]++;
                            break;
                        case 'G':
                        case 'g':
                            coverage[4]++;
                            break;
                        case 'T':
                        case 't':
                            coverage[5]++;
                            break;
                        default:
                            coverage[6]++;
                            break;
                    } // end of switch
                }
            }// end of iMap
            
            if((coverage[0]>=2) && (coverage[1] >= 2)) {
                validCs.push_back(coordinate);
            }
            
        }//end of potentialCs for
       
        // ### pairing information ###
        // // only one or zero C, no need further processing
        if(validCs.size() <= 1) {
            return;
        }

        int nStat = 4;//int nStatSNP=37;
        int pairStat[nStat]; 
        map<string,int> pairStatSNP;//c++ will initial as 0
        //int pairStatSNP[nStatSNP];
        // pairStat[0]: MM  ||  MR  RM  ||  MV1  V1M  ||  VV
        // pairStat[1]: MU   ||  MV  RU   ||  MV2  V1U   ||  VR
        // pairStat[2]: UM   ||  UR   VM  ||  UV1   V2M  ||  RV
        // pairStat[3]: UU    ||  UV   VU   ||  UV2   V2U   ||  RR
        int pairedState=0;  //  1     2   3  //1  MV meth--ref/var //2  VM  ref/var--meth  //3  VV  var--var
        int processed_snp_validC=0;process_snp=false;
        while(HeterSNPs.size()>0 && validCs[0] > HeterSNPs[0].pos )
        {
        	HeterSNPs.pop_front();
        }
        // just test a valid C with its next valid C
        for(int i=0; i < validCs.size()-1; i++) {
            int coordinate1 = validCs[i];
            int coordinate2 = validCs[i+1];
	     pairedState=0;
	     if(HeterSNPs.size() > processed_snp_validC)
	     {
	     	 if(HeterSNPs[processed_snp_validC].pos == coordinate1) // 2 3
	     	 {
	     	 	 pairedState=2;
	     	 	 if(processed_snp >= 2 && HeterSNPs.size() > processed_snp_validC+1 && HeterSNPs[processed_snp_validC+1].pos == coordinate2)
	     	 	 {
	     	 	 	 pairedState=3;
	     	 	}
	     	 }else if(HeterSNPs[processed_snp_validC].pos == coordinate2)// 1 
	     	 {
	     	 	 pairedState=1;
	     	 }
	     }
	     
            for(int j=0; j<nStat; j++) {
                pairStat[j] = 0; // initialize the pairStat
            }
//
            for(int iMap=0; iMap < vDnaMethylMap.size(); iMap++) {
            	alignedread read = vDnaMethylMap[iMap];
                int index1 = coordinate1 - read.position;
                if(index1 < 0) {
                    break;
                }
                int index2 = coordinate2 - read.position;
                if(index2 < 0) {
                    break;
                }
                
                int RealLen1=0,RealLen2=0;//genome move
                int ReadStep1=0,ReadStep2=0;//read move
                if((index1 < Get_read_valid_Length(read.cigar)  )&&(index2 < Get_read_valid_Length(read.cigar)  )) {
	                for(;ReadStep1<strlen(read.methState);ReadStep1++)
	                {
	                	if(RealLen1==index1) break;
			        	if(read.methState[ReadStep1] == 'S' || read.methState[ReadStep1] == 'I')
			        	{
			        		continue;
			        	}else if(read.methState[ReadStep1] == 'D')
			        	{
			        		RealLen1++;
			        		continue;
			        	}else{
			        		RealLen1++;
			        	}
	                }
	                //ReadStep1--;
                	if(ReadStep2<ReadStep1) 
                	{
                		ReadStep2=ReadStep1;
                		RealLen2=RealLen1-1;
                	}
	                for(;ReadStep2<strlen(read.methState);ReadStep2++)
	                {
			        	if(RealLen2==index2) break;
			        	if(read.methState[ReadStep2] == 'S' || read.methState[ReadStep2] == 'I')
			        	{
			        		continue;
			        	}else if(read.methState[ReadStep2] == 'D')
			        	{
			        		RealLen2++;
			        		continue;
			        	}else{
			        		RealLen2++;
			        	}
	                }
                	//ReadStep2--;
                    char methylationLabel_1 = Char2Trans[ read.methState[ReadStep1] ]; 
                    char methylationLabel_2 = Char2Trans[ read.methState[ReadStep2] ]; 
                    // switch 
                    if(pairedState==0)  //MM meth-meth
                    	processSwicthMM(methylationLabel_1,methylationLabel_2,pairStat);
                    else //1  MV meth--ref/var //2  VM  ref/var--meth  //3  VV  var--var
                    {
                    	if(pairedState==1 && methylationLabel_2=='=') methylationLabel_2=HeterSNPs[processed_snp_validC].RefBase;
                    	else  if(pairedState==2 && methylationLabel_1=='=') methylationLabel_1 =HeterSNPs[processed_snp_validC].RefBase;
                    	else if(pairedState==3)
                    	{
                    		if(methylationLabel_1=='=')
                    			methylationLabel_1 =HeterSNPs[processed_snp_validC].RefBase;
                    		if(methylationLabel_2=='=')
                    			methylationLabel_2 =HeterSNPs[processed_snp_validC+1].RefBase;
                    	}
                    	processSwicthMV(methylationLabel_1,methylationLabel_2,pairStatSNP);
                    }
                }
            } // end of iMap
            
            if(pairedState==0)
            {
	            int count=0;
	            for(int c=0;c<nStat;c++) count+=pairStat[c];
	            if(count > 0)
	            {
			     fprintf(OUTFILE,"%s\t%d\t%d",chrom,coordinate1,coordinate2 );
			     for(int j=0; j<nStat; j++) {
			           fprintf(OUTFILE,"\t%d",pairStat[j]);
			     }
			     double pvalue = fishers_exact(pairStat[0],pairStat[1],pairStat[2],pairStat[3]);
			     fprintf(OUTFILE,"\t%f\n",pvalue);
			     Pvals pval;chrom_offset++;
			     pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
			     Pvalues.push_back(pval);
			     isprint=1;
		    }
		    else
	 	    {
				if( isprint==1) 
				{
					fprintf(OUTFILE,"******\n");
					isprint=0;
				}
	  	   }
	    }else
	    {
				int count_rmM=0,count=0;
				map<string,int>::iterator it;
    			for(it=pairStatSNP.begin();it!=pairStatSNP.end();++it)
    			{
    				string key=it->first;
        			if( ! ((key[0]=='M' || key[0]=='U') && (key[1]=='M' || key[1]=='U')) )
        			{
        				count_rmM+=it->second;
        			}
        			count+=it->second;
        		}
	    	   
	    	   char refBase1=HeterSNPs[processed_snp_validC].RefBase;
	    	   char var1=HeterSNPs[processed_snp_validC].variantBase1;
	    	   char var2=HeterSNPs[processed_snp_validC].variantBase2;
	    	   char refBase2;char var3;char var4;
	    	   if(pairedState==3)
	    	   {
	    	   	   refBase2=HeterSNPs[processed_snp_validC+1].RefBase;
	    	   	   var3=HeterSNPs[processed_snp_validC+1].variantBase1;
	    	   	   var4=HeterSNPs[processed_snp_validC+1].variantBase2;
	    	   }
				if(count_rmM==0 && count > 0)
				{
					fprintf(OUTFILE,"%s\t%d\t%d",chrom,coordinate1,coordinate2 );
					fprintf(OUTFILE,"\t%d",pairStatSNP["MM"]);
					fprintf(OUTFILE,"\t%d",pairStatSNP["MU"]);
					fprintf(OUTFILE,"\t%d",pairStatSNP["UM"]);
					fprintf(OUTFILE,"\t%d",pairStatSNP["UU"]);
					double pvalue = fishers_exact(pairStatSNP["MM"],pairStatSNP["MU"],pairStatSNP["UM"],pairStatSNP["UU"]);
					fprintf(OUTFILE,"\t%f\n",pvalue);
					Pvals pval;
					chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
					Pvalues.push_back(pval);
					isprint=1;
				}else if(count > 0)
	            {
			     fprintf(OUTFILE,"%s\t%d\t%d",chrom,coordinate1,coordinate2 );
			     map<string,int>::iterator itp;
			     for(itp=pairStatSNP.begin();itp!=pairStatSNP.end();++itp)
			     {
			     	 if(itp->second>0) fprintf(OUTFILE,"\t%s:%d",itp->first.c_str(),itp->second);
			     }
			     if(pairedState!=3) //doublevar
			     {
			     		if(HeterSNPs[processed_snp_validC].doubleVar) 
			     		{
		     				char stat_tmp[2];string stat_string; double pvalue;
			     			if(pairedState==1)
			     			{
				     			int MV1,MV2,UV1,UV2;
								stat_tmp[0]='M';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			MV1=pairStatSNP[stat_string];
				     			stat_tmp[0]='M';stat_tmp[1]=var2;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			MV2=pairStatSNP[stat_string];
								stat_tmp[0]='U';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			UV1=pairStatSNP[stat_string];
				     			stat_tmp[0]='U';stat_tmp[1]=var2;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			UV2=pairStatSNP[stat_string];
				     			pvalue = fishers_exact(MV1,MV2,UV1,UV2);
				     			
						     	Pvals pval;
						     	chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
						     	Pvalues.push_back(pval);
			     			}else
			     			{
			     				int V1M,V1U,V2M,V2U;
								stat_tmp[0]=var1;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V1M=pairStatSNP[stat_string];
								stat_tmp[0]=var2;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V2M=pairStatSNP[stat_string];
								stat_tmp[0]=var1;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V1U=pairStatSNP[stat_string];
								stat_tmp[0]=var2;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V2U=pairStatSNP[stat_string];
				     			pvalue = fishers_exact(V1M,V1U,V2M,V2U);
							
							     Pvals pval;
							     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
							     Pvalues.push_back(pval);
			     			}
			     			
			     			fprintf(OUTFILE,"\t%c|%c%c\t%f\n",refBase1,var1,var2,pvalue);
			     		}//end doubleVar
			     		else //single var 
			     		{
		     				char stat_tmp[2];string stat_string; double pvalue;
			     			if(pairedState==1)
			     			{
			     				int MR,MV,UR,UV;
							stat_tmp[0]='M';stat_tmp[1]=refBase1;stat_tmp[2]='\0';stat_string=stat_tmp; 
							MR=pairStatSNP[stat_string];
							stat_tmp[0]='M';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			MV=pairStatSNP[stat_string];
							stat_tmp[0]='U';stat_tmp[1]=refBase1;stat_tmp[2]='\0';stat_string=stat_tmp; 
							UR=pairStatSNP[stat_string];
							stat_tmp[0]='U';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			UV=pairStatSNP[stat_string];
				     			pvalue = fishers_exact(MR,MV,UR,UV);
						     Pvals pval;
						     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
						     Pvalues.push_back(pval);
			     			}else
			     			{
			     				int RM,RU,VM,VU;
							stat_tmp[0]=refBase1;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
							RM=pairStatSNP[stat_string];
							stat_tmp[0]=refBase1;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
							RU=pairStatSNP[stat_string];
							stat_tmp[0]=var1;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
							VM=pairStatSNP[stat_string];
							stat_tmp[0]=var1;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
							VU=pairStatSNP[stat_string];
				     			pvalue = fishers_exact(RM,RU,VM,VU);
						     Pvals pval;
						     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
						     Pvalues.push_back(pval);
			     			}
			     			
			     			fprintf(OUTFILE,"\t%c|%c\t%f\n",refBase1,var1,pvalue);
			     		}
			     }
			     else // VV
			     {
			      	char stat_tmp[2];string stat_string; 
			           if(HeterSNPs[processed_snp_validC].doubleVar && HeterSNPs[processed_snp_validC+1].doubleVar)
			           {
			           	   //V1 V2 || V3 V4
			           	   int V1V3,V1V4,V2V3,V2V4;
						//VV
						stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V1V3=pairStatSNP[stat_string];
						stat_tmp[0]=var1;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V1V4=pairStatSNP[stat_string];
						stat_tmp[0]=var2;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V2V3=pairStatSNP[stat_string];
						stat_tmp[0]=var2;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V2V4+=pairStatSNP[stat_string];
				     		double pvalue = fishers_exact(V1V3,V1V4,V2V3,V2V4);
			    	      	fprintf(OUTFILE,"\t%c|%c%c\t%c|%c%c\t%f\n",refBase1,var1,var2,refBase2,var3,var4,pvalue);
					     Pvals pval;
					     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
					     Pvalues.push_back(pval);
			    	  }else if(!HeterSNPs[processed_snp_validC].doubleVar && !HeterSNPs[processed_snp_validC+1].doubleVar)
			    	  {
			    	  	  int RR,RV,VR,VV;
						stat_tmp[0]=refBase1;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
						RR=pairStatSNP[stat_string];
						stat_tmp[0]=refBase1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						RV=pairStatSNP[stat_string];
						stat_tmp[0]=var1;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
						VR=pairStatSNP[stat_string];
						stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						VV=pairStatSNP[stat_string];
						double pvalue = fishers_exact(RR,RV,VR,VV);
			    	      	  fprintf(OUTFILE,"\t%c|%c\t%c|%c\t%f\n",refBase1,var1,refBase2,var3,pvalue);
					     Pvals pval;
					     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
					     Pvalues.push_back(pval);
			    	      }else if(HeterSNPs[processed_snp_validC].doubleVar && !HeterSNPs[processed_snp_validC+1].doubleVar)
			    	      {
			    	      	  int V1R,V1V3,V2R,V2V3;
						stat_tmp[0]=var1;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V1R=pairStatSNP[ stat_string];
						stat_tmp[0]=var2;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V2R=pairStatSNP[stat_string];
						stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V1V3=pairStatSNP[stat_string];
						stat_tmp[0]=var2;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V2V3=pairStatSNP[stat_string];
						double pvalue = fishers_exact(V1R,V1V3,V2R,V2V3);
			    	      	  fprintf(OUTFILE,"\t%c|%c%c\t%c|%c\t%f\n",refBase1,var1,var2,refBase2,var3,pvalue);
					     Pvals pval;
					     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
					     Pvalues.push_back(pval);
			    	      }else if(!HeterSNPs[processed_snp_validC].doubleVar && HeterSNPs[processed_snp_validC+1].doubleVar)
			    	      {
			    	      	  int RV3,RV4,V1V3,V1V4;
						stat_tmp[0]=refBase1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						RV3=pairStatSNP[stat_string];
						stat_tmp[0]=refBase1;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
						RV4=pairStatSNP[stat_string];
						stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V1V3=pairStatSNP[stat_string];
						stat_tmp[0]=var1;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
						V1V4=pairStatSNP[stat_string];
						double pvalue = fishers_exact(RV3,RV4,V1V3,V1V4);
			    	      	  fprintf(OUTFILE,"\t%c|%c\t%c|%c%c\t%f\n",refBase1,var1,refBase2,var3,var4,pvalue);
					     Pvals pval;
					     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
					     Pvalues.push_back(pval);
			    	      }
			     }
			     isprint=1;
		    }else
		    {
				if( isprint==1) 
				{
					fprintf(OUTFILE,"******\n");
					isprint=0;
				}
		    }
	    }
    	 if(pairedState==2)
	    	     processed_snp_validC++;
	     if(pairedState==3)
	     	     processed_snp_validC++;
	     
	}//end of validCs
	
	//remove processed SNP
	int removeI=0;
        while(HeterSNPs.size()>0) {
            if(HeterSNPs.size()>1 && removeI < processed_snp) { // keep one HeterSNP -1
                HeterSNPs.pop_front();
                removeI++;
            }
            else {
                break;
            }
        }
        //delete[] coverage;
}


void init_pairStatMap(map <string,int>& pairStatMap,map <int,string> & intpairStatMap)
{
	int snp_stat=1;//remove 0(default)
	string pairState="MA";//1-16
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="AM";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="UA";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="AU";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	
	pairState="MT";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="TM";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="UT";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="TU";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	
	pairState="MG";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="GM";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="UG";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="GU";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	
	pairState="MC";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="CM";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="UC";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="CU";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	//---------------------------------17-32
	pairState="AA";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="AT";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="AG";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="AC";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	
	pairState="TA";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="TT";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="TG";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="TC";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	
	pairState="GA";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="GT";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="GG";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="GC";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	
	pairState="CA";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="CT";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="CG";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="CC";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	
	//33-36
	pairState="MM";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="MU";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="UM";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
	pairState="UU";
	intpairStatMap[snp_stat]=pairState; pairStatMap[pairState]=snp_stat++;
}

//processed one reads of sam file
void processOneRead(char strand,vector<Pvals> &Pvalues, deque<alignedread> & vDnaMethylMap, int & processed_coordinate,deque<int> & potentialCs,deque<int> & validCs,int& isprint,FILE* OUTFILE) {
	// no potentail Cs, no need further processing
        if(potentialCs.size() <= 0) {
            return;
        }
	
        // coverage at potential Cs
        //validCs.clear(); // clear the valid Cs
        while(validCs.size()>0) {
            if(validCs.size()>1) { //remove = //because 
                validCs.pop_front();
            }
            else {
                break;
            }
        }

        // region to process: from start to end
        int start = processed_coordinate + 1;
        if(start < vDnaMethylMap[0].position ) { // vDnaMethylMap.get(0).getCoordinate() is the smallest available coordinate
            start = vDnaMethylMap[0].position;
        }
        int end = vDnaMethylMap[0].position + Get_read_valid_Length(vDnaMethylMap[0].cigar);

        int Nstat = 7;
        // index 0: M, m - methylated
        // index 1: U, u - un-methylated
        // index 2: A, a - A
        // index 3: C, c - C
        // index 4: G, g - G
        // index 5: T, t - T
        // index 6: others

        char* chrom = vDnaMethylMap[0].chrom;
        //int* coverage = new int[Nstat];
        int coverage[7];
        for(int i=0; i < potentialCs.size(); i++) {
            int coordinate = potentialCs[i];
            if(coordinate < start) {
                continue; // skip already processed Cs
            }
            if(coordinate >= end) {
                break; // stop if over the end of the first read
            }

            // initialize coverage
            for(int j=0; j<Nstat; j++) {
                coverage[j] = 0;
            }
        
            for(int iMap=0; iMap < vDnaMethylMap.size(); iMap++) {
                alignedread read = vDnaMethylMap[iMap];
                int index = coordinate - read.position;
                if(index < 0) {
                    break;
                }
		int RealLen=0;//genome move
		int ReadStep=0;//read move
                if(index < Get_read_valid_Length(read.cigar) ) {
	                for(;ReadStep<strlen(read.methState);ReadStep++)
	                {
		        	if(read.methState[ReadStep] == 'S' || read.methState[ReadStep] == 'I')
		        	{
		        		continue;
		        	}else if(read.methState[ReadStep] == 'D')
		        	{
		        		RealLen++;
		        		continue;
		        	}else{
		        		RealLen++;
		        	}
		        	if(RealLen-1==index) break;
	                }
	                //ReadStep--;
                    char methylationLabel = read.methState[ReadStep];
                    switch(methylationLabel) {
                        case 'M':
                        case 'm':
                        case 'X':
                        case 'H':
                        case 'Z':
                            coverage[0]++;
                            break;
                        case 'U':
                        case 'u':
                        case 'x':
                        case 'h':
                        case 'z':
                            coverage[1]++;
                            break;
                        case 'A':
                        case 'a':
                            coverage[2]++;
                            break;
                        case 'C':
                        case 'c':
                            coverage[3]++;
                            break;
                        case 'G':
                        case 'g':
                            coverage[4]++;
                            break;
                        case 'T':
                        case 't':
                            coverage[5]++;
                            break;
                        default:
                            coverage[6]++;
                            break;
                    } // end of switch
                }
            }// end of iMap

  /*
            printf("%s %d",chrom,coordinate);
            for(int j=0; j<Nstat; j++) {
                printf("\t%d",coverage[j]);
            }
            printf("\n");
*/
            if((coverage[0]>=2) && (coverage[1] >= 2)) {
                validCs.push_back(coordinate);
            }
            
        }//end of potentialCs for
       
        // ### pairing information ###
        // // only one or zero C, no need further processing
        if(validCs.size() <= 1) {
            return;
        }

        int nStat = 4;
        //int* pairStat = new int[nStat];
        int pairStat[nStat];
        // pairStat[0]: MM
        // pairStat[1]: MU
        // pairStat[2]: UM
        // pairStat[3]: UU

        // just test a valid C with its next valid C
        for(int i=0; i < validCs.size()-1; i++) {
            int coordinate1 = validCs[i];
            int coordinate2 = validCs[i+1];

            for(int j=0; j<nStat; j++) {
                pairStat[j] = 0; // initialize the pairStat
            }
            for(int iMap=0; iMap < vDnaMethylMap.size(); iMap++) {
            	alignedread read = vDnaMethylMap[iMap];
                int index1 = coordinate1 - read.position;
                if(index1 < 0) {
                    break;
                }
                int index2 = coordinate2 - read.position;
                if(index2 < 0) {
                    break;
                }
                
		int RealLen1=0,RealLen2=0;//genome move
		int ReadStep1=0,ReadStep2=0;//read move
              if((index1 < Get_read_valid_Length(read.cigar)  )&&(index2 < Get_read_valid_Length(read.cigar)  )) {
	                for(;ReadStep1<strlen(read.methState);ReadStep1++)
	                {
		        	if(read.methState[ReadStep1] == 'S' || read.methState[ReadStep1] == 'I')
		        	{
		        		continue;
		        	}else if(read.methState[ReadStep1] == 'D')
		        	{
		        		RealLen1++;
		        		continue;
		        	}else{
		        		RealLen1++;
		        	}
		        	if(RealLen1-1==index1) break;
	                }
	                //ReadStep1--;
                	if(ReadStep2<ReadStep1) 
                	{
                		ReadStep2=ReadStep1;
                		RealLen2=RealLen1-1;
                	}
	                for(;ReadStep2<strlen(read.methState);ReadStep2++)
	                {
		        	if(read.methState[ReadStep2] == 'S' || read.methState[ReadStep2] == 'I')
		        	{
		        		continue;
		        	}else if(read.methState[ReadStep2] == 'D')
		        	{
		        		RealLen2++;
		        		continue;
		        	}else{
		        		RealLen2++;
		        	}
		        	if(RealLen2-1==index2) break;
	                }
                	//ReadStep2--;
                    char methylationLabel_1 = read.methState[ReadStep1]; 
                    char methylationLabel_2 = read.methState[ReadStep2]; 

                    switch(methylationLabel_1) {
                        case 'M':
                        case 'm':
                        case 'X':
                        case 'H':
                        case 'Z':
                            switch(methylationLabel_2) {
                                case 'M':
                                case 'm':
	                        case 'X':
	                        case 'H':
	                        case 'Z':
                                    pairStat[0]++;
                                    break;
                                case 'U':
                                case 'u':
	                        case 'x':
	                        case 'h':
	                        case 'z':
                                    pairStat[1]++;
                                    break;
                            }
                            break;
                        case 'U':
                        case 'u':
                        case 'x':
                        case 'h':
                        case 'z':
                            switch(methylationLabel_2) {
                                case 'M':
                                case 'm':
	                        case 'X':
	                        case 'H':
	                        case 'Z':
                                    pairStat[2]++;
                                    break;
                                case 'U':
                                case 'u':
	                        case 'x':
	                        case 'h':
	                        case 'z':
                                    pairStat[3]++;
                                    break;
                            }
                            break;
                        default:
                            break;
                    } // end of switch
                }
            } // end of iMap
            int count=0;
            for(int c=0;c<nStat;c++) count+=pairStat[c];
            if(count > 0)
            {
		     fprintf(OUTFILE,"%s\t%d\t%d",chrom,coordinate1,coordinate2 );
		     for(int j=0; j<nStat; j++) {
		           fprintf(OUTFILE,"\t%d",pairStat[j]);
		     }
		     double pvalue = fishers_exact(pairStat[0],pairStat[1],pairStat[2],pairStat[3]); //MM+UU  ----  MU+UM
		     fprintf(OUTFILE,"\t%f\n",pvalue);
		     Pvals pval;
		     chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
		     Pvalues.push_back(pval);
		     isprint=1;
	    }else
	    {
			if( isprint==1) 
			{
				fprintf(OUTFILE,"******\n");
				isprint=0;
			}
	    }
	}//end of validCs
	//delete[] coverage;
}


//{----------------------------------- FILE HANDLING ---------------------------------------------------------

//sort meth and variant
bool sort_methvar(const heterSNP & h1, const heterSNP & h2){
	return h1.pos < h2.pos;
}
//remove same pos and modify methv
bool uniq_and_mod(std::deque< heterSNP > & v){
    if(v.size() <= 1)
        return 0;
    std::deque< heterSNP >::iterator it;
    std::deque< heterSNP >::iterator last = v.begin();
    for( it = v.begin()+1; it != v.end(); )
    {
        if( it->pos == last->pos ){
        	last -> methv=2;
            it = v.erase( it );
        }
        else
        {
        	last = it;
            ++it;
        }
    }
}
//-------qvalue
static bool sort_raw_pval(const Pvals &r1, const Pvals &r2) {
  return r1.pvalue < r2.pvalue;
}

static bool sorted_by_coordinate(const Pvals &r1, const Pvals &r2) {
  return r1.coordinate < r2.coordinate;
}
////BH
void adjust(vector<Pvals> &Pvalues) {

      sort(Pvalues.begin(), Pvalues.end(), sort_raw_pval);

      for (size_t ind = 0; ind < Pvalues.size(); ++ind) 
      {
        const double current_score = Pvalues[ind].pvalue;
        //Assign a new one.
        const double adjust_pval = (double(Pvalues.size()))*current_score/(double(ind + 1));
        Pvalues[ind].adjust_pval = adjust_pval;
      }

      for (vector<Pvals>::reverse_iterator it = Pvalues.rbegin() + 1; it != Pvalues.rend(); ++it) 
      {
        const Pvals &prev_locus = *(it - 1);
        Pvals &cur_locus = *(it);
        cur_locus.adjust_pval = std::min(prev_locus.adjust_pval, cur_locus.adjust_pval);
      }

      for (vector<Pvals>::iterator it = Pvalues.begin();it != Pvalues.end(); ++it) 
      {
        Pvals &cur_locus = *(it);
        if (cur_locus.adjust_pval > 1.0)
          cur_locus.adjust_pval = 1.0;
      }

      // Restore original order
      sort(Pvalues.begin(), Pvalues.end(), sorted_by_coordinate);

}
//}----------------------------------- FILE HANDLING ---------------------------------------------------------
void extractPotentialCs_single(alignedread read,deque<int> & potentialCs,int & last_coordinate) {
        int coordinate = read.position;
        char* methylationStatus = read.methState;
        int length = strlen(methylationStatus);
        int start = last_coordinate;
        if(start < coordinate) {
            start = coordinate;
        }
        int RealLen=0;//genome move steps
        for(int i=0; i<length; i++) {
        	if(methylationStatus[i] == 'S' || methylationStatus[i] == 'I')
        	{
        		continue;
        	}else if(methylationStatus[i] == 'D')
        	{
        		RealLen++;
        		continue;
        	}else{
        		RealLen++;
        	}
        	if(coordinate + RealLen -1< start) continue;
        	if( methylationStatus[i] != '=' && Char2Trans[methylationStatus[i]] !='=') 
        	{
        		potentialCs.push_back(coordinate+RealLen-1);
        	}
        }
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
//process meth file
int process_meth(char* buffer,heterSNP* var_t, char* chrom, int NMETH, int NCOVER, float MFloat){
    int i = 0, j = 0, s = 0, e = 0;

    var_t->methv = 1;

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    //chrom = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) chrom[j - s] = buffer[j];
    chrom[j - s] = '\0';

    char* tempstring;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    var_t->pos = atoi(tempstring);
    free(tempstring);

    tempstring = NULL; //strand + -
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
  /*
    if(meth_strand != '.' && (meth_strand == '+' || meth_strand == '-' ) ) {
        if(meth_strand != tempstring[0]) {
          free(tempstring);
          return -1;
        }
    }
  */
    var_t->strand = tempstring[0];

//    var_t->RefBase = 'M'; 
//    var_t->variantBase1 = 'U'; 
    free(tempstring);
    
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;

    tempstring = NULL;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1); 
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    int nmeth = atoi(tempstring);
    free(tempstring);
    if( nmeth < NMETH)
      return -1;

    tempstring = NULL;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    int ncover = atoi(tempstring);
    free(tempstring);
    if( ncover < NCOVER) 
        return -1;

    if(!( ( (float)nmeth/ncover >= MFloat && (float)nmeth/ncover <= (float) (1.0 - MFloat) ) && ncover >= NCOVER && nmeth >= NMETH )) {
        return -1;
    }

    return 1;
}
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<map>
#include <vector>
#include <deque>
#ifdef __cplusplus
extern "C"
{
#endif
	#include "./samtools-0.1.18/sam_header.h"
	#include "./samtools-0.1.18/sam.h"
#ifdef __cplusplus
}
#endif
#include "processPairedBlock.hpp"
#define memSize 50000
#define MAX_Len_ALLOWED 100
using std::map;

void extractPotentialCs(alignedread read,deque<int> & potentialCs,int & last_coordinate);
bool twopoint_Comp(twopoint a,twopoint b);
void paiedend_methyhaplo(bool bamformat,char* Align_fileName,char* Output_Name);
void adjust(vector<Pvals> &Pvalues);
int processbamread(const bam_header_t *header, const bam1_t *b, char* readid,int &Flag,char* Chrom,int &pos,int &mapQuality,char* CIG,char* Chrom_P,int &pos_P,int &Insert_Size,char* forReadString,char* forQuality, char* methState);
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

void processOneRead_heterSNPs_paired(char strand,vector<Pvals> &Pvalues,deque<heterSNP>& HeterSNPs,deque<alignedread> & vDnaMethylMap, int & processed_coordinate,deque<int> & validCs,int& isprint,FILE* OUTFILE);
void processOneRead_paired(deque<twopoint> jointposition,char strand,vector<Pvals> &Pvalues, deque<alignedread> & vDnaMethylMap, int & processed_coordinate,deque<int> & potentialCs,deque<int> & validCs,int& isprint,deque <validC> &paired_validCs);
void removePotentialCs(deque<int> & potentialCs, int & processed_coordinate);
void repalceReadmeth(char* methState);

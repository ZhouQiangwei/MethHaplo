#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include <vector>
#include <deque>
#define CHROM_NAME_LEN 40
#define MAX_Len_ALLOWED 100

using std::vector;using std::string;using std::deque;

struct heterSNP
{
	char RefBase;
	char variantBase1;
	char variantBase2;
	bool doubleVar;
	int pos;
//	unsigned int methv; //0 , 1 -- meth, 2 -- meth and variant
//	char strand;
};
struct GenomeSNP
{
	char chrom[CHROM_NAME_LEN];
	deque<struct heterSNP> HeterSNPs;
	int index;
};
struct GenomeMETH
{
        char chrom[CHROM_NAME_LEN];
        deque<int> NegHeterMeths;
	deque<int> PlusHeterMeths;
        int index;
};
struct Pvals {
	unsigned coordinate;
	double pvalue;
	double adjust_pval;
	char strand;
};

struct pairedBLOCK
{
	int start;int end;
	int index;//sort paired Block
	std::vector<std::string> coordinate_expand;
	
};

struct twopoint
{
	char chrom[MAX_Len_ALLOWED];
	int coordinate1;int coordinate2;
	std::string jointstate;
	double pvalue;
	char refBase;
	char var1;char var2;
	double isdoubleVar;
	char refBase2;
	char var3;char var4;
	double isdoubleVar2;
	int pairedState;//1  MV meth--ref/var //2  VM  ref/var--meth  //3  VV  var--var // 0 MM
};

struct pairedDNAmethyMap
{
	//vector<alignread> paired_reads;//methstate = left"+"right
	int left_start;int left_end;
	int right_start;int right_end;
	std::deque<twopoint> jointposition;
	bool mateBlock;
};


struct OnlypairedMap
{
	//vector<alignread> paired_reads;//methstate = left"+"right
	int start;int end;
	std::deque<twopoint> jointposition;
	
};


struct validC
{
	int coordinate;
	bool isSNP;
	bool pairedValid;//only paired-end reads
	char RefBase;
	char variantBase1;
	char variantBase2;
	bool doubleVar;
};


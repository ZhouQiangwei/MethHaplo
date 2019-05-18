#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <map>
#include <vector>
#include <deque>
#include "paired.hpp"


#define BATBUF 3000

#define memSize 50000

using std::vector;using std::string;using std::map;
using std::deque;

//{-----------------------------  FUNCTION PRTOTYPES  -----------------------

void processOneRead(char strand,vector<Pvals> &Pvalues, deque<alignedread> & vDnaMethylMap, deque<int> & validCs,int& isprint,FILE* OUTFILE);
void processOneRead_heterSNPs(char strand,vector<Pvals> &Pvalues, deque<heterSNP>& HeterSNPs,deque<alignedread> & vDnaMethylMap, deque<int> & validCs,int& isprint,FILE* OUTFILE);
void removePotentialCs(deque<int> & potentialCs, int & processed_coordinate);
void processSwicthMM(char& methylationLabel_1,char& methylationLabel_2,int* pairStat);
void processSwicthMV(char& methylationLabel_1,char& methylationLabel_2,std::map<string,int> & pairStatSNP);
void init_pairStatMap(map <string,int>& pairStatMap,std::map <int,string> & intpairStatMap);
void adjust(vector<Pvals> &Pvalues);

//--------------------- hypergeomtric distribution
/* gsl_cdf_hypergeometric_Q(k, n1, n2, t); //lower.tail = FALSE  || gsl_cdf_hypergeometric_P(k, n1, n2, t);  //lower.tail = TRUE
* fishers_exact(MM,MU,UM,UU)==gsl_cdf_hypergeometric_Q(MM-1,MM+MU,UM+UU,MM+UM)
* 
*/
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

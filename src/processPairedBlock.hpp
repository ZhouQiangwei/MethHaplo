#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <map>
#include <deque>
#include <time.h>
#include <math.h>
#include <algorithm>
#include "common.hpp"
#include "alignRead.h"
#include <gsl/gsl_cdf.h>
#include <tr1/cmath>
#include <gsl/gsl_sf_gamma.h>
#include<sstream>

using std::deque;
using std::string;using std::map;
void print_block(deque<pairedDNAmethyMap> & paired_DnaMethylBlock,deque<twopoint> Onlypairedjointposition,int processed_coordinate,FILE* outFile);
void process_paired_block(deque<alignedread> & paired_vDnaMethylMap,deque<pairedDNAmethyMap>::iterator cur_iter,deque<pairedDNAmethyMap> &pos_paired_DnaMethylBlock\
		,deque<validC> & paired_validCs,deque<twopoint>& Onlypairedjointposition);

int Get_read_valid_Length(char* cig);
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_log(const char *log);
void Show_Progress(long & nReads);
void Show_Progress_SNP(long & nReads);
void Show_Progress_Meth(long & nReads);
void Show_Progress_float(long & nReads);
bool twopoint_Comp(twopoint a,twopoint b);
void processSwicthMM(char& methylationLabel_1,char& methylationLabel_2,int* pairStat);
void processSwicthMV(char& methylationLabel_1,char& methylationLabel_2, map<string,int> & pairStatSNP);
std::string int2str(int &int_temp);

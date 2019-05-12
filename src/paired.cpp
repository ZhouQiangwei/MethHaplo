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
#include <time.h>
#include <math.h>
#include <algorithm>
#include <limits.h>
#include "paired.hpp"

#define BATBUF 3000
using std::vector;using std::string;using std::map;using std::deque;
extern int test;
extern map <string,int> pairStatMap;
extern map <int,string> intpairStatMap;
extern bool hetero_SNP;
extern int Quality_cutoff;
extern int minIS;
extern int maxIS;
extern char Char2Trans[];
extern map <string,int> String_Hash;
extern GenomeSNP chromSNPs[];
extern GenomeMETH chromMeths[];
extern unsigned chrom_offset;
extern int readDBSIZE;

int filterRead_pair(deque<heterSNP> & HeterSNPs, deque<int> & HeterMeths, alignedread & read, int read_valid_len);
void sortdequebyvec(deque<alignedread> & pos_vDnaMethylMap);
int Get_read_valid_Length_and_is(char* cig, int & containIS);

void removePotentialCs(deque<int> & potentialCs, int & processed_coordinate) {
        while(potentialCs.size()>0) {
            if(potentialCs[0] <= processed_coordinate) { //remove = //because 
                potentialCs.pop_front();
            }
            else {
                break;
            }
        }
}

void extractPotentialCs(alignedread read,deque<int> & potentialCs,int & last_coordinate) {
        int coordinate = read.position;
        char* methylationStatus = read.methState;
        int length = strlen(methylationStatus);
        int start = last_coordinate;
        if(start < coordinate) {
            start = coordinate;
        }
        int RealLen=0;//genome move steps
        bool thisRight=false;
        if(length>0)
        {
        	if(methylationStatus[length-1]!='+')
        	{
        		thisRight=true;
        	}
        }else 
        	return;
        bool Right_begin=false;
        for(int i=0; i<length; i++) {
        	if(thisRight && !Right_begin)
        	{
        		if(methylationStatus[i] == '+')
        		{
        				Right_begin=true;
        		}
        		continue;
        	}
    		if(!thisRight && methylationStatus[i] == '+')
    		{
    				break;
    		}
    		else if(methylationStatus[i] == 'S' || methylationStatus[i] == 'I')
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

bool Paired_align_Comp(alignedread a,alignedread b)
{
	int mina=a.position>a.mateposition?a.mateposition:a.position;
	int minb=b.position>b.mateposition?b.mateposition:b.position;
    return mina < minb;
}

void FreeQ_twopint(std::deque <twopoint>  & t ) 
{
	t.clear();
	//std::vector <twopoint> tmp; 
	//std::swap(t,tmp);
}
void FreeQ_paired_vDnaMethylMap(std::deque <alignedread>  & t ) 
{
	t.clear();
	//std::vector <alignedread> tmp; 
	//std::swap(t,tmp);
}
void FreeQ_paired_validCs(std::deque <validC>  & t ) 
{
	t.clear();
	//std::vector <validC> tmp; 
	//std::swap(t,tmp);
}
void FreeQ_DnaMethylBlock(std::deque <pairedDNAmethyMap> & t)
{
	t.clear();
	//std::vector <pairedDNAmethyMap> tmp; 
	//std::swap(t,tmp);
}
static bool vDna_compare(const alignedread& a,const alignedread& b){
    return a.position < b.position;

    int mina;
    int minb;
    if(strcmp(a.matechrom, "*") == 0) mina = a.position;
    else mina =a.position>a.mateposition?a.mateposition:a.position;
    if(strcmp(b.matechrom, "*") == 0) minb = b.position;
    else minb=b.position>b.mateposition?b.mateposition:b.position;
    return mina < minb;
}

void store_s_read(deque<alignedread> & pos_vDnaMethylMap, deque<alignedread> & neg_vDnaMethylMap, std::map<std::string, alignedread> & pos_readDB,  std::map<std::string, alignedread> & neg_readDB, alignedread & read, int & pos_smallerrightpos, int & neg_smallerrightpos, bool newchr){
    if(pos_readDB.size() == 0 && neg_readDB.size() == 0) return;

    if( ((read.flag & 0x40) && (read.flag & 0x20)) || ((read.flag & 0x80) && (read.flag & 0x10 )) ){ //plus
        if(read.position > pos_smallerrightpos && pos_smallerrightpos!=0 ){
            for (std::map<string, alignedread>::iterator iter = pos_readDB.begin(); iter != pos_readDB.end(); ) {
                if ((iter->second).mateposition < read.position && (iter->second).position < read.position) {
                    strcpy((iter->second).matechrom, "*");
                    pos_vDnaMethylMap.push_back(iter->second);
                    pos_readDB.erase(iter++);
                }else break; //++iter;
            }
            if(pos_readDB.size()>0 && strcpy((pos_readDB.begin()->second).matechrom, "=")==0) pos_smallerrightpos = pos_readDB.begin()->second.mateposition;
            else pos_smallerrightpos=0;
        }
    }else if( ((read.flag & 0x40) && (read.flag & 0x10)) || ((read.flag & 0x80) && (read.flag & 0x20 )))//4 2 //neg 
    {
        if(neg_smallerrightpos!=0 && read.position > neg_smallerrightpos){
            for (std::map<string, alignedread>::iterator iter = neg_readDB.begin(); iter != neg_readDB.end(); ) {
                if ((iter->second).mateposition < read.position && (iter->second).position < read.position) {
                    strcpy((iter->second).matechrom, "*");
                    neg_vDnaMethylMap.push_back(iter->second);
                    neg_readDB.erase(iter++);
                }else break; //++iter;
            }
            if(neg_readDB.size()>0 && strcpy((neg_readDB.begin()->second).matechrom, "=")==0) neg_smallerrightpos = neg_readDB.begin()->second.mateposition;
            else neg_smallerrightpos=0;
        }
    }
    if(newchr){
        for (std::map<string, alignedread>::iterator iter = pos_readDB.begin(); iter != pos_readDB.end(); ) {
                strcpy((iter->second).matechrom, "*");
                pos_vDnaMethylMap.push_back(iter->second);
                pos_readDB.erase(iter++);
        }
	pos_smallerrightpos = 0;
        for (std::map<string, alignedread>::iterator iter = neg_readDB.begin(); iter != neg_readDB.end(); ) {
                strcpy((iter->second).matechrom, "*");
                neg_vDnaMethylMap.push_back(iter->second);
                neg_readDB.erase(iter++);
        }
	neg_smallerrightpos = 0;
    }

}

void copy_read(alignedread & ra, alignedread & rb){
    ra.position = rb.position;
    //ra.chrom = rb.chrom; //same chrom
    ra.real_len = rb.real_len;
    strcpy(ra.cigar, rb.cigar);
    //strcpy(ra.methState, rb.methState);
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

void paiedend_methyhaplo(bool bamformat,char* Align_fileName,char* Output_Name)
{
	char pos_Output_tmp[100];
	char neg_Output_tmp[100];
	
	long nReads=0;
	FILE* SAM_File;
	samfile_t *bamin = 0;bam1_t *b;bam_header_t *header;
	if(bamformat)
	{
		if ((bamin = samopen(Align_fileName, "rb", 0)) == 0) {
			fprintf(stderr, "fail to open \"%s\" for reading.\n", Align_fileName);
		}
		b = bam_init1();
		header=bam_header_dup((const bam_header_t*)bamin->header);
		Show_log("Processing bam file ");
	}
	else
	{
		SAM_File=File_Open(Align_fileName,"r");
		Show_log("Processing sam file ");
	}
	
	static deque<int> pos_validCs;
	static deque<int> Neg_validCs;
	
	//paired-end
	static deque<validC> plus_paired_validCs;
	static deque<validC> Neg_paired_validCs;
	
	static deque<alignedread> pos_vDnaMethylMap;
	static deque<alignedread> Neg_vDnaMethylMap;
	
	//pvalues
	vector<Pvals> Pvalues;
	//read sam read mapping
    std::map <std::string,alignedread> pos_readDB;
    std::map <std::string,alignedread> neg_readDB;

	int pos_first_start_coordinate = -1; //+ strand
	int pos_first_end_coordinate   = -1;
	int pos_processed_coordinate   = -1; // the last coordinate already processed
	int neg_first_start_coordinate = -1; // - strand
	int neg_first_end_coordinate   = -1;
	int neg_processed_coordinate   = -1; // the last coordinate already processed
	
	int withFirstResult = 0; // indicator whether the first mapped resul is read
	//processed the read number
	nReads = 1;
	int pos_print=0,neg_print=0;
	//open output file and write the results
	sprintf(pos_Output_tmp,"%s.tmp.plus.txt",Output_Name);
	sprintf(neg_Output_tmp,"%s.tmp.neg.txt",Output_Name);
	FILE* pos_OUTFILE=File_Open(pos_Output_tmp,"w");
	FILE* neg_OUTFILE=File_Open(neg_Output_tmp,"w");
	fprintf(pos_OUTFILE,"#chrom\tcoordinate1\tcoordinate2\t{ #MM\t#MU\t#UM\t#UU }/{#MV\tRef|Var}\tpvalue\n");
	fprintf(neg_OUTFILE,"#chrom\tcoordinate1\tcoordinate2\t{ #MM\t#MU\t#UM\t#UU }/{#MV\tRef|Var}\tpvalue\n");
	char s2t[BATBUF],Dummy[BATBUF];
	int position,Quality;string old_chr;string now_chr;string printedchr;
	int H=0;//index to genome
	static deque<heterSNP> plus_HeterSNPs;
	static deque<heterSNP> neg_HeterSNPs;
	bool heterSNP_chrom=false;int r; int prevloc=0;int nowloc = 0;int containIS = 0;
	int pos_smallerrightpos = 0;
	int neg_smallerrightpos = 0;
	char printinf[100];
	alignedread read;
	while( (!bamformat && fgets(s2t,BATBUF,SAM_File)!=0) || (bamformat && (r = samread(bamin, b)) >= 0 ))
	{
		if(s2t[0]=='@') continue;
		if(bamformat) 
		{
			r = processbamread(header, b, read.readid,read.flag,read.chrom,read.position,read.Quality,read.cigar,read.matechrom,read.mateposition,read.IS,read.sequence,read.quality, read.methState);
			if(r == -1) continue;
		}
		
		//process bar
		nReads++;
		if (nReads % 1000000 == 0) {
			printf("\ncoord: %d %d ||size: %d %d ||readpos: %d %d ||DB: %d %d %d %d\n", pos_first_end_coordinate, neg_first_end_coordinate, pos_vDnaMethylMap.size(), Neg_vDnaMethylMap.size(), read.position, read.flag, pos_readDB.size(), neg_readDB.size(), pos_smallerrightpos, neg_smallerrightpos);
			Show_Progress(nReads);
		}
		if(!bamformat) sscanf(s2t,"%s%d%s%d%d%s%s%d%d%s%s%*[^0-9]%i\tMD:Z:%s",read.readid,&(read.flag),read.chrom,&(read.position),&(read.Quality),read.cigar,\
				read.matechrom,&(read.mateposition),&(read.IS),read.sequence,read.quality,&(read.mismatches),read.methState);
        
		if(read.Quality < Quality_cutoff || (read.flag & 0x4) || (read.flag & 0x200)) continue;
		if(String_Hash.find(read.chrom)==String_Hash.end()) continue;
	
		if(read.chrom!=old_chr) H=String_Hash[read.chrom];
		if( ((read.flag & 0x40) && (read.flag & 0x20)) || ((read.flag & 0x80) && (read.flag & 0x10 ))){
			if(chromMeths[H].PlusHeterMeths.size() + chromSNPs[H].HeterSNPs.size() <2) continue;
		}else if( ((read.flag & 0x40) && (read.flag & 0x10)) || ((read.flag & 0x80) && (read.flag & 0x20 ))){
			if(chromMeths[H].NegHeterMeths.size() + chromSNPs[H].HeterSNPs.size() <2) continue;
		}

		nowloc=read.position;
		read.span = read.IS>0?read.IS:-1*read.IS;
        int read_valid_len=Get_read_valid_Length_and_is(read.cigar, containIS);
        if(containIS)
    		repalceReadmeth(read.methState);
		//if(read.span < read.real_len) ... need merge
		read.real_len=read_valid_len;
		now_chr=read.chrom;
        if(strcmp(read.matechrom, "=") != 0 || read.span > maxIS ) strcpy(read.matechrom, "*");
		if((read.flag & 0x2)){
            if(now_chr != old_chr)
                    store_s_read(pos_vDnaMethylMap, Neg_vDnaMethylMap, pos_readDB, neg_readDB, read, pos_smallerrightpos, neg_smallerrightpos, true);
            else
                    store_s_read(pos_vDnaMethylMap, Neg_vDnaMethylMap, pos_readDB, neg_readDB, read, pos_smallerrightpos, neg_smallerrightpos, false);
        	}

		if( (read.flag & 0x2) &&  !strcmp("=",read.matechrom) && read.span <=maxIS && read.span>=minIS && !(read.flag & 0x200) && !(read.flag & 0x100) && !(read.flag & 0x800) )//proper paired
		{
			if((read.flag & 0x40) && (read.flag & 0x20)) /*first and plus, map:1*/   //the read is the first read(1) in a pair or second reads(2),smaller position
			{
                strcat(read.methState, "+");
				pos_readDB[read.readid]= read;
                if(pos_smallerrightpos == 0) pos_smallerrightpos = read.mateposition;
                continue;
			}else if((read.flag & 0x80) && (read.flag & 0x20 )){ //if((read.flag & 0x40) && (read.flag & 0x10 )){ /*second and plus, map:2*/
                strcat(read.methState, "+");
                neg_readDB[read.readid] = read;
                if(neg_smallerrightpos == 0) neg_smallerrightpos = read.mateposition;
                continue;
            }
			else //the read align position is bigger
			{
				map<string ,alignedread >::iterator l_it;
				l_it=pos_readDB.find(read.readid);
			    if(l_it!=pos_readDB.end())
			    {
			    	string tmp=(l_it->second).methState;
                    tmp = tmp + read.methState;
			    	tmp.copy(read.methState,tmp.length(),0);
			    	*(read.methState+tmp.length())='\0';
					
                    read.mateposition = read.position;
                    read.position = (l_it->second).position;
                    read.mate_real_len = read.real_len;
                    read.real_len = (l_it->second).real_len;
                    pos_readDB.erase(l_it);
			    }else if( (l_it=neg_readDB.find(read.readid))!= neg_readDB.end()){
                    string tmp=(l_it->second).methState;
                    tmp = tmp + read.methState;
                    tmp.copy(read.methState,tmp.length(),0);
                    *(read.methState+tmp.length())='\0';

                    read.mateposition = read.position;
                    read.position = (l_it->second).position;
                    read.mate_real_len = read.real_len;
                    read.real_len = (l_it->second).real_len;
                    neg_readDB.erase(l_it);
                }else
                    strcpy(read.matechrom, "*");
			    
			}
		}

    if(pos_readDB.size() > readDBSIZE){
        for (std::map<string, alignedread>::iterator iter = pos_readDB.begin(); iter != pos_readDB.end(); ) {
                strcpy((iter->second).matechrom, "*");
                pos_vDnaMethylMap.push_back(iter->second);
                pos_readDB.erase(iter++);
        }
        pos_smallerrightpos=0;
    }

    if(neg_readDB.size() > readDBSIZE){
        for (std::map<string, alignedread>::iterator iter = neg_readDB.begin(); iter != neg_readDB.end(); ) {
                strcpy((iter->second).matechrom, "*");
                Neg_vDnaMethylMap.push_back(iter->second);
                neg_readDB.erase(iter++);
        }
        neg_smallerrightpos = 0;
    }

	        if(now_chr!=old_chr) {
			if(printedchr != now_chr){
				sprintf(printinf, "procesing chromosome, %s \n", now_chr.c_str());
				Show_log(printinf);
				printedchr = now_chr;
			}
	                if(hetero_SNP) {
                		plus_HeterSNPs.clear();
        	        	neg_HeterSNPs.clear();
	                	if(chromSNPs[H].HeterSNPs.size()>0) heterSNP_chrom=true;
                		else heterSNP_chrom=false;
        	    	}
		    	prevloc = nowloc;
        	}else{
	        	if(prevloc != 0 && nowloc < prevloc) {
				fprintf(stderr, "\nalignment file is not sorted, please sort the sam or bam file first. %s %d %d\n", read.readid, read.position, prevloc);
				exit(0);
        	        }else prevloc = nowloc;
		}

		if(!(read.flag & 0x4)) //match
		{
			if(withFirstResult==0)
			{
		    		old_chr=read.chrom;
		    		withFirstResult=-1;
		   	}
		   	if( ( ((read.flag & 0x40) && (read.flag & 0x20)) || ((read.flag & 0x80) && (read.flag & 0x10 ))) ) //+ strand 96 144
		   	{
                if(chromMeths[H].PlusHeterMeths.size() + chromSNPs[H].HeterSNPs.size() <2) continue;
                if(now_chr==old_chr)// the same chromsome
                {
                    if( pos_first_start_coordinate==-1 || pos_first_end_coordinate==-1) 
                        if(hetero_SNP && heterSNP_chrom) plus_HeterSNPs=chromSNPs[H].HeterSNPs;
                        
                    int filter=filterRead_pair(plus_HeterSNPs, chromMeths[H].PlusHeterMeths, read, read_valid_len);
                    if(filter == -1) continue;

                    while(chromMeths[H].PlusHeterMeths.size() > 0 && read.position > chromMeths[H].PlusHeterMeths[0]){
                        if(pos_validCs.size()==0 || (pos_validCs.size()>0 && pos_validCs.back() < chromMeths[H].PlusHeterMeths[0]))
                            pos_validCs.push_back(chromMeths[H].PlusHeterMeths[0]);
                        chromMeths[H].PlusHeterMeths.pop_front();
                    }

                    //add new reads
                    //read.strand='+';
                    pos_vDnaMethylMap.push_back(read);
                    if( pos_first_start_coordinate==-1 || pos_first_end_coordinate==-1)//diff chromsome
                    {
                        pos_first_start_coordinate = read.position;
			if(strcmp(read.matechrom, "=") == 0)
                        	pos_first_end_coordinate   = read.mateposition + read_valid_len;
			else pos_first_end_coordinate   = read.position + read_valid_len;
                    }
                }
				//deal reads
                if(pos_readDB.size()==0 && (pos_first_end_coordinate < read.position || now_chr!=old_chr || pos_vDnaMethylMap.size() > memSize))
                    sortdequebyvec(pos_vDnaMethylMap);
                    //sort(pos_vDnaMethylMap.begin(), pos_vDnaMethylMap.end());
				while(pos_readDB.size()==0 && (pos_first_end_coordinate < read.position || now_chr!=old_chr || pos_vDnaMethylMap.size() > memSize) )
	            {
					if(pos_vDnaMethylMap.size() < 1) break;
					//if(hetero_SNP && heterSNP_chrom && plus_HeterSNPs.size() >0)
						processOneRead_heterSNPs_paired('+',Pvalues, plus_HeterSNPs,pos_vDnaMethylMap, pos_processed_coordinate, pos_validCs,pos_print,pos_OUTFILE);
					//else
						//processOneRead_paired(plus_jointposition,'+',Pvalues,pos_vDnaMethylMap, pos_processed_coordinate, pos_potentialCs,pos_validCs,pos_print,plus_paired_validCs);
					pos_processed_coordinate = pos_first_end_coordinate - 1;
					pos_vDnaMethylMap.pop_front();
					if(pos_vDnaMethylMap.size()<=0) break;
					pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
					if(strcmp(pos_vDnaMethylMap[0].matechrom, "=") == 0)
						pos_first_end_coordinate = pos_vDnaMethylMap[0].mateposition + pos_vDnaMethylMap[0].mate_real_len;
					else
						pos_first_end_coordinate = pos_vDnaMethylMap[0].position + pos_vDnaMethylMap[0].real_len;
					// skip the reads already processed
					while(pos_first_end_coordinate - 1 <= pos_processed_coordinate) {
						pos_vDnaMethylMap.pop_front();
						if(pos_vDnaMethylMap.size()<=0) break;
						pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
					if(strcmp(pos_vDnaMethylMap[0].matechrom, "=") == 0)
						pos_first_end_coordinate   = pos_vDnaMethylMap[0].mateposition + pos_vDnaMethylMap[0].mate_real_len;
					else
						pos_first_end_coordinate = pos_vDnaMethylMap[0].position + pos_vDnaMethylMap[0].real_len;
					
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
			}else if( ( ((read.flag & 0x40) && (read.flag & 0x10)) || ((read.flag & 0x80) && (read.flag & 0x20 ))) ) //- strand 80 160
			{
//printf("\nYYYsss %d", chromMeths[H].NegHeterMeths.size() + chromSNPs[H].HeterSNPs.size());

                if(chromMeths[H].NegHeterMeths.size() + chromSNPs[H].HeterSNPs.size() <2) continue;
                if(neg_first_start_coordinate==-1 || neg_first_end_coordinate==-1)
                {
                    if(hetero_SNP && heterSNP_chrom) neg_HeterSNPs=chromSNPs[H].HeterSNPs;
                }
                if(now_chr==old_chr)
                {
                    int filter=filterRead_pair(neg_HeterSNPs, chromMeths[H].NegHeterMeths, read, read_valid_len);
                    if(filter== -1) continue;
                 // meth valid C

                    while(chromMeths[H].NegHeterMeths.size() > 0 && read.position > chromMeths[H].NegHeterMeths[0]){
                        if(Neg_validCs.size()==0 || (Neg_validCs.back() < chromMeths[H].NegHeterMeths[0]))
                            Neg_validCs.push_back(chromMeths[H].NegHeterMeths[0]);
                        chromMeths[H].NegHeterMeths.pop_front();
                    }

                    //add new reads
                    if(neg_first_start_coordinate==-1 || neg_first_end_coordinate==-1) {
                        neg_first_start_coordinate = read.position;
			if(strcmp(read.matechrom, "=") == 0)
                        	neg_first_end_coordinate   = read.mateposition + read_valid_len;
			else neg_first_end_coordinate   = read.position + read_valid_len;
                    }
                    Neg_vDnaMethylMap.push_back(read);
                }
                if(neg_readDB.size()==0 && (neg_first_end_coordinate < read.position || now_chr!=old_chr || Neg_vDnaMethylMap.size() > memSize))
                    //sort(Neg_vDnaMethylMap.begin(), Neg_vDnaMethylMap.end());
                    sortdequebyvec(Neg_vDnaMethylMap);
				//deal reads
				while(neg_readDB.size()==0 && (neg_first_end_coordinate < read.position || now_chr!=old_chr || Neg_vDnaMethylMap.size() > memSize) )
	            {
					if(Neg_vDnaMethylMap.size() < 1) break;
	            	//if(hetero_SNP && heterSNP_chrom && neg_HeterSNPs.size()>0)
	            		processOneRead_heterSNPs_paired('-',Pvalues,neg_HeterSNPs,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_validCs,neg_print,neg_OUTFILE);
	            	//else
	            		//processOneRead_paired(neg_jointposition,'-',Pvalues,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_potentialCs,Neg_validCs,neg_print,Neg_paired_validCs);
	            	neg_processed_coordinate = neg_first_end_coordinate - 1;
	            	Neg_vDnaMethylMap.pop_front();
	            	if(Neg_vDnaMethylMap.size()<=0) break;
				    neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
				    if(strcmp(Neg_vDnaMethylMap[0].matechrom, "=") == 0)
    				        neg_first_end_coordinate   = Neg_vDnaMethylMap[0].mateposition + Neg_vDnaMethylMap[0].mate_real_len;
				    else
					neg_first_end_coordinate   = Neg_vDnaMethylMap[0].position + Neg_vDnaMethylMap[0].real_len;
				    // skip the reads already processed
			        while(neg_first_end_coordinate - 1 <= neg_processed_coordinate) {
			        	Neg_vDnaMethylMap.pop_front();
			        	if(Neg_vDnaMethylMap.size()<=0) break;
				        neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
					if(strcmp(Neg_vDnaMethylMap[0].matechrom, "=") == 0)
				            neg_first_end_coordinate   = Neg_vDnaMethylMap[0].mateposition + Neg_vDnaMethylMap[0].mate_real_len;
					else
					    neg_first_end_coordinate   = Neg_vDnaMethylMap[0].position + Neg_vDnaMethylMap[0].real_len;
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
				pos_validCs.clear();
				Neg_vDnaMethylMap.clear();
				Neg_validCs.clear();
				//add new potentialCs
	        	//add new reads and init 
				pos_first_start_coordinate = -1;
				pos_first_end_coordinate   = -1;
				pos_processed_coordinate   = -1;
				neg_first_start_coordinate = -1; 
				neg_first_end_coordinate   = -1;
				neg_processed_coordinate   = -1; 

				if( ( ((read.flag & 0x40) && !(read.flag & 0x10)) || ((read.flag & 0x80) && (read.flag & 0x10 ))) ) //+ strand
				{
					if(hetero_SNP && heterSNP_chrom)
					{
						plus_HeterSNPs=chromSNPs[H].HeterSNPs;
					}
                    int filter=0;
                    if((filter=filterRead_pair(plus_HeterSNPs, chromMeths[H].PlusHeterMeths, read, read_valid_len)) == -1)
                        continue;
                            // meth valid C

                    while(chromMeths[H].PlusHeterMeths.size() > 0 && read.position > chromMeths[H].PlusHeterMeths[0]){
                        if(pos_validCs.size()==0 || (pos_validCs.back() < chromMeths[H].PlusHeterMeths[0]))
                            pos_validCs.push_back(chromMeths[H].PlusHeterMeths[0]);
                        chromMeths[H].PlusHeterMeths.pop_front();
                    }

                    while(pos_validCs.size()>0 && pos_first_start_coordinate > pos_validCs[0])
                        pos_validCs.pop_front();

					pos_first_start_coordinate = read.position;
					if(read.flag & 0x2 && strcmp(read.matechrom, "=") == 0) pos_first_end_coordinate   = read.mateposition + read_valid_len;
					else pos_first_end_coordinate   = read.position + read_valid_len;
					//read.strand='+';
					pos_vDnaMethylMap.push_back(read);
                    for (map<string, alignedread>::iterator l_it=pos_readDB.begin(); l_it!=pos_readDB.end(); /*i++*/){
                        pos_vDnaMethylMap.push_back(l_it->second);
                        map<string, alignedread>::iterator iterTemp = l_it;
                        ++l_it;
                        pos_readDB.erase(iterTemp);
                    }
				}else if( ( ((read.flag & 0x40) && (read.flag & 0x10)) || ((read.flag & 0x80) && !(read.flag & 0x10 ))) ) //- strand
				{
					if(hetero_SNP && heterSNP_chrom)
					{
						neg_HeterSNPs=chromSNPs[H].HeterSNPs;
					}
                    int filter=0;
                    if((filter=filterRead_pair(neg_HeterSNPs, chromMeths[H].NegHeterMeths, read, read_valid_len)) == -1)
                        continue;
                            // meth valid C

                    while(chromMeths[H].NegHeterMeths.size() > 0 && read.position > chromMeths[H].NegHeterMeths[0]){
                        if(Neg_validCs.size()==0 || (Neg_validCs.back() < chromMeths[H].NegHeterMeths[0]))
                            Neg_validCs.push_back(chromMeths[H].NegHeterMeths[0]);
                        chromMeths[H].NegHeterMeths.pop_front();
                    }

                    while(Neg_validCs.size()>0 && neg_first_start_coordinate > Neg_validCs[0])
                        Neg_validCs.pop_front();
	            			
					neg_first_start_coordinate = read.position;
					if(read.flag & 0x2 && strcmp(read.matechrom, "=") == 0) neg_first_end_coordinate = read.mateposition + read_valid_len;
					else neg_first_end_coordinate = read.position + read_valid_len;
					//read.strand='-';
                    for (map<string, alignedread>::iterator l_it=neg_readDB.begin(); l_it!=neg_readDB.end(); /*i++*/){
                        Neg_vDnaMethylMap.push_back(l_it->second);
                        map<string, alignedread>::iterator iterTemp = l_it;
                        ++l_it;
                        neg_readDB.erase(iterTemp);
                    }
					Neg_vDnaMethylMap.push_back(read);
				 }
			}
		}// end if
		old_chr=now_chr;
	}//end sam
	Show_Progress_float(nReads);
	//sam done!
	//foreach the last  +strand
    
    for (map<string, alignedread>::iterator l_it=pos_readDB.begin(); l_it!=pos_readDB.end(); /*i++*/){
        pos_vDnaMethylMap.push_back(l_it->second);
        map<string, alignedread>::iterator iterTemp = l_it;
        ++l_it;
        pos_readDB.erase(iterTemp);
    }
    sortdequebyvec(pos_vDnaMethylMap);
    //sort(pos_vDnaMethylMap.begin(), pos_vDnaMethylMap.end());
    int lastpos;alignedread pa;
    if(pos_vDnaMethylMap.size()>0) pa = pos_vDnaMethylMap.back();
    if(strcmp(pa.matechrom, "=") != 0) lastpos = pa.position + pa.real_len;
    else lastpos = pa.mateposition + pa.mate_real_len;
    
    while(chromMeths[H].PlusHeterMeths.size() > 0 && lastpos >= chromMeths[H].PlusHeterMeths[0] ){
        if(pos_validCs.size()==0 || (pos_validCs.back() < chromMeths[H].PlusHeterMeths[0]))
            pos_validCs.push_back(chromMeths[H].PlusHeterMeths[0]);
        chromMeths[H].PlusHeterMeths.pop_front();
    }

	while(true) 
    {
		if(pos_vDnaMethylMap.size() < 1) break;

	    //if(hetero_SNP && heterSNP_chrom && plus_HeterSNPs.size() >0 )
	    	processOneRead_heterSNPs_paired('+',Pvalues,plus_HeterSNPs,pos_vDnaMethylMap, pos_processed_coordinate, pos_validCs,pos_print,pos_OUTFILE);
	    //else
	    	//processOneRead_paired(plus_jointposition,'+',Pvalues,pos_vDnaMethylMap, pos_processed_coordinate, pos_potentialCs,pos_validCs,pos_print,plus_paired_validCs);
 
	    pos_processed_coordinate = pos_first_end_coordinate - 1;
	    pos_vDnaMethylMap.pop_front();    
	    if(pos_vDnaMethylMap.size()<=0) break;
	    pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
	    if(strcmp(pos_vDnaMethylMap[0].matechrom, "=") == 0)
	        pos_first_end_coordinate   = pos_vDnaMethylMap[0].mateposition + pos_vDnaMethylMap[0].mate_real_len;
	    else
		pos_first_end_coordinate   = pos_vDnaMethylMap[0].position + pos_vDnaMethylMap[0].real_len;
			    // skip the reads already processed

	    while(pos_first_end_coordinate - 1 <= pos_processed_coordinate) {
	    	pos_vDnaMethylMap.pop_front();
	    	if(pos_vDnaMethylMap.size()<=0) break;
	    	pos_first_start_coordinate = pos_vDnaMethylMap[0].position;
		if(strcmp(pos_vDnaMethylMap[0].matechrom, "=") == 0)
	    	    pos_first_end_coordinate   = pos_vDnaMethylMap[0].mateposition + pos_vDnaMethylMap[0].mate_real_len;
		else
		    pos_first_end_coordinate   = pos_vDnaMethylMap[0].position + pos_vDnaMethylMap[0].real_len;
	    }
        if(pos_vDnaMethylMap.size() <= 1 && pos_print==1) 
        {
            fprintf(pos_OUTFILE,"******\n");
            pos_print=0;
            break;
        }
    }
    //-strand
    int nlastpos;alignedread na;
    for (map<string, alignedread>::iterator l_it=neg_readDB.begin(); l_it!=neg_readDB.end(); /*i++*/){
        Neg_vDnaMethylMap.push_back(l_it->second);
        map<string, alignedread>::iterator iterTemp = l_it;
        ++l_it;
        neg_readDB.erase(iterTemp);
    }
    sortdequebyvec(Neg_vDnaMethylMap);
//    sort(Neg_vDnaMethylMap.begin(), Neg_vDnaMethylMap.end());
    if(Neg_vDnaMethylMap.size()>0) na = Neg_vDnaMethylMap.back();
    if(strcmp(na.matechrom, "=") != 0) nlastpos = na.position + na.real_len;
    else nlastpos = na.mateposition + na.mate_real_len;
    while(chromMeths[H].NegHeterMeths.size() > 0 && nlastpos >= chromMeths[H].NegHeterMeths[0]){
        if(Neg_validCs.size()==0 || (Neg_validCs.back() < chromMeths[H].NegHeterMeths[0]))
            Neg_validCs.push_back(chromMeths[H].NegHeterMeths[0]);
        chromMeths[H].NegHeterMeths.pop_front();
    }

    while(true) 
    {
    	if(Neg_vDnaMethylMap.size() < 1) break;
    	//if(hetero_SNP && heterSNP_chrom && neg_HeterSNPs.size() >0)
    		processOneRead_heterSNPs_paired('-',Pvalues,neg_HeterSNPs,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_validCs,neg_print,neg_OUTFILE);
      	//else
      		//processOneRead_paired(neg_jointposition,'-',Pvalues,Neg_vDnaMethylMap, neg_processed_coordinate, Neg_potentialCs,Neg_validCs,neg_print,Neg_paired_validCs);
        neg_processed_coordinate = neg_first_end_coordinate - 1;
        Neg_vDnaMethylMap.pop_front();
        if(Neg_vDnaMethylMap.size()<=0) break;
		neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
		if(strcmp(Neg_vDnaMethylMap[0].matechrom, "=") == 0)
		    neg_first_end_coordinate   = Neg_vDnaMethylMap[0].mateposition + Neg_vDnaMethylMap[0].mate_real_len;
		else
		    neg_first_end_coordinate   = Neg_vDnaMethylMap[0].position + Neg_vDnaMethylMap[0].real_len;
			    // skip the reads already processed
		while(neg_first_end_coordinate - 1 <= neg_processed_coordinate) {
			Neg_vDnaMethylMap.pop_front();
			if(Neg_vDnaMethylMap.size()<=0) break;
		    neg_first_start_coordinate = Neg_vDnaMethylMap[0].position;
		    if(strcmp(Neg_vDnaMethylMap[0].matechrom, "=") == 0)
		        neg_first_end_coordinate   = Neg_vDnaMethylMap[0].mateposition + Neg_vDnaMethylMap[0].mate_real_len;
		    else
			neg_first_end_coordinate   = Neg_vDnaMethylMap[0].position + Neg_vDnaMethylMap[0].real_len;
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
    exit(0);
            		
}

void repalceReadmeth_old(char* methState){
    int len=strlen(methState);
    int i;
    char* q = NULL;
    q = methState; int qi=0;
    for(i=0; i<len; i++)
    {
        if(methState[i] == 'S' || methState[i] == 'I')
        {
            continue;
        }else
        {
            q[qi] = methState[i];
            qi++;
            continue;
        }
    }
    q[qi]='\0';
    methState = q;
}

void repalceReadmeth(char *methState)  
{
    if(methState == NULL) {
        fprintf(stderr, "\nAlign file dont contain MD tag, please run calmeth in BatMeth2 program first.\n");
        exit(0);
    }
    char *p=methState,*q=methState;  
    for(;*q!='\0';q++)  
    {  
        if(*q!='S' && *q!='I')  
            *p++=*q;   
    }  
    *p='\0';
    methState = p; 
    return;  
}  

int processbamread(const bam_header_t *header, const bam1_t *b, char* readid,int &Flag,char* Chrom,int &pos,int &mapQuality,char* CIG,char* Chrom_P,int &pos_P,int &Insert_Size,char* forReadString,char* forQuality, char* methState)
{
        uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
        int i;
        const bam1_core_t *c = &b->core;
        strcpy(readid,  bam1_qname(b));
        Flag = c->flag;

        if( (Flag & 0x100) || (Flag & 0x200) || (Flag & 0x400) || (Flag & 0x800) || (Flag & 0x4))
                return -1;
	if((Flag & 0x10) && (Flag & 0x20)) return -1;
	if(c->qual < Quality_cutoff || (int) c->pos + 1 <= 0 ) return -1;

        if (c->tid < 0) return -1;
        else {
                if (header) strcpy(Chrom, header->target_name[c->tid]); 
                else sprintf(Chrom, "%d", c->tid); 
        }
        pos = c->pos + 1;
        mapQuality = c->qual;

        //define
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
        
        s = bam1_aux(b);i=0;int md=0;
        while (s < b->data + b->data_len) {
                uint8_t type, key[2];
            key[0] = s[0]; key[1] = s[1];
            s += 2; type = *s; ++s;
            if(key[0] == 'M' && key[1] == 'D'){
                if(type == 'Z') 
                {
			md =1;
                    while (*s) {
			if(*s >='0' && *s<='9') {
				md=0;
				break;
			}
			methState[i++] = *s++; 
		    }
                    ++s;
		    break;
                }
            }else ++s;
        }
        methState[i]='\0';
        if(md == 0) {
		printf("\nError: no methstate MD:Z: tag in sam or bam file, please run calmeth with -f paramater to generate the new file contain MD:Z: tag.\n");
		exit(0);
	}
        return 1;        
}

//processed one reads of sam file while with heter SNPs
void processOneRead_heterSNPs_paired(char strand,vector<Pvals> &Pvalues,deque<heterSNP>& HeterSNPs,deque<alignedread> & vDnaMethylMap, int & processed_coordinate, deque<int> & validCs, int& isprint, FILE* OUTFILE) 
{

    if(validCs.size()>0){
        while(HeterSNPs.size() > 0 && validCs[0] > HeterSNPs[0].pos)
            HeterSNPs.pop_front();  
    }
    // coverage at potential Cs
    //validCs.clear(); // clear the valid Cs
    int i;
    for(i=0; i< HeterSNPs.size(); i++){
        if(HeterSNPs.size() == 0 || HeterSNPs[i].pos >= vDnaMethylMap.back().position ) // validCs.back())
            break;
        else if(find(validCs.begin(), validCs.end(), HeterSNPs[i].pos) == validCs.end())
            validCs.push_back(HeterSNPs[i].pos);
    }

    sort(validCs.begin(), validCs.end() );

    if(validCs.size() <= 1) {
        return;
    }
    
    char* chrom = vDnaMethylMap[0].chrom;

    int nStat = 4;int nStatSNP=37;
    int pairStat[nStat]; 
	map<string,int> pairStatSNP;
	//map<int ,int> pairStatSNP;
        // pairStat[0]: MM  ||  MR  RM  ||  MV1  V1M  ||  VV
        // pairStat[1]: MU   ||  MV  RU   ||  MV2  V1U   ||  VR
        // pairStat[2]: UM   ||  UR   VM  ||  UV1   V2M  ||  RV
        // pairStat[3]: UU    ||  UV   VU   ||  UV2   V2U   ||  RR
        int pairedState=0;  //  1   2 3
        int processed_snp_validC=0;
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
                 if(HeterSNPs.size() > processed_snp_validC+1 && HeterSNPs.size() > processed_snp_validC+1 && HeterSNPs[processed_snp_validC+1].pos == coordinate2)
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
            int pos2i = 2;
            for(int iMap=0; iMap < vDnaMethylMap.size(); iMap++) {
            	alignedread read = vDnaMethylMap[iMap];
                if(strcmp(read.matechrom, "=")==0){
                    if(coordinate1 > read.position + read.real_len && coordinate1 < read.mateposition)
                        continue;
                    while(read.mateposition > coordinate2 && read.position + read.real_len < coordinate2) {
                        if(i < validCs.size()-pos2i) {
                            coordinate2 = validCs[i+pos2i];
                            pos2i++;
                        }else break;
                    }

                }
                int index1;
                if(coordinate1 < read.position + read.real_len)
                    index1 = coordinate1 - read.position;
                else index1 = coordinate1 - read.mateposition + read.real_len + 1;
                if(index1 < 0) {
                    break;
                }
                int index2;
                if(coordinate2 < read.position + read.real_len)
                    index2 = coordinate2 - read.position;
                else index2 = coordinate2 - read.mateposition + read.real_len + 1;

                if(index2 < 0) {
                    break;
                }

                int RealLen1=0,RealLen2=0;//genome move
                int ReadStep1=0,ReadStep2=0;//read move
                if((index1 < read.real_len + read.mate_real_len  )&&(index2 < read.real_len + read.mate_real_len )) {

                	//ReadStep2--;
                    char methylationLabel_1 = Char2Trans[ read.methState[index1] ]; 
                    char methylationLabel_2 = Char2Trans[ read.methState[index2] ]; 
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
            }else// 1 2 3
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
	    }//end pairedState!=0
        //printf("\nTes\n");
      
    	if(pairedState==2)
    		processed_snp_validC++;
	    if(pairedState==3)
	    	processed_snp_validC++;
	     
	}//end of validCs
	
    while(validCs.size()>0) {
        if(validCs.size()>1) { //remove = //because 
            validCs.pop_front();
        }
        else {
            break;
        }
    }
    if(validCs.size()>0){
        while(HeterSNPs.size() > 0 && validCs[0] > HeterSNPs[0].pos)
            HeterSNPs.pop_front();
    }
    //delete[] coverage;
}

//processed one reads of sam file
void processOneRead_paired(deque<twopoint> jointposition,char strand,vector<Pvals> &Pvalues, deque<alignedread> & vDnaMethylMap, int & processed_coordinate,deque<int> & potentialCs,deque<int> & validCs,int& isprint,deque <validC> &paired_validCs) {
        FILE* OUTFILE;
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
        int end = vDnaMethylMap[0].position + vDnaMethylMap[0].real_len;

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
                if(index < read.real_len) {
                	int len=strlen(read.methState);
                	if(len>0 && read.methState[len-1] != '+')
                	{
                		while(read.methState[ReadStep] != '+')
                		{
                			ReadStep++;
                		}
                		ReadStep++;
                	}
	                for(;ReadStep<len;ReadStep++)
	                {
	                	if(read.methState[ReadStep] == '+')
	                		break;
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
              
            	  int len=strlen(read.methState);
            	  if(len>0 && read.methState[len-1] != '+')
            	  {
            		  while(read.methState[ReadStep1] != '+')
            		  {
            			  ReadStep1++;
            		  }
            		  ReadStep1++;
            	  }
	                for(;ReadStep1<len;ReadStep1++)
	                {
	                	if(read.methState[ReadStep1] == '+')
	                		break;
	                	else if(read.methState[ReadStep1] == 'S' || read.methState[ReadStep1] == 'I')
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
	                	if(read.methState[ReadStep1] == '+')
	                		break;
	                	else if(read.methState[ReadStep2] == 'S' || read.methState[ReadStep2] == 'I')
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
		     fprintf(OUTFILE,"%s\t%d\t%d\n",chrom,coordinate1,coordinate2);
		     	 for(int j=0; j<nStat; j++) {
		     		 fprintf(OUTFILE,"\t%d",pairStat[j]);
		     	 }
		     double pvalue = fishers_exact(pairStat[0],pairStat[1],pairStat[2],pairStat[3]); //MM+UU  ----  MU+UM
		     fprintf(OUTFILE,"\t%f\n",pvalue);
		     //Pvals pval;
		     //chrom_offset++;pval.coordinate=chrom_offset;pval.pvalue=pvalue;pval.strand=strand;
		     //Pvalues.push_back(pval);
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

//filter reads
int filterRead_pair(deque<heterSNP> & HeterSNPs, deque<int> & HeterMeths, alignedread & read, int read_valid_len)
{
    if(strcmp(read.matechrom, "=") != 0){
        if(HeterSNPs.size()==0){
            if((HeterMeths[0] <= read.position + read_valid_len && HeterMeths[1] > read.position + read_valid_len) || HeterMeths[0] > read.position + read_valid_len)
            return -1;
        }else if(HeterMeths.size()==0){
            if( (HeterSNPs[0].pos <= read.position + read_valid_len && HeterSNPs[1].pos > read.position + read_valid_len) ||  HeterSNPs[0].pos > read.position + read_valid_len)
            return -1;
        }
        else{
            int minpos=0, secpos=0;
            deque<int> tmp;
            if(HeterSNPs.size()>0) tmp.push_back(HeterMeths[0]);
            if(HeterSNPs.size()>1) tmp.push_back(HeterMeths[1]);
            if(HeterSNPs.size() > 0) tmp.push_back(HeterSNPs[0].pos);
            if(HeterSNPs.size() > 1) tmp.push_back(HeterSNPs[1].pos);
            sort(tmp.begin(), tmp.end());
            if(tmp.size()<2) {fprintf(stderr, "\nUnexpected bug!\n");}
            minpos = tmp[0];
            secpos = tmp[1];
            if((minpos <= read.position + read_valid_len && secpos > read.position + read_valid_len) || minpos > read.position + read_valid_len)
                return -1;
        }
        return 1;
    }else{ // paired-end
        if(HeterSNPs.size()==0){
        if(HeterSNPs.size()==0){
            if((HeterMeths[0] <= read.mateposition + read.mate_real_len && HeterMeths[1] > read.mateposition + read.mate_real_len) || HeterMeths[0] > read.mateposition + read.mate_real_len)
                return -1;
        }else if(HeterMeths.size()==0){
            if( (HeterSNPs[0].pos <= read.mateposition + read.mate_real_len && HeterSNPs[1].pos > read.mateposition + read.mate_real_len ) ||  HeterSNPs[0].pos > read.mateposition + read.mate_real_len)
                return -1;
        }
        else{
            int minpos=0, secpos=0;
            deque<int> tmp;
            if(HeterSNPs.size()>0) tmp.push_back(HeterMeths[0]);
            if(HeterSNPs.size()>1) tmp.push_back(HeterMeths[1]);
            if(HeterSNPs.size() > 0) tmp.push_back(HeterSNPs[0].pos);
            if(HeterSNPs.size() > 1) tmp.push_back(HeterSNPs[1].pos);
            sort(tmp.begin(), tmp.end());
            if(tmp.size()<2) {fprintf(stderr, "\nUnexpected bug!\n");}
            minpos = tmp[0];
            secpos = tmp[1];
            if((minpos <= read.mateposition + read.mate_real_len && secpos > read.mateposition + read.mate_real_len) || minpos > read.mateposition + read.mate_real_len)
                return -1;
        }
        return 1; 
        }
    }
}
void sortdequebyvec(deque<alignedread> & vDnaMethylMap){
    vector<alignedread> tempsort;
    for(int i = 0; i < vDnaMethylMap.size(); i++) {
        tempsort.push_back(vDnaMethylMap.at(i));
    }
    sort(tempsort.begin(), tempsort.end(), vDna_compare);
    vDnaMethylMap.clear();
    for(int i = 0; i < tempsort.size(); i++) {
        vDnaMethylMap.push_back(tempsort.at(i));
    }
    vector<alignedread>().swap(tempsort); 
}

int Get_read_valid_Length_and_is(char* cig, int & containIS)
{
    containIS=0;
        if(cig==NULL)
        {
                printf("\nCIGAR is null\n");
                exit(0);
        }
        char temp[8];unsigned n=0;int readLength=0;
        while(*cig!='\0')
        {
                if(*cig>='0' && *cig<='9')
                {
                                temp[n]=*cig;
                                cig++;n++;
                }else if(*cig=='M')
                {
                        temp[n]='\0';int length=atoi(temp);
                        readLength+=length;
                        cig++;n=0;
                }else if(*cig=='D')
                {
                        temp[n]='\0';int length=atoi(temp);
                        cig++;n=0;
                        readLength+= length;
                }else
                {
                    if(containIS == 0 && (*cig=='I' || *cig=='S'))
                        containIS = 1;
                    cig++;n=0;
                    continue;
                }
        }
        return readLength;
}

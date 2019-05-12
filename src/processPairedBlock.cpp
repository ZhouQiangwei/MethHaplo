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
#include "processPairedBlock.hpp"
extern map <string,int> pairStatMap;
extern map <int,string> intpairStatMap;

extern char Char2Trans[];

using std::vector;using std::string;using std::deque;

bool twopoint_Comp(twopoint a,twopoint b)
{
    return a.coordinate1 < b.coordinate1;
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

//------MV meth--variant
void processSwicthMV(char& methylationLabel_1,char& methylationLabel_2,map<string,int> & pairStatSNP)
{
	char stat_tmp[2];
	stat_tmp[0]=methylationLabel_1;stat_tmp[1]=methylationLabel_2;stat_tmp[2]='\0';
	string stat_string(stat_tmp);
	pairStatSNP[stat_string]++;
}

//MM --- 0 meth--meth
void processSwicthMM(char& methylationLabel_1,char& methylationLabel_2,int* pairStat)
{
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
bool paired_block_Comp(pairedDNAmethyMap a, pairedDNAmethyMap b)
{
	return a.left_start < b.left_start;
}

void Show_log(const char *log)
{
        struct tm *current_date;
        time_t seconds;
        time(&seconds);
        current_date = localtime(&seconds);
        printf("\n[%d-%d-%d %2d:%2d] %s",current_date->tm_mon+1,current_date->tm_mday, 1900+current_date->tm_year,current_date->tm_hour,current_date->tm_min, log);
}

void Show_Progress_SNP(long & nReads)
{
	struct tm *current_date;
	time_t seconds;
	time(&seconds);
	current_date = localtime(&seconds);
	printf("\r[%d-%d-%d %2d:%2d] %ld hetero SNPs processed ...",current_date->tm_mon+1,current_date->tm_mday, 1900+current_date->tm_year,current_date->tm_hour,current_date->tm_min,nReads);
	fflush(stdout);
}

void Show_Progress_Meth(long & nReads)
{
        struct tm *current_date;
        time_t seconds;
        time(&seconds);
        current_date = localtime(&seconds);
        printf("\r[%d-%d-%d %2d:%2d] %ld hetero methy site processed ...",current_date->tm_mon+1,current_date->tm_mday, 1900+current_date->tm_year,current_date->tm_hour,current_date->tm_min,nReads);
        fflush(stdout);
}

void Show_Progress(long & nReads)
{
	struct tm *current_date;
	time_t seconds;
	time(&seconds);
	current_date = localtime(&seconds);
	printf("\r[%d-%d-%d %2d:%2d] %d Million reads processed ...",current_date->tm_mon+1,current_date->tm_mday, 1900+current_date->tm_year,current_date->tm_hour,current_date->tm_min,\
									(int)(nReads / 1000000) );
	fflush(stdout);
}

void Show_Progress_float(long & nReads)
{
	struct tm *current_date;
	time_t seconds;
	time(&seconds);
	current_date = localtime(&seconds);
	printf("\r[%d-%d-%d %2d:%2d] %.1f Million reads processed ...",current_date->tm_mon+1,current_date->tm_mday, 1900+current_date->tm_year,current_date->tm_hour,current_date->tm_min,\
									(float)((float)nReads / 1000000) );
	fflush(stdout);
}
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

std::string int2str(int &int_temp)
{
    std::stringstream stream;
    stream << int_temp;
    return stream.str();   
}

int Get_read_valid_Length(char* cig)
{
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
			cig++;n=0;
			continue;
		}
	}
	return readLength;
}

void write_block(deque<twopoint> &jointposition,FILE* outFile)
{
//	if(jointposition.size()>0)
	//	fprintf(outFile,"******\n");
	int up_coordinate1=0,up_coordinate2=0;
	while(jointposition.size()>0)
	{
		twopoint tpTMP;
		tpTMP=jointposition[0];
		jointposition.pop_front();
		if(tpTMP.coordinate1 > up_coordinate1 && tpTMP.coordinate1 > up_coordinate2 )
			fprintf(outFile,"******\n");
		if(tpTMP.pairedState==0)
		{
			fprintf(outFile,"%s\t%d\t%d\t%s\t%f\n",tpTMP.chrom,tpTMP.coordinate1,tpTMP.coordinate2,\
					tpTMP.jointstate.c_str(),tpTMP.pvalue);
		}
		else if(tpTMP.pairedState==1 || tpTMP.pairedState==2)
		{	if(tpTMP.isdoubleVar)
				fprintf(outFile,"%s\t%d\t%d\t%s\t%c|%c%c\t%f\n",tpTMP.chrom,tpTMP.coordinate1,tpTMP.coordinate2,\
						tpTMP.jointstate.c_str(),tpTMP.refBase,tpTMP.var1,tpTMP.var2,tpTMP.pvalue);
			else
				fprintf(outFile,"%s\t%d\t%d\t%s\t%c|%c\t%f\n",tpTMP.chrom,tpTMP.coordinate1,tpTMP.coordinate2,\
						tpTMP.jointstate.c_str(),tpTMP.refBase,tpTMP.var1,tpTMP.pvalue);
		}
		else if(tpTMP.pairedState==3)
		{	if(tpTMP.isdoubleVar && tpTMP.isdoubleVar2)
				fprintf(outFile,"%s\t%d\t%d\t%s\t%c|%c%c\t%c|%c%c\t%f\n",tpTMP.chrom,tpTMP.coordinate1,tpTMP.coordinate2,\
						tpTMP.jointstate.c_str(),tpTMP.refBase,tpTMP.var1,tpTMP.var2,\
						tpTMP.refBase2,tpTMP.var3,tpTMP.var4,tpTMP.pvalue);
			else if(!tpTMP.isdoubleVar && tpTMP.isdoubleVar2)
				fprintf(outFile,"%s\t%d\t%d\t%s\t%c|%c\t%c|%c%c\t%f\n",tpTMP.chrom,tpTMP.coordinate1,tpTMP.coordinate2,\
						tpTMP.jointstate.c_str(),tpTMP.refBase,tpTMP.var1,\
						tpTMP.refBase2,tpTMP.var3,tpTMP.var4,tpTMP.pvalue);
			else if(tpTMP.isdoubleVar && !tpTMP.isdoubleVar2)
				fprintf(outFile,"%s\t%d\t%d\t%s\t%c|%c%c\t%c|%c\t%f\n",tpTMP.chrom,tpTMP.coordinate1,tpTMP.coordinate2,\
						tpTMP.jointstate.c_str(),tpTMP.refBase,tpTMP.var1,tpTMP.var2,\
						tpTMP.refBase2,tpTMP.var3,tpTMP.pvalue);
			else if(!tpTMP.isdoubleVar && !tpTMP.isdoubleVar2)
				fprintf(outFile,"%s\t%d\t%d\t%s\t%c|%c\t%c|%c\t%f\n",tpTMP.chrom,tpTMP.coordinate1,tpTMP.coordinate2,\
						tpTMP.jointstate.c_str(),tpTMP.refBase,tpTMP.var1,\
						tpTMP.refBase2,tpTMP.var3,tpTMP.pvalue);
		}
		up_coordinate1=tpTMP.coordinate1;
		up_coordinate2=tpTMP.coordinate2;
	}
}
bool Onlypairedjointposition_uniq(twopoint a,twopoint b)
{
	return (a.coordinate1==b.coordinate1 && a.coordinate2==b.coordinate2 );
}
void print_block(deque<pairedDNAmethyMap> & paired_DnaMethylBlock,deque<twopoint> Onlypairedjointposition,int processed_coordinate,FILE* outFile)
{
	if( paired_DnaMethylBlock.size()==0) return;
	if(Onlypairedjointposition.size()>1) std::sort(Onlypairedjointposition.begin(),Onlypairedjointposition.end(),twopoint_Comp);
	if(Onlypairedjointposition.size()>1) 
	{
		Onlypairedjointposition.erase( unique(Onlypairedjointposition.begin(), Onlypairedjointposition.end(), Onlypairedjointposition_uniq), Onlypairedjointposition.end() );
	}
	
	while(paired_DnaMethylBlock.size()>0)
	{
		while(Onlypairedjointposition.size() >0 && (Onlypairedjointposition[0].coordinate2 < paired_DnaMethylBlock[0].left_start || \
				Onlypairedjointposition[0].coordinate1 >= paired_DnaMethylBlock[0].left_start && Onlypairedjointposition[0].coordinate2 <= paired_DnaMethylBlock[0].left_end))
		{
			Onlypairedjointposition.pop_front();
		}
		pairedDNAmethyMap MapTmp;
		MapTmp=paired_DnaMethylBlock[0];
		
		pairedDNAmethyMap restore_TMP = MapTmp;
		
		deque<twopoint> jointposition=MapTmp.jointposition;
		deque<twopoint> restore_jointposition=jointposition;
		if( MapTmp.right_end > processed_coordinate || \
				(Onlypairedjointposition.size()>0 && Onlypairedjointposition[0].coordinate2 > processed_coordinate ) )
		{
			return;
		}
		else
		{
			if(Onlypairedjointposition.size()==0 || \
					(Onlypairedjointposition.size()>0 && Onlypairedjointposition[0].coordinate1 > MapTmp.left_end )) //l1 l2 c1
			{
				//write
				write_block(paired_DnaMethylBlock[0].jointposition,outFile);
				
				paired_DnaMethylBlock.pop_front();
			}else if(Onlypairedjointposition.size()>0 && Onlypairedjointposition[0].coordinate1 >= MapTmp.left_start && \
					Onlypairedjointposition[0].coordinate1 <= MapTmp.left_end &&\
					Onlypairedjointposition[0].coordinate2 > MapTmp.left_end) // l1 c1 l2 < c2
			{
				if(paired_DnaMethylBlock.size()>1)
				{
					if(Onlypairedjointposition[0].coordinate2 < paired_DnaMethylBlock[1].left_start )//l2 c2 l2.1
					{
						restore_jointposition.push_back(Onlypairedjointposition[0]);
						restore_TMP.jointposition=restore_jointposition;
						restore_TMP.left_end=Onlypairedjointposition[0].coordinate2;
						paired_DnaMethylBlock.pop_front();
						paired_DnaMethylBlock.insert(paired_DnaMethylBlock.begin(),restore_TMP);
						//sort(paired_DnaMethylBlock.begin(), paired_DnaMethylBlock.end(),paired_block_Comp);
						
						
						Onlypairedjointposition.pop_front();
						
					}else if(Onlypairedjointposition[0].coordinate2 >= paired_DnaMethylBlock[1].left_start && \
							Onlypairedjointposition[0].coordinate2 <= paired_DnaMethylBlock[1].left_end )//l2 l2.1 c2 l2.2
					{
						restore_jointposition.push_back(Onlypairedjointposition[0]);
						restore_jointposition.insert(restore_jointposition.end(),paired_DnaMethylBlock[1].jointposition.begin(),paired_DnaMethylBlock[1].jointposition.end());
						restore_TMP.left_end=paired_DnaMethylBlock[1].left_end;
						restore_TMP.right_start=paired_DnaMethylBlock[1].right_start;restore_TMP.right_end=paired_DnaMethylBlock[1].right_end;
						restore_TMP.jointposition=restore_jointposition;
						paired_DnaMethylBlock.pop_front();paired_DnaMethylBlock.pop_front();
						paired_DnaMethylBlock.insert(paired_DnaMethylBlock.begin(),restore_TMP);
						//sort(paired_DnaMethylBlock.begin(), paired_DnaMethylBlock.end(),paired_block_Comp);
						
						Onlypairedjointposition.pop_front();
					}
					
				}else
				{
					restore_jointposition.push_back(Onlypairedjointposition[0]);
					restore_TMP.jointposition=restore_jointposition;
					restore_TMP.left_end=Onlypairedjointposition[0].coordinate2;
					paired_DnaMethylBlock.pop_front();
					paired_DnaMethylBlock.insert(paired_DnaMethylBlock.begin(),restore_TMP);
					Onlypairedjointposition.pop_front();
				}
			}else if(Onlypairedjointposition[0].coordinate1 < MapTmp.left_start && Onlypairedjointposition[0].coordinate2 >= MapTmp.left_start && \
					Onlypairedjointposition[0].coordinate2 <= MapTmp.left_end)//c1 l1 c2 l2
			{
				restore_jointposition.insert(restore_jointposition.begin(),Onlypairedjointposition[0]);
				restore_TMP.left_start=Onlypairedjointposition[0].coordinate1;
				restore_TMP.jointposition=restore_jointposition;
				paired_DnaMethylBlock.pop_front();
				paired_DnaMethylBlock.insert(paired_DnaMethylBlock.begin(),restore_TMP);
				//sort(paired_DnaMethylBlock.begin(), paired_DnaMethylBlock.end(),paired_block_Comp);
				
				
				Onlypairedjointposition.pop_front();
			}else if(Onlypairedjointposition[0].coordinate1 < MapTmp.left_start && Onlypairedjointposition[0].coordinate2 > MapTmp.left_end \
					) //c1 l1 l2 c2
			{
				if(Onlypairedjointposition.size() > 1) // c1 l1 l2 l2.1 c2 l2.2 
				{
					if( (Onlypairedjointposition[1].coordinate1 > Onlypairedjointposition[0].coordinate2) || \
							(MapTmp.left_end < Onlypairedjointposition[1].coordinate1) ) // c1 c2 c3 c4 || c1 l1 l2 c3 c2
					{//write
						write_block(paired_DnaMethylBlock[0].jointposition,outFile);				
						paired_DnaMethylBlock.pop_front();
					}else //if(Onlypairedjointposition[1].coordinate1 ) // c1 c3 c2 c4 or c1 c3 c4 c2
					{//write
						write_block(paired_DnaMethylBlock[0].jointposition,outFile);				
						paired_DnaMethylBlock.pop_front();
					}
				}else
				{//write
					write_block(paired_DnaMethylBlock[0].jointposition,outFile);				
					paired_DnaMethylBlock.pop_front();
				}
				
			}
			
			
		}

	}//end while paired_DnaMethylBlock
}

void process_paired_block(deque<alignedread> & paired_vDnaMethylMap,deque<pairedDNAmethyMap>::iterator cur_iter,deque<pairedDNAmethyMap> &pos_paired_DnaMethylBlock\
		,deque<validC> & paired_validCs,deque<twopoint>& Onlypairedjointposition)
{
	pairedDNAmethyMap Tmp;
	int left_start=cur_iter->left_start;
	int left_end=cur_iter->left_end;//end edge
	int right_start=cur_iter->right_start;
	int right_end=cur_iter->right_end;
	int left_blcok_coord=0,right_blcok_coord=0;
	int left_blcok_iter=0,right_blcok_iter=0;
	//int left_positive_validC_coord;right_positive_validC_coord;
	if(paired_validCs.size()<=1) return;
	int pairedState=0;//0 MM 1  MV meth--ref/var //2  VM  ref/var--meth  //3  VV  var--var
	
	while(paired_validCs.size()>0 && paired_validCs[0].coordinate < left_start )
	{
		paired_validCs.pop_front();
	}
	if(paired_validCs.size()<=1) return;
	
	//get the two point 
	int coordinate1=0;int coordinate2=0;
	for(int i=0; i < paired_validCs.size()-1; i++) {
		coordinate1 = paired_validCs[i].coordinate;
        coordinate2 = paired_validCs[i+1].coordinate;
        if( coordinate1 <= left_end && coordinate2>left_end)
        {
        	left_blcok_coord = coordinate1;
        	left_blcok_iter=i;
        }
        if( coordinate1 < right_start && coordinate2 >= right_start)
        {
        	right_blcok_coord = coordinate2;
        	right_blcok_iter=i;
        	break;
        }
	}//end validC
	if(left_blcok_iter==0 || right_blcok_iter==0) 
	{
		paired_vDnaMethylMap.clear();
		return;
	}
	if(coordinate1==0 || coordinate2==0) 
	{
		paired_vDnaMethylMap.clear();
		return;
	}
	//get pairedState
	if(paired_validCs[left_blcok_iter].isSNP && paired_validCs[right_blcok_iter].isSNP)
		pairedState=3;
	else if(!paired_validCs[left_blcok_iter].isSNP && paired_validCs[right_blcok_iter].isSNP)
		pairedState=1;
	else if(paired_validCs[left_blcok_iter].isSNP && !paired_validCs[right_blcok_iter].isSNP)
		pairedState=2;
	
    int nStat = 4;int nStatSNP=37;
    int pairStat[nStat]; map<string,int> pairStatSNP;
    for(int j=0; j<nStat; j++) {
        pairStat[j] = 0; 
    }

    twopoint thispoint;

	if(paired_vDnaMethylMap.size()>0 && paired_vDnaMethylMap[paired_vDnaMethylMap.size()-1].position < left_start || paired_vDnaMethylMap[paired_vDnaMethylMap.size()-1].position < left_blcok_coord)
	{
		paired_vDnaMethylMap.clear();
	}
	deque<alignedread>::iterator iter=paired_vDnaMethylMap.begin(); 
    while(iter!=paired_vDnaMethylMap.end())
    {
    	if( iter->position < left_start || iter->position < left_blcok_coord)//.position is the right mate start
    		iter++;
    	else 
    		break;
    }
    paired_vDnaMethylMap.erase(paired_vDnaMethylMap.begin(),iter);

    if(paired_vDnaMethylMap.size()<3)
    	return;
    
    char* chrom = paired_vDnaMethylMap[0].chrom;
    //process paired read
	for(int iMap=0; iMap < paired_vDnaMethylMap.size(); iMap++)
	{
		alignedread read = paired_vDnaMethylMap[iMap];
        int index1 = left_blcok_coord - read.mateposition;
        if(index1 < 0) {
            break;
        }
        int index2 = right_blcok_coord - read.position;
        if(index2 < 0) {
            break;
        }
        
        int readStep1=0,readStep2=0;//read move
        int RealLen1=0,RealLen2=0;//genome move
        while(true)//left methstate
        {
        	if(readStep1 >= strlen(read.methState) ) break;
        	if(RealLen1==index1) break;
        	if(read.methState[readStep1]=='S' || read.methState[readStep1]=='I')
        	{
        		readStep1++;
        		continue;
        	}
        	else if(read.methState[readStep1]=='D')
        	{
        		RealLen1++;readStep1++;
        		continue;
        	}
        	else
        	{
        		RealLen1++;
        	}
        	
        	readStep1++;	
        }
        bool right_end_begain=false;
        while(true)//right
        {
        	if(readStep2 >= strlen(read.methState) ) break;
        	if(read.methState[readStep2]=='+') 
        	{
        		right_end_begain=true;
        		readStep2++;
        		continue;
        	}
        	if(!right_end_begain) 
        	{
        		readStep2++;
        		continue;
        	}
        	
        	if(RealLen2==index2) break;
        	if(read.methState[readStep2]=='S' || read.methState[readStep2]=='I')
        	{
        		readStep2++;
        		continue;
        	}
        	else if(read.methState[readStep2]=='D')
        	{
        		RealLen2++;readStep2++;
        		continue;
        	}
        	else
        	{
        		RealLen2++;
        	}
        	
        	readStep2++;	
        }
        
        
        char methylationLabel_1 = Char2Trans[ read.methState[readStep1] ]; 
        char methylationLabel_2 = Char2Trans[ read.methState[readStep2] ];
        
	     if(pairedState==0)  //MM meth-meth
		     processSwicthMM(methylationLabel_1,methylationLabel_2,pairStat);
         else //1  MV meth--ref/var //2  VM  ref/var--meth  //3  VV  var--var
         {
         	if(pairedState==1 && methylationLabel_2=='=') methylationLabel_2= paired_validCs[right_blcok_iter].RefBase;
         	else  if(pairedState==2 && methylationLabel_1=='=') methylationLabel_1 = paired_validCs[left_blcok_iter].RefBase;
         	else if(pairedState==3)
         	{
         		if(methylationLabel_1=='=')
         			methylationLabel_1 =paired_validCs[left_blcok_iter].RefBase;
         		if(methylationLabel_2=='=')
         			methylationLabel_2 = paired_validCs[right_blcok_iter].RefBase;
         	}
         	processSwicthMV(methylationLabel_1,methylationLabel_2,pairStatSNP);
         }
	     
         
         
	}//end vDNAmethy
	
	     //store
         if(pairedState==0)
         {
	            int count=0;
	            for(int c=0;c<nStat;c++) count+=pairStat[c];
	            if(count > 0)
	            {
	            	strcpy(thispoint.chrom,chrom);
	            	thispoint.coordinate1=coordinate1;
	            	thispoint.coordinate2=coordinate2;
	            	for(int j=0; j<nStat; j++) {
	            		if(j==0)
	            			thispoint.jointstate=int2str(pairStat[j]);
	            		else
	            			thispoint.jointstate= (thispoint.jointstate + "\t" +int2str(pairStat[j]));
	            	}
	            	double pvalue = fishers_exact(pairStat[0],pairStat[1],pairStat[2],pairStat[3]);
	            	thispoint.pvalue=pvalue;
	            	thispoint.pairedState=0;
	            }
         }else //pairedState!=0
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
	    	   
	    	   if(count_rmM==0 && count > 0)
	    	   {
		            strcpy(thispoint.chrom,chrom);
		            thispoint.coordinate1=coordinate1;
		            thispoint.coordinate2=coordinate2;
		            
		    		thispoint.jointstate=int2str(pairStatSNP["MM"]);
		    		thispoint.jointstate= (thispoint.jointstate + "\t" +int2str(pairStatSNP["MU"]));
		    		thispoint.jointstate= (thispoint.jointstate + "\t" +int2str(pairStatSNP["UM"]));
		    		thispoint.jointstate= (thispoint.jointstate + "\t" +int2str(pairStatSNP["UU"]));
		            double pvalue = fishers_exact(pairStatSNP["MM"],pairStatSNP["MU"],pairStatSNP["UM"],pairStatSNP["UU"]);
			     	thispoint.pvalue=pvalue;
	    	    }else if(count > 0)
	            {
		            strcpy(thispoint.chrom,chrom);
		            thispoint.coordinate1=coordinate1;
		            thispoint.coordinate2=coordinate2;
		            map<string,int>::iterator itp;
	    			for(itp=pairStatSNP.begin();itp!=pairStatSNP.end();++itp)
	    			{
	    				if(itp->second>0)
	    				{
	    					if(thispoint.jointstate.empty())
					    		 thispoint.jointstate=int2str(itp->second);
					    	 else
					    		 thispoint.jointstate= (thispoint.jointstate + "\t" +int2str(itp->second));
	    				}
	    			}
	    			
	    	    	if(pairedState!=3) 
	    	    	{
		     			char stat_tmp[2];string stat_string; double pvalue;
			     		if(pairedState==1)
			     		{
			     			if(thispoint.refBase=paired_validCs[right_blcok_iter].doubleVar)
			     			{
			     				char var1=paired_validCs[right_blcok_iter].variantBase1;
			     				char var2=paired_validCs[right_blcok_iter].variantBase2;
				     			int MV1,MV2,UV1,UV2;
								stat_tmp[0]='M';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			MV1=pairStatSNP[ stat_string];
				     			stat_tmp[0]='M';stat_tmp[1]=var2;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			MV2=pairStatSNP[stat_string];
								stat_tmp[0]='U';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			UV1=pairStatSNP[stat_string];
				     			stat_tmp[0]='U';stat_tmp[1]=var2;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			UV2=pairStatSNP[ stat_string ];
				     			pvalue = fishers_exact(MV1,MV2,UV1,UV2);
				     			thispoint.pairedState=1;
				     			thispoint.refBase=paired_validCs[right_blcok_iter].RefBase;thispoint.isdoubleVar=true;
				     			thispoint.var1=paired_validCs[right_blcok_iter].variantBase1;
				     			thispoint.var2=paired_validCs[right_blcok_iter].variantBase2;
			     			}else
			     			{
			     				char refBase1=paired_validCs[right_blcok_iter].RefBase;
								char var1=paired_validCs[right_blcok_iter].variantBase1;
			     				int MR,MV,UR,UV;
			     				stat_tmp[0]='M';stat_tmp[1]=refBase1;stat_tmp[2]='\0';stat_string=stat_tmp; 
								MR=pairStatSNP[ stat_string ];
								stat_tmp[0]='M';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			MV=pairStatSNP[ stat_string ];
				     			stat_tmp[0]='U';stat_tmp[1]=refBase1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			UR=pairStatSNP[ stat_string ];
								stat_tmp[0]='U';stat_tmp[1]=var1;stat_tmp[2]='\0';stat_string=stat_tmp; 
				     			UV=pairStatSNP[ stat_string ];
				     			pvalue = fishers_exact(MR,MV,UR,UV);
				     			thispoint.pairedState=1;
				     			thispoint.refBase=paired_validCs[right_blcok_iter].RefBase;thispoint.isdoubleVar=false;
				     			thispoint.var1=paired_validCs[right_blcok_iter].variantBase1;
			     			}
			     		}else//2
			     		{
			     			if(thispoint.refBase=paired_validCs[left_blcok_iter].doubleVar)
			     			{
			     				char var1=paired_validCs[left_blcok_iter].variantBase1;
			     				char var2=paired_validCs[left_blcok_iter].variantBase2;
			     				int V1M,V1U,V2M,V2U;
								stat_tmp[0]=var1;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V1M=pairStatSNP[ stat_string ];
								stat_tmp[0]=var2;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V2M=pairStatSNP[ stat_string ];
								stat_tmp[0]=var1;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V1U=pairStatSNP[ stat_string ];
								stat_tmp[0]=var2;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
								V2U=pairStatSNP[ stat_string ];
				     			pvalue = fishers_exact(V1M,V1U,V2M,V2U);
				     			thispoint.pairedState=2;
				     			thispoint.refBase=paired_validCs[left_blcok_iter].RefBase;thispoint.isdoubleVar=true;
				     			thispoint.var1=paired_validCs[left_blcok_iter].variantBase1;
				     			thispoint.var2=paired_validCs[left_blcok_iter].variantBase2;
			     			}else
			     			{
			     				char refBase1=paired_validCs[left_blcok_iter].RefBase;
								char var1=paired_validCs[left_blcok_iter].variantBase1;
			     				int RM,RU,VM,VU;
			     				stat_tmp[0]=refBase1;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
			     				RM=pairStatSNP[ stat_string ];
			     				stat_tmp[0]=refBase1;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
			     				RU=pairStatSNP[ stat_string ];
			     				stat_tmp[0]=var1;stat_tmp[1]='M';stat_tmp[2]='\0';stat_string=stat_tmp; 
			     				VM=pairStatSNP[ stat_string ];
			     				stat_tmp[0]=var1;stat_tmp[1]='U';stat_tmp[2]='\0';stat_string=stat_tmp; 
			     				VU=pairStatSNP[ stat_string ];
				     			pvalue = fishers_exact(RM,RU,VM,VU);
				     			thispoint.pairedState=2;
				     			thispoint.refBase=paired_validCs[left_blcok_iter].RefBase;thispoint.isdoubleVar=false;
				     			thispoint.var1=paired_validCs[left_blcok_iter].variantBase1;
			     			}
			     		}
			     			
			     		thispoint.pvalue=pvalue;

			     }
			     else // VV
			     {
			    	 thispoint.pairedState=3;
			      	 char stat_tmp[2];string stat_string; 
			         if(paired_validCs[left_blcok_iter].doubleVar && paired_validCs[right_blcok_iter].doubleVar)
			         {
			        	 char var1=paired_validCs[left_blcok_iter].variantBase1;
			        	 char var2=paired_validCs[left_blcok_iter].variantBase2;
			        	 char var3=paired_validCs[right_blcok_iter].variantBase1;
			        	 char var4=paired_validCs[right_blcok_iter].variantBase2;	 
			           	   //V1 V2 || V3 V4
			           	   int V1V3,V1V4,V2V3,V2V4;
			           	   //VV
			           	   stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			           	   V1V3=pairStatSNP[ stat_string ];
			           	   stat_tmp[0]=var1;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
			           	   V1V4=pairStatSNP[ stat_string ];
			           	   stat_tmp[0]=var2;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			           	   V2V3=pairStatSNP[ stat_string ];
			           	   stat_tmp[0]=var2;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
			           	   V2V4+=pairStatSNP[ stat_string ];
			           	   double pvalue = fishers_exact(V1V3,V1V4,V2V3,V2V4);
			           	   thispoint.refBase=paired_validCs[left_blcok_iter].RefBase;thispoint.isdoubleVar=true;
			           	   thispoint.var1=paired_validCs[left_blcok_iter].variantBase1;
			           	   thispoint.var2=paired_validCs[left_blcok_iter].variantBase2;
			           	   
			           	   thispoint.refBase2=paired_validCs[right_blcok_iter].RefBase;thispoint.isdoubleVar2=true;
			           	   thispoint.var3=paired_validCs[right_blcok_iter].variantBase1;
			           	   thispoint.var4=paired_validCs[right_blcok_iter].variantBase2;
			           	   thispoint.pvalue=pvalue;			     			
			    	  }else if(!paired_validCs[left_blcok_iter].doubleVar && !paired_validCs[right_blcok_iter].doubleVar)
			    	  {
				         char refBase1=paired_validCs[left_blcok_iter].RefBase;
				         char var1=paired_validCs[left_blcok_iter].variantBase1;
				         char refBase2=paired_validCs[right_blcok_iter].RefBase;
				         char var3=paired_validCs[right_blcok_iter].variantBase1;
				         
			    	  	  int RR,RV,VR,VV;
			    	  	  stat_tmp[0]=refBase1;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	  	  RR=pairStatSNP[ stat_string ];
			    	  	  stat_tmp[0]=refBase1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	  	  RV=pairStatSNP[ stat_string ];
			    	  	  stat_tmp[0]=var1;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	  	  VR=pairStatSNP[ stat_string ];
			    	  	  stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	  	  VV=pairStatSNP[ stat_string ];
			    	  	  double pvalue = fishers_exact(RR,RV,VR,VV);
			           	  thispoint.refBase=paired_validCs[left_blcok_iter].RefBase;thispoint.isdoubleVar=false;
				          thispoint.var1=paired_validCs[left_blcok_iter].variantBase1;
				           	   
				          thispoint.refBase2=paired_validCs[right_blcok_iter].RefBase;thispoint.isdoubleVar2=false;
				          thispoint.var3=paired_validCs[right_blcok_iter].variantBase1;
			    	  	  thispoint.pvalue=pvalue;
			    	  }else if(paired_validCs[left_blcok_iter].doubleVar && !paired_validCs[right_blcok_iter].doubleVar)
			    	  {
					      char var1=paired_validCs[left_blcok_iter].variantBase1;
					      char var2=paired_validCs[left_blcok_iter].variantBase2;
					      char refBase2=paired_validCs[right_blcok_iter].RefBase;
					      char var3=paired_validCs[right_blcok_iter].variantBase1;
			    	     //
					      int V1R,V1V3,V2R,V2V3;
			    	      stat_tmp[0]=var1;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	      V1R=pairStatSNP[ stat_string ];
			    	      stat_tmp[0]=var2;stat_tmp[1]=refBase2;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	      V2R=pairStatSNP[ stat_string ];
			    	      stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	      V1V3=pairStatSNP[ stat_string ];
			    	      stat_tmp[0]=var2;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    	      V2V3=pairStatSNP[ stat_string ];
			    	      double pvalue = fishers_exact(V1R,V1V3,V2R,V2V3);
			           	  thispoint.refBase=paired_validCs[left_blcok_iter].RefBase;thispoint.isdoubleVar=true;
				          thispoint.var1=paired_validCs[left_blcok_iter].variantBase1;
				          thispoint.var2=paired_validCs[left_blcok_iter].variantBase2;
				           	   
				          thispoint.refBase2=paired_validCs[right_blcok_iter].RefBase;thispoint.isdoubleVar2=false;
				          thispoint.var3=paired_validCs[right_blcok_iter].variantBase1;
			    	    
				          thispoint.pvalue=pvalue;
			    	 }else if(!paired_validCs[left_blcok_iter].doubleVar && paired_validCs[right_blcok_iter].doubleVar)
			    	 {
				         char refBase1=paired_validCs[left_blcok_iter].RefBase;
				         char var1=paired_validCs[left_blcok_iter].variantBase1;
				         char var3=paired_validCs[right_blcok_iter].variantBase1;
				         char var4=paired_validCs[right_blcok_iter].variantBase2;
			    		 //
			    		 int RV3,RV4,V1V3,V1V4;
			    		 stat_tmp[0]=refBase1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    		 RV3=pairStatSNP[ stat_string ];
			    		 stat_tmp[0]=refBase1;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    		 RV4=pairStatSNP[ stat_string ];
			    		 stat_tmp[0]=var1;stat_tmp[1]=var3;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    		 V1V3=pairStatSNP[ stat_string ];
			    		 stat_tmp[0]=var1;stat_tmp[1]=var4;stat_tmp[2]='\0';stat_string=stat_tmp; 
			    		 V1V4=pairStatSNP[ stat_string ];
			    		 double pvalue = fishers_exact(RV3,RV4,V1V3,V1V4);
			           	 thispoint.refBase=paired_validCs[left_blcok_iter].RefBase;thispoint.isdoubleVar=false;
				         thispoint.var1=paired_validCs[left_blcok_iter].variantBase1;
				           	   
				         thispoint.refBase2=paired_validCs[right_blcok_iter].RefBase;thispoint.isdoubleVar2=true;
				         thispoint.var3=paired_validCs[right_blcok_iter].variantBase1;
				         thispoint.var4=paired_validCs[right_blcok_iter].variantBase2;
			    		 thispoint.pvalue=pvalue;
			    	 }
			     }//VV

		    }//end if count>0
	    }//end pairedState!=0
     
        Onlypairedjointposition.push_back(thispoint);
        paired_validCs.pop_front();

	
	
}

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

#define MAX_Len_ALLOWED 100
#define METHSTAT 2000
#define MEDIA 50

struct alignedread 
{
	//int readlength;
	int real_len;
	int mate_real_len;
	char readid[MAX_Len_ALLOWED];
	int flag; 
	char chrom[MAX_Len_ALLOWED]; 
	int Quality;
	//int mateQuality;
	char cigar[MEDIA];
	//char matecigar[MEDIA];
	char matechrom[MAX_Len_ALLOWED];  
	char methState[METHSTAT];
	//char matemethState[METHSTAT];
	//int matech; char matestrand; 
	int position;int mateposition; 
	//int mquality; 
	int IS;//span distance
	char sequence[METHSTAT]; 
	char quality[METHSTAT];
	char strand;
	int mismatches;
	//int indels; // no of mismatches and no of insertions/deletions
	
	//int cflag; 
	int span;  
	//int tid; int mtid; // matetid
	//int minposition;//paired end reads mapping min position
	bool operator< (const alignedread& rt) const
	{
		return position < rt.position;
	}
        bool operator== (const alignedread& rt) const
        {
                return position == rt.position;
        }
};

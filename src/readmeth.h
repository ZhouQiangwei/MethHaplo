#ifndef _READVARIANT_H
#define _READVARIANT_H 
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

int count_variants(char* methfile, char* vcffile, char* sampleid, int* samplecol, int NMETH, int NCOVER, float MFloat);

int parse_variant(VARIANT* variant, char* buffer, int samplecol);

int read_variantfile(char* methfile,char* vcffile, VARIANT* varlist, HASHTABLE* ht, int* hetvariants, int samplecol, int NMETH, int NCOVER, float NFloat,int* variants);

#endif

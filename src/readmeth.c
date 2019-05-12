#include <limits.h>
#include <stdio.h>
#include "hashmap.h"
#include "assert.h"

extern char meth_strand;
typedef struct _CHROMMETHS
{
    unsigned long mf_p;//file pos 
    int start;
    char* chrom;
} CHROMMETHS;

int count_methfl(char* buffer, int NMETH, int NCOVER, float MFloat ) {
    int i = 0, j = 0, s = 0, e = 0;

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

    char* tempstring = NULL;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    int nmeth = atoi(tempstring);
    if( nmeth < NMETH){
        free(tempstring);
        return 0;
    }
    free(tempstring);

    tempstring = NULL;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    int ncover = atoi(tempstring);
    if( ncover < NCOVER) {
       free(tempstring);
       return 0;
    }
    free(tempstring);

    if(!( (float)nmeth/ncover >= MFloat && (float)nmeth/ncover <= (float)(1.0 - MFloat) )) {
        return 0;
    }

    return 1;
}

// count the # of variants in VCF file to allocate space for VCF variant array

int count_variants(char* methfile, char* vcffile, char* sampleid, int* samplecol, int NMETH, int NCOVER, float MFloat) {
    FILE* fp = fopen(vcffile, "r");
    if (fp == NULL) {
        fprintf(stderr, "could not open file %s\n\n", vcffile);
        return -1;
    }
    int variants = 0;
    char buffer[100000];
    int i = 0, j = 0, cols = 0;

    while (fgets(buffer, 100000, fp)) {
        if (buffer[0] != '#') variants++; // this should work for non-VCF files as well.
        else if (buffer[0] == '#' && buffer[1] == '#') continue;
        else if (buffer[0] == '#' && buffer[1] == 'C' && buffer[2] == 'H' && buffer[3] == 'R' && buffer[4] == 'O' && buffer[5] == 'M') {
            // find the column of the sample we want to phase 
            j = 0;
            while (buffer[i++] != '\n') {
                if ((buffer[i] == ' ' || buffer[i] == '\t') && j == 1) j = 0;
                else if (buffer[i] != ' ' && buffer[i] != '\t' && buffer[i] != '\n' && j == 0) {
                    j = 1;
                    cols++;
                }
            }
        }
    }
    fclose(fp);
    fprintf(stderr, "VCF file %s has %d variants\n", vcffile, variants);

    FILE* mfp = fopen(methfile, "r");
    if (mfp == NULL) {
        fprintf(stderr, "could not open file %s\n\n", methfile);
        return -1; 
    }   
    int Mvariants = 0;

    while (fgets(buffer, 100000, mfp)) {
        if (buffer[0] != '#') Mvariants += count_methfl(buffer, NMETH, NCOVER, MFloat); // this should work for non-VCF files as well.
    }   
    fclose(mfp);
    fprintf(stderr, "methratio file %s has %d meths\n", methfile, Mvariants);

    return (variants+Mvariants);
}


//get meth file chromosome
void get_methchrom(char* buffer,char* prevchrom,unsigned long mf_p, map_t chrom_map){
    int i = 0, j = 0, s = 0, e = 0;
    int position = 0;
    char* chrom;
    
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    chrom = (char*) malloc(e - s + 1); 
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
    position = atoi(tempstring);
    free(tempstring);

    if( strcmp(prevchrom,chrom) != 0 ) {
        CHROMMETHS* chromms;
	chromms = malloc(sizeof(CHROMMETHS));
        chromms->mf_p = mf_p;
	chromms->start = position;
	chromms->chrom = malloc(strlen(chrom)*sizeof(char));
	strcpy(chromms->chrom, chrom);
	int error = hashmap_put(chrom_map, chromms->chrom, chromms);
	assert(error==MAP_OK);
	//free(chromms);
    }
    strcpy(prevchrom,chrom);
    free(chrom);
    return;
}

//process meth line
int process_meth(char* buffer,VARIANT* var_t,int vf_position,int* position,char* prevchrom,int NMETH, int NCOVER, float MFloat){
    int i = 0, j = 0, s = 0, e = 0;

    var_t->depth = 0;
    var_t->A1 = 0;
    var_t->A2 = 0;
    var_t->H1 = 0;
    var_t->H2 = 0;

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    //var_t->chrom = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) var_t->chrom[j - s] = buffer[j];
    var_t->chrom[j - s] = '\0';
    if( strcmp(var_t->chrom,prevchrom) != 0 ) 
	return -2; //not same chrom
    char* tempstring;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    (*position) = atoi(tempstring);
    var_t->position = *position;
    free(tempstring);

    tempstring = NULL; //strand + -
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    if(meth_strand != '.' && (meth_strand == '+' || meth_strand == '-' ) ) {
        if(meth_strand != tempstring[0]) {
	    free(tempstring);
	    return -1;
	}
    }

    if(strcmp(tempstring, "+") == 0){
	//var_t->RA = (char*) malloc(5);
	var_t->RA[0] = 'M'; var_t->RA[1] = '_'; var_t->RA[2] = '+'; var_t->RA[3] = '\0';
	//var_t->AA = (char*) malloc(5);
        var_t->AA[0] = 'U'; var_t->AA[1] = '_'; var_t->AA[2] = '+'; var_t->AA[3] = '\0';

    }else if(strcmp(tempstring, "-") == 0){
    	//var_t->RA = (char*) malloc(5);
        var_t->RA[0] = 'M'; var_t->RA[1] = '_'; var_t->RA[2] = '-'; var_t->RA[3] = '\0';
	//var_t->AA = (char*) malloc(5);
        var_t->AA[0] = 'U'; var_t->AA[1] = '_'; var_t->AA[2] = '-'; var_t->AA[3] = '\0';
    }
    var_t->type = 0;
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
    if( nmeth < NMETH)
	return -1;
    free(tempstring);

    tempstring = NULL;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    int ncover = atoi(tempstring);
    if( ncover < NCOVER) 
        return -1;
    free(tempstring);


    if(!( (float)nmeth/ncover >= MFloat && (float)nmeth/ncover <= (float)(1.0 - MFloat) )) {
	return -1;
    }

    //var_t->allele1 = (char*) malloc(strlen(var_t->RA) + 1);
    strcpy(var_t->allele1, var_t->RA);
    //var_t->allele2 = (char*) malloc(strlen(var_t->AA) + 1);
    strcpy(var_t->allele2, var_t->AA);
    //var_t->genotype = (char*) malloc(4);
    
    if( ( (float)nmeth/ncover >= MFloat && (float)nmeth/ncover <= (float) (1.0 - MFloat) ) && ncover >= NCOVER && nmeth >= NMETH ) {
        var_t->genotype[0] = '0'; var_t->genotype[1] = '/'; var_t->genotype[2] = '1'; var_t->genotype[3] = '\0';
        var_t->heterozygous = '1';
    }else if( (float)nmeth/ncover <= MFloat ) {
        var_t->genotype[0] = '0'; var_t->genotype[1] = '/'; var_t->genotype[2] = '0'; var_t->genotype[3] = '\0';
	var_t->heterozygous = '0';
    }else if(  (float)nmeth/ncover >= 1.0 - MFloat ) {
        var_t->genotype[0] = '1'; var_t->genotype[1] = '/'; var_t->genotype[2] = '1'; var_t->genotype[3] = '\0';
	var_t->heterozygous = '0';
    }else {
        var_t->genotype[0] = '1'; var_t->genotype[1] = '/'; var_t->genotype[2] = '1'; var_t->genotype[3] = '\0';
	var_t->heterozygous = '0';
    }
    
    if( *position < vf_position)
	return 1;
    else if( *position > vf_position)
        return 2;
    else { //equal
        var_t->AA[3] = 'V'; var_t->AA[4] = '\0';
        return 0;
    }
}


///////////////////////////////////////////////////////////////////////////////////
void init_var(VARIANT* var) {
	var->chrom = malloc(1000);
	strcpy(var->chrom, "-----");
    var->RA = malloc(1000);
    var->AA = malloc(1000);
    var->genotype = malloc(1000);
    var->allele1 = malloc(1000);
    var->allele2 = malloc(1000);
}
void free_var(VARIANT* var) {
	free(var->chrom);
    free(var->RA);
    free(var->AA);
    free(var->genotype);
    free(var->allele1);
    free(var->allele2);
}
void copy_var(VARIANT* snpfrag, VARIANT* snp_t) {
    snpfrag->chrom = (char*)malloc(strlen(snp_t->chrom)+1);
    strcpy(snpfrag->chrom, snp_t->chrom);
    snpfrag->position = snp_t->position;
    snpfrag->RA = (char*)malloc(strlen(snp_t->RA)+1);
    strcpy(snpfrag->RA, snp_t->RA);
    snpfrag->AA = (char*)malloc(strlen(snp_t->AA)+1);
    strcpy(snpfrag->AA, snp_t->AA);
    snpfrag->allele1 = (char*)malloc(strlen(snp_t->allele1)+1);
    strcpy(snpfrag->allele1, snp_t->allele1);
    snpfrag->allele2 = (char*)malloc(strlen(snp_t->allele2)+1);
    strcpy(snpfrag->allele2, snp_t->allele2);
    snpfrag->genotype = (char*)malloc(strlen(snp_t->genotype)+1);
    strcpy(snpfrag->genotype, snp_t->genotype);
    snpfrag->heterozygous = snp_t->heterozygous;
    snpfrag->type = snp_t->type;
}
void copy_pvar(VARIANT* snpfrag, VARIANT* snp_t) {
    strcpy(snpfrag->chrom, snp_t->chrom);
    snpfrag->position = snp_t->position;
    strcpy(snpfrag->RA, snp_t->RA);
    strcpy(snpfrag->AA, snp_t->AA);
    strcpy(snpfrag->allele1, snp_t->allele1);
    strcpy(snpfrag->allele2, snp_t->allele2);
    strcpy(snpfrag->genotype, snp_t->genotype);
    snpfrag->heterozygous = snp_t->heterozygous;
    snpfrag->type = snp_t->type;
}
//
int parse_meth(FILE* mfp,char* prevchrom,int position,int* i,VARIANT* varlist,int* Start, VARIANT* var_b, int NMETH, int NCOVER, float NFloat, int* hetvariants){
    char buffer[100000];
    VARIANT var_t;
    init_var(&var_t);
    int in_t = -2;
    int ret = -1;
    while (fgets(buffer, 100000, mfp)) {
	in_t = process_meth(buffer, &var_t, position, Start, prevchrom, NMETH, NCOVER, NFloat);
	if(in_t == -1) continue;
	else if(in_t == 2) {
	    copy_pvar(var_b, &var_t);// current pos big than var
	    free_var(&var_t);
	    return ret; 
	}
	else if(in_t == 0 || in_t == 1) {
            if(var_t.chrom == NULL) printf("\n NULL snp pointer\n");
	    copy_var(&varlist[(*i)], &var_t);
	    (*i)++;
	    if(in_t == 0) ret = 1;//pos same, and ignore var
	    (*hetvariants)++;
	}else {
	    ret = -2;
	    free_var(&var_t);
	    return ret;//finished this chromosome
	}
    }

    ret = -2;
    free_var(&var_t);
    return ret;
}

int read_variantfile(char* methfile,char* vcffile, VARIANT* varlist, HASHTABLE* ht, int* hetvariants, int samplecol, int NMETH, int NCOVER, float NFloat, int* variants) {
    FILE* fp = fopen(vcffile, "r");
    FILE* mfp = fopen(methfile,"r");
    char buffer[100000];
    int i = 0;
    //	char allele1[256]; char allele2[256]; char genotype[256]; int quality; 
    char prevchrom[256];
    strcpy(prevchrom, "----");
    int chromosomes = 0; //int blocks=0;
    *hetvariants = 0;
    int het = 0;

    map_t chrom_map;//char* - struct
    chrom_map = hashmap_new();
    
    unsigned long mf_p = ftell(mfp);
    while (fgets(buffer, 100000, mfp)) {
	if (buffer[0] == '#') {
	    mf_p = ftell(mfp);
	    continue;
	}
	else {
	    get_methchrom(buffer,prevchrom,mf_p,chrom_map);
	    mf_p = ftell(mfp);
	}
    }
    rewind(mfp);
    strcpy(prevchrom, "----");
    int same_chr = 0, Start = 0;
    VARIANT var_V, var_M;
    init_var(&var_V);
    init_var(&var_M);
    int ret = -1;
    while (fgets(buffer, 100000, fp)) {
        if (buffer[0] == '#') continue;
        else {
            het = parse_variant(&var_V, buffer, samplecol);
 //           (*hetvariants) += het;
	    if (strcmp(var_V.chrom, prevchrom) != 0) {
                // insert chromname into hashtable 
                insert_keyvalue(ht, var_V.chrom, strlen(var_V.chrom), chromosomes);
                chromosomes++;
            }
	    
	    if (strcmp(var_V.chrom, prevchrom) != 0) { // diff chromsome
		//store old chromsome meth
		if(same_chr == 1)
		{
		    parse_meth(mfp,prevchrom,INT_MAX,&i,varlist,&Start,&var_M, NMETH, NCOVER, NFloat, hetvariants);
		}
		
		CHROMMETHS *chrommeth;
		int error = hashmap_get(chrom_map, var_V.chrom, (void**)(&chrommeth));
		if( error == MAP_OK && chrommeth != NULL ) {
		    same_chr = 1;
	            Start = chrommeth->start;
		    rewind(mfp);
		    fseek(mfp, chrommeth->mf_p, SEEK_SET);
		    if(Start <= var_V.position){
		        ret = parse_meth(mfp,var_V.chrom,var_V.position,&i,varlist,&Start,&var_M, NMETH, NCOVER, NFloat, hetvariants);
		    }
		}else {
		    same_chr = 0;
		    Start = 0;
		}
		strcpy(prevchrom, var_V.chrom);
	    }else { // same chrom
		if(same_chr == 1 && Start <= var_V.position) {
                    if( strcmp(var_M.chrom, "-----") !=0  && var_M.position != varlist[i-1].position) {
                        copy_var(&varlist[i],&var_M);
                        i++;
                        (*hetvariants)++;
                        if(var_M.position == var_V.position) {
                                ret = -1;
                                continue;
                        }
                    }
	            ret = parse_meth(mfp,var_V.chrom,var_V.position,&i,varlist,&Start,&var_M, NMETH, NCOVER, NFloat, hetvariants);
		}
	    }
	    
	    if(ret == -1 ) { //&& het != 0 ) {
	        if(var_V.chrom == NULL) printf("\n NULL snp pointer\n");
	        if (var_V.type != 0) var_V.position++;
		copy_var(&varlist[i], &var_V);
		i++;
		ret = -1;
	    }
	    else if(ret == 1) {
	        ret = -1;
	        continue;
	    }
	    else if(ret == -2) {
		same_chr = 0;
		Start = 0;
		if(var_V.chrom == NULL) printf("\n NULL snp pointer\n");
		if (var_V.type != 0) var_V.position++;
		copy_var(&varlist[i], &var_V);
		i++;
		int error = hashmap_remove(chrom_map, var_V.chrom);
		assert(error==MAP_OK);
		ret = -1;
		continue;
	    }
            //if (het ==0) continue; else (*hetvariants)++;
            //	fprintf(stdout,"%s %d %s %s %s %s\n",varlist[i].chrom,varlist[i].position,varlist[i].RA,varlist[i].AA,varlist[i].genotype,prevchrom); 
        }
    }
    fclose(fp); //chromosomes--;
    if(ret != -2 && Start > var_V.position) {
        if(var_M.chrom != NULL) {
            copy_var(&varlist[i], &var_M);
            i++;
	    (*hetvariants)++;
        }   
        parse_meth(mfp,prevchrom,INT_MAX,&i,varlist,&Start,&var_M, NMETH, NCOVER, NFloat, hetvariants);
        int error = hashmap_remove(chrom_map, var_V.chrom);
//        assert(error==MAP_OK);
    }
    // deal remain meth
    int size_of_hapmap = hashmap_length(chrom_map);
    if( size_of_hapmap > 0){
        CHROMMETHS** chrommeths;
	chrommeths = (CHROMMETHS**) malloc(sizeof (CHROMMETHS*) * size_of_hapmap);
	hashmap_values(chrom_map, (any_t*)chrommeths);
	if(size_of_hapmap == 1) {
	    CHROMMETHS* meth_t = (CHROMMETHS*)chrommeths[0];
	    insert_keyvalue(ht, meth_t->chrom, strlen(meth_t->chrom), chromosomes);
	    chromosomes++;
	    fseek(mfp, meth_t->mf_p, SEEK_SET);
	    parse_meth(mfp,meth_t->chrom,INT_MAX,&i,varlist,&meth_t->start,&var_M, NMETH, NCOVER, NFloat, hetvariants);
	}else {
	    rewind(mfp);
 	    int bi = 0, j = 0, s = 0, e = 0;
            int position = 0;
            char* chrom;strcpy(prevchrom,"-------");
	    while(fgets(buffer, 100000, mfp)) {
	        while (buffer[bi] == ' ' || buffer[bi] == '\t') bi++;
	        s = bi;
	        while (buffer[bi] != ' ' && buffer[bi] != '\t') bi++;
	        e = bi;
	        chrom = (char*) malloc(e - s + 1);
	        for (j = s; j < e; j++) chrom[j - s] = buffer[j];
	        chrom[j - s] = '\0';
	        if( strcmp(prevchrom,chrom) != 0 ) {
		    CHROMMETHS *meth_t;
                    int error = hashmap_get(chrom_map, chrom, (void**)(&meth_t));
                    if( error == MAP_OK && meth_t != NULL ) {
			insert_keyvalue(ht, meth_t->chrom, strlen(meth_t->chrom), chromosomes);
			chromosomes++;
			parse_meth(mfp,meth_t->chrom,INT_MAX,&i,varlist,&meth_t->start,&var_M, NMETH, NCOVER, NFloat, hetvariants);
			error = hashmap_remove(chrom_map, chrom);
	                assert(error==MAP_OK);
		    }
    	 }
	        strcpy(prevchrom,chrom);
		free(chrom);
	}
	free_var(&var_V);
	free_var(&var_M);
/*
	    int is;
	    for(is = 0; is< size_of_hapmap; is++){
            //if(chrommeths[i] != NULL) {
                CHROMMETHS* meth_t = (CHROMMETHS*)chrommeths[is];
                // insert chromname into hashtable 
                insert_keyvalue(ht, meth_t->chrom, strlen(meth_t->chrom), chromosomes);
                chromosomes++;
		fseek(mfp, meth_t->mf_p, SEEK_SET);
	        parse_meth(mfp,meth_t->chrom,INT_MAX,&i,varlist,&meth_t->start,&var_M, NMETH, NCOVER, NFloat, hetvariants);
            //}
	    }
*/
	}
    }
    hashmap_free(chrom_map);
    fclose(mfp);
    fprintf(stderr, "vcffile %s chromosomes %d hetvariants %d %d\n", vcffile, chromosomes, i, *hetvariants);
    (*variants) = i;
    return chromosomes;

}

// build a physical map that maps  intervals on chromosomes to the first variant that precedes the start of that interval

void build_intervalmap(CHROMVARS* chromvars, int chromosomes, VARIANT* varlist, int snps) {
    int i = 0, j = 0, k = 0, blocks = 0;
    chromvars[j].first = 0;
    j = 0;
    for (i = 0; i < snps - 1; i++) {
        if (strcmp(varlist[i].chrom, varlist[i + 1].chrom) != 0) {
            chromvars[j].last = i;
            chromvars[j].variants = chromvars[j].last - chromvars[j].first + 1;
            //fprintf(stderr,"new chrom %d %d %s %s\n",j,chromvars[j].variants,varlist[i].chrom,varlist[i+1].chrom);
            j++;
            chromvars[j].first = i + 1;
        }
    }
    chromvars[j].last = i;
    //	int** intervalmap; // map 1000bp of chromosome to first snp in that region indexed by snp_array 
    // first SNP to the right of the given base position including that position
    for (j = 0; j < chromosomes; j++) {
        blocks = (int) (varlist[chromvars[j].last].position / BSIZE) + 2;
        chromvars[j].blocks = blocks;
        //	fprintf(stderr,"chromosomes %d blocks %d \n",j,blocks);
        chromvars[j].intervalmap = (int*) malloc(sizeof (int)*blocks);
        for (i = 0; i < blocks; i++) chromvars[j].intervalmap[i] = -1;
        //fprintf(stderr,"blocks for chrom %d: %d \n",j,blocks);
        k = chromvars[j].first;
        for (i = 0; i < blocks; i++) {
            while (varlist[k].position <= BSIZE * i && k < chromvars[j].last) k++;
            if (k == chromvars[j].last) break;
            if (varlist[k].position > BSIZE * i && chromvars[j].intervalmap[i] == -1) chromvars[j].intervalmap[i] = k;
            //		if (chromvars[j].intervalmap[i] != -1) printf("FSNPtoright chrom %d block %d: %d %d \n",j,BSIZE*i,chromvars[j].intervalmap[i],varlist[chromvars[j].intervalmap[i]].position);
            //			else printf("FSNPtoright chrom %d block %d: %d \n",j,BSIZE*i,intervalmap[j][i]);
        }
    }
}

// this will only work for pure insertions and deletions, not for block substitutions

int calculate_rightshift(VARIANT* varlist, int ss, REFLIST* reflist) {
    int i = 0, j = 0;
    int a1 = 0, a2 = 0;
    int shift = 0;
    a1 = strlen(varlist[ss].allele1);
    a2 = strlen(varlist[ss].allele2);
    if (a1 > a2 && a2 == 1) {
        i = varlist[ss].position; // first base of deletion assuming position is +1 and not previous base 
        j = varlist[ss].position + a1 - a2;
        while (i - 1 < reflist->lengths[reflist->current] && j - 1 < reflist->lengths[reflist->current] && reflist->sequences[reflist->current][i - 1] == reflist->sequences[reflist->current][j - 1]) {
            i++;
            j++;
            shift++;
        }
        return shift;
    } else if (a1 == 1 && a2 > a1) {
        i = 1;
        j = varlist[ss].position;
        while (j - 1 < reflist->lengths[reflist->current] && varlist[ss].allele2[i] == reflist->sequences[reflist->current][j - 1] && i < a2) {
            i++;
            j++;
            shift++;
        }
        if (i == a2) // covered the full length of the inserted bases
        {
            i = varlist[ss].position;
            while (i - 1 < reflist->lengths[reflist->current] && j - 1 < reflist->lengths[reflist->current] && reflist->sequences[reflist->current][i - 1] == reflist->sequences[reflist->current][j - 1]) {
                i++;
                j++;
                shift++;
            }
        }
        return shift;
    } else return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int count_variants_oldformat(char* snpfile) {
    FILE* sf;
    char ch = '0';
    int snps = 0;
    sf = fopen(snpfile, "r");
    if (sf == NULL) {
        fprintf(stdout, "error opening file \n");
        exit(0);
    }
    while (1) {
        ch = fgetc(sf);
        if (ch == EOF) break;
        if (ch == '\n') snps++;
    }
    fclose(sf);
    fprintf(stderr, "read %d variants from file %s\n", snps, snpfile);
    return snps;
}

int read_variantfile_oldformat(char* snpfile, VARIANT* varlist, HASHTABLE* ht, int snps) {
    //fprintf(stderr,"old variant format is no longer supported, use --VCF variantfile.VCF option \n"); return -1;

    time_t ts;
    time(&ts);
    srand48((long int) ts);
    FILE* sf; //  char ch ='0'; 

    char allele1[256];
    char allele2[256];
    char genotype[256];
    int quality;
    char type[64];
    char prevchrom[256];
    char chrom[256];
    strcpy(prevchrom, "----");
    int i = 0;
    int chromosomes = 0;
    //char buffer[10000];

    sf = fopen(snpfile, "r");
    for (i = 0; i < snps; i++) {
        fscanf(sf, "%s %s %d %s %s %s %d\n", type, chrom, &varlist[i].position, allele1, allele2, genotype, &quality);
        varlist[i].allele1 = (char*) malloc(strlen(allele1) + 1);
        strcpy(varlist[i].allele1, allele1);
        varlist[i].allele2 = (char*) malloc(strlen(allele2) + 1);
        strcpy(varlist[i].allele2, allele2);
        varlist[i].chrom = (char*) malloc(strlen(chrom) + 1);
        strcpy(varlist[i].chrom, chrom);
        if (strcmp(type, "SNP") == 0 || strstr(type, "SNP") != NULL || strstr(type, "SNV") != NULL) {
            varlist[i].type = 0;
            if (genotype[0] == varlist[i].allele1[0] && genotype[2] == varlist[i].allele2[0]) varlist[i].heterozygous = '1';
            else if (genotype[0] == varlist[i].allele2[0] && genotype[2] == varlist[i].allele1[0]) varlist[i].heterozygous = '1';
            else varlist[i].heterozygous = '0';
            //printf("variant SNP %s %d \n",varlist[i].chrom,varlist[i].position);
        } else if (strcmp(type, "DNM") == 0) {
            varlist[i].type = 0;
        } else if (strcmp(type, "DEL") == 0) {
            varlist[i].type = -1 * strlen(allele1);
        } else if (strcmp(type, "INS") == 0) {
            varlist[i].type = strlen(allele2);
        }
        if (strcmp(varlist[i].chrom, prevchrom) != 0) {
            //fprintf(stdout,"%d %s %d %s %s %s %d %s\n",varlist[i].type,varlist[i].chrom,varlist[i].position,allele1,allele2,genotype,quality,prevchrom); 
            //			fprintf(stderr,"chromosomes %d %d\n",chromosomes,i);
            // insert chromname into hashtable 
            insert_keyvalue(ht, varlist[i].chrom, strlen(varlist[i].chrom), chromosomes);
            strcpy(prevchrom, varlist[i].chrom);
            chromosomes++;
        }
    }
    fclose(sf); //chromosomes--;

    fprintf(stderr, "read %d variants from file %s chromosomes %d\n", snps, snpfile, chromosomes);
    return chromosomes;

}



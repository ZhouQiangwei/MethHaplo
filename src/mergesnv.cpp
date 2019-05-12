/*

methyhap output result is 4 haplotype, we need merge to 2 haplotype

 */
#include <stdio.h>
#include <string.h>
//#include <iostream>
#include <cstdlib>
#include <ctype.h>

char Char_To_Code[256];

struct allele{
	int pos;
	char* ref;
	char* var;
	char* allele0;
	char* allele1;
	int hap1;
	int hap2;
	char* type;
};

int count_haps(char* hapfile);
bool process_hapfile(char* buffer, struct hapfrag* phapfrag, int nphase);
void readhapfile(char* phapfile, char* nhapfile, struct hapfrag* phapfrags, struct hapfrag* nhapfrags);
void merge_hap(struct hapfrag* plushaps, struct hapfrag* neghaps, int phaplen, int nhaplen);
void get_score(struct hapfrag &plushap, struct hapfrag &neghap, int start, int end, int* f13, int* f14);
void print_hap(struct hapfrag* plushaps, struct hapfrag* neghaps, int *ps, int *pe, int *ns, int *ne, int* currentprint);
void cal_score(struct hapfrag* plushaps, struct hapfrag* neghaps, int* ps, int* pe, int* ns, int* ne);
void process_hap(struct hapfrag* plushaps, struct hapfrag* neghaps, int* ps, int* pe, int* ns, int* ne);
int count_phased(char* buffer);

int main(int argc, char** argv) {
	char phapfile[1024];
	char nhapfile[1024];
	strcpy(phapfile, "None");
	strcpy(nhapfile, "None");
	int i = 0;
	for(i = 1; i < argc; i++){
		if(strcmp(argv[i], "--phap") == 0){
			strcpy(phapfile, argv[i+1]);
		}else if(strcmp(argv[i], "--nhap") == 0){
			strcpy(nhapfile, argv[i+1]);
		}
	}
	if(strcmp(phapfile, "None") == 0 || strcmp(nhapfile, "None") == 0) {
		fprintf(stderr, "\nmergehap --phap plushap --nhap neghap.\n");
		fprintf(stderr, "\nphapfile or nhapfile is null.\n");
		exit(0);
	}
	struct hapfrag* phapfrags;
	struct hapfrag* nhapfrags;
	int phaplen = 0, nhaplen = 0;
	phaplen = count_haps(phapfile);
	nhaplen = count_haps(nhapfile);
	phapfrags = (struct hapfrag*) malloc(sizeof(struct hapfrag)*phaplen);
	nhapfrags = (struct hapfrag*) malloc(sizeof(struct hapfrag)*nhaplen);
	if(phaplen > 0 && nhaplen > 0)
		readhapfile(phapfile, nhapfile, phapfrags, nhapfrags);
	else {
		fprintf(stderr, "\nplus haplotype or negative haplotype size is 0, no need to merge .\n");
		exit(0);
	}
	merge_hap(phapfrags, nhapfrags, phaplen, nhaplen);
	//clean struct
	for(i = 0; i < phaplen; i++){
		free(phapfrags[i].alist);
	}
	free(phapfrags);
	for(i = 0; i < nhaplen; i++){
		free(nhapfrags[i].alist);
	}
	free(nhapfrags);
	return 0;
}

int count_haps(char* hapfile) {
	FILE* fp = fopen(hapfile, "r");
	int count_blocks = 0;
	char buffer[10000];
	while(fgets(buffer, 10000, fp)){
		if(buffer[0] != '#')
			count_blocks++;
	}
	fclose(fp);
	return count_blocks;
}

//fprintf(fp, "BLOCK: offset: %d len: %d phased: %d ", clist[i].offset + 1, clist[i].length, clist[i].phased);
//fprintf(fp, "SPAN: %d fragments %d %s:%d-%d\n", span, clist[i].frags, snpfrag[clist[i].offset].chromosome, snpfrag[clist[i].offset].position , snpfrag[clist[i].lastvar].position);
//BLOCK: offset: 234 len: 17 phased: 4 SPAN: 536 fragments 63 chr1:142858382-142858918
//234     1       0       chr1    142858382       C       T       0/1:189:88:88:38:49:55.68%:1.1808E-19:20:20:16:22:27:22 0       0.000000        0.000000

bool process_hapfile(char* buffer, struct hapfrag* hapfrag, int nphase) {
	int i = 0, j = 0, s = 0, e = 0;

	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;

	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	if(buffer[s] == '-') return false;
	hapfrag->alist[nphase].hap1 = buffer[s] - '0';

	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	hapfrag->alist[nphase].hap2 = buffer[s] - '0';

	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	if(nphase == 0){
		hapfrag->chrom = (char*) malloc(e-s+1);
		for(j = s; j < e; j++) hapfrag->chrom[j-s] = buffer[j];
		hapfrag->chrom[j-s] = '\0';	
	}

	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	char* tempstring;
	tempstring = (char*) malloc(e-s+1);
	for(j = s; j < e; j++) tempstring[j-s] = buffer[j];
	tempstring[j-s] = '\0';
	hapfrag->alist[nphase].pos = atoi(tempstring);
	if(nphase == 0)
		hapfrag->firstvar = atoi(tempstring);
	else if(nphase == hapfrag->snps - 1)
		hapfrag->lastvar = atoi(tempstring);
	free(tempstring);

	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	hapfrag->alist[nphase].ref = (char*) malloc(e-s+1);
	for(j = s; j < e; j++) hapfrag->alist[nphase].ref[j-s] = buffer[j];
	hapfrag->alist[nphase].ref[j-s] = '\0';

	if(hapfrag->alist[nphase].hap1 == 0) {
		hapfrag->alist[nphase].allele0 = (char*) malloc(e-s+1);
		for(j = s; j < e; j++) hapfrag->alist[nphase].allele0[j-s] = buffer[j];
		hapfrag->alist[nphase].allele0[j-s] = '\0';
	}else {
		hapfrag->alist[nphase].allele1 = (char*) malloc(e-s+1);
		for(j = s; j < e; j++) hapfrag->alist[nphase].allele1[j-s] = buffer[j];
		hapfrag->alist[nphase].allele1[j-s] = '\0';
	}


	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	hapfrag->alist[nphase].var = (char*) malloc(e-s+1);
	for(j = s; j < e; j++) hapfrag->alist[nphase].var[j-s] = buffer[j];
	hapfrag->alist[nphase].var[j-s] = '\0';

	if(hapfrag->alist[nphase].hap2 == 1) {
		hapfrag->alist[nphase].allele1 = (char*) malloc(e-s+1);
		for(j = s; j < e; j++) hapfrag->alist[nphase].allele1[j-s] = buffer[j];
		hapfrag->alist[nphase].allele1[j-s] = '\0';
	}else {
		hapfrag->alist[nphase].allele0 = (char*) malloc(e-s+1);
		for(j = s; j < e; j++) hapfrag->alist[nphase].allele0[j-s] = buffer[j];
		hapfrag->alist[nphase].allele0[j-s] = '\0';
	}
        hapfrag->alist[nphase].type = (char*) malloc(10);
        strcpy(hapfrag->alist[nphase].type, "SNV");
	if(hapfrag->alist[nphase].hap1 == 0) {
	    if((strcmp(hapfrag->alist[nphase].allele0,"C") == 0 && strcmp(hapfrag->alist[nphase].allele1, "T") == 0) || (strcmp(hapfrag->alist[nphase].allele0,"G") == 0 && strcmp(hapfrag->alist[nphase].allele1, "A") == 0) )
	        strcpy(hapfrag->alist[nphase].type, "METH");
	}else if(hapfrag->alist[nphase].hap1 == 1 ) {
		if( (strcmp(hapfrag->alist[nphase].allele0,"T") == 0 && strcmp(hapfrag->alist[nphase].allele1, "C") == 0) || (strcmp(hapfrag->alist[nphase].allele0,"A") == 0 && strcmp(hapfrag->alist[nphase].allele1, "G") == 0))
                strcpy(hapfrag->alist[nphase].type, "METH");
	}

//	printf("\nEEE %d %d %s %s\n", hapfrag->alist[nphase].hap1, hapfrag->alist[nphase].hap2,
//		hapfrag->alist[nphase].allele0, hapfrag->alist[nphase].allele1);
	return true;

}

//BLOCK: offset: 234 len: 17 phased: 4 SPAN: 536 fragments 63 chr1:142858382-142858918
int count_phased(char* buffer) {
	// the 7 col
	int i = 0, j = 0, s = 0, e = 0;
	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;

	while(buffer[i] == ' ' || buffer[i] == '\t') i++;
	s = i;
	while(buffer[i] != ' ' && buffer[i] != '\t') i++;
	e = i;
	char* tempstring;
	tempstring = (char*) malloc(e-s+1);
	for(j = s; j < e; j++) tempstring[j-s] = buffer[j];
	tempstring[j-s] = '\0';
	int countphases = atoi(tempstring);
	free(tempstring);

	return countphases;
}

void readhapfile(char* phapfile, char* nhapfile, struct hapfrag* phapfrags, struct hapfrag* nhapfrags){
	char buffer[10000];
	FILE* fp_ph = fopen(phapfile, "r");
	int snps = 0; int i = -1; int nphase = 0;
	while(fgets(buffer, 10000, fp_ph)){
		if(buffer[0] == '#')
			continue;
		if(buffer[0] == 'B' && buffer[1] == 'L' && buffer[2] == 'O'){
			snps = count_phased(buffer);
			phapfrags[i].snps = nphase;
			i++;
			phapfrags[i].snps = snps;
			phapfrags[i].printstate = 0;
			phapfrags[i].alist = (struct allele*) malloc( sizeof(struct allele) *snps);
			nphase = 0;
			continue;
		}
		if(process_hapfile(buffer, &phapfrags[i], nphase))
			nphase++;
	}
	phapfrags[i].snps = nphase;
	fclose(fp_ph);

	FILE* fp_nh = fopen(nhapfile, "r");
	i = -1; nphase = 0; snps = 0;
	while(fgets(buffer, 10000, fp_nh)){
		if(buffer[0] == '*')
			continue;
		if(buffer[0] == 'B' && buffer[1] == 'L' && buffer[2] == 'O'){
			snps = count_phased(buffer);
			nhapfrags[i].snps = nphase;
			i++;
			nhapfrags[i].snps = snps;
			nhapfrags[i].printstate = 0;
			nhapfrags[i].alist = (struct allele*) malloc(sizeof(struct allele)*snps);
			nphase = 0;
			continue;
		}
		
		if(process_hapfile(buffer, &nhapfrags[i], nphase))
			nphase++;
	}
	nhapfrags[i].snps = nphase;
	fclose(fp_nh);
}

void merge_hap(struct hapfrag* plushaps, struct hapfrag* neghaps, int phaplen, int nhaplen) {
    int ps = 0, pe = -1;
    int ns = 0, ne = -1;
    int i = 0, j = 0;
    int currentprint = 0;
    char pchrom[1024];
    char nchrom[1024];
    strcpy(pchrom, plushaps[i].chrom); strcpy(nchrom, neghaps[j].chrom);
    while (i < phaplen && j < nhaplen ) {

    	if (plushaps[i].lastvar < neghaps[j].firstvar) {
        	ne = j-1; pe = i;
     		process_hap(plushaps, neghaps, &ps, &pe, &ns, &ne);
     		//print out until plus i && neg j-1;
     		printf("#Block\n");
     		print_hap(plushaps, neghaps, &ps, &pe, &ns, &ne, &currentprint);
        	ps = i+1; ns = j;
        	i++;
	}else if (neghaps[j].lastvar < plushaps[i].firstvar){
        	pe = i-1; ne = j;
        	process_hap(plushaps, neghaps, &ps, &pe, &ns, &ne);
          	//print out until plus i && neg j;
          	printf("#Block\n");
        	print_hap(plushaps, neghaps, &ps, &pe, &ns, &ne, &currentprint);
        	ps = i; ns = j+1;
        	j++;
	    }else if(plushaps[i].lastvar == neghaps[j].lastvar) {
	    	pe = i; ne = j;
	    	process_hap(plushaps, neghaps, &ps, &pe, &ns, &ne);
	    	printf("#Block\n");
	    	print_hap(plushaps, neghaps, &ps, &pe, &ns, &ne, &currentprint);
	    	ps = i; ns = j;
	    	i++; j++;
	    }else if(plushaps[i].lastvar < neghaps[j].lastvar){
        	if(i == phaplen-1) {
        		pe = i; ne = j;
//        		int ops = ps, ons = ns, ope = pe, one = ne;
        		process_hap(plushaps, neghaps, &ps, &pe, &ns, &ne);
			int ops = ps, ons = ns, ope = pe, one = ne;
        		printf("#Blockt %d %d %d %d\n", ops, ope, ons, one);
     			print_hap(plushaps, neghaps, &ops, &ope, &ons, &one, &currentprint);
     			j++;
        	}
        	i++;
	    }else {
        	if(j == nhaplen-1) {
        		pe = i; ne = j;
        		process_hap(plushaps, neghaps, &ps, &pe, &ns, &ne);
        		printf("#Blockf\n");
     			print_hap(plushaps, neghaps, &ps, &pe, &ns, &ne, &currentprint);
     			i++;
        	}
	    	j++;
	    }
    }
    if (j < nhaplen ){
        ns = j; ne = nhaplen -1;
        ps = pe+1; 
        printf("#Block\n");
        print_hap(plushaps, neghaps, &ps, &pe, &ns, &ne, &currentprint);
    }else if (i < phaplen ){
        ps = i; pe = phaplen -1;
        ns = ne+1;
        printf("#Block\n");
        print_hap(plushaps, neghaps, &ps, &pe, &ns, &ne, &currentprint);
    }

}


//deal with the first overlap region
void process_hap(struct hapfrag* plushaps, struct hapfrag* neghaps, int* ps, int* pe, int* ns, int* ne){
	if ( (*ps < 0 || *pe - *ps <0) && ( *ns < 0 || *ne - *ns < 0) ){ //{plushaps overlap region size && neghaps overlap region size all is 0}
		return;
	}else if (*ps < 0 || *pe - *ps <0){ //{plushaps overlap region size is 0}
		//printf("#BLOCK\n");
		//print_hap(plushaps, neghaps, ps, pe, ns, ne);
	    //print neghaps[ns to ne]
	    return;
	}else if (*ns < 0 || *ne - *ns<0){   //{neghaps overlap region size is 0}
	    //print plushaps[ps to pe]
	    //printf("#BLOCK\n");
	    //print_hap(plushaps, neghaps, ps, pe, ns, ne);
	    return;
	}
	plushaps[*ps].printstate = 1;
	neghaps[*ns].printstate = 1;
        cal_score(plushaps, neghaps, ps, pe, ns, ne);
}



void cal_score(struct hapfrag* plushaps, struct hapfrag* neghaps, int* ps, int* pe, int* ns, int* ne) {
//	printf("\nGG %d %d %d %d\n", *ps, *pe, *ns, *ne);
	if (*ns > *ne && *ps > *pe)
    	return;
	int f13 = 0, f14 = 0;
	int start = 0, end = 0; // overlap region
	int next = 0; // means next hap is plushap (-1) or neghap (1)
	hapfrag plushap = plushaps[*ps];
	hapfrag neghap = neghaps[*ns];
	if (plushap.firstvar >= neghap.firstvar && plushap.lastvar >= neghap.lastvar){ //{case 2}
    	start = plushap.firstvar;
    	end = neghap.lastvar;
    	next = 1;
	}else if (plushap.firstvar >= neghap.firstvar && plushap.lastvar <= neghap.lastvar){ //{case 3}
    	start = plushap.firstvar;
    	end = plushap.lastvar;
    	next = -1;
	}else if (plushap.firstvar <= neghap.firstvar && plushap.lastvar >= neghap.lastvar){ //{case 4}
	    start = neghap.firstvar;
	    end = neghap.lastvar;
	    next = 1;
	}else if (plushap.firstvar <= neghap.firstvar && plushap.lastvar <= neghap.lastvar){ //{case1}
	    start = neghap.firstvar;
	    end = plushap.lastvar;
		next = -1;
	}
    // nexthap == 1 means 13 && 24, 
    // nexthap == 2 means 24 && 13
    get_score(plushaps[*ps], neghaps[*ns], start-1, end+1, &f13, &f14);
//    printf("\nScore %d %d\n", f13, f14);
    if (f13 >= 2 && f14 <= -2){
        if (next == -1){
        	plushaps[*ps].nexthap = 1;
        	if(plushaps[*ps].printstate == 0) plushaps[*ps].printstate =1;
         	neghaps[*ns].printstate = plushaps[*ps].printstate;
    	}else if (next == 1){
        	neghaps[*ns].nexthap = 1;
        	if(neghaps[*ns].printstate == 0) neghaps[*ns].printstate =1;
        	plushaps[*ps].printstate = neghaps[*ns].printstate;
        }
    }else if (f13 <= -2 && f14 >= 2){
        if (next == -1){
			plushaps[*ps].nexthap = 2;
			if(plushaps[*ps].printstate == 0) plushaps[*ps].printstate = 1;
			neghaps[*ns].printstate =  0 - plushaps[*ps].printstate;
    	}else if (next == 1){
        	neghaps[*ns].nexthap = 2;
        	if(neghaps[*ns].printstate == 0) neghaps[*ns].printstate =1;
        	plushaps[*ps].printstate = 0- neghaps[*ns].printstate;
        }
    }else {
    	plushaps[*ps].printstate = 0;
    	neghaps[*ns].printstate = 0;
    }
//    printf("\nPRINT %d %d\n", plushaps[*ps].printstate, neghaps[*ns].printstate);

    if (next == -1){
    	(*ps)++;
 //   	printf("\nTTT %d %d\n", *ps, *pe);
    	if (*ps <= *pe)
			cal_score(plushaps, neghaps, ps, pe, ns, ne);
        else 
        	(*ps)--;
        return;
    }else if (next == 1){
        (*ns)++;
//        printf("\nCCC %d %d\n", *ns, *ne);
    	if (*ns <= *ne)
        	cal_score(plushaps, neghaps, ps, pe, ns, ne);
    	else 
    		(*ns)--;
        return;
    }
}

struct Score_Matrix
{
	unsigned* Score;
};

int SCORE_M(char c1, char c2, char c3, char c4){
	int s = 2;
	if(c1 == 'C' && c2 == 'G' && c3 == 'T' && c4 == 'A')
		return 2*s;
	if(c1 == 'T' && c2 == 'A' && c3 == 'C' && c4 == 'G')
		return 2*s;
	return -2*s;
}

int SCORE(char c1, char c2, char c3, char c4){
	int s = 2;
	if(c1 == c2 && c3 == c4)
		return 3*s;
	else if(c1=='T' && c2=='C') {
		if(c3 == c4 && c3 == 'T') 
			return 0;
		else if(c3 == c4)
			return 3*s;
		else if(c3 == 'G' && c4 == 'A')
			return 2*s;
	}else if(c1 == 'C' && c2 == 'C') {
		if(c3 == 'G' && c4 == 'A')
			return 3*s;
	}else if(c1 == 'T' && c2 == 'T') {
		if(c3 == 'T' && c4 =='C')
			return 0;
		if(c1 == 'G' && c2 == 'A')
			return 3*s;
	}else if(c1 == 'G' && c2 == 'G') {
		if(c3 == 'T' && c4 == 'C')
			return 3*s;
	}else if(c1=='G' && c2=='A') {
		if(c3 == c4)
			if(c3 == 'A') return 0;
			else return 3*s;
		if(c3 == 'T' && c4 == 'C')
			return 2*s;
	}else if(c1 == 'A' && c2 == 'A') {
		if(c3 == 'T' && c4 == 'C')
			return 3*s;
		if(c3 == 'G' && c4 == 'A')
			return 3*s;
	}

	return -3*s;
}

void get_score(struct hapfrag &plushap, struct hapfrag &neghap, int start, int end, int* f13, int* f14) {
	int pi = 0, ni = 0; 
	char p1 = ' '; char p2 = ' '; char n3 = ' '; char n4 = ' ';
	while (plushap.alist[pi].pos < start)
    	pi++;
	while (neghap.alist[ni].pos < start)
    	ni++;

	while (plushap.alist[pi].pos <= end && neghap.alist[ni].pos <= end){
//printf("\nPOS %d %d\n", plushap.alist[pi].pos, neghap.alist[ni].pos);
	    if(pi >= plushap.snps || ni >= neghap.snps) break;
    	    if (plushap.alist[pi].pos + 1  < neghap.alist[ni].pos){
        	pi++;
        	continue;
    	    }else if (neghap.alist[ni].pos < plushap.alist[pi].pos){
        	ni++;
        	continue;
    	    }
    	    p1 = plushap.alist[pi].allele0[0];
    	    p2 = plushap.alist[pi].allele1[0];
    	    n3 = neghap.alist[ni].allele0[0];
    	    n4 = neghap.alist[ni].allele1[0];
//    	    pintf("\nFFF %d %d %c %c %c %c\n", plushap.alist[pi].pos, neghap.alist[ni].pos, p1, p2, n3, n4);
	    if(plushap.alist[pi].pos + 1 == neghap.alist[ni].pos) {
		if(plushap.alist[pi].ref[0] == 'C') {
		    (*f13) += SCORE_M(p1, n3, p2, n4);
                    (*f14) += SCORE_M(p1, n4, p2, n3);
//printf("sCORE %d %d", *f13, *f14);
    		    pi++;
		    continue;
		}else {
		    pi++;
		    continue;
		}
	    }
	    strcpy(plushap.alist[pi].type, "SNV"); strcpy(neghap.alist[ni].type, "SNV");
    	    (*f13) += SCORE(p1, n3, p2, n4);
    	    (*f14) += SCORE(p1, n4, p2, n3);
	    pi++; ni++;
    }

}
/*
#plus neg 
C2T null plus_methylation
C2T C2T  snv+(methylation?)
null C2T neg_methylation
*/
void print_hap(struct hapfrag* plushaps, struct hapfrag* neghaps, int *ps, int *pe, int *ns, int *ne, int* currentprint) {
	int pi = 0, ni = 0;
	int currentloci = *currentprint;
//	printf("\n=== %d %d %d %d %d %d\n", *ps, *pe, *ns, *ne, plushaps[*ps].printstate, neghaps[*ns].printstate);
	if (*ps <= *pe && (*ns > *ne || plushaps[*ps].printstate == 0) ){
    	while (*ps <= *pe && *ps >= 0 && (plushaps[*ps].printstate == 0 || *ns > *ne) ){
        	pi = 0;
        	while (pi < plushaps[*ps].snps){
            	if (plushaps[*ps].printstate == -1){
            		if(plushaps[*ps].alist[pi].pos > *currentprint)
            			printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap2, plushaps[*ps].alist[pi].hap1,
            				plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
            	}
                else if(plushaps[*ps].printstate == 1){//we need more work here.
                	if(plushaps[*ps].alist[pi].pos > *currentprint)
           				printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap1, plushaps[*ps].alist[pi].hap2,
           					plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
                }else{
                	printf("%s\t%d\t%d\t%d\t%s\t%s\t+\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap1, plushaps[*ps].alist[pi].hap2,
           				plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
                }
                currentloci = plushaps[*ps].alist[pi].pos;
                pi++;
            }
        	(*ps)++;
        }
        if(*ns > *ne)
    		return;
	}
	if (*ns <= *ne && (*ps > *pe || neghaps[*ns].printstate==0) ){
    	while (*ns <= *ne && *ns >= 0 && (neghaps[*ns].printstate==0 || *ps > *pe)){
        	ni = 0;
        	while (ni < neghaps[*ns].snps){
            	if (neghaps[*ns].printstate == -1) {
            		if(neghaps[*ns].alist[ni].pos > *currentprint)
	            		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap2, neghaps[*ns].alist[ni].hap1,
    	        			neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
                }else if (neghaps[*ns].printstate == 1){
                	if(neghaps[*ns].alist[ni].pos > *currentprint)
	            		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap1, neghaps[*ns].alist[ni].hap2,
    	        			neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
                }else{
                	printf("%s\t%d\t%d\t%d\t%s\t%s\t-\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap1, neghaps[*ns].alist[ni].hap2,
            			neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
                }
                currentloci = neghaps[*ns].alist[ni].pos;
                ni++;
            }
        	(*ns)++;
        }
    	return;
    }

    //overlap and merged
	while (pi < plushaps[*ps].snps || ni < neghaps[*ns].snps) {
		while (pi >= plushaps[*ps].snps && ni < neghaps[*ns].snps) {
        	if (neghaps[*ns].printstate == -1) {
        		if(neghaps[*ns].alist[ni].pos > *currentprint)
	           		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap2, neghaps[*ns].alist[ni].hap1,
    	       			neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
       		}else{ 
       			if(neghaps[*ns].alist[ni].pos > *currentprint)
	           		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap1, neghaps[*ns].alist[ni].hap2,
    	   				neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
            }
            currentloci = neghaps[*ns].alist[ni].pos;
        	ni++;
	    }
    	while (ni >= neghaps[*ns].snps && pi < plushaps[*ps].snps) {
        	if (plushaps[*ps].printstate == -1) {
        		if(plushaps[*ps].alist[pi].pos > *currentprint)
	           		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap2, plushaps[*ps].alist[pi].hap1,
    	       			plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
        	}else{ 
        		if(plushaps[*ps].alist[pi].pos > *currentprint)
	           		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap1, plushaps[*ps].alist[pi].hap2,
    	       			plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
            }
            currentloci = plushaps[*ps].alist[pi].pos;
        	pi++;
        }
        if(ni >= neghaps[*ns].snps && pi >= plushaps[*ps].snps)
        	break;
// printf("\nCCC %d %d %d %d\n", pi, ni, plushaps[*ps].alist[pi].pos, neghaps[*ns].alist[ni].pos);

    	while (pi < plushaps[*ps].snps && plushaps[*ps].alist[pi].pos < neghaps[*ns].alist[ni].pos) {
        	if (plushaps[*ps].printstate == -1) {
        		if(plushaps[*ps].alist[pi].pos > *currentprint)
	           		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap2, plushaps[*ps].alist[pi].hap1,
    	       			plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
        	}else{ 
        		if(plushaps[*ps].alist[pi].pos > *currentprint)
	           		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap1, plushaps[*ps].alist[pi].hap2,
    	       			plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
            }
            currentloci = plushaps[*ps].alist[pi].pos;
        	pi++;
        }
    	while (ni < neghaps[*ns].snps  && plushaps[*ps].alist[pi].pos > neghaps[*ns].alist[ni].pos) {
        	if (neghaps[*ns].printstate == -1) {
        		if(neghaps[*ns].alist[ni].pos > *currentprint)
	           		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap2, neghaps[*ns].alist[ni].hap1,
    	       			neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
       		}else{ 
       			if(neghaps[*ns].alist[ni].pos > *currentprint)
           			printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap1, neghaps[*ns].alist[ni].hap2,
       					neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var,neghaps[*ns].alist[ni].type);
            }
            currentloci = neghaps[*ns].alist[ni].pos;
        	ni++;
        }
    	while (pi < plushaps[*ps].snps && ni < neghaps[*ns].snps && plushaps[*ps].alist[pi].pos == neghaps[*ns].alist[ni].pos) {
        	if (plushaps[*ps].alist[pi].ref[0] == 'C') {
	        	if (neghaps[*ns].printstate == -1) {
	        		if(neghaps[*ns].alist[ni].pos > *currentprint)
	        	    	printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap2, neghaps[*ns].alist[ni].hap1,
    	        			neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
   	    		}else{ 
   	    			if(neghaps[*ns].alist[ni].pos > *currentprint)
	        	    	printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", neghaps[*ns].chrom, neghaps[*ns].alist[ni].pos, neghaps[*ns].alist[ni].hap1, neghaps[*ns].alist[ni].hap2,
    	       				neghaps[*ns].alist[ni].ref, neghaps[*ns].alist[ni].var, neghaps[*ns].alist[ni].type);
            	}
            	currentloci = neghaps[*ns].alist[ni].pos;
        	}else{ 
            	if (plushaps[*ps].printstate == -1){
            		if(plushaps[*ps].alist[pi].pos > *currentprint)
	            		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap2, plushaps[*ps].alist[pi].hap1,
    	        			plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
            	}
                else {//we need more work here.
                	if(plushaps[*ps].alist[pi].pos > *currentprint)
	            		printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n", plushaps[*ps].chrom, plushaps[*ps].alist[pi].pos, plushaps[*ps].alist[pi].hap1, plushaps[*ps].alist[pi].hap2,
    	        			plushaps[*ps].alist[pi].ref, plushaps[*ps].alist[pi].var, plushaps[*ps].alist[pi].type);
                }
                currentloci = plushaps[*ps].alist[pi].pos;
        	}
        	pi++; ni++;
        }
    }

    *currentprint = currentloci;
	(*ps)++; (*ns)++;
	if (*ps <= *pe || *ns <= *ne)
        print_hap(plushaps, neghaps, ps, pe, ns, ne, currentprint);
}



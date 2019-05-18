#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <map>
#include <vector>
#include <algorithm>
#include <limits.h>
#define BATBUF 3000

FILE* File_Open(const char* File_Name,const char* Mode)
{
        FILE* Handle;
        Handle=fopen(File_Name,Mode);
        if (Handle==NULL)
        {
                printf("File %s Cannot be opened ....",File_Name);
                exit(1);
        }
        else return Handle;
}
struct Link{
    unsigned pos;
    char hap1;
    char hap2;
    int ismr; //0 no 1 yes 2 possible
    int datatype;
};
struct Block{
    char chr[100];
    std::vector<Link> adjanct;
    unsigned start;
    unsigned end;
    int datatype; //1 wgbs, 2 hic
};
bool comp(const Block &a, const Block &b){
    if( strcmp(a.chr, b.chr) == 0 ) return a.start< b.start;
    const char *a1=a.chr;
    const char *a2=b.chr;
    while(*a1==*a2)
    {
        a1++;
        a2++;
    }
    return *a1<*a2;
}
bool countscore(Block &a, Block &b, int &score, int& nn);

int main(int argc,char *argv[])
{
    if(argc<3){
        fprintf(stderr, "Usage:\n    program plushaplo neghaplo > merge.haplo\n");
        exit(0);
    }
    printf("[Methyhaplo] %s\n", argv[1]);
    printf("[Methyhaplo] %s\n", argv[2]);
    char methyhaploIN[100];strcpy(methyhaploIN, argv[1]);
    char hicIN[100];strcpy(hicIN, argv[2]);
    std::vector<Block> Hap;
    Block tempblk;tempblk.start=0;tempblk.end=0;
    char Dummy[BATBUF];
    unsigned pos=0;int blk1,blk2;
    Link linktmp;char chrom[100];char blockcase[100];
    char haps1[10];char haps2[10];
    fprintf(stderr, "Read hap1 %s\n", methyhaploIN);
    FILE* BSfile = File_Open(methyhaploIN,"r");
    int varline=0,perline=0;
    while(fgets(Dummy,BATBUF,BSfile)!=0){
        sscanf(Dummy, "%*s%s",blockcase);
        if(Dummy[0]=='*' || blockcase[0]=='-'){
            if(tempblk.start!=0 && tempblk.end!=0){
                tempblk.datatype=1;
                Hap.push_back(tempblk);
                tempblk.adjanct.clear();
                tempblk.start=0;
            }
            continue;
        }
        if(Dummy[0]=='B') continue;
        varline++;perline++;
        if(perline>=100000){
            perline=0;
            fprintf(stderr, "Processed %d\n", varline);
        }
        //1       1       0       chrM    327     C       T       0/1
        sscanf(Dummy,"%*s%d%d%s%d%s%s", &blk1, &blk2, chrom, &pos, haps1, haps2);
        if(haps1==NULL || haps2==NULL){
            fprintf(stderr, "%s %d %s %s\n", chrom, pos, haps1, haps2);
        }
        if(blk1==0) 
        {
            linktmp.pos=pos;
            linktmp.hap1=haps1[0];
            linktmp.hap2=haps2[0];
        }
        else if(blk1==1){
            linktmp.pos=pos;
            linktmp.hap1=haps2[0];
            linktmp.hap2=haps1[0];
        }else continue;
        linktmp.datatype=1;
        if(haps1[0] == 'C' && haps2[0]=='T'){
            linktmp.ismr=2;
        }
        tempblk.adjanct.push_back(linktmp);

        if(tempblk.start==0){
            tempblk.start=pos;
            strcpy(tempblk.chr, chrom);
        }else{
            tempblk.end=pos;
        }

    }
    if(tempblk.start!=0 && tempblk.end!=0){
        tempblk.datatype=1;
        Hap.push_back(tempblk);
        tempblk.adjanct.clear();
        tempblk.start=0;
    }
    fclose(BSfile);

    fprintf(stderr, "Read hap2 %s\n", hicIN);
    //hic hap
    tempblk.start=0;
    varline=0;perline=0;
    FILE* hicfile = File_Open(hicIN,"r");
    while(fgets(Dummy,BATBUF,hicfile)!=0){
        sscanf(Dummy, "%*s%s",blockcase);
        if(Dummy[0]=='*' || blockcase[0]=='-'){
            if(tempblk.start!=0 && tempblk.end!=0){
                tempblk.datatype=2;
                Hap.push_back(tempblk);
                tempblk.adjanct.clear();
                tempblk.start=0;
            }
            continue;
        }
        if(Dummy[0]=='B') continue;
        varline++;perline++;
        if(perline>=100000){
            perline=0;
            fprintf(stderr, "Processed %d\n", varline);
        }
        //1       0       1       1       642     C       T       0/1     0       .       100.00
        sscanf(Dummy,"%*s%d%d%s%d%s%s", &blk1, &blk2, chrom, &pos, haps1, haps2);
        if(blk1==0) 
        {
            linktmp.pos=pos;
            linktmp.hap1=haps1[0];
            linktmp.hap2=haps2[0];
        }
        else if(blk1==1){
            linktmp.pos=pos;
            linktmp.hap1=haps2[0];
            linktmp.hap2=haps1[0];
        }else continue;
        if(haps1[0] == 'G' && haps2[0]=='A'){
            linktmp.ismr=2;
        }
        linktmp.datatype=2;
        tempblk.adjanct.push_back(linktmp);

        if(tempblk.start==0){
            tempblk.start=pos;
            strcpy(tempblk.chr, chrom);
        }else{
            tempblk.end=pos;
        }

    }
    if(tempblk.start!=0 && tempblk.end!=0){
        tempblk.datatype=2;
        Hap.push_back(tempblk);
        tempblk.adjanct.clear();
        tempblk.start=0;
    }
    fclose(hicfile);

    fprintf(stderr, "Process hap\n");
    sort(Hap.begin(), Hap.end(), comp);
    fprintf(stderr, "Sorted hap\n");
    fprintf(stderr, "merged hap\n");
    int count = Hap.size();
    Block blkmerge;blkmerge.start=0;blkmerge.end=0;
    Block blktmp;int score=0, nn=0;
    for(int i=0; i<count; i++){
        blktmp = Hap[i];
        if(blkmerge.start==0) {
            blkmerge=blktmp;
            continue;
        }
        if(blkmerge.end < blktmp.start){
            if(blkmerge.start<blkmerge.end){
                printf("Block\t%s\t%d\t%d\n", blkmerge.chr, blkmerge.start, blkmerge.end);
                std::vector<Link> linkout = blkmerge.adjanct;
                for(int i=0; i<linkout.size(); i++){
                    printf("%s\t%d\t%c\t%c\t%d\n", blkmerge.chr, linkout[i].pos, linkout[i].hap1, linkout[i].hap2, linkout[i].datatype);
                }
            }
            blkmerge.adjanct.clear();
            blkmerge=blktmp;
            score=0;nn=0;
            continue;
        }else{
            bool mergett = countscore(blkmerge, blktmp,score, nn);
            if(mergett){
                score=0; nn=0;
                continue;
            }
            else{
                if(blkmerge.start<blkmerge.end){
                    printf("Block\t%s\t%d\t%d\n", blkmerge.chr, blkmerge.start, blkmerge.end);
                    std::vector<Link> linkout = blkmerge.adjanct;
                    for(int i=0; i<linkout.size(); i++){
                        printf("%s\t%d\t%c\t%c\t%d\n", blkmerge.chr, linkout[i].pos, linkout[i].hap1, linkout[i].hap2, linkout[i].datatype);
                    }
                }
                blkmerge.adjanct.clear();
                blkmerge=blktmp;
                score=0;nn=0;
                continue;
            }
        }
    }
    //print end blkmerge
    if(blkmerge.start<blkmerge.end){
        printf("Block\t%s\t%d\t%d\n", blkmerge.chr, blkmerge.start, blkmerge.end);
        std::vector<Link> linkout = blkmerge.adjanct;
        for(int i=0; i<linkout.size(); i++){
            printf("%s\t%d\t%c\t%c\t%d\n", blkmerge.chr, linkout[i].pos, linkout[i].hap1, linkout[i].hap2, linkout[i].datatype);
        }
    }
    blkmerge.adjanct.clear();
    Hap.clear();
}

bool compLink(const Link &a, const Link &b){
    return a.pos< b.pos;
}
void removedup(Block & a){
    int count = a.adjanct.size();
    sort(a.adjanct.begin(), a.adjanct.end(), compLink);
    Block tmp;
    unsigned pos=0;
    tmp.start=a.start; tmp.end=a.end; strcpy(tmp.chr, a.chr);
    tmp.datatype = 3;

    for(int i=0; i<count; i++){
        if(pos!=a.adjanct[i].pos){
            tmp.adjanct.push_back(a.adjanct[i]);
        }
        pos=a.adjanct[i].pos;
    }
    std::swap(a, tmp);
}
bool countscore(Block &a, Block &b, int &score, int &nn){
    int tt=2;int ff=-2;int ss=0;
    int mrs=2;int mrd=-2;
    int start = a.start;
    int end = a.end;
    if(start < b.start) start = b.start;
    if(end > b.end) end = b.end;
    int size_a=a.adjanct.size(); 
    int size_b=b.adjanct.size();
    int last_score=0;
    std::vector <Link> Linkmerge;
    for(int i=0; i<size_b; i++){
        Linkmerge.push_back(b.adjanct[i]);
    }
    for(int i=0; i<size_a; i++){
        Linkmerge.push_back(a.adjanct[i]);
    }
    sort(Linkmerge.begin(), Linkmerge.end(), compLink);
    Link tmp_last,tmp_now; tmp_last.pos=0;tmp_now.pos=0;
    for(int i=0;i<size_a+size_b;i++){
        tmp_now = Linkmerge[i];
        if(tmp_last.pos==0){
            tmp_last = Linkmerge[i];
            continue;
        }
        if(tmp_last.pos < start-1 || tmp_last.pos > end+1 ) {
            tmp_last = Linkmerge[i];
            continue;
        }
        //meth
        int mrpp=1;
        if(tmp_last.pos + 1 == tmp_now.pos && tmp_last.datatype == 1 && tmp_now.datatype == 2){
            if(tmp_last.hap1 == 'C' && tmp_last.hap2 == 'T'){
                if(tmp_now.hap1 == 'G' && tmp_now.hap2 == 'A'){
                    ss += mrs;
                    nn++; 
                }else if(tmp_now.hap1 == 'A' && tmp_now.hap2 == 'G'){
                    ss += mrd;
                    nn++;
                }
            }
            else if(tmp_last.hap1 == 'T' && tmp_last.hap2 == 'C'){
                if(tmp_now.hap1 == 'A' && tmp_now.hap2 == 'G'){
                    ss += mrs;
                    nn++; 
                }else if(tmp_now.hap1 == 'G' && tmp_now.hap2 == 'A'){
                    ss += mrd;
                    nn++;
                }
            }else
                mrpp=-1;
            /* 可以放开这部分
            if(mrpp == 1 && last_score == ff){
                ss -= ff;
                nn--;
            }
            */
        }
        //snv
        if(tmp_last.datatype != tmp_now.datatype && tmp_last.pos == tmp_now.pos){
            nn++;
            if(tmp_last.hap1 == tmp_now.hap1 && tmp_last.hap2 == tmp_now.hap2){
                ss += tt;
                last_score = tt;
            }else if(tmp_last.hap1 == tmp_now.hap2 && tmp_last.hap2 == tmp_now.hap1){
                ss += ff;
                last_score = ff;
            }
        }else
            last_score = 0;
        tmp_last = Linkmerge[i];
    }
    score = ss;
    if(ss != 0 ){
        for(int i=0; i<size_b; i++){
            if(ss<0){
                int tmphap = b.adjanct[i].hap1;
                b.adjanct[i].hap1 = b.adjanct[i].hap2;
                b.adjanct[i].hap2 = tmphap;
            }
            a.adjanct.push_back(b.adjanct[i]);
            removedup(a);
            if(a.start > b.start) a.start = b.start;
            if(a.end < b.end) a.end = b.end;
        }
        return true;
    }
    return false;
}

/*sort默认从小到大排序.
priority_queue默认大的优先,所以队首的是最大的元素.
因为默认优先级相反,所以定义的方式也刚好相反.*/

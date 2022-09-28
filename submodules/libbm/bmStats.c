#include "binaMeth.h"
#include "bmCommon.h"
#include <errno.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>

//Returns -1 if there are no applicable levels, otherwise an integer indicating the most appropriate level.
//Like Kent's library, this divides the desired bin size by 2 to minimize the effect of blocks overlapping multiple bins
static int32_t determineZoomLevel(binaMethFile_t *fp, int basesPerBin) {
    int32_t out = -1;
    int64_t diff;
    //uint32_t bestDiff = -1;
    uint32_t bestDiff = 1000000;
    uint16_t i;

    basesPerBin/=2;
    for(i=0; i<fp->hdr->nLevels; i++) {
        diff = basesPerBin - (int64_t) fp->hdr->zoomHdrs->level[i];
        if(diff >= 0 && diff < bestDiff) {
            bestDiff = diff;
            out = i;
        }
    }
    return out;
}

/// @cond SKIP
struct val_t {
    uint32_t nBases;
    float min, max, sum, sumsq;
    double scalar;
};

struct vals_t {
    uint32_t n;
    struct val_t **vals;
};
/// @endcond

void destroyVals_t(struct vals_t *v) {
    uint32_t i;
    if(!v) return;
    for(i=0; i<v->n; i++) free(v->vals[i]);
    if(v->vals) free(v->vals);
    free(v);
}

//Determine the base-pair overlap between an interval and a block
double getScalar(uint32_t i_start, uint32_t i_end, uint32_t b_start, uint32_t b_end) {
    double rv = 0.0;
    if(b_start <= i_start) {
        if(b_end > i_start) rv = ((double)(b_end - i_start))/(b_end-b_start);
    } else if(b_start < i_end) {
        if(b_end < i_end) rv = ((double)(b_end - b_start))/(b_end-b_start);
        else rv = ((double)(i_end - b_start))/(b_end-b_start);
    }

    return rv;
}

//Returns NULL on error
static struct vals_t *getVals(binaMethFile_t *fp, bmOverlapBlock_t *o, int i, uint32_t tid, uint32_t start, uint32_t end) {
    void *buf = NULL, *compBuf = NULL;
    uLongf sz = fp->hdr->bufSize;
    int compressed = 0, rv;
    uint32_t *p, vtid, vstart, vend;
    struct vals_t *vals = NULL;
    struct val_t *v = NULL;

    if(sz) {
        compressed = 1;
        buf = malloc(sz); 
    }
    sz = 0; //This is now the size of the compressed buffer

    if(bmSetPos(fp, o->offset[i])) goto error;

    vals = calloc(1,sizeof(struct vals_t));
    if(!vals) goto error;

    v = malloc(sizeof(struct val_t));
    if(!v) goto error;

    if(sz < o->size[i]) compBuf = malloc(o->size[i]);
    if(!compBuf) goto error;

    if(bmRead(compBuf, o->size[i], 1, fp) != 1) goto error;
    if(compressed) {
        sz = fp->hdr->bufSize;
        rv = uncompress(buf, &sz, compBuf, o->size[i]);
        if(rv != Z_OK) goto error;
    } else {
        buf = compBuf;
        sz = o->size[i];
    }

    p = buf;
    while(((uLongf) ((void*)p-buf)) < sz) {
        vtid = p[0];
        vstart = p[1];
        vend = p[2];
        v->nBases = p[3];
        v->min = ((float*) p)[4];
        v->max = ((float*) p)[5];
        v->sum = ((float*) p)[6];
        v->sumsq = ((float*) p)[7];
        v->scalar = getScalar(start, end, vstart, vend);

        if(tid == vtid) {
            if((start <= vstart && end > vstart) || (start < vend && start >= vstart)) {
                vals->vals = realloc(vals->vals, sizeof(struct val_t*)*(vals->n+1));
                if(!vals->vals) goto error;
                vals->vals[vals->n++] = v;
                v = malloc(sizeof(struct val_t));
                if(!v) goto error;
            }
            if(vstart > end) break;
        } else if(vtid > tid) {
            break;
        }
        p+=8;
    }

    free(v);
    free(buf);
    if(compressed) free(compBuf);
    return vals;

error:
    if(buf) free(buf);
    if(compBuf && compressed) free(compBuf);
    if(v) free(v);
    destroyVals_t(vals);
    return NULL;
}

//On error, errno is set to ENOMEM and NaN is returned (though NaN can be returned normally)
static double blockMean(binaMethFile_t *fp, bmOverlapBlock_t *blocks, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, j;
    double output = 0.0, coverage = 0.0;
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            output += v->vals[j]->sum * v->vals[j]->scalar;
            coverage += v->vals[j]->nBases * v->vals[j]->scalar;
        }
        destroyVals_t(v);
    }


    if(!coverage) return strtod("NaN", NULL);

    return output/coverage;

error:
    if(v) free(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intMean(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version, uint8_t strand, uint8_t context) {
    double sum = 0.0;
    uint32_t nBases = 0, i, start_use, end_use;

    if(!ints->l) return strtod("NaN", NULL);

    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        if(version & BM_END){
            end_use = ints->end[i];
            if(ints->start[i] < start) start_use = start;
            if(ints->end[i] > end) end_use = end;
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if((ints->strand[i] == strand || strand == 2) 
                    && (ints->context[i] == context || context == 0) ){
                    nBases += end_use-start_use;
                    sum += (end_use-start_use)*((double) ints->value[i]); 
                }
            }else if(version & BM_STRAND){
                if(ints->strand[i] == strand || strand == 2){
                    nBases += end_use-start_use;
                    sum += (end_use-start_use)*((double) ints->value[i]); 
                }
            }else if(version & BM_CONTEXT){
                if(ints->context[i] == context || context == 0){
                    nBases += end_use-start_use;
                    sum += (end_use-start_use)*((double) ints->value[i]); 
                }
            }else{
                nBases += end_use-start_use;
                sum += (end_use-start_use)*((double) ints->value[i]);
            }
        }else{
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if((ints->strand[i] == strand || strand == 2) 
                    && (ints->context[i] == context || context == 0) ){
                    nBases += 1;
                    sum += ((double) ints->value[i]); 
                }
            }else if(version & BM_STRAND){
                if(ints->strand[i] == strand || strand == 2){
                    nBases += 1;
                    sum += ((double) ints->value[i]);
                }
            }else if(version & BM_CONTEXT){
                if(ints->context[i] == context || context == 0){
                    nBases += 1;
                    sum += ((double) ints->value[i]);
                }
            }else{
                nBases += 1;
                sum += ((double) ints->value[i]);
            }
        }
    }

    return sum/nBases;
}

double *intMean_array(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version, uint8_t strand) {
    uint32_t i, start_use, end_use;
    //0 c, 1 cg, 2 chg, 3 chh
    unsigned int Tsize = 4;
    double *sum = malloc(sizeof(double)*Tsize);
    uint32_t *nBases = malloc(sizeof(uint32_t)*Tsize);
    for(i=0;i<Tsize;i++) nBases[i] = 0;
    for(i=0;i<Tsize;i++) sum[i] = 0;
    double *output = malloc(sizeof(double)*Tsize);

    if(!ints->l) {
        //return strtod("NaN", NULL);
        for(i=0;i<Tsize;i++){
            output[i] = strtod("NaN", NULL);
        }
        return output;
    }

    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        if(version & BM_END){
            end_use = ints->end[i];
            if(ints->start[i] < start) start_use = start;
            if(ints->end[i] > end) end_use = end;
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if(ints->strand[i] == strand || strand == 2){
                    nBases[ints->context[i]] += end_use-start_use;
                    sum[ints->context[i]] += (end_use-start_use)*((double) ints->value[i]);
                    if(ints->context[i]!=0){
                        nBases[0] += end_use-start_use;
                        sum[0] += (end_use-start_use)*((double) ints->value[i]);
                    }
                }
            }else if(version & BM_STRAND){
                if(ints->strand[i] == strand || strand == 2){
                    nBases[0] += end_use-start_use;
                    sum[0] += (end_use-start_use)*((double) ints->value[i]); 
                }
            }else if(version & BM_CONTEXT){
                nBases[ints->context[i]] += end_use-start_use;
                sum[ints->context[i]] += (end_use-start_use)*((double) ints->value[i]);
                if(ints->context[i]!=0){
                    nBases[0] += end_use-start_use;
                    sum[0] += (end_use-start_use)*((double) ints->value[i]);
                }
            }else{
                nBases[0] += end_use-start_use;
                sum[0] += (end_use-start_use)*((double) ints->value[i]);
            }
        }else{
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if(ints->strand[i] == strand || strand == 2){
                    nBases[ints->context[i]] += 1;
                    sum[ints->context[i]] += ((double) ints->value[i]);
                    if(ints->context[i]!=0){
                        nBases[0] += 1;
                        sum[0] += ((double) ints->value[i]); 
                    }
                }
            }else if(version & BM_STRAND){
                if(ints->strand[i] == strand || strand == 2){
                    nBases[0] += 1;
                    sum[0] += ((double) ints->value[i]);
                }
            }else if(version & BM_CONTEXT){
                nBases[ints->context[i]] += 1;
                sum[ints->context[i]] += ((double) ints->value[i]);
                if(ints->context[i]!=0){
                    nBases[0] += 1;
                    sum[0] += ((double) ints->value[i]); 
                }
            }else{
                nBases[0] += 1;
                sum[0] += ((double) ints->value[i]);
            }
        }
    }

    for(i=0;i<Tsize;i++){
        if(nBases[i] <= 0) output[i] = strtod("NaN", NULL);
        else output[i] = sum[i]/nBases[i];
    }
    return output;
}



double *intweightedMean_array(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version, uint8_t strand) {
    uint32_t i, start_use, end_use;
    unsigned int Tsize = 4;
    double *sum = malloc(sizeof(double)*Tsize);
    uint32_t *nBases = malloc(sizeof(uint32_t)*Tsize);
    for(i=0;i<Tsize;i++) nBases[i] = 0;
    for(i=0;i<Tsize;i++) sum[i] = 0;
    uint32_t coverC = 0;
    double *output = malloc(sizeof(double)*Tsize);
    if(!(version & BM_COVER)){
        fprintf(stderr, "Error: bm file without coverage information!!!\n mean average will replace weight coverage.");
        double *Amean = intMean_array(ints, start, end, version, strand);
        return Amean;
    }
    if(!ints->l) {
        //return strtod("NaN", NULL);
        for(i=0;i<Tsize;i++){
            output[i] = strtod("NaN", NULL);
        }
        return output;
    }
    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        coverC = (int)(ints->value[i] * ints->coverage[i]+0.5);
        if(version & BM_END){
            end_use = ints->end[i];
            if(ints->start[i] < start) start_use = start;
            if(ints->end[i] > end) end_use = end;
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if(ints->strand[i] == strand || strand == 2){
                    nBases[ints->context[i]] += (end_use-start_use)*coverC;
                    sum[ints->context[i]] += (end_use-start_use)*ints->coverage[i];
                    if(ints->context[i]!=0){
                        nBases[0] += (end_use-start_use)*coverC;
                        sum[0] += (end_use-start_use)*ints->coverage[i];
                    }
                }
            }else if(version & BM_STRAND){
                if(strand == 2 || ints->strand[i] == strand){
                    nBases[0] += (end_use-start_use)*coverC;
                    sum[0] += (end_use-start_use)*ints->coverage[i];
                }
            }else if(version & BM_CONTEXT){
                nBases[ints->context[i]] += (end_use-start_use)*coverC;
                sum[ints->context[i]] += (end_use-start_use)*ints->coverage[i];
                if(ints->context[i]!=0){
                    nBases[0] += (end_use-start_use)*coverC;
                    sum[0] += (end_use-start_use)*ints->coverage[i];
                }
            }else{
                nBases[0] += (end_use-start_use)*coverC;
                sum[0] += (end_use-start_use)*ints->coverage[i];
            }
        }else{
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if(ints->strand[i] == strand || strand == 2){
                    nBases[ints->context[i]] += coverC;
                    sum[ints->context[i]] += ints->coverage[i];
                    if(ints->context[i]!=0){
                        nBases[0] += coverC;
                        sum[0] += ints->coverage[i];
                    }
                }
            }else if(version & BM_STRAND){
                if(strand == 2 || ints->strand[i] == strand){
                    nBases[0] += coverC;
                    sum[0] += ints->coverage[i];
                }
            }else if(version & BM_CONTEXT){
                nBases[ints->context[i]] += coverC;
                sum[ints->context[i]] += ints->coverage[i];
                if(ints->context[i]!=0){
                    nBases[0] += coverC;
                    sum[0] += ints->coverage[i];
                }
            }else{
                nBases[0] += coverC;
                sum[0] += ints->coverage[i];
            }
        }
    }

    for(i=0;i<Tsize;i++){
        //printf("%f %d %d %d\n", sum[i], nBases[i], ints->l, strand );
        if(sum[i] <= 0) output[i] = strtod("NaN", NULL);
        else output[i] = nBases[i]/sum[i];
    }
    return output;
}

void intweightedMean_array_count(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version, uint8_t strand, uint16_t *countC, uint16_t *countCT) {
    uint32_t i, start_use, end_use;
    unsigned int Tsize = 4;
    uint32_t coverC = 0;
    if(!(version & BM_COVER)){
        fprintf(stderr, "Error: bm file without coverage information!!!\n");
        return;
    }
    for(i=0;i<Tsize;i++){
        countC[i] = 0;
        countCT[i] = 0;
    }
    if(!ints->l) {
        return;
    }
    
    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        coverC = (int)(ints->value[i] * ints->coverage[i]+0.5);
        if(version & BM_END){
            end_use = ints->end[i];
            if(ints->start[i] < start) start_use = start;
            if(ints->end[i] > end) end_use = end;
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if(ints->strand[i] == strand || strand == 2){
                    countC[ints->context[i]] += (end_use-start_use)*coverC;
                    countCT[ints->context[i]] += (end_use-start_use)*ints->coverage[i];
                    if(ints->context[i]!=0){
                        countC[0] += (end_use-start_use)*coverC;
                        countCT[0] += (end_use-start_use)*ints->coverage[i];
                    }
                }
            }else if(version & BM_STRAND){
                if(strand == 2 || ints->strand[i] == strand){
                    countC[0] += (end_use-start_use)*coverC;
                    countCT[0] += (end_use-start_use)*ints->coverage[i];
                }
            }else if(version & BM_CONTEXT){
                countC[ints->context[i]] += (end_use-start_use)*coverC;
                countCT[ints->context[i]] += (end_use-start_use)*ints->coverage[i];
                if(ints->context[i]!=0){
                    countC[0] += (end_use-start_use)*coverC;
                    countCT[0] += (end_use-start_use)*ints->coverage[i];
                }
            }else{
                countC[0] += (end_use-start_use)*coverC;
                countCT[0] += (end_use-start_use)*ints->coverage[i];
            }
        }else{
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if(ints->strand[i] == strand || strand == 2){
                    countC[ints->context[i]] += coverC;
                    countCT[ints->context[i]] += ints->coverage[i];
                    if(ints->context[i]!=0){
                        countC[0] += coverC;
                        countCT[0] += ints->coverage[i];
                    }
                }
            }else if(version & BM_STRAND){
                if(strand == 2 || ints->strand[i] == strand){
                    countC[0] += coverC;
                    countCT[0] += ints->coverage[i];
                }
            }else if(version & BM_CONTEXT){
                countC[ints->context[i]] += coverC;
                countCT[ints->context[i]] += ints->coverage[i];
                if(ints->context[i]!=0){
                    countC[0] += coverC;
                    countCT[0] += ints->coverage[i];
                }
            }else{
                countC[0] += coverC;
                countCT[0] += ints->coverage[i];
            }
        }
    }

    return;
}

static double intweightedMean(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version, uint8_t strand, uint8_t context) {
    double sum = 0.0;
    uint32_t nBases = 0, i, start_use, end_use;
    if(!(version & BM_COVER)){
        fprintf(stderr, "Error: bm file without coverage information!!!\n mean average will replace weight coverage.");
        float Amean = intMean(ints, start, end, version, strand, context);
        return Amean;
    }
    if(!ints->l) return strtod("NaN", NULL);
    for(i=0; i<ints->l; i++) {
        start_use = ints->start[i];
        if(version & BM_END){
            end_use = ints->end[i];
            if(ints->start[i] < start) start_use = start;
            if(ints->end[i] > end) end_use = end;
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if((ints->strand[i] == strand || strand == 2) 
                    && (ints->context[i] == context || context == 0) ){
                    nBases += (end_use-start_use)*(uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                    sum += (end_use-start_use)*ints->coverage[i]; 
                }
            }else if(version & BM_STRAND){
                if(strand == 2 || ints->strand[i] == strand){
                    nBases += (end_use-start_use)*(uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                    sum += (end_use-start_use)*ints->coverage[i];
                }
            }else if(version & BM_CONTEXT){
                if(ints->context[i] == context || context == 0){
                    nBases += (end_use-start_use)*(uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                    sum += (end_use-start_use)*ints->coverage[i];
                }
            }else{
                nBases += (end_use-start_use)*(uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                sum += (end_use-start_use)*ints->coverage[i];
            }
        }else{
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if((ints->strand[i] == strand || strand == 2) 
                    && (ints->context[i] == context || context == 0) ){
                    nBases += (uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                    sum += ints->coverage[i];
                }
            }else if(version & BM_STRAND){
                if(strand == 2 || ints->strand[i] == strand){
                    nBases += (uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                    sum += ints->coverage[i];
                }
            }else if(version & BM_CONTEXT){
                if(ints->context[i] == context || context == 0){
                    nBases += (uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                    sum += ints->coverage[i];
                }
            }else{
                nBases += (uint32_t)((double) ints->value[i] * ints->coverage[i]+0.5);
                sum += ints->coverage[i];
            }
        }
    }

    return nBases/sum;
}

//Does UCSC compensate for partial block/range overlap?
static double blockDev(binaMethFile_t *fp, bmOverlapBlock_t *blocks, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, j;
    double mean = 0.0, ssq = 0.0, coverage = 0.0, diff;
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            coverage += v->vals[j]->nBases * v->vals[j]->scalar;
            mean += v->vals[j]->sum * v->vals[j]->scalar;
            ssq += v->vals[j]->sumsq * v->vals[j]->scalar;
        }
        destroyVals_t(v);
        v = NULL;
    }

    if(coverage<=1.0) return strtod("NaN", NULL);
    diff = ssq-mean*mean/coverage;
    if(coverage > 1.0) diff /= coverage-1;
    if(fabs(diff) > 1e-8) { //Ignore floating point differences
        return sqrt(diff);
    } else {
        return 0.0;
    }

error:
    if(v) destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

//This uses compensated summation to account for finite precision math
static double intDev(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version, uint8_t strand, uint8_t context) {
    double v1 = 0.0, mean, rv;
    uint32_t nBases = 0, i, start_use, end_use;

    if(!ints->l) return strtod("NaN", NULL);
    mean = intMean(ints, start, end, version, strand, context);

    if(version & BM_END){
        for(i=0; i<ints->l; i++) {
            start_use = ints->start[i];
            end_use = ints->end[i];
            if(ints->start[i] < start) start_use = start;
            if(ints->end[i] > end) end_use = end;
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if((ints->strand[i] == strand || strand == 2) 
                    && (ints->context[i] == context || context == 0) ){
                    nBases += end_use-start_use;
                    v1 += (end_use-start_use) * pow(ints->value[i]-mean, 2.0);
                }
            }else if(version & BM_STRAND){
                if(ints->strand[i] == strand || strand == 2){
                    nBases += end_use-start_use;
                    v1 += (end_use-start_use) * pow(ints->value[i]-mean, 2.0); //running sum of squared difference
                }
            }else if(version & BM_CONTEXT){
                if(ints->context[i] == context || context == 0){
                    nBases += end_use-start_use;
                    v1 += (end_use-start_use) * pow(ints->value[i]-mean, 2.0); //running sum of squared difference
                }
            }else{
                nBases += end_use-start_use;
                v1 += (end_use-start_use) * pow(ints->value[i]-mean, 2.0); //running sum of squared difference
            }
        }
    }else{
        for(i=0; i<ints->l; i++) {
            if((version & BM_STRAND) && (version & BM_CONTEXT)){
                if((ints->strand[i] == strand || strand == 2) 
                    && (ints->context[i] == context || context == 0) ){
                    nBases += 1;
                    v1 += pow(ints->value[i]-mean, 2.0);
                }
            }else if(version & BM_STRAND){
                if(ints->strand[i] == strand || strand == 2){
                    nBases += 1;
                    v1 += pow(ints->value[i]-mean, 2.0); //running sum of squared difference
                }
            }else if(version & BM_CONTEXT){
                if(ints->context[i] == context || context == 0){
                    nBases += 1;
                    v1 += pow(ints->value[i]-mean, 2.0); //running sum of squared difference
                }
            }else{
                nBases += 1;
                v1 += pow(ints->value[i]-mean, 2.0); //running sum of squared difference
            }
        }
    }

    if(nBases>=2) rv = sqrt(v1/(nBases-1));
    else if(nBases==1) rv = sqrt(v1);
    else rv = strtod("NaN", NULL);

    return rv;
}

static double blockMax(binaMethFile_t *fp, bmOverlapBlock_t *blocks, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, j, isNA = 1;
    double o = strtod("NaN", NULL);
    struct vals_t *v = NULL;

    if(!blocks->n) return o;

    //Iterate the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            if(isNA) {
                o = v->vals[j]->max;
                isNA = 0;
            } else if(v->vals[j]->max > o) {
                o = v->vals[j]->max;
            }
        }
        destroyVals_t(v);
    }

    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intMax(bmOverlappingIntervals_t* ints) {
    uint32_t i;
    double o;

    if(ints->l < 1) return strtod("NaN", NULL);

    o = ints->value[0];
    for(i=1; i<ints->l; i++) {
        if(ints->value[i] > o) o = ints->value[i];
    }

    return o;
}

static double blockMin(binaMethFile_t *fp, bmOverlapBlock_t *blocks, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, j, isNA = 1;
    double o = strtod("NaN", NULL);
    struct vals_t *v = NULL;

    if(!blocks->n) return o;

    //Iterate the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            if(isNA) {
                o = v->vals[j]->min;
                isNA = 0;
            } else if(v->vals[j]->min < o) o = v->vals[j]->min;
        }
        destroyVals_t(v);
    }

    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intMin(bmOverlappingIntervals_t* ints) {
    uint32_t i;
    double o;

    if(ints->l < 1) return strtod("NaN", NULL);

    o = ints->value[0];
    for(i=1; i<ints->l; i++) {
        if(ints->value[i] < o) o = ints->value[i];
    }

    return o;
}

//Does UCSC compensate for only partial block/interval overlap?
static double blockCoverage(binaMethFile_t *fp, bmOverlapBlock_t *blocks, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, j;
    double o = 0.0;
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            o+= v->vals[j]->nBases * v->vals[j]->scalar;
        }
        destroyVals_t(v);
    }

    if(o == 0.0) return strtod("NaN", NULL);
    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intCoverage(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version) {
    uint32_t i, start_use, end_use;
    double o = 0.0;

    if(!ints->l) return strtod("NaN", NULL);

    if(version & BM_END){
        for(i=0; i<ints->l; i++) {
            start_use = ints->start[i];
            end_use = ints->end[i];
            if(start_use < start) start_use = start;
            if(end_use > end) end_use = end;
            o += end_use - start_use;
        }
    }else{
        for(i=0; i<ints->l; i++) {
            o += 1;
        }
    }

    return o/(end-start);
}

static double blockSum(binaMethFile_t *fp, bmOverlapBlock_t *blocks, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, j, sizeUse;
    double o = 0.0;
    struct vals_t *v = NULL;

    if(!blocks->n) return strtod("NaN", NULL);

    //Iterate over the blocks
    for(i=0; i<blocks->n; i++) {
        v = getVals(fp, blocks, i, tid, start, end);
        if(!v) goto error;
        for(j=0; j<v->n; j++) {
            //Multiply the block average by min(bases covered, block overlap with interval)
            sizeUse = v->vals[j]->scalar;
            if(sizeUse > v->vals[j]->nBases) sizeUse = v->vals[j]->nBases;
            o+= (v->vals[j]->sum * sizeUse) / v->vals[j]->nBases;
        }
        destroyVals_t(v);
    }

    if(o == 0.0) return strtod("NaN", NULL);
    return o;

error:
    destroyVals_t(v);
    errno = ENOMEM;
    return strtod("NaN", NULL);
}

static double intSum(bmOverlappingIntervals_t* ints, uint32_t start, uint32_t end, uint16_t version) {
    uint32_t i, start_use, end_use;
    double o = 0.0;

    if(!ints->l) return strtod("NaN", NULL);

    if(version & BM_END){
        for(i=0; i<ints->l; i++) {
            start_use = ints->start[i];
            end_use = ints->end[i];
            if(start_use < start) start_use = start;
            if(end_use > end) end_use = end;
            o += (end_use - start_use) * ints->value[i];
        }
    }else{
        for(i=0; i<ints->l; i++) {
            o += ints->value[i];
        }
    }

    return o;
}

//Returns NULL on error, otherwise a double* that needs to be free()d
double *bmStatsFromZoom(binaMethFile_t *fp, int32_t level, uint32_t tid, uint32_t start, uint32_t end, uint32_t nBins, enum bmStatsType type) {
    bmOverlapBlock_t *blocks = NULL;
    double *output = NULL;
    uint32_t pos = start, i, end2;

    if(!fp->hdr->zoomHdrs->idx[level]) {
        fp->hdr->zoomHdrs->idx[level] = bmReadIndex(fp, fp->hdr->zoomHdrs->indexOffset[level]);
        if(!fp->hdr->zoomHdrs->idx[level]) return NULL;
    }
    errno = 0; //Sometimes libCurls sets and then doesn't unset errno on errors

    output = malloc(sizeof(double)*nBins);
    if(!output) return NULL;

    for(i=0, pos=start; i<nBins; i++) {
        end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        blocks = walkRTreeNodes(fp, fp->hdr->zoomHdrs->idx[level]->root, tid, pos, end2);
        if(!blocks) goto error;

        switch(type) {
        case 0:
            //mean
            output[i] = blockMean(fp, blocks, tid, pos, end2);
            break;
        case 1:
            //stdev
            output[i] = blockDev(fp, blocks, tid, pos, end2);
            break;
        case 2:
            //max
            output[i] = blockMax(fp, blocks, tid, pos, end2);
            break;
        case 3:
            //min
            output[i] = blockMin(fp, blocks, tid, pos, end2);
            break;
        case 4:
            //cov
            output[i] = blockCoverage(fp, blocks, tid, pos, end2)/(end2-pos);
            break;
        case 5:
            //sum
            output[i] = blockSum(fp, blocks, tid, pos, end2);
            break;
        default:
            goto error;
            break;
        }
        if(errno) goto error;
        destroyBWOverlapBlock(blocks);
        pos = end2;
    }

    return output;

error:
    fprintf(stderr, "got an error in bmStatsFromZoom in the range %"PRIu32"-%"PRIu32": %s\n", pos, end2, strerror(errno));
    if(blocks) destroyBWOverlapBlock(blocks);
    if(output) free(output);
    return NULL;
}

double *bmStatsFromFull(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand, uint8_t context) {
    bmOverlappingIntervals_t *ints = NULL;
    double *output = malloc(sizeof(double)*nBins);
    uint32_t i, pos = start, end2;
    if(!output) return NULL;

    for(i=0; i<nBins; i++) {
        //end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        if(i==0){
            end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        }else{
            end2 = start + movestep*(i+1);
        }
        ints = bmGetOverlappingIntervals(fp, chrom, pos, end2);

        if(!ints) {
            output[i] = strtod("NaN", NULL);
            continue;
        }

        switch(type) {
        default :
        case 0:
            output[i] = intMean(ints, pos, end2, fp->hdr->version, strand, context);
            break;
        case 6:
            output[i] = intweightedMean(ints, pos, end2, fp->hdr->version, strand, context);
            break;
        case 1:
            output[i] = intDev(ints, pos, end2, fp->hdr->version, strand, context);
            break;
        case 2:
            output[i] = intMax(ints);
            break;
        case 3:
            output[i] = intMin(ints);
            break;
        case 4:
            output[i] = intCoverage(ints, pos, end2, fp->hdr->version);
            break;
        case 5:
            output[i] = intSum(ints, pos, end2, fp->hdr->version);
            break;
        }
        bmDestroyOverlappingIntervals(ints);
        //pos = end2;
        pos = pos + movestep;
    }

    return output;
}

double *bmStatsFromFull_array(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand) {
    bmOverlappingIntervals_t *ints = NULL;
    unsigned int Tsize = 4;
    double *output = malloc(sizeof(double)*nBins*Tsize);
    uint32_t i, pos = start, end2, j, k;
    if(!output) return NULL;
    double *outmean = malloc(sizeof(double)*Tsize);

    for(i=0,k=0; i<nBins; i++) {
        //end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        if(i==0){
            end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        }else{
            end2 = start + movestep*(i+1);
        }
        ints = bmGetOverlappingIntervals(fp, chrom, pos, end2);

        if(!ints) {
            for(j=0;j<Tsize;j++){
                output[k+j] = strtod("NaN", NULL);
            }
            continue;
        }

        switch(type) {
            default :
            case 0:
                outmean = intMean_array(ints, pos, end2, fp->hdr->version, strand);
                for(j=0;j<Tsize;j++){
                    output[k+j] = outmean[j];
                }
                break;
            case 6:
                outmean = intweightedMean_array(ints, pos, end2, fp->hdr->version, strand);
                for(j=0;j<Tsize;j++){
                    output[k+j] = outmean[j];
                }
                break;
        }
        k += Tsize;
        bmDestroyOverlappingIntervals(ints);
        //pos = end2;
        pos = pos + movestep;
    }
    free(outmean);
    return output;
}

void bmStatsFromFull_array_count(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand, uint16_t *countC, uint16_t *countCT) {
    bmOverlappingIntervals_t *ints = NULL;
    unsigned int Tsize = 4;
    //double *output = malloc(sizeof(double)*nBins*Tsize);
    uint32_t i, pos = start, end2, j, k;
    if(!countC || !countCT) return;
    uint16_t *out_C = malloc(sizeof(uint16_t)*Tsize);
    uint16_t *out_CT = malloc(sizeof(uint16_t)*Tsize);

    for(i=0,k=0; i<nBins; i++) {
        //end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        if(i==0){
            end2 = start + ((double)(end-start)*(i+1))/((int) nBins);
        }else{
            end2 = start + movestep*(i+1);
        }
        ints = bmGetOverlappingIntervals(fp, chrom, pos, end2);

        if(!ints) {
            for(j=0;j<Tsize;j++){
                countC[k+j] = 0;
                countCT[k+j] = 0;
            }
            continue;
        }

        switch(type) {
            default :
            case 6:
                intweightedMean_array_count(ints, pos, end2, fp->hdr->version, strand, out_C, out_CT);
                for(j=0;j<Tsize;j++){
                    countC[k+j] = out_C[j];
                    countCT[k+j] = out_CT[j];
                }
                break;
        }
        k += Tsize;
        bmDestroyOverlappingIntervals(ints);
        //pos = end2;
        pos = pos + movestep;
    }
    free(out_C); free(out_CT);
    return;
}

//Returns a list of floats of length nBins that must be free()d
//On error, NULL is returned
double *bmStats(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand, uint8_t context) {
    //int32_t level = determineZoomLevel(fp, ((double)(end-start))/((int) nBins));
    uint32_t tid = bmGetTid(fp, chrom);
    if(tid == (uint32_t) -1) return NULL;

    //if(level == -1 || level == 0)
    return bmStatsFromFull(fp, chrom, start, end, nBins, movestep, type, strand, context);
    //return bmStatsFromZoom(fp, level, tid, start, end, nBins, type);
}

double *bmStats_array(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand) {
    //int32_t level = determineZoomLevel(fp, ((double)(end-start))/((int) nBins));
    uint32_t tid = bmGetTid(fp, chrom);
    if(tid == (uint32_t) -1) return NULL;

    //if(level == -1 || level == 0)
    return bmStatsFromFull_array(fp, chrom, start, end, nBins, movestep, type, strand);
    //return bmStatsFromZoom(fp, level, tid, start, end, nBins, type);
}

void bmStats_array_count(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand, uint16_t *countC, uint16_t *countCT) {
    //int32_t level = determineZoomLevel(fp, ((double)(end-start))/((int) nBins));
    uint32_t tid = bmGetTid(fp, chrom);
    if(tid == (uint32_t) -1) return;

    //if(level == -1 || level == 0)
    bmStatsFromFull_array_count(fp, chrom, start, end, nBins, movestep, type, strand, countC, countCT);
    //return bmStatsFromZoom(fp, level, tid, start, end, nBins, type);
    return;
}

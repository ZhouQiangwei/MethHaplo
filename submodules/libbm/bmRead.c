#include "binaMeth.h"
#include "bmCommon.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

static uint64_t readChromBlock(binaMethFile_t *bm, chromList_t *cl, uint32_t keySize);

//Return the position in the file
long bmTell(binaMethFile_t *fp) {
    if(fp->URL->type == BWG_FILE) return ftell(fp->URL->x.fp);
    return (long) (fp->URL->filePos + fp->URL->bufPos);
}

//Seek to a given position, always from the beginning of the file
//Return 0 on success and -1 on error
//To do, use the return code of urlSeek() in a more useful way.
int bmSetPos(binaMethFile_t *fp, size_t pos) {
    CURLcode rv = urlSeek(fp->URL, pos);
    if(rv == CURLE_OK) return 0;
    return -1;
}

//returns the number of full members read (nmemb on success, something less on error)
size_t bmRead(void *data, size_t sz, size_t nmemb, binaMethFile_t *fp) {
    size_t i, rv;
    for(i=0; i<nmemb; i++) {
        rv = urlRead(fp->URL, data+i*sz, sz);
        if(rv != sz) return i;
    }
    return nmemb;
}

//Initializes curl and sets global variables
//Returns 0 on success and 1 on error
//This should be called only once and bmCleanup() must be called when finished.
int bmInit(size_t defaultBufSize) {
    //set the buffer size, number of iterations, sleep time between iterations, etc.
    GLOBAL_DEFAULTBUFFERSIZE = defaultBufSize;

    //call curl_global_init()
#ifndef NOCURL
    CURLcode rv;
    rv = curl_global_init(CURL_GLOBAL_ALL);
    if(rv != CURLE_OK) return 1;
#endif
    return 0;
}

//This should be called before quiting, to release memory acquired by curl
void bmCleanup() {
#ifndef NOCURL
    curl_global_cleanup();
#endif
}

static bmZoomHdr_t *bmReadZoomHdrs(binaMethFile_t *bm) {
    if(bm->isWrite) return NULL;
    uint16_t i;
    bmZoomHdr_t *zhdr = malloc(sizeof(bmZoomHdr_t));
    if(!zhdr) return NULL;
    uint32_t *level = malloc(bm->hdr->nLevels * sizeof(uint64_t));
    if(!level) {
        free(zhdr);
        return NULL;
    }
    uint32_t padding = 0;
    uint64_t *dataOffset = malloc(sizeof(uint64_t) * bm->hdr->nLevels);
    if(!dataOffset) {
        free(zhdr);
        free(level);
        return NULL;
    }
    uint64_t *indexOffset = malloc(sizeof(uint64_t) * bm->hdr->nLevels);
    if(!indexOffset) {
        free(zhdr);
        free(level);
        free(dataOffset);
        return NULL;
    }

    for(i=0; i<bm->hdr->nLevels; i++) {
        if(bmRead((void*) &(level[i]), sizeof(uint32_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &padding, sizeof(uint32_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &(dataOffset[i]), sizeof(uint64_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &(indexOffset[i]), sizeof(uint64_t), 1, bm) != 1) goto error;
    }

    zhdr->level = level;
    zhdr->dataOffset = dataOffset;
    zhdr->indexOffset = indexOffset;
    zhdr->idx = calloc(bm->hdr->nLevels, sizeof(bmRTree_t*));
    if(!zhdr->idx) goto error;

    return zhdr;

error:
    for(i=0; i<bm->hdr->nLevels; i++) {
        if(zhdr->idx[i]) bmDestroyIndex(zhdr->idx[i]);
    }
    free(zhdr);
    free(level);
    free(dataOffset);
    free(indexOffset);
    return NULL;
}

static void bmHdrDestroy(binaMethHdr_t *hdr) {
    int i;
    if(hdr->zoomHdrs) {
        free(hdr->zoomHdrs->level);
        free(hdr->zoomHdrs->dataOffset);
        free(hdr->zoomHdrs->indexOffset);
        for(i=0; i<hdr->nLevels; i++) {
            if(hdr->zoomHdrs->idx[i]) bmDestroyIndex(hdr->zoomHdrs->idx[i]);
        }
        free(hdr->zoomHdrs->idx);
        free(hdr->zoomHdrs);
    }
    free(hdr);
}

static void bmHdrRead(binaMethFile_t *bm) {
    uint32_t magic;
    if(bm->isWrite) return;
    bm->hdr = calloc(1, sizeof(binaMethHdr_t));
    if(!bm->hdr) return;

    if(bmRead((void*) &magic, sizeof(uint32_t), 1, bm) != 1) goto error; //0x0
    if(!(magic & BM_MAGIC) && magic != BIGWIG_MAGIC && magic != BIGBED_MAGIC) goto error;

    if(bmRead((void*) &(bm->hdr->version), sizeof(uint16_t), 1, bm) != 1) goto error; //0x4
    if(bmRead((void*) &(bm->hdr->nLevels), sizeof(uint16_t), 1, bm) != 1) goto error; //0x6
    if(bmRead((void*) &(bm->hdr->ctOffset), sizeof(uint64_t), 1, bm) != 1) goto error; //0x8
    if(bmRead((void*) &(bm->hdr->dataOffset), sizeof(uint64_t), 1, bm) != 1) goto error; //0x10
    if(bmRead((void*) &(bm->hdr->indexOffset), sizeof(uint64_t), 1, bm) != 1) goto error; //0x18
    if(bmRead((void*) &(bm->hdr->fieldCount), sizeof(uint16_t), 1, bm) != 1) goto error; //0x20
    if(bmRead((void*) &(bm->hdr->definedFieldCount), sizeof(uint16_t), 1, bm) != 1) goto error; //0x22
    if(bmRead((void*) &(bm->hdr->sqlOffset), sizeof(uint64_t), 1, bm) != 1) goto error; //0x24
    if(bmRead((void*) &(bm->hdr->summaryOffset), sizeof(uint64_t), 1, bm) != 1) goto error; //0x2c
    if(bmRead((void*) &(bm->hdr->bufSize), sizeof(uint32_t), 1, bm) != 1) goto error; //0x34
    if(bmRead((void*) &(bm->hdr->extensionOffset), sizeof(uint64_t), 1, bm) != 1) goto error; //0x38

    //zoom headers
    if(bm->hdr->nLevels) {
        if(!(bm->hdr->zoomHdrs = bmReadZoomHdrs(bm))) goto error;
    }

    //File summary information
    if(bm->hdr->summaryOffset) {
        if(urlSeek(bm->URL, bm->hdr->summaryOffset) != CURLE_OK) goto error;
        if(bmRead((void*) &(bm->hdr->nBasesCovered), sizeof(uint64_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &(bm->hdr->minVal), sizeof(uint64_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &(bm->hdr->maxVal), sizeof(uint64_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &(bm->hdr->sumData), sizeof(uint64_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &(bm->hdr->sumSquared), sizeof(uint64_t), 1, bm) != 1) goto error;
    }

    //In case of uncompressed remote files, let the IO functions know to request larger chunks
    bm->URL->isCompressed = (bm->hdr->bufSize > 0)?1:0;

    return;

error:
    bmHdrDestroy(bm->hdr);
    fprintf(stderr, "[bmHdrRead] There was an error while reading in the header!\n");
    bm->hdr = NULL;
}

static void destroyChromList(chromList_t *cl) {
    uint32_t i;
    if(!cl) return;
    if(cl->nKeys && cl->chrom) {
        for(i=0; i<cl->nKeys; i++) {
            if(cl->chrom[i]) free(cl->chrom[i]);
        }
    }
    if(cl->chrom) free(cl->chrom);
    if(cl->len) free(cl->len);
    free(cl);
}

static uint64_t readChromLeaf(binaMethFile_t *bm, chromList_t *cl, uint32_t valueSize) {
    uint16_t nVals, i;
    uint32_t idx;
    char *chrom = NULL;

    if(bmRead((void*) &nVals, sizeof(uint16_t), 1, bm) != 1) return -1;
    chrom = calloc(valueSize+1, sizeof(char));
    if(!chrom) return -1;

    for(i=0; i<nVals; i++) {
        if(bmRead((void*) chrom, sizeof(char), valueSize, bm) != valueSize) goto error;
        if(bmRead((void*) &idx, sizeof(uint32_t), 1, bm) != 1) goto error;
        if(bmRead((void*) &(cl->len[idx]), sizeof(uint32_t), 1, bm) != 1) goto error;
        cl->chrom[idx] = strdup(chrom);
        if(!(cl->chrom[idx])) goto error;
    }

    free(chrom);
    return nVals;

error:
    free(chrom);
    return -1;
}

static uint64_t readChromNonLeaf(binaMethFile_t *bm, chromList_t *cl, uint32_t keySize) {
    uint64_t offset , rv = 0, previous;
    uint16_t nVals, i;

    if(bmRead((void*) &nVals, sizeof(uint16_t), 1, bm) != 1) return -1;

    previous = bmTell(bm) + keySize;
    for(i=0; i<nVals; i++) {
        if(bmSetPos(bm, previous)) return -1;
        if(bmRead((void*) &offset, sizeof(uint64_t), 1, bm) != 1) return -1;
        if(bmSetPos(bm, offset)) return -1;
        rv += readChromBlock(bm, cl, keySize);
        previous += 8 + keySize;
    }

    return rv;
}

static uint64_t readChromBlock(binaMethFile_t *bm, chromList_t *cl, uint32_t keySize) {
    uint8_t isLeaf, padding;

    if(bmRead((void*) &isLeaf, sizeof(uint8_t), 1, bm) != 1) return -1;
    if(bmRead((void*) &padding, sizeof(uint8_t), 1, bm) != 1) return -1;

    if(isLeaf) {
        return readChromLeaf(bm, cl, keySize);
    } else { //I've never actually observed one of these, which is good since they're pointless
        return readChromNonLeaf(bm, cl, keySize);
    }
}

static chromList_t *bmReadChromList(binaMethFile_t *bm) {
    chromList_t *cl = NULL;
    uint32_t magic, keySize, valueSize, itemsPerBlock;
    uint64_t rv, itemCount;
    if(bm->isWrite) return NULL;
    if(bmSetPos(bm, bm->hdr->ctOffset)) return NULL;

    cl = calloc(1, sizeof(chromList_t));
    if(!cl) return NULL;

    if(bmRead((void*) &magic, sizeof(uint32_t), 1, bm) != 1) goto error;
    if(magic != CIRTREE_MAGIC) goto error;

    if(bmRead((void*) &itemsPerBlock, sizeof(uint32_t), 1, bm) != 1) goto error;
    if(bmRead((void*) &keySize, sizeof(uint32_t), 1, bm) != 1) goto error;
    if(bmRead((void*) &valueSize, sizeof(uint32_t), 1, bm) != 1) goto error;
    if(bmRead((void*) &itemCount, sizeof(uint64_t), 1, bm) != 1) goto error;

    cl->nKeys = itemCount;
    cl->chrom = calloc(itemCount, sizeof(char*));
    cl->len = calloc(itemCount, sizeof(uint32_t));
    if(!cl->chrom) goto error;
    if(!cl->len) goto error;

    if(bmRead((void*) &magic, sizeof(uint32_t), 1, bm) != 1) goto error;
    if(bmRead((void*) &magic, sizeof(uint32_t), 1, bm) != 1) goto error;

    //Read in the blocks
    rv = readChromBlock(bm, cl, keySize);
    if(rv == (uint64_t) -1) goto error;
    if(rv != itemCount) goto error;

    return cl;

error:
    destroyChromList(cl);
    return NULL;
}

//This is here mostly for convenience
static void bmDestroyWriteBuffer(bmWriteBuffer_t *wb) {
    if(wb->p) free(wb->p);
    if(wb->compressP) free(wb->compressP);
    if(wb->firstZoomBuffer) free(wb->firstZoomBuffer);
    if(wb->lastZoomBuffer) free(wb->lastZoomBuffer);
    if(wb->nNodes) free(wb->nNodes);
    free(wb);
}

void bmClose(binaMethFile_t *fp) {
    if(DEBUG>1) fprintf(stderr, "kkkkxx\n");
    if(!fp) return;
    if(DEBUG>1) fprintf(stderr, "lllllxxxx\n");
    if(bmFinalize(fp)) {
        fprintf(stderr, "[bmClose] There was an error while finishing writing a binaMeth file! The output is likely truncated.\n");
    }
    if(DEBUG>1) fprintf(stderr, "222222222\n");
    if(fp->URL) urlClose(fp->URL);
    if(fp->hdr) bmHdrDestroy(fp->hdr);
    if(fp->cl) destroyChromList(fp->cl);
    if(fp->idx) bmDestroyIndex(fp->idx);
    if(fp->writeBuffer) bmDestroyWriteBuffer(fp->writeBuffer);
    free(fp);
}

int bmIsBinaMeth(char *fname, CURLcode (*callBack) (CURL*)) {
    uint32_t magic = 0;
    URL_t *URL = NULL;

    URL = urlOpen(fname, *callBack, NULL);

    if(!URL) return 0;
    if(urlRead(URL, (void*) &magic, sizeof(uint32_t)) != sizeof(uint32_t)) magic = 0;
    urlClose(URL);
    if(magic == BIGWIG_MAGIC) return 1;
    return 0;
}

uint32_t BMtype(char *fname, CURLcode (*callBack) (CURL*)) {
    uint32_t magic = 0;
    URL_t *URL = NULL;

    URL = urlOpen(fname, *callBack, NULL);

    if(!URL) return 0;
    if(urlRead(URL, (void*) &magic, sizeof(uint32_t)) != sizeof(uint32_t)) magic = 0;
    urlClose(URL);
    //if(magic == BIGWIG_MAGIC) return 1;
    return magic;
}

char *bbGetSQL(binaMethFile_t *bm) {
    char *o = NULL;
    uint64_t len;
    if(!bm->hdr->sqlOffset) return NULL;
    len = bm->hdr->summaryOffset - bm->hdr->sqlOffset; //This includes the NULL terminator
    o = malloc(sizeof(char) * len);
    if(!o) goto error;
    if(bmSetPos(bm, bm->hdr->sqlOffset)) goto error;
    if(bmRead((void*) o, len, 1, bm) != 1) goto error;
    return o;

error:
    if(o) free(o);
    printf("Got an error in bbGetSQL!\n");
    return NULL;
}

int bbIsBigBed(char *fname, CURLcode (*callBack) (CURL*)) {
    uint32_t magic = 0;
    URL_t *URL = NULL;

    URL = urlOpen(fname, *callBack, NULL);

    if(!URL) return 0;
    if(urlRead(URL, (void*) &magic, sizeof(uint32_t)) != sizeof(uint32_t)) magic = 0;
    urlClose(URL);
    if(magic == BIGBED_MAGIC) return 1;
    return 0;
}

binaMethFile_t *bmOpen(char *fname, CURLcode (*callBack) (CURL*), const char *mode) {
    binaMethFile_t *bmg = calloc(1, sizeof(binaMethFile_t));
    if(!bmg) {
        fprintf(stderr, "[bmOpen] Couldn't allocate space to create the output object!\n");
        return NULL;
    }
    if((!mode) || (strchr(mode, 'w') == NULL)) {
        bmg->isWrite = 0;
        bmg->URL = urlOpen(fname, *callBack, NULL);
        if(!bmg->URL) {
            fprintf(stderr, "[bmOpen] urlOpen is NULL!\n");
            goto error;
        }

        //Attempt to read in the fixed header
        bmHdrRead(bmg);
        if(!bmg->hdr) {
            fprintf(stderr, "[bmOpen] bmg->hdr is NULL!\n");
            goto error;
        }

        //Read in the chromosome list
        bmg->cl = bmReadChromList(bmg);
        if(!bmg->cl) {
            fprintf(stderr, "[bmOpen] bmg->cl is NULL (%s)!\n", fname);
            goto error;
        }

        //Read in the index
        if(bmg->hdr->indexOffset) {
            bmg->idx = bmReadIndex(bmg, 0);
            if(!bmg->idx) {
                fprintf(stderr, "[bmOpen] bmg->idx is NULL bmg->hdr->dataOffset 0x%"PRIx64"!\n", bmg->hdr->dataOffset);
                goto error;
            }
        }
        //type used for meth type
        bmg->type = bmg->hdr->version;
    } else {
        bmg->isWrite = 1;
        bmg->URL = urlOpen(fname, NULL, "w+");
        if(!bmg->URL) goto error;
        bmg->writeBuffer = calloc(1,sizeof(bmWriteBuffer_t));
        if(!bmg->writeBuffer) goto error;
        bmg->writeBuffer->l = 24;
    }

    return bmg;

error:
    bmClose(bmg);
    return NULL;
}

binaMethFile_t *bbOpen(char *fname, CURLcode (*callBack) (CURL*)) {
    binaMethFile_t *bb = calloc(1, sizeof(binaMethFile_t));
    if(!bb) {
        fprintf(stderr, "[bbOpen] Couldn't allocate space to create the output object!\n");
        return NULL;
    }

    //Set the type to 1 for bigBed
    bb->type = 1;

    bb->URL = urlOpen(fname, *callBack, NULL);
    if(!bb->URL) goto error;

    //Attempt to read in the fixed header
    bmHdrRead(bb);
    if(!bb->hdr) goto error;

    //Read in the chromosome list
    bb->cl = bmReadChromList(bb);
    if(!bb->cl) goto error;

    //Read in the index
    bb->idx = bmReadIndex(bb, 0);
    if(!bb->idx) goto error;

    return bb;

error:
    bmClose(bb);
    return NULL;
}

#ifndef LIBBIGWIG_H
#define LIBBIGWIG_H

#include "binaMethIO.h"
#include "bmValues.h"
#include <inttypes.h>
#include <zlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \mainpage libBinaMeth
 *
 * \section Introduction
 *
 * libBinaMeth is a C library for parsing local/remote binaMeth and bigBed files. This is similar to Kent's library from UCSC, except 
 *  * The license is much more liberal
 *  * This code doesn't call `exit()` on error, thereby killing the calling application.
 *
 * External files are accessed using [curl](http://curl.haxx.se/).
 *
 * Please submit issues and pull requests [here](https://github.com/dpryan79/libBinaMeth).
 *
 * \section Compilation
 *
 * Assuming you already have the curl libraries installed (not just the curl binary!):
 *
 *     make install prefix=/some/path
 *
 * \section Writing binaMeth files
 *
 * There are three methods for storing values in a binaMeth file, further described in the [wiggle format](http://genome.ucsc.edu/goldenpath/help/wiggle.html). The entries within the file are grouped into "blocks" and each such block is limited to storing entries of a single type. So, it is unwise to use a single bedGraph-like endtry followed by a single fixed-step entry followed by a variable-step entry, as that would require three separate blocks, with additional space required for each.
 *
 * \section Testing file types
 *
 * As of version 0.3.0, libBinaMeth supports reading bigBed files. If an application needs to support both bigBed and binaMeth input, then the `bmIsBinaMeth` and `bbIsBigBed` functions can be used to determine the file type. These both use the "magic" number at the beginning of the file to determine the file type.
 *
 * \section Interval and entry iterators
 *
 * As of version 0.3.0, libBinaMeth supports iterating over intervals in binaMeth files and entries in bigBed files. The number of intervals/entries returned with each iteration can be controlled by setting the number of blocks processed in each iteration (intervals and entries are group inside of binaMeth and bigBed files into blocks of entries). See `test/testIterator.c` for an example.
 *
 * \section Examples
 * 
 * Please see [README.md](README.md) and the files under `test/` for examples.
 */
 

/*! \file binaMeth.h
 *
 * These are the functions and structured that should be used by external users. While I don't particularly recommend dealing with some of the structures (e.g., a binaMethHdr_t), they're described here in case you need them.
 *
 * BTW, this library doesn't switch endianness as appropriate, since I kind of assume that there's only one type produced these days.
 */

/*!
 * The library version number
 */
#define LIBBIGWIG_VERSION 0.4.6

/*!
 * If 1, then this library was compiled with remote file support.
 */
#ifdef NOCURL
#define LIBBIGWIG_CURL 0
#ifndef CURLTYPE_DEFINED
#define CURLTYPE_DEFINED
typedef int CURLcode;
typedef void CURL;
#endif
#else
#define LIBBIGWIG_CURL 1
#endif

/*
* 0, without; 1, important; 2, all
*/
#define DEBUG 0
/*!
 * chrom 111
 */
#define BM_COVER 0x7
/*!
 * strand 111000
 */
#define BM_STRAND 0x38
/*!
 * context 111000000
 */
#define BM_CONTEXT 0x1c0
/*!
 * context 111000000000
 */
#define BM_ID 0xe00
/*!
 * context 111000000000000
 */
#define BM_END 0x7000
/*
* bigmeth file, 1000000000000000
*/
#define BM_MAGIC 0x8000

/*!
 * The magic number of a binaMeth file.
 */
#define BIGWIG_MAGIC 0x888FFC26
/*!
 * The magic number of a bigBed file.
 */
#define BIGBED_MAGIC 0x8789F2EB
/*!
 * The magic number of a "cirTree" block in a file.
 */
#define CIRTREE_MAGIC 0x78ca8c91
/*!
 * The magic number of an index block in a file.
 */
#define IDX_MAGIC 0x2468ace0
/*!
 * The default number of children per block.
 */
#define DEFAULT_nCHILDREN 64
/*!
 * The default decompression buffer size in bytes. This is used to determin
 */
#define DEFAULT_BLOCKSIZE 32768

/*!
 * An enum that dictates the type of statistic to fetch for a given interval
 */
enum bmStatsType {
    doesNotExist = -1, /*!< This does nothing */
    mean = 0, /*!< The mean value */
    average = 0, /*!< The mean value */
    weighted = 6,
    stdev = 1, /*!< The standard deviation of the values */
    dev = 1, /*!< The standard deviation of the values */
    max = 2, /*!< The maximum value */
    min = 3, /*!< The minimum value */
    cov = 4, /*!< The number of bases covered */
    coverage = 4, /*!<The number of bases covered */ 
    sum = 5 /*!< The sum of per-base values */
};

//Should hide this from end users
/*!
 * @brief BinaMeth files have multiple "zoom" levels, each of which has its own header. This hold those headers
 *
 * N.B., there's 4 bytes of padding in the on disk representation of level and dataOffset.
 */
typedef struct {
    uint32_t *level; /**<The zoom level, which is an integer starting with 0.*/
    //There's 4 bytes of padding between these
    uint64_t *dataOffset; /**<The offset to the on-disk start of the data. This isn't used currently.*/
    uint64_t *indexOffset; /**<The offset to the on-disk start of the index. This *is* used.*/
    bmRTree_t **idx; /**<Index for each zoom level. Represented as a tree*/
} bmZoomHdr_t;

/*!
 * @brief The header section of a binaMeth file.
 *
 * Some of the values aren't currently used for anything. Others may optionally not exist.
 */
typedef struct {
    uint16_t version; /**<The version information of the file.*/
    uint16_t nLevels; /**<The number of "zoom" levels.*/
    uint64_t ctOffset; /**<The offset to the on-disk chromosome tree list.*/
    uint64_t dataOffset; /**<The on-disk offset to the first block of data.*/
    uint64_t indexOffset; /**<The on-disk offset to the data index.*/
    uint16_t fieldCount; /**<Total number of fields.*/
    uint16_t definedFieldCount; /**<Number of fixed-format BED fields.*/
    uint64_t sqlOffset; /**<The on-disk offset to an SQL string. This is unused.*/
    uint64_t summaryOffset; /**<If there's a summary, this is the offset to it on the disk.*/
    uint32_t bufSize; /**<The compression buffer size (if the data is compressed).*/
    uint64_t extensionOffset; /**<Unused*/
    bmZoomHdr_t *zoomHdrs; /**<Pointers to the header for each zoom level.*/
    //total Summary
    uint64_t nBasesCovered; /**<The total bases covered in the file.*/
    double minVal; /**<The minimum value in the file.*/
    double maxVal; /**<The maximum value in the file.*/
    double sumData; /**<The sum of all values in the file.*/
    double sumSquared; /**<The sum of the squared values in the file.*/
//    uint32_t *pcoverage; /*The coverage of all postive strand cytosines.*/
//    uint32_t *ncoverage; /*The coverage of all negative strand cytosines.*/
} binaMethHdr_t;

//Should probably replace this with a hash
/*!
 * @brief Holds the chromosomes and their lengths
 */
typedef struct {
    int64_t nKeys; /**<The number of chromosomes */
    char **chrom; /**<A list of null terminated chromosomes */
    uint32_t *len; /**<The lengths of each chromosome */
//    uint32_t *NmethC; /*The number of C of each chromosome*/
//    uint32_t *NmethCover; /*The number of cover of each chromosome*/
} chromList_t;

//TODO remove from binaMeth.h
/// @cond SKIP
typedef struct bmLL bmLL;
struct bmLL {
    bmRTreeNode_t *node;
    struct bmLL *next;
};
typedef struct bmZoomBuffer_t bmZoomBuffer_t;
struct bmZoomBuffer_t { //each individual entry takes 32 bytes
    void *p;
    uint32_t l, m;
    struct bmZoomBuffer_t *next;
};
/// @endcond

/*!
 * @brief This is only needed for writing binaMeth files (and won't be created otherwise)
 * This should be removed from binaMeth.h
 */
typedef struct {
    uint64_t nBlocks; /**<The number of blocks written*/
    uint32_t blockSize; /**<The maximum number of children*/
    uint64_t nEntries; /**<The number of entries processed. This is used for the first contig and determining how the zoom levels are computed*/
    uint64_t runningWidthSum; /**<The running sum of the entry widths for the first contig (again, used for the first contig and computing zoom levels)*/
    uint32_t tid; /**<The current TID that's being processed*/
    uint32_t start; /**<The start position of the block*/
    uint32_t end; /**<The end position of the block*/
    uint32_t span; /**<The span of each entry, if applicable*/
    uint32_t step; /**<The step size, if applicable*/
    uint8_t ltype; /**<The type of the last entry added*/
    uint32_t l; /**<The current size of p. This and the type determine the number of items held*/
    uint16_t nItems;
    void *p; /**<A buffer of size hdr->bufSize*/
    bmLL *firstIndexNode; /**<The first index node in the linked list*/
    bmLL *currentIndexNode; /**<The last index node in a linked list*/
    bmZoomBuffer_t **firstZoomBuffer; /**<The first node in a linked list of leaf nodes*/
    bmZoomBuffer_t **lastZoomBuffer; /**<The last node in a linked list of leaf nodes*/
    uint64_t *nNodes; /**<The number of leaf nodes per zoom level, useful for determining duplicate levels*/
    uLongf compressPsz; /**<The size of the compression buffer*/
    void *compressP; /**<A compressed buffer of size compressPsz*/
} bmWriteBuffer_t;

/*!
 * @brief A structure that holds everything needed to access a binaMeth file.
 */
typedef struct {
    URL_t *URL; /**<A pointer that can handle both local and remote files (including a buffer if needed).*/
    binaMethHdr_t *hdr; /**<The file header.*/
    chromList_t *cl; /**<A list of chromosome names (the order is the ID).*/
    bmRTree_t *idx; /**<The index for the full dataset.*/
    bmWriteBuffer_t *writeBuffer; /**<The buffer used for writing.*/
    int isWrite; /**<0: Opened for reading, 1: Opened for writing.*/
    int type; /**<0: binaMeth, 1: bigBed.*/
    int mtype;
    /* 0: bigmeth, 1: bigmeth for region with ID, 2: bigmeth for region without ID, 3: ID with strand */
    /**<0: binaMeth, 1: bigBed.*/
} binaMethFile_t;

/*!
 * @brief Holds interval:value associations
 * 0. chrom start end value coverage strand context
 * 4    4   4   4   2   1   1
 * 12 + 4 = 16
 * type 1. chrom start end value coverage geneID
 * type 2. chrom start end value coverage 
 */
typedef struct {
    uint32_t l; /**<Number of intervals held*/
    uint32_t m; /**<Maximum number of values/intervals the struct can hold*/

    uint32_t *start; /**<The start positions (0-based half open)*/
    uint32_t *end; /**<The end positions (0-based half open)*/
    float *value; /**<The value associated with each position*/

    uint16_t *coverage;
    uint8_t *strand;
    uint8_t *context;
    char **entryid;
} bmOverlappingIntervals_t;

/*!
 * @brief Holds interval:str associations
 */
typedef struct {
    uint32_t l; /**<Number of intervals held*/
    uint32_t m; /**<Maximum number of values/intervals the struct can hold*/
    uint32_t *start; /**<The start positions (0-based half open)*/
    uint32_t *end; /**<The end positions (0-based half open)*/
    char **str; /**<The strings associated with a given entry.*/
} bbOverlappingEntries_t;

/*!
 * @brief A structure to hold iterations
 * One of intervals and entries should be used to access records from binaMeth or bigBed files, respectively.
 */
typedef struct {
    binaMethFile_t *bm; /**<Pointer to the binaMeth/bigBed file.*/
    uint32_t tid; /**<The contig/chromosome ID.*/
    uint32_t start; /**<Start position of the query interval.*/
    uint32_t end; /**<End position of the query interval.*/
    uint64_t offset; /**<Offset into the blocks.*/
    uint32_t blocksPerIteration; /**<Number of blocks to use per iteration.*/
    int withString; /**<For bigBed entries, whether to return the string with the entries.*/
    void *blocks; /**<Overlapping blocks.*/
    bmOverlappingIntervals_t *intervals; /**<Overlapping intervals (or NULL).*/
    bbOverlappingEntries_t *entries; /**<Overlapping entries (or NULL).*/
    void *data; /**<Points to either intervals or entries. If there are no further intervals/entries, then this is NULL. Use this to test for whether to continue iterating.*/
} bmOverlapIterator_t;

/*!
 * @brief Initializes curl and global variables. This *MUST* be called before other functions (at least if you want to connect to remote files).
 * For remote file, curl must be initialized and regions of a file read into an internal buffer. If the buffer is too small then an excessive number of connections will be made. If the buffer is too large than more data than required is fetched. 128KiB is likely sufficient for most needs.
 * @param bufSize The internal buffer size used for remote connection.
 * @see bmCleanup
 * @return 0 on success and 1 on error.
 */
int bmInit(size_t bufSize);

/*!
 * @brief The counterpart to bmInit, this cleans up curl.
 * @see bmInit
 */
void bmCleanup(void);

/*!
 * @brief Determine if a file is a binaMeth file.
 * This function will quickly check either local or remote files to determine if they appear to be valid binaMeth files. This can be determined by reading the first 4 bytes of the file.
 * @param fname The file name or URL (http, https, and ftp are supported)
 * @param callBack An optional user-supplied function. This is applied to remote connections so users can specify things like proxy and password information. See `test/testRemote` for an example.
 * @return 1 if the file appears to be binaMeth, otherwise 0.
 */
int bmIsBinaMeth(char *fname, CURLcode (*callBack)(CURL*));

/*!
 * @brief Determine is a file is a bigBed file.
 * This function will quickly check either local or remote files to determine if they appear to be valid binaMeth files. This can be determined by reading the first 4 bytes of the file.
 * @param fname The file name or URL (http, https, and ftp are supported)
 * @param callBack An optional user-supplied function. This is applied to remote connections so users can specify things like proxy and password information. See `test/testRemote` for an example.
 * @return 1 if the file appears to be binaMeth, otherwise 0.
 */
int bbIsBigBed(char *fname, CURLcode (*callBack)(CURL*));

/*!
 * @brief Opens a local or remote binaMeth file.
 * This will open a local or remote binaMeth file. Writing of local binaMeth files is also supported.
 * @param fname The file name or URL (http, https, and ftp are supported)
 * @param callBack An optional user-supplied function. This is applied to remote connections so users can specify things like proxy and password information. See `test/testRemote` for an example.
 * @param mode The mode, by default "r". Both local and remote files can be read, but only local files can be written. For files being written the callback function is ignored. If and only if the mode contains "w" will the file be opened for writing (in all other cases the file will be opened for reading.
 * @return A binaMethFile_t * on success and NULL on error.
 */
binaMethFile_t *bmOpen(char *fname, CURLcode (*callBack)(CURL*), const char* mode);

/*!
 * @brief Opens a local or remote bigBed file.
 * This will open a local or remote bigBed file. Note that this file format can only be read and NOT written!
 * @param fname The file name or URL (http, https, and ftp are supported)
 * @param callBack An optional user-supplied function. This is applied to remote connections so users can specify things like proxy and password information. See `test/testRemote` for an example.
 * @return A binaMethFile_t * on success and NULL on error.
 */
binaMethFile_t *bbOpen(char *fname, CURLcode (*callBack)(CURL*));

/*!
 * @brief Returns a string containing the SQL entry (or NULL).
 * The "auto SQL" field contains the names and value types of the entries in
 * each bigBed entry. If you need to parse a particular value out of each entry,
 * then you'll need to first parse this.
 * @param fp The file pointer to a valid binaMethFile_t
 * @return A char *, which you MUST free!
 */
char *bbGetSQL(binaMethFile_t *fp);

/*!
 * @brief Closes a binaMethFile_t and frees up allocated memory
 * This closes both binaMeth and bigBed files.
 * @param fp The file pointer.
 */
void bmClose(binaMethFile_t *fp);

/*******************************************************************************
*
* The following are in bmStats.c
*
*******************************************************************************/

/*!
 * @brief Converts between chromosome name and ID
 *
 * @param fp A valid binaMethFile_t pointer
 * @param chrom A chromosome name
 * @return An ID, -1 will be returned on error (note that this is an unsigned value, so that's ~4 billion. binaMeth/bigBed files can't store that many chromosomes anyway.
 */
uint32_t bmGetTid(binaMethFile_t *fp, char *chrom);

/*!
 * @brief Frees space allocated by `bmGetOverlappingIntervals`
 * @param o A valid `bmOverlappingIntervals_t` pointer.
 * @see bmGetOverlappingIntervals
 */
void bmDestroyOverlappingIntervals(bmOverlappingIntervals_t *o);

/*!
 * @brief Frees space allocated by `bbGetOverlappingEntries`
 * @param o A valid `bbOverlappingEntries_t` pointer.
 * @see bbGetOverlappingEntries
 */
void bbDestroyOverlappingEntries(bbOverlappingEntries_t *o);

/*!
 * @brief Return binaMeth entries overlapping an interval.
 * Find all binaMeth entries overlapping a range and returns them, including their associated values.
 * @param fp A valid binaMethFile_t pointer. This MUST be for a binaMeth file!
 * @param chrom A valid chromosome name.
 * @param start The start position of the interval. This is 0-based half open, so 0 is the first base.
 * @param end The end position of the interval. Again, this is 0-based half open, so 100 will include the 100th base...which is at position 99.
 * @return NULL on error or no overlapping values, otherwise a `bmOverlappingIntervals_t *` holding the values and intervals.
 * @see bmOverlappingIntervals_t
 * @see bmDestroyOverlappingIntervals
 * @see bmGetValues
 */
bmOverlappingIntervals_t *bmGetOverlappingIntervals(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end);

/*!
 * @brief Return bigBed entries overlapping an interval.
 * Find all bigBed entries overlapping a range and returns them.
 * @param fp A valid binaMethFile_t pointer. This MUST be for a bigBed file!
 * @param chrom A valid chromosome name.
 * @param start The start position of the interval. This is 0-based half open, so 0 is the first base.
 * @param end The end position of the interval. Again, this is 0-based half open, so 100 will include the 100th base...which is at position 99.
 * @param withString If not 0, return the string associated with each entry in the output. If 0, there are no associated strings returned. This is useful if the only information needed are the locations of the entries, which require significantly less memory.
 * @return NULL on error or no overlapping values, otherwise a `bbOverlappingEntries_t *` holding the intervals and (optionally) the associated string.
 * @see bbOverlappingEntries_t
 * @see bbDestroyOverlappingEntries
 */
bbOverlappingEntries_t *bbGetOverlappingEntries(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, int withString);

/*!
 * @brief Creates an iterator over intervals in a binaMeth file
 * Iterators can be traversed with `bmIteratorNext()` and destroyed with `bmIteratorDestroy()`.
 * Intervals are in the `intervals` member and `data` can be used to determine when to end iteration.
 * @param fp A valid binaMethFile_t pointer. This MUST be for a binaMeth file!
 * @param chrom A valid chromosome name.
 * @param start The start position of the interval. This is 0-based half open, so 0 is the first base.
 * @param end The end position of the interval. Again, this is 0-based half open, so 100 will include the 100th base...which is at position 99.
 * @param blocksPerIteration The number of blocks (internal groupings of intervals in binaMeth files) to return per iteration.
 * @return NULL on error, otherwise a bmOverlapIterator_t pointer
 * @see bmOverlapIterator_t
 * @see bmIteratorNext
 * @see bmIteratorDestroy
 */ 
bmOverlapIterator_t *bmOverlappingIntervalsIterator(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t blocksPerIteration);

/*!
 * @brief Creates an iterator over entries in a bigBed file
 * Iterators can be traversed with `bmIteratorNext()` and destroyed with `bmIteratorDestroy()`.
 * Entries are in the `entries` member and `data` can be used to determine when to end iteration.
 * @param fp A valid binaMethFile_t pointer. This MUST be for a bigBed file!
 * @param chrom A valid chromosome name.
 * @param start The start position of the interval. This is 0-based half open, so 0 is the first base.
 * @param end The end position of the interval. Again, this is 0-based half open, so 100 will include the 100th base...which is at position 99.
 * @param withString Whether the returned entries should include their associated strings.
 * @param blocksPerIteration The number of blocks (internal groupings of entries in bigBed files) to return per iteration.
 * @return NULL on error, otherwise a bmOverlapIterator_t pointer
 * @see bbGetOverlappingEntries
 * @see bmOverlapIterator_t
 * @see bmIteratorNext
 * @see bmIteratorDestroy
 */ 
bmOverlapIterator_t *bbOverlappingEntriesIterator(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, int withString, uint32_t blocksPerIteration);

/*!
 * @brief Traverses to the entries/intervals in the next group of blocks.
 * @param iter A bmOverlapIterator_t pointer that is updated (or destroyed on error)
 * @return NULL on error, otherwise a bmOverlapIterator_t pointer with the intervals or entries from the next set of blocks.
 * @see bmOverlapIterator_t
 * @see bmIteratorDestroy
 */ 
bmOverlapIterator_t *bmIteratorNext(bmOverlapIterator_t *iter);

/*!
 * @brief Destroys a bmOverlapIterator_t
 * @param iter The bmOverlapIterator_t that should be destroyed
 */
void bmIteratorDestroy(bmOverlapIterator_t *iter);

/*!
 * @brief Return all per-base binaMeth values in a given interval.
 * Given an interval (e.g., chr1:0-100), return the value at each position in a binaMeth file. Positions without associated values are suppressed by default, but may be returned if `includeNA` is not 0.
 * @param fp A valid binaMethFile_t pointer.
 * @param chrom A valid chromosome name.
 * @param start The start position of the interval. This is 0-based half open, so 0 is the first base.
 * @param end The end position of the interval. Again, this is 0-based half open, so 100 will include the 100th base...which is at position 99.
 * @param includeNA If not 0, report NA values as well (as NA).
 * @return NULL on error or no overlapping values, otherwise a `bmOverlappingIntervals_t *` holding the values and positions.
 * @see bmOverlappingIntervals_t
 * @see bmDestroyOverlappingIntervals
 * @see bmGetOverlappingIntervals
 */
bmOverlappingIntervals_t *bmGetValues(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, int includeNA);

/*!
 * @brief Determines per-interval binaMeth statistics
 * Can determine mean/min/max/coverage/standard deviation of values in one or more intervals in a binaMeth file. You can optionally give it an interval and ask for values from X number of sub-intervals.
 * @param fp The file from which to extract statistics.
 * @param chrom A valid chromosome name.
 * @param start The start position of the interval. This is 0-based half open, so 0 is the first base.
 * @param end The end position of the interval. Again, this is 0-based half open, so 100 will include the 100th base...which is at position 99.
 * @param nBins The number of bins within the interval to calculate statistics for.
 * @param type The type of statistic.
 * @see bmStatsType
 * @return A pointer to an array of double precission floating point values. Note that binaMeth files only hold 32-bit values, so this is done to help prevent overflows.
 */
double *bmStats(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand, uint8_t context);
double *bmStats_array(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand);
void bmStats_array_count(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand, uint16_t *countC, uint16_t *countCT);

/*!
 * @brief Determines per-interval binaMeth statistics
 * Can determine mean/min/max/coverage/standard deviation of values in one or more intervals in a binaMeth file. You can optionally give it an interval and ask for values from X number of sub-intervals. The difference with bmStats is that zoom levels are never used.
 * @param fp The file from which to extract statistics.
 * @param chrom A valid chromosome name.
 * @param start The start position of the interval. This is 0-based half open, so 0 is the first base.
 * @param end The end position of the interval. Again, this is 0-based half open, so 100 will include the 100th base...which is at position 99.
 * @param nBins The number of bins within the interval to calculate statistics for.
 * @param type The type of statistic.
 * @see bmStatsType
 * @return A pointer to an array of double precission floating point values. Note that binaMeth files only hold 32-bit values, so this is done to help prevent overflows.
*/
double *bmStatsFromFull(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand, uint8_t context);
double *bmStatsFromFull_array(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, uint32_t movestep, enum bmStatsType type, uint8_t strand);

//Writer functions

/*!
 * @brief Create a largely empty binaMeth header
 * Every binaMeth file has a header, this creates the template for one. It also takes care of space allocation in the output write buffer.
 * @param fp The binaMethFile_t* that you want to write to.
 * @param maxZooms The maximum number of zoom levels. If you specify 0 then there will be no zoom levels. A value <0 or > 65535 will result in a maximum of 10.
 * @return 0 on success.
 */
int bmCreateHdr(binaMethFile_t *fp, int32_t maxZooms);

/*!
 * @brief Take a list of chromosome names and lengths and return a pointer to a chromList_t
 * This MUST be run before `bmWriteHdr()`. Note that the input is NOT free()d!
 * @param chroms A list of chromosomes.
 * @param lengths The length of each chromosome.
 * @param n The number of chromosomes (thus, the length of `chroms` and `lengths`)
 * @return A pointer to a chromList_t or NULL on error.
 */
chromList_t *bmCreateChromList(char **chroms, uint32_t *lengths, int64_t n);
chromList_t *bmCreateChromList_ifp(binaMethFile_t *ibm);
/*!
 * @brief Write a the header to a binaMeth file.
 * You must have already opened the output file, created a header and a chromosome list.
 * @param bm The output binaMethFile_t pointer.
 * @see bmCreateHdr
 * @see bmCreateChromList
 */
int bmWriteHdr(binaMethFile_t *bm);

/*!
 * @brief Write a new block of bedGraph-like intervals to a binaMeth file
 * Adds entries of the form:
 * chromosome	start	end	value
 * to the file. These will always be added in a new block, so you may have previously used a different storage type.
 * 
 * In general it's more efficient to use the bmAppend* functions, but then you MUST know that the previously written block is of the same type. In other words, you can only use bmAppendIntervals() after bmAddIntervals() or a previous bmAppendIntervals().
 * @param fp The output file pointer.
 * @param chrom A list of chromosomes, of length `n`.
 * @param start A list of start positions of length`n`.
 * @param end A list of end positions of length`n`.
 * @param values A list of values of length`n`.
 * @param n The length of the aforementioned lists.
 * @return 0 on success and another value on error.
 * @see bmAppendIntervals
 */
int bmAddIntervals(binaMethFile_t *fp, char **chrom, uint32_t *start, uint32_t *end, float *values, uint16_t *coverage, uint8_t *strand,
    uint8_t *context, char **entryid, uint32_t n);

/*!
 * @brief Append bedGraph-like intervals to a previous block of bedGraph-like intervals in a binaMeth file.
 * If you have previously used bmAddIntervals() then this will append additional entries into the previous block (or start a new one if needed).
 * @param fp The output file pointer.
 * @param start A list of start positions of length`n`.
 * @param end A list of end positions of length`n`.
 * @param values A list of values of length`n`.
 * @param n The length of the aforementioned lists.
 * @return 0 on success and another value on error.
 * @warning Do NOT use this after `bmAddIntervalSpanSteps()`, `bmAppendIntervalSpanSteps()`, `bmAddIntervalSpanSteps()`, or `bmAppendIntervalSpanSteps()`.
 * @see bmAddIntervals
 */
int bmAppendIntervals(binaMethFile_t *fp, uint32_t *start, uint32_t *end, float *values, uint16_t *coverage, uint8_t *strand,
    uint8_t *context, char **entryid, uint32_t n);

/*!
 * @brief Add a new block of variable-step entries to a binaMeth file
 * Adds entries for the form
 * chromosome	start	value
 * to the file. Each block of such entries has an associated "span", so each value describes the region chromosome:start-(start+span)
 *
 * This will always start a new block of values.
 * @param fp The output file pointer.
 * @param chrom A list of chromosomes, of length `n`.
 * @param start A list of start positions of length`n`.
 * @param span The span of each entry (the must all be the same).
 * @param values A list of values of length`n`.
 * @param n The length of the aforementioned lists.
 * @return 0 on success and another value on error.
 * @see bmAppendIntervalSpans
 */
int bmAddIntervalSpans(binaMethFile_t *fp, char *chrom, uint32_t *start, uint32_t span, float *values, uint32_t n);

/*!
 * @brief Append to a previous block of variable-step entries.
 * If you previously used `bmAddIntervalSpans()`, this will continue appending more values to the block(s) it created.
 * @param fp The output file pointer.
 * @param start A list of start positions of length`n`.
 * @param values A list of values of length`n`.
 * @param n The length of the aforementioned lists.
 * @return 0 on success and another value on error.
 * @warning Do NOT use this after `bmAddIntervals()`, `bmAppendIntervals()`, `bmAddIntervalSpanSteps()` or `bmAppendIntervalSpanSteps()`
 * @see bmAddIntervalSpans
 */
int bmAppendIntervalSpans(binaMethFile_t *fp, uint32_t *start, float *values, uint32_t n);

/*!
 * @brief Add a new block of fixed-step entries to a binaMeth file
 * Adds entries for the form
 * value
 * to the file. Each block of such entries has an associated "span", "step", chromosome and start position. See the wiggle format for more details.
 *
 * This will always start a new block of values.
 * @param fp The output file pointer.
 * @param chrom The chromosome that the entries describe.
 * @param start The starting position of the block of entries.
 * @param span The span of each entry (i.e., the number of bases it describes).
 * @param step The step between entry start positions.
 * @param values A list of values of length`n`.
 * @param n The length of the aforementioned lists.
 * @return 0 on success and another value on error.
 * @see bmAddIntervalSpanSteps
 */
int bmAddIntervalSpanSteps(binaMethFile_t *fp, char *chrom, uint32_t start, uint32_t span, uint32_t step, float *values, uint32_t n);

/*!
 * @brief Append to a previous block of fixed-step entries.
 * If you previously used `bmAddIntervalSpanSteps()`, this will continue appending more values to the block(s) it created.
 * @param fp The output file pointer.
 * @param values A list of values of length`n`.
 * @param n The length of the aforementioned lists.
 * @return 0 on success and another value on error.
 * @warning Do NOT use this after `bmAddIntervals()`, `bmAppendIntervals()`, `bmAddIntervalSpans()` or `bmAppendIntervalSpans()`
 * @see bmAddIntervalSpanSteps
 */
int bmAppendIntervalSpanSteps(binaMethFile_t *fp, float *values, uint32_t n);

#ifdef __cplusplus
}
#endif

#endif // LIBBIGWIG_H

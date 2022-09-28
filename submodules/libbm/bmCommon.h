/*! \file bmCommon.h
 *
 * You have no reason to use these functions. They may change without warning because there's no reason for them to be used outside of libBigWig's internals.
 *
 * These are structures and functions from a variety of files that are used across files internally but don't need to be see by libBigWig users.
 */

/*!
 * @brief Like fsetpos, but for local or remote binaMeth files.
 * This will set the file position indicator to the specified point. For local files this literally is `fsetpos`, while for remote files it fills a memory buffer with data starting at the desired position.
 * @param fp A valid opened binaMethFile_t.
 * @param pos The position within the file to seek to.
 * @return 0 on success and -1 on error.
 */
int bmSetPos(binaMethFile_t *fp, size_t pos);

/*!
 * @brief A local/remote version of `fread`.
 * Reads data from either local or remote binaMeth files.
 * @param data An allocated memory block big enough to hold the data.
 * @param sz The size of each member that should be copied.
 * @param nmemb The number of members to copy.
 * @param fp The binaMethFile_t * from which to copy the data.
 * @see bmSetPos
 * @return For nmemb==1, the size of the copied data. For nmemb>1, the number of members fully copied (this is equivalent to `fread`).
 */
size_t bmRead(void *data, size_t sz, size_t nmemb, binaMethFile_t *fp);

/*!
 * @brief Determine what the file position indicator say.
 * This is equivalent to `ftell` for local or remote files.
 * @param fp The file.
 * @return The position in the file.
 */
long bmTell(binaMethFile_t *fp);

/*!
 * @brief Reads a data index (either full data or a zoom level) from a binaMeth file.
 * There is little reason for end users to use this function. This must be freed with `bmDestroyIndex`
 * @param fp A valid binaMethFile_t pointer
 * @param offset The file offset where the index begins
 * @return A bmRTree_t pointer or NULL on error.
 */
bmRTree_t *bmReadIndex(binaMethFile_t *fp, uint64_t offset);

/*!
 * @brief Destroy an bmRTreeNode_t and all of its children.
 * @param node The node to destroy.
 */
void bmDestroyIndexNode(bmRTreeNode_t *node);

/*!
 * @brief Frees space allocated by `bmReadIndex`
 * There is generally little reason to use this, since end users should typically not need to run `bmReadIndex` themselves.
 * @param idx A bmRTree_t pointer allocated by `bmReadIndex`.
 */
void bmDestroyIndex(bmRTree_t *idx);

/// @cond SKIP
bmOverlapBlock_t *walkRTreeNodes(binaMethFile_t *bm, bmRTreeNode_t *root, uint32_t tid, uint32_t start, uint32_t end);
void destroyBWOverlapBlock(bmOverlapBlock_t *b);
/// @endcond

/*!
 * @brief Finishes what's needed to write a binaMethFile
 * Flushes the buffer, converts the index linked list to a tree, writes that to disk, handles zoom level stuff, writes magic at the end
 * @param fp A valid binaMethFile_t pointer
 * @return 0 on success
 */
int bmFinalize(binaMethFile_t *fp);

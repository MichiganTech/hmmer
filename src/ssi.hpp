/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

#pragma once

/* ssi.h
 * Database indexing (SSI format support)
 *
 * See: ssi_format.tex in Docs/
 */

#include <stdio.h>
#include <stdint.h>

/* Limits
 */
#define SSI_MAXFILES 32767	  /* 2^15-1 */
#define SSI_MAXKEYS  2147483647L  /* 2^31-1 */
#define SSI_MAXRAM   200	  /* allow 200MB indexes before external sort mode */

/* typedef: SSIOFFSET
 * Use the union to save space, since the two offset types are
 * mutually exclusive, controlled by "mode"
 */
struct ssioffset_s {
  char mode;			/* GSI_OFFSET_I32, for example */
  union {
    uint32_t   i32;           /* an offset that fseek() can use         */
    uint64_t   i64;           /* an offset that e.g. fseeko64() can use */
  } off;
};
typedef struct ssioffset_s SSIOFFSET;
#define SSI_OFFSET_I32    0
#define SSI_OFFSET_I64    1


/* Structure: SSIFILE
 * xref:      SSI API documentation in ssi-format.tex
 */
struct ssifile_s {
  FILE        *fp;		/* open SSI index file                 */
  uint32_t   flags;		/* optional behavior flags             */
  uint16_t   nfiles;		/* number of files = 16 bit int        */
  uint32_t   nprimary;	/* number of primary keys              */
  uint32_t   nsecondary;	/* number of secondary keys            */
  uint32_t   flen;		/* length of filenames (inc '\0')      */
  uint32_t   plen;		/* length of primary keys (inc '\0')   */
  uint32_t   slen;		/* length of secondary keys (inc '\0') */
  uint32_t   frecsize;	/* # bytes in a file record            */
  uint32_t   precsize;	/* # bytes in a primary key record     */
  uint32_t   srecsize;	/* # bytes in a secondary key record   */
  SSIOFFSET    foffset;		/* disk offset, start of file records  */
  SSIOFFSET    poffset;		/* disk offset, start of pri key recs  */
  SSIOFFSET    soffset;		/* disk offset, start of sec key recs  */
 
  char imode;			/* mode for index file offsets, 32 v. 64 bit    */
  char smode;			/* mode for sequence file offsets, 32 v. 64 bit */

  /* File information:
   */
  char       **filename;	/* list of file names [0..nfiles-1]    */
  uint32_t  *fileformat;	/* file formats                        */
  uint32_t  *fileflags;       /* optional per-file behavior flags    */
  uint32_t  *bpl;     	/* bytes per line in file              */
  uint32_t  *rpl;     	/* residues per line in file           */
};
typedef struct ssifile_s SSIFILE;

/* optional per-index behavior flags in SSIFILE structure's flags:
 */
#define SSI_USE64        1<<0	/* seq offsets are 64-bit        */
#define SSI_USE64_INDEX  1<<1	/* index file offsets are 64-bit */

/* optional per-file behavior flags in fileflags
 */
#define SSI_FAST_SUBSEQ  1<<0	/* can do subseq lookup in this file */

/* Structure: SSIINDEX
 *
 * Used when building up an index and writing it to disk
 */
struct ssipkey_s {		/* Primary key data: */
  char        *key;             /* key name          */
  uint16_t   fnum;		/* file number       */
  SSIOFFSET    r_off;		/* record offset     */
  SSIOFFSET    d_off;		/* data offset       */
  uint32_t   len;		/* sequence length   */
};


struct ssiskey_s {		/* Secondary key data: */
  char        *key;             /* secondary key name  */
  char        *pkey;            /* primary key name    */
};


struct ssiindex_s {
  int           smode;		/* sequence mode: SSI_OFFSET_I32 or _I64 */
  int           imode;		/* index mode:    SSI_OFFSET_I32 or _I64 */
  int           external;	/* true if pkeys and skeys are on disk   */
  int           max_ram;	/* maximum RAM in MB before switching to external */

  char        **filenames;
  uint32_t   *fileformat;
  uint32_t   *bpl;
  uint32_t   *rpl;
  uint32_t    flen;		/* length of longest filename, inc '\0' */
  uint16_t    nfiles;
 
  struct ssipkey_s *pkeys;
  uint32_t         plen;	/* length of longest pkey, including '\0' */
  uint32_t         nprimary;
  char              *ptmpfile;	/* name of tmp file, for external sort mode */
  FILE              *ptmp;	/* handle on open ptmpfile */

  struct ssiskey_s *skeys;
  uint32_t         slen;	/* length of longest skey, including '\0' */
  uint32_t         nsecondary;
  char              *stmpfile;	/* name of tmp file, for external sort mode */
  FILE              *stmp;	/* handle on open ptmpfile */
};
typedef struct ssiindex_s SSIINDEX;

/* These control malloc and realloc chunk sizes in the index
 * construction code.
 */
#define SSI_FILE_BLOCK    10
#define SSI_KEY_BLOCK     100

/* Error codes set by the API
 */
#define SSI_ERR_OK           0
#define SSI_ERR_NODATA       1	/* no data? an fread() failed */
#define SSI_ERR_NO_SUCH_KEY  2 /* that key's not in the index */
#define SSI_ERR_MALLOC       3
#define SSI_ERR_NOFILE       4	/* no such file? an fopen() failed */
#define SSI_ERR_BADMAGIC     5	/* magic number mismatch in GSIOpen() */
#define SSI_ERR_BADFORMAT    6	/* didn't read what I expected to fread() */
#define SSI_ERR_NO64BIT      7	/* needed 64-bit support and didn't have it */
#define SSI_ERR_SEEK_FAILED  8 /* an fseek() (or similar) failed */
#define SSI_ERR_TELL_FAILED  9 /* an ftell() (or similar) failed */
#define SSI_ERR_NO_SUBSEQS   10 /* fast subseq is disallowed */
#define SSI_ERR_RANGE        11 /* subseq requested is out of range */
#define SSI_ERR_BADARG       12 /* something wrong with a function argument */
#define SSI_ERR_TOOMANY_FILES 13 /* ran out of range for files in an index */
#define SSI_ERR_TOOMANY_KEYS  14 /* ran out of range for keys in an index */
#define SSI_ERR_FWRITE        15
#define SSI_ERR_EXTERNAL_SORT 16 /* external sort failed */

/* The SSI file reading API:
 */


/* Function: SSIOpen()
 *
 * Purpose:  Opens the SSI index file {filename} and returns
 *           a SSIFILE * stream thru {ret_sfp}.
 *           The caller must eventually close this stream using
 *           SSIClose(). More than one index file can be open
 *           at once.
 *
 * Args:     filename - full path to a SSI index file
 *
 * Returns:  Returns 0 on success, nonzero on failure.
 */
int
SSIOpen(
  char *filename,
  SSIFILE **ret_sfp
);


/* load_indexfile(): given a SSIFILE structure with an open and positioned
 *    stream (fp) -- but no other data loaded -- read the next SSIFILE
 *    in from disk. We use this routine without its SSIOpen() wrapper
 *    as part of the external mergesort when creating large indices.
 */
int
load_indexfile(
  SSIFILE *sfp
);


/* Function: SSIGetOffsetByName()
 *
 * Purpose:  Looks up the string {key} in the open index {sfp}.
 *           {key} can be either a primary or secondary key. If {key}
 *           is found, {*ret_fh} contains a unique handle on
 *           the file that contains {key} (suitable for an SSIFileInfo()
 *           call, or for comparison to the handle of the last file
 *           that was opened for retrieval), and {offset} is filled
 *           in with the offset in that file.
 *          
 * Args:     sfp         - open index file
 *           key         - string to search for
 *           ret_fh      - RETURN: handle on file that key is in
 *           ret_offset  - RETURN: offset of the start of that key's record
 *
 * Returns:  0 on success.
 *           non-zero on error.
 */
int
SSIGetOffsetByName(
  SSIFILE *sfp,
  char *key,
  int *ret_fh,
  SSIOFFSET *ret_offset
);


/* Function: SSIGetOffsetByNumber()
 *
 * Purpose:  Looks up primary key #{n} in the open index {sfp}.
 *           {n} ranges from 0..nprimary-1. When key #{n}
 *           is found, {*ret_fh} contains a unique
 *           handle on the file that contains {key} (suitable
 *           for an SSIFileInfo() call, or for comparison to
 *           the handle of the last file that was opened for retrieval),
 *           and {offset} is filled in with the offset in that file.
 *          
 * Args:     sfp        - open index file
 *           n          - primary key number to retrieve.
 *           ret_fh     - RETURN: handle on file that key is in
 *           ret_offset - RETURN: offset of the start of that key's record
 *
 * Returns:  0 on success.
 *           non-zero on error.
 */
int
SSIGetOffsetByNumber(
  SSIFILE *sfp,
  int n,
  int *ret_fh,
  SSIOFFSET *ret_offset
);


/* Function: SSIGetSubseqOffset()
 *
 * Purpose:  Implements SSI_FAST_SUBSEQ.
 *
 *           Looks up a primary or secondary {key} in the open
 *           index {sfp}. Asks for the nearest offset to a
 *           subsequence starting at position {requested_start}
 *           in the sequence (numbering the sequence 1..L).
 *           If {key} is found, on return, {ret_fh}
 *           contains a unique handle on the file that contains
 *           {key} (suitable for an SSIFileInfo() call, or for
 *           comparison to the handle of the last file that was
 *           opened for retrieval); {record_offset} contains the
 *           disk offset to the start of the record; {data_offset}
 *           contains the disk offset either exactly at the requested
 *           residue, or at the start of the line containing the
 *           requested residue; {ret_actual_start} contains the
 *           coordinate (1..L) of the first valid residue at or
 *           after {data_offset}. {ret_actual_start} is <=
 *           {requested_start}.
 *
 * Args:     sfp             - open index file
 *           key             - primary or secondary key to find
 *           requested_start - residue we'd like to start at (1..L)
 *           ret_fh          - RETURN: handle for file the key is in
 *           record_offset   - RETURN: offset of entire record
 *           data_offset     - RETURN: offset of subseq (see above)
 *           ret_actual_start- RETURN: coord (1..L) of residue at data_offset
 *
 * Returns:  0 on success, non-zero on failure.
 */
int
SSIGetSubseqOffset(
  SSIFILE *sfp,
  char *key,
  int requested_start,
  int *ret_fh,
  SSIOFFSET *record_offset,
  SSIOFFSET *data_offset,
  int *ret_actual_start
);


/* Function: SSISetFilePosition()
 *
 * Purpose:  Uses {offset} to sets the file position for {fp}, usually an
 *           open sequence file, relative to the start of the file.
 *           Hides the details of system-dependent shenanigans necessary for
 *           file positioning in large (>2 GB) files.
 *          
 *           Behaves just like fseek(fp, offset, SEEK_SET) for 32 bit
 *           offsets and <2 GB files.
 *          
 *           Warning: if all else fails, in desperation, it will try to
 *           use fsetpos(). This requires making assumptions about fpos_t
 *           that may be unwarranted... assumptions that ANSI C prohibits
 *           me from making... though I believe the ./configure
 *           script robustly tests whether I can play with fpos_t like this.
 *
 * Args:     fp      - file to position.
 *           offset  - SSI offset relative to file start.
 *                
 * Returns:  0 on success, nonzero on error.
 */
int
SSISetFilePosition(
  FILE *fp, SSIOFFSET *offset
);


/* Function: SSIFileInfo()
 *
 * Purpose:  Given a file number {fh} in an open index file
 *           {sfp}, retrieve file name {ret_filename} and
 *           the file format {ret_format}.
 *          
 *           {ret_filename} is a pointer to a string maintained
 *           internally by {sfp}. It should not be free'd;
 *           SSIClose(sfp) takes care of it.
 *
 * Args:     sfp          - open index file
 *           fh           - handle on file to look up
 *           ret_filename - RETURN: name of file n
 *           ret_format   - RETURN: format of file n
 *
 * Returns:  0 on success, nonzero on failure.
 */
int
SSIFileInfo(
  SSIFILE *sfp,
  int fh,
  char **ret_filename,
  int *ret_format
);


/* Function: SSIClose()
 *
 * Purpose:  Close an open {SSIFILE *}.
 *
 * Args:     sfp - index file to close.
 *
 * Returns:  (void)
 */
void
SSIClose(
  SSIFILE *sfp
);


/* clear_ssifile(): free the innards of SSIFILE, without
 * destroying the structure or closing the stream.
 */
void
clear_ssifile(
  SSIFILE *sfp
);
 

/* Function: SSIRecommendMode()
 *
 * Purpose:  Examines the file and determines whether it should be
 *           indexed with large file support or not; returns
 *           SSI_OFFSET_I32 for most files, SSI_OFFSET_I64 for large
 *           files, or -1 on failure.
 *
 * Args:     file - name of file to check for size
 *
 * Returns:  -1 on failure (including case where file is too big)
 *           SSI_OFFSET_I32 for most files (<= 2^31-1 bytes)
 *           SSI_OFFSET_I64 for large files (> 2^31-1 bytes)
 */
int
SSIRecommendMode(
  char *file
);


/* Function: SSICreateIndex()
 *
 * Purpose:  Creates and initializes a SSI index structure.
 *           Sequence file offset type is specified by {mode}.
 *
 * Args:     mode    - SSI_OFFSET_I32 or SSI_OFFSET_I64, sequence file index mode.
 *
 * Returns:  ptr to new index structure, or NULL on failure.
 *           Caller is responsible for free'ing the returned
 *           structure with SSIFreeIndex().
 */
SSIINDEX*
SSICreateIndex(
  int mode
);


/* Function: SSIGetFilePosition()
 *
 * Purpose:  Fills {ret_offset} with the current disk
 *           offset of {fp}, relative to the start of the file.
 *           {mode} is set to either SSI_OFFSET_I32 or
 *           SSI_OFFSET_I64. If {mode} is _I32 (32 bit), just wraps
 *           a call to ftell(); otherwise deals with system-dependent
 *           details of 64-bit file offsets.
 *
 * Args:     fp         - open stream
 *           mode       - SSI_OFFSET_I32 or SSI_OFFSET_I64
 *           ret_offset - RETURN: file position      
 *
 * Returns:  0 on success. nonzero on error.
 */
int
SSIGetFilePosition(
  FILE *fp,
  int mode,
  SSIOFFSET *ret_offset
);


/* Function: SSIAddFileToIndex()
 *
 * Purpose:  Adds the sequence file {filename}, which is known to
 *           be in format {fmt}, to the index {g}. Creates and returns
 *           a unique filehandle {fh} for then associating primary keys
 *           with this file using SSIAddPrimaryKeyToIndex().
 *
 * Args:     g         - active index
 *           filename  - file to add
 *           fmt       - format code for this file (e.g. SQFILE_FASTA)
 *           ret_fh    - RETURN: unique handle for this file
 *
 * Returns:  0 on success; nonzero on error.
 */
int
SSIAddFileToIndex(
  SSIINDEX *g,
  char *filename,
  int fmt,
  int *ret_fh
);


/* Function: SSISetFileForSubseq()
 *
 * Purpose:  Set SSI_FAST_SUBSEQ for the file indicated by
 *           filehandle {fh} in the index {g}, setting
 *           parameters {bpl} and {rpl} to the values given.
 *           {bpl} is the number of bytes per sequence data line.
 *           {rpl} is the number of residues per sequence data line.
 *           Caller must be sure that {bpl} and {rpl} do not change
 *           on any line of any sequence record in the file
 *           (except for the last data line of each record). If
 *           this is not the case in this file, SSI_FAST_SUBSEQ
 *           will not work, and this routine should not be
 *           called.
 *
 * Args:     g    - the active index
 *           fh   - handle for file to set SSI_FAST_SUBSEQ on
 *           bpl  - bytes per data line
 *           rpl  - residues per data line
 *
 * Returns:  0 on success; 1 on error.
 */
int
SSISetFileForSubseq(
  SSIINDEX *g,
  int fh,
  int bpl,
  int rpl
);


/* Function: SSIAddPrimaryKeyToIndex()
 *
 * Purpose:  Put primary key {key} in the index {g}, while telling
 *           the index this primary key is in the file associated
 *           with filehandle {fh} (returned by a previous call
 *           to SSIAddFileToIndex()), and its record starts at
 *           position {r_off} in the file.
 *          
 *           {d_off} and {L} are optional; they may be left unset
 *           by passing NULL and 0, respectively. (If one is
 *           provided, both must be provided.) If they are provided,
 *           {d_off} gives the position of the first line of sequence
 *           data in the record, and {L} gives the length of
 *           the sequence in residues. They are used when
 *           SSI_FAST_SUBSEQ is set for this file. If SSI_FAST_SUBSEQ
 *           is not set for the file, {d_off} and {L} will be
 *           ignored by the index reading API even if they are stored
 *           by the index writing API, so it doesn't hurt for the
 *           indexing program to provide them; typically they
 *           won't know whether it's safe to set SSI_FAST_SUBSEQ
 *           for the whole file until the whole file has been
 *           read and every key has already been added to the index.
 *          
 * Args:     g      - active index
 *           key    - primary key to add
 *           fh     - handle on file that this key's in
 *           r_off  - offset to start of record
 *           d_off  - offset to start of sequence data
 *           L      - length of sequence, or 0
 *
 * Returns:  0 on success, nonzero on error.
 */
int
SSIAddPrimaryKeyToIndex(
  SSIINDEX *g,
  char *key,
  int fh,
  SSIOFFSET *r_off,
  SSIOFFSET *d_off,
  int L
);


/* Function: SSIAddSecondaryKeyToIndex()
 *
 * Purpose:  Puts secondary key {key} in the index {g}, associating
 *           it with primary key {pkey} that was previously
 *           registered by SSIAddPrimaryKeyToIndex().
 *
 * Args:     g    - active index
 *           key  - secondary key to add            
 *           pkey - primary key to associate this key with
 *
 * Returns:  0 on success, nonzero on failure.
 */
int
SSIAddSecondaryKeyToIndex(
  SSIINDEX *g,
  char *key,
  char *pkey
);


/* Function: SSIWriteIndex()
 *
 * Purpose:  Writes complete index {g} in SSI format to a
 *           binary file {file}. Does all          
 *           the overhead of sorting the primary and secondary keys,
 *           and maintaining the association of secondary keys
 *           with primary keys during and after the sort.
 *
 * Args:     file  - file to write to
 *           g     - index to sort & write out.     
 *
 * Returns:  0 on success, nonzero on error.
 */
/* needed for qsort() */
int
pkeysort(
  const void *k1,
  const void *k2
);


int
skeysort(
  const void *k1,
  const void *k2
);


int
SSIWriteIndex(
  char *file,
  SSIINDEX *g
);


/* Function: SSIFreeIndex()
 *
 * Purpose:  Free an index structure {g}.
 *
 * Args:     g  - ptr to an open index.
 *
 * Returns:  (void)
 */
void
SSIFreeIndex(
  SSIINDEX *g);


/* Function: SSIErrorString()
 *
 * Purpose:  Returns a ptr to an internal string corresponding
 *           to error {n}, a code returned from any of the
 *           functions in the API that return non-zero on error.
 *
 * Args:     n - error code
 *
 * Returns:  ptr to an internal string.
 */
char*
SSIErrorString(
  int n);


int
read_i16(
  FILE *fp,
  uint16_t *ret_result);


int
write_i16(
  FILE *fp,
  uint16_t n);


int
read_i32(
  FILE *fp,
  uint32_t *ret_result);


int
write_i32(
  FILE *fp,
  uint32_t n);


int
read_i64(
  FILE *fp,
  uint64_t *ret_result);


int
write_i64(
  FILE *fp,
  uint64_t n);


int     
read_offset(
  FILE *fp,
  char mode,
  SSIOFFSET *ret_offset);


int
write_offset(
  FILE *fp,
  SSIOFFSET *offset);


int
parse_pkey_info(
  char *buf,
  char mode,
  struct ssipkey_s *pkey);


int
parse_skey_info(
  char *buf,
  struct ssiskey_s *skey);


/* Function: binary_search()
 *
 * Purpose:  Find a key in a SSI index, by a binary search
 *           in an alphabetically sorted list of keys. If successful,
 *           return 0, and the index file is positioned to read
 *           the rest of the data for that key. Else returns nonzero.
 *
 * Args:     sfp    - an open SSIFILE
 *           key    - key to find
 *           klen   - key length to allocate (plen or slen from sfp)
 *           base   - base offset (poffset or soffset)
 *           recsize - size of each key record in bytes (precsize or srecsize)
 *           maxidx  - # of keys (nprimary or nsecondary)
 *
 * Returns:  0 on success, and leaves file positioned for reading remaining
 *           data for the key.
 *           Nonzero on failure:
 *                SSI_ERR_NO_SUCH_KEY  - that key's not in the index
 *                SSI_ERR_MALLOC       - a memory allocation failure
 *                SSI_ERR_NODATA       - an fread() failed
 */
int
binary_search(
  SSIFILE *sfp,
  char *key,
  int klen,
  SSIOFFSET *base,
  uint32_t recsize,
  uint32_t maxidx);


/* Function: indexfile_position()
 *
 * Purpose:  Position the open index file {sfp} at the start
 *           of record {n} in a list of records that starts at
 *           base offset {base}, where each record takes up {l}
 *           bytes. (e.g. the position is byte (base + n*l)).
 *
 * Args:     sfp - open SSIFILE
 *           base  - offset of record 0 (e.g. sfp->foffset)
 *           len   - size of each record in bytes (e.g. sfp->frecsize)
 *           n     - which record to get (e.g. 0..sfp->nfiles)
 *
 * Returns:  0 on success, non-zero on failure.
 */
int
indexfile_position(
  SSIFILE *sfp,
  SSIOFFSET *base,
  uint32_t len,
  uint32_t n);


/* Function: current_index_size()
 *
 * Purpose:  Calculates the size of the current index,
 *           in megabytes.
 */
uint64_t
current_index_size(
  SSIINDEX *g);


/* Function: activate_external_sort()
 *
 * Purpose:  Switch to external sort mode.
 *           Open file handles for external index files (ptmp, stmp).
 *           Flush current index information to these files.
 *           Free current memory, turn over control to the tmpfiles.
 *          
 * Return:   0 on success; non-zero on failure.
 */
int
activate_external_sort(
  SSIINDEX *g);


/*****************************************************************
 * Debugging API
 *****************************************************************/
void
SSIForceExternalSort(
  SSIINDEX *g);

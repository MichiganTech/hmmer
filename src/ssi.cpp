/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>

#include "endian.h"
#include "squid.h"
#include "ssi.h"
#include "file.h"

static uint32_t v20magic = 0xf3f3e9b1; /* SSI 1.0: "ssi1" + 0x80808080 */
static uint32_t v20swap  = 0xb1e9f3f3; /* byteswapped */


int
SSIOpen(
  char *filename, 
  SSIFILE **ret_sfp
){
  SSIFILE  *sfp = NULL;
  int       status;
  if ((sfp = malloc(sizeof(SSIFILE))) == NULL)   return SSI_ERR_MALLOC;
  if ((sfp->fp = fopen(filename, "rb")) == NULL) {
    free(sfp);
    return SSI_ERR_NOFILE;    
  }
  status = load_indexfile(sfp);
  *ret_sfp = sfp;
  return status;
}


int
load_indexfile(
  SSIFILE *sfp
){
  uint32_t   magic;
  uint16_t   i;		/* counter over files */
  int          status;		/* overall return status if an error is thrown */

  status = SSI_ERR_BADFORMAT; /* default: almost every kind of error is a bad format error */

  sfp->filename   = NULL;
  sfp->fileformat = NULL;
  sfp->fileflags  = NULL;
  sfp->bpl        = NULL;
  sfp->rpl        = NULL;
  sfp->nfiles     = 0;          
  if (! read_i32(sfp->fp, &magic))               {status = SSI_ERR_BADMAGIC;  goto FAILURE; }
  if (magic != v20magic && magic != v20swap)     {status = SSI_ERR_BADMAGIC;  goto FAILURE; }
  if (! read_i32(sfp->fp, &(sfp->flags))) goto FAILURE; 


  sfp->imode = (sfp->flags & SSI_USE64_INDEX) ? SSI_OFFSET_I64 : SSI_OFFSET_I32;
  sfp->smode = (sfp->flags & SSI_USE64) ?       SSI_OFFSET_I64 : SSI_OFFSET_I32;

  if (! read_i16(sfp->fp, &(sfp->nfiles)))     goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->nprimary)))   goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->nsecondary))) goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->flen)))       goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->plen)))       goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->slen)))       goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->frecsize)))   goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->precsize)))   goto FAILURE;
  if (! read_i32(sfp->fp, &(sfp->srecsize)))   goto FAILURE;
  
  if (! read_offset(sfp->fp, sfp->imode, &(sfp->foffset))) goto FAILURE;
  if (! read_offset(sfp->fp, sfp->imode, &(sfp->poffset))) goto FAILURE;
  if (! read_offset(sfp->fp, sfp->imode, &(sfp->soffset))) goto FAILURE;

  /* Read the file information and keep it.
   * We expect the number of files to be small, so reading it
   * once should be advantageous overall. If SSI ever had to
   * deal with large numbers of files, you'd probably want to
   * read file information on demand.
   */
  if (sfp->nfiles == 0)                                                   goto FAILURE;
  if ((sfp->filename=malloc(sizeof(char *)    *sfp->nfiles)) == NULL)   {status = SSI_ERR_MALLOC; goto FAILURE; }
  for (i = 0; i < sfp->nfiles; i++) sfp->filename[i] = NULL; 
  if ((sfp->fileformat=malloc(sizeof(uint32_t)*sfp->nfiles)) == NULL) {status = SSI_ERR_MALLOC; goto FAILURE; }
  if ((sfp->fileflags =malloc(sizeof(uint32_t)*sfp->nfiles)) == NULL) {status = SSI_ERR_MALLOC; goto FAILURE; }
  if ((sfp->bpl     =malloc(sizeof(uint32_t)*sfp->nfiles)) == NULL)   {status = SSI_ERR_MALLOC; goto FAILURE; }
  if ((sfp->rpl     =malloc(sizeof(uint32_t)*sfp->nfiles)) == NULL)   {status = SSI_ERR_MALLOC; goto FAILURE; }

  for (i = 0; i < sfp->nfiles; i++) 
    {
      /* We have to explicitly position, because header and file 
       * records may expand in the future; frecsize and foffset 
       * give us forwards compatibility. 
       */ 
      if (indexfile_position(sfp, &(sfp->foffset), sfp->frecsize, i) !=0)  goto FAILURE;
      if ((sfp->filename[i] =malloc(sizeof(char)*sfp->flen)) == NULL)        {status = SSI_ERR_MALLOC; goto FAILURE; }
      if (fread(sfp->filename[i],sizeof(char),sfp->flen, sfp->fp)!=sfp->flen) goto FAILURE;
      if (! read_i32(sfp->fp, &(sfp->fileformat[i])))                             goto FAILURE;
      if (! read_i32(sfp->fp, &(sfp->fileflags[i])))                              goto FAILURE;
      if (! read_i32(sfp->fp, &(sfp->bpl[i])))                                    goto FAILURE;
      if (! read_i32(sfp->fp, &(sfp->rpl[i])))                                    goto FAILURE;
    }
  
  /* Success. Return 0.
   */
  return 0;			

 FAILURE:
  /* Failure: free the damaged structure, return status code.
   */
  SSIClose(sfp);
  return status;
}


int
SSIGetOffsetByName(
  SSIFILE *sfp, 
  char *key, 
  int *ret_fh,
  SSIOFFSET *ret_offset
){
  int         status;
  uint16_t  fnum;

  /* Look in the primary keys.
   */
  status = binary_search(sfp, key, sfp->plen, &(sfp->poffset), sfp->precsize,
			 sfp->nprimary);
  if (status == 0) {		
    /* We found it as a primary key; get our data & return.
     */
    if (! read_i16(sfp->fp, &fnum)) return SSI_ERR_NODATA;
    *ret_fh = (int) fnum;
    if (! read_offset(sfp->fp, sfp->smode, ret_offset))  return SSI_ERR_NODATA;

    return 0;	/* success! (we don't need the other key data) */
  } else if (status == SSI_ERR_NO_SUCH_KEY) {
    /* Not in the primary keys? OK, try the secondary keys.
     */
    if (sfp->nsecondary > 0) {
      char *pkey;
      status = binary_search(sfp, key, sfp->slen, &(sfp->soffset), sfp->srecsize,
			     sfp->nsecondary);
      if (status != 0) return status;
      if ((pkey = malloc(sizeof(char) * sfp->plen)) == NULL) return SSI_ERR_MALLOC;
      if (fread(pkey, sizeof(char), sfp->plen, sfp->fp) != sfp->plen) return SSI_ERR_NODATA;

      status = SSIGetOffsetByName(sfp, pkey, ret_fh, ret_offset);
      free(pkey);
    }
    return status;

  } else return status;		
  /*NOTREACHED*/
}


int
SSIGetOffsetByNumber(
  SSIFILE *sfp, 
  int n, 
  int *ret_fh, 
  SSIOFFSET *ret_offset
){
  uint16_t fnum;
  char      *pkey;

  if (n >= sfp->nprimary) return SSI_ERR_NO_SUCH_KEY;
  if (indexfile_position(sfp, &(sfp->poffset), sfp->precsize, n) != 0) 
    return SSI_ERR_SEEK_FAILED;

  if ((pkey = malloc(sizeof(char) * sfp->plen)) == NULL) return SSI_ERR_MALLOC;
  if (fread(pkey, sizeof(char), sfp->plen, sfp->fp) != sfp->plen) return SSI_ERR_NODATA;
  if (! read_i16(sfp->fp, &fnum))                      return SSI_ERR_NODATA;
  if (! read_offset(sfp->fp, sfp->smode, ret_offset))  return SSI_ERR_NODATA;  
  *ret_fh = fnum;
  free(pkey);
  return 0;
}


int
SSIGetSubseqOffset(
  SSIFILE *sfp, 
  char *key, 
  int requested_start,
  int *ret_fh, 
  SSIOFFSET *record_offset,
  SSIOFFSET *data_offset, 
  int *ret_actual_start
){
  int        status;
  uint32_t len;
  int        r, b, i, l;	/* tmp variables for "clarity", to match docs */
  
  /* Look up the key. Rely on the fact that SSIGetOffsetByName()
   * leaves the index file positioned at the rest of the data for this key.
   */
  status = SSIGetOffsetByName(sfp, key, ret_fh, record_offset);
  if (status != 0) return status;

  /* Check that we're allowed to do subseq lookup on that file.
   */
  if (! (sfp->fileflags[*ret_fh] & SSI_FAST_SUBSEQ))
    return SSI_ERR_NO_SUBSEQS;

  /* Read the data we need for subseq lookup
   */
  if (! read_offset(sfp->fp, sfp->smode, data_offset)) return SSI_ERR_NODATA;
  if (! read_i32(sfp->fp, &len))                         return SSI_ERR_NODATA;

  /* Set up tmp variables for clarity of equations below,
   * and to make them match documentation (ssi-format.tex).
   */
  r = sfp->rpl[*ret_fh];    /* residues per line */
  b = sfp->bpl[*ret_fh];    /* bytes per line    */
  i = requested_start;	    /* start position 1..L */
  l = (i-1)/r;		    /* data line # (0..) that the residue is on */
  if (r == 0 || b == 0) return SSI_ERR_NO_SUBSEQS;
  if (i < 0 || i > len) return SSI_ERR_RANGE;
  
  /* When b = r+1, there's nothing but sequence on each data line (and the \0),
   * and we can find each residue precisely.
   */
  if (b == r+1) {
    if (sfp->smode == SSI_OFFSET_I32) {
      data_offset->mode    = SSI_OFFSET_I32;
      data_offset->off.i32 = data_offset->off.i32 + l*b + (i-1)%r;
    } else if (sfp->smode == SSI_OFFSET_I64) {
      data_offset->mode    = SSI_OFFSET_I64;
      data_offset->off.i64 = data_offset->off.i64 + l*b + (i-1)%r;
    } 
    *ret_actual_start = requested_start;
  } else { 
    /* else, there's other stuff on seq lines, so the best
     * we can do easily is to position at start of relevant line.
     */
    if (sfp->smode == SSI_OFFSET_I32) {
      data_offset->mode    = SSI_OFFSET_I32;
      data_offset->off.i32 = data_offset->off.i32 + l*b;
    } else if (sfp->smode == SSI_OFFSET_I64) {
      data_offset->mode    = SSI_OFFSET_I64;
      data_offset->off.i64 = data_offset->off.i64 + l*b;
    } 
    /* yes, the eq below is = 1 + (i-1)/r*r but it's not = i. that's an integer /. */
    *ret_actual_start = 1 + l*r;
  }
  return 0;
}


int
SSISetFilePosition(
  FILE *fp, 
  SSIOFFSET *offset
){
  if (offset->mode == SSI_OFFSET_I32) {
    if (fseek(fp, offset->off.i32, SEEK_SET) != 0)       return SSI_ERR_SEEK_FAILED;
  }
#ifndef HAS_64BIT_FILE_OFFSETS
  else return SSI_ERR_NO64BIT;
#elif defined HAVE_FSEEKO && SIZEOF_OFF_T == 8
  else if (fseeko(fp, offset->off.i64, SEEK_SET) != 0)   return SSI_ERR_SEEK_FAILED;
#elif defined HAVE_FSEEKO64 && SIZEOF_OFF64_T == 8
  else if (fseeko64(fp, offset->off.i64, SEEK_SET) != 0) return SSI_ERR_SEEK_FAILED;
#elif defined HAVE_FSEEK64
  else if (fseek64(fp, offset->off.i64, SEEK_SET) != 0)  return SSI_ERR_SEEK_FAILED;
#elif defined ARITHMETIC_FPOS_T && SIZEOF_FPOS_T == 8
  else if (fsetpos(fp, &(offset->off.i64)) != 0)         return SSI_ERR_SEEK_FAILED;
#endif
  return 0;
}


int
SSIFileInfo(
  SSIFILE *sfp, 
  int fh, 
  char **ret_filename, 
  int *ret_format
){
  if (fh < 0 || fh >= sfp->nfiles) return SSI_ERR_BADARG;
  *ret_filename = sfp->filename[fh];
  *ret_format   = sfp->fileformat[fh];
  return 0;
}


void
SSIClose(
  SSIFILE *sfp
){
  if (sfp != NULL) {
    clear_ssifile(sfp);
    if (sfp->fp       != NULL) fclose(sfp->fp);
    free(sfp);
  }
}  


void
clear_ssifile(
  SSIFILE *sfp
){
  int i;

  if (sfp->filename != NULL) {
    for (i = 0; i < sfp->nfiles; i++) 
      if (sfp->filename[i] != NULL) free(sfp->filename[i]);
    free(sfp->filename);
  }
  if (sfp->fileformat   != NULL) free(sfp->fileformat);
  if (sfp->fileflags    != NULL) free(sfp->fileflags);
  if (sfp->bpl          != NULL) free(sfp->bpl);
  if (sfp->rpl          != NULL) free(sfp->rpl);
}
  

int
SSIRecommendMode(
  char *file
){
#if HAVE_STAT64
  struct stat64 s1;
  if (stat64(file, &s1) == 0) {
    if (s1.st_size <= 2146483647L) return SSI_OFFSET_I32;
    else                           return SSI_OFFSET_I64;
  }
#else 
  struct stat s2;
  if (stat(file, &s2) == 0) {
    if (s2.st_size <= 2146483647L) return SSI_OFFSET_I32;
    else                           return SSI_OFFSET_I64;
  }
#endif
  return -1;
}
 

SSIINDEX*
SSICreateIndex(
  int mode
){
  SSIINDEX *g;

  g = NULL;
  if ((g = malloc(sizeof(SSIINDEX))) == NULL)                               goto FAILURE;
  g->smode    = mode;
  g->imode    = SSI_OFFSET_I32;	/* index always starts as 32-bit; may get upgraded later */
  g->external = false;
  g->max_ram  = SSI_MAXRAM;

#ifndef HAS_64BIT_FILE_OFFSETS
  if (mode == SSI_OFFSET_I64) 
    Die("\
Can't create a 64-bit SSI index on this system, sorry;\n\
I don't have 64-bit file offset functions available.\n");
#endif

  g->filenames  = NULL;
  g->fileformat = NULL;
  g->bpl        = NULL;
  g->rpl        = NULL;
  g->flen       = 0;
  g->nfiles     = 0;

  g->pkeys         = NULL;
  g->plen          = 0;
  g->nprimary      = 0;
  g->ptmpfile      = "tmp.ssi.1"; /* hardcoded, for now. */
  g->ptmp          = NULL;
  
  g->skeys         = NULL;
  g->slen          = 0;
  g->nsecondary    = 0;
  g->stmpfile      = "tmp.ssi.2"; /* hardcoded, for now. */
  g->stmp          = NULL;

  /* All mallocs must go after NULL initializations, because of the cleanup strategy;
   * we'll try to free anything non-NULL if a malloc fails.
   */
  if ((g->filenames = malloc(sizeof(char *)     * SSI_FILE_BLOCK)) == NULL) goto FAILURE;
  if ((g->fileformat= malloc(sizeof(uint32_t) * SSI_FILE_BLOCK)) == NULL) goto FAILURE; 
  if ((g->bpl       = malloc(sizeof(uint32_t) * SSI_FILE_BLOCK)) == NULL) goto FAILURE; 
  if ((g->rpl       = malloc(sizeof(uint32_t) * SSI_FILE_BLOCK)) == NULL) goto FAILURE; 
  
  if ((g->pkeys = malloc(sizeof(struct ssipkey_s)* SSI_KEY_BLOCK))== NULL)  goto FAILURE;
  if ((g->skeys = malloc(sizeof(struct ssipkey_s)* SSI_KEY_BLOCK))== NULL)  goto FAILURE;

  return g;

 FAILURE:
  SSIFreeIndex(g);		/* free the damaged structure */
  return NULL;
}


int 
SSIGetFilePosition(
  FILE *fp, 
  int mode, 
  SSIOFFSET *ret_offset
){
  if (mode == SSI_OFFSET_I32) 
    {
      ret_offset->mode    = SSI_OFFSET_I32;
      ret_offset->off.i32 = ftell(fp);
      if (ret_offset->off.i32 == -1) return SSI_ERR_TELL_FAILED;
    }
  else if (mode != SSI_OFFSET_I64) abort(); /* only happens on a coding error */
  else {
    ret_offset->mode    = SSI_OFFSET_I64;
#ifndef HAS_64BIT_FILE_OFFSETS
    return SSI_ERR_NO64BIT;
#elif defined HAVE_FTELLO && SIZEOF_OFF_T == 8
    if ((ret_offset->off.i64 = ftello(fp)) == -1)   return SSI_ERR_TELL_FAILED;
#elif defined HAVE_FTELLO64 && SIZEOF_OFF64_T == 8
    if ((ret_offset->off.i64 = ftello64(fp)) == -1) return SSI_ERR_TELL_FAILED;
#elif defined HAVE_FTELL64
    if ((ret_offset->off.i64 = ftell64(fp)) == -1)  return SSI_ERR_TELL_FAILED;
#elif defined ARITHMETIC_FPOS_T && SIZEOF_FPOS_T == 8
    if (fgetpos(fp, &(ret_offset->off.i64)) != 0)   return SSI_ERR_TELL_FAILED;
#endif
  }
  return 0;
}


int
SSIAddFileToIndex(
  SSIINDEX *g, 
  char *filename, 
  int fmt, 
  int *ret_fh
){
  int n;
  
  if (g->nfiles >= SSI_MAXFILES) return SSI_ERR_TOOMANY_FILES;

  n = strlen(filename);
  if ((n+1) > g->flen) g->flen = n+1;

  g->filenames[g->nfiles]  = FileTail(filename, false);
  g->fileformat[g->nfiles] = fmt;
  g->bpl[g->nfiles]        = 0;
  g->rpl[g->nfiles]        = 0;
  *ret_fh                  = g->nfiles;   /* handle is simply = file number */
  g->nfiles++;

  if (g->nfiles % SSI_FILE_BLOCK == 0) {
    g->filenames = realloc(g->filenames,  sizeof(char *) * (g->nfiles+SSI_FILE_BLOCK));
    if (g->filenames == NULL) return SSI_ERR_MALLOC;
    g->fileformat= realloc(g->fileformat, sizeof(uint32_t) * (g->nfiles+SSI_FILE_BLOCK));
    if (g->fileformat == NULL) return SSI_ERR_MALLOC;
    g->bpl       = realloc(g->bpl,        sizeof(uint32_t) * (g->nfiles+SSI_FILE_BLOCK));
    if (g->bpl == NULL) return SSI_ERR_MALLOC;
    g->rpl       = realloc(g->rpl,        sizeof(uint32_t) * (g->nfiles+SSI_FILE_BLOCK));
    if (g->rpl == NULL) return SSI_ERR_MALLOC;
  }
  return 0;
}


int
SSISetFileForSubseq(
  SSIINDEX *g, 
  int fh, 
  int bpl, 
  int rpl
){
  if (fh < 0 || fh >= g->nfiles) return SSI_ERR_BADARG;
  if (bpl <= 0 || rpl <= 0)      return SSI_ERR_BADARG;
  g->bpl[fh] = bpl;
  g->rpl[fh] = rpl;
  return 0;
}


int
SSIAddPrimaryKeyToIndex(
  SSIINDEX *g, 
  char *key, 
  int fh,
  SSIOFFSET *r_off, 
  SSIOFFSET *d_off, 
  int L
){
  int n;			/* a string length */
  
  if (fh >= SSI_MAXFILES)         return SSI_ERR_TOOMANY_FILES;
  if (g->nprimary >= SSI_MAXKEYS) return SSI_ERR_TOOMANY_KEYS;
  if (L > 0 && d_off == NULL) abort(); /* need both. */

  /* Before adding the key: check how big our index is.
   * If it's getting too large, switch to external mode.
   */
  if (!g->external && current_index_size(g) >= g->max_ram) 
    if (activate_external_sort(g) != 0)  return SSI_ERR_NOFILE;

  /* Update maximum pkey length, if needed.
   */
  n = strlen(key);
  if ((n+1) > g->plen) g->plen = n+1;

  /* External mode? Simply append to disk...
   */
  if (g->external) {
    if (g->smode == SSI_OFFSET_I32) {
      fprintf(g->ptmp, "%s\t%d\t%lu\t%lu\t%lu\n", 
	      key, fh, (unsigned long) r_off->off.i32, 
	      (unsigned long) (d_off == NULL? 0 : d_off->off.i32),
	      (unsigned long) L);
    } else {
      fprintf(g->ptmp, "%s\t%d\t%lu\t%lu\t%lu\n", 
	      key, fh, r_off->off.i64, 
	      d_off == NULL? 0 : d_off->off.i64, 
	      (unsigned long) L);
    }
    g->nprimary++;
    return 0;
  }

  /* Else: internal mode, keep keys in memory...
   */
  if ((g->pkeys[g->nprimary].key = strdup(key)) == NULL) return SSI_ERR_MALLOC;
  g->pkeys[g->nprimary].fnum  = (uint16_t) fh;
  g->pkeys[g->nprimary].r_off = *r_off;
  if (d_off != NULL && L > 0) {
    g->pkeys[g->nprimary].d_off = *d_off;
    g->pkeys[g->nprimary].len   = L;
  } else {
	/* yeah, this looks stupid, but look: we have to give a valid
           looking, non-NULL d_off of some sort, or writes will fail. 
           It's going to be unused anyway. */
    g->pkeys[g->nprimary].d_off = *r_off;
    g->pkeys[g->nprimary].len   = 0;
  }
  g->nprimary++;

  if (g->nprimary % SSI_KEY_BLOCK == 0) {
    g->pkeys = realloc(g->pkeys, sizeof(struct ssipkey_s) * (g->nprimary+SSI_KEY_BLOCK));
    if (g->pkeys == NULL) return SSI_ERR_MALLOC;
  }
  return 0;
}


int
SSIAddSecondaryKeyToIndex(
  SSIINDEX *g, 
  char *key, 
  char *pkey
){
  int n;			/* a string length */
  
  if (g->nsecondary >= SSI_MAXKEYS) return SSI_ERR_TOOMANY_KEYS;

  /* Before adding the key: check how big our index is.
   * If it's getting too large, switch to external mode.
   */
  if (!g->external && current_index_size(g) >= g->max_ram) 
    if (activate_external_sort(g) != 0)  return SSI_ERR_NOFILE;

  /* Update maximum secondary key length, if necessary.
   */
  n = strlen(key);
  if ((n+1) > g->slen) g->slen = n+1;

  /* if external mode: write info to disk.
   */
  if (g->external) {
    fprintf(g->stmp, "%s\t%s\n", key, pkey);
    g->nsecondary++;
    return 0;
  }

  /* else, internal mode... store info in memory.
   */
  if ((g->skeys[g->nsecondary].key  = strdup(key))   == NULL) return SSI_ERR_MALLOC;
  if ((g->skeys[g->nsecondary].pkey = strdup(pkey)) == NULL) return SSI_ERR_MALLOC;
  g->nsecondary++;

  if (g->nsecondary % SSI_KEY_BLOCK == 0) {
    g->skeys = realloc(g->skeys, sizeof(struct ssiskey_s) * (g->nsecondary+SSI_KEY_BLOCK));
    if (g->skeys == NULL) return SSI_ERR_MALLOC;
  }
  return 0;
}


int 
pkeysort(
  const void *k1, 
  const void *k2
){
  struct ssipkey_s *key1;
  struct ssipkey_s *key2;
  key1 = (struct ssipkey_s *) k1;
  key2 = (struct ssipkey_s *) k2;
  return strcmp(key1->key, key2->key);
}


int 
skeysort(
  const void *k1, 
  const void *k2
){
  struct ssiskey_s *key1;
  struct ssiskey_s *key2;
  key1 = (struct ssiskey_s *) k1;
  key2 = (struct ssiskey_s *) k2;
  return strcmp(key1->key, key2->key);
}


int
SSIWriteIndex(
  char *file, 
  SSIINDEX *g
){
  FILE      *fp;
  int        status;
  int        i;
  uint32_t header_flags, file_flags;
  uint32_t frecsize, precsize, srecsize;
  uint64_t foffset, poffset, soffset;
  char       *s, *s2;

  if ((fp = fopen(file,"wb")) == NULL) return SSI_ERR_NOFILE;
  status = 0;

  /* How big is the index? If it's going to be > 2GB, we need
   * to flip to 64-bit index mode. 2047 (instead of 2048) gives us
   * some slop room.
   * die'ing here is pretty brutal - if we flip to 64-bit index
   * mode, we hve 100's of millions of keys, so we've processed
   * a long time before reaching this point. Ah well.
   */
  if (current_index_size(g) >= 2047) {
    g->imode = SSI_OFFSET_I64;
#ifndef HAS_64BIT_FILE_OFFSETS
    Die("\
Can't switch to 64-bit SSI index mode on this system, sorry;\n\
I don't have 64-bit file offset functions available.\n");
#endif
  }

  /* Magic-looking numbers come from adding up sizes 
   * of things in bytes
   */
  frecsize = 16 + g->flen;
  precsize = (g->smode == SSI_OFFSET_I64) ? 22+g->plen : 14+g->plen;
  srecsize = g->slen + g->plen;

  header_flags = 0;
  if (g->smode == SSI_OFFSET_I64) header_flags |= SSI_USE64;
  if (g->imode == SSI_OFFSET_I64) header_flags |= SSI_USE64_INDEX;

  /* Magic-looking numbers again come from adding up sizes 
   * of things in bytes
   */
  foffset = (header_flags & SSI_USE64_INDEX) ? 66 : 54;
  poffset = foffset + frecsize*g->nfiles;
  soffset = poffset + precsize*g->nprimary;
  
  /* Sort the keys
   * If external mode, make system calls to UNIX/POSIX "sort" in place, then
   * open new sorted files for reading thru ptmp and stmp handles.
   * If internal mode, call qsort.
   * 
   * Note that you'd better force a POSIX locale for the sort; else,
   * some silly distro (e.g. Mandrake Linux >=8.1) may have specified
   * LC_COLLATE=en_US, and this'll give a sort "bug" in which it doesn't
   * sort by byte order.
   */
  if (g->external) {
    char cmd[1024];

    fclose(g->ptmp);
    g->ptmp = NULL;
    sprintf(cmd, "env LC_ALL=POSIX sort -o %s %s\n", g->ptmpfile, g->ptmpfile);
    if ((status = system(cmd)) != 0) return SSI_ERR_EXTERNAL_SORT;
    if ((g->ptmp = fopen(g->ptmpfile, "r")) == NULL) return SSI_ERR_EXTERNAL_SORT;

    fclose(g->stmp);
    g->stmp = NULL;
    sprintf(cmd, "env LC_ALL=POSIX sort -o %s %s\n", g->stmpfile, g->stmpfile);
    if ((status = system(cmd)) != 0) return SSI_ERR_EXTERNAL_SORT;
    if ((g->stmp = fopen(g->stmpfile, "r")) == NULL) return SSI_ERR_EXTERNAL_SORT;
  } else {
    qsort((void *) g->pkeys, g->nprimary,   sizeof(struct ssipkey_s), pkeysort); 
    qsort((void *) g->skeys, g->nsecondary, sizeof(struct ssiskey_s), skeysort); 
  }

  /* Write the header
   */
  if (! write_i32(fp, v20magic))      return SSI_ERR_FWRITE;
  if (! write_i32(fp, header_flags))  return SSI_ERR_FWRITE;
  if (! write_i16(fp, g->nfiles))     return SSI_ERR_FWRITE;
  if (! write_i32(fp, g->nprimary))   return SSI_ERR_FWRITE;
  if (! write_i32(fp, g->nsecondary)) return SSI_ERR_FWRITE;
  if (! write_i32(fp, g->flen))       return SSI_ERR_FWRITE;
  if (! write_i32(fp, g->plen))       return SSI_ERR_FWRITE;
  if (! write_i32(fp, g->slen))       return SSI_ERR_FWRITE;
  if (! write_i32(fp, frecsize))      return SSI_ERR_FWRITE;
  if (! write_i32(fp, precsize))      return SSI_ERR_FWRITE;
  if (! write_i32(fp, srecsize))      return SSI_ERR_FWRITE;
  if (g->imode == SSI_OFFSET_I32) {
    if (! write_i32(fp, foffset))     return SSI_ERR_FWRITE;
    if (! write_i32(fp, poffset))     return SSI_ERR_FWRITE;
    if (! write_i32(fp, soffset))     return SSI_ERR_FWRITE;
  } else {
    if (! write_i64(fp, foffset))     return SSI_ERR_FWRITE;
    if (! write_i64(fp, poffset))     return SSI_ERR_FWRITE;
    if (! write_i64(fp, soffset))     return SSI_ERR_FWRITE;
  }

  /* The file section
   */
  if ((s = malloc(sizeof(char) * g->flen)) == NULL) return SSI_ERR_MALLOC;
  for (i = 0; i < g->nfiles; i++)
    {
      file_flags = 0;
      if (g->bpl[i] > 0 && g->rpl[i] > 0) file_flags |= SSI_FAST_SUBSEQ;
      
      strcpy(s, g->filenames[i]);
      if (fwrite(s, sizeof(char), g->flen, fp) != g->flen) return SSI_ERR_FWRITE;
      if (! write_i32(fp, g->fileformat[i]))               return SSI_ERR_FWRITE;
      if (! write_i32(fp, file_flags))                     return SSI_ERR_FWRITE;
      if (! write_i32(fp, g->bpl[i]))                      return SSI_ERR_FWRITE;
      if (! write_i32(fp, g->rpl[i]))                      return SSI_ERR_FWRITE;
    }
  free(s);

  /* The primary key section
   */
  if ((s = malloc(sizeof(char) * g->plen)) == NULL) return SSI_ERR_MALLOC;
  if (g->external) {
    char *buf    = NULL;
    size_t buflen = 0;
    struct ssipkey_s pkey;
    for (i = 0; i < g->nprimary; i++) 
      {
	if (getline(&buf, &buflen, g->ptmp) == NULL)       return SSI_ERR_NODATA;
	if (parse_pkey_info(buf, g->smode, &pkey) != 0)      return SSI_ERR_BADFORMAT;
	strcpy(s, pkey.key);
	if (fwrite(s, sizeof(char), g->plen, fp) != g->plen) return SSI_ERR_FWRITE;
	if (! write_i16(   fp, pkey.fnum))                   return SSI_ERR_FWRITE;
	if (! write_offset(fp, &(pkey.r_off)))               return SSI_ERR_FWRITE;
	if (! write_offset(fp, &(pkey.d_off)))               return SSI_ERR_FWRITE;
	if (! write_i32(   fp, pkey.len))                    return SSI_ERR_FWRITE;
      }
    free(buf);
  } else {
    for (i = 0; i < g->nprimary; i++)
      {
	strcpy(s, g->pkeys[i].key);
	if (fwrite(s, sizeof(char), g->plen, fp) != g->plen) return SSI_ERR_FWRITE;
	if (! write_i16(   fp, g->pkeys[i].fnum))            return SSI_ERR_FWRITE;
	if (! write_offset(fp, &(g->pkeys[i].r_off)))        return SSI_ERR_FWRITE;
	if (! write_offset(fp, &(g->pkeys[i].d_off)))        return SSI_ERR_FWRITE;
	if (! write_i32(   fp, g->pkeys[i].len))             return SSI_ERR_FWRITE;
      }
  }

  /* The secondary key section
   */
  if (g->nsecondary > 0) {
    if ((s2  = malloc(sizeof(char) * g->slen)) == NULL) return SSI_ERR_MALLOC;

    if (g->external) {
      struct ssiskey_s skey;
      char *buf  = NULL;
      size_t n    = 0;

      for (i = 0; i < g->nsecondary; i++)
	{
	  if (getline(&buf, &n, g->stmp) == NULL)  return SSI_ERR_NODATA;
	  if (parse_skey_info(buf, &skey) != 0)           return SSI_ERR_BADFORMAT;
	  strcpy(s2, skey.key);
	  strcpy(s,  skey.pkey);
	  if (fwrite(s2, sizeof(char), g->slen, fp) != g->slen) return SSI_ERR_FWRITE;
	  if (fwrite(s,  sizeof(char), g->plen, fp) != g->plen) return SSI_ERR_FWRITE;
	}
      free(buf);
    } else {
      for (i = 0; i < g->nsecondary; i++)
	{
	  strcpy(s2, g->skeys[i].key);
	  strcpy(s,  g->skeys[i].pkey);
	  if (fwrite(s2, sizeof(char), g->slen, fp) != g->slen) return SSI_ERR_FWRITE;
	  if (fwrite(s,  sizeof(char), g->plen, fp) != g->plen) return SSI_ERR_FWRITE;
	} 
    }
    free(s2);
  }

  free(s);
  fclose(fp);
  return status;
}


void
SSIFreeIndex(
  SSIINDEX *g
){
  int i;
  if (g != NULL) 
    {
      if (g->external == false) {
	for (i = 0; i < g->nprimary;   i++) free(g->pkeys[i].key);
	for (i = 0; i < g->nsecondary; i++) free(g->skeys[i].key);
	for (i = 0; i < g->nsecondary; i++) free(g->skeys[i].pkey);
	if (g->pkeys       != NULL)         free(g->pkeys);       	
	if (g->skeys       != NULL)         free(g->skeys);       
      } else {
	if (g->ptmp        != NULL)         fclose(g->ptmp);
	if (g->stmp        != NULL)         fclose(g->stmp);       
#if DEBUGLEVEL == 0
	remove(g->ptmpfile);
	remove(g->stmpfile);
#endif
      }
      for (i = 0; i < g->nfiles;     i++) free(g->filenames[i]);
      if (g->filenames   != NULL)         free(g->filenames);
      if (g->fileformat  != NULL)         free(g->fileformat);
      if (g->bpl         != NULL)         free(g->bpl);       
      if (g->rpl         != NULL)         free(g->rpl);       
      free(g);
    }
}


char*
SSIErrorString(
  int n
){
  switch (n) {
  case SSI_ERR_OK:            return "ok (no error)"; 
  case SSI_ERR_NODATA:        return "no data, fread() failed";
  case SSI_ERR_NO_SUCH_KEY:   return "no such key";
  case SSI_ERR_MALLOC:        return "out of memory, malloc() failed";
  case SSI_ERR_NOFILE:        return "file not found, fopen() failed";
  case SSI_ERR_BADMAGIC:      return "not a SSI file? (bad magic)"; 
  case SSI_ERR_BADFORMAT:     return "corrupt format? unexpected data";
  case SSI_ERR_NO64BIT:       return "no large file support for this system";
  case SSI_ERR_SEEK_FAILED:   return "failed to reposition on disk";
  case SSI_ERR_TELL_FAILED:   return "failed to get file position on disk";
  case SSI_ERR_NO_SUBSEQS:    return "no fast subseq support for this seqfile";
  case SSI_ERR_RANGE:         return "subseq start is out of range";
  case SSI_ERR_BADARG:        return "an argument is out of range";
  case SSI_ERR_TOOMANY_FILES: return "number of files exceeds limit";
  case SSI_ERR_TOOMANY_KEYS:  return "number of keys exceeds limit";
  case SSI_ERR_FWRITE:        return "an fwrite() failed";
  case SSI_ERR_EXTERNAL_SORT: return "some problem with external sorting";
  default:                    return "unrecognized code";
  }
  /*NOTREACHED*/
}


int
read_i16(
  FILE *fp, 
  uint16_t *ret_result
){
  uint16_t result;
  if (fread(&result, sizeof(uint16_t), 1, fp) != 1) return 0;
  *ret_result = ntohs(result);
  return 1;
}


int
write_i16(
  FILE *fp, 
  uint16_t n
){
  n = htons(n);
  if (fwrite(&n, sizeof(uint16_t), 1, fp) != 1) return 0;
  return 1;
}


int
read_i32(
  FILE *fp, 
  uint32_t *ret_result
){
  uint32_t result;
  if (fread(&result, sizeof(uint32_t), 1, fp) != 1) return 0;
  *ret_result = ntohl(result);
  return 1;
}


int
write_i32(
  FILE *fp, 
  uint32_t n
){
  n = htonl(n);
  if (fwrite(&n, sizeof(uint32_t), 1, fp) != 1) return 0;
  return 1;
}


int
read_i64(
  FILE *fp, 
  uint64_t *ret_result
){
  uint64_t result;
  if (fread(&result, sizeof(uint64_t), 1, fp) != 1) return 0;
  *ret_result = le64toh(result);
  return 1;
}


int
write_i64(
  FILE *fp, 
  uint64_t n
){
  n = htole64(n);
  if (fwrite(&n, sizeof(uint64_t), 1, fp) != 1) return 0;
  return 1;
}


int			
read_offset(
  FILE *fp, 
  char mode, 
  SSIOFFSET *ret_offset
){
  if (mode == SSI_OFFSET_I32) {
    ret_offset->mode = SSI_OFFSET_I32;
    if (! read_i32(fp, &(ret_offset->off.i32))) return 0;
  } else if (mode == SSI_OFFSET_I64) {
    ret_offset->mode = SSI_OFFSET_I64;
    if (! read_i64(fp, &(ret_offset->off.i64))) return 0;
  } else return 0;

  return 1;
}


int
write_offset(
  FILE *fp, 
  SSIOFFSET *offset
){
  if      (offset->mode == SSI_OFFSET_I32) return write_i32(fp, offset->off.i32);
  else if (offset->mode == SSI_OFFSET_I64) return write_i64(fp, offset->off.i64);
  else abort();
  /*UNREACHED*/
  return 1; /* silence bitchy compilers */
}


int
parse_pkey_info(
  char *buf, 
  char mode, 
  struct ssipkey_s *pkey
){
  char *s, *tok;
  
  s = buf;
  if ((tok = strtok(s, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;  
  pkey->key  = tok;
  if ((tok = strtok(NULL, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;  
  pkey->fnum = (uint16_t) atoi(tok);

  if (mode == SSI_OFFSET_I32) {
    if ((tok = strtok(NULL, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;  
    pkey->r_off.mode = mode;
    pkey->r_off.off.i32  = (uint32_t) strtoul(tok, NULL, 10);
    if ((tok = strtok(NULL, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;  
    pkey->d_off.mode = mode;
    pkey->d_off.off.i32  = (uint32_t) strtoul(tok, NULL, 10);
  }else {
    if ((tok = strtok(NULL, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;  
    pkey->r_off.mode = mode;
    pkey->r_off.off.i64  = (uint64_t) strtoull(tok, NULL, 10);
    if ((tok = strtok(NULL, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;  
    pkey->d_off.mode = mode;
    pkey->d_off.off.i64  = (uint64_t) strtoull(tok, NULL, 10);
  }
  if ((tok = strtok(NULL, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;
  pkey->len = (uint32_t) strtoul(tok, NULL, 10);

  return 0;
}


int
parse_skey_info(
  char *buf, 
  struct ssiskey_s *skey
){
  char *s, *tok;
  
  s = buf;
  if ((tok = strtok(s, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;
  skey->key = tok;
  if ((tok = strtok(NULL, "\t\n")) == NULL) return SSI_ERR_BADFORMAT;
  skey->pkey = tok;
  return 0;
}


int
binary_search(
  SSIFILE *sfp, 
  char *key, 
  int klen, 
  SSIOFFSET *base,
  uint32_t recsize, 
  uint32_t maxidx
){
  char        *name;
  uint32_t   left, right, mid;
  int          cmp;
  int          status;
  
  if (maxidx == 0) return SSI_ERR_NO_SUCH_KEY; /* special case: empty index */
  if ((name = malloc (sizeof(char)*klen)) == NULL) return SSI_ERR_MALLOC;
  left  = 0;
  right = maxidx-1;
  while (1) {			/* A binary search: */
    mid   = (left+right) / 2;	/* careful here. only works because
				   we limit unsigned vars to signed ranges. */
    if ((status = indexfile_position(sfp, base, recsize, mid)) != 0)
      { free(name); return status; }
    if (fread(name, sizeof(char), klen, sfp->fp) != klen) 
      { free(name); return SSI_ERR_NODATA; }
    cmp = strcmp(name, key);
    if      (cmp == 0) break;	          /* found it!              */
    else if (left >= right)	          /* oops, missed it; fail  */
      { free(name); return SSI_ERR_NO_SUCH_KEY; }
    else if (cmp < 0)       left  = mid+1; /* it's right of mid     */
    else if (cmp > 0) {
      if (mid == 0) { free(name); return SSI_ERR_NO_SUCH_KEY; } /* special case, beware */
      else right = mid-1;                  /* it's left of mid      */
    }
  }
  free(name);
  return 0;			/* and sfp->fp is positioned... */
}


int
indexfile_position(
  SSIFILE *sfp, 
  SSIOFFSET *base, 
  uint32_t len, 
  uint32_t n
){
  SSIOFFSET pos;
  int       status;

  if (base->mode == SSI_OFFSET_I32) {
    pos.mode    = SSI_OFFSET_I32;
    pos.off.i32 = base->off.i32 + n*len;
  } else if (base->mode == SSI_OFFSET_I64) {
    pos.mode    = SSI_OFFSET_I64;
    pos.off.i64 = base->off.i64 + n*len;
  } else return 0;
  if ((status = SSISetFilePosition(sfp->fp, &pos)) != 0) return status;
  return 0;
}


uint64_t 
current_index_size(
  SSIINDEX *g
){
  uint64_t frecsize, precsize, srecsize;
  uint64_t total;

  /* Magic-looking numbers come from adding up sizes 
   * of things in bytes
   */
  frecsize = 16 + g->flen;
  precsize = (g->smode == SSI_OFFSET_I64) ? 22+g->plen : 14+g->plen;
  srecsize = g->plen+g->slen;
  total = (66L +		       /* header size, if 64bit index offsets */
	   frecsize * g->nfiles +      /* file section size                   */
	   precsize * g->nprimary +    /* primary key section size            */
	   srecsize * g->nsecondary) / /* secondary key section size          */
          1048576L;
  return total;
}


int
activate_external_sort(
  SSIINDEX *g
){
  int i;
				/* it's a bit late to be checking this, but... */
  if (g->external)             return 0; /* we already are external, fool */
  if (FileExists(g->ptmpfile)) return 1;	 
  if (FileExists(g->stmpfile)) return 1;
  if ((g->ptmp = fopen(g->ptmpfile, "w")) == NULL) return 1;
  if ((g->stmp = fopen(g->stmpfile, "w")) == NULL) return 1;

  /* Flush the current indices.
   */
  //SQD_DPRINTF1(("Switching to external sort - flushing ssiindex to disk...\n"));
  for (i = 0; i < g->nprimary; i++) {
    if (g->smode == SSI_OFFSET_I32) {
      fprintf(g->ptmp, "%s\t%u\t%lu\t%lu\t%lu\n", 
	      g->pkeys[i].key, g->pkeys[i].fnum,
	      (unsigned long) g->pkeys[i].r_off.off.i32, 
	      (unsigned long) g->pkeys[i].d_off.off.i32, 
	      (unsigned long) g->pkeys[i].len);
    } else {
      fprintf(g->ptmp, "%s\t%u\t%llu\t%llu\t%lu\n", 
	      g->pkeys[i].key, g->pkeys[i].fnum,
	      (unsigned long long) g->pkeys[i].r_off.off.i64, 
	      (unsigned long long) g->pkeys[i].d_off.off.i64, 
	      (unsigned long) g->pkeys[i].len);
    }
  }
  for (i = 0; i < g->nsecondary; i++)
    fprintf(g->stmp, "%s\t%s\n", g->skeys[i].key, g->skeys[i].pkey);
  
  /* Free the memory now that we've flushed our lists to disk
   */
  for (i = 0; i < g->nprimary;   i++) free(g->pkeys[i].key);
  for (i = 0; i < g->nsecondary; i++) free(g->skeys[i].key);
  for (i = 0; i < g->nsecondary; i++) free(g->skeys[i].pkey);
  if (g->pkeys       != NULL)         free(g->pkeys);       	
  if (g->skeys       != NULL)         free(g->skeys);       
  g->pkeys = NULL;
  g->skeys = NULL;

  /* Turn control over to external accumulation mode.
   */
  g->external = true;
  return 0;
}


/*****************************************************************
 * Debugging API
 *****************************************************************/
void
SSIForceExternalSort(SSIINDEX *g)
{
  if (activate_external_sort(g) != 0)
    Die("failed to turn external sorting on.");
}


/*****************************************************************
 * Test driving mode
 *****************************************************************/
#ifdef MUGGINS_LETS_ME_SLEEP 
/* Minimally: 
   cc -g -Wall -o shiva -DDEBUGLEVEL=1 -DMUGGINS_LETS_ME_SLEEP ssi.c sqerror.c sre_string.c types.c sre_ctype.c sre_math.c file.c -lm 
*/

int
main(int argc, char **argv)
{
  char      name[32], accession[32];
  SSIINDEX *ssi;
  int       mode;
  SSIOFFSET r_off, d_off;
  FILE     *ofp;
  int       i;
  int       fh;			/* a file handle */
  int       status;		/* return status from a SSI call */
  
  mode = SSI_OFFSET_I32;
  if ((ssi = SSICreateIndex(mode)) == NULL)
    Die("Failed to allocate SSI index");

  /* Generate two FASTA files, tmp.0 and tmp.1, and index them.
   */
  if ((ofp = fopen("tmp.0", "w")) == NULL) 
    Die("failed to open tmp.0");
  if ((status = SSIAddFileToIndex(ssi, "tmp.0", SQFILE_FASTA, &fh)) != 0)
    Die("SSIAddFileToIndex() failed: %s", SSIErrorString(status));
  for (i = 0; i < 10; i++) {
    if ((status = SSIGetFilePosition(ofp, mode, &r_off)) != 0)
      Die("SSIGetFilePosition() failed: %s", SSIErrorString(status));
    sprintf(name, "seq%d", i);
    sprintf(accession, "ac%d", i);
    fprintf(ofp, ">%s [%s] Description? we don't need no steenking description.\n", 
	    name, accession);
    if ((status = SSIGetFilePosition(ofp, mode, &d_off)) != 0) 
      Die("SSIGetFilePosition() failed: %s", SSIErrorString(status));
    fprintf(ofp, "AAAAAAAAAA\n");
    fprintf(ofp, "CCCCCCCCCC\n");
    fprintf(ofp, "GGGGGGGGGG\n");
    fprintf(ofp, "TTTTTTTTTT\n");

    if ((status = SSIAddPrimaryKeyToIndex(ssi, name, fh, &r_off, &d_off, 40)) != 0)
      Die("SSIAddPrimaryKeyToIndex() failed: %s", SSIErrorString(status));
    if ((status = SSIAddSecondaryKeyToIndex(ssi, accession, name)) != 0)
      Die("SSIAddSecondaryKeyToIndex() failed: %s", SSIErrorString(status));
  }
  SSISetFileForSubseq(ssi, fh, 11, 10);
  fclose(ofp);
  
  if ((ofp = fopen("tmp.1", "w")) == NULL) 
    Die("failed to open tmp.1");
  if ((status = SSIAddFileToIndex(ssi, "tmp.1", SQFILE_FASTA, &fh)) != 0)
    Die("SSIAddFileToIndex() failed: %s", SSIErrorString(status));
  for (i = 10; i < 20; i++) {
    if ((status = SSIGetFilePosition(ofp, mode, &r_off)) != 0)
      Die("SSIGetFilePosition() failed: %s", SSIErrorString(status));
    sprintf(name, "seq%d", i);
    sprintf(accession, "ac%d", i);
    fprintf(ofp, ">%s [%s] i/o, i/o, it's off to disk we go.\n", 
	    name, accession);
    if ((status = SSIGetFilePosition(ofp, mode, &d_off)) != 0)
      Die("SSIGetFilePosition() failed: %s", SSIErrorString(status));
    fprintf(ofp, "AAAAAAAAAA 10\n");
    fprintf(ofp, "CCCCCCCCCC 20\n");
    fprintf(ofp, "GGGGGGGGGG 30\n");
    fprintf(ofp, "TTTTTTTTTT 40\n");

    if ((status = SSIAddPrimaryKeyToIndex(ssi, name, fh, &r_off, &d_off, 40)) != 0)
      Die("SSIAddPrimaryKeyToIndex() failed: %s", SSIErrorString(status));
    if ((status = SSIAddSecondaryKeyToIndex(ssi, accession, name)) != 0)
      Die("SSIAddSecondaryKeyToIndex() failed: %s", SSIErrorString(status));
  }
  SSISetFileForSubseq(ssi, fh, 14, 10);
  fclose(ofp);
  
  /* Write the index to tmp.ssi
   */  
  if ((status = SSIWriteIndex("tmp.ssi", ssi)) != 0) 
    Die("SSIWriteIndex() failed: %s", SSIErrorString(status));
  SSIFreeIndex(ssi);

  /* Now reopen the index and run some tests.
   */
  exit(0);
}


#endif /* test driving code */




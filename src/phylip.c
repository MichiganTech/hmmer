/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* phylip.c
 * 
 * Import/export of PHYLIP interleaved multiple sequence alignment
 * format files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"
#include "msa.h"


MSA *
ReadPhylip(
  MSAFILE *afp
){
  MSA  *msa;
  char *s, *s1, *s2;
  char  name[11];   /* seq name max len = 10 char */
  int   nseq;
  int   idx;      /* index of current sequence */
  int   nblock;
  
  if (feof(afp->f)) return NULL;

  /* Skip until we see a nonblank line; it's the header,
   * containing nseq/alen
   */
  nseq = 0;
  while ((s = MSAFileGetLine(afp)) != NULL)
    {
      if ((s1 = strtok(s, WHITESPACE)) == NULL) continue;
      if ((s2 = strtok(s, WHITESPACE)) == NULL)
  Die("Failed to parse nseq/alen from first line of PHYLIP file %s\n", afp->fname);
      if (! IsInt(s1) || ! IsInt(s2))
  Die("nseq and/or alen not an integer in first line of PHYLIP file %s\n", afp->fname);
      nseq = atoi(s1);
      break;
    }

  msa = MSAAlloc(nseq, 0);
  idx    = 0;
  nblock = 0;
  while ((s = MSAFileGetLine(afp)) != NULL) 
    {
      /* ignore blank lines. nonblank lines start w/ nonblank char */
      if (isspace(*s)) continue;
        /* First block has seq names */
      if (nblock == 0) {
  strncpy(name, s, 10);
  name[10] = '\0';
  GKIStoreKey(msa->index, name);
  msa->sqname[idx] = strdup(name);
  s += 10;    
      }
        /* be careful of trailing whitespace on lines */
      if ((s1 = strtok(s, WHITESPACE)) == NULL)
  Die("Failed to parse sequence at line %d of PHYLIP file %s\n", 
      afp->linenumber, afp->fname);
      msa->sqlen[idx] = strlen(strcat(msa->aseq[idx], s1));

      idx++;
      if (idx == nseq) { idx = 0; nblock++; }
    }
  msa->nseq = nseq;
  MSAVerifyParse(msa);    /* verifies; sets alen, wgt; frees sqlen[] */
  return msa;
}


void
WritePhylip(
  FILE *fp, 
  MSA *msa
){
  int    idx;     /* counter for sequences         */
  int    cpl = 50;    /* 50 seq char per line          */
  char   buf[51];   /* buffer for writing seq        */
  int    pos;

  /* First line has nseq, alen
   */
  fprintf(fp, " %ld  %ld\n", msa->nseq, msa->alen);

  /* Alignment section.
   * PHYLIP is a multiblock format, blocks (optionally) separated
   * by blanks; names only attached to first block. Names are
   * restricted to ten char; we achieve this by simple truncation (!).
   * (Do we need to convert gap characters from our ./- convention?)
   */
  for (pos = 0; pos < msa->alen; pos += cpl)
    {
      if (pos > 0) fprintf(fp, "\n");

      for (idx = 0; idx < msa->nseq; idx++)
  {
    strncpy(buf, msa->aseq[idx] + pos, cpl);
    buf[cpl] = '\0';        
    if (pos > 0) fprintf(fp, "%s\n", buf);
    else         fprintf(fp, "%-10.10s%s\n", msa->sqname[idx], buf);
  }
    }
  return;
}



#ifdef TESTDRIVE_PHYLIP
/*****************************************************************
 * phylip.c test driver:
 * 
 */
int 
main(
  int argc, 
  char **argv
){
  MSAFILE *afp;
  MSA     *msa;
  char    *file;
  
  file = argv[1];

  if ((afp = MSAFileOpen(file, MSAFILE_UNKNOWN, NULL)) == NULL)
    Die("Couldn't open %s\n", file);

  printf("format %d\n", afp->format);

  while ((msa = ReadPhylip(afp)) != NULL)
    {
      WritePhylip(stdout, msa);
      MSAFree(msa); 
    }
  
  MSAFileClose(afp);
  exit(0);
}
/******************************************************************/
#endif /* testdrive_phylip */
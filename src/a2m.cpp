/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* a2m.c
 * 
 * reading/writing A2M (aligned FASTA) files.
 * 
 */
#include "a2m.h"


MSA *
ReadA2M(
  MSAFILE *afp
){
  MSA  *msa;
  char *buf;
  char *name;
  char *desc;
  char *seq;
  int   idx;
  
  if (feof(afp->f)) return NULL;

  name = NULL;
  msa  = MSAAlloc(10, 0);
  idx  = 0;
  while ((buf = MSAFileGetLine(afp)) != NULL) 
    {
      if (*buf == '>') 
  {
    buf++;    /* skip the '>' */
    if ((name = strtok(buf, WHITESPACE)) == NULL)
      Die("Blank name in A2M file %s (line %d)\n", afp->fname, afp->linenumber);
    desc = strtok(buf, "\n");
  
    idx = GKIStoreKey(msa->index, name);
    if (idx >= msa->nseqalloc) MSAExpand(msa);

    msa->sqname[idx] = strdup(name);
    if (desc != NULL) MSASetSeqDescription(msa, idx, desc);
    msa->nseq++;
  } 
      else if (name != NULL) 
  {
    if ((seq = strtok(buf, WHITESPACE)) == NULL) continue; 
    msa->sqlen[idx] = strlen(strcat(msa->aseq[idx], seq));
  }
    } 
  if (name == NULL) { MSAFree(msa); return NULL; }

  MSAVerifyParse(msa);
  return msa;
}



void
WriteA2M(
  FILE *fp, 
  MSA *msa
){
  int  idx;     /* sequence index */
  int  pos;     /* position in sequence */
  char buf[64];     /* buffer for individual lines */
  int  cpl = 60;    /* char per line; must be < 64 unless buf is bigger */

  buf[cpl] = '\0';
  for (idx = 0; idx < msa->nseq; idx++)
    {
      fprintf(fp, ">%s %s\n", 
        msa->sqname[idx],
        (msa->sqdesc != NULL && msa->sqdesc[idx] != NULL) ? msa->sqdesc[idx] : "");
      for (pos = 0; pos < msa->alen; pos+=cpl)
  {
    strncpy(buf, &(msa->aseq[idx][pos]), cpl);
    fprintf(fp, "%s\n", buf);
  }
    }
}

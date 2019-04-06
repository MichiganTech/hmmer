/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* clustal.c
 * Import/export of ClustalV/W multiple sequence alignment
 * formatted files. Derivative of msf.c; MSF is a pretty
 * generic interleaved format.  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "squid.hpp"
#include "msa.hpp"



MSA *
ReadClustal(
  MSAFILE *afp
){
  MSA    *msa;
  char   *s;
  int     sqidx;
  char   *name;
  char   *seq;
  char   *s2;

  if (feof(afp->f)) return NULL;

  /* Skip until we see the CLUSTAL header
   */
  while ((s = MSAFileGetLine(afp)) != NULL)
    {
      if (strncmp(s, "CLUSTAL", 7) == 0 &&
    strstr(s, "multiple sequence alignment") != NULL)
  break;
    }
  if (s == NULL) return NULL;

  msa = MSAAlloc(10, 0);

  /* Now we're in the sequence section.
   * As discussed above, if we haven't seen a sequence name, then we
   * don't include the sequence in the alignment.
   * Watch out for conservation markup lines that contain *.: chars
   */
  while ((s = MSAFileGetLine(afp)) != NULL)
    {
      if ((name = strtok(s, WHITESPACE))  == NULL) continue;
      if ((seq  = strtok(NULL, WHITESPACE)) == NULL) continue;
      s2 = strtok(NULL, "\n");

      /* The test for a conservation markup line
       */
      if (strpbrk(name, ".*:") != NULL && strpbrk(seq, ".*:") != NULL)
  continue;
      if (s2 != NULL)
  Die("Parse failed at line %d, file %s: possibly using spaces as gaps",
      afp->linenumber, afp->fname);
 
      /* It's not blank, and it's not a coord line: must be sequence
       */
      sqidx = MSAGetSeqidx(msa, name, msa->lastidx+1);
      msa->lastidx = sqidx;
      msa->aseq[sqidx] = realloc(msa->aseq[sqidx], strlen(msa->aseq[sqidx]) + strlen(seq) + 1);
      msa->sqlen[sqidx] = strlen(strcat(msa->aseq[sqidx], seq));
    }

  MSAVerifyParse(msa);    /* verifies, and also sets alen and wgt. */
  return msa;
}



void
WriteClustal(
  FILE *fp,
  MSA *msa
){
  int    idx;     /* counter for sequences         */
  int    len;     /* tmp variable for name lengths */
  int    namelen;   /* maximum name length used      */
  int    pos;     /* position counter              */
  char   buf[64];         /* buffer for writing seq        */
  int    cpl = 50;    /* char per line (< 64)          */

        /* calculate max namelen used */
  namelen = 0;
  for (idx = 0; idx < msa->nseq; idx++)
    if ((len = strlen(msa->sqname[idx])) > namelen)
      namelen = len;

  fprintf(fp, "CLUSTAL W(1.5) multiple sequence alignment\n");

  /*****************************************************
   * Write the sequences
   *****************************************************/

  for (pos = 0; pos < msa->alen; pos += cpl)
    {
      fprintf(fp, "\n");  /* Blank line between sequence blocks */
      for (idx = 0; idx < msa->nseq; idx++)
  {
    strncpy(buf, msa->aseq[idx] + pos, cpl);
    buf[cpl] = '\0';
    fprintf(fp, "%*s %s\n", namelen, msa->sqname[idx], buf);
  }
    }

  return;
}



#ifdef TESTDRIVE_CLUSTAL
/*****************************************************************
 * msf.c test driver:
 * cc -DTESTDRIVE_CLUSTAL -g -O2 -Wall -o test clustal.c msa.c gki.c sqerror.c sre_string.c file.c hsregex.c sre_math.c sre_ctype.c -lm
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

  if ((afp = MSAFileOpen(file, MSAFILE_CLUSTAL, NULL)) == NULL)
    Die("Couldn't open %s\n", file);

  while ((msa = ReadClustal(afp)) != NULL)
    {
      WriteClustal(stdout, msa);
      MSAFree(msa);
    }
 
  MSAFileClose(afp);
  exit(0);
}
/******************************************************************/
#endif /* testdrive_clustal */
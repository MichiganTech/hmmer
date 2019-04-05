/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* alphabet.c
 * Configuration of the global symbol alphabet information.
 */

#include "config.h"
//#include "squidconf.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>

#include "structs.h"
#include "funcs.h"
#include "squid.h"
#include "alphabet.h"


void
DetermineAlphabet(
  char **rseqs, 
  int  nseq
){
  int idx;
  int other, nucleic, amino;
  int type;

  /* Autodetection of alphabet type.
   */
  type = hmmNOTSETYET;
  other = nucleic = amino = 0;
  for (idx = 0; idx < nseq; idx++) {
    switch (Seqtype(rseqs[idx])) {
    case kRNA:
      nucleic++;
      break;
    case kDNA:
      nucleic++;
      break;
    case kAmino:
      amino++;
      break;
    case kOtherSeq:
      other++;
      break;
    default:
      Die("No such alphabet type");
    }
  }

  if      (nucleic == nseq) type = hmmNUCLEIC;
  else if (amino   == nseq) type = hmmAMINO;
  else if (nucleic > amino && nucleic > other) {
    Warn("Looks like nucleic acid sequence, hope that's right");
    type = hmmNUCLEIC;
  } else if (amino > nucleic && amino > other) {
    Warn("Looks like amino acid sequence, hope that's right");
    type = hmmAMINO;
  } else Die("Sorry, I can't tell if that's protein or DNA");

  /* Now set up the alphabet.
   */
  SetAlphabet(type);
}


void
SetAlphabet(
  int type
){
  int x;
  pthread_mutex_t  alphabet_lock; /* alphabet is global; must protect to be threadsafe */
  int              rtn;      /* return code from pthreads */

  if ((rtn = pthread_mutex_init(&alphabet_lock, NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));
  if ((rtn = pthread_mutex_lock(&alphabet_lock)) != 0)
    Die("pthread_mutex_lock FAILED: %s\n", strerror(rtn));

  /* Because the alphabet information is global, we must
   * be careful to make this a thread-safe function. The mutex
   * (above) takes care of that. But, indeed, it's also
   * just good sense (and more efficient) to simply never
   * allow resetting the alphabet. If type is Alphabet_type,
   * silently return; else die with an alphabet mismatch
   * warning.
   */
  if (Alphabet_type != hmmNOTSETYET) {
    if (type != Alphabet_type)
      Die("An alphabet type conflict occurred.\nYou probably mixed a DNA seq file with a protein model, or vice versa.");

    if ((rtn = pthread_mutex_unlock(&alphabet_lock)) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
    return;
  }

  switch(type) {
  case hmmAMINO:
    Alphabet_type     = type;
    strcpy(Alphabet, "ACDEFGHIKLMNPQRSTVWYUBZX");
    Alphabet_size     = 20;
    Alphabet_iupac    = 24;
    for (x = 0; x < Alphabet_iupac; x++) {
      memset(Degenerate[x], 0, Alphabet_size);
    }
    for (x = 0; x < Alphabet_size; x++) {
      Degenerate[x][x] = 1;
      DegenCount[x] = 1;
    }
    set_degenerate('U', "S");  /* selenocysteine is treated as serine */
    set_degenerate('B', "ND");
    set_degenerate('Z', "QE");
    set_degenerate('X', "ACDEFGHIKLMNPQRSTVWY");
    break;
  case hmmNUCLEIC:
    Alphabet_type     = type;
    strcpy(Alphabet, "ACGTUNRYMKSWHBVDX");
    Alphabet_size     = 4;
    Alphabet_iupac    = 17;
    for (x = 0; x < Alphabet_iupac; x++) {
      memset(Degenerate[x], 0, Alphabet_size);
    }
    for (x = 0; x < Alphabet_size; x++) {
      Degenerate[x][x] = 1;
      DegenCount[x] = 1;
    }
    set_degenerate('U', "T");
    set_degenerate('N', "ACGT");
    set_degenerate('X', "ACGT");
    set_degenerate('R', "AG");
    set_degenerate('Y', "CT");
    set_degenerate('M', "AC");
    set_degenerate('K', "GT");
    set_degenerate('S', "CG");
    set_degenerate('W', "AT");
    set_degenerate('H', "ACT");
    set_degenerate('B', "CGT");
    set_degenerate('V', "ACG");
    set_degenerate('D', "AGT");
    break;
  default:
    Die("No support for non-nucleic or protein alphabets");
  }

  if ((rtn = pthread_mutex_unlock(&alphabet_lock)) != 0)
    Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
}


unsigned char
SymbolIndex(
  char sym
){
  char *s;
  return ((s = strchr(Alphabet, (char) toupper((int) sym))) == NULL) ?
         Alphabet_iupac-1 : s - Alphabet;
}


unsigned char*
DigitizeSequence(
  char *seq, 
  int L
){
  unsigned char *dsq;
  int i;

  dsq = MallocOrDie (sizeof(unsigned char) * (L+2));
  dsq[0] = dsq[L+1] = (unsigned char) Alphabet_iupac;
  for (i = 1; i <= L; i++)
    dsq[i] = SymbolIndex(seq[i-1]);
  return dsq;
}


char*
DedigitizeSequence(
  unsigned char *dsq, 
  int L
){
  char *seq;
  int i;

  seq = MallocOrDie(sizeof(char) * (L+1));
  for (i = 0; i < L; i++)
    seq[i] = Alphabet[dsq[i+1]];
  seq[L] = '\0';
  return seq;
}


void
DigitizeAlignment(
  MSA *msa, 
  unsigned char ***ret_dsqs
){
  unsigned char **dsq;
  // idx      /* counter for sequences     */
  // apos;      /* position in aligned seq   */

  dsq = MallocOrDie (sizeof(unsigned char *) * msa->nseq);
  for (size_t idx = 0; idx < msa->nseq; idx++) {
    dsq[idx] = MallocOrDie (sizeof(unsigned char) * (msa->alen+2));

    dsq[idx][0] = (unsigned char) Alphabet_iupac; /* sentinel byte at start */

    int dpos = 1;      /* position in digitized seq */
    for (size_t apos = 0; apos < msa->alen; apos++) {
      if (! isgap(msa->aseq[idx][apos]))  /* skip gaps */
        dsq[idx][dpos++] = SymbolIndex(msa->aseq[idx][apos]);
    }
    dsq[idx][dpos] = (unsigned char) Alphabet_iupac; /* sentinel byte at end */
  }
  *ret_dsqs = dsq;
}


void
P7CountSymbol(
  float *counters, 
  unsigned char symidx, 
  float wt
){
  if (symidx < Alphabet_size)
    counters[symidx] += wt;
  else
    for (int x = 0; x < Alphabet_size; x++) {
      if (Degenerate[symidx][x])
        counters[x] += wt / (float) DegenCount[symidx];
    }
}


void
set_degenerate(char iupac, char *syms) {
  DegenCount[strchr(Alphabet,iupac)-Alphabet] = strlen(syms);
  while (*syms) {
    Degenerate[strchr(Alphabet,iupac)-Alphabet]
    [strchr(Alphabet,*syms)-Alphabet] = 1;
    syms++;
  }
}

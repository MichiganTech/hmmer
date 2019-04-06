/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* debug.c
 *
 * Printing out or naming various useful things from HMMER
 * innards.
 */

#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "structs.hpp"
#include "debug.hpp"


char*
Statetype(
  char st
){
  switch (st) {
  case STS:
    return "S";
  case STN:
    return "N";
  case STB:
    return "B";
  case STM:
    return "M";
  case STD:
    return "D";
  case STI:
    return "I";
  case STE:
    return "E";
  case STJ:
    return "J";
  case STC:
    return "C";
  case STT:
    return "T";
  default:
    return "BOGUS";
  }
}


char*
AlphabetType2String(
  int type
){
  switch (type) {
  case hmmAMINO:
    return "protein";
  case hmmNUCLEIC:
    return "nucleic acid";
  case hmmNOTSETYET:
    return "unknown";
  default:
    return "BOGUS";
  }
}


void
P7PrintTrace(
  FILE *fp,
  struct p7trace_s *tr,
  struct plan7_s *hmm,
  unsigned char *dsq
){
  int          tpos;    /* counter for trace position */

  if (tr == NULL) {
    fprintf(fp, " [ trace is NULL ]\n");
    return;
  }

  if (hmm == NULL) {
    fprintf(fp, "st  node   rpos  - traceback len %d\n", tr->tlen);
    fprintf(fp, "--  ---- ------\n");
    for (tpos = 0; tpos < tr->tlen; tpos++) {
      fprintf(fp, "%1s  %4d %6d\n",
              Statetype(tr->statetype[tpos]),
              tr->nodeidx[tpos],
              tr->pos[tpos]);
    }
  } else {
    if (!(hmm->flags & PLAN7_HASBITS))
      Die("oi, you can't print scores from that hmm, it's not ready.");

    unsigned int sym;
    int          sc;
    sc = 0;
    fprintf(fp, "st  node   rpos  transit emission - traceback len %d\n", tr->tlen);
    fprintf(fp, "--  ---- ------  ------- --------\n");
    for (tpos = 0; tpos < tr->tlen; tpos++) {
      if (dsq != NULL) sym = dsq[tr->pos[tpos]];

      fprintf(fp, "%1s  %4d %6d  %7d",
              Statetype(tr->statetype[tpos]),
              tr->nodeidx[tpos],
              tr->pos[tpos],
              (tpos < tr->tlen-1) ?
              TransitionScoreLookup(hmm, tr->statetype[tpos], tr->nodeidx[tpos],
                                    tr->statetype[tpos+1], tr->nodeidx[tpos+1]) : 0);

      if (tpos < tr->tlen-1)
        sc += TransitionScoreLookup(hmm, tr->statetype[tpos], tr->nodeidx[tpos],
                                    tr->statetype[tpos+1], tr->nodeidx[tpos+1]);

      if (dsq != NULL) {
        if (tr->statetype[tpos] == STM) {
          fprintf(fp, " %8d %c", hmm->msc[sym][tr->nodeidx[tpos]], Alphabet[sym]);
          sc += hmm->msc[sym][tr->nodeidx[tpos]];
        } else if (tr->statetype[tpos] == STI) {
          fprintf(fp, " %8d %c", hmm->isc[sym][tr->nodeidx[tpos]],
                  (char) tolower((int) Alphabet[sym]));
          sc += hmm->isc[sym][tr->nodeidx[tpos]];
        } else if ((tr->statetype[tpos] == STN && tr->statetype[tpos-1] == STN) ||
                   (tr->statetype[tpos] == STC && tr->statetype[tpos-1] == STC) ||
                   (tr->statetype[tpos] == STJ && tr->statetype[tpos-1] == STJ)) {
          fprintf(fp, " %8d %c", 0, (char) tolower((int) Alphabet[sym]));
        }
      } else {
        fprintf(fp, " %8s %c", "-", '-');
      }

      fputs("\n", fp);
    }
    fprintf(fp, "                 ------- --------\n");
    fprintf(fp, "           total: %6d\n\n", sc);
  }
}


bool
TraceVerify(
  struct p7trace_s *tr,
  int M,
  int N
){
  int tpos;      // position in trace
  int k;      // current position in HMM nodes 1..M
  int i;      // current position in seq 1..N
  int nn, nc, nj;    // number of STN's, STC's, STJ's seen
  int nm;      // number of STM's seen

  // Basic checks on ends.

  if (tr->statetype[0] != STS)          return false;
  if (tr->statetype[1] != STN)          return false;
  if (tr->statetype[tr->tlen-2] != STC) return false;
  if (tr->statetype[tr->tlen-1] != STT) return false;
  if (tr->pos[1] != 0)                  return false;

  // Check for consistency throughout trace

  k = i = nn = nc = nj = nm = 0;
  for (tpos = 0; tpos < tr->tlen; tpos++) {
    switch (tr->statetype[tpos]) {
    case STS:
      if (tr->nodeidx[tpos] != 0) return false;
      if (tr->pos[tpos]     != 0) return false;
      if (k != 0)                 return false;
      if (i != 0)                 return false;
      if (tpos != 0)              return false;
      break;

    case STN:      // first N doesn't emit.
      if (tr->nodeidx[tpos] != 0) return false;
      if (k != 0)                 return false;
      if (nn > 0) {
        if (tr->pos[tpos] != i+1) return false;
        i++;
      } else {
        if (tr->pos[tpos] != 0) return false;
        if (i != 0)             return false;
      }
      nn++;
      break;

    case STB:
      if (tr->nodeidx[tpos] != 0) return false;
      if (tr->pos[tpos]     != 0) return false;
      nm = 0;
      break;

    case STM:      // can enter anywhere on first M
      if (tr->pos[tpos] != i+1) return false;
      if (tr->nodeidx[tpos] < 1 || tr->nodeidx[tpos] > M) return false;
      i++;
      if (nm == 0)  k = tr->nodeidx[tpos];
      else {
        if (tr->nodeidx[tpos] != k+1) return false;
        k++;
      }
      nm++;
      break;

    case STI:
      if (tr->pos[tpos] != i+1)   return false;
      if (tr->nodeidx[tpos] != k) return false;
      if (tr->nodeidx[tpos] < 1 || tr->nodeidx[tpos] > M-1) return false;
      if (k >= M)                 return false;
      i++;
      break;

    case STD:
      if (tr->pos[tpos] != 0)       return false;
      if (tr->nodeidx[tpos] != k+1) return false;
      if (tr->nodeidx[tpos] < 1 || tr->nodeidx[tpos] > M) return false;
      k++;
      break;

    case STE:
      if (tr->nodeidx[tpos] != 0) return false;
      if (tr->pos[tpos]     != 0) return false;
      nj = 0;
      break;

    case STJ:
      if (tr->nodeidx[tpos] != 0) return false;
      if (nj > 0) {
        if (tr->pos[tpos] != i+1) return false;
        i++;
      } else if (tr->pos[tpos] != 0) return false;
      nj++;
      break;

    case STC:
      if (tr->nodeidx[tpos] != 0) return false;
      if (nc > 0) {
        if (tr->pos[tpos] != i+1) return false;
        i++;
      } else if (tr->pos[tpos] != 0)  return false;
      nc++;
      break;

    case STT:
      if (tpos != tr->tlen - 1)   return false;
      if (tr->nodeidx[tpos] != 0) return false;
      if (tr->pos[tpos]     != 0) return false;
      if (i != N)                 return false;
      break;

    case STBOGUS:
    default:
      return false;
    }  // end switch over statetypes
  } // end loop over trace positions

  return true;
}


bool
TraceCompare(
  struct p7trace_s *t1,
  struct p7trace_s *t2
){
  int tpos;

  if (t1->tlen != t2->tlen) return false;

  for (tpos = 0; tpos < t1->tlen; tpos++) {
    if (t1->statetype[tpos] != t2->statetype[tpos]) return false;
    if (t1->nodeidx[tpos]   != t2->nodeidx[tpos])   return false;
    if (t1->pos[tpos]       != t2->pos[tpos])       return false;
  }
  return true;
}

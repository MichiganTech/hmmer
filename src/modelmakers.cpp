/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* modelmakers.c
 *
 * Construction of models from multiple alignments. Three versions:
 *    Handmodelmaker() -- use #=RF annotation to indicate match columns
 *    Fastmodelmaker() -- Krogh/Haussler heuristic
 *    Maxmodelmaker()  -- MAP model construction algorithm (Eddy,
 *                          unpublished)
 *
 * The meat of the model construction code is in matassign2hmm().
 * The three model construction strategies simply label which columns
 * are supposed to be match states, and then hand this info to
 * matassign2hmm().
 *
 * Two wrinkles to watch for:
 * 1) The alignment is assumed to contain sequence fragments. Look in
 *    fake_tracebacks() for how internal entry/exit points are handled.
 * 2) Plan7 disallows DI and ID transitions, but an alignment may
 *    imply these. Look in trace_doctor() for how DI and ID transitions
 *    are removed.
 */

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include "structs.hpp"
#include "vectorops.hpp"
#include "modelmakers.hpp"



/* flags used for matassign[] arrays --
 *   assignment of aligned columns to match/insert states
 */
#define ASSIGN_MATCH      (1<<0)
#define FIRST_MATCH       (1<<1)
#define LAST_MATCH        (1<<2)
#define ASSIGN_INSERT     (1<<3)
#define EXTERNAL_INSERT_N (1<<4)
#define EXTERNAL_INSERT_C (1<<5)


void
P7Handmodelmaker(
  MSA *msa,
  unsigned char **dsq,
  char *isfrag,
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr
){
  int     *matassign;           /* MAT state assignments if 1; 1..alen */
  // apos                /* counter for aligned columns         */

  /* Make sure we have all the info about the alignment that we need */
  if (msa->rf == NULL)
    Die("Alignment must have RF annotation to hand-build an HMM");

  /* Allocation */
  matassign = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));

  /* Determine match assignment from optional annotation
   */
  matassign[0] = 0;
  for (size_t apos = 0; apos < msa->alen; apos++) {
    matassign[apos+1] = 0;
    if (!isgap(msa->rf[apos]))
      matassign[apos+1] |= ASSIGN_MATCH;
    else
      matassign[apos+1] |= ASSIGN_INSERT;
  }

  /* Hand matassign off for remainder of model construction
   */
  /*   print_matassign(matassign, msa->alen); */
  matassign2hmm(msa, dsq, isfrag, matassign, ret_hmm, ret_tr);

  free(matassign);
  return;
}


void
P7Fastmodelmaker(
  MSA *msa,
  unsigned char **dsq,
  char *isfrag,
  float symfrac,
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr
){
  int     *matassign;           /* MAT state assignments if 1; 1..alen */
  // idx                 /* counter over sequences              */
  // apos                /* counter for aligned columns         */
  float   *r;            /* weighted frac of gaps in column     */
  float    totwgt;          /* total non-fragment seq weight       */
  float    maxR;               /* maximum r_i                         */
  int      incfrags;    /* true to include candidate frags  */

  /* Allocations: matassign is 1..alen array of bit flags;
   *              gapfrac is 1..alen array of fractions 0<=gapfrac[i]<=1
   */
  matassign = MallocOrDie (sizeof(int)   * (msa->alen+1));
  r         = MallocOrDie (sizeof(float) * (msa->alen+1));

  /* Determine total non-frag weight, just once.
   */
  incfrags = false;
  totwgt   = 0.;
  for (size_t idx = 0; idx < msa->nseq; idx++)
    if (! isfrag[idx]) totwgt += msa->wgt[idx];

  /* Fallback position, if we don't have any non-candidate frags:
   * count all seqs.
   */
  if (totwgt == 0.) { /* yes, this fp compare is safe */
    totwgt = FSum(msa->wgt, msa->nseq);
    incfrags = true;
  }

  /* Determine weighted sym freq in each column; keep track of max.
   * Mind the off-by-one (r is 1..alen, msa->aseq is 0..alen-1)
   */
  for (size_t apos = 0; apos < msa->alen; apos++) {
    r[apos+1] = 0.;
    for (size_t idx = 0; idx < msa->nseq; idx++)
      if ((incfrags || ! isfrag[idx])
          && ! isgap(msa->aseq[idx][apos]))
        r[apos+1] += msa->wgt[idx];
    r[apos+1] /= totwgt;
  }
  maxR = FMax(r+1, msa->alen);

  /* Determine match assignment. (Both matassign and r are 1..alen)
   */
  matassign[0] = 0;
  for (size_t apos = 1; apos <= msa->alen; apos++) {
    matassign[apos] = 0;
    if (r[apos] >= symfrac * maxR) matassign[apos] |= ASSIGN_MATCH;
    else                           matassign[apos] |= ASSIGN_INSERT;
  }

  /* Once we have matassign calculated, all modelmakers behave
   * the same; matassign2hmm() does this stuff (traceback construction,
   * trace counting) and sets up ret_hmm and ret_tr.
   */
  matassign2hmm(msa, dsq, isfrag, matassign, ret_hmm, ret_tr);

  free(matassign);
  free(r);
  return;
}


void
matassign2hmm(
  MSA *msa,
  unsigned char **dsq,
  char *isfrag,
  int *matassign,
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr
){
  struct plan7_s    *hmm;       /* RETURN: new hmm                     */
  struct p7trace_s **tr;        /* fake tracebacks for each seq        */
  int      M;                   /* length of new model in match states */
  // idx                 /* counter over sequences              */
  size_t apos;                /* counter for aligned columns         */

  /* how many match states in the HMM? */
  M = 0;
  for (apos = 1; apos <= msa->alen; apos++) {
    if (matassign[apos] & ASSIGN_MATCH)
      M++;
  }

  if (M == 0)
    Die("No conserved consensus columns found; aborting construction!\n\
This is an unusual situation. Reexamine your sequence alignment. It is\n\
probably unusually full of gaps, or lots of sequence fragments. You may be\n\
able to force HMMER to model it; see the --fast (and --gapmax), or --hand\n\
options to hmmbuild.");

  /* delimit N-terminal tail */
  for (apos=1; apos <= msa->alen && matassign[apos] & ASSIGN_INSERT; apos++)
    matassign[apos] |= EXTERNAL_INSERT_N;
  if (apos <= msa->alen) matassign[apos] |= FIRST_MATCH;

  /* delimit C-terminal tail */
  for (apos=msa->alen; matassign[apos] & ASSIGN_INSERT && apos > 0; apos--)
    matassign[apos] |= EXTERNAL_INSERT_C;
  if (apos > 0) matassign[apos] |= LAST_MATCH;

  /* print_matassign(matassign, msa->alen);  */

  /* make fake tracebacks for each seq */
  fake_tracebacks(msa->aseq, msa->nseq, msa->alen, isfrag, matassign, &tr);
  /* build model from tracebacks */
  hmm = AllocPlan7(M);
  ZeroPlan7(hmm);
  for (size_t idx = 0; idx < msa->nseq; idx++) {
    /* P7PrintTrace(stdout, tr[idx], NULL, NULL);   */
    P7TraceCount(hmm, dsq[idx], msa->wgt[idx], tr[idx]);
  }
  /* annotate new model */
  annotate_model(hmm, matassign, msa);

  /* Set #=RF line of alignment to reflect our assignment
   * of match, delete. matassign is valid from 1..alen and is off
   * by one from msa->rf.
   */
  if (msa->rf != NULL) free(msa->rf);
  msa->rf = (char *) MallocOrDie (sizeof(char) * (msa->alen + 1));
  for (apos = 0; apos < msa->alen; apos++)
    msa->rf[apos] = (matassign[apos+1] & ASSIGN_MATCH) ? 'x' : '.';
  msa->rf[msa->alen] = '\0';

  /* Cleanup and return. */
  if (ret_tr != NULL) *ret_tr = tr;
  else   {
    for (size_t idx = 0; idx < msa->nseq; idx++) P7FreeTrace(tr[idx]);
    free(tr);
  }
  if (ret_hmm != NULL) *ret_hmm = hmm;
  else FreePlan7(hmm);
  return;
}


void
fake_tracebacks(
  char **aseq,
  int nseq,
  int alen,
  char *isfrag,
  int *matassign,
  struct p7trace_s ***ret_tr
){
  struct p7trace_s **tr;
  int  idx;                     /* counter over sequences          */
  int  apos;                    /* position in alignment columns   */

  tr = (struct p7trace_s **) MallocOrDie (sizeof(struct p7trace_s *) * nseq);

  for (idx = 0; idx < nseq; idx++) {
    P7AllocTrace(alen+6, &tr[idx]); /* allow room for S,N,B,E,C,T */

    /* all traces start with S state... */
    tr[idx]->statetype[0] = STS;
    tr[idx]->nodeidx[0]   = 0;
    tr[idx]->pos[0]       = 0;
    /* ...and transit to N state; N-term tail
       is emitted on N->N transitions */
    tr[idx]->statetype[1] = STN;
    tr[idx]->nodeidx[1]   = 0;
    tr[idx]->pos[1]       = 0;

    int  i; //position in raw sequence (1..L)
    int  k; //position in HMM
    int  tpos; //position in traceback
    i = 1;
    k = 0;
    tpos = 2;
    for (apos = 0; apos < alen; apos++) {
      tr[idx]->statetype[tpos] = STBOGUS; /* bogus, deliberately, to debug */

      if (matassign[apos+1] & FIRST_MATCH) {
        /* BEGIN */
        tr[idx]->statetype[tpos] = STB;
        tr[idx]->nodeidx[tpos]   = 0;
        tr[idx]->pos[tpos]       = 0;
        tpos++;
      }

      if (matassign[apos+1] & ASSIGN_MATCH && ! isgap(aseq[idx][apos])) {
        /* MATCH */
        k++;    /* move to next model pos */
        tr[idx]->statetype[tpos] = STM;
        tr[idx]->nodeidx[tpos]   = k;
        tr[idx]->pos[tpos]       = i;
        i++;
        tpos++;
      } else if (matassign[apos+1] & ASSIGN_MATCH) {
        /* DELETE */
        /* being careful about S/W transitions.
                     * Count B->D1 transitions only if we're not
                     * a candidate fragment seq.
         */
        k++;    /* *always* move on model when ASSIGN_MATCH */
        if (tr[idx]->statetype[tpos-1] != STB || ! isfrag[idx]) {
          tr[idx]->statetype[tpos] = STD;
          tr[idx]->nodeidx[tpos]   = k;
          tr[idx]->pos[tpos]       = 0;
          tpos++;
        }
      } else if (matassign[apos+1] & EXTERNAL_INSERT_N &&
                 ! isgap(aseq[idx][apos])) {
        /* N-TERMINAL TAIL */
        tr[idx]->statetype[tpos] = STN;
        tr[idx]->nodeidx[tpos]   = 0;
        tr[idx]->pos[tpos]       = i;
        i++;
        tpos++;
      } else if (matassign[apos+1] & EXTERNAL_INSERT_C &&
                 ! isgap(aseq[idx][apos])) {
        /* C-TERMINAL TAIL */
        tr[idx]->statetype[tpos] = STC;
        tr[idx]->nodeidx[tpos]   = 0;
        tr[idx]->pos[tpos]       = i;
        i++;
        tpos++;
      } else if (! isgap(aseq[idx][apos])) {
        /* INSERT */
        tr[idx]->statetype[tpos] = STI;
        tr[idx]->nodeidx[tpos]   = k;
        tr[idx]->pos[tpos]       = i;
        i++;
        tpos++;
      }

      if (matassign[apos+1] & LAST_MATCH) {
        /* END */
        /* be careful about S/W transitions; if we're
               * a candidate sequence fragment, we need to roll
         * back over some D's because there's no D->E transition
               * for fragments. k==M right now; don't change it.
         */
        if (isfrag[idx])
          while (tr[idx]->statetype[tpos-1] == STD)
            tpos--;
        tr[idx]->statetype[tpos] = STE;
        tr[idx]->nodeidx[tpos]   = 0;
        tr[idx]->pos[tpos]       = 0;
        tpos++;
        /* and then transit E->C;
           alignments that use J are undefined;
           C-term tail is emitted on C->C transitions */
        tr[idx]->statetype[tpos] = STC;
        tr[idx]->nodeidx[tpos]   = 0;
        tr[idx]->pos[tpos]       = 0;
        tpos++;
      }
    }
    /* all traces end with T state */
    tr[idx]->statetype[tpos] = STT;
    tr[idx]->nodeidx[tpos]   = 0;
    tr[idx]->pos[tpos]       = 0;
    tr[idx]->tlen = ++tpos;
    /* deal with DI, ID transitions */
    /* k == M here */
    trace_doctor(tr[idx], k, NULL, NULL);

  }    /* end for sequence # idx */

  *ret_tr = tr;
  return;
}


void
trace_doctor(
  struct p7trace_s *tr,
  int mlen,
  int *ret_ndi,
  int *ret_nid
){
  int opos;      /* position in old trace                 */
  int npos;      /* position in new trace (<= opos)       */
  int ndi, nid;      /* number of DI, ID transitions doctored */

  /* overwrite the trace from left to right */
  ndi  = nid  = 0;
  opos = npos = 0;
  while (opos < tr->tlen) {
    /* fix implied D->I transitions; D transforms to M, I pulled in */
    if (tr->statetype[opos] == STD && tr->statetype[opos+1] == STI) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]; /* D transforms to M      */
      tr->pos[npos]       = tr->pos[opos+1];   /* insert char moves back */
      opos += 2;
      npos += 1;
      ndi++;
    } /* fix implied I->D transitions; D transforms to M, I is pushed in */
    else if (tr->statetype[opos]== STI && tr->statetype[opos+1]== STD) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos+1];/* D transforms to M    */
      tr->pos[npos]       = tr->pos[opos];      /* insert char moves up */
      opos += 2;
      npos += 1;
      nid++;
    } /* fix implied B->I transitions; pull I back to its M */
    else if (tr->statetype[opos]== STI && tr->statetype[opos-1]== STB) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]; /* offending I transforms to M */
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    } /* fix implied I->E transitions; push I to next M */
    else if (tr->statetype[opos]== STI && tr->statetype[opos+1]== STE) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]+1;/* offending I transforms to M */
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    } /* rare: N-N-B-E becomes N-B-M_1-E (swap B,N) */
    else if (tr->statetype[opos]==STB && tr->statetype[opos+1]==STE
             && tr->statetype[opos-1]==STN && tr->pos[opos-1] > 0) {
      tr->statetype[npos]   = STM;
      tr->nodeidx[npos]     = 1;
      tr->pos[npos]         = tr->pos[opos-1];
      tr->statetype[npos-1] = STB;
      tr->nodeidx[npos-1]   = 0;
      tr->pos[npos-1]       = 0;
      opos++;
      npos++;
    } /* rare: B-E-C-C-x becomes B-M_M-E-C-x (swap E,C) */
    else if (tr->statetype[opos]==STE && tr->statetype[opos-1]==STB
             && tr->statetype[opos+1]==STC
             && tr->statetype[opos+2]==STC) {
      tr->statetype[npos]   = STM;
      tr->nodeidx[npos]     = mlen;
      tr->pos[npos]         = tr->pos[opos+2];
      tr->statetype[npos+1] = STE;
      tr->nodeidx[npos+1]   = 0;
      tr->pos[npos+1]       = 0;
      tr->statetype[npos+2] = STC; /* first C must be a nonemitter  */
      tr->nodeidx[npos+2]   = 0;
      tr->pos[npos+2]       = 0;
      opos+=3;
      npos+=3;
    } /* everything else is just copied */
    else {
      tr->statetype[npos] = tr->statetype[opos];
      tr->nodeidx[npos]   = tr->nodeidx[opos];
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    }
  }
  tr->tlen = npos;

  if (ret_ndi != NULL) *ret_ndi = ndi;
  if (ret_nid != NULL) *ret_nid = nid;
  return;
}


void
annotate_model(
  struct plan7_s *hmm,
  int *matassign,
  MSA *msa
){
  // apos      /* position in matassign, 1.alen  */
  size_t k;      /* position in model, 1.M         */
  char *pri;      /* X-PRM, X-PRI, X-PRT annotation */

  /* Transfer reference coord annotation from alignment,
   * if available
   */
  if (msa->rf != NULL) {
    hmm->rf[0] = ' ';
    for (size_t apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH) /* ainfo is off by one from HMM */
        hmm->rf[k++] = (msa->rf[apos-1] == ' ') ? '.' : msa->rf[apos-1];
    hmm->rf[k] = '\0';
    hmm->flags |= PLAN7_RF;
  }

  /* Transfer consensus structure annotation from alignment,
   * if available
   */
  if (msa->ss_cons != NULL) {
    hmm->cs[0] = ' ';
    for (size_t apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
        hmm->cs[k++] = (msa->ss_cons[apos-1] == ' ') ? '.' : msa->ss_cons[apos-1];
    hmm->cs[k] = '\0';
    hmm->flags |= PLAN7_CS;
  }

  /* Transfer surface accessibility annotation from alignment,
   * if available
   */
  if (msa->sa_cons != NULL) {
    hmm->ca[0] = ' ';
    for (size_t apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
        hmm->ca[k++] = (msa->sa_cons[apos-1] == ' ') ? '.' : msa->sa_cons[apos-1];
    hmm->ca[k] = '\0';
    hmm->flags |= PLAN7_CA;
  }

  /* Store the alignment map
   */
  for (size_t apos = k = 1; apos <= msa->alen; apos++)
    if (matassign[apos] & ASSIGN_MATCH)
      hmm->map[k++] = apos;
  hmm->flags |= PLAN7_MAP;

  /* Translate and transfer X-PRM annotation.
   * 0-9,[a-zA-Z] are legal; translate as 0-9,10-35 into hmm->mpri.
   * Any other char is translated as -1, and this will be interpreted
   * as a flag that means "unknown", e.g. use the normal mixture Dirichlet
   * procedure for this column.
   */
  if ((pri = MSAGetGC(msa, "X-PRM")) != NULL) {
    hmm->mpri = MallocOrDie(sizeof(int) * (hmm->M+1));
    for (size_t apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH) {
        if      (isdigit((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - '0';
        else if (islower((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - 'a' + 10;
        else if (isupper((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - 'A' + 10;
        else hmm->mpri[k] = -1;
        k++;
      }
  }
  /* And again for X-PRI annotation on insert priors:
   */
  if ((pri = MSAGetGC(msa, "X-PRI")) != NULL) {
    hmm->ipri = MallocOrDie(sizeof(int) * (hmm->M+1));
    for (size_t apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH) {
        if      (isdigit((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - '0';
        else if (islower((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - 'a' + 10;
        else if (isupper((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - 'A' + 10;
        else hmm->ipri[k] = -1;
        k++;
      }
  }
  /* And one last time for X-PRT annotation on transition priors:
   */
  if ((pri = MSAGetGC(msa, "X-PRT")) != NULL) {
    hmm->tpri = MallocOrDie(sizeof(int) * (hmm->M+1));
    for (size_t apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH) {
        if      (isdigit((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - '0';
        else if (islower((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - 'a' + 10;
        else if (isupper((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - 'A' + 10;
        else hmm->tpri[k] = -1;
        k++;
      }
  }
}
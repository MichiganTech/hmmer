/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/


/* plan7.c
 *
 * Support for Plan 7 HMM data structure, plan7_s.
 */

//#include "squidconf.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "config.hpp"
#include "squid.hpp"
#include "structs.hpp"
#include "vectorops.hpp"



struct plan7_s *
AllocPlan7(
  int M
){
  struct plan7_s *hmm;

  hmm = AllocPlan7Shell();
  AllocPlan7Body(hmm, M);
  return hmm;
}


struct plan7_s *
AllocPlan7Shell(
){
  struct plan7_s *hmm;

  hmm    = (struct plan7_s *) MallocOrDie (sizeof(struct plan7_s));
  hmm->M = 0;

  hmm->name     = NULL;
  hmm->acc      = NULL;
  hmm->desc     = NULL;
  hmm->rf       = NULL;
  hmm->cs       = NULL;
  hmm->ca       = NULL;
  hmm->comlog   = NULL;
  hmm->nseq     = 0;
  hmm->ctime    = NULL;
  hmm->map      = NULL;
  hmm->checksum = 0;

  hmm->tpri = NULL;
  hmm->mpri = NULL;
  hmm->ipri = NULL;

  hmm->ga1 = hmm->ga2 = 0.0;
  hmm->tc1 = hmm->tc2 = 0.0;
  hmm->nc1 = hmm->nc2 = 0.0;

  hmm->t      = NULL;
  hmm->mat    = NULL;
  hmm->ins    = NULL;

  hmm->tsc     = hmm->msc     = hmm->isc     = NULL;
  hmm->tsc_mem = hmm->msc_mem = NULL;

  hmm->begin  = NULL;
  hmm->end    = NULL;

  hmm->bsc = hmm->bsc_mem = NULL;
  hmm->esc = hmm->esc_mem = NULL;

  /* DNA translation is not enabled by default */
  hmm->dnam   = NULL;
  hmm->dnai   = NULL;
  hmm->dna2   = -INFTY;
  hmm->dna4   = -INFTY;
  /* statistical parameters set to innocuous empty values */
  hmm->mu     = 0.;
  hmm->lambda = 0.;

  hmm->flags = 0;
  return hmm;
}


void
FreePlan7(
  struct plan7_s *hmm
){
  if (hmm->name    != NULL) free(hmm->name);
  if (hmm->acc     != NULL) free(hmm->acc);
  if (hmm->desc    != NULL) free(hmm->desc);
  if (hmm->rf      != NULL) free(hmm->rf);
  if (hmm->cs      != NULL) free(hmm->cs);
  if (hmm->ca      != NULL) free(hmm->ca);
  if (hmm->comlog  != NULL) free(hmm->comlog);
  if (hmm->ctime   != NULL) free(hmm->ctime);
  if (hmm->map     != NULL) free(hmm->map);
  if (hmm->tpri    != NULL) free(hmm->tpri);
  if (hmm->mpri    != NULL) free(hmm->mpri);
  if (hmm->ipri    != NULL) free(hmm->ipri);
  if (hmm->bsc_mem != NULL) free(hmm->bsc_mem);
  if (hmm->begin   != NULL) free(hmm->begin);
  if (hmm->esc_mem != NULL) free(hmm->esc_mem);
  if (hmm->end     != NULL) free(hmm->end);
  if (hmm->msc_mem != NULL) free(hmm->msc_mem);
  if (hmm->isc_mem != NULL) free(hmm->isc_mem);
  if (hmm->tsc_mem != NULL) free(hmm->tsc_mem);
  if (hmm->mat     != NULL) free(hmm->mat[0]);
  if (hmm->ins     != NULL) free(hmm->ins[0]);
  if (hmm->t       != NULL) free(hmm->t[0]);
  if (hmm->msc     != NULL) free(hmm->msc);
  if (hmm->isc     != NULL) free(hmm->isc);
  if (hmm->tsc     != NULL) free(hmm->tsc);
  if (hmm->mat     != NULL) free(hmm->mat);
  if (hmm->ins     != NULL) free(hmm->ins);
  if (hmm->t       != NULL) free(hmm->t);
  if (hmm->dnam    != NULL) free(hmm->dnam);
  if (hmm->dnai    != NULL) free(hmm->dnai);
  free(hmm);
}


void
ZeroPlan7(
  struct plan7_s *hmm
){
  int k;
  for (k = 1; k < hmm->M; k++) {
    FSet(hmm->t[k], 7, 0.);
    FSet(hmm->mat[k], Alphabet_size, 0.);
    FSet(hmm->ins[k], Alphabet_size, 0.);
  }
  FSet(hmm->mat[hmm->M], Alphabet_size, 0.);
  hmm->tbd1 = 0.;
  FSet(hmm->begin+1, hmm->M, 0.);
  FSet(hmm->end+1, hmm->M, 0.);
  for (k = 0; k < 4; k++)
    FSet(hmm->xt[k], 2, 0.);
  hmm->flags &= ~PLAN7_HASBITS;  /* invalidates scores */
  hmm->flags &= ~PLAN7_HASPROB;  /* invalidates probabilities */
}


void
Plan7SetName(
  struct plan7_s *hmm,
  char *name
){
  if (hmm->name != NULL) free(hmm->name);
  hmm->name = strdup(name);
  StringChop(hmm->name);
}


void
Plan7SetAccession(
  struct plan7_s *hmm,
  char *acc
){
  if (hmm->acc != NULL) free(hmm->acc);
  hmm->acc = strdup(acc);
  StringChop(hmm->acc);
  hmm->flags |= PLAN7_ACC;
}


void
Plan7SetDescription(
  struct plan7_s *hmm,
  char *desc
){
  if (hmm->desc != NULL) free(hmm->desc);
  hmm->desc = strdup(desc);
  StringChop(hmm->desc);
  hmm->flags |= PLAN7_DESC;
}


void
Plan7ComlogAppend(
  struct plan7_s *hmm,
  int argc,
  char **argv
){
  int len;
  int i;

  /* figure out length of command line, w/ spaces and \n */
  len = argc;
  for (i = 0; i < argc; i++)
    len += strlen(argv[i]);

  /* allocate */
  if (hmm->comlog != NULL) {
    len += strlen(hmm->comlog);
    hmm->comlog = ReallocOrDie(hmm->comlog, sizeof(char)* (len+1));
  } else {
    hmm->comlog = MallocOrDie(sizeof(char)* (len+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
  }

  /* append */
  strcat(hmm->comlog, "\n");
  for (i = 0; i < argc; i++) {
    strcat(hmm->comlog, argv[i]);
    if (i < argc-1) strcat(hmm->comlog, " ");
  }
}


void
Plan7SetCtime(
  struct plan7_s *hmm
){
  time_t date = time(NULL);
  if (hmm->ctime != NULL) free(hmm->ctime);
  hmm->ctime = strdup(ctime(&date));
  StringChop(hmm->ctime);
}


void
Plan7SetNullModel(
  struct plan7_s *hmm,
  float null[MAXABET],
  float p1
){
  int x;
  for (x = 0; x < Alphabet_size; x++)
    hmm->null[x] = null[x];
  hmm->p1 = p1;
}


void
P7Logoddsify(
  struct plan7_s *hmm,
  int viterbi_mode
){
  int k;      /* counter for model position */
  int x;      /* counter for symbols        */
  float accum;

  if (hmm->flags & PLAN7_HASBITS) return;

  /* Symbol emission scores
   */
  for (k = 1; k <= hmm->M; k++) {
    /* match/insert emissions in main model */
    for (x = 0; x < Alphabet_size; x++) {
      hmm->msc[x][k] = Prob2Score(hmm->mat[k][x], hmm->null[x]);
      if (k < hmm->M)
        hmm->isc[x][k] =  Prob2Score(hmm->ins[k][x], hmm->null[x]);
    }
    /* degenerate match/insert emissions */
    for (x = Alphabet_size; x < Alphabet_iupac; x++) {
      hmm->msc[x][k] = DegenerateSymbolScore(hmm->mat[k], hmm->null, x);
      if (k < hmm->M)
        hmm->isc[x][k] = DegenerateSymbolScore(hmm->ins[k], hmm->null, x);
    }
  }

  /* State transitions.
   *
   * A note on "folding" of D_1 and D_M.
   * These two delete states are folded out of search form models
   * in order to prevent null cycles in the dynamic programming
   * algorithms (see code below). However, we use their log transitions
   * when we save the model! So the following log transition probs
   * are used *only* in save files, *never* in search algorithms:
   *    log (tbd1), D1 -> M2, D1 -> D2
   *    Mm-1 -> Dm, Dm-1 -> Dm
   *
   * In a search algorithm, these have to be interpreted as -INFTY
   * because their contributions are folded into bsc[] and esc[]
   * entry/exit scores. They can't be set to -INFTY here because
   * we need them in save files.
   */
  for (k = 1; k < hmm->M; k++) {
    hmm->tsc[TMM][k] = Prob2Score(hmm->t[k][TMM], hmm->p1);
    hmm->tsc[TMI][k] = Prob2Score(hmm->t[k][TMI], hmm->p1);
    hmm->tsc[TMD][k] = Prob2Score(hmm->t[k][TMD], 1.0);
    hmm->tsc[TIM][k] = Prob2Score(hmm->t[k][TIM], hmm->p1);
    hmm->tsc[TII][k] = Prob2Score(hmm->t[k][TII], hmm->p1);
    hmm->tsc[TDM][k] = Prob2Score(hmm->t[k][TDM], hmm->p1);
    hmm->tsc[TDD][k] = Prob2Score(hmm->t[k][TDD], 1.0);
  }

  /* B->M entry transitions. Note how D_1 is folded out.
   * M1 is just B->M1
   * M2 is sum (or max) of B->M2 and B->D1->M2
   * M_k is sum (or max) of B->M_k and B->D1...D_k-1->M_k
   * These have to be done in log space, else you'll get
   * underflow errors; and we also have to watch for log(0).
   * A little sloppier than it probably has to be; historically,
   * doing in this in log space was in response to a bug report.
   */
  accum = hmm->tbd1 > 0.0 ? log(hmm->tbd1) : -9999.;
  for (k = 1; k <= hmm->M; k++) {
    float tbm;
    tbm = hmm->begin[k] > 0. ? log(hmm->begin[k]) : -9999.;  /* B->M_k part */

    /* B->D1...D_k-1->M_k part we get from accum*/
    if (k > 1 && accum > -9999.) {
      if (hmm->t[k-1][TDM] > 0.0) {
        if (viterbi_mode) tbm =  MAX(tbm, accum + log(hmm->t[k-1][TDM]));
        else              tbm =  LogSum(tbm, accum + log(hmm->t[k-1][TDM]));
      }

      accum = (hmm->t[k-1][TDD] > 0.0) ? accum + log(hmm->t[k-1][TDD]) : -9999.;
    }
    /* Convert from log_e to scaled integer log_2 odds. */
    if (tbm > -9999.)
      hmm->bsc[k] = (int) floor(0.5 + INTSCALE * 1.44269504 * (tbm - log(hmm->p1)));
    else
      hmm->bsc[k] = -INFTY;
  }

  /* M->E exit transitions. Note how D_M is folded out.
   * M_M is 1 by definition
   * M_M-1 is sum of M_M-1->E and M_M-1->D_M->E, where D_M->E is 1 by definition
   * M_k is sum of M_k->E and M_k->D_k+1...D_M->E
   * Must be done in log space to avoid underflow errors.
   * A little sloppier than it probably has to be; historically,
   * doing in this in log space was in response to a bug report.
   */
  hmm->esc[hmm->M] = 0;
  accum = 0.;
  for (k = hmm->M-1; k >= 1; k--) {
    float tme;
    tme = hmm->end[k] > 0. ? log(hmm->end[k]) : -9999.;
    if (accum > -9999.) {
      if (hmm->t[k][TMD] > 0.0) {
        if (viterbi_mode) tme = MAX(tme, accum + log(hmm->t[k][TMD]));
        else              tme = LogSum(tme, accum + log(hmm->t[k][TMD]));
      }
      accum = (hmm->t[k][TDD] > 0.0) ? accum + log(hmm->t[k][TDD]) : -9999.;
    }
    /* convert from log_e to scaled integer log odds. */
    hmm->esc[k] = (tme > -9999.) ? (int) floor(0.5 + INTSCALE * 1.44269504 * tme) : -INFTY;
  }

  /* special transitions */
  hmm->xsc[XTN][LOOP] = Prob2Score(hmm->xt[XTN][LOOP], hmm->p1);
  hmm->xsc[XTN][MOVE] = Prob2Score(hmm->xt[XTN][MOVE], 1.0);
  hmm->xsc[XTE][LOOP] = Prob2Score(hmm->xt[XTE][LOOP], 1.0);
  hmm->xsc[XTE][MOVE] = Prob2Score(hmm->xt[XTE][MOVE], 1.0);
  hmm->xsc[XTC][LOOP] = Prob2Score(hmm->xt[XTC][LOOP], hmm->p1);
  hmm->xsc[XTC][MOVE] = Prob2Score(hmm->xt[XTC][MOVE], 1.-hmm->p1);
  hmm->xsc[XTJ][LOOP] = Prob2Score(hmm->xt[XTJ][LOOP], hmm->p1);
  hmm->xsc[XTJ][MOVE] = Prob2Score(hmm->xt[XTJ][MOVE], 1.0);

  hmm->flags |= PLAN7_HASBITS;  /* raise the log-odds ready flag */
}


void
Plan7Rescale(
  struct plan7_s *hmm,
  float scale
){
  int k;
  int st;

  /* emissions and transitions in the main model.
   * Note that match states are 1..M, insert states are 1..M-1,
   * and only nodes 1..M-1 have a valid array of transitions.
   */
  for(k = 1; k <= hmm->M; k++)
    FScale(hmm->mat[k], Alphabet_size, scale);
  for(k = 1; k <  hmm->M; k++)
    FScale(hmm->ins[k], Alphabet_size, scale);
  for(k = 1; k <  hmm->M; k++)
    FScale(hmm->t[k],   7,             scale);

  /* begin, end transitions; only valid [1..M] */
  FScale(hmm->begin+1, hmm->M, scale);
  FScale(hmm->end+1,   hmm->M, scale);

  /* B->D1 transition */
  hmm->tbd1 *= scale;

  /* special transitions */
  for (st = 0; st < 4; st++)
    FScale(hmm->xt[st], 2, scale);

  return;
}


void
Plan7Renormalize(
  struct plan7_s *hmm
){
  int   k;      /* counter for model position */
  int   st;      /* counter for special states */
  float d;      /* denominator */

  /* match emissions */
  for (k = 1; k <= hmm->M; k++)
    FNorm(hmm->mat[k], Alphabet_size);
  /* insert emissions */
  for (k = 1; k < hmm->M; k++)
    FNorm(hmm->ins[k], Alphabet_size);
  /* begin transitions */
  d = FSum(hmm->begin+1, hmm->M) + hmm->tbd1;
  FScale(hmm->begin+1, hmm->M, 1./d);
  hmm->tbd1 /= d;
  /* main model transitions */
  for (k = 1; k < hmm->M; k++) {
    d = FSum(hmm->t[k], 3) + hmm->end[k];
    FScale(hmm->t[k], 3, 1./d);
    hmm->end[k] /= d;

    FNorm(hmm->t[k]+3, 2);  /* insert */
    FNorm(hmm->t[k]+5, 2);  /* delete */
  }
  /* null model emissions */
  FNorm(hmm->null, Alphabet_size);
  /* special transitions  */
  for (st = 0; st < 4; st++)
    FNorm(hmm->xt[st], 2);
  /* enforce nonexistent transitions */
  /* (is this necessary?) */
  hmm->t[0][TDM] = hmm->t[0][TDD] = 0.0;

  hmm->flags &= ~PLAN7_HASBITS;  /* clear the log-odds ready flag */
  hmm->flags |= PLAN7_HASPROB;  /* set the probabilities OK flag */
}


void
Plan7RenormalizeExits(
  struct plan7_s *hmm
){
  int   k;

  for (k = 1; k < hmm->M; k++) {
    float d = FSum(hmm->t[k], 3);
    FScale(hmm->t[k], 3, 1./(d + d*hmm->end[k]));
  }
}


void
Plan7NakedConfig(
  struct plan7_s *hmm
){
  hmm->xt[XTN][MOVE] = 1.;        /* disallow N-terminal tail */
  hmm->xt[XTN][LOOP] = 0.;
  hmm->xt[XTE][MOVE] = 1.;        /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1.;        /* disallow C-terminal tail */
  hmm->xt[XTC][LOOP] = 0.;
  hmm->xt[XTJ][MOVE] = 0.;        /* J state unused */
  hmm->xt[XTJ][LOOP] = 1.;
  FSet(hmm->begin+2, hmm->M-1, 0.);   /* disallow internal entries. */
  hmm->begin[1]    = 1. - hmm->tbd1;
  FSet(hmm->end+1,   hmm->M-1, 0.);   /* disallow internal exits. */
  hmm->end[hmm->M] = 1.;
  Plan7RenormalizeExits(hmm);
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}


void
Plan7GlobalConfig(
  struct plan7_s *hmm
){
  hmm->xt[XTN][MOVE] = 1. - hmm->p1;  /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 1.;        /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1. - hmm->p1;  /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 0.;        /* J state unused */
  hmm->xt[XTJ][LOOP] = 1.;
  FSet(hmm->begin+2, hmm->M-1, 0.);   /* disallow internal entries. */
  hmm->begin[1]    = 1. - hmm->tbd1;
  FSet(hmm->end+1,   hmm->M-1, 0.);   /* disallow internal exits. */
  hmm->end[hmm->M] = 1.;
  Plan7RenormalizeExits(hmm);
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}


void
Plan7LSConfig(
  struct plan7_s *hmm
){
  hmm->xt[XTN][MOVE] = 1.-hmm->p1;    /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 0.5;       /* expectation 2 domains/seq */
  hmm->xt[XTE][LOOP] = 0.5;
  hmm->xt[XTC][MOVE] = 1.-hmm->p1;    /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.-hmm->p1;   /* allow J junction state */
  hmm->xt[XTJ][LOOP] = hmm->p1;
  FSet(hmm->begin+2, hmm->M-1, 0.);  /* start at M1/D1 */
  hmm->begin[1]    = 1. - hmm->tbd1;
  FSet(hmm->end+1,   hmm->M-1, 0.);  /* end at M_m/D_m */
  hmm->end[hmm->M] = 1.;
  Plan7RenormalizeExits(hmm);
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}


void
Plan7SWConfig(
  struct plan7_s *hmm,
  float pentry,
  float pexit
){
  float basep;      /* p1 for exits: the base p */
  int   k;      /* counter over states      */

  /* Configure special states.
   */
  hmm->xt[XTN][MOVE] = 1-hmm->p1;    /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 1.;       /* disallow jump state   */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1-hmm->p1;    /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.;           /* J is unused */
  hmm->xt[XTJ][LOOP] = 0.;

  /* Configure entry.
   */
  hmm->begin[1] = (1. - pentry) * (1. - hmm->tbd1);
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.- hmm->tbd1)) / (float)(hmm->M-1));

  /* Configure exit.
   */
  hmm->end[hmm->M] = 1.0;
  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  Plan7RenormalizeExits(hmm);
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}


void
Plan7FSConfig(
  struct plan7_s *hmm,
  float pentry,
  float pexit
){
  float basep;      /* p1 for exits: the base p */
  int   k;      /* counter over states      */

  /* Configure special states.
   */
  hmm->xt[XTN][MOVE] = 1-hmm->p1;    /* allow N-terminal tail     */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 0.5;       /* allow loops / multihits   */
  hmm->xt[XTE][LOOP] = 0.5;
  hmm->xt[XTC][MOVE] = 1-hmm->p1;    /* allow C-terminal tail     */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.-hmm->p1;   /* allow J junction between domains */
  hmm->xt[XTJ][LOOP] = hmm->p1;

  /* Configure entry.
   */
  hmm->begin[1] = (1. - pentry) * (1. - hmm->tbd1);
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.-hmm->tbd1)) / (float)(hmm->M-1));

  /* Configure exit.
   */
  hmm->end[hmm->M] = 1.0;
  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  Plan7RenormalizeExits(hmm);
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}


int
DegenerateSymbolScore(
  float *p,
  float *null,
  int ambig
){
  int x;
  float numer = 0.;
  float denom = 0.;

  for (x = 0; x < Alphabet_size; x++) {
    if (Degenerate[ambig][x]) {
      numer += null[x] * sreLOG2(p[x] / null[x]);
      denom += null[x];
    }
  }
  return (int) (INTSCALE * numer / denom);
}


void
Plan9toPlan7(
  struct plan9_s *hmm,
  struct plan7_s **ret_plan7
){
  struct plan7_s *plan7;
  int k, x;

  plan7 = AllocPlan7(hmm->M);

  for (k = 1; k < hmm->M; k++) {
    plan7->t[k][TMM] = hmm->mat[k].t[MATCH];
    plan7->t[k][TMD] = hmm->mat[k].t[DELETE];
    plan7->t[k][TMI] = hmm->mat[k].t[INSERT];
    plan7->t[k][TDM] = hmm->del[k].t[MATCH];
    plan7->t[k][TDD] = hmm->del[k].t[DELETE];
    plan7->t[k][TIM] = hmm->ins[k].t[MATCH];
    plan7->t[k][TII] = hmm->ins[k].t[INSERT];
  }

  for (k = 1; k <= hmm->M; k++)
    for (x = 0; x < Alphabet_size; x++)
      plan7->mat[k][x] = hmm->mat[k].p[x];

  for (k = 1; k < hmm->M; k++)
    for (x = 0; x < Alphabet_size; x++)
      plan7->ins[k][x] = hmm->ins[k].p[x];

  plan7->tbd1 = hmm->mat[0].t[DELETE] / (hmm->mat[0].t[DELETE] + hmm->mat[0].t[MATCH]);

  /* We have to make up the null transition p1; use default */
  P7DefaultNullModel(plan7->null, &(plan7->p1));
  for (x = 0; x < Alphabet_size; x++)
    plan7->null[x] = hmm->null[x];

  if (hmm->name != NULL)
    Plan7SetName(plan7, hmm->name);
  if (hmm->flags & HMM_REF) {
    strcpy(plan7->rf, hmm->ref);
    plan7->flags |= PLAN7_RF;
  }
  if (hmm->flags & HMM_CS) {
    strcpy(plan7->cs, hmm->cs);
    plan7->flags |= PLAN7_CS;
  }

  Plan7LSConfig(plan7);    /* configure specials for ls-style alignment */
  Plan7Renormalize(plan7);  /* mainly to correct for missing ID and DI */
  plan7->flags |= PLAN7_HASPROB;  /* probabilities are valid */
  plan7->flags &= ~PLAN7_HASBITS;  /* scores are not valid    */
  *ret_plan7 = plan7;
}

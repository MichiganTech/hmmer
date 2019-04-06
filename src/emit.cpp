/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* emit.c
 *
 * Generation of sequences/traces from an HMM.
 */

#include "emit.hpp"
#include "vectorops.hpp"


void
EmitSequence(
  struct plan7_s *hmm,
  unsigned char **ret_dsq,
  int *ret_L,
  struct p7trace_s **ret_tr
){
  struct p7trace_s *tr;
  unsigned char    *dsq;        /* generated sequence, digitized */
  char  type;                 /* current state type */
  int   k;      /* current node index */
  int   L;      /* length of sequence */
  int   alloc_tlen;    /* allocated space for traceback */
  int   alloc_L;    /* allocated space for sequence  */
  int   tpos;      /* position in traceback */
  float t[4];      /* little array for choosing M transition from */

  /* Initialize; allocations
   */
  P7AllocTrace(64, &tr);
  alloc_tlen = 64;
  dsq = MallocOrDie(sizeof(unsigned char) * 64);
  alloc_L = 64;

  TraceSet(tr, 0, STS, 0, 0);
  TraceSet(tr, 1, STN, 0, 0);
  dsq[0] = Alphabet_iupac;
  L      = 1;
  k      = 0;
  type   = STN;
  tpos   = 2;

  while (type != STT) {
    /* Deal with state transition
     */
    switch (type) {
    case STB:
      hmm->begin[0] = hmm->tbd1; /* begin[0] hack (documented in structs.h) */
      k = FChoose(hmm->begin, hmm->M+1);
      if (k == 0) {
        type = STD;
        k = 1;
      } else {
        type = STM;
      }
      break;

    case STI:
      type = (FChoose(hmm->t[k]+TIM, 2) == 0)    ? STM : STI;
      if (type == STM) k++;
      break;
    case STN:
      type = (FChoose(hmm->xt[XTN], 2)  == LOOP) ? STN : STB;
      k = 0;
      break;
    case STE:
      type = (FChoose(hmm->xt[XTE], 2)  == LOOP) ? STJ : STC;
      k = 0;
      break;
    case STC:
      type = (FChoose(hmm->xt[XTC], 2)  == LOOP) ? STC : STT;
      k = 0;
      break;
    case STJ:
      type = (FChoose(hmm->xt[XTJ], 2)  == LOOP) ? STJ : STB;
      k = 0;
      break;

    case STD:
      if (k < hmm->M) {
        type = (FChoose(hmm->t[k]+TDM, 2) == 0) ? STM : STD;
        k++;
      } else {
        type = STE;
        k = 0;
      }
      break;

    case STM:
      if (k < hmm->M) {
        FCopy(t, hmm->t[k], 3);
        t[3] = hmm->end[k];
        switch (FChoose(t,4)) {
        case 0:
          k++;
          type = STM;
          break;
        case 1:
          type = STI;
          break;
        case 2:
          k++;
          type = STD;
          break;
        case 3:
          k=0;
          type = STE;
          break;
        default:
          Die("never happens");
        }
      } else {
        k    = 0;
        type = STE;
      }
      break;

    case STT:
    case STBOGUS:
    default:
      Die("can't happen.");
    }

    /* Choose a symbol emission, if necessary
     */

    unsigned char sym;    /* a generated symbol index */
    sym = Alphabet_iupac;  /* use sentinel byte to mean "not set yet" */
    if      (type == STM) sym = FChoose(hmm->mat[k], Alphabet_size);
    else if (type == STI) sym = FChoose(hmm->ins[k], Alphabet_size);
    else if ((type == STN && tr->statetype[tpos-1] == STN) ||
             (type == STC && tr->statetype[tpos-1] == STC) ||
             (type == STJ && tr->statetype[tpos-1] == STJ))
      sym = FChoose(hmm->null, Alphabet_size);

    /* Add to the traceback; deal with realloc if necessary
     */
    TraceSet(tr, tpos, type, k, (sym != Alphabet_iupac) ? L : 0);
    tpos++;
    if (tpos == alloc_tlen) {
      alloc_tlen += 64;
      P7ReallocTrace(tr, alloc_tlen);
    }

    /* Add to the digitized seq; deal with realloc, if necessary
     */
    if (sym != Alphabet_iupac) {
      dsq[L] = sym;
      L++;
      if (L == alloc_L) {
        alloc_L += 64;
        dsq = ReallocOrDie(dsq, sizeof(unsigned char) * alloc_L);
      }
    }
  }

  /* Finish off the trace
   */
  tr->tlen = tpos;

  /* Finish off the dsq with sentinel byte and null terminator.
   * Emitted sequence length is L-1. No, this isn't a bug - L is
   * the index of the next position, and since we're ending, the
   * sequence is L-1 and the sentinel goes in position L.
   */
  dsq[L]   = Alphabet_iupac;
  L--;

  /* Return
   */
  if (ret_dsq != NULL) *ret_dsq = dsq;
  else free(dsq);
  if (ret_L   != NULL) *ret_L   = L;
  if (ret_tr  != NULL) *ret_tr  = tr;
  else P7FreeTrace(tr);
  return;
}


void
EmitConsensusSequence(
  struct plan7_s *hmm,
  char **ret_seq,
  unsigned char **ret_dsq,
  int *ret_L,
  struct p7trace_s **ret_tr
){
  struct p7trace_s *tr;         /* RETURN: traceback */
  unsigned char    *dsq;        /* RETURN: sequence in digitized form */
  char             *seq;        /* RETURN: sequence in undigitized form */
  float *mp, *ip, *dp;          /* state occupancies from StateOccupancy() */
  int    nmat, ndel, nins;  /* number of matches, deletes, inserts used */
  int    k;      /* counter for nodes */
  int    tpos;      /* position in trace */
  int    i;                     /* position in seq (equiv pos in dsq is i+1 */
  int    x;      /* symbol choice (M) or # symbols (I) */
  float  mthresh;    /* >= this, show symbol as upper case */

  if (Alphabet_type == hmmAMINO) mthresh = 0.5;
  else                           mthresh = 0.9;

  StateOccupancy(hmm, &mp, &ip, &dp);

  /* First pass: how many states do we need in the trace?
   *             how long will the sequence be?
   */
  nmat = ndel = nins = 0;
  for (k = 1; k <= hmm->M; k++) {
    if (mp[k] >= 0.5) nmat++;
    else ndel++;
    if (k < hmm->M && ip[k] >= 0.5)
      nins += (int) (1.f / (1.f - hmm->t[k][TII]));
  }

  /* Allocations
   */
  P7AllocTrace(6 + nmat + ndel + nins, &tr);
  dsq = MallocOrDie(sizeof(unsigned char) * (nmat+nins+3));
  seq = MallocOrDie(sizeof(char) * (nmat+nins+1));

  /* Main pass.
   * Construct consensus trace, seq, and dsq.
   */
  TraceSet(tr, 0, STS, 0, 0);
  TraceSet(tr, 1, STN, 0, 0);
  TraceSet(tr, 2, STB, 0, 0);
  dsq[0] = Alphabet_iupac;  /* guard byte */
  tpos = 3;
  i    = 0;
  for (k = 1; k <= hmm->M; k++) {
    if (mp[k] >= 0.5) {
      x = FArgMax(hmm->mat[k], Alphabet_size);
      TraceSet(tr, tpos, STM, k, i+1);
      seq[i]   = Alphabet[x];
      dsq[i+1] = x;
      if (hmm->mat[k][x] < mthresh)
        seq[i] = tolower((int) seq[i]);
      i++;
      tpos++;
    } else {
      TraceSet(tr, tpos, STD, k, 0);
      tpos++;
    }

    if (k < hmm->M && ip[k] >= 0.5) {
      x = (int) (1.f / (1.f - hmm->t[k][TII]));
      while (x--) {
        TraceSet(tr, tpos, STI, k, i+1);
        seq[i]   = 'x';
        dsq[i+1] = Alphabet_iupac - 1;
        i++;
        tpos++;
      }
    }
  }
  TraceSet(tr, tpos, STE, 0, 0);
  tpos++;
  TraceSet(tr, tpos, STC, 0, 0);
  tpos++;
  TraceSet(tr, tpos, STT, 0, 0);
  tpos++;
  dsq[i+1] = Alphabet_iupac;

  free(mp);
  free(ip);
  free(dp);
  if (ret_seq != NULL) *ret_seq = seq;
  else free(seq);
  if (ret_dsq != NULL) *ret_dsq = dsq;
  else free(dsq);
  if (ret_L   != NULL) *ret_L   = i;
  if (ret_tr  != NULL) *ret_tr  = tr;
  else P7FreeTrace(tr);
}


void
StateOccupancy(
  struct plan7_s *hmm,
  float **ret_mp,
  float **ret_ip,
  float **ret_dp
){
  float *fmp, *fip, *fdp;       /* forward probabilities  */
  int k;      /* counter for nodes      */

  /* Initial allocations
   */
  fmp = MallocOrDie (sizeof(float) * (hmm->M+1));
  fip = MallocOrDie (sizeof(float) * (hmm->M));
  fdp = MallocOrDie (sizeof(float) * (hmm->M+1));

  /* Forward pass.
   */
  fdp[1] = hmm->tbd1;
  fmp[1] = hmm->begin[1];
  fip[1] = fmp[1] * hmm->t[1][TMI];
  for (k = 2; k <= hmm->M; k++) {
    /* M: from M,D,I at k-1, or B; count t_II as 1.0 */
    fmp[k] = fmp[k-1] * hmm->t[k-1][TMM] +
             fip[k-1] +
             fdp[k-1] * hmm->t[k-1][TDM] +
             hmm->begin[k];
    /* D: from M,D at k-1 */
    fdp[k] = fmp[k-1] * hmm->t[k-1][TMD] +
             fdp[k-1] * hmm->t[k-1][TDD];
    /* I: from M at k; don't count II */
    if (k < hmm->M) {
      fip[k] = fmp[k] * hmm->t[k][TMI];
    }

    assert((fabs(1.0f - fmp[k] - fdp[k]) < 1e-6f));
    fmp[k] /= fmp[k]+fdp[k];  /* prevent propagating fp errors */
    fdp[k] /= fmp[k]+fdp[k];
  }
  /* We don't need a backward pass; all backwards P's are 1.0
   * by definition (you can always get out of a state with P=1).
   * The only situation where this might not be true is for
   * a TII of 1.0, when TIM = 0 -- but in that case, if there's
   * a finite chance of getting into that insert state, the model
   * generates infinitely long sequences, so we can consider this
   * situation "perverse" and disallow it elsewhere in building
   * profile HMMs.
   */

  /* Return.
   */
  *ret_mp = fmp;
  *ret_dp = fdp;
  *ret_ip = fip;
}
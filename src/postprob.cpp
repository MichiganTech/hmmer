/************************************************************
 * Copyright (C) 1998 Ian Holmes (ihh@sanger.ac.uk)
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

#include "postprob.hpp"


float
P7Backward(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s **ret_mx
){
  struct dpmatrix_s *mx;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */
  mx = AllocPlan7Matrix(L+1, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the L row.
   * Note that xmx[i][stS] = xmx[i][stN] by definition for all i,
   *    so stS need not be calculated in backward DP matrices.
   */
  xmx[L][XMC] = hmm->xsc[XTC][MOVE];                 /* C<-T                          */
  xmx[L][XME] = xmx[L][XMC] + hmm->xsc[XTE][MOVE];   /* E<-C, no C-tail               */
  xmx[L][XMJ] = xmx[L][XMB] = xmx[L][XMN] = -INFTY;  /* need seq to get out from here */
  for (k = hmm->M; k >= 1; k--) {
    mmx[L][k] = xmx[L][XME] + hmm->esc[k];           /* M<-E ...                      */
    mmx[L][k] += hmm->msc[dsq[L]][k];                /* ... + emitted match symbol    */
    imx[L][k] = dmx[L][k] = -INFTY;                  /* need seq to get out from here */
  }

  /* Recursion. Done as a pull.
   * Note slightly wasteful boundary conditions:
   *    M_M precalculated, D_M set to -INFTY,
   *    D_1 wastefully calculated.
   * Scores for transitions to D_M also have to be hacked to -INFTY,
   * as Plan7Logoddsify does not do this for us (I think? - ihh).
   */
  hmm->tsc[TDD][hmm->M-1] = hmm->tsc[TMD][hmm->M-1] = -INFTY;    /* no D_M state -- HACK -- should be in Plan7Logoddsify */
  for (i = L-1; i >= 0; i--) {
    /* Do the special states first.
     * remember, C, N and J emissions are zero score by definition
     */
    xmx[i][XMC] = xmx[i+1][XMC] + hmm->xsc[XTC][LOOP];

    xmx[i][XMB] = -INFTY;
    /* The following section has been hacked to fit a bug in core_algorithms.c
     * The "correct" code is:
     * for (k = hmm->M; k >= 1; k--)
     * xmx[i][XMB] = ILogsum(xmx[i][XMB], mmx[i+1][k] + hmm->bsc[k];
     *
     * The following code gives the same results as core_algorithms.c:
     */
    xmx[i][XMB] = ILogsum(xmx[i][XMB], mmx[i+1][hmm->M] + hmm->bsc[hmm->M-1]);
    for (k = hmm->M-1; k >= 1; k--)
      xmx[i][XMB] = ILogsum(xmx[i][XMB], mmx[i+1][k] + hmm->bsc[k]);

    xmx[i][XMJ] = ILogsum(xmx[i][XMB] + hmm->xsc[XTJ][MOVE],
                          xmx[i+1][XMJ] + hmm->xsc[XTJ][LOOP]);

    xmx[i][XME] = ILogsum(xmx[i][XMC] + hmm->xsc[XTE][MOVE],
                          xmx[i][XMJ] + hmm->xsc[XTE][LOOP]);

    xmx[i][XMN] = ILogsum(xmx[i][XMB] + hmm->xsc[XTN][MOVE],
                          xmx[i+1][XMN] + hmm->xsc[XTN][LOOP]);

    /* Now the main states. Note the boundary conditions at M.
     */

    if (i>0) {
      mmx[i][hmm->M] = xmx[i][XME] + hmm->esc[hmm->M] + hmm->msc[dsq[i]][hmm->M];
      dmx[i][hmm->M] = -INFTY;
      for (k = hmm->M-1; k >= 1; k--) {
        mmx[i][k]  = ILogsum(ILogsum(xmx[i][XME] + hmm->esc[k],
                                     mmx[i+1][k+1] + hmm->tsc[TMM][k]),
                             ILogsum(imx[i+1][k] + hmm->tsc[TMI][k],
                                     dmx[i][k+1] + hmm->tsc[TMD][k]));
        mmx[i][k] += hmm->msc[dsq[i]][k];

        imx[i][k]  = ILogsum(imx[i+1][k] + hmm->tsc[TII][k],
                             mmx[i+1][k+1] + hmm->tsc[TIM][k]);
        imx[i][k] += hmm->isc[dsq[i]][k];

        dmx[i][k]  = ILogsum(dmx[i][k+1] + hmm->tsc[TDD][k],
                             mmx[i+1][k+1] + hmm->tsc[TDM][k]);

      }
    }

  }

  sc = xmx[0][XMN];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreePlan7Matrix(mx);

  return Scorify(sc);    /* the total Backward score. */
}



void
P7EmitterPosterior(
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s *forward,
  struct dpmatrix_s *backward,
  struct dpmatrix_s *mx
){
  int i;
  int k;
  int sc;

  sc = backward->xmx[0][XMN];

  for (i = L; i >= 1; i--) {
    mx->xmx[i][XMC] = forward->xmx[i-1][XMC] + hmm->xsc[XTC][LOOP] + backward->xmx[i][XMC] - sc;

    mx->xmx[i][XMJ] = forward->xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP] + backward->xmx[i][XMJ] - sc;

    mx->xmx[i][XMN] = forward->xmx[i-1][XMN] + hmm->xsc[XTN][LOOP] + backward->xmx[i][XMN] - sc;

    mx->xmx[i][XMB] = mx->xmx[i][XME] = -INFTY;

    for (k = 1; k < hmm->M; k++) {
      mx->mmx[i][k]  = backward->mmx[i][k];
      mx->mmx[i][k] += ILogsum(ILogsum(forward->mmx[i-1][k-1] + hmm->tsc[TMM][k-1],
                                       forward->imx[i-1][k-1] + hmm->tsc[TIM][k-1]),
                               ILogsum(forward->xmx[i-1][XMB] + hmm->bsc[k],
                                       forward->dmx[i-1][k-1] + hmm->tsc[TDM][k-1]));
      mx->mmx[i][k] -= sc;

      mx->imx[i][k]  = backward->imx[i][k];
      mx->imx[i][k] += ILogsum(forward->mmx[i-1][k] + hmm->tsc[TMI][k],
                               forward->imx[i-1][k] + hmm->tsc[TII][k]);
      mx->imx[i][k] -= sc;

      mx->dmx[i][k] = -INFTY;
    }
    mx->mmx[i][hmm->M]  = backward->mmx[i][hmm->M];
    mx->mmx[i][hmm->M] += ILogsum(ILogsum(forward->mmx[i-1][hmm->M-1] + hmm->tsc[TMM][hmm->M-1],
                                          forward->imx[i-1][hmm->M-1] + hmm->tsc[TIM][hmm->M-1]),
                                  ILogsum(forward->xmx[i-1][XMB] + hmm->bsc[hmm->M],
                                          forward->dmx[i-1][hmm->M-1] + hmm->tsc[TDM][hmm->M-1]));
    mx->mmx[i][hmm->M] -= sc;

    mx->imx[i][hmm->M] = mx->dmx[i][hmm->M] = mx->dmx[i][0] = -INFTY;

  }
}


float P7FillOptimalAccuracy(
  int L,
  int M,
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr
){
  struct p7trace_s  *tr;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization of the zero row.
   * Each cell in the optimal accuracy matrix holds the log of the expected
   * of correctly assigned symbols up to that point.
   * To begin with, everything is log(0) = -INFTY.
   */
  xmx[0][XMN] = xmx[0][XMB] = xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;
  for (k = 0; k <= M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

    for (k = 1; k <= M; k++) {
      /* match state */
      mmx[i][k]  = -INFTY;
      if ((sc = mmx[i-1][k-1]) > mmx[i][k])
        mmx[i][k] = sc;
      if ((sc = imx[i-1][k-1]) > mmx[i][k])
        mmx[i][k] = sc;
      if ((sc = dmx[i-1][k-1]) > mmx[i][k])
        mmx[i][k] = sc;
      if ((sc = xmx[i-1][XMB]) > mmx[i][k])
        mmx[i][k] = sc;
      mmx[i][k] = ILogsum(mmx[i][k], posterior->mmx[i][k]);

      /* delete state */
      dmx[i][k] = -INFTY;
      if ((sc = mmx[i][k-1]) > dmx[i][k])
        dmx[i][k] = sc;
      if ((sc = dmx[i][k-1]) > dmx[i][k])
        dmx[i][k] = sc;

      /* insert state */
      imx[i][k] = -INFTY;
      if ((sc = mmx[i-1][k]) > imx[i][k])
        imx[i][k] = sc;
      if ((sc = imx[i-1][k]) > imx[i][k])
        imx[i][k] = sc;
      imx[i][k] = ILogsum(imx[i][k], posterior->imx[i][k]);
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */

    /* N state */
    xmx[i][XMN] = -INFTY;
    if ((sc = ILogsum(xmx[i-1][XMN], posterior->xmx[i][XMN])) > -INFTY)
      xmx[i][XMN] = sc;

    /* E state */
    xmx[i][XME] = -INFTY;
    for (k = 1; k <= M; k++)
      if ((sc =  mmx[i][k]) > xmx[i][XME])
        xmx[i][XME] = sc;

    /* J state */
    xmx[i][XMJ] = -INFTY;
    if ((sc = ILogsum(xmx[i-1][XMJ], posterior->xmx[i][XMJ])) > -INFTY)
      xmx[i][XMJ] = sc;
    if ((sc = xmx[i][XME]) > xmx[i][XMJ])      /* no E->J emission */
      xmx[i][XMJ] = sc;

    /* B state */
    xmx[i][XMB] = -INFTY;
    if ((sc = xmx[i][XMN]) > -INFTY)
      xmx[i][XMB] = sc;
    if ((sc = xmx[i][XMJ]) > xmx[i][XMB])
      xmx[i][XMB] = sc;

    /* C state */
    xmx[i][XMC] = -INFTY;
    if ((sc = ILogsum(xmx[i-1][XMC], posterior->xmx[i][XMC])) > -INFTY)
      xmx[i][XMC] = sc;
    if ((sc = xmx[i][XME]) > xmx[i][XMC])      /* no E->C emission */
      xmx[i][XMC] = sc;
  }

  /* T state (not stored) */
  sc = xmx[L][XMC];

  if (ret_tr != NULL) {
    P7OptimalAccuracyTrace(L, M, posterior, mx, &tr);
    *ret_tr = tr;
  }

  return Score2Prob(sc,1);  /* the log of the expected accuracy. */
}


void
P7OptimalAccuracyTrace(
  int L,
  int M,
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr
){
  struct p7trace_s *tr;
  int curralloc;    /* current allocated length of trace */
  int tpos;      /* position in trace */
  int i;      /* position in seq (1..L) */
  int k = 0;      /* position in model (1..M) */
  int **xmx, **mmx, **imx, **dmx;
  int sc;      /* temp var for pre-emission score */

  /* Overallocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + L is minimum trace;
   * add L more as buffer.
   */
  curralloc = L * 2 + 6;
  P7AllocTrace(curralloc, &tr);

  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization of trace
   * We do it back to front; ReverseTrace() is called later.
   */
  tr->statetype[0] = STT;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STC;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
  i    = L;      /* current i (seq pos) we're trying to assign */

  /* Traceback
   */
  while (tr->statetype[tpos-1] != STS) {
    switch (tr->statetype[tpos-1]) {
    case STM:      /* M connects from i-1,k-1, or B */
      sc = mmx[i+1][k+1];
      if (sc == ILogsum(mmx[i][k], posterior->mmx[i+1][k+1]) && i > 0 && k > 0) {
        tr->statetype[tpos] = STM;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = i--;
      } else if (sc == ILogsum(imx[i][k], posterior->mmx[i+1][k+1]) && i > 0 && k > 0) {
        tr->statetype[tpos] = STI;
        tr->nodeidx[tpos]   = k;
        tr->pos[tpos]       = i--;
      } else if (sc == ILogsum(dmx[i][k], posterior->mmx[i+1][k+1]) && i > 0 && k > 1) {
        tr->statetype[tpos] = STD;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = 0;
      } else if (sc == ILogsum(xmx[i][XMB], posterior->mmx[i+1][k+1])) {
        tr->statetype[tpos] = STB;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      } else Die("traceback failed");
      break;

    case STD:      /* D connects from M,D */
      if (dmx[i][k+1] == mmx[i][k] && i > 0 && k > 0) {
        tr->statetype[tpos] = STM;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = i--;
      } else if (dmx[i][k+1] == dmx[i][k] && k > 1) {
        tr->statetype[tpos] = STD;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = 0;
      } else Die("traceback failed");
      break;

    case STI:      /* I connects from M,I */
      sc = imx[i+1][k];
      if (sc == ILogsum(mmx[i][k], posterior->imx[i+1][k]) && i > 0 && k > 0) {
        tr->statetype[tpos] = STM;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = i--;
      } else if (sc == ILogsum(imx[i][k], posterior->imx[i+1][k]) && i > 0 && k > 0) {
        tr->statetype[tpos] = STI;
        tr->nodeidx[tpos]   = k;
        tr->pos[tpos]       = i--;
      } else Die("traceback failed");
      break;

    case STN:      /* N connects from S, N */
      if (i == 0 && xmx[i][XMN] == -INFTY) {
        tr->statetype[tpos] = STS;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      } else if (i > 0 && xmx[i+1][XMN] == ILogsum(xmx[i][XMN], posterior->xmx[i+1][XMN]) && i > 0) {
        tr->statetype[tpos] = STN;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;    /* note convention adherence:  */
        tr->pos[tpos-1]     = i--;  /* first N doesn't emit        */
      } else Die("traceback failed");
      break;

    case STB:      /* B connects from N, J */
      if (xmx[i][XMB] == xmx[i][XMN]) {
        tr->statetype[tpos] = STN;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      } else if (xmx[i][XMB] == xmx[i][XMJ]) {
        tr->statetype[tpos] = STJ;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      } else Die("traceback failed");
      break;

    case STE:      /* E connects from any M state. k set here */
      for (k = M; k >= 1; k--)
        if (xmx[i][XME] == mmx[i][k] && i > 0) {
          tr->statetype[tpos] = STM;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = i--;
          break;
        }
      if (k <= 0) Die("traceback failed");
      break;

    case STC:      /* C comes from C, E */
      if (xmx[i][XMC] == ILogsum(xmx[i-1][XMC], posterior->xmx[i][XMC]) && i > 0) {
        tr->statetype[tpos] = STC;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;    /* note convention adherence: */
        tr->pos[tpos-1]     = i--;  /* first C doesn't emit       */
      } else if (xmx[i][XMC] == xmx[i][XME]) {
        tr->statetype[tpos] = STE;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0; /* E is a nonemitter */
      } else Die("Traceback failed.");
      break;

    case STJ:      /* J connects from E, J */
      if (xmx[i][XMJ] == ILogsum(xmx[i-1][XMJ], posterior->xmx[i][XMJ]) && i > 0) {
        tr->statetype[tpos] = STJ;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;    /* note convention adherence: */
        tr->pos[tpos-1]     = i--;  /* first J doesn't emit       */
      } else if (xmx[i][XMJ] == xmx[i][XME]) {
        tr->statetype[tpos] = STE;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0; /* E is a nonemitter */
      } else Die("Traceback failed.");
      break;

    default:
      Die("traceback failed");

    } /* end switch over statetype[tpos-1] */

    tpos++;
    if (tpos == curralloc) {
      /* grow trace if necessary  */
      curralloc += L;
      P7ReallocTrace(tr, curralloc);
    }

  } /* end traceback, at S state; tpos == tlen now */
  tr->tlen = tpos;
  P7ReverseTrace(tr);
  *ret_tr = tr;
}
/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* core_algorithms.c
 *
 * Simple and robust "research" implementations of Forward, Backward,
 * and Viterbi for Plan7. For optimized replacements for some of these functions,
 * see fast_algorithms.c.
 */

#include <string.h>
#include <assert.h>

#include "algorithms.hpp"
#include "config.hpp"
#include "structs.hpp"
#include "vectorops.hpp"
#include "structs.hpp"


struct dpmatrix_s*
AllocPlan7Matrix(
  int rows,
  int M,
  int ***xmx,
  int ***mmx,
  int ***imx,
  int ***dmx
){
  struct dpmatrix_s *mx;
  mx = CreatePlan7Matrix(rows-1, M, 0, 0);
  if (xmx != NULL) *xmx = mx->xmx;
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  return mx;
}


void
FreePlan7Matrix(
  struct dpmatrix_s *mx
){
  free (mx->xmx_mem);
  free (mx->mmx_mem);
  free (mx->imx_mem);
  free (mx->dmx_mem);
  free (mx->xmx);
  free (mx->mmx);
  free (mx->imx);
  free (mx->dmx);
  free (mx);
}


struct dpshadow_s*
AllocShadowMatrix(
  int rows,
  int M,
  char ***xtb,
  char ***mtb,
  char ***itb,
  char ***dtb
){
  struct dpshadow_s *tb;
  int i;

  tb         = (struct dpshadow_s *) MallocOrDie (sizeof(struct dpshadow_s));
  tb->xtb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->mtb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->itb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->dtb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->esrc   = (int *)   MallocOrDie (sizeof(int)  * rows);
  tb->xtb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*5));
  tb->mtb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*(M+2)));
  tb->itb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*(M+2)));
  tb->dtb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*(M+2)));
  for (i = 1; i < rows; i++) {
    tb->xtb[i] = tb->xtb[0] + (i*5);
    tb->mtb[i] = tb->mtb[0] + (i*(M+2));
    tb->itb[i] = tb->itb[0] + (i*(M+2));
    tb->dtb[i] = tb->dtb[0] + (i*(M+2));
  }

  if (xtb != NULL) *xtb = tb->xtb;
  if (mtb != NULL) *mtb = tb->mtb;
  if (itb != NULL) *itb = tb->itb;
  if (dtb != NULL) *dtb = tb->dtb;
  return tb;
}


void
FreeShadowMatrix(
  struct dpshadow_s *tb
){
  free (tb->xtb[0]);
  free (tb->mtb[0]);
  free (tb->itb[0]);
  free (tb->dtb[0]);
  free (tb->esrc);
  free (tb->xtb);
  free (tb->mtb);
  free (tb->itb);
  free (tb->dtb);
  free (tb);
}


int
P7ViterbiSpaceOK(
  int L,
  int M,
  struct dpmatrix_s *mx
){
  int newM;
  int newN;

  if (M <= mx->maxM && L <= mx->maxN) return true;

  if (M > mx->maxM) newM = M + mx->padM;
  else newM = mx->maxM;
  if (L > mx->maxN) newN = L + mx->padN;
  else newN = mx->maxN;

  if (P7ViterbiSize(newN, newM) <= RAMLIMIT)
    return true;
  else
    return false;
}


int
P7ViterbiSize(
  int L,
  int M
){
  float Mbytes;

  /* We're excessively precise here, but it doesn't cost
   * us anything to be pedantic. The four terms are:
   *   1. the matrix structure itself;
   *   2. the O(NM) main matrix (this dominates!)
   *   3. ptrs into the rows of the matrix
   *   4. storage for 5 special states. (xmx)
   */
  Mbytes =  (float) sizeof(struct dpmatrix_s);
  Mbytes += 3. * (float) (L+1) * (float) (M+2) * (float) sizeof(int);
  Mbytes += 4. * (float) (L+1) * (float) sizeof(int *);
  Mbytes += 5. * (float) (L+1) * (float) sizeof(int);
  Mbytes /= 1048576.;
  return (int) Mbytes;
}


float
P7Forward(
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

  /* Initialization of the zero row.
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  xmx[0][XMN] = 0;                         /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   */
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;
    for (k = 1; k < hmm->M; k++) {
      mmx[i][k]  = ILogsum(ILogsum(mmx[i-1][k-1] + hmm->tsc[TMM][k-1],
                                   imx[i-1][k-1] + hmm->tsc[TIM][k-1]),
                           ILogsum(xmx[i-1][XMB] + hmm->bsc[k],
                                   dmx[i-1][k-1] + hmm->tsc[TDM][k-1]));
      mmx[i][k] += hmm->msc[dsq[i]][k];

      dmx[i][k]  = ILogsum(mmx[i][k-1] + hmm->tsc[TMD][k-1],
                           dmx[i][k-1] + hmm->tsc[TDD][k-1]);
      imx[i][k]  = ILogsum(mmx[i-1][k] + hmm->tsc[TMI][k],
                           imx[i-1][k] + hmm->tsc[TII][k]);
      imx[i][k] += hmm->isc[dsq[i]][k];
    }
    mmx[i][hmm->M] = ILogsum(ILogsum(mmx[i-1][hmm->M-1] + hmm->tsc[TMM][hmm->M-1],
                                     imx[i-1][hmm->M-1] + hmm->tsc[TIM][hmm->M-1]),
                             ILogsum(xmx[i-1][XMB] + hmm->bsc[hmm->M],
                                     dmx[i-1][hmm->M-1] + hmm->tsc[TDM][hmm->M-1]));
    mmx[i][hmm->M] += hmm->msc[dsq[i]][hmm->M];

    /* Now the special states.
     * remember, C and J emissions are zero score by definition
     */
    xmx[i][XMN] = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP];

    xmx[i][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      xmx[i][XME] = ILogsum(xmx[i][XME], mmx[i][k] + hmm->esc[k]);

    xmx[i][XMJ] = ILogsum(xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP],
                          xmx[i][XME]   + hmm->xsc[XTE][LOOP]);

    xmx[i][XMB] = ILogsum(xmx[i][XMN] + hmm->xsc[XTN][MOVE],
                          xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]);

    xmx[i][XMC] = ILogsum(xmx[i-1][XMC] + hmm->xsc[XTC][LOOP],
                          xmx[i][XME] + hmm->xsc[XTE][MOVE]);
  }

  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreePlan7Matrix(mx);

  return Scorify(sc);    /* the total Forward score. */
}


// float
// P7Viterbi(
//   unsigned char *dsq,
//   int L,
//   struct plan7_s *hmm,
//   struct dpmatrix_s *mx,
//   struct p7trace_s **ret_tr
// ){
//   struct p7trace_s  *tr;
//   int **xmx;
//   int **mmx;
//   int **imx;
//   int **dmx;
//   int   i,k;
//   int   sc;

//   /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
//    */
//   ResizePlan7Matrix(mx, L, hmm->M, &xmx, &mmx, &imx, &dmx);

//   /* Initialization of the zero row.
//    */
//   xmx[0][XMN] = 0;                         /* S->N, p=1            */
//   xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
//   xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
//   for (k = 0; k <= hmm->M; k++)
//     mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

//   /* Recursion. Done as a pull.
//    * Note some slightly wasteful boundary conditions:
//    *    tsc[0] = -INFTY for all eight transitions (no node 0)
//    *    D_M and I_M are wastefully calculated (they don't exist)
//    */
//   for (i = 1; i <= L; i++) {
//     mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

//     for (k = 1; k <= hmm->M; k++) {
//       /* match state */
//       mmx[i][k]  = -INFTY;
//       if ((sc = mmx[i-1][k-1] + hmm->tsc[TMM][k-1]) > mmx[i][k])
//         mmx[i][k] = sc;
//       if ((sc = imx[i-1][k-1] + hmm->tsc[TIM][k-1]) > mmx[i][k])
//         mmx[i][k] = sc;
//       if ((sc = xmx[i-1][XMB] + hmm->bsc[k]) > mmx[i][k])
//         mmx[i][k] = sc;
//       if ((sc = dmx[i-1][k-1] + hmm->tsc[TDM][k-1]) > mmx[i][k])
//         mmx[i][k] = sc;
//       if (hmm->msc[dsq[i]][k] != -INFTY) mmx[i][k] += hmm->msc[dsq[i]][k];
//       else                                     mmx[i][k] = -INFTY;

//       /* delete state */
//       dmx[i][k] = -INFTY;
//       if ((sc = mmx[i][k-1] + hmm->tsc[TMD][k-1]) > dmx[i][k])
//         dmx[i][k] = sc;
//       if ((sc = dmx[i][k-1] + hmm->tsc[TDD][k-1]) > dmx[i][k])
//         dmx[i][k] = sc;

//       /* insert state */
//       if (k < hmm->M) {
//         imx[i][k] = -INFTY;
//         if ((sc = mmx[i-1][k] + hmm->tsc[TMI][k]) > imx[i][k])
//           imx[i][k] = sc;
//         if ((sc = imx[i-1][k] + hmm->tsc[TII][k]) > imx[i][k])
//           imx[i][k] = sc;
//         if (hmm->isc[dsq[i]][k] != -INFTY) imx[i][k] += hmm->isc[dsq[i]][k];
//         else                                    imx[i][k] = -INFTY;
//       }
//     }

//     /* Now the special states. Order is important here.
//      * remember, C and J emissions are zero score by definition,
//      */
//     /* N state */
//     xmx[i][XMN] = -INFTY;
//     if ((sc = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
//       xmx[i][XMN] = sc;

//     /* E state */
//     xmx[i][XME] = -INFTY;
//     for (k = 1; k <= hmm->M; k++)
//       if ((sc =  mmx[i][k] + hmm->esc[k]) > xmx[i][XME])
//         xmx[i][XME] = sc;
//     /* J state */
//     xmx[i][XMJ] = -INFTY;
//     if ((sc = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
//       xmx[i][XMJ] = sc;
//     if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
//       xmx[i][XMJ] = sc;

//     /* B state */
//     xmx[i][XMB] = -INFTY;
//     if ((sc = xmx[i][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
//       xmx[i][XMB] = sc;
//     if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
//       xmx[i][XMB] = sc;

//     /* C state */
//     xmx[i][XMC] = -INFTY;
//     if ((sc = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
//       xmx[i][XMC] = sc;
//     if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
//       xmx[i][XMC] = sc;
//   }
//   /* T state (not stored) */
//   sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

//   if (ret_tr != NULL) {
//     P7ViterbiTrace(hmm, dsq, L, mx, &tr);
//     *ret_tr = tr;
//   }

//   return Scorify(sc);    /* the total Viterbi score. */
// }


void
P7ViterbiTrace(
  struct plan7_s *hmm,
  unsigned char *dsq,
  int N,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr
){
  struct p7trace_s *tr;
  int curralloc;    /* current allocated length of trace */
  int tpos;      /* position in trace */
  int i;      /* position in seq (1..N) */
  int k = 0;      /* position in model (1..M) */
  int **xmx, **mmx, **imx, **dmx;
  int sc;      /* temp var for pre-emission score */

  /* Overallocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + N is minimum trace;
   * add N more as buffer.
   */
  curralloc = N * 2 + 6;
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
  i    = N;      /* current i (seq pos) we're trying to assign */

  /* Traceback
   */
  while (tr->statetype[tpos-1] != STS) {
    switch (tr->statetype[tpos-1]) {
    case STM:      /* M connects from i-1,k-1, or B */
      sc = mmx[i+1][k+1] - hmm->msc[dsq[i+1]][k+1];
      if (sc <= -INFTY) {
        P7FreeTrace(tr);
        *ret_tr = NULL;
        return;
      } else if (sc == xmx[i][XMB] + hmm->bsc[k+1]) {
        /* Check for wing unfolding */
        if (Prob2Score(hmm->begin[k+1], hmm->p1) + 1 * INTSCALE <= hmm->bsc[k+1])
          while (k > 0) {
            tr->statetype[tpos] = STD;
            tr->nodeidx[tpos]   = k--;
            tr->pos[tpos]       = 0;
            tpos++;
            if (tpos == curralloc) {
              /* grow trace if necessary  */
              curralloc += N;
              P7ReallocTrace(tr, curralloc);
            }
          }

        tr->statetype[tpos] = STB;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      } else if (sc == mmx[i][k] + hmm->tsc[TMM][k]) {
        tr->statetype[tpos] = STM;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = i--;
      } else if (sc == imx[i][k] + hmm->tsc[TIM][k]) {
        tr->statetype[tpos] = STI;
        tr->nodeidx[tpos]   = k;
        tr->pos[tpos]       = i--;
      } else if (sc == dmx[i][k] + hmm->tsc[TDM][k]) {
        tr->statetype[tpos] = STD;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = 0;
      } else
        Die("traceback failed");
      break;

    case STD:      /* D connects from M,D */
      if (dmx[i][k+1] <= -INFTY) {
        P7FreeTrace(tr);
        *ret_tr = NULL;
        return;
      } else if (dmx[i][k+1] == mmx[i][k] + hmm->tsc[TMD][k]) {
        tr->statetype[tpos] = STM;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = i--;
      } else if (dmx[i][k+1] == dmx[i][k] + hmm->tsc[TDD][k]) {
        tr->statetype[tpos] = STD;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = 0;
      } else Die("traceback failed");
      break;

    case STI:      /* I connects from M,I */
      sc = imx[i+1][k] - hmm->isc[dsq[i+1]][k];
      if (sc <= -INFTY) {
        P7FreeTrace(tr);
        *ret_tr = NULL;
        return;
      } else if (sc == mmx[i][k] + hmm->tsc[TMI][k]) {
        tr->statetype[tpos] = STM;
        tr->nodeidx[tpos]   = k--;
        tr->pos[tpos]       = i--;
      } else if (sc == imx[i][k] + hmm->tsc[TII][k]) {
        tr->statetype[tpos] = STI;
        tr->nodeidx[tpos]   = k;
        tr->pos[tpos]       = i--;
      } else Die("traceback failed");
      break;

    case STN:      /* N connects from S, N */
      if (i == 0 && xmx[i][XMN] == 0) {
        tr->statetype[tpos] = STS;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      } else if (i > 0 && xmx[i+1][XMN] == xmx[i][XMN] + hmm->xsc[XTN][LOOP]) {
        tr->statetype[tpos] = STN;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;    /* note convention adherence:  */
        tr->pos[tpos-1]     = i--;  /* first N doesn't emit        */
      } else Die("traceback failed");
      break;

    case STB:      /* B connects from N, J */
      if (xmx[i][XMB] <= -INFTY) {
        P7FreeTrace(tr);
        *ret_tr = NULL;
        return;
      } else if (xmx[i][XMB] == xmx[i][XMN] + hmm->xsc[XTN][MOVE]) {
        tr->statetype[tpos] = STN;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      } else if (xmx[i][XMB] == xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) {
        tr->statetype[tpos] = STJ;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
      }

      else Die("traceback failed");
      break;

    case STE:      /* E connects from any M state. k set here */
      if (xmx[i][XME] <= -INFTY) {
        P7FreeTrace(tr);
        *ret_tr = NULL;
        return;
      }
      for (k = hmm->M; k >= 1; k--)
        if (xmx[i][XME] == mmx[i][k] + hmm->esc[k]) {
          /* check for wing unfolding */
          if (Prob2Score(hmm->end[k], 1.) + 1*INTSCALE <=  hmm->esc[k]) {
            int dk;    /* need a tmp k while moving thru delete wing */
            for (dk = hmm->M; dk > k; dk--) {
              tr->statetype[tpos] = STD;
              tr->nodeidx[tpos]   = dk;
              tr->pos[tpos]       = 0;
              tpos++;
              if (tpos == curralloc) {
                /* grow trace if necessary  */
                curralloc += N;
                P7ReallocTrace(tr, curralloc);
              }
            }
          }

          tr->statetype[tpos] = STM;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = i--;
          break;
        }
      if (k < 0) Die("traceback failed");
      break;

    case STC:      /* C comes from C, E */
      if (xmx[i][XMC] <= -INFTY) {
        P7FreeTrace(tr);
        *ret_tr = NULL;
        return;
      } else if (xmx[i][XMC] == xmx[i-1][XMC] + hmm->xsc[XTC][LOOP]) {
        tr->statetype[tpos] = STC;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;    /* note convention adherence: */
        tr->pos[tpos-1]     = i--;  /* first C doesn't emit       */
      } else if (xmx[i][XMC] == xmx[i][XME] + hmm->xsc[XTE][MOVE]) {
        tr->statetype[tpos] = STE;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0; /* E is a nonemitter */
      }

      else Die("Traceback failed.");
      break;

    case STJ:      /* J connects from E, J */
      if (xmx[i][XMJ] <= -INFTY) {
        P7FreeTrace(tr);
        *ret_tr = NULL;
        return;
      } else if (xmx[i][XMJ] == xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP]) {
        tr->statetype[tpos] = STJ;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;    /* note convention adherence: */
        tr->pos[tpos-1]     = i--;  /* first J doesn't emit       */
      } else if (xmx[i][XMJ] == xmx[i][XME] + hmm->xsc[XTE][LOOP]) {
        tr->statetype[tpos] = STE;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0; /* E is a nonemitter */
      }

      else Die("Traceback failed.");
      break;

    default:
      Die("traceback failed");

    } /* end switch over statetype[tpos-1] */

    tpos++;
    if (tpos == curralloc) {
      /* grow trace if necessary  */
      curralloc += N;
      P7ReallocTrace(tr, curralloc);
    }

  } /* end traceback, at S state; tpos == tlen now */
  tr->tlen = tpos;
  P7ReverseTrace(tr);
  *ret_tr = tr;
}


float
P7SmallViterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr
){
  struct p7trace_s *ctr;        /* collapsed trace of optimal parse */
  struct p7trace_s *tr;         /* full trace of optimal alignment */
  struct p7trace_s **tarr;      /* trace array */
  int   ndom;      /* number of subsequences */
  int   i;      /* counter over domains   */
  int   pos;      /* position in sequence */
  int   tpos;      /* position in trace */
  int   tlen;      /* length of full trace   */
  int   totlen;                 /* length of L matched by model (as opposed to N/C/J) */
  float sc;      /* score of optimal alignment */
  int   t2;      /* position in a subtrace */

  /* Step 1. Call P7ParsingViterbi to calculate an optimal parse
   *         of the sequence into single-hit subsequences; this parse
   *         is returned in a "collapsed" trace
   */
  sc = P7ParsingViterbi(dsq, L, hmm, &ctr);

  /* If we don't want full trace, we're done;
   * also, if parsing viterbi returned a NULL trace we're done. */
  if (ctr == NULL || ret_tr == NULL) {
    P7FreeTrace(ctr);
    return sc;
  }

  /* Step 2. Call either P7Viterbi or P7WeeViterbi on each subsequence
   *         to recover a full traceback of each, collecting them in
   *         an array.
   */
  ndom = ctr->tlen/2 - 1;
  tarr = MallocOrDie(sizeof(struct p7trace_s *) * ndom);
  tlen = totlen = 0;
  for (i = 0; i < ndom; i++) {
    int   sqlen;      /* length of a subsequence */
    sqlen = ctr->pos[i*2+2] - ctr->pos[i*2+1];   /* length of subseq */

    if (P7ViterbiSpaceOK(sqlen, hmm->M, mx)) {
      P7Viterbi(dsq + ctr->pos[i*2+1], sqlen, hmm, mx, &(tarr[i]));
    } else if (sqlen == 1) {
      /* xref bug#h30. P7WeeViterbi() can't take L=1. This
       is a hack to work around the problem, which is rare.
       Attempts to use our main dp mx will violate our
       RAMLIMIT guarantee, so allocate a tiny linear one. */
      struct dpmatrix_s *tiny;
      tiny = CreatePlan7Matrix(1, hmm->M, 0, 0);
      P7Viterbi(dsq + ctr->pos[i*2+1], sqlen, hmm, tiny, &(tarr[i]));
      FreePlan7Matrix(tiny);
    } else {
      P7WeeViterbi(dsq + ctr->pos[i*2+1], sqlen, hmm, &(tarr[i]));
    }

    tlen  += tarr[i]->tlen - 4; /* not counting S->N,...,C->T */
    totlen += sqlen;
  }

  /* Step 3. Compose the subtraces into one big final trace.
   *         This is wasteful because we're going to TraceDecompose()
   *         it again in both hmmsearch and hmmpfam to look at
   *         individual domains; but we do it anyway so the P7SmallViterbi
   *         interface looks exactly like the P7Viterbi interface. Maybe
   *         long traces shouldn't include all the N/J/C states anyway,
   *         since they're unambiguously implied.
   */

  /* Calculate total trace len and alloc;
   * nonemitting SNCT + nonemitting J's + emitting NJC
   */
  tlen += 4 + (ndom-1) + (L-totlen);
  P7AllocTrace(tlen, &tr);
  tr->tlen = tlen;

  /* Add N-terminal trace framework
   */
  tr->statetype[0] = STS;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STN;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
  /* add implied N's */
  for (pos = 1; pos <= ctr->pos[1]; pos++) {
    tr->statetype[tpos] = STN;
    tr->nodeidx[tpos]   = 0;
    tr->pos[tpos]       = pos;
    tpos++;
  }

  /* Add each subseq trace in, with its appropriate
   * sequence offset derived from the collapsed trace
   */
  for (i = 0; i < ndom; i++) {
    /* skip SN, CT framework at ends */
    for (t2 = 2; t2 < tarr[i]->tlen-2; t2++) {
      tr->statetype[tpos] = tarr[i]->statetype[t2];
      tr->nodeidx[tpos]   = tarr[i]->nodeidx[t2];
      if (tarr[i]->pos[t2] > 0)
        tr->pos[tpos]       = tarr[i]->pos[t2] + ctr->pos[i*2+1];
      else
        tr->pos[tpos]       = 0;
      tpos++;
    }
    /* add nonemitting J or C */
    tr->statetype[tpos] = (i == ndom-1) ? STC : STJ;
    tr->nodeidx[tpos]   = 0;
    tr->pos[tpos]       = 0;
    tpos++;
    /* add implied emitting J's */
    if (i != ndom-1)
      for (pos = ctr->pos[i*2+2]+1; pos <= ctr->pos[(i+1)*2+1]; pos++) {
        tr->statetype[tpos] = STJ;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = pos;
        tpos++;
      }
  }

  /* add implied C's */
  for (pos = ctr->pos[ndom*2]+1; pos <= L; pos++) {
    tr->statetype[tpos] = STC;
    tr->nodeidx[tpos]   = 0;
    tr->pos[tpos]       = pos;
    tpos++;
  }
  /* add terminal T */
  tr->statetype[tpos] = STT;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0;
  tpos++;

  for (i = 0; i < ndom; i++) P7FreeTrace(tarr[i]);
  free(tarr);
  P7FreeTrace(ctr);

  *ret_tr = tr;
  return sc;
}


float
P7ParsingViterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr
){
  struct dpmatrix_s *mx;        /* two rows of score matrix */
  struct dpmatrix_s *tmx;       /* two rows of misused score matrix: traceback ptrs */
  struct p7trace_s  *tr;        /* RETURN: collapsed traceback */
  int  **xmx, **mmx, **dmx, **imx; /* convenience ptrs to score matrix */
  int  **xtr, **mtr, **dtr, **itr; /* convenience ptrs to traceback pointers */
  int   *btr, *etr;             /* O(L) trace ptrs for B, E state pts in seq */
  int    sc;      /* integer score of optimal alignment  */
  int    i,k,tpos;    /* index for seq, model, trace position */
  int    cur;    /* indices for rolling dp matrix */
  int    curralloc;    /* size of allocation for tr */


  /* Alloc a DP matrix and traceback pointers, two rows each, O(M).
   * Alloc two O(L) arrays to trace back through the sequence thru B and E.
   */
  mx  = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  tmx = AllocPlan7Matrix(2, hmm->M, &xtr, &mtr, &itr, &dtr);
  btr = MallocOrDie(sizeof(int) * (L+1));
  etr = MallocOrDie(sizeof(int) * (L+1));

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;                         /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  btr[0]      = 0;
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  etr[0]      = -1;
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull. Rolling index trick. Trace ptr propagation trick.
   * Note some slightly wasteful boundary conditions:
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   *
   * Notes on traceback pointer propagation.
   *  - In the path B->E, we propagate the i that B was aligned to in the optimal
   *    alignment, via mtr, dtr, and itr.
   *  - When we reach an E, we record the i of the B it started from in etr.
   *  - In a looping path E->J...->B or terminal path E->C...->T, we propagate
   *    the i that E was aligned to in the optimal alignment via xtr[][XMC]
   *    and xtr[][XMJ].
   *  - When we enter B, we record the i of the best previous E, or 0 if there
   *    isn't one, in btr.
   */
  cur = 1;
  for (i = 1; i <= L; i++) {
    int    prv;    /* indices for rolling dp matrix */
    cur = i % 2;
    prv = !cur;//TODO: BUG: this was never properly initialized

    mmx[cur][0] = imx[cur][0] = dmx[cur][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
      /* match state */
      mmx[cur][k] = -INFTY;
      if ((sc = mmx[prv][k-1] + hmm->tsc[TMM][k-1]) > -INFTY) {
        mmx[cur][k] = sc;
        mtr[cur][k] = mtr[prv][k-1];
      }
      if ((sc = imx[prv][k-1] + hmm->tsc[TIM][k-1]) > mmx[cur][k]) {
        mmx[cur][k] = sc;
        mtr[cur][k] = itr[prv][k-1];
      }
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k]) {
        mmx[cur][k] = sc;
        mtr[cur][k] = i-1;
      }
      if ((sc = dmx[prv][k-1] + hmm->tsc[TDM][k-1]) > mmx[cur][k]) {
        mmx[cur][k] = sc;
        mtr[cur][k] = dtr[prv][k-1];
      }
      if (hmm->msc[dsq[i]][k] != -INFTY)
        mmx[cur][k] += hmm->msc[dsq[i]][k];
      else
        mmx[cur][k] = -INFTY;

      /* delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > -INFTY) {
        dmx[cur][k] = sc;
        dtr[cur][k] = mtr[cur][k-1];
      }
      if ((sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k]) {
        dmx[cur][k] = sc;
        dtr[cur][k] = dtr[cur][k-1];
      }

      /* insert state */
      if (k < hmm->M) {
        imx[cur][k] = -INFTY;
        if ((sc = mmx[prv][k] + hmm->tsc[TMI][k]) > -INFTY) {
          imx[cur][k] = sc;
          itr[cur][k] = mtr[prv][k];
        }
        if ((sc = imx[prv][k] + hmm->tsc[TII][k]) > imx[cur][k]) {
          imx[cur][k] = sc;
          itr[cur][k] = itr[prv][k];
        }
        if (hmm->isc[dsq[i]][k] != -INFTY)
          imx[cur][k] += hmm->isc[dsq[i]][k];
        else
          imx[cur][k] = -INFTY;
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
    /* N state */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[cur][XMN] = sc;
    /* E state */
    xmx[cur][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME]) {
        xmx[cur][XME] = sc;
        etr[i] = mtr[cur][k];
      }
    /* J state */
    xmx[cur][XMJ] = -INFTY;
    if ((sc = xmx[prv][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY) {
      xmx[cur][XMJ] = sc;
      xtr[cur][XMJ] = xtr[prv][XMJ];
    }
    if ((sc = xmx[cur][XME]   + hmm->xsc[XTE][LOOP]) > xmx[cur][XMJ]) {
      xmx[cur][XMJ] = sc;
      xtr[cur][XMJ] = i;
    }
    /* B state */
    xmx[cur][XMB] = -INFTY;
    if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY) {
      xmx[cur][XMB] = sc;
      btr[i] = 0;
    }
    if ((sc = xmx[cur][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[cur][XMB]) {
      xmx[cur][XMB] = sc;
      btr[i] = xtr[cur][XMJ];
    }
    /* C state */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY) {
      xmx[cur][XMC] = sc;
      xtr[cur][XMC] = xtr[prv][XMC];
    }
    if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC]) {
      xmx[cur][XMC] = sc;
      xtr[cur][XMC] = i;
    }
  }
  /* T state (not stored) */
  sc = xmx[cur][XMC] + hmm->xsc[XTC][MOVE];

  /*****************************************************************
   * Collapsed traceback stage.
   * xtr[L%2][XMC] contains the position j of the previous E
   * etr[j]        contains the position i of the previous B
   * btr[i]        contains the position j of the previous E, or 0
   * continue until btr[i] = 0.
   *****************************************************************/

  curralloc = 2;    /* minimum: no hits */
  P7AllocTrace(curralloc, &tr);

  /* Init of collapsed trace. Back to front; we ReverseTrace() later.
   */
  tpos = 0;
  tr->statetype[tpos] = STT;
  tr->pos[tpos]       = 0;
  i                   = xtr[L%2][XMC];
  while (i > 0) {
    curralloc += 2;
    P7ReallocTrace(tr, curralloc);

    tpos++;
    tr->statetype[tpos] = STE;
    tr->pos[tpos] = i;
    i = etr[i];

    tpos++;
    tr->statetype[tpos] = STB;
    tr->pos[tpos] = i;
    i = btr[i];
  }

  tpos++;
  tr->statetype[tpos] = STS;
  tr->pos[tpos]       = 0;
  tr->tlen = tpos + 1;
  P7ReverseTrace(tr);

  FreePlan7Matrix(mx);
  FreePlan7Matrix(tmx);
  free(btr);
  free(etr);

  *ret_tr = tr;
  return Scorify(sc);
}


float
P7WeeViterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr
){
  struct p7trace_s *tr;         /* RETURN: traceback */
  int          *kassign;        /* 0..L+1, alignment of seq positions to model nodes */
  char         *tassign;        /* 0..L+1, alignment of seq positions to state types */
  int          *endlist;        /* stack of end points on sequence to work on */
  int          *startlist;      /* stack of start points on sequence to work on */
  int          lpos;            /* position in endlist, startlist */
  float        ret_sc = 0.0;    /* optimal score over complete seq */
  int          tlen;    /* length needed for trace */
  int          i, k, tpos;  /* index in sequence, model, trace */

  /* Someday, reexamine impl of get_wee_midpoint, and remove this
   * L>1 limitation. (xref bug #h30).
   */
  if (L==1) Die("P7WeeViterbi() cannot accept L=1 subsequence.\n");

  /* Initialize.
   */
  kassign   = MallocOrDie (sizeof(int) * (L+1));
  tassign   = MallocOrDie (sizeof(char)* (L+1));
  endlist   = MallocOrDie (sizeof(int) * (L+1));
  startlist = MallocOrDie (sizeof(int) * (L+1));

  lpos = 0;
  startlist[lpos] = 1;
  endlist[lpos]   = L;
  kassign[1]      = 1;
  kassign[L]      = hmm->M;
  tassign[1]      = STS;  /* temporary boundary condition! will become N or M */
  tassign[L]      = STT;  /* temporary boundary condition! will become M or C */


  /* Recursive divide-and-conquer alignment.
   */
  while (lpos >= 0) {
    //NOTE: k2, t2, and s2 are initialized by reference.
    int k1, k2, k3; /* start, middle, and end in model      */
    char t1, t2, t3;  /* start, middle, and end in state type */
    int s1, s2, s3;  /* start, middle, and end in sequence   */
    /* Pop a segment off the stack */
    s1 = startlist[lpos];
    k1 = kassign[s1];
    t1 = tassign[s1];
    s2 = 0;
    k2 = 0;
    t2 = 0;
    s3 = endlist[lpos];
    k3 = kassign[s3];
    t3 = tassign[s3];
    lpos--;

    /* find optimal midpoint of segment */
    float sc;    /* score of segment optimal alignment */
    sc = get_wee_midpt(hmm, dsq, k1, t1, s1, k3, t3, s3, &k2, &t2, &s2);
    kassign[s2] = k2;
    tassign[s2] = t2;
    /* score is valid on first pass */
    if (t1 == STS && t3 == STT) ret_sc = sc;

    /* push N-terminal segment on stack */
    if (t2 != STN && (s2 - s1 > 1 || (s2 - s1 == 1 && t1 == STS))) {
      lpos++;
      startlist[lpos] = s1;
      endlist[lpos]   = s2;
    }
    /* push C-terminal segment on stack */
    if (t2 != STC && (s3 - s2 > 1 || (s3 - s2 == 1 && t3 == STT))) {
      lpos++;
      startlist[lpos] = s2;
      endlist[lpos]   = s3;
    }

    if (t2 == STN) { //if we see STN midpoint, we know the whole N-term is STN
      for (; s2 >= s1; s2--) {
        kassign[s2] = 1;
        tassign[s2] = STN;
      }
    }
    if (t2 == STC) { //if we see STC midpoint, we know whole C-term is STC
      for (; s2 <= s3; s2++) {
        kassign[s2] = hmm->M;
        tassign[s2] = STC;
      }
    }
  }

  /*****************************************************************
   * Construct a traceback structure from kassign/tassign by interpolating
   * necessary states.
   * Trace allocation is as follows. We clearly need L emitting states.
   * We also need nonemitting states as follows:
   * STS,STN,STB,STE,STC,STT = 6
   * STD: count k2-k1-1 in kassign M->M's
   * Also, count N->M's and M->C's (potential wing unfoldings)...
   *   ...and be careful to check wing unfoldings when there aren't
   *      any emitting N or C flanks! (bugfix, 2.1.1b)
   *****************************************************************/

  tlen = L + 6;
  for (i = 1; i < L; i++) {
    if (tassign[i] == STM && tassign[i+1] == STM)
      tlen += kassign[i+1] - kassign[i] - 1;
    if (tassign[i] == STN && tassign[i+1] == STM)
      tlen += kassign[i+1] - 1;
    if (tassign[i] == STM && tassign[i+1] == STC)
      tlen += hmm->M - kassign[i];
  }
  if (tassign[1] == STM) tlen += kassign[1] - 1;
  if (tassign[L] == STM) tlen += hmm->M - kassign[L];
  P7AllocTrace(tlen, &tr);

  tr->statetype[0] = STS;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STN;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;

  for (i = 1; i <= L; i++) {
    switch(tassign[i]) {
    case STM:
      /* check for first match state */
      if (tr->statetype[tpos-1] == STN) {
        tr->statetype[tpos] = STB;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
        tpos++;
        /* check for wing unfolding */
        if (Prob2Score(hmm->begin[kassign[i]], hmm->p1) + INTSCALE <= hmm->bsc[kassign[i]])
          for (k = 1; k < kassign[i]; k++) {
            tr->statetype[tpos] = STD;
            tr->nodeidx[tpos]   = k;
            tr->pos[tpos]       = 0;
            tpos++;
          }
      }
      /* do the match state itself */
      tr->statetype[tpos] = STM;
      tr->nodeidx[tpos]   = kassign[i];
      tr->pos[tpos]       = i;
      tpos++;
      /* do any deletes necessary 'til next match */
      if (i < L && tassign[i+1] == STM && kassign[i+1] - kassign[i] > 1)
        for (k = kassign[i] + 1; k < kassign[i+1]; k++) {
          tr->statetype[tpos] = STD;
          tr->nodeidx[tpos]   = k;
          tr->pos[tpos]       = 0;
          tpos++;
        }
      /* check for last match state */
      if (i == L || tassign[i+1] == STC) {
        /* check for wing unfolding */
        if (Prob2Score(hmm->end[kassign[i-1]], 1.) + INTSCALE <=  hmm->esc[kassign[i-1]])
          for (k = kassign[i]+1; k <= hmm->M; k++) {
            tr->statetype[tpos] = STD;
            tr->nodeidx[tpos]   = k;
            tr->pos[tpos]       = 0;
            tpos++;
          }
        /* add on the end state */
        tr->statetype[tpos] = STE;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
        tpos++;
        /* and a nonemitting C state */
        tr->statetype[tpos] = STC;
        tr->nodeidx[tpos]   = 0;
        tr->pos[tpos]       = 0;
        tpos++;
      }
      break;

    case STI:
      tr->statetype[tpos] = STI;
      tr->nodeidx[tpos]   = kassign[i];
      tr->pos[tpos]       = i;
      tpos++;
      break;

    case STN:
      tr->statetype[tpos] = STN;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = i;
      tpos++;
      break;

    case STC:
      tr->statetype[tpos] = STC;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = i;
      tpos++;
      break;

    default:
      Die("Bogus state %s", Statetype(tassign[i]));
    }
  }
  /* terminate the trace */
  tr->statetype[tpos] = STT;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0;
  tr->tlen = tpos+1;

  *ret_tr = tr;

  free(kassign);
  free(tassign);
  free(startlist);
  free(endlist);
  return ret_sc;
}


float
get_wee_midpt(
  struct plan7_s *hmm,
  unsigned char *dsq,
  int k1,
  char t1,
  int s1,
  int k3,
  char t3,
  int s3,
  int *ret_k2,
  char *ret_t2,
  int *ret_s2
){
  struct dpmatrix_s *fwd;
  struct dpmatrix_s *bck;
  int        **xmx;             /* convenience ptr into special states */
  int        **mmx;             /* convenience ptr into match states   */
  int        **imx;             /* convenience ptr into insert states  */
  int        **dmx;             /* convenience ptr into delete states  */
  int          k2;
  char         t2;
  int          s2;
  int          cur, nxt;  /* current, next row index (0 or 1)*/
  int          i,k;    /* indices for seq, model */
  int          sc;    /* integer score */
  int          max;    /* maximum integer score */
  int          start;    /* s1 to start at (need, for STS special case) */


  /* Choose our midpoint.
   * Special cases: s1, s3 adjacent and t1 == STS: s2 = s1
   *                s1, s3 adjacent and t3 == STT: s2 = s3
   *                (where we must replace STS, STT eventually)
   */
  s2 = s1 + (s3-s1) / 2;
  if (s3-s1 == 1 && t1 == STS) s2 = s1;
  if (s3-s1 == 1 && t3 == STT) s2 = s3;

  /* STS is a special case. STS aligns to row zero by convention,
   * but we'll be passed s1=1, t1=STS. We have to init on row
   * zero then start DP on row 1.
   */
  start = (t1 == STS) ? 0 : s1;

  /* Allocate our forward two rows.
   * Initialize row zero.
   */
  fwd = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  cur = start%2;
  xmx[cur][XMN] = xmx[cur][XMB] = -INFTY;
  xmx[cur][XME] = xmx[cur][XMC] = -INFTY;
  for (k = k1; k <= k3; k++)
    mmx[cur][k] = imx[cur][k] = dmx[cur][k] = -INFTY;

  /* Where to put our zero for our start point...
   * (only possible to start on an emitting state; J disallowed)
   */
  switch (t1) {
  case STM:
    mmx[cur][k1]  = 0;
    break;
  case STI:
    imx[cur][k1]  = 0;
    break;
  case STN:
    xmx[cur][XMN] = 0;
    break;
  case STC:
    xmx[cur][XMC] = 0;
    break;
  case STS:
    xmx[cur][XMN] = 0;
    break;
  default:
    Die("you can't init get_wee_midpt with a %s\n", Statetype(t1));
  }

  /* Still initializing.
   * Deal with pulling horizontal matrix moves in initial row.
   * These are any transitions to nonemitters:
   *    STM-> E, D
   *    STI-> none
   *    STN-> B
   *    STC-> (T, but we never observe this in the forward pass of a d&c)
   *    STE-> C
   *    STS-> (N, already implied by setting xmx[cur][XMN] = 0)
   *    STB-> M
   */
  if (t1 == STM) {
    for (k = k1+1; k <= k3; k++) { /* transits into STD */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > -INFTY)
        dmx[cur][k] = sc;
      if ((sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k])
        dmx[cur][k] = sc;
    }
    /* transit into STE */
    xmx[cur][XME] = -INFTY;
    if ((sc = mmx[cur][k1] + hmm->esc[k1]) > -INFTY)
      xmx[cur][XME] = sc;
  }

  /* transit into STB from STN */
  xmx[cur][XMB] = -INFTY;
  if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
    xmx[cur][XMB] = sc;
  /* transit into STC from STE */
  xmx[cur][XMC] = -INFTY;
  if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > -INFTY)
    xmx[cur][XMC] = sc;

  /* Done initializing.
   * Start recursive DP; sweep forward to chosen s2 midpoint. Done as a pull.
   */
  for (i = start+1; i <= s2; i++) {
    int prv;  /* previous row index (0 or 1)*/
    cur = i % 2;
    prv = !cur;

    mmx[cur][k1] = imx[cur][k1] = dmx[cur][k1] = -INFTY;

    /* Insert state in column k1, and B->M transition in k1.
     */
    if (k1 < hmm->M) {
      imx[cur][k1] = -INFTY;
      if ((sc = mmx[prv][k1] + hmm->tsc[TMI][k1]) > -INFTY)
        imx[cur][k1] = sc;
      if ((sc = imx[prv][k1] + hmm->tsc[TII][k1]) > imx[cur][k1])
        imx[cur][k1] = sc;
      if (hmm->isc[dsq[i]][k1] != -INFTY)
        imx[cur][k1] += hmm->isc[dsq[i]][k1];
      else
        imx[cur][k1] = -INFTY;
    }
    if ((sc = xmx[prv][XMB] + hmm->bsc[k1]) > -INFTY)
      mmx[cur][k1] = sc;
    if (hmm->msc[dsq[i]][k1] != -INFTY)
      mmx[cur][k1] += hmm->msc[dsq[i]][k1];
    else
      mmx[cur][k1] = -INFTY;

    /* Main chunk of recursion across model positions
     */
    for (k = k1+1; k <= k3; k++) {
      /* match state */
      mmx[cur][k]  = -INFTY;
      if ((sc = mmx[prv][k-1] + hmm->tsc[TMM][k-1]) > -INFTY)
        mmx[cur][k] = sc;
      if ((sc = imx[prv][k-1] + hmm->tsc[TIM][k-1]) > mmx[cur][k])
        mmx[cur][k] = sc;
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k])
        mmx[cur][k] = sc;
      if ((sc = dmx[prv][k-1] + hmm->tsc[TDM][k-1]) > mmx[cur][k])
        mmx[cur][k] = sc;
      if (hmm->msc[dsq[i]][k] != -INFTY)
        mmx[cur][k] += hmm->msc[dsq[i]][k];
      else
        mmx[cur][k] = -INFTY;

      /* delete state */
      dmx[cur][k] = -INFTY;
      if (k < hmm->M) {
        if ((sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > -INFTY)
          dmx[cur][k] = sc;
        if ((sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k])
          dmx[cur][k] = sc;
      }

      /* insert state */
      imx[cur][k] = -INFTY;
      if (k < hmm->M) {
        if ((sc = mmx[prv][k] + hmm->tsc[TMI][k]) > -INFTY)
          imx[cur][k] = sc;
        if ((sc = imx[prv][k] + hmm->tsc[TII][k]) > imx[cur][k])
          imx[cur][k] = sc;
        if (hmm->isc[dsq[i]][k] != -INFTY)
          imx[cur][k] += hmm->isc[dsq[i]][k];
        else
          imx[cur][k] = -INFTY;
      }
    }
    /* N state */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[cur][XMN] = sc;
    /* E state */
    xmx[cur][XME] = -INFTY;
    for (k = k1; k <= k3 && k <= hmm->M; k++)
      if ((sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME])
        xmx[cur][XME] = sc;
    /* B state */
    xmx[cur][XMB] = -INFTY;
    if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[cur][XMB] = sc;
    /* C state */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[cur][XMC] = sc;
    if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC])
      xmx[cur][XMC] = sc;
  }

  /* Row s2%2 in fwd matrix now contains valid scores from s1 (start) to s2,
   * with J transitions disallowed (no cycles through model).
   */

  /*****************************************************************
   * Backwards pass.
   *****************************************************************/

  /* Allocate our backwards two rows. Init last row.
   */
  bck = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  nxt = s3%2;
  xmx[nxt][XMN] = xmx[nxt][XMB] = -INFTY;
  xmx[nxt][XME] = xmx[nxt][XMC] = -INFTY;
  for (k = k1; k <= k3 + 1; k++)
    mmx[nxt][k] = imx[nxt][k] = dmx[nxt][k] = -INFTY;
  cur = !nxt;
  mmx[cur][k3+1] = imx[cur][k3+1] = dmx[cur][k3+1] = -INFTY;

  /* Where to put the zero for our end point on last row.
   */
  switch (t3) {
  case STM:
    mmx[nxt][k3]  = 0;
    break;
  case STI:
    imx[nxt][k3]  = 0;
    break;
  case STN:
    xmx[nxt][XMN] = 0;
    break;
  case STC:
    xmx[nxt][XMC] = 0;
    break;   /* must be an emitting C */
  case STT:
    xmx[nxt][XMC] = hmm->xsc[XTC][MOVE];
    break; /* C->T implied */
  default:
    Die("you can't init get_wee_midpt with a %s\n", Statetype(t3));
  }

  /* Still initializing.
   * In the case t3==STT, there are a few horizontal moves possible
   * on row s3, because STT isn't an emitter. All other states are
   * emitters, so their connections have to be to the previous row s3-1.
   */
  if (t3 == STT) { // E->C
    xmx[nxt][XME] = xmx[nxt][XMC] + hmm->xsc[XTE][MOVE];
    // M->E
    for (k = k3; k >= k1; k--) {
      mmx[nxt][k] = xmx[nxt][XME] + hmm->esc[k];
      if (s3 != s2)
        mmx[nxt][k] += hmm->msc[dsq[s3]][k];
    }
  }

  /* Start recursive DP; sweep backwards to chosen s2 midpoint.
   * Done as a pull. M, I scores at current row do /not/ include
   * emission scores. Be careful of integer underflow.
   */
  for (i = s3-1; i >= s2; i--) {
    /* note i < L, so i+1 is always a legal index */
    cur = i%2;
    nxt = !cur;
    /* C pulls from C (T is special cased) */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[nxt][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[cur][XMC] = sc;
    /* B pulls from M's */
    xmx[cur][XMB] = -INFTY;
    for (k = k1; k <= k3; k++)
      if ((sc = mmx[nxt][k] + hmm->bsc[k]) > xmx[cur][XMB])
        xmx[cur][XMB] = sc;
    /* E pulls from C (J disallowed) */
    xmx[cur][XME] = -INFTY;
    if ((sc = xmx[cur][XMC] + hmm->xsc[XTE][MOVE]) > -INFTY)
      xmx[cur][XME] = sc;
    /* N pulls from B, N */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[cur][XMB] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[cur][XMN] = sc;
    if ((sc = xmx[nxt][XMN] + hmm->xsc[XTN][LOOP]) > xmx[cur][XMN])
      xmx[cur][XMN] = sc;

    /* Main recursion across model
     */
    for (k = k3; k >= k1; k--)  {
      /* special case k == M */
      if (k == hmm->M) {
        mmx[cur][k] = xmx[cur][XME]; /* p=1 transition to E by definition */
        dmx[cur][k] = -INFTY;  /* doesn't exist */
        imx[cur][k] = -INFTY;  /* doesn't exist */
        if (i != s2)
          mmx[cur][k] += hmm->msc[dsq[i]][k];
        continue;
      }      /* below this k < M, so k+1 is a legal index */

      /* pull into match state */
      mmx[cur][k] = -INFTY;
      if ((sc = xmx[cur][XME] + hmm->esc[k]) > -INFTY)
        mmx[cur][k] = sc;
      if ((sc = mmx[nxt][k+1] + hmm->tsc[TMM][k]) > mmx[cur][k])
        mmx[cur][k] = sc;
      if ((sc = imx[nxt][k] + hmm->tsc[TMI][k]) > mmx[cur][k])
        mmx[cur][k] = sc;
      if ((sc = dmx[cur][k+1] + hmm->tsc[TMD][k]) > mmx[cur][k])
        mmx[cur][k] = sc;
      if (i != s2)
        mmx[cur][k] += hmm->msc[dsq[i]][k];

      /* pull into delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[nxt][k+1] + hmm->tsc[TDM][k]) > -INFTY)
        dmx[cur][k] = sc;
      if ((sc = dmx[cur][k+1] + hmm->tsc[TDD][k]) > dmx[cur][k])
        dmx[cur][k] = sc;
      /* pull into insert state */
      imx[cur][k] = -INFTY;
      if ((sc = mmx[nxt][k+1] + hmm->tsc[TIM][k]) > -INFTY)
        imx[cur][k] = sc;
      if ((sc = imx[nxt][k] + hmm->tsc[TII][k]) > imx[cur][k])
        imx[cur][k] = sc;
      if (i != s2)
        imx[cur][k] += hmm->isc[dsq[i]][k];
    }
  }

  /*****************************************************************
   * DP complete; we have both forward and backward passes. Now we
   * look across the s2 row and find the optimal emitting state.
   *****************************************************************/

  cur = s2%2;
  max = -INFTY;
  for (k = k1; k <= k3; k++) {
    if ((sc = fwd->mmx[cur][k] + bck->mmx[cur][k]) > max) {
      k2 = k;
      t2 = STM;
      max = sc;
    }
    if ((sc = fwd->imx[cur][k] + bck->imx[cur][k]) > max) {
      k2 = k;
      t2 = STI;
      max = sc;
    }
  }
  if ((sc = fwd->xmx[cur][XMN] + bck->xmx[cur][XMN]) > max) {
    k2 = 1;
    t2 = STN;
    max = sc;
  }
  if ((sc = fwd->xmx[cur][XMC] + bck->xmx[cur][XMC]) > max) {
    k2 = hmm->M;
    t2 = STC;
    max = sc;
  }

  /*****************************************************************
   * Garbage collection, return.
   *****************************************************************/

  FreePlan7Matrix(fwd);
  FreePlan7Matrix(bck);
  *ret_k2 = k2;
  *ret_t2 = t2;
  *ret_s2 = s2;
  return Scorify(max);
}


struct p7trace_s*
P7ViterbiAlignAlignment(
  MSA *msa,
  struct plan7_s *hmm
){
  struct dpmatrix_s *mx;        /* Viterbi calculation lattice (two rows) */
  struct dpshadow_s *tb;        /* shadow matrix of traceback pointers */
  struct p7trace_s  *tr;        /* RETURN: traceback */
  int  **xmx, **mmx, **imx, **dmx;
  char **xtb, **mtb, **itb, **dtb;
  float **con;                  /* [1..alen][0..Alphabet_size-1], consensus counts */
  float  *mocc;                 /* fractional occupancy of a column; used to weight transitions */
  size_t  i;      /* counter for columns */
  int     k;      /* counter for model positions */
  // idx      /* counter for seqs */
  int     sym;      /* counter for alphabet symbols */
  int     sc;      /* temp variable for holding score */
  float   denom;    /* total weight of seqs; used to "normalize" counts */

  /* The "consensus" is a counts matrix, [1..alen][0..Alphabet_size-1].
   * Gaps are not counted explicitly, but columns with lots of gaps get
   * less total weight because they have fewer counts.
   */
  /* allocation */
  con  = MallocOrDie(sizeof(float *) * (msa->alen+1));
  mocc = MallocOrDie(sizeof(float)   * (msa->alen+1));
  for (i = 1; i <= msa->alen; i++) {
    con[i] = MallocOrDie(sizeof(float) * Alphabet_size);
    FSet(con[i], Alphabet_size, 0.0);
  }
  mocc[0] = -9999.;
  /* initialization */
  /* note: aseq is off by one, 0..alen-1 */
  /* "normalized" to have a max total count of 1 per col */
  denom = FSum(msa->wgt, msa->nseq);
  for (i = 1; i <= msa->alen; i++) {
    for (size_t idx = 0; idx < msa->nseq; idx++)
      if (! isgap(msa->aseq[idx][i-1]))
        P7CountSymbol(con[i], SYMIDX(msa->aseq[idx][i-1]), msa->wgt[idx]);
    FScale(con[i], Alphabet_size, 1./denom);
    mocc[i] = FSum(con[i], Alphabet_size);
  }

  /* Allocate a DP matrix with 2 rows, 0..M columns,
   * and a shadow matrix with 0,1..alen rows, 0..M columns.
   */
  mx = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  tb = AllocShadowMatrix(msa->alen+1, hmm->M, &xtb, &mtb, &itb, &dtb);

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;                         /* S->N, p=1            */
  xtb[0][XMN] = STS;
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xtb[0][XMB] = STN;
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  tb->esrc[0] = 0;
  xtb[0][XMC] = xtb[0][XMJ] = STBOGUS;
  for (k = 0; k <= hmm->M; k++) {
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */
    mtb[0][k] = itb[0][k] = dtb[0][k] = STBOGUS;
  }

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= msa->alen; i++) {
    int     cur, prv;
    cur = i % 2;
    prv = ! cur;

    mmx[cur][0] = imx[cur][0] = dmx[cur][0] = -INFTY;
    mtb[i][0]   = itb[i][0]   = dtb[i][0]   = STBOGUS;

    for (k = 1; k <= hmm->M; k++) {
      /* match state */
      mmx[cur][k]  = -INFTY;
      mtb[i][k]    = STBOGUS;
      if (mmx[prv][k-1] > -INFTY && hmm->tsc[TMM][k-1] > -INFTY &&
          (sc = mmx[prv][k-1] + hmm->tsc[TMM][k-1]) > mmx[cur][k]) {
        mmx[cur][k] = sc;
        mtb[i][k] = STM;
      }
      if (imx[prv][k-1] > -INFTY && hmm->tsc[TIM][k-1] > -INFTY &&
          (sc = imx[prv][k-1] + hmm->tsc[TIM][k-1] * mocc[i-1]) > mmx[cur][k]) {
        mmx[cur][k] = sc;
        mtb[i][k] = STI;
      }
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k]) {
        mmx[cur][k] = sc;
        mtb[i][k] = STB;
      }
      if (dmx[prv][k-1] > -INFTY && hmm->tsc[TDM][k-1] > -INFTY &&
          (sc = dmx[prv][k-1] + hmm->tsc[TDM][k-1]) > mmx[cur][k]) {
        mmx[cur][k] = sc;
        mtb[i][k] = STD;
      }
      /* average over "consensus" sequence */
      for (sym = 0; sym < Alphabet_size; sym++) {
        if (con[i][sym] > 0 && hmm->msc[sym][k] == -INFTY) {
          mmx[cur][k] = -INFTY;
          break;
        }
        mmx[cur][k] += hmm->msc[sym][k] * con[i][sym];
      }

      /* delete state */
      dmx[cur][k] = -INFTY;
      dtb[i][k]   = STBOGUS;
      if (mmx[cur][k-1] > -INFTY && hmm->tsc[TMD][k-1] > -INFTY &&
          (sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > dmx[cur][k]) {
        dmx[cur][k] = sc;
        dtb[i][k] = STM;
      }
      if (dmx[cur][k-1] > -INFTY && hmm->tsc[TDD][k-1] > -INFTY &&
          (sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k]) {
        dmx[cur][k] = sc;
        dtb[i][k] = STD;
      }

      /* insert state */
      if (k < hmm->M) {
        imx[cur][k] = -INFTY;
        itb[i][k]   = STBOGUS;
        if (mmx[prv][k] > -INFTY && hmm->tsc[TMI][k] > -INFTY &&
            (sc = mmx[prv][k] + hmm->tsc[TMI][k] * mocc[i]) > imx[cur][k]) {
          imx[cur][k] = sc;
          itb[i][k] = STM;
        }
        if (imx[prv][k] > -INFTY && hmm->tsc[TII][k] > -INFTY &&
            (sc = imx[prv][k] + hmm->tsc[TII][k] * mocc[i-1] * mocc[i]) > imx[cur][k]) {
          imx[cur][k] = sc;
          itb[i][k] = STI;
        }
        /* average over "consensus" sequence */
        for (sym = 0; sym < Alphabet_size; sym++) {
          if (con[i][sym] > 0 && hmm->isc[sym][k] == -INFTY) {
            imx[cur][k] = -INFTY;
            break;
          }
          imx[cur][k] += hmm->isc[sym][k] * con[i][sym];
        }
      }
    }

    /* Now the special states. Order is important here.
     * remember, N, C, and J emissions are zero score by definition.
     */
    /* N state */
    xmx[cur][XMN] = -INFTY;
    xtb[i][XMN]   = STBOGUS;
    if (xmx[prv][XMN] > -INFTY && hmm->xsc[XTN][LOOP] > -INFTY &&
        (sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP] * mocc[i]) > -INFTY) {
      xmx[cur][XMN] = sc;
      xtb[i][XMN] = STN;
    }
    /* E state */
    xmx[cur][XME] = -INFTY;
    xtb[i][XME]   = STBOGUS;
    for (k = 1; k <= hmm->M; k++)
      if (mmx[cur][k] > -INFTY && hmm->esc[k] > -INFTY &&
          (sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME]) {
        xmx[cur][XME] = sc;
        tb->esrc[i] = k;
      }

    /* we don't check J state */
    /* B state; don't connect from J */
    xmx[cur][XMB] = -INFTY;
    xtb[i][XMB]   = STBOGUS;
    if (xmx[cur][XMN] > -INFTY && hmm->xsc[XTN][MOVE] > -INFTY &&
        (sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > xmx[cur][XMB]) {
      xmx[cur][XMB] = sc;
      xtb[i][XMB] = STN;
    }

    /* C state */
    xmx[cur][XMC] = -INFTY;
    xtb[i][XMC]   = STBOGUS;
    if (xmx[prv][XMC] > -INFTY && hmm->xsc[XTC][LOOP] > -INFTY &&
        (sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP] * mocc[i]) > -INFTY) {
      xmx[cur][XMC] = sc;
      xtb[i][XMC] = STC;
    }
    if (xmx[cur][XME] > -INFTY && hmm->xsc[XTE][MOVE] > -INFTY &&
        (sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC]) {
      xmx[cur][XMC] = sc;
      xtb[i][XMC] = STE;
    }
  }
  /* T state (not stored in mx) */
  sc = xmx[msa->alen%2][XMC] + hmm->xsc[XTC][MOVE];

  /* do the traceback */
  tr = ShadowTrace(tb, hmm, msa->alen);
  /* cleanup and return */
  FreePlan7Matrix(mx);
  FreeShadowMatrix(tb);
  for (i = 1; i <= msa->alen; i++)
    free(con[i]);
  free(con);
  free(mocc);

  return tr;
}


struct p7trace_s*
ShadowTrace(
  struct dpshadow_s *tb,
  struct plan7_s *hmm,
  int L
){
  struct p7trace_s *tr;
  int curralloc;    /* current allocated length of trace */
  int tpos;      /* position in trace */
  int i;      /* position in seq (1..N) */
  int k;      /* position in model (1..M) */
  char nxtstate;          /* next state to assign in traceback */

  /* Overallocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + L is minimum trace;
   * add L more as buffer.
   */
  curralloc = L * 2 + 6;
  P7AllocTrace(curralloc, &tr);

  /* Initialization of trace
   * We do it back to front; ReverseTrace() is called later.
   */
  tr->statetype[0] = STT;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tpos     = 1;
  i        = L;      /* current i (seq pos) we're trying to assign   */
  k        = 0;      /* current k (model pos) we're trying to assign */
  nxtstate = STC;    /* assign the C state first, for C->T */

  /* Traceback
   */
  while (nxtstate != STS) {
    switch (nxtstate) {
    case STM:
      tr->statetype[tpos] = STM;
      nxtstate            = tb->mtb[i][k];
      tr->nodeidx[tpos]   = k--;
      tr->pos[tpos]       = i--;
      tpos++;
      break;

    case STI:
      tr->statetype[tpos] = STI;
      nxtstate            = tb->itb[i][k];
      tr->nodeidx[tpos]   = k;
      tr->pos[tpos]       = i--;
      tpos++;
      break;

    case STD:
      tr->statetype[tpos] = STD;
      nxtstate            = tb->dtb[i][k];
      tr->nodeidx[tpos]   = k--;
      tr->pos[tpos]       = 0;
      tpos++;
      break;

    case STN:
      tr->statetype[tpos] = STN;
      nxtstate            = tb->xtb[i][XMN];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]  = (nxtstate == STN) ? i-- : 0; /* N->N; 2nd one emits. */
      tpos++;
      break;

    case STB:
      /* Check for wing unfolding */
      if (Prob2Score(hmm->begin[k+1], hmm->p1) + 1 * INTSCALE <= hmm->bsc[k+1])
        while (k > 0) {
          tr->statetype[tpos] = STD;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = 0;
          tpos++;
          if (tpos == curralloc) {
            /* grow trace if necessary  */
            curralloc += L;
            P7ReallocTrace(tr, curralloc);
          }
        }

      tr->statetype[tpos] = STB;
      nxtstate            = tb->xtb[i][XMB];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = 0;
      tpos++;
      break;

    case STJ:
      tr->statetype[tpos] = STJ;
      nxtstate            = tb->xtb[i][XMJ];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]  = (nxtstate == STJ) ? i-- : 0; /* J->J; 2nd one emits. */
      tpos++;
      break;

    case STE:
      tr->statetype[tpos] = STE;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = 0;
      k                   = tb->esrc[i];
      nxtstate            = STM;
      tpos++;
      /* check for wing unfolding */
      if (Prob2Score(hmm->end[k], 1.) + 1*INTSCALE <=  hmm->esc[k]) {
        int dk;    /* need a tmp k while moving thru delete wing */
        for (dk = hmm->M; dk > k; dk--) {
          tr->statetype[tpos] = STD;
          tr->nodeidx[tpos]   = dk;
          tr->pos[tpos]       = 0;
          tpos++;
          if (tpos == curralloc) {
            /* grow trace if necessary  */
            curralloc += L;
            P7ReallocTrace(tr, curralloc);
          }
        }
      }
      break;

    case STC:
      tr->statetype[tpos] = STC;
      nxtstate            = tb->xtb[i][XMC];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]  = (nxtstate == STC) ? i-- : 0; /* C->C; 2nd one emits. */
      tpos++;
      break;

    default:
      Die("HMMER: Bad state (%s) in ShadowTrace()\n", Statetype(nxtstate));

    } /* end switch over nxtstate */

    if (tpos == curralloc) {
      /* grow trace if necessary  */
      curralloc += L;
      P7ReallocTrace(tr, curralloc);
    }

  } /* end traceback, just before assigning S state */

  tr->statetype[tpos] = STS;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0;
  tr->tlen            = tpos + 1;

  P7ReverseTrace(tr);
  return tr;
}


float
PostprocessSignificantHit(
  struct tophit_s *ghit,
  struct tophit_s *dhit,
  struct p7trace_s *tr,
  struct plan7_s *hmm,
  unsigned char *dsq,
  int L,
  char *seqname,
  char *seqacc,
  char *seqdesc,
  int do_forward,
  float sc_override,
  int do_null2,
  struct threshold_s *thresh,
  int hmmpfam_mode
){
  struct p7trace_s **tarr;      /* array of per-domain traces */
  struct fancyali_s *ali;       /* alignment of a domain      */
  int ntr;      /* number of domain traces from Viterbi */
  int tidx;      /* index for traces (0..ntr-1) */
  int ndom;      /* # of domains accepted in sequence */
  int didx;      /* index for domains (1..ndom) */
  int k1, k2;      /* start, stop coord in model */
  int i1, i2;      /* start, stop in sequence    */
  float   whole_sc;    /* whole sequence score = \sum domain scores */
  float  *score;                /* array of raw scores for each domain */
  int    *usedomain;            /* true if this domain is accepted */
  double  whole_pval;
  double  pvalue;
  double  sortkey;

  /* Special case: rarely, the alignment was totally impossible
   * and tr is NULL.
   */
  if (tr == NULL) return sc_override;

  /* Break the trace into one or more individual domains.
   */
  TraceDecompose(tr, &tarr, &ntr);
  if (ntr == 0) Die("TraceDecompose() screwup"); /* "can't happen" (!) */

  /* Rescore each domain, apply null2 correction if asked.
   * Mark positive-scoring ones (we'll definitely report those),
   * and include their score in the whole sequence score.
   */
  score     = MallocOrDie(sizeof(float) * ntr);
  usedomain = MallocOrDie(sizeof(int)   * ntr);
  ndom      = 0;
  whole_sc  = 0.;
  for (tidx = 0; tidx < ntr; tidx++) {
    score[tidx]  = P7TraceScore(hmm, dsq, tarr[tidx]);
    if (do_null2) score[tidx] -= TraceScoreCorrection(hmm, tarr[tidx], dsq);
    if (score[tidx] > 0.0) {
      usedomain[tidx] = true;
      ndom++;
      whole_sc += score[tidx];
    } else
      usedomain[tidx] = false;
  }

  /* Make sure at least one positive scoring domain is in
   * the trace. If not, invoke "weak single domain" rules:
   * we will always report at least one domain per sequence, even
   * if it has a negative score. (HMMER's Plan7 architecture can report
   * one negative scoring domain but not more.)
   */
  if (ndom == 0) {
    tidx            = FArgMax(score, ntr);
    usedomain[tidx] = true;
    whole_sc        = score[tidx];
    ndom            = 1;
  }

  /* Implement --do_forward: override the trace-dependent sum-of-domain
   * whole score, use the P7Forward() score that the called passed
   * us instead. This is a hack; null2 is trace-dependent and
   * thus undefined for P7Forward() scoring; see commentary in hmmpfam.c.
   */
  if (do_forward) whole_sc = sc_override;

  /* Go through and put all the accepted domains into the hit list.
   */
  whole_pval = PValue(hmm, whole_sc);
  for (tidx = 0, didx = 1; tidx < ntr; tidx++) {
    if (! usedomain[tidx]) continue;

    TraceSimpleBounds(tarr[tidx], &i1, &i2, &k1, &k2);
    pvalue = PValue(hmm, score[tidx]);

    if (pvalue <= thresh->domE && score[tidx] >= thresh->domT) {
      ali     = CreateFancyAli(tarr[tidx], hmm, dsq, seqname);

      if (hmmpfam_mode)
        sortkey = -1.*(double)i1; /* hmmpfam: sort on position in seq    */
      else
        sortkey = score[tidx];    /* hmmsearch: sort on E (monotonic w/ sc) */

      RegisterHit(dhit, sortkey,
                  pvalue,     score[tidx],
                  whole_pval, whole_sc,
                  hmmpfam_mode ? hmm->name : seqname,
                  hmmpfam_mode ? hmm->acc  : seqacc,
                  hmmpfam_mode ? hmm->desc : seqdesc,
                  i1,i2, L,
                  k1,k2, hmm->M,
                  didx,ndom,ali);
    }
    didx++;
  }

  /* Now register the global hit, with the domain-derived score.
   */

  /* sorting:
   * hmmpfam has to worry that score and E-value are not monotonic
   * when multiple HMMs (with different EVD parameters) are potential
   * targets. Therefore in hmmpfam_mode we apply a weird hack
   * to sort primarily on E-value, but on score
   * for really good hits with E=0.0... works because we can
   * assume 100000. > -log(DBL_MIN).
   * hmmsearch simply sorts on score (which for a single HMM, we
   * know is monotonic with E-value).
   */
  if (hmmpfam_mode)
    sortkey = (whole_pval > 0.0) ? -1.*log(whole_pval) : 100000. + whole_sc;
  else
    sortkey = whole_sc;

  /* Note: we've recalculated whole_sc and it may have decreased
   *       after the null2 correction was applied. For Pfam GA, TC,
   *       or NC cutoffs, we have to be sure that everything on the
   *       hitlist is correct (the hmmpfam output routine assumes it,
   *       otherwise it would have to reload each HMM to get its
   *       cutoffs). In all other cases, though, we don't care if
   *       the hit list has a bit too many things on it, because the
   *       output routine in hmmsearch or hmmpfam will check against
   *       the cutoffs. Hence we only need to check against globT
   *       (it may be set by GA, TC, or NC) but not globE.
   *                 - SRE, CSHL genome mtg May 2001
   */
  if (whole_sc >= thresh->globT) {
    RegisterHit(ghit, sortkey,
                whole_pval, whole_sc,
                0., 0.,                    /* no mother seq */
                hmmpfam_mode ? hmm->name : seqname,
                hmmpfam_mode ? hmm->acc  : seqacc,
                hmmpfam_mode ? hmm->desc : seqdesc,
                0,0,0,                    /* seq positions  */
                0,0,0,                    /* HMM positions  */
                0, ndom,            /* # domains info    */
                NULL);                    /* alignment info */
  }

  /* Clean up and return.
   */
  for (tidx = 0; tidx < ntr; tidx++)
    P7FreeTrace(tarr[tidx]);
  free(tarr);
  free(score);
  free(usedomain);
  return whole_sc;
}


float
P7Viterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
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
  int  *mc, *dc;        /* pointers to rows of mmx, dmx */
  int  *mpp, *mpc, *ip;      /* ptrs to mmx[i-1], mmx[i], imx[i-1] */
  int  *bp;         /* ptr into bsc[] */
  int  *dpp;                 /* ptr into dmx[i-1] (previous row) */
  int  *tpmm, *tpmi, *tpmd, *tpim, *tpii, *tpdm, *tpdd; /* ptrs into tsc */
  int   M;

  /* Make sure we have space for a DP matrix with 0..L rows, 0..M-1 columns.
   */
  ResizePlan7Matrix(mx, L, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;                         /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Initializations that help icc vectorize.
   */
  M = hmm->M;

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */

  tpmm  = hmm->tsc[TMM];
  tpim  = hmm->tsc[TIM];
  tpdm  = hmm->tsc[TDM];
  tpmd  = hmm->tsc[TMD];
  tpdd  = hmm->tsc[TDD];
  tpmi  = hmm->tsc[TMI];
  tpii  = hmm->tsc[TII];
  bp    = hmm->bsc;
  for (i = 1; i <= L; i++) {
    int  *ic; /* pointers to rows of imx */
    int  *ms; /* pointers to msc[i] */
    int  *is; /* pointers to isc[i] */
    int  *ep; /* ptr into esc[] */
    int   xmb;/* value of xmx[i-1][XMB] */
    int   xme;/* max for xmx[i][XME] */
    mc    = mmx[i];
    dc    = dmx[i];
    ic    = imx[i];
    mpp   = mmx[i-1];
    dpp   = dmx[i-1];
    ip    = imx[i-1];
    xmb   = xmx[i-1][XMB];
    ms    = hmm->msc[dsq[i]];
    is    = hmm->isc[dsq[i]];
    mc[0] = -INFTY;
    dc[0] = -INFTY;
    ic[0] = -INFTY;

    for (k = 1; k <= M; k++) {
      mc[k] = mpp[k-1]   + tpmm[k-1];
      if ((sc = ip[k-1]  + tpim[k-1]) > mc[k])  mc[k] = sc;
      if ((sc = dpp[k-1] + tpdm[k-1]) > mc[k])  mc[k] = sc;
      if ((sc = xmb  + bp[k])         > mc[k])  mc[k] = sc;
      mc[k] += ms[k];
      if (mc[k] < -INFTY) mc[k] = -INFTY;

      dc[k] = dc[k-1] + tpdd[k-1];
      if ((sc = mc[k-1] + tpmd[k-1]) > dc[k]) dc[k] = sc;
      if (dc[k] < -INFTY) dc[k] = -INFTY;

      if (k < M) {
        ic[k] = mpp[k] + tpmi[k];
        if ((sc = ip[k] + tpii[k]) > ic[k]) ic[k] = sc;
        ic[k] += is[k];
        if (ic[k] < -INFTY) ic[k] = -INFTY;
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
    /* N state */
    xmx[i][XMN] = -INFTY;
    if ((sc = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[i][XMN] = sc;

    /* E state */
    xme = -INFTY;
    mpc = mmx[i];
    ep  = hmm->esc;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mpc[k] + ep[k]) > xme) xme = sc;
    xmx[i][XME] = xme;
    /* J state */
    xmx[i][XMJ] = -INFTY;
    if ((sc = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
      xmx[i][XMJ] = sc;
    if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
      xmx[i][XMJ] = sc;

    /* B state */
    xmx[i][XMB] = -INFTY;
    if ((sc = xmx[i][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[i][XMB] = sc;
    if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
      xmx[i][XMB] = sc;

    /* C state */
    xmx[i][XMC] = -INFTY;
    if ((sc = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[i][XMC] = sc;
    if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
      xmx[i][XMC] = sc;
  }
  /* T state (not stored) */
  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_tr != NULL) {
    P7ViterbiTrace(hmm, dsq, L, mx, &tr);
    *ret_tr = tr;
  }

  return Scorify(sc);    /* the total Viterbi score. */
}


void
AllocPlan7Body(
  struct plan7_s *hmm,
  int M
){
  int k, x;

  hmm->M = M;

  hmm->rf     = MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = MallocOrDie ((M+2) * sizeof(char));
  hmm->ca     = MallocOrDie ((M+2) * sizeof(char));
  hmm->map    = MallocOrDie ((M+1) * sizeof(int));

  hmm->t      = MallocOrDie (M     *           sizeof(float *));
  hmm->tsc    = MallocOrDie (7     *           sizeof(int *));
  hmm->mat    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = MallocOrDie (M     *           sizeof(float *));
  hmm->msc    = MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc    = MallocOrDie (MAXCODE   *       sizeof(int *));

  hmm->t[0]   = MallocOrDie ((7*M)     *       sizeof(float));
  /* Allocate extra memory so tsc[TMM,TIM,TDM,TMD,TDD] start on the
   * 16-byte cache boundary, and tsc[TMI,TII] start
   * 12 bytes offset from the boundary.
   */
  hmm->tsc_mem = MallocOrDie (((7*(M+16)))  *   sizeof(int));
  hmm->mat[0] = MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((MAXABET*M) *     sizeof(float));
  /* Allocate extra mem. to make sure all members of msc,isc start
   * on 12-byte offsets from cache boundary.
   */
  hmm->msc_mem = MallocOrDie ((MAXCODE*(M+1+16)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXCODE*(M+16)) *   sizeof(int));

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * MAXABET;
      hmm->t[k]   = hmm->t[0]   + k * 7;
    }
  }

  /* align tsc pointers */
  hmm->tsc[TMM] = (int *) (((((size_t) hmm->tsc_mem) + 15) & (~0xf)));
  hmm->tsc[TMI] = (int *) (((((size_t) hmm->tsc_mem) + (M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TMD] = (int *) (((((size_t) hmm->tsc_mem) + 2*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TIM] = (int *) (((((size_t) hmm->tsc_mem) + 3*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TII] = (int *) (((((size_t) hmm->tsc_mem) + 4*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TDM] = (int *) (((((size_t) hmm->tsc_mem) + 5*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TDD] = (int *) (((((size_t) hmm->tsc_mem) + 6*(M+12)*sizeof(int) + 15) & (~0xf)));

  for (x = 0; x < MAXCODE; x++) {
    hmm->msc[x] = (int *) (((((size_t)hmm->msc_mem) + x*(M+1+12)*sizeof(int) + 15) & (~0xf)) + 12);
    hmm->isc[x] = (int *) (((((size_t)hmm->isc_mem) + x*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  }
  /* tsc[0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < 7; x++)
    hmm->tsc[x][0] = -INFTY;

  hmm->begin  = MallocOrDie  ((M+1) * sizeof(float));
  hmm->bsc_mem= MallocOrDie  ((M+1+12) * sizeof(int));
  hmm->end    = MallocOrDie  ((M+1) * sizeof(float));
  hmm->esc_mem= MallocOrDie  ((M+1+12) * sizeof(int));

  hmm->bsc = (int *) (((((size_t) hmm->bsc_mem) + 15) & (~0xf)) + 12);
  hmm->esc = (int *) (((((size_t) hmm->esc_mem) + 15) & (~0xf)) + 12);

  return;
}


struct dpmatrix_s *
CreatePlan7Matrix(
  int N,
  int M,
  int padN,
  int padM
){
  struct dpmatrix_s *mx;

  mx         = (struct dpmatrix_s *) MallocOrDie (sizeof(struct dpmatrix_s));
  mx->xmx    = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->mmx    = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->imx    = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->dmx    = (int **) MallocOrDie (sizeof(int *) * (N+1));

  /* For the memory accessed by the altivec routines, we want to have
   * accesses aligned to 16-byte boundaries as far as possible.
   * To accomplish this, we align extra memory and then set the first
   * pointer on each row to point to 4 bytes before a boundary.
   * This means element 1, which is the first one we work on, will be
   * on a 16-byte boundary. We still make sure we 'own' the three bytes
   * before, though, so we can load the vector with element 0 cache-aligned too.
   * The real pointers to memory are kept in xmx_mem,mmx_mem,imx_mem,dmx_mem.
   */
  mx->xmx_mem = (void *) MallocOrDie (sizeof(int) * (N+1)*(5 + 16));
  mx->mmx_mem = (void *) MallocOrDie (sizeof(int) * (N+1)*(M+2+16));
  mx->imx_mem = (void *) MallocOrDie (sizeof(int) * (N+1)*(M+2+16));
  mx->dmx_mem = (void *) MallocOrDie (sizeof(int) * (N+1)*(M+2+16));

  mx->xmx[0] = (int *) (((((size_t) mx->xmx_mem) + 15) & (~0xf)) + 12);
  mx->mmx[0] = (int *) (((((size_t) mx->mmx_mem) + 15) & (~0xf)) + 12);
  mx->imx[0] = (int *) (((((size_t) mx->imx_mem) + 15) & (~0xf)) + 12);
  mx->dmx[0] = (int *) (((((size_t) mx->dmx_mem) + 15) & (~0xf)) + 12);

  /* And make sure the beginning of each row is aligned the same way */
  for (int i = 1; i <= N; i++) {
    mx->xmx[i] = mx->xmx[0] + i*(5+11) ; /* add 11 bytes per row, making it divisible by 4 */
    int n = 12 - (M+2)%4;
    mx->mmx[i] = mx->mmx[0] + i*(M+2+n);
    mx->imx[i] = mx->imx[0] + i*(M+2+n);
    mx->dmx[i] = mx->dmx[0] + i*(M+2+n);
  }

  mx->workspace_mem = (int *) MallocOrDie (sizeof(int) * (M+4) * 12 + 16);
  mx->workspace     = (int *) ((((size_t) mx->workspace_mem) + 15) & (~0xf));

  mx->maxN = N;
  mx->maxM = M;
  mx->padN = padN;
  mx->padM = padM;

  return mx;
}


void
ResizePlan7Matrix(
  struct dpmatrix_s *mx,
  int N,
  int M,
  int ***xmx,
  int ***mmx,
  int ***imx,
  int ***dmx
){

  if (N <= mx->maxN && M <= mx->maxM) {
    if (xmx != NULL) *xmx = mx->xmx;
    if (mmx != NULL) *mmx = mx->mmx;
    if (imx != NULL) *imx = mx->imx;
    if (dmx != NULL) *dmx = mx->dmx;
    return;
  }

  if (N > mx->maxN) {
    N          += mx->padN;
    mx->maxN    = N;
    mx->xmx     = (int **) ReallocOrDie (mx->xmx, sizeof(int *) * (mx->maxN+1));
    mx->mmx     = (int **) ReallocOrDie (mx->mmx, sizeof(int *) * (mx->maxN+1));
    mx->imx     = (int **) ReallocOrDie (mx->imx, sizeof(int *) * (mx->maxN+1));
    mx->dmx     = (int **) ReallocOrDie (mx->dmx, sizeof(int *) * (mx->maxN+1));
  }

  if (M > mx->maxM) {
    M += mx->padM;
    mx->maxM = M;
  }

  mx->xmx_mem = ReallocOrDie (mx->xmx_mem, sizeof(int) * (mx->maxN+1)*(5 + 16));
  mx->mmx_mem = ReallocOrDie (mx->mmx_mem, sizeof(int) * (mx->maxN+1)*(mx->maxM+2+16));
  mx->imx_mem = ReallocOrDie (mx->imx_mem, sizeof(int) * (mx->maxN+1)*(mx->maxM+2+16));
  mx->dmx_mem = ReallocOrDie (mx->dmx_mem, sizeof(int) * (mx->maxN+1)*(mx->maxM+2+16));

  mx->xmx[0] = (int *) (((((size_t) mx->xmx_mem) + 15) & (~0xf)) + 12);
  mx->mmx[0] = (int *) (((((size_t) mx->mmx_mem) + 15) & (~0xf)) + 12);
  mx->imx[0] = (int *) (((((size_t) mx->imx_mem) + 15) & (~0xf)) + 12);
  mx->dmx[0] = (int *) (((((size_t) mx->dmx_mem) + 15) & (~0xf)) + 12);

  /* And make sure the beginning of each row is aligned the same way */
  for (int i = 1; i <= mx->maxN; i++) {
    mx->xmx[i] = mx->xmx[0] + i*(5+11) ; /* add 11 bytes per row, making it divisible by 4 */
    int n = 12 - (mx->maxM+2)%4;
    mx->mmx[i] = mx->mmx[0] + i*(mx->maxM+2+n);
    mx->imx[i] = mx->imx[0] + i*(mx->maxM+2+n);
    mx->dmx[i] = mx->dmx[0] + i*(mx->maxM+2+n);
  }

  if (xmx != NULL) *xmx = mx->xmx;
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
}

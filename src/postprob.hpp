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

/* postprob.h
 *
 * Functions for working with posterior probabilities,
 * including unfussed "backwards" and "optimal accuracy"
 * implementations.
 */
/* postprob.c
 * Author: Ian Holmes (ihh@sanger.ac.uk, Jun 5 1998)
 * Derived from core_algorithms.c (SRE, Nov 11 1996)
 *
 *****************************************************************
 * IHH's notes:
 *
 * Functions for working with posterior probabilities,
 * including unfussed "backwards" and "optimal accuracy"
 * implementations.
 *****************************************************************
 * SRE's notes:
 *
 * Simple API example:
 *     struct p7trace_s  *tr;
 *     struct dpmatrix_s *fwd;
 *     struct dpmatrix_s *bck;
 *     struct dpmatrix_s *posterior;
 *     char *postcode;
 *
 *     (get a traceback from somewhere: P7Viterbi() or a modelmaker)
 *     (get an HMM from somewhere: read file or construct it)
 *     P7Forward (dsq, len, hmm, &fwd);
 *     P7Backward(dsq, len, hmm, &bck);
 *     posterior = bck;              -- can alloc posterior, but also can re-use bck --
 *     P7EmitterPosterior(len, hmm, fwd, bck, posterior);
 *     postcode = PostalCode(len, posterior, tr);
 *
 *     MSAAppendGR(msa, "POST", seqidx, postcode);  -- or a similar annotation call --
 *
 *     free(postcode);
 *     FreePlan7Matrix(fwd);
 *     FreePlan7Matrix(bck);
 *
 * P7OptimalAccuracy() - the Durbin/Holmes optimal accuracy
 *                       alignment algorithm. Takes a sequence
 *                       and an HMM, returns an alignment as
 *                       a trace structure.
 *
 * P7Backward()        - The Backward() algorithm, counterpart
 *                       of P7Forward() in core_algorithms.c.
 *
 * P7EmitterPosterior()- The heart of postprob.c: given a Forward
 *                       and a Backward matrix, calculate a new matrix
 *                       that contains the posterior probabilities
 *                       for each symbol i being emitted by
 *                       state k (so, \sum_k p(k | x_i) = 1.0).
 *
 * P7FillOptimalAccuracy() - The core DP algorithm called by
 *                       P7OptimalAccuracy().
 *
 * P7OptimalAccuracyTrace() - the traceback algorithm called by
 *                       P7FillOptimalAccuracy().
 *
 * PostalCode()        - Create a character string for annotating
 *                       an alignment.
 *
 * No small memory variants of these algorithms are available
 * right now.
 */

#pragma once

#include "structs.hpp"

/* Extra algorithms to work with posterior probabilities.
 */

float
P7OptimalAccuracy(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr);


/* Function: P7Backward()
 *
 * Purpose:  The Backward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
P7Backward(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s **ret_mx);


/* Function: P7EmitterPosterior()
 *
 * Purpose:  Combines Forward and Backward matrices into a posterior
 *           probability matrix.
 *           The entries in row i of this matrix are the logs of the
 *           posterior probabilities of each state emitting symbol i of
 *           the sequence, i.e. all entries for non-emitting states are -INFTY.
 *           The caller must allocate space for the matrix, although the
 *           backward matrix can be used instead (overwriting it will not
 *           compromise the algorithm).
 *
 * Args:     L        - length of sequence
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           mx       - pre-allocated dynamic programming matrix
 *
 * Return:   void
 */
void 
P7EmitterPosterior(
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s *forward,
  struct dpmatrix_s *backward,
  struct dpmatrix_s *mx);


/* Function: P7FillOptimalAccuracy()
 *
 * Purpose:  The core of the optimal accuracy dynamic programming algorithm.
 *           Identical to Viterbi() except that scores are given by a
 *           posterior matrix (that the caller must pre-calculate).
 *           Also, the caller must pre-allocate the optimal accuracy matrix
 *           (this allows the forward matrix to be re-used).
 *           P7OptimalAccuracy() does all this for you and cleans up.
 *
 *
 * Args:     L         - length of sequence
 *           M         - length of model
 *           posterior - pre-calculated emitter posterior matrix
 *           mx        - pre-allocated dynamic programming matrix
 *           ret_tr    - RETURN: traceback; pass NULL if it's not wanted
 *
 * Return:   log ( sum_{residues} P(label|M,D) ), as a bit score
 *           (i.e. log of expected accuracy)
 */
float
P7FillOptimalAccuracy(
  int L,
  int M,
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);


/* Function: P7OptimalAccuracyTrace()
 *
 * Purpose:  Traceback of an optimal accuracy matrix: i.e. retrieval
 *           of optimum alignment.
 *
 * Args:     L         - length of sequence
 *           M         - length of HMM
 *           posterior - the posterior matrix
 *           mx        - the matrix to trace back in, (L+1) x M
 *           ret_tr    - RETURN: traceback.
 *
 * Return:   (void)
 *           ret_tr is allocated here. Free using P7FreeTrace().
 */
void 
P7OptimalAccuracyTrace(
  int L,
  int M,
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);

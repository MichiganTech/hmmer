/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* infocontent.c
 *
 * for evolving models to a specified information content
 */


#include "structs.hpp"


/* Function: AdjustAveInfoContent()
 *
 * Purpose:  evolve match emissions to specified average info content
 *           over the length of the sequence
 *
 * Args:     hmm   - HMM with match emissions to evolve
 *       desired   - desired ave information content for match emissions
 *           matrixfile - file containing rate matrix(ces)
 *        and corresponding amino acid freqs
 *
 * Return:   (void)
 */
void
AdjustAveInfoContent (
  struct plan7_s *hmm,
  float desired,
  char *matrixfile);


/* Function: CalculateBackgroundEntropy()
 *
 * Purpose: determine the background entropy of an alignment based on some
 *          expected amino acid frequency
 *
 * Args:  none
 *
 * Return: (entropy)
 */
float
CalculateBackgroundEntropy();


/* Function: CalculateEmitsEntropy()
 *
 * Purpose: determine entropy of match emits
 *
 * Args:  emits - 20 x L array holding match emits for each node of model
 *    L - length of the hmm
 *
 *
 * Return: (entropy)
 */
float
CalculateEmitsEntropy(
  double *emits,
  int L);


/* Function: NormalizeEmits()
 *
 * Purpose: normalize emission probabilities of an hmm
 *
 * Args:  emits - 20 x L array holding match emits for each node of model
 *    L - length of the hmm
 *
 * Return: (void)
 */
void
NormalizeEmits(
  double *emits,
  int L);


/* Function: EvolveEmits()
 *
 * Purpose: evolve emission probabilities of an hmm
 *
 * Args:  emits - 20 x L array holding match emits for each node of model
 *    L - length of the hmm
 *        P - conditional matrices representing the evolutionary model(s)
 *    nmodels - number of evolutionary models
 *
 * Return: (void)
 */
void
EvolveEmits(
  double *emits,
  double *P,
  int L);

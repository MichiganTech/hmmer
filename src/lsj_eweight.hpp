/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* lsjfuncs.c
 *
 * entropy targeting:
 * Code for setting effective sequence number (in hmmbuild) by
 * achieving a certain target entropy loss, relative to background
 * null distribution.
 */

#include "structs.hpp"
#include "lsjfuncs.hpp"

/* Function: Eweight()
 *
 * Purpose:  Main entropy-based weighting function.
 *
 * Args:
 *           **mat       - Current model match state counts.
 *           **pri       - Model priors.
 *       numb_seqs       - Number of sequences in alignment.
 *       targetent       - Target mean match state entropy.
 *
 * Return: eff_no        - New effective sequence number.
 */
float
Eweight(
  struct plan7_s *hmm, 
  struct p7prior_s *pri,
  float numb_seqs,
  float targetent);
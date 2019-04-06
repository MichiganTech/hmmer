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

#include <ctype.h>
#include <stdlib.h>

#include "structs.hpp"


/* Function: EmitSequence()
 *
 * Purpose:  Given a model, sample a sequence and/or traceback.
 *
 * Args:     hmm     - the model
 *           ret_dsq - RETURN: generated digitized sequence (pass NULL if unwanted)
 *           ret_L   - RETURN: length of generated sequence
 *           ret_tr  - RETURN: generated trace (pass NULL if unwanted)
 *
 * Returns:  void
 */
void
EmitSequence(
  struct plan7_s *hmm,
  unsigned char **ret_dsq,
  int *ret_L,
  struct p7trace_s **ret_tr);



/* Function: EmitConsensusSequence()
 *
 * Purpose:  Generate a "consensus sequence". For the purposes
 *           of a profile HMM, this is defined as:
 *              - for each node:
 *                 - if StateOccupancy() says that M is used
 *                     with probability >= 0.5, this M is "consensus".
 *                     Then, choose maximally likely residue.
 *                     if P>0.5 (protein) or P>0.9 (DNA), make
 *                     it upper case; else make it lower case.
 *                 - if StateOccupancy() says that I
 *                     is used with P >= 0.5, this I is "consensus";
 *                     use it 1/(1-TII) times (its expectation value).
 *                     Generate an "x" from each I.
 *
 *           The function expects that the model is config'ed
 *           by Plan7NakedConfig(): that is, for a single global pass
 *           with no N,C,J involvement.
 *
 *
 * Args:     hmm     - the model
 *           ret_seq - RETURN: consensus sequence (pass NULL if unwanted)
 *           ret_dsq - RETURN: digitized consensus sequence (pass NULL if unwanted)
 *           ret_L   - RETURN: length of generated sequence
 *           ret_tr  - RETURN: generated trace (pass NULL if unwanted)
 *
 * Returns:  void
 */
void
EmitConsensusSequence(
  struct plan7_s *hmm,
  char **ret_seq,
  unsigned char **ret_dsq,
  int *ret_L,
  struct p7trace_s **ret_tr);


/* Function: StateOccupancy()
 *
 * Purpose:  Calculate the expected state occupancy for
 *           a given HMM in generated traces.
 *
 *           Note that expected prob of getting into
 *           any special state in a trace is trivial:
 *              S,N,B,E,C,T = 1.0
 *              J = E->J transition prob
 *
 * Args:     hmm    - the model
 *           ret_mp - RETURN: [1..M] prob's of occupying M
 *           ret_ip - RETURN: [1..M-1] prob's of occupying I
 *           ret_dp - RETURN: [1..M] prob's of occupying D
 *
 * Returns:  void
 *           mp, ip, dp are malloc'ed here. Caller must free().
 */
void
StateOccupancy(
  struct plan7_s *hmm,
  float **ret_mp,
  float **ret_ip,
  float **ret_dp);

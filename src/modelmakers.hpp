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

#include "msa.hpp"
#include "structs.hpp"


/* Function: P7Handmodelmaker()
 *
 * Purpose:  Manual model construction:
 *           Construct an HMM from an alignment, where the #=RF line
 *           of a HMMER alignment file is given to indicate
 *           the columns assigned to matches vs. inserts.
 *
 *           NOTE: Handmodelmaker() will slightly revise the alignment
 *           if necessary, if the assignment of columns implies
 *           DI and ID transitions.
 *
 *           Returns both the HMM in counts form (ready for applying
 *           Dirichlet priors as the next step), and fake tracebacks
 *           for each aligned sequence.
 *
 * Args:     msa  - multiple sequence alignment
 *           dsq  - digitized unaligned aseq's
 *           isfrag  - [0..nseq-1] flags for candidate seq frags
 *           ret_hmm - RETURN: counts-form HMM
 *           ret_tr  - RETURN: array of tracebacks for aseq's
 *
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here; FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm).
 */
void
P7Handmodelmaker(
  MSA *msa,
  unsigned char **dsq,
  char *isfrag,
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr
);


/* Function: P7Fastmodelmaker()
 *
 * Purpose:  Heuristic model construction.
 *
 *           Construct an HMM from an alignment using a heuristic,
 *           based on the fractional occupancy of each columns w/
 *           residues vs gaps. Roughly, any column w/ a fractional
 *           occupancy of $\geq$ <symfrac> is assigned as a MATCH column;
 *           for instance, if thresh = 0.5, columns w/ $\geq$ 50\%
 *           residues are assigned to match... roughly speaking.
 *
 *           "Roughly speaking" because sequences are weighted; we
 *           guard against fragments counting against us by not
 *           counting frags at all unless we have to; and we guard
 *           against highly gappy alignments by reweighting the threshold
 *           according to the most occupied column in the alignment.
 *
 *           The detailed calculation is:
 *
 *           Count the weighted fraction r_i of symbols in column i
 *           (weighted by relative sequence weights): 0<=r_i<=1, with
 *           r_i = 1.0 in fully occupied columns. Only sequences that
 *           are not flagged as candidate fragments are counted toward
 *           the r_i calculation. (If *all* sequences are candidate
 *           fragments, then revert to counting them, rather than
 *           giving up w/ undefined r_i's.)
 *
 *           Determine the fraction of symbols in the most occupied
 *           column: R = max_i r_i. Normally R=1.0, but if we're given
 *           a gappy alignment, R may be less than that.
 *
 *           Then given threshold t:
 *
 *           if r_i >= tR, column is assigned as a MATCH column;
 *           else it's an INSERT column.
 *
 *           NOTE: p7_Fastmodelmaker() will slightly revise the
 *           alignment if the assignment of columns implies
 *           DI and ID transitions.
 *
 *           Returns the HMM in counts form (ready for applying Dirichlet
 *           priors as the next step). Also returns fake traceback
 *           for each training sequence.
 *
 * Args:     msa       - multiple sequence alignment
 *           abc       - symbol alphabet to use
 *           dsq       - digitized unaligned aseq's
 *           isfrag    - [0..nseq-1] T/F flags for candidate seq frags
 *           symfrac   - threshold for residue occupancy
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here; FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm).
 */
void
P7Fastmodelmaker(
  MSA *msa,
  unsigned char **dsq,
  char *isfrag,
  float symfrac,
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr);


/* Function: matassign2hmm()
 *
 * Purpose:  Given an assignment of alignment columns to match vs.
 *           insert, finish the final part of the model construction
 *           calculation that is constant between model construction
 *           algorithms.
 *
 * Args:     msa       - multiple sequence alignment
 *           dsq       - digitized unaligned aseq's
 *           matassign - 1..alen bit flags for column assignments
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here for the calling
 *           modelmaker function.
 */
void
matassign2hmm(
  MSA *msa,
  unsigned char **dsq,
  char *isfrag,
  int *matassign,
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr);


/* Function: fake_tracebacks()
 *
 * Purpose:  From a consensus assignment of columns to MAT/INS, construct fake
 *           tracebacks for each individual sequence.
 *
 * Note:     Fragment tolerant by default. Internal entries are
 *           B->M_x, instead of B->D1->D2->...->M_x; analogously
 *           for internal exits.
 *
 * Args:     aseqs     - alignment [0..nseq-1][0..alen-1]
 *           nseq      - number of seqs in alignment
 *           alen      - length of alignment in columns
 *           isfrag    - T/F flags for candidate fragments
 *           matassign - assignment of column; [1..alen] (off one from aseqs)
 *           ret_tr    - RETURN: array of tracebacks
 *
 * Return:   (void)
 *           ret_tr is alloc'ed here. Caller must free.
 */
void
fake_tracebacks(
  char **aseq,
  int nseq,
  int alen,
  char *isfrag,
  int *matassign,
  struct p7trace_s ***ret_tr);


/* Function: trace_doctor()
 *
 * Purpose:  Plan 7 disallows D->I and I->D "chatter" transitions.
 *           However, these transitions may be implied by many
 *           alignments for hand- or heuristic- built HMMs.
 *           trace_doctor() collapses I->D or D->I into a
 *           single M position in the trace.
 *           Similarly, B->I and I->E transitions may be implied
 *           by an alignment.
 *
 *           trace_doctor does not examine any scores when it does
 *           this. In ambiguous situations (D->I->D) the symbol
 *           will be pulled arbitrarily to the left, regardless
 *           of whether that's the best column to put it in or not.
 *
 * Args:     tr      - trace to doctor
 *           M       - length of model that traces are for
 *           ret_ndi - number of DI transitions doctored
 *           ret_nid - number of ID transitions doctored
 *
 * Return:   (void)
 *           tr is modified
 */
void
trace_doctor(
  struct p7trace_s *tr,
  int mlen,
  int *ret_ndi,
  int *ret_nid);


/* Function: annotate_model()
 *
 * Purpose:  Add rf, cs optional annotation to a new model.
 *
 * Args:     hmm       - new model
 *           matassign - which alignment columns are MAT; [1..alen]
 *           msa       - alignment, including annotation to transfer
 *
 * Return:   (void)
 */
void
annotate_model(
  struct plan7_s *hmm,
  int *matassign,
  MSA *msa);

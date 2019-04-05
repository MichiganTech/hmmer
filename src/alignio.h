/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* alignio.c
 * 
 * Input/output of sequence alignments.
 */

#include "squid.h"

/* Function: AllocAlignment()
 * 
 * Purpose:  Allocate space for an alignment, given the number
 *           of sequences and the alignment length in columns.
 *           
 * Args:     nseq     - number of sequences
 *           alen     - width of alignment
 *           ret_aseq - RETURN: alignment itself
 *           ainfo    - RETURN: other info associated with alignment
 *           
 * Return:   (void)
 *           aseq, ainfo free'd by caller: FreeAlignment(aseq, &ainfo).
 *           note that ainfo itself is alloc'ed in caller, usually
 *           just by a "AINFO ainfo" definition.
 */
void
AllocAlignment(
  int nseq, 
  int alen, 
  char ***ret_aseq, 
  AINFO *ainfo);
 

/* Function: InitAinfo()
 *
 * Purpose:  Initialize the fields in ainfo structure to
 *           default (null) values. Does nothing with 
 *           fields that are dependent on nseq or alen.
 *
 * Args:     ainfo  - optional info structure for an alignment
 *
 * Returns:  (void). ainfo is modified.
 */
void
InitAinfo(
  AINFO *ainfo);


/* Function: FreeAlignment()
 * 
 * Purpose:  Free the space allocated to alignment, names, and optional
 *           information. 
 *           
 * Args:     aseqs - sequence alignment
 *           ainfo - associated alignment data.
 */                  
void
FreeAlignment(
  char **aseqs, 
  AINFO *ainfo);


/* Function: SAMizeAlignment()
 *
 * Purpose:  Make a "best effort" attempt to convert an alignment
 *           to SAM gap format: - in delete col, . in insert col.
 *           Only works if alignment adheres to SAM's upper/lower
 *           case convention, which is true for instance of old
 *           HMMER alignments.
 *
 * Args:     aseq  - alignment to convert
 *           nseq  - number of seqs in alignment
 *           alen  - length of alignment
 *
 * Returns:  (void)
 */
void
SAMizeAlignment(
  char **aseq, 
  int nseq, 
  int alen);


/* Function: SAMizeAlignmentByGapFrac()
 *
 * Purpose:  Convert an alignment to SAM's gap and case
 *           conventions, using gap fraction in a column
 *           to choose match versus insert columns. In match columns,
 *           residues are upper case and gaps are '-'.
 *           In insert columns, residues are lower case and
 *           gaps are '.'
 *
 * Args:     aseq   - aligned sequences
 *           nseq   - number of sequences
 *           alen   - length of alignment
 *           maxgap - if more gaps than this fraction, column is insert.
 *
 * Returns:  (void) Characters in aseq may be altered.
 */
void
SAMizeAlignmentByGapFrac(
  char **aseq, 
  int nseq, 
  int alen, 
  float maxgap);


/* Function: MakeAlignedString()
 * 
 * Purpose:  Given a raw string of some type (secondary structure, say),
 *           align it to a given aseq by putting gaps wherever the
 *           aseq has gaps. 
 *           
 * Args:     aseq:  template for alignment
 *           alen:  length of aseq
 *           ss:    raw string to align to aseq
 *           ret_s: RETURN: aligned ss
 *           
 * Return:   1 on success, 0 on failure (and squid_errno is set.)
 *           ret_ss is malloc'ed here and must be free'd by caller.
 */
int
MakeAlignedString(
  char *aseq, 
  int alen, 
  char *ss, 
  char **ret_s);


/* Function: MakeDealignedString()
 * 
 * Purpose:  Given an aligned string of some type (either sequence or 
 *           secondary structure, for instance), dealign it relative
 *           to a given aseq. Return a ptr to the new string.
 *           
 * Args:     aseq  : template alignment 
 *           alen  : length of aseq
 *           ss:   : string to make dealigned copy of; same length as aseq
 *           ret_s : RETURN: dealigned copy of ss
 *           
 * Return:   1 on success, 0 on failure (and squid_errno is set)
 *           ret_s is alloc'ed here and must be freed by caller
 */
int
MakeDealignedString(
  char *aseq, 
  int alen, 
  char *ss, 
  char **ret_s);


/* Function: DealignedLength()
 *
 * Purpose:  Count the number of non-gap symbols in seq.
 *           (i.e. find the length of the unaligned sequence)
 * 
 * Args:     aseq - aligned sequence to count symbols in, \0 terminated
 * 
 * Return:   raw length of seq.
 */
int
DealignedLength(
  char *aseq);


/* Function: WritePairwiseAlignment()
 * 
 * Purpose:  Write a nice formatted pairwise alignment out,
 *           with a BLAST-style middle line showing identities
 *           as themselves (single letter) and conservative
 *           changes as '+'.
 *           
 * Args:     ofp          - open fp to write to (stdout, perhaps)
 *           aseq1, aseq2 - alignments to write (not necessarily 
 *                          flushed right with gaps)
 *           name1, name2 - names of sequences
 *           spos1, spos2 - starting position in each (raw) sequence
 *           pam          - PAM matrix; positive values define
 *                          conservative changes
 *           indent       - how many extra spaces to print on left
 *           
 * Return:  1 on success, 0 on failure
 */
int
WritePairwiseAlignment(
  FILE *ofp,
  char *aseq1, 
  char *name1, 
  int spos1,
  char *aseq2, 
  char *name2, 
  int spos2,
  int **pam, 
  int indent);


/* Function: MingapAlignment()
 * 
 * Purpose:  Remove all-gap columns from a multiple sequence alignment
 *           and its associated data. The alignment is assumed to be
 *           flushed (all aseqs the same length).
 */
int
MingapAlignment(
  char **aseqs, 
  AINFO *ainfo);


/* Function: RandomAlignment()
 * 
 * Purpose:  Create a random alignment from raw sequences.
 * 
 *           Ideally, we would like to sample an alignment from the
 *           space of possible alignments according to its probability,
 *           given a prior probability distribution for alignments.
 *           I don't see how to describe such a distribution, let alone
 *           sample it.
 *           
 *           This is a rough approximation that tries to capture some
 *           desired properties. We assume the alignment is generated
 *           by a simple HMM composed of match and insert states.
 *           Given parameters (pop, pex) for the probability of opening
 *           and extending an insertion, we can find the expected number
 *           of match states, M, in the underlying model for each sequence.
 *           We use an average M taken over all the sequences (this is
 *           an approximation. The expectation of M given all the sequence
 *           lengths is a nasty-looking summation.)
 *           
 *           M = len / ( 1 + pop ( 1 + 1/ (1-pex) ) ) 
 *           
 *           Then, we assign positions in each raw sequence onto the M match
 *           states and M+1 insert states of this "HMM", by rolling random
 *           numbers and inserting the (rlen-M) inserted positions randomly
 *           into the insert slots, taking into account the relative probability
 *           of open vs. extend.
 *           
 *           The resulting alignment has two desired properties: insertions
 *           tend to follow the HMM-like exponential distribution, and
 *           the "sparseness" of the alignment is controllable through
 *           pop and pex.
 *           
 * Args:     rseqs     - raw sequences to "align", 0..nseq-1
 *           sqinfo    - array of 0..nseq-1 info structures for the sequences
 *           nseq      - number of sequences   
 *           pop       - probability to open insertion (0<pop<1)
 *           pex       - probability to extend insertion (0<pex<1)
 *           ret_aseqs - RETURN: alignment (flushed)
 *           ainfo     - fill in: alignment info
 * 
 * Return:   1 on success, 0 on failure. Sets squid_errno to indicate cause
 *           of failure.
 */                      
int
RandomAlignment(
  char **rseqs, 
  SQINFO *sqinfo, 
  int nseq, 
  float pop, 
  float pex,
  char ***ret_aseqs, 
  AINFO *ainfo);


/* Function: AlignmentHomogenousGapsym()
 *
 * Purpose:  Sometimes we've got to convert alignments to 
 *           a lowest common denominator, and we need
 *           a single specific gap character -- for example,
 *           PSI-BLAST blastpgp -B takes a very simplistic
 *           alignment input format which appears to only
 *           allow '-' as a gap symbol.
 *           
 *           Anything matching the isgap() macro is
 *           converted.
 *
 * Args:     aseq   - aligned character strings, [0..nseq-1][0..alen-1]
 *           nseq   - number of aligned strings
 *           alen   - length of alignment
 *           gapsym - character to use for gaps.         
 *
 * Returns:  void ("never fails")
 */
void
AlignmentHomogenousGapsym(
  char **aseq, 
  int nseq, 
  int alen, 
  char gapsym);

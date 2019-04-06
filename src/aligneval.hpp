#pragma once
/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* aligneval.c
 *
 * Comparison of multiple alignments. Three functions are
 * provided, using subtly different scoring schemes:
 *    CompareMultAlignments()    - basic scoring scheme
 *    CompareRefMultAlignments() - only certain "canonical" columns
 *                                 are scored
 *                                
 * The similarity measure is a fractional alignment identity averaged
 * over all sequence pairs. The score for all pairs is:
 *      (identically aligned symbols) / (total aligned columns in
 *      known alignment)
 *     
 * A column c is identically aligned for sequences i, j if:
 *    1) both i,j have a symbol aligned in column c, and the
 *       same pair of symbols is aligned somewhere in the test
 *       alignment
 *    2) S[i][c] is aligned to a gap in sequence j, and that symbol
 *       is aligned to a gap in the test alignment
 *    3) converse of 2)
 *   
 *   
 * The algorithm is as follows:
 *    1) For each known/test aligned pair of sequences (k1,k2 and t1,t2)
 *        construct a list for each sequence, in which for every
 *        counted symbol we record the raw index of the symbol in
 *        the other sequence that it aligns to, or -1 if it aligns
 *        to a gap or uncounted symbol.
 *       
 *    2)  Compare the list for k1 to the list for t1 and count an identity
 *        for each correct alignment.
 *       
 *    3) Repeat 2) for comparing k2 to t2. Note that this means correct sym/sym
 *       alignments count for 2; correct sym/gap alignments count for 1.
 *   
 *    4) The score is (identities from 2 + identities from 3) /
 *       (totals from 2 + totals from 3).
 *
 * Written originally for koala's ss2 pairwise alignment package.
 *
 */


#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "squid.hpp"


/* Function: ComparePairAlignments
 *
 * Purpose:  Calculate and return a number representing how well two different alignments
 *           of a pair of sequences compare. The number is, roughly speaking,
 *           the fraction of columns which are identically aligned.
 *
 *           For all columns c in which either known1[c] or known2[c]
 *           is a non-gap, count an identity if those same symbols are
 *           aligned somewhere in calc1/calc2. The score is identities/total
 *           columns examined. (i.e. fully gapped columns don't count)
 *
 *           more explicitly, identities come from:
 *             both known and test aligned pairs have the same symbol in the first sequence aligned to
 *               a gap in the second sequence;
 *             both known and test aligned pairs have the same symbol in the second sequence
 *               aligned to a gap in the first sequence;
 *             the known alignment has symbols aligned at this column, and the test
 *               alignment aligns the same two symbols.
 *
 * Args:     known1, known2: trusted alignment of two sequences
 *           calc1, calc2:   test alignment of two sequences
 * 
 * Return:   Returns -1.0 on internal failure.
 */
float
ComparePairAlignments(
  char *known1,
  char *known2,
  char *calc1,
  char *calc2);


/* Function: CompareRefPairAlignments()
 *
 * Same as above, but the only columns that count are the ones
 * with indices in *refcoord. *refcoord and the known1, known2
 * pair must be in sync with each other (come from the same
 * multiple sequence alignment)
 *
 * Args:     ref           - 0..alen-1 array of 1 or 0
 *           known1,known2 - trusted alignment
 *           calc1, calc2  - test alignment          
 *
 * Return:  the fractional alignment identity on success, -1.0 on failure.
 */
float
CompareRefPairAlignments(
  int  *ref,
  char *known1,
  char *known2,
  char *calc1,
  char *calc2);


/* Function: make_alilist()
 *
 * Purpose:  Construct a list (array) mapping the raw symbols of s1
 *           onto the indexes of the aligned symbols in s2 (or -1
 *           for gaps in s2). The list (s1_list) will be of the
 *           length of s1's raw sequence.
 *          
 * Args:     s1          - sequence to construct the list for
 *           s2          - sequence s1 is aligned to
 *           ret_s1_list - RETURN: the constructed list (caller must free)
 *           ret_listlen - RETURN: length of the list
 *          
 * Returns:  1 on success, 0 on failure
 */
int
make_alilist(
  char *s1,
  char *s2,
  int **ret_s1_list,
  int *ret_listlen);


/* Function: make_ref_alilist()
 *
 * Purpose:  Construct a list (array) mapping the raw symbols of s1
 *           which are under canonical columns of the ref alignment
 *           onto the indexes of the aligned symbols in s2 (or -1
 *           for gaps in s2 or noncanonical symbols in s2).
 *          
 * Args:     ref:        - array of indices of canonical coords (1 canonical, 0 non)
 *           k1          - s1's known alignment (w/ respect to refcoords)
 *           k2          - s2's known alignment (w/ respect to refcoords)
 *           s1          - sequence to construct the list for
 *           s2          - sequence s1 is aligned to
 *           ret_s1_list - RETURN: the constructed list (caller must free)
 *           ret_listlen - RETURN: length of the list
 *          
 * Returns:  1 on success, 0 on failure
 */
/*ARGSUSED*/
int
make_ref_alilist(
  int *ref,
  char *k1,
  //char *k2,
  char *s1,
  char *s2,
  int **ret_s1_list,
  int *ret_listlen);


/* Function: compare_lists()
 *
 * Purpose:  Given four alignment lists (k1,k2, t1,t2), calculate the
 *           alignment score.
 *          
 * Args:     k1   - list of k1's alignment to k2
 *           k2   - list of k2's alignment to k1
 *           t1   - list of t1's alignment to t2
 *           t2   - list of t2's alignment to t2
 *           len1 - length of k1, t1 lists (same by definition)
 *           len2 - length of k2, t2 lists (same by definition)
 *           ret_sc - RETURN: identity score of alignment
 *
 * Return:   1 on success, 0 on failure.
 */          
bool
compare_lists(
  int *k1,
  int *k2,
  int *t1,
  int *t2,
  int len1,
  int len2,
  float *ret_sc);


/* Function: CompareMultAlignments
 *
 * Purpose:  Invokes pairwise alignment comparison for every possible pair,
 *           and returns the average score over all N(N-1) of them or -1.0
 *           on an internal failure.
 *
 *           Can be slow for large N, since it's quadratic.
 *
 * Args:     kseqs  - trusted multiple alignment
 *           tseqs  - test multiple alignment
 *           N      - number of sequences
 *          
 * Return:   average identity score, or -1.0 on failure.         
 */
float
CompareMultAlignments(
  char **kseqs,
  char **tseqs,
  int N);


/* Function: CompareRefMultAlignments()
 *
 * Purpose:  Same as above, except an array of reference coords for
 *           the canonical positions of the known alignment is also
 *           provided.
 *
 * Args:     ref      : 0..alen-1 array of 1/0 flags, 1 if canon
 *           kseqs    : trusted alignment
 *           tseqs    : test alignment
 *           N        : number of sequences
 *
 * Return:   average identity score, or -1.0 on failure
 */
float
CompareRefMultAlignments(
  int   *ref,
  char **kseqs,
  char **tseqs,
  int N);


/* Function: PairwiseIdentity()
 *
 * Purpose:  Calculate the pairwise fractional identity between
 *           two aligned sequences s1 and s2. This is simply
 *           (idents / MIN(len1, len2)).
 *
 *           Note how many ways there are to calculate pairwise identity,
 *           because of the variety of choices for the denominator:
 *           idents/(idents+mismat) has the disadvantage that artifactual
 *             gappy alignments would have high "identities".
 *           idents/(AVG|MAX)(len1,len2) both have the disadvantage that
 *             alignments of fragments to longer sequences would have
 *             artifactually low "identities".
 *          
 *           Case sensitive; also, watch out in nucleic acid alignments;
 *           U/T RNA/DNA alignments will be counted as mismatches!
 */
float
PairwiseIdentity(
  char *s1,
  char *s2);


/* Function: AlignmentIdentityBySampling()
 *
 * Purpose:  Estimate and return the average pairwise
 *           fractional identity of an alignment,
 *           using sampling.
 *          
 *           For use when there's so many sequences that
 *           an all vs. all rigorous calculation will
 *           take too long.
 *          
 *           Case sensitive!
 *
 * Args:     aseq       - aligned sequences
 *           L          - length of alignment
 *           N          - number of seqs in alignment
 *           nsample    - number of samples                    
 *
 * Returns:  average fractional identity, 0..1.
 */
float
AlignmentIdentityBySampling(
  char **aseq,
  //int L,
  int N,
  int nsample);


/* Function: MajorityRuleConsensus()
 *
 * Purpose:  Given a set of aligned sequences, produce a
 *           majority rule consensus sequence. If >50% nonalphabetic
 *           (usually meaning gaps) in the column, ignore the column.
 *
 * Args:     aseq  - aligned sequences, [0..nseq-1][0..alen-1]
 *           nseq  - number of sequences
 *           alen  - length of alignment       
 *
 * Returns:  ptr to allocated consensus sequence.
 *           Caller is responsible for free'ing this.
 */
char *
MajorityRuleConsensus(
  char **aseq,
  int nseq,
  int alen);

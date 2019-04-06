#pragma once

/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* weight.c
 *
 * Calculate weights for sequences in an alignment.
 */

#include <ctype.h>
#include <string.h>

#include "squid.hpp"


/* Function: GSCWeights()
 *
 * Purpose:  Use Erik's tree-based algorithm to set weights for
 *           sequences in an alignment. upweight() and downweight()
 *           are derived from Graeme Mitchison's code.
 *          
 * Args:     aseq        - array of (0..nseq-1) aligned sequences
 *           nseq        - number of seqs in alignment 
 *           alen        - length of alignment
 *           wgt         - allocated [0..nseq-1] array of weights to be returned
 *          
 * Return:   (void)
 *           wgt is filled in.
 */
void
GSCWeights(
  char **aseq,
  int nseq,
//  int alen,
  float *wgt);


void
upweight(
  struct phylo_s *tree,
  int nseq,
  float *lwt,
  float *rwt,
  int node);


void
downweight(
  struct phylo_s *tree,
  int nseq,
  float *lwt,
  float *rwt,
  float *fwt,
  int node);


/* Function: VoronoiWeights()
 *
 * Purpose:  Calculate weights using the scheme of Sibbald &
 *           Argos (JMB 216:813-818 1990). The scheme is
 *           slightly modified because the original algorithm
 *           actually doesn't work on gapped alignments.
 *           The sequences are assumed to be protein.
 *          
 * Args:     aseq  - array of (0..nseq-1) aligned sequences
 *           nseq  - number of sequences
 *           alen  - length of alignment
 *           wgt   - allocated [0..nseq-1] array of weights to be returned
 *
 * Return:   void
 *           wgt is filled in.
 */
void
VoronoiWeights(
  char **aseq,
  int nseq,
  int alen,
  float *wgt);


/* Function: simple_distance()
 *
 * Purpose:  For two identical-length null-terminated strings, return
 *           the fractional difference between them. (0..1)
 *           (Gaps don't count toward anything.)
 */
float
simple_distance(
  char *s1,
  char *s2);


/* Function: simple_diffmx()
 *
 * Purpose:  Given a set of flushed, aligned sequences, construct
 *           an NxN fractional difference matrix using the
 *           simple_distance rule.
 *          
 * Args:     aseqs        - flushed, aligned sequences
 *           num          - number of aseqs
 *           ret_dmx      - RETURN: difference matrix (caller must free)
 *          
 * Return:   1 on success, 0 on failure.
 */
int
simple_diffmx(
  char **aseqs,
  int num,
  float ***ret_dmx);


/* Function: BlosumWeights()
 *
 * Purpose:  Assign weights to a set of aligned sequences
 *           using the BLOSUM rule:
 *             - do single linkage clustering at some pairwise identity
 *             - in each cluster, give each sequence 1/clustsize
 *               total weight.
 *
 *           The clusters have no pairwise link >= maxid.
 *          
 *           O(N) in memory. Probably ~O(NlogN) in time; O(N^2)
 *           in worst case, which is no links between sequences
 *           (e.g., values of maxid near 1.0).
 *
 * Args:     aseqs - alignment
 *           nseq  - number of seqs in alignment
 *           alen  - # of columns in alignment
 *           maxid - fractional identity (e.g. 0.62 for BLOSUM62)
 *           wgt   - [0..nseq-1] array of weights to be returned
 */              
void
BlosumWeights(
  char **aseqs,
  int nseq,
  //int alen,
  float maxid,
  float *wgt);


/* Function: PositionBasedWeights()
 *
 * Purpose:  Implementation of Henikoff and Henikoff position-based
 *           weights (JMB 243:574-578, 1994) [Henikoff94b].
 *          
 *           A significant advantage of this approach that Steve and Jorja
 *           don't point out is that it is O(N) in memory, unlike
 *           many other approaches like GSC weights or Voronoi.
 *          
 *           A potential disadvantage that they don't point out
 *           is that in the theoretical limit of infinite sequences
 *           in the alignment, weights go flat: eventually every
 *           column has at least one representative of each of 20 aa (or 4 nt)
 *           in it.
 *          
 *           They also don't give a rule for how to handle gaps.
 *           The rule used here seems the obvious and sensible one
 *           (ignore them). This means that longer sequences
 *           initially get more weight; hence a "double
 *           normalization" in which the weights are first divided
 *           by sequence length (to compensate for that effect),
 *           then normalized to sum to nseq.
 *          
 * Limitations:
 *           Implemented in a way that's alphabet-independent:
 *           it uses the 26 upper case letters as "residues".
 *           Any alphabetic character in aseq is interpreted as
 *           a unique "residue" (case insensitively; lower case
 *           mapped to upper case). All other characters are
 *           interpreted as gaps.
 *          
 *           This way, we don't have to pass around any alphabet
 *           type info (DNA vs. RNA vs. protein) and don't have
 *           to deal with remapping IUPAC degenerate codes
 *           probabilistically. However, on the down side,
 *           a sequence with a lot of degenerate IUPAC characters
 *           will get an artifactually high PB weight.
 *          
 * Args:     aseq   - sequence alignment to weight
 *           nseq   - number of sequences in alignment
 *           alen   - length of alignment
 *           wgt    - RETURN: weights filled in (pre-allocated 0..nseq-1)
 *
 * Returns:  (void)
 *           wgt is allocated (0..nseq-1) by caller, and filled in here.
 */
void
PositionBasedWeights(
  char **aseq,
  int nseq,
  int alen,
  float *wgt);



/* Function: FilterAlignment()
 *
 * Purpose:  Constructs a new alignment by removing near-identical
 *           sequences from a given alignment (where identity is
 *           calculated *based on the alignment*).
 *           Does not affect the given alignment.
 *           Keeps earlier sequence, discards later one.
 *          
 *           Usually called as an ad hoc sequence "weighting" mechanism.
 *          
 * Limitations:
 *           Unparsed Stockholm markup is not propagated into the
 *           new alignment.
 *          
 * Args:     msa      -- original alignment
 *           cutoff   -- fraction identity cutoff. 0.8 removes sequences > 80% id.
 *           ret_new  -- RETURN: new MSA, usually w/ fewer sequences
 *                        
 * Return:   (void)
 *           ret_new must be free'd by caller: MSAFree().
 */
void
FilterAlignment(
  MSA *msa,
  float cutoff,
  MSA **ret_new);


/* Function: SampleAlignment()
 *
 * Purpose:  Constructs a new, smaller alignment by sampling a given
 *           number of sequences at random. Does not change the
 *           alignment nor the order of the sequences.
 *          
 *           If you ask for a sample that is larger than nseqs,
 *           it silently returns the original alignment.
 *          
 *           Not really a weighting method, but this is as good
 *           a place as any to keep it, since it's similar in
 *           construction to FilterAlignment().
 *          
 * Args:     msa      -- original alignment
 *           sample   -- number of sequences in new alignment (0 < sample <= nseq)
 *           ret_new  -- RETURN: new MSA
 *                        
 * Return:   (void)
 *           ret_new must be free'd by caller: MSAFree().
 */
void
SampleAlignment(
  MSA *msa,
  int sample,
  MSA **ret_new);


/* Function: SingleLinkCluster()
 *
 * Purpose:  Perform simple single link clustering of seqs in a
 *           sequence alignment. A pairwise identity threshold
 *           defines whether two sequences are linked or not.
 *          
 *           Important: runs in O(N) memory, unlike standard
 *           graph decomposition algorithms that use O(N^2)
 *           adjacency matrices or adjacency lists. Requires
 *           O(N^2) time in worst case (which is when you have
 *           no links at all), O(NlogN) in "average"
 *           case, and O(N) in best case (when there is just
 *           one cluster in a completely connected graph.
 *          
 *           (Developed because hmmbuild could no longer deal
 *           with GP120, a 16,013 sequence alignment.)
 *          
 * Limitations:
 *           CASE-SENSITIVE. Assumes aseq have been put into
 *           either all lower or all upper case; or at least,
 *           within a column, there's no mixed case.
 *          
 * Algorithm:
 *           I don't know if this algorithm is published. I
 *           haven't seen it in graph theory books, but that might
 *           be because it's so obvious that nobody's bothered.
 *          
 *           In brief, we're going to do a breadth-first search
 *           of the graph, and we're going to calculate links
 *           on the fly rather than precalculating them into
 *           some sort of standard adjacency structure.
 *          
 *           While working, we keep two stacks of maximum length N:
 *                a : list of vertices that are still unconnected.
 *                b : list of vertices that we've connected to
 *                    in our current breadth level, but we haven't
 *                    yet tested for other connections to a.
 *           The current length (number of elements in) a and b are
 *           kept in na, nb.
 *                   
 *           We store our results in an array of length N:
 *                c : assigns each vertex to a component. for example
 *                    c[4] = 1 means that vertex 4 is in component 1.
 *                    nc is the number of components. Components
 *                    are numbered from 0 to nc-1. We return c and nc
 *                    to our caller.
 *                   
 *           The algorithm is:
 *          
 *           Initialisation:
 *                a  <-- all the vertices
 *                na <-- N
 *                b  <-- empty set
 *                nb <-- 0
 *                nc <-- 0
 *               
 *           Then:
 *                while (a is not empty)
 *                  pop a vertex off a, push onto b
 *                  while (b is not empty)
 *                    pop vertex v off b
 *                    assign c[v] = nc
 *                    for each vertex w in a:
 *                       compare v,w. If w is linked to v, remove w
 *                       from a, push onto b.
 *                  nc++    
 *           q.e.d. :)      
 *
 * Args:     aseq   - aligned sequences
 *           nseq   - number of sequences in aseq
 *           alen   - alignment length
 *           maxid  - fractional identity threshold 0..1. if id >= maxid, seqs linked
 *           ret_c  - RETURN: 0..nseq-1 assignments of seqs to components (clusters)
 *           ret_nc - RETURN: number of components
 *
 * Returns:  void.
 *           ret_c is allocated here. Caller free's with free(*ret_c)
 */
void
SingleLinkCluster(
  char **aseq,
  int nseq,
  //int alen,
  float maxid,
  int **ret_c,
  int *ret_nc);
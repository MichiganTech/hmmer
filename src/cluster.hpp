#pragma once
/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* cluster.c
 *
 * almost identical to bord.c, from fd
 * also now contains routines for constructing difference matrices
 * from alignments
 * 
 * "branch ordering": Input a symmetric or upper-right-diagonal 
 * NxN difference matrix (usually constructed by pairwise alignment 
 * and similarity calculations for N sequences). Use the simple 
 * cluster analysis part of the Fitch/Margoliash tree-building algorithm
 * (as described by Fitch and Margoliash 1967 as well as Feng
 * and Doolittle 1987) to calculate the topology of an "evolutionary
 * tree" consistent with the difference matrix. Returns an array
 * which represents the tree.
 * 
 * The input difference matrix is just an NxN matrix of floats.
 * A good match is a small difference score (the algorithm is going
 * to search for minima among the difference scores). The original difference
 * matrix remains unchanged by the calculations.
 * 
 * The output requires some explanation. A phylogenetic
 * tree is a binary tree, with N "leaves" and N-1 "nodes". The
 * topology of the tree may be completely described by N-1 structures
 * containing two pointers; each pointer points to either a leaf
 * or another node. Here, this is implemented with integer indices
 * rather than pointers. An array of N-1 pairs of ints is returned.
 * If the index is in the range (0..N-1), it is a "leaf" -- the
 * number of one of the sequences. If the index is in the range
 * (N..2N-2), it is another "node" -- (index-N) is the index
 * of the node in the returned array.
 * 
 * If both indices of a member of the returned array point to
 * nodes, the tree is "compound": composed of more than one
 * cluster of related sequences.
 * 
 * The higher-numbered elements of the returned array were the
 * first constructed, and hence represent the distal tips
 * of the tree -- the most similar sequences. The root
 * is node 0.
 ******************************************************************
 *
 * Algorithm
 * 
 * INITIALIZATIONS:
 *  - copy the difference matrix (otherwise the caller's copy would
 *       get destroyed by the operations of this algorithm). If
 *       it's asymmetric, make it symmetric.
 *  - make a (0..N-1) array of ints to keep track of the indices in
 *       the difference matrix as they get swapped around. Initialize
 *       this matrix to 0..N-1.
 *  - make a (0..N-2) array of int[2] to store the results (the tree
 *       topology). Doesn't need to be initialized.
 *  - keep track of a "N'", the current size of the difference
 *       matrix being operated on.
 *
 * PROCESSING THE DIFFERENCE MATRIX:
 *  - for N' = N down to N' = 2  (N-1 steps):
 *    - in the half-diagonal N'xN' matrix, find the indices i,j at which
 *      there's the minimum difference score
 *      
 *     Store the results:
 *    - at position N'-2 of the result array, store coords[i] and 
 *         coords[j].
 *    
 *     Move i,j rows, cols to the outside edges of the matrix:
 *    - swap row i and row N'-2
 *    - swap row j and row N'-1   
 *    - swap column i and column N'-2
 *    - swap column j and column N'-1
 *    - swap indices i, N'-2 in the index array
 *    - swap indices j, N'-1 in the index array
 *    
 *     Build a average difference score for differences to i,j:
 *    - for all columns, find avg difference between rows i and j and store in row i: 
 *       row[i][col] = (row[i][col] + row[j][col]) / 2.0
 *    - copy the contents of row i to column i (it's a symmetric
 *       matrix, no need to recalculate)
 *    - store an index N'+N-2 at position N'-2 of the index array: means
 *       that this row/column is now a node rather than a leaf, and
 *       contains minimum values
 *       
 *     Continue:
 *    - go to the next N'
 *    
 * GARBAGE COLLECTION & RETURN.
 * 
 **********************************************************************
 *
 * References:
 * 
 * Feng D-F and R.F. Doolittle. "Progressive sequence alignment as a
 *    prerequisite to correct phylogenetic trees." J. Mol. Evol. 
 *    25:351-360, 1987.
 *    
 * Fitch W.M. and Margoliash E. "Construction of phylogenetic trees."
 *    Science 155:279-284, 1967.
 *    
 **********************************************************************
 */


#include <stdio.h>
#include <string.h>
#include <math.h>

#include "squid.h"
//#include "sqfuncs.h"


/* Function: Cluster()
 * 
 * Purpose:  Cluster analysis on a distance matrix. Constructs a
 *           phylogenetic tree which contains the topology
 *           and info for each node: branch lengths, how many
 *           sequences are included under the node, and which
 *           sequences are included under the node.
 *           
 * Args:     dmx     - the NxN distance matrix ( >= 0.0, larger means more diverged)
 *           N       - size of mx (number of sequences)
 *           mode    - CLUSTER_MEAN, CLUSTER_MAX, or CLUSTER_MIN
 *           ret_tree- RETURN: the tree 
 *
 * Return:   1 on success, 0 on failure.          
 *           The caller is responsible for freeing the tree's memory,
 *           by calling FreePhylo(tree, N).
 */
int
Cluster(
  float **dmx, 
  int N, 
  enum clust_strategy mode, 
  struct phylo_s **ret_tree);


/* Function: AllocPhylo()
 * 
 * Purpose:  Allocate space for a phylo_s array. N-1 structures
 *           are allocated, one for each node; in each node, a 0..N
 *           is_in flag array is also allocated and initialized to
 *           all zeros.
 *           
 * Args:     N  - size; number of sequences being clustered
 *                
 * Return:   pointer to the allocated array
 *           
 */
struct phylo_s*
AllocPhylo(
  int N);


/* Function: FreePhylo()
 * 
 * Purpose:  Free a clustree array that was built to cluster N sequences.
 * 
 * Args:     tree - phylogenetic tree to free
 *           N    - size of clustree; number of sequences it clustered
 *
 * Return:   (void)               
 */
void
FreePhylo(
  struct phylo_s *tree, 
  int N);


/* Function: MakeDiffMx()
 * 
 * Purpose:  Given a set of aligned sequences, construct
 *           an NxN fractional difference matrix. (i.e. 1.0 is
 *           completely different, 0.0 is exactly identical).
 *           
 * Args:     aseqs        - flushed, aligned sequences
 *           num          - number of aseqs
 *           ret_dmx      - RETURN: difference matrix 
 *           
 * Return:   1 on success, 0 on failure.
 *           Caller must free diff matrix with FMX2Free(dmx)
 */
void
MakeDiffMx(
  char **aseqs, 
  int num, 
  float ***ret_dmx);


/* Function: MakeIdentityMx()
 * 
 * Purpose:  Given a set of aligned sequences, construct
 *           an NxN fractional identity matrix. (i.e. 1.0 is
 *           completely identical, 0.0 is completely different).
 *           Virtually identical to MakeDiffMx(). It's
 *           less confusing to have two distinct functions, I find.
 *           
 * Args:     aseqs        - flushed, aligned sequences
 *           num          - number of aseqs
 *           ret_imx      - RETURN: identity matrix (caller must free)
 *           
 * Return:   1 on success, 0 on failure.
 *           Caller must free imx using FMX2Free(imx)
 */
void
MakeIdentityMx(
  char **aseqs, 
  int num, 
  float ***ret_imx);


/* Function: PrintNewHampshireTree()
 * 
 * Purpose:  Print out a tree in the "New Hampshire" standard
 *           format. See PHYLIP's draw.doc for a definition of
 *           the New Hampshire format.
 *
 *           Like a CFG, we generate the format string left to
 *           right by a preorder tree traversal.
 *           
 * Args:     fp   - file to print to
 *           ainfo- alignment info, including sequence names 
 *           tree - tree to print
 *           N    - number of leaves
 *           
 */
void
PrintNewHampshireTree(
  FILE *fp, 
  AINFO *ainfo, 
  struct phylo_s *tree, 
  int N);


/* Function: PrintPhylo()
 * 
 * Purpose:  Debugging output of a phylogenetic tree structure.
 */
void
PrintPhylo(
  FILE *fp, 
  AINFO *ainfo, 
  struct phylo_s *tree, 
  int N);
#pragma once

/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* hmmalign.c
 *
 * main() for aligning a set of sequences to an HMM.
 */

#include "config.h"    /* compile-time configuration constants */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "structs.h"    /* data structures, macros, #define's   */
#include "funcs.h"    /* function declarations                */
#include "globals.h"    /* alphabet global variables            */
#include "squid.h"    /* general sequence analysis library    */
#include "msa.h"    /* squid's multiple alignment i/o       */
#include "vectorops.h"
#include "getopt.h"






int
DealignedLength(
  char *aseq);


/* Function: GCGchecksum()
 *
 * Purpose:  Calculate a GCG checksum for a sequence.
 *           Code provided by Steve Smith of Genetics
 *           Computer Group.
 *
 * Args:     seq -  sequence to calculate checksum for.
 *                  may contain gap symbols.
 *           len -  length of sequence (usually known,
 *                  so save a strlen() call)       
 *
 * Returns:  GCG checksum.
 */
int
GCGchecksum(
	char *seq, 
	int len);


/* Function: GCGMultchecksum()
 * 
 * Purpose:  GCG checksum for a multiple alignment: sum of
 *           individual sequence checksums (including their
 *           gap characters) modulo 10000.
 *
 *           Implemented using spec provided by Steve Smith of
 *           Genetics Computer Group.
 *           
 * Args:     seqs - sequences to be checksummed; aligned or not
 *           nseq - number of sequences
 *           
 * Return:   the checksum, a number between 0 and 9999
 */                      
int
GCGMultchecksum(
	char **seqs, 
	int nseq);


/* Function: MSAToSqinfo()
 * Purpose:  Take an MSA and generate a SQINFO array suitable
 *           for use in annotating the unaligned sequences.
 *           Return the array.
 *           
 *           Permanent temporary code. sqinfo was poorly designed.
 *           it must eventually be replaced, but the odds
 *           of this happening soon are nil, so I have to deal.
 *
 * Args:     msa   - the alignment
 *
 * Returns:  ptr to allocated sqinfo array.
 *           Freeing is ghastly: free in each individual sqinfo[i] 
 *           with FreeSequence(NULL, &(sqinfo[i])), then
 *           free(sqinfo).
 */
SQINFO *
MSAToSqinfo(
  MSA *msa);


void
SeqinfoCopy(
  SQINFO *sq1, 
  SQINFO *sq2
);



/* Function: include_alignment()
 * Purpose:  Given the name of a multiple alignment file,
 *           align that alignment to the HMM, and add traces
 *           to an existing array of traces. If do_mapped
 *           is true, we use the HMM's map file. If not,
 *           we use P7ViterbiAlignAlignment().
 *
 * Args:     seqfile  - name of alignment file
 *           hmm      - model to align to
 *           do_mapped- true if we're to use the HMM's alignment map
 *           rsq      - RETURN: array of rseqs to add to
 *           dsq      - RETURN: array of dsq to add to
 *           sqinfo   - RETURN: array of SQINFO to add to
 *           tr       - RETURN: array of traces to add to
 *           nseq     - RETURN: number of seqs
 *
 * Returns:  new, realloc'ed arrays for rsq, dsq, sqinfo, tr; nseq is
 *           increased to nseq+ainfo.nseq.
 */
void
include_alignment(
	char *seqfile, 
	struct plan7_s *hmm, 
	int do_mapped,
	char ***rsq, 
	unsigned char ***dsq, 
	SQINFO **sqinfo,
	struct p7trace_s ***tr, 
	int *nseq);

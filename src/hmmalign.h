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

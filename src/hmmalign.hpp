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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "config.hpp"
#include "structs.hpp"    /* data structures, macros, #define's   */
#include "globals.hpp"    /* alphabet global variables            */
#include "squid.hpp"    /* general sequence analysis library    */
#include "msa.hpp"    /* squid's multiple alignment i/o       */
#include "vectorops.hpp"
#include "getopt.hpp"



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

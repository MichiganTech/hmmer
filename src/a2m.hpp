/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* a2m.c
 *
 * reading/writing A2M (aligned FASTA) files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.hpp"
#include "msa.hpp"

/* Function: ReadA2M()
 * Purpose:  Parse an alignment read from an open A2M format
 *           alignment file. A2M is a single alignment format.
 *           Return the alignment, or NULL if we've already
 *           read the alignment.
 *
 * Args:     afp - open alignment file
 *
 * Returns:  MSA *  - an alignment object.
 *                    Caller responsible for an MSAFree()
 */
MSA *
ReadA2M(
  MSAFILE *afp);


/* Function: WriteA2M()
 * Purpose:  Write an "aligned FASTA" (aka a2m, to UCSC) formatted
 *           alignment.
 *
 * Args:     fp    - open FILE to write to.
 *           msa   - alignment to write
 *
 * Returns:  void
 */
void
WriteA2M(
  FILE *fp,
  MSA *msa);
/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* phylip.c
 *
 * Import/export of PHYLIP interleaved multiple sequence alignment
 * format files.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "squid.hpp"
#include "msa.hpp"


/* Function: ReadPhylip()
 * Purpose:  Parse an alignment from an open Phylip format
 *           alignment file. Phylip is a single-alignment format.
 *           Return the alignment, or NULL if we have no data.
 *
 * Args:     afp - open alignment file
 *
 * Returns:  MSA * - an alignment object
 *                   Caller responsible for an MSAFree()
 *           NULL if no more alignments       
 */
MSA *
ReadPhylip(
  MSAFILE *afp);



/* Function: WritePhylip()
 * Purpose:  Write an alignment in Phylip format to an open file.
 *
 * Args:     fp    - file that's open for writing.
 *           msa   - alignment to write.
 *
 * Returns:  (void)
 */
void
WritePhylip(
  FILE *fp,
  MSA *msa);
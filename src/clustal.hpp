/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* clustal.c
 *
 * Import/export of ClustalV/W multiple sequence alignment
 * formatted files. Derivative of msf.c; MSF is a pretty
 * generic interleaved format.  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "squid.hpp"
#include "msa.hpp"



/* Function: ReadClustal()
 * Purpose:  Parse an alignment read from an open Clustal format
 *           alignment file. Clustal is a single-alignment format.
 *           Return the alignment, or NULL if we have no data.
 *          
 * Args:     afp  - open alignment file
 *
 * Returns:  MSA * - an alignment object
 *                   caller responsible for an MSAFree()
 *           NULL if no more alignments
 *
 * Diagnostics:
 *           Will Die() here with a (potentially) useful message
 *           if a parsing error occurs.
 */
MSA *
ReadClustal(
  MSAFILE *afp);


/* Function: WriteClustal()
 * Purpose:  Write an alignment in Clustal format to an open file.
 *
 * Args:     fp    - file that's open for writing.
 *           msa   - alignment to write.
 *
 * Returns:  (void)
 */
void
WriteClustal(
  FILE *fp,
  MSA *msa);
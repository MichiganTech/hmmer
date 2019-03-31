/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* msf.c
 * Import/export of GCG MSF multiple sequence alignment
 * formatted files. Designed using format specifications
 * kindly provided by Steve Smith of Genetics Computer Group.
 */

#include <stdio.h>

#include "msa.h"


/* Function: ReadMSF()
 * Purpose:  Parse an alignment read from an open MSF format
 *           alignment file. (MSF is a single-alignment format.)
 *           Return the alignment, or NULL if we've already
 *           read the alignment.
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
MSA*
ReadMSF(
  MSAFILE *afp);


/* Function: WriteMSF()
 * Purpose:  Write an alignment in MSF format to an open file.
 *
 * Args:     fp    - file that's open for writing.
 *           msa   - alignment to write. 
 *
 *                   Note that msa->type, usually optional, must be
 *                   set for WriteMSF to work. If it isn't, a fatal
 *                   error is generated.
 *
 * Returns:  (void)
 */
void
WriteMSF(
  FILE *fp, 
  MSA *msa);
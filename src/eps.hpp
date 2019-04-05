/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* eps.c
 * Some crude support for Encapsulated PostScript (EPS) output,
 * DSC compliant.
 */

#pragma once

#include <stdlib.h>

#include "msa.h"       

/* Function: EPSWriteSmallMSA()
 * Purpose:  Write an alignment in singleblock, Stockholm/SELEX like
 *           format to an open file. Very crude.
 *           Currently fails if the alignment is >50 columns long, because
 *           it doesn't think it will fit on a single page.
 *
 * Args:     fp  - open file for writing
 *           msa - alignment to write     
 *
 * Returns:  (void)
 */
void
EPSWriteSmallMSA(
  FILE *fp, 
  MSA *msa);
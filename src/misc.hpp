/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* misc.c
 *
 * Functions that I don't know quite where to put yet.
 */


#include <stdio.h>

#include "config.hpp"
#include "getopt.hpp"
#include "squid.hpp"
#include "structs.hpp"


/* Function: HMMERBanner()
 *
 * Purpose:  Print a package version and copyright banner.
 *           Used by all the main()'s.
 *
 *    Expects to be able to pick up defined preprocessor variables:
 *    variable          example
 *    --------           --------------
 *    PACKAGE_NAME      "HMMER"
 *    PACKAGE_VERSION   "2.0.42"
 *    PACKAGE_DATE      "April 1999"
 *    PACKAGE_COPYRIGHT "Copyright (C) 1992-1999 Washington University School of Medicine"
 *    PACKAGE_LICENSE   "Freely distributed under the GNU General Public License (GPL)."
 *
 *    This gives us a general mechanism to update release information
 *    without changing multiple points in the code.
 *
 * Args:     fp     - where to print it
 *           banner - one-line program description, e.g.:
 *                    "foobar - make bars from foo with elan"
 * Returns:  (void)
 */
void
HMMERBanner(FILE *fp, char *banner);


/* Function: Getword()
 *
 * Purpose:  little function used by ReadPrior() and ReadHMM() to parse
 *           next valid field out of an open file, ignoring
 *           comments. '#' marks the beginning of a comment.
 *
 * Arg:      fp   - open file for reading
 *           type - sqdARG_INT, sqdARG_FLOAT, or sqdARG_STRING from squid.h
 */
char*
Getword(
  FILE *fp,
  enum arg_type type);


/* Function: SetAutocuts()
 *
 * Purpose:  Set score thresholds using the GA, TC, or NC information
 *           in an HMM.
 *
 * Args:     thresh - score threshold structure. autocut must be set
 *                    properly (CUT_GA, CUT_NC, or CUT_TC).
 *           hmm    - HMM containing appropriate score cutoff info
 *
 * Returns:  1 on success.
 *           0 if HMM does not have the score cutoffs available -- caller
 *             will have to decide on a fallback plan.
 *           Has no effect (and returns success) if autocut is
 *           CUT_NONE.
 */
int
SetAutocuts(
  struct threshold_s *thresh,
  struct plan7_s *hmm);

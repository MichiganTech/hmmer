/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* debug.c
 *
 * Printing out or naming various useful things from HMMER
 * innards.
 */

#include <stdio.h>

#include "structs.hpp"


/* Function: Statetype()
 *
 * Purpose:  Returns the state type in text.
 * Example:  Statetype(S) = "S"
 */
char*
Statetype(char st);


/* Function: AlphabetType2String()
 * Purpose:  Returns a string "protein" for hmmAMINO,
 *           "nucleic acid" for hmmNUCLEIC, etc... used
 *           for formatting diagnostics.
 *
 * Args:     type - Alphabet type, e.g. hmmAMINO
 *
 * Returns:  char *
 */
char*
AlphabetType2String(
  int type);


/* Function: P7PrintTrace()
 *
 * Purpose:  Print out a traceback structure.
 *           If hmm is non-NULL, also print transition and emission scores.
 *
 * Args:     fp  - stderr or stdout, often
 *           tr  - trace structure to print
 *           hmm - NULL or hmm containing scores to print
 *           dsq - NULL or digitized sequence trace refers to.
 */
void
P7PrintTrace(
  FILE *fp,
  struct p7trace_s *tr,
  struct plan7_s *hmm,
  unsigned char *dsq);


/* Function: TraceVerify()
 * Purpose:  Check a traceback structure for internal consistency.
 *           Used in Shiva testsuite, for example.
 *
 * Args:     tr  - traceback to verify
 *           M   - length of HMM
 *           N   - length of sequence
 *
 * Returns:  1 if OK. 0 if not.
 */
bool
TraceVerify(
  struct p7trace_s *tr,
  int M,
  int N);


/* Function: TraceCompare()
 * Purpose:  Compare two tracebacks; return 1 if they're
 *           identical, else 0. Written for Shiva testsuite.
 *
 * Args:     t1 - first trace
 *           t2 - second trace
 *
 * Returns:  1 if identical; 0 elsewise
 */
bool
TraceCompare(
  struct p7trace_s *t1,
  struct p7trace_s *t2);
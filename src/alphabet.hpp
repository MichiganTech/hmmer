/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* alphabet.c
 * Configuration of the global symbol alphabet information.
 */

//#include "config.h"
//#include "squidconf.h"

//#include <stdlib.h>
//#include <string.h>
//#include <ctype.h>
//#include <pthread.h>

//#include "structs.h"
//#include "squid.h"
#include "msa.hpp"


/* Function: DetermineAlphabet()
 *
 * Purpose:  From a set of sequences (raw or aligned), make a good
 *           guess whether they're Nucleic, Amino, or something
 *           else, and set alphabet accordingly.
 *
 *           If Alphabet_type is already set, that means our
 *           autodetection was overridden from the command line,
 *           and we just set the other globals accordingly.
 */
void
DetermineAlphabet(
  char **rseqs,
  int  nseq);


/* Function: SetAlphabet()
 *
 * Purpose:  Set the alphabet globals, given an alphabet type
 *           of either hmmAMINO or hmmNUCLEIC.
 */
void
SetAlphabet(
  int type);


/* Function: SymbolIndex()
 *
 * Purpose:  Convert a symbol to its index in Alphabet[].
 *           Bogus characters are converted to 'X'.
 *           More robust than the SYMIDX() macro but
 *           presumably slower.
 */
unsigned char
SymbolIndex(
  char sym);


/* Function: DigitizeSequence()
 *
 * Purpose:  Internal representation of a sequence in HMMER is
 *           as a char array. 1..L are the indices
 *           of seq symbols in Alphabet[]. 0,L+1 are sentinel
 *           bytes, set to be Alphabet_iupac -- i.e. one more
 *           than the maximum allowed index.
 *
 *           Assumes that 'X', the fully degenerate character,
 *           is the last character in the allowed alphabet.
 *
 * Args:     seq - sequence to be digitized (0..L-1)
 *           L   - length of sequence
 *
 * Return:   digitized sequence, dsq.
 *           dsq is allocated here and must be free'd by caller.
 */
unsigned char *
DigitizeSequence(
  char *seq,
  int L);


/* Function: DedigitizeSequence()
 * Purpose:  Returns a 0..L-1 character string, converting the
 *           dsq back to the real alphabet.
 */
char *
DedigitizeSequence(
  unsigned char *dsq,
  int L);


/* Function: DigitizeAlignment()
 *
 * Purpose:  Given an alignment, return digitized unaligned
 *           sequence array. (Tracebacks are always relative
 *           to digitized unaligned seqs, even if they are
 *           faked from an existing alignment in modelmakers.c.)
 *
 * Args:     msa      - alignment to digitize
 *           ret_dsqs - RETURN: array of digitized unaligned sequences
 *
 * Return:   (void)
 *           dsqs is alloced here. Free2DArray(dseqs, nseq).
 */
void
DigitizeAlignment(
  MSA *msa,
  unsigned char ***ret_dsqs);


/* Function: P7CountSymbol()
 *
 * Purpose:  Given a possibly degenerate symbol code, increment
 *           a symbol counter array (generally an emission
 *           probability vector in counts form) appropriately.
 *
 * Args:     counters:  vector to count into. [0..Alphabet_size-1]
 *           symidx:    symbol index to count: [0..Alphabet_iupac-1]
 *           wt:        weight to use for the count; often 1.0
 *
 * Return:   (void)
 */
void
P7CountSymbol(
  float *counters,
  unsigned char symidx,
  float wt);


/* Function: set_degenerate()
 *
 * Purpose:  convenience function for setting up
 *           Degenerate[][] global for the alphabet.
 */
void
set_degenerate(
  char iupac,
  char *syms);

/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* sre_string.c
 * 
 * my library of extra string functions. Some for portability
 * across UNIXes
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include "squid.h"


/* Function: StringChop()
 * 
 * Purpose:  Chop trailing whitespace off of a string.
 */
void
StringChop(
  char *s);

int
Strinsert(
  char  *s1,            /* string to insert a char into  */
  char   c,		/* char to insert                */
  int    pos);		/* position in s1 to insert c at */


int
Strdelete(
  char *s1,             /* string to delete a char from       */
  int   pos);		/* position of char to delete 0..n-1  */


void
s2lower(
  char *s);


void
s2upper(
  char *s);


/* Function: RandomSequence()
 * 
 * Purpose:  Generate an iid symbol sequence according
 *           to some alphabet, alphabet_size, probability
 *           distribution, and length. Return the
 *           sequence.
 *           
 * Args:     alphabet  - e.g. "ACGT"
 *           p         - probability distribution [0..n-1]
 *           n         - number of symbols in alphabet
 *           len       - length of generated sequence 
 *           
 * Return:   ptr to random sequence, or NULL on failure.
 */
char *
RandomSequence(
  char *alphabet, 
  float *p, 
  int n, 
  int len);

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
 * RCS $Id: sre_string.c,v 1.11 2001/06/07 16:59:37 eddy Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include "sre_string.hpp"
#include "vectorops.hpp"


void
StringChop(
  char *s
){
  int   i;

  i = strlen(s) - 1;		         /* set i at last char in string     */
  while (i >= 0 && isspace((int) s[i])) i--;   /* i now at last non-whitespace char, or -1 */
  s[i+1] = '\0';
}


int
Strinsert(
  char  *s1,            /* string to insert a char into  */
  char   c,		/* char to insert                */
  int    pos /* position in s1 to insert c at */
){		
  char    oldc;
  char   *s;

  for (s = s1 + pos; c; s++)
    {
				/* swap current char for inserted one */
      oldc = *s;		/* pick up current */
      *s   = c;   		/* put down inserted one    */
      c    = oldc;		/* old becomes next to insert */
    }
  *s = '\0';

  return 1;
}


int
Strdelete(
  char *s1,             /* string to delete a char from       */
  int   pos /* position of char to delete 0..n-1  */
){		
  char *s;                     

  for (s = s1 + pos; *s; s++)
    *s = *(s + 1);

  return 1;
}


void
s2lower(
  char *s
){
  for (; *s != '\0'; s++)
    *s = tolower((int) *s);
}


void
s2upper(
  char *s
){
  for (; *s != '\0'; s++)
    *s = toupper((int) *s);
}


char *
RandomSequence(
  char *alphabet,
  float *p,
  int n,
  int len
){
  char *s;
  int   x;

  s = (char *) MallocOrDie (sizeof(char) * (len+1));
  for (x = 0; x < len; x++)
    s[x] = alphabet[FChoose(p,n)];
  s[x] = '\0';
  return s;
}


#ifdef CUBS_WIN
/* A timing test for sre_strcat()
 * cc -O2 -g sre_string.c sre_ctype.c sqerror.c sre_math.c hsregex.c -lm
 * 15.200u - 5.360u = 9.84u if sre_strcat() with no length info passed
 * 13.660u - 5.360u = 8.30u if strcat(), with a single malloc().
 * 11.370u - 5.360u = 6.01u if sre_strcat() with length info passed.
 */
int main(
){
  float p[4] = {0.25, 0.25, 0.25, 0.25};
  int   buflen;
  int   len;
  int   nappends;
  int   nstrings;
  char *s1 = NULL;
  char *s2;
  int   i;

  nappends = 100;
  nstrings = 1000;
  while (nstrings--)
    {
      /* s1 = malloc(sizeof(char) * (255*nappends+1));
	 s1[0] = '\0';
      */

      s1 = NULL;
      len = 0;
      for (i = 0; i < nappends; i++)
	{
	  buflen = CHOOSE(255) + 1;
	  s2 = RandomSequence("ACGT", p, 4, buflen);
     
	  /* strcat(s1,s2); */
	  if ((len = sre_strcat(&s1, len, s2, buflen)) < 0) exit(1);
	  free(s2);
	}
      free(s1);
    }
  exit(0);
}
#endif /*CUBS_WIN*/

/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/


#include "squid.h"


/* Library version info is made available as a global to
 * any interested program. These are defined in iupac.c
 * with the other globals.
  */
char* squid_version;  /* version number  */
char* squid_date; /* date of release */
int  squid_errno = 0;  /* error codes     */


struct iupactype *iupac;

char**  stdcode1; /* 1-letter amino acid translation code */
char**  stdcode3; /* 3-letter amino acid translation code */
float*  dnafq;        /* nucleotide occurrence frequencies    */
float*  aafq;   /* amino acid occurrence frequencies    */
char*   aa_alphabet;  /* amino acid alphabet                  */
int*    aa_index;     /* convert 0..19 indices to 0..26       */


/* Strparse() defines and manages these. 
 * sqd_parse[0] contains the substring that matched the pattern.
 * sqd_parse[1-9] contain substrings matched with ()'s.
 */
char *sqd_parse[10];


void
Free2DArray(
  void **p, 
  int dim1
){
  int i;
  
  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) free(p[i]);
    free(p);
  }
}


void
Free3DArray(
  void ***p, 
  int dim1, 
  int dim2
){
  int i, j;

  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) {
  for (j = 0; j < dim2; j++)
    if (p[i][j] != NULL) free(p[i][j]);
  free(p[i]);
      }
    free(p);
  }
}



bool
FileExists(
  char *filename
){
  FILE *fp;
  if ((fp = fopen(filename, "r"))) { fclose(fp); return true; }
  return false;
}
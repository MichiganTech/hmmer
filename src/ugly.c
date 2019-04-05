/* Generated automatically from squid.h.in by configure. */
/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/


#include "squid.h"
#include "alignio.h"

 char *aminos      = "ABCDEFGHIKLMNPQRSTVWXYZ*";
 char *primenuc    = "ACGTUN";
 char *protonly    = "EFIPQZ";


char *sqd_parse[10];


int
Strparse(
  char *rexp, 
  char *s, 
  int ntok
){
  sqd_regexp *pat;
  int         code;
  int         len;
  int         i;
        /* sanity check */
  if (ntok >= NSUBEXP )  Die("Strparse(): ntok must be <= %d", NSUBEXP-1); 

  /* Free previous global substring buffers
   */
  for (i = 0; i <= ntok; i++)
    if (sqd_parse[i] != NULL) 
      { 
  free(sqd_parse[i]);
  sqd_parse[i] = NULL;
      }

  /* Compile and match the pattern, using our modified 
   * copy of Henry Spencer's regexp library
   */
  if ((pat = sqd_regcomp(rexp)) == NULL) 
    Die("regexp compilation failed.");
  code = sqd_regexec(pat, s);

  /* Fill the global substring buffers
   */
  if (code == 1) 
    for (i = 0; i <= ntok; i++)
      if (pat->startp[i] != NULL && pat->endp[i] != NULL)
  {
    len = pat->endp[i] - pat->startp[i];
    sqd_parse[i] = (char *) MallocOrDie(sizeof(char) * (len+1));
    strncpy(sqd_parse[i], pat->startp[i], len);
    sqd_parse[i][len] = '\0';
  }

  free(pat);
  return code;
}


void *
sre_malloc(
  char *file, 
  int line, 
  size_t size
){
  void *ptr;

  if ((ptr = malloc (size)) == NULL)
    Die("malloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}


void *
sre_realloc(
  char *file, 
  int line, 
  void *p, 
  size_t size
){
  void *ptr;

  if ((ptr = realloc(p, size)) == NULL)
    Die("realloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}


bool
IsBlankline(
  char *s
){
  for (; *s != '\0'; s++)
    if (! isspace(*s)) return false;
  return true;
}


bool
IsInt(
  char *s
){
  int hex = 0;

  if (s == NULL) {squid_errno = SQERR_PARAMETER; return false; }

        /* skip whitespace */
  while (isspace((int) (*s))) s++;      
        /* skip leading sign */
  if (*s == '-' || *s == '+') s++;
        /* skip leading conversion signals */
  if ((strncmp(s, "0x", 2) == 0 && (int) strlen(s) > 2) ||
      (strncmp(s, "0X", 2) == 0 && (int) strlen(s) > 2))
    {
      s += 2;
      hex = 1;
    }
  else if (*s == '0' && (int) strlen(s) > 1)
    s++;
        /* examine remainder for garbage chars */
  if (!hex)
    while (*s != '\0')
      {
  if (!isdigit((int) (*s))) return false;
  s++;
      }
  else
    while (*s != '\0')
      {
  if (!isxdigit((int) (*s))) return false;
  s++;
      }

  return true;
}


int
ParsePAMFile(
  FILE *fp, 
  int ***ret_pam, 
  float *ret_scale
){
  int    **pam;
  char     buffer[512];   /* input buffer from fp                  */
  int      order[27];   /* order of fields, obtained from header */
  int      nsymbols;    /* total number of symbols in matrix     */
  char    *sptr;
  int      idx;
  int      row, col;
  float    scale;
  int      gotscale = false;
  
  scale = 0.0;    /* just to silence gcc uninit warnings */
  if (fp == NULL) { squid_errno = SQERR_NODATA; return 0; }
  
  /* Look at the first non-blank, non-comment line in the file.
   * It gives single-letter codes in the order the PAM matrix
   * is arrayed in the file. 
   */
  do {
    if (fgets(buffer, 512, fp) == NULL) 
      { squid_errno = SQERR_NODATA; return 0; }

    /* Get the scale factor from the header.
     * For BLOSUM files, we assume the line looks like:
     *     BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
     * and we assume that the fraction is always 1/x;
     * 
     * For PAM files, we assume the line looks like:
     *     PAM 120 substitution matrix, scale = ln(2)/2 = 0.346574
     * and we assume that the number following the final '=' is our scale
     */
    if (strstr(buffer, "BLOSUM Clustered Scoring Matrix") != NULL &&
  (sptr = strchr(buffer, '/')) != NULL)
      {
  sptr++;
  if (! isdigit((int) (*sptr))) { squid_errno = SQERR_FORMAT; return 0; }
  scale = (float) (log(2.0) / atof(sptr));
  gotscale = true;
      }
    else if (strstr(buffer, "substitution matrix,") != NULL)
      {
  while ((sptr = strrchr(buffer, '=')) != NULL) {
    sptr += 2;
    float tmpFloat;
    if (1 == sscanf(sptr, "%f", &tmpFloat)) {
      scale = atof(sptr);
      gotscale = true;
      break;
    }
  }
      }
  } while ((sptr = strtok(buffer, " \t\n")) == NULL || *sptr == '#');

  idx = 0;
  do {
    order[idx] = (int) *sptr - (int) 'A';
    if (order[idx] < 0 || order[idx] > 25) order[idx] = 26;
    idx++;
  } while ((sptr = strtok(NULL, " \t\n")) != NULL);
  nsymbols = idx;
  
  /* Allocate a pam matrix. For speed of indexing, we use
   * a 27x27 matrix so we can do lookups using the ASCII codes
   * of amino acid single-letter representations, plus one
   * extra field to deal with the "*" (terminators).
   */
  if ((pam = (int **) calloc (27, sizeof(int *))) == NULL)
    Die("calloc failed");
  for (idx = 0; idx < 27; idx++)
    if ((pam[idx] = (int *) calloc (27, sizeof(int))) == NULL)
      Die("calloc failed");

  /* Parse the rest of the file.
   */
  for (row = 0; row < nsymbols; row++)
    {
      if (fgets(buffer, 512, fp) == NULL) 
  { squid_errno = SQERR_NODATA; return 0; }

      if ((sptr = strtok(buffer, " \t\n")) == NULL)
  { squid_errno = SQERR_NODATA; return 0; }
      for (col = 0; col < nsymbols; col++)
  {
    if (sptr == NULL) { squid_errno = SQERR_NODATA; return 0; }

    /* Watch out for new BLAST format, with leading characters
     */
    if (*sptr == '*' || isalpha((int) *sptr))
      col--;  /* hack hack */
    else
      pam [order[row]] [order[col]] = atoi(sptr);

    sptr = strtok(NULL, " \t\n");
  }
    }
  
  /* Return
   */
  if (ret_scale != NULL)
    {
      if (gotscale) *ret_scale = scale;
      else
  {
    Warn("Failed to parse PAM matrix scale factor. Defaulting to ln(2)/2!");
    *ret_scale = log(2.0) / 2.0;
  }
    }
  *ret_pam = pam;
  return 1;
}


void
SqdClean(
){
  int i;

  // Free global substring buffers that Strparse() uses
  for (i = 0; i <= 9; i++)
    if (sqd_parse[i] != NULL) {
      free(sqd_parse[i]);
      sqd_parse[i] = NULL;
    }
}
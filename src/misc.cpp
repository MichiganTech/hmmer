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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>

#include "misc.hpp"


void
HMMERBanner(
  FILE *fp,
  char *banner
){
  fprintf(fp, "%s\n", banner);
  fprintf(fp, "%s %s (%s)\n", PACKAGE_NAME, PACKAGE_VERSION, PACKAGE_DATE);
  fprintf(fp, "%s\n", PACKAGE_COPYRIGHT);
  fprintf(fp, "%s\n", PACKAGE_LICENSE);
  fprintf(fp, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
}


char*
Getword(
  FILE *fp,
  enum arg_type type
){
  static char buffer[512];
  static char *sptr = NULL;


  if ((sptr = fgets(buffer, 512, fp)) == NULL) return NULL;
  if ((sptr = strchr(buffer, '#')) != NULL) return NULL;

  char tmpStr[512];
  float tmpFloat;
  int tmpInt;

  switch (type) {
  case sqdARG_STRING:
    if (sscanf(sptr, "%s", tmpStr) == 0) {
      Warn("Parse failed: expected string, got nothing");
      return NULL;
    }
    break;
  case sqdARG_INT:
    if (sscanf(sptr, "%d", &tmpInt) == 0) {
      Warn("Parse failed: expected integer, got %s", sptr);
      return NULL;
    }
    break;
  case sqdARG_FLOAT:
    if (sscanf(sptr, "%f", &tmpFloat) == 0) {
      Warn("Parse failed: expected real value, got %s", sptr);
      return NULL;
    }
    break;
  default: break;
  }

  return sptr;
}


int
SetAutocuts(
  struct threshold_s *thresh,
  struct plan7_s *hmm
){
  if (thresh->autocut == CUT_GA) {
    if (! (hmm->flags & PLAN7_GA)) return 0;
    thresh->globT = hmm->ga1;
    thresh->domT  = hmm->ga2;
    thresh->globE = thresh->domE = FLT_MAX;
  } else if (thresh->autocut == CUT_NC) {
    if (! (hmm->flags & PLAN7_NC)) return 0;
    thresh->globT = hmm->nc1;
    thresh->domT  = hmm->nc2;
    thresh->globE = thresh->domE = FLT_MAX;
  } else if (thresh->autocut == CUT_TC) {
    if (! (hmm->flags & PLAN7_TC)) return 0;
    thresh->globT = hmm->tc1;
    thresh->domT  = hmm->tc2;
    thresh->globE = thresh->domE = FLT_MAX;
  }
  return 1;
}

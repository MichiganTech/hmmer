/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* hmmfetch.c
 *
 * Recover a specific HMM file from an HMM database, using
 * an SSI index (created with hmmindex).
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.hpp"
#include "structs.hpp"
#include "globals.hpp"
#include "getopt.hpp"


static char banner[] = "hmmfetch -- retrieve specific HMM from an HMM database";

static char usage[] = "\
Usage: hmmfetch [-options] <hmmfile> <HMM name>\n\
Available options are:\n\
  -h   : print short usage and version info, then exit\n\
  -n   : interpret <HMM name> instead as an HMM number (0..nhmm-1)\n\
";

static char experts[] = "\
";

static struct opt_s OPTIONS[] = {
  { "-h", true, sqdARG_NONE  },
  { "-n", true, sqdARG_NONE  },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


int
main(int argc, char **argv) {
  char    *hmmfile;             /* HMM file to open                */
  char    *key;      /* HMM name to retrieve            */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  struct plan7_s *hmm;    /* a hidden Markov model           */

  char *optname;    /* name of option found by Getopt() */
  char *optarg;      /* argument found by Getopt()       */
  int   optind;            /* index in argv[]                  */

  int   by_number;    /* fetch by number, not name        */

  /***********************************************
   * Parse the command line
   ***********************************************/

  by_number = false;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg)) {
    if      (strcmp(optname, "-n") == 0) by_number = true;
    else if (strcmp(optname, "-h") == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];
  key     = argv[optind++];

  /***********************************************
   * Open HMM file, make sure SSI index exists
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);
  if (hmmfp->ssi == NULL)
    Die("There is no SSI index for %s; you need to use hmmindex on it.", hmmfile);

  /***********************************************
   * find key in hmmfile; get HMM; show as ASCII
   ***********************************************/

  if (by_number) {
    if (! IsInt(key)) Die("%s does not appear to be a number.", key);
    int   nhmm;      /* hmm number */
    nhmm = atoi(key);
    if (! HMMFilePositionByIndex(hmmfp, nhmm))
      Die("failed to position %s to HMM #%d", hmmfile, nhmm);
  } else {
    if (! HMMFilePositionByName(hmmfp, key))
      Die("No such hmm %s in HMM file %s\n", key, hmmfile);
  }

  if (! HMMFileRead(hmmfp, &hmm))
    Die("Unexpected end of HMM file");
  if (hmm == NULL)
    Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);

  WriteAscHMM(stdout, hmm);

  FreePlan7(hmm);
  HMMFileClose(hmmfp);

  /***********************************************
   * Exit
   ***********************************************/

  SqdClean();
  return 0;
}
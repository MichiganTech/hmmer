/* masks_test.c
 *
 * Test driver for sequence masking routines in masks.c
 *
 */


#include <stdio.h>
#include <string.h>

#include "config.hpp"
#include "structs.hpp"
#include "globals.hpp"

static char banner[] = "\
masks_test : testing of repeat masking code in masks.c";

static char usage[] = "\
Usage: testdriver [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -v              : verbose output\n\
";

static char experts[] = "\
  --xnu     <file>: apply xnu to seqs in <file>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       true,  sqdARG_NONE  },
  { "-v",       true,  sqdARG_NONE  },
  { "--xnu",    false, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

/* The test sequence and result from the XNU software distribution
 */
static char *test1 = "\
ACDEFGHIKLMNPQRQRQRQRQRQRQRQRQRSTVWYACDEFGHIKLMNPQRQRQRQRQRQ\
RQRQRQRSTVWYACDEFGHIKLMNPQRQRQRQRQRQRQRQRQRSTVWYACDEFGHIKLMN\
PQRQRQRQRQRQRQRQRQRSTVWYACDEFGHIKLMNPQRQRQRQRQRQRQRQRQRSTVWY\
ACDEFGHIKLMNPQRQRQRQRQRQRQRQRQRSTVWY";

static char *answer1 = "\
ACDEFGHIKLMNPXXXXXXXXXXXXXXXXXXSTVWYACDEFGHIKLMNPXXXXXXXXXXX\
XXXXXXXSTVWYACDEFGHIKLMNPXXXXXXXXXXXXXXXXXXSTVWYACDEFGHIKLMN\
PXXXXXXXXXXXXXXXXXXSTVWYACDEFGHIKLMNPXXXXXXXXXXXXXXXXXXSTVWY\
ACDEFGHIKLMNPXXXXXXXXXXXXXXXXXXSTVWY";

int
main(int argc, char **argv) {
  char          *seq;
  unsigned char *dsq;
  int   len;
  int   i,j;
  char *result;

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   be_verbose;
  char *xnufile;    /* NULL, or file to run xnu on     */


  /***********************************************
   * Parse command line
   ***********************************************/

  be_verbose = false;
  xnufile    = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-v")       == 0) {
      be_verbose = true;
    } else if (strcmp(optname, "--xnu")    == 0) {
      xnufile    = optarg;
    } else if (strcmp(optname, "-h")       == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 0)
    Die("Incorrect number of arguments.\n%s\n", usage);

  SetAlphabet(hmmAMINO);

  /* XNU test
   */
  seq = test1;
  len = (int) strlen(seq);
  dsq = DigitizeSequence(seq, len);
  XNU(dsq, len);
  result = MallocOrDie(sizeof(char) * (len+1));

  for (i = 0; i < len; i++)
    result[i] = Alphabet[dsq[i+1]];
  result[len] = '\0';

  if (be_verbose) {
    printf("XNU test:\n");
    for (i = 1; i <= len; i+=60) {
      for (j = i; j < i+60 && j <= len; j++)
        putc(Alphabet[(int) dsq[j]], stdout);
      putc('\n', stdout);
    }
    if (strcmp(answer1, result) == 0)
      printf("-- OK; Identical to expected\n");
  }

  if (strcmp(answer1, result) != 0)
    Die("XNU test failed.");
  free(result);
  free(dsq);

  /* On demand XNU test.
   */
  if (xnufile != NULL) {
    SQFILE *sqfp;
    SQINFO  sqinfo;
    int     xnum;

    if ((sqfp = SeqfileOpen(xnufile, SQFILE_UNKNOWN, NULL)) == NULL)
      Die("Failed to open sequence database file %s\n%s\n", xnufile, usage);
    while (ReadSeq(sqfp, &seq, &sqinfo)) {
      dsq = DigitizeSequence(seq, sqinfo.len);
      xnum = XNU(dsq, sqinfo.len);
      result = DedigitizeSequence(dsq, sqinfo.len);

      printf("%-20s\t%5d\n", sqinfo.name, xnum);
      if (be_verbose)
        WriteSeq(stdout, SQFILE_FASTA, result, &sqinfo);

      free(dsq);
      FreeSequence(seq, &sqinfo);
      free(result);
    }
    SeqfileClose(sqfp);
  }

  return EXIT_SUCCESS;
}
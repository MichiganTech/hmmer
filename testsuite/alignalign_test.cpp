/* alignalign_test.c
 *
 * Test driver for P7ViterbiAlignAlignment().
 *
 * The test is to
 *    1) read an alignment and a corresponding HMM
 *    2) align the alignment to the HMM to get a master trace
 *    3) map the alignment to the HMM to get another master trace
 *    4) Test that the two master traces are identical; if not, fail.
 *         This doesn't have to be true always, but it's true for the
 *         fn3 test example.
 *    5) Get imposed traces for each sequence
 *    6) Viterbi align individual seqs to the model;
 *           compare the imposed trace with the Viterbi trace;
 *    7) If an excessive number of individual traces differ from
 *          those imposed by master, fail.
 *
 */


#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "config.hpp"
#include "structs.hpp"
#include "globals.hpp"
#include "squid.hpp"
#include "msa.hpp"
#include "getopt.hpp"

static char banner[] = "\
alignalign_test : testing of P7ViterbiAlignAlignment() code";

static char usage[] = "\
Usage: alignalign_test [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -v              : be verbose\n\
";

static char experts[] = "\
  --ali <f>       : read alignment from <f>\n\
  --hmm <f>       : read HMM from <f>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       true,  sqdARG_NONE },
  { "-v",       true,  sqdARG_NONE },
  { "--ali",    false, sqdARG_STRING },
  { "--hmm",    false, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


/* Function: GCGchecksum()
 * Date:     SRE, Mon May 31 11:13:21 1999 [St. Louis]
 *
 * Purpose:  Calculate a GCG checksum for a sequence.
 *           Code provided by Steve Smith of Genetics
 *           Computer Group.
 *
 * Args:     seq -  sequence to calculate checksum for.
 *                  may contain gap symbols.
 *           len -  length of sequence (usually known,
 *                  so save a strlen() call)      
 *
 * Returns:  GCG checksum.
 */
int
GCGchecksum(char *seq, int len)
{
  int i;      /* position in sequence */
  int chk = 0;      /* calculated checksum  */

  for (i = 0; i < len; i++)
    chk = (chk + (i % 57 + 1) * (toupper((int) seq[i]))) % 10000;
  return chk;
}


/* Function: GCGMultchecksum()
 *
 * Purpose:  GCG checksum for a multiple alignment: sum of
 *           individual sequence checksums (including their
 *           gap characters) modulo 10000.
 *
 *           Implemented using spec provided by Steve Smith of
 *           Genetics Computer Group.
 *          
 * Args:     seqs - sequences to be checksummed; aligned or not
 *           nseq - number of sequences
 *          
 * Return:   the checksum, a number between 0 and 9999
 */                     
int
GCGMultchecksum(char **seqs, int nseq)
{
  int chk = 0;
  int idx;

  for (idx = 0; idx < nseq; idx++)
    chk = (chk + GCGchecksum(seqs[idx], strlen(seqs[idx]))) % 10000;
  return chk;
}


/* Function: DealignAseqs()
 *
 * Given an array of (num) aligned sequences aseqs,
 * strip the gaps. Store the raw sequences in a new allocated array.
 *
 * Caller is responsible for free'ing the memory allocated to
 * rseqs.
 *
 * Returns 1 on success. Returns 0 and sets squid_errno on
 * failure.
 */
int
DealignAseqs(char **aseqs, int num, char ***ret_rseqs)
{
  char **rseqs;                 /* de-aligned sequence array   */
  int    idx;     /* counter for sequences       */
  int    depos;     /* position counter for dealigned seq*/
  int    apos;      /* position counter for aligned seq */
  int    seqlen;    /* length of aligned seq */

        /* alloc space */
  rseqs = (char **) MallocOrDie (num * sizeof(char *));
        /* main loop */
  for (idx = 0; idx < num; idx++)
    {
      seqlen = strlen(aseqs[idx]);
        /* alloc space */
      rseqs[idx] = (char *) MallocOrDie ((seqlen + 1) * sizeof(char));

        /* strip gaps */
      depos = 0;
      for (apos = 0; aseqs[idx][apos] != '\0'; apos++)
  if (!isgap(aseqs[idx][apos]))
    {
      rseqs[idx][depos] = aseqs[idx][apos];
      depos++;
    }
      rseqs[idx][depos] = '\0';
    }
  *ret_rseqs = rseqs;
  return 1;
}


int
main(int argc, char **argv) {
  char    *hmmfile;         /* file to read HMM(s) from                */
  HMMFILE *hmmfp;               /* opened hmmfile for reading              */
  struct plan7_s  *hmm;         /* HMM to search with                      */
  char    *afile;               /* file to read alignment from             */
  int      format;              /* format determined for afile             */
  MSAFILE *afp;                 /* afile, open for reading                 */
  MSA     *msa;     /* multiple sequence alignment from afile  */
  char   **rseq;                /* raw, dealigned aseq                     */
  unsigned char     *dsq; /* digitized target sequence               */
  struct dpmatrix_s *mx;        /* reused DP alignment matrix              */
  struct p7trace_s  *mtr; /* master traceback from alignment         */
  struct p7trace_s  *maptr;     /* master traceback from mapping           */
  struct p7trace_s **tr;        /* individual tracebacks imposed by mtr    */
  struct p7trace_s **itr;       /* individual trace from P7Viterbi()       */
  int       idx;    /* counter for seqs                        */
  int       ndiff;    /* number of differing traces              */
  int       rlen;   /* length of an unaligned sequence         */

  int be_verbose;
  int be_standard;    /* true when running standard test */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

  /***********************************************
   * Parse command line
   ***********************************************/

  hmmfile     = "fn3.hmm";
  afile       = "fn3.seed";
  format      = MSAFILE_STOCKHOLM;
  be_verbose  = false;
  be_standard = true;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-v")       == 0) be_verbose = true;
    else if (strcmp(optname, "--ali")    == 0) {
      afile   = optarg;
      be_standard = false;
    } else if (strcmp(optname, "--hmm")    == 0) {
      hmmfile = optarg;
      be_standard = false;
    } else if (strcmp(optname, "-h")       == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 0)
    Die("Incorrect number of arguments.\n%s\n", usage);

  /***********************************************
   * Get one alignment from test file: must be Stockholm format.
   ***********************************************/

  if ((afp = MSAFileOpen(afile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", afile);
  if ((msa = MSAFileRead(afp)) == NULL)
    Die("Didn't read an alignment from %s", afile);
  MSAFileClose(afp);

  for (idx = 0; idx < msa->nseq; idx++)
    s2upper(msa->aseq[idx]);
  DealignAseqs(msa->aseq, msa->nseq, &rseq);

  /***********************************************
   * Open HMM file
   * Read a single HMM from it.
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("Failed to open HMM file %s\n", hmmfile);
  if (!HMMFileRead(hmmfp, &hmm))
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL)
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  HMMFileClose(hmmfp);
  P7Logoddsify(hmm, true);

  if (! (hmm->flags & PLAN7_MAP))
    Die("HMM in %s has no map", hmmfile);
  if (GCGMultchecksum(msa->aseq, msa->nseq) != hmm->checksum)
    Die("Checksum for alignment in %s does not match that in HMM (%d != %d)",
        afile, GCGMultchecksum(msa->aseq, msa->nseq), hmm->checksum);

  /***********************************************
   * First test:
   * mapped alignment should match re-aligned alignment:
   * obtain and compare the two master traces
   ***********************************************/

  mtr   = P7ViterbiAlignAlignment(msa, hmm);
  maptr = MasterTraceFromMap(hmm->map, hmm->M, msa->alen);
  if (! TraceVerify(mtr, hmm->M, msa->alen))
    Die("Trace verify on P7ViterbiAlignAlignment() result failed\n");
  if (! TraceVerify(maptr, hmm->M, msa->alen))
    Die("Trace verify on MasterTraceFromMap() result failed\n");
  if (! TraceCompare(mtr, maptr))
    Die("Master traces differ for alignment versus map\n");

  /**************************************************
   * Second test:
   * seq traces implied by mapped alignment should generally match
   * re-aligned individual sequences.
   ***************************************************/

  ImposeMasterTrace(msa->aseq, msa->nseq, mtr, &tr);

  itr = MallocOrDie(sizeof(struct p7trace_s *) * msa->nseq);

  /* Create a DP matrix; initially only two rows big, but growable;
   * we overalloc by 25 rows (L dimension) when we grow; not growable
   * in model dimension, since we know the hmm size
   */
  mx = CreatePlan7Matrix(1, hmm->M, 25, 0);

  /* align individuals, compare traces */
  ndiff = 0;
  for (idx = 0; idx < msa->nseq; idx++) {
    rlen = strlen(rseq[idx]);
    dsq  = DigitizeSequence(rseq[idx], rlen);
    P7Viterbi(dsq, rlen, hmm, mx, &(itr[idx]));
    if (! TraceCompare(itr[idx], tr[idx]))
      ndiff++;
    free(dsq);
  }
  FreePlan7Matrix(mx);

  /* Determine success/failure.
   */
  if (ndiff > msa->nseq / 2)
    Die("alignalign: Test FAILED; %d/%d differ\n", ndiff, msa->nseq);

  if (be_standard) {
    if (ndiff != 12)
      Die("alignalign: Test FAILED; %d traces differ, should be 12\n", ndiff);
    if (msa->nseq != 109)
      Die("alignalign: Test FAILED; %d seqs read, should be 109\n", msa->nseq);
  }

  if (be_verbose) printf("alignalign: Test passed; %d/%ld differ, as expected\n",
                           ndiff, msa->nseq);

  /* Cleanup.
   */
  P7FreeTrace(mtr);
  P7FreeTrace(maptr);
  for (idx = 0; idx < msa->nseq; idx++) {
    P7FreeTrace(tr[idx]);
    P7FreeTrace(itr[idx]);
  }
  free(tr);
  free(itr);
  Free2DArray((void **) rseq, msa->nseq);
  MSAFree(msa);
  FreePlan7(hmm);
  SqdClean();

  return 0;
}

/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* hmmpfam.c
 *
 * Search a single sequence against an HMM database.
 */

#include "config.h"    /* compile-time configuration constants */
//#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>

//#include "squid.h"    /* general sequence analysis library    */
#include "structs.h"    /* data structures, macros, #define's   */
#include "funcs.h"    /* function declarations                */
#include "globals.h"    /* alphabet global variables            */

static char banner[] = "hmmpfam - search one or more sequences against HMM database";

static char usage[]  = "\
Usage: hmmpfam [-options] <hmm database> <sequence file or database>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -n        : nucleic acid models/sequence (default protein)\n\
   -A <n>    : sets alignment output limit to <n> best domain alignments\n\
   -E <x>    : sets E value cutoff (globE) to <x>; default 10\n\
   -T <x>    : sets T bit threshold (globT) to <x>; no threshold by default\n\
   -Z <n>    : sets Z (# models) for E-value calculation\n\
";

static char experts[] = "\
   --acc         : use HMM accession numbers instead of names in output\n\
   --compat      : make best effort to use last version's output style\n\
   --cpu <n>     : run <n> threads in parallel (if threaded)\n\
   --cut_ga      : use Pfam GA gathering threshold cutoffs\n\
   --cut_nc      : use Pfam NC noise threshold cutoffs\n\
   --cut_tc      : use Pfam TC trusted threshold cutoffs\n\
   --domE <x>    : sets domain Eval cutoff (2nd threshold) to <x>\n\
   --domT <x>    : sets domain T bit thresh (2nd threshold) to <x>\n\
   --forward     : use the full Forward() algorithm instead of Viterbi\n\
   --informat <s>: sequence file is in format <s>\n\
   --null2       : turn OFF the post hoc second null model\n\
   --xnu         : turn ON XNU filtering of query protein sequence\n\
\n";


static struct opt_s OPTIONS[] = {
  { "-h",        true,  sqdARG_NONE },
  { "-n",        true,  sqdARG_NONE },
  { "-A",        true,  sqdARG_INT  },
  { "-E",        true,  sqdARG_FLOAT},
  { "-T",        true,  sqdARG_FLOAT},
  { "-Z",        true,  sqdARG_INT  },
  { "--acc",     false, sqdARG_NONE },
  { "--compat",  false, sqdARG_NONE },
  { "--cpu",     false, sqdARG_INT  },
  { "--cut_ga",  false, sqdARG_NONE },
  { "--cut_nc",  false, sqdARG_NONE },
  { "--cut_tc",  false, sqdARG_NONE },
  { "--domE",    false, sqdARG_FLOAT},
  { "--domT",    false, sqdARG_FLOAT},
  { "--forward", false, sqdARG_NONE },
  { "--informat",false, sqdARG_STRING},
  { "--null2",   false, sqdARG_NONE },
  { "--xnu",     false, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void main_loop_serial(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo,
                             struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
                             struct tophit_s *ghit, struct tophit_s *dhit, int *nhmm);
static void main_loop_threaded(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo,
                               struct threshold_s *thresh, int do_xnu,
                               int do_forward, int do_null2, int num_threads,
                               struct tophit_s *ghit, struct tophit_s *dhit, int *nhmm);


/* POSIX threads version:
 * the threads share a workpool_s structure amongst themselves,
 * for obtaining locks on input HMM file and output histogram and
 * tophits structures.
 */
struct workpool_s {
  /* Shared configuration resources that don't change:
   */
  char          *hmmfile;       /* name of HMM file                */
  unsigned char *dsq;           /* digitized query sequence        */
  char  *seqname;               /* sequence name                   */
  int    L;      /* length of dsq                   */
  int    do_forward;    /* true to score using Forward     */
  int    do_null2;    /* true to apply null2 correction  */
  struct threshold_s *thresh;   /* score/evalue cutoff information */

  /* Shared (mutex-protected) input resources:
   */
  HMMFILE *hmmfp;               /* ptr to open HMM file           */
  int nhmm;      /* number of HMMs searched so far */
  pthread_mutex_t input_lock;   /* mutex for locking input        */

  /* Shared (mutex-protected) output resources:
   */
  struct tophit_s *ghit;        /* per-sequence top hits */
  struct tophit_s *dhit;        /* per-domain top hits */
  pthread_mutex_t output_lock;  /* mutex for locking output */

  /* Thread pool information
   */
  pthread_t *thread;            /* our pool of threads */
  int        num_threads;       /* number of threads   */
};

static struct workpool_s *workpool_start(char *hmmfile, HMMFILE *hmmfp,
    unsigned char *dsq, char *seqname, int L,
    int do_forward, int do_null2,
    struct threshold_s *thresh,
    struct tophit_s *ghit, struct tophit_s *dhit,
    int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);





int
main(int argc, char **argv) {
  char            *hmmfile;  /* file to read HMMs from                  */
  HMMFILE         *hmmfp;       /* opened hmmfile for reading              */
  char            *seqfile;     /* file to read target sequence from       */
  SQFILE          *sqfp;        /* opened seqfile for reading              */
  int              format;  /* format of seqfile                       */
  char              *seq;  /* target sequence                         */
  SQINFO             sqinfo;  /* optional info for seq                   */
  struct fancyali_s *ali;  /* an alignment for display                */

  float   sc;      /* log-odds score in bits                  */
  double  pvalue;    /* pvalue of an HMM score                  */
  double  evalue;    /* evalue of an HMM score                  */
  double  motherp;    /* pvalue of a whole seq HMM score         */
  float   mothersc;    /* score of a whole seq parent of domain   */
  int     sqfrom, sqto;    /* coordinates in sequence                 */
  int     hmmfrom, hmmto;  /* coordinate in HMM                       */
  char   *name, *acc, *desc;    /* hit HMM name, accession, description            */
  int     hmmlen;    /* length of HMM hit                       */
  int     nhmm;      /* number of HMMs searched                 */
  int     domidx;    /* number of this domain                   */
  int     ndom;      /* total # of domains in this seq          */

  int    Alimit;    /* A parameter limiting output alignments   */
  struct threshold_s thresh;    /* contains all threshold (cutoff) info     */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   do_forward;    /* true to use Forward() not Viterbi()      */
  int   do_nucleic;    /* true to do DNA/RNA instead of protein    */
  int   do_null2;    /* true to adjust scores with null model #2 */
  int   do_xnu;      /* true to do XNU filtering                 */
  int   be_backwards;    /* true to be backwards-compatible in output*/
  int   show_acc;    /* true to sub HMM accessions for names     */
  int   i;
  int   nreported;

  int   num_threads;            /* number of worker threads */

  /***********************************************
   * Parse command line
   ***********************************************/

  format      = SQFILE_UNKNOWN;  /* default: autodetect format w/ Babelfish */
  do_forward  = false;
  do_nucleic  = false;
  do_null2    = true;
  do_xnu      = false;
  be_backwards= false;
  show_acc    = false;

  num_threads     = sysconf(_SC_NPROCESSORS_ONLN);

  Alimit         = INT_MAX;  /* no limit on alignment output     */
  thresh.globE   = 10.0;  /* use a reasonable Eval threshold; */
  thresh.globT   = -FLT_MAX;  /*   but no bit threshold,          */
  thresh.domT    = -FLT_MAX;  /*   no domain bit threshold,       */
  thresh.domE    = FLT_MAX;     /*   and no domain Eval threshold.  */
  thresh.autocut = CUT_NONE;  /*   and no Pfam cutoffs used.      */
  thresh.Z       = 0;    /* Z not preset, so determined by # of HMMs */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n")        == 0) do_nucleic     = true;
    else if (strcmp(optname, "-A")        == 0) Alimit         = atoi(optarg);
    else if (strcmp(optname, "-E")        == 0) thresh.globE   = atof(optarg);
    else if (strcmp(optname, "-T")        == 0) thresh.globT   = atof(optarg);
    else if (strcmp(optname, "-Z")        == 0) thresh.Z       = atoi(optarg);
    else if (strcmp(optname, "--acc")     == 0) show_acc       = true;
    else if (strcmp(optname, "--compat")  == 0) be_backwards   = true;
    else if (strcmp(optname, "--cpu")     == 0) num_threads    = atoi(optarg);
    else if (strcmp(optname, "--cut_ga")  == 0) thresh.autocut = CUT_GA;
    else if (strcmp(optname, "--cut_nc")  == 0) thresh.autocut = CUT_NC;
    else if (strcmp(optname, "--cut_tc")  == 0) thresh.autocut = CUT_TC;
    else if (strcmp(optname, "--domE")    == 0) thresh.domE    = atof(optarg);
    else if (strcmp(optname, "--domT")    == 0) thresh.domT    = atof(optarg);
    else if (strcmp(optname, "--forward") == 0) do_forward     = true;
    else if (strcmp(optname, "--null2")   == 0) do_null2       = false;
    else if (strcmp(optname, "--xnu")     == 0) do_xnu         = true;
    else if (strcmp(optname, "--informat") == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == SQFILE_UNKNOWN)
        Die("unrecognized sequence file format \"%s\"", optarg);
    } else if (strcmp(optname, "-h")      == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);

  hmmfile = argv[optind++];
  seqfile = argv[optind++];

  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (format == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    format = SQFILE_FASTA;

  /***********************************************
   * Open sequence database (must be in curr directory);
   * get target sequence.
   ***********************************************/

  if (do_nucleic) SetAlphabet(hmmNUCLEIC);
  else            SetAlphabet(hmmAMINO);

  if (do_nucleic && do_xnu)
    Die("You can't use -n and --xnu together: I can't xnu DNA data.");

  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence file %s\n%s\n", seqfile, usage);

  /***********************************************
   * Open HMM database (might be in HMMERDB or current directory)
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM database %s\n%s", hmmfile, usage);

  /***********************************************
   * Show the banner
   ***********************************************/

  HMMERBanner(stdout, banner);
  printf(   "HMM file:                 %s\n", hmmfile);
  printf(   "Sequence file:            %s\n", seqfile);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /***********************************************
   * Search each HMM against each sequence
   ***********************************************/

  while (ReadSeq(sqfp, format, &seq, &sqinfo)) {
    struct tophit_s   *ghit;      /* list of top hits and alignments for seq  */
    struct tophit_s   *dhit;  /* list of top hits/alignments for domains  */
    ghit = AllocTophits(20);   /* keeps full seq scores */
    dhit = AllocTophits(20);   /* keeps domain scores   */

    /* 1. Search sequence against all HMMs.
     *    Significant scores+alignments accumulate in ghit, dhit.
     */
    if (num_threads > 1)
      main_loop_threaded(hmmfile, hmmfp, seq, &sqinfo,
                         &thresh, do_xnu, do_forward, do_null2, num_threads,
                         ghit, dhit, &nhmm);
    else
      main_loop_serial(hmmfile, hmmfp, seq, &sqinfo,
                       &thresh, do_xnu, do_forward, do_null2,
                       ghit, dhit, &nhmm);

    /* set Z for good now that we're done */
    if (!thresh.Z) thresh.Z = nhmm;

    /* 2. (Done searching all HMMs for this query seq; start output)
     *    Report the overall sequence hits, sorted by significance.
     */
    if (be_backwards) {
      printf("Query:  %s  %s\n", sqinfo.name,
             (sqinfo.flags & SQINFO_DESC) ? sqinfo.desc : "");
    } else {
      printf("\nQuery sequence: %s\n", sqinfo.name);
      printf("Accession:      %s\n", (sqinfo.flags &SQINFO_ACC)  ? sqinfo.acc  : "[none]");
      printf("Description:    %s\n", (sqinfo.flags &SQINFO_DESC) ? sqinfo.desc : "[none]");
    }
    /* We'll now sort the global hit list by evalue...
     * (not score! that was bug #12. in hmmpfam, score and evalue are not
     *  monotonic.)
     */
    FullSortTophits(ghit);
    int     namewidth;    /* max width of printed HMM name           */
    int     descwidth;    /* max width of printed description        */
    namewidth = MAX(8, TophitsMaxName(ghit)); /* must print whole name, no truncation  */
    descwidth = MAX(52-namewidth, 11);        /* may truncate desc, but avoid neg len! */

    printf("\nScores for sequence family classification (score includes all domains):\n");
    printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Model", descwidth, "Description", "Score", "E-value", " N ");
    printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", descwidth, "-----------", "-----", "-------", "---");
    for (i = 0, nreported = 0; i < ghit->num; i++) {
      char *safedesc;
      GetRankedHit(ghit, i,
                   &pvalue, &sc, NULL, NULL,
                   &name, &acc, &desc,
                   NULL, NULL, NULL,           /* seq positions */
                   NULL, NULL, NULL,           /* HMM positions */
                   NULL, &ndom,                /* domain info   */
                   NULL);                     /* alignment info*/

      evalue = pvalue * (double) thresh.Z;

      /* safedesc is a workaround for an apparent Linux printf()
       * bug with the *.*s format. dbmalloc crashes with a memchr() ptr out of bounds
       * flaw if the malloc'ed space for desc is short. The workaround
             * is to make sure the ptr for *.* has a big malloc space.
       */
      if (desc != NULL && strlen(desc) < 80) {
        safedesc = MallocOrDie(sizeof(char) * 80);
        strcpy(safedesc, desc);
      } else safedesc = strdup(desc);

      /* sneaky trick warning:
       * if we're using dynamic Pfam score cutoffs (GA, TC, NC),
       * then the list of hits is already correct and does not
       * need any score cutoffs. Unset the thresholds. They'll
             * be reset in the main_loop if we still have sequences
             * to process.
       */
      if (thresh.autocut != CUT_NONE) {
        thresh.globE = thresh.domE = FLT_MAX;
        thresh.globT = thresh.domT = -FLT_MAX;
      }

      if (evalue <= thresh.globE && sc >= thresh.globT) {
        printf("%-*s %-*.*s %7.1f %10.2g %3d\n",
               namewidth,
               (show_acc && acc != NULL) ?  acc : name,
               descwidth, descwidth, safedesc != NULL ? safedesc : "",
               sc, evalue, ndom);
        nreported++;
      }
      free(safedesc);
    }
    if (nreported == 0) printf("\t[no hits above thresholds]\n");

    /* 3. Report domain hits (sorted on sqto coordinate)
     */
    FullSortTophits(dhit);
    namewidth = MAX(8, TophitsMaxName(dhit)); /* must print whole name, no truncation  */

    printf("\nParsed for domains:\n");
    printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
           namewidth, "Model", "Domain ", "seq-f", "seq-t", "hmm-f", "hmm-t", "score", "E-value");
    printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
           namewidth, "--------", "-------", "-----", "-----", "-----", "-----", "-----", "-------");

    for (i = 0, nreported = 0; i < dhit->num; i++) {
      GetRankedHit(dhit, i,
                   &pvalue, &sc, &motherp, &mothersc,
                   &name, &acc, NULL,
                   &sqfrom, &sqto, NULL,
                   &hmmfrom, &hmmto, &hmmlen,
                   &domidx, &ndom,
                   NULL);
      evalue = pvalue * (double) thresh.Z;

      /* Does the "mother" (complete) sequence satisfy global thresholds? */
      if (motherp * (double)thresh. Z > thresh.globE || mothersc < thresh.globT)
        continue;
      else if (evalue <= thresh.domE && sc >= thresh.domT) {
        printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
               namewidth,
               (show_acc && acc != NULL) ?  acc : name,
               domidx, ndom,
               sqfrom, sqto,
               sqfrom == 1 ? '[' : '.', sqto == sqinfo.len ? ']' : '.',
               hmmfrom, hmmto,
               hmmfrom == 1 ? '[':'.',  hmmto == hmmlen ? ']' : '.',
               sc, evalue);
        nreported++;
      }
    }
    if (nreported == 0) printf("\t[no hits above thresholds]\n");


    /* 3. Alignment output, also by domain.
     *    dhits is already sorted and namewidth is set, from above code.
     *    Number of displayed alignments is limited by Alimit parameter;
     *    also by domE (evalue threshold), domT (score theshold).
     */
    if (Alimit != 0) {
      printf("\nAlignments of top-scoring domains:\n");
      for (i = 0, nreported = 0; i < dhit->num; i++) {
        if (nreported == Alimit) break; /* limit to Alimit output alignments */
        GetRankedHit(dhit, i,
                     &pvalue, &sc, &motherp, &mothersc,
                     &name, &acc, NULL,
                     &sqfrom, &sqto, NULL,              /* seq position info  */
                     &hmmfrom, &hmmto, &hmmlen,         /* HMM position info  */
                     &domidx, &ndom,                    /* domain info        */
                     &ali);                        /* alignment info     */
        evalue = pvalue * (double) thresh.Z;

        if (motherp * (double) thresh.Z > thresh.globE || mothersc < thresh.globT)
          continue;
        else if (evalue <= thresh.domE && sc >= thresh.domT) {
          printf("%s: domain %d of %d, from %d to %d: score %.1f, E = %.2g\n",
                 (show_acc && acc != NULL) ?  acc : name,
                 domidx, ndom, sqfrom, sqto, sc, evalue);
          PrintFancyAli(stdout, ali);
          nreported++;
        }
      }
      if (nreported == 0)      printf("\t[no hits above thresholds]\n");
      if (nreported == Alimit) printf("\t[output cut off at A = %d top alignments]\n", Alimit);
    }


    printf("//\n");
    FreeSequence(seq, &sqinfo);
    FreeTophits(ghit);
    FreeTophits(dhit);

    HMMFileRewind(hmmfp);
  }

  /***********************************************
   * Clean-up and exit.
   ***********************************************/
  SeqfileClose(sqfp);
  HMMFileClose(hmmfp);
  SqdClean();

  return 0;
}


/* Function: main_loop_serial()
 * Purpose:  Search a sequence against an HMM database;
 *           main loop for the serial version.
 *
 *           On return, ghit and dhit contain info for all hits
 *           that satisfy the set thresholds. If an evalue
 *           cutoff is used at all, the lists will be overestimated --
 *           because the evalue will be underestimated until
 *           we know the final Z. (Thus the main program must recheck
 *           thresholds before printing any results.) If only
 *           score cutoffs are used, then the lists are correct,
 *           and may be printed exactly as they come (after
 *           appropriate sorting, anyway). This is especially
 *           important for dynamic thresholding using Pfam
 *           score cutoffs -- the main caller cannot afford to
 *           rescan the HMM file just to get the GA/TC/NC cutoffs
 *           back out for each HMM, and neither do I want to
 *           burn the space to store them as I make a pass thru
 *           Pfam.
 *
 * Args:     hmmfile - name of HMM file
 *           hmmfp   - open HMM file (and at start of file)
 *           dsq     - digitized sequence
 *           sqinfo  - ptr to SQINFO optional info for dsq
 *           thresh  - score/evalue threshold information
 *           do_xnu     - true to apply XNU filter to sequence
 *           do_forward - true to use Forward() scores
 *           do_null2   - true to adjust scores w/ ad hoc null2 model
 *           num_threads- number of threads, if threaded
 *           ghit    - global hits list
 *           dhit    - domain hits list
 *           ret_nhmm    - number of HMMs searched.
 *
 * Returns:  (void)
 */
static void
main_loop_serial(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo,
                 struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
                 struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm) {
  unsigned char     *dsq;       /* digitized sequence                      */
  int                nhmm;  /* number of HMMs searched                 */
  struct plan7_s    *hmm;       /* current HMM to search with              */
  struct p7trace_s  *tr;  /* traceback of alignment                  */
  struct dpmatrix_s *mx;        /* growable DP matrix                      */
  float   sc;                   /* an alignment score                      */

  tr = NULL;

  /* Prepare sequence.
   */
  SQD_DPRINTF1(("main_loop_serial:\n"));

  dsq = DigitizeSequence(seq, sqinfo->len);
  if (do_xnu && Alphabet_type == hmmAMINO) XNU(dsq, sqinfo->len);

  /*
   * We'll create for at least N=300xM=300, and thus consume at least 1 MB,
   * regardless of RAMLIMIT -- this helps us avoid reallocating some weird
   * asymmetric matrices.
   *
   * We're growable in both M and N, because inside of P7SmallViterbi,
   * we're going to be calling P7Viterbi on subpieces that vary in size,
   * and for different models.
   */
  mx = CreatePlan7Matrix(300, 300, 25, 25);

  nhmm = 0;
  while (HMMFileRead(hmmfp, &hmm)) {
    if (hmm == NULL)
      Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
    P7Logoddsify(hmm, !(do_forward));
    SQD_DPRINTF1(("   ... working on HMM %s\n", hmm->name));
    SQD_DPRINTF1(("   ... mx is now M=%d by N=%d\n", mx->maxM, mx->maxN));

    if (! SetAutocuts(thresh, hmm))
      Die("HMM %s did not contain the GA, TC, or NC cutoffs you needed",
          hmm->name);

    /* Score sequence, do alignment (Viterbi), recover trace
     */

#ifdef ALTIVEC

    /* By default we call an Altivec routine that doesn't save the full
     * trace. This also means the memory usage is just proportional to the
     * model (hmm->M), so we don't need to check if space is OK unless
     * we need the trace.
     */

    if(do_forward && do_null2) {
      /* Need the trace - first check space */
      if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx)) {
        /* Slower altivec version */
        sc = P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
      } else {
        /* Low-memory C version */
        sc = P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
      }
    } else {
      /* Fastest altivec version */
      sc = P7ViterbiNoTrace(dsq, sqinfo->len, hmm, mx);
      tr = NULL;
    }

#else /* not altivec */

    if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx)) {
      SQD_DPRINTF1(("   ... using normal P7Viterbi(); size ~%d MB\n",
                    P7ViterbiSize(sqinfo->len, hmm->M)));
      sc = P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
    } else {
      SQD_DPRINTF1(("   ... using P7SmallViterbi(); size ~%d MB\n",
                    P7ViterbiSize(sqinfo->len, hmm->M)));
      sc = P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
    }

#endif

    /* Implement do_forward; we'll override the whole_sc with a P7Forward()
     * calculation.
     * HMMER is so trace- (alignment-) dependent that this gets a bit hacky.
     * Some important implications:
     *   1) if --do_forward is selected, the domain (Viterbi) scores do not
     *      necessarily add up to the whole sequence (Forward) score.
     *   2) The implementation of null2 for a Forward score is undefined,
     *      since the null2 correction is trace-dependent. As a total hack,
     *      we use a null2 correction derived from the whole trace
     *      (which was the behavior of HMMER 2.1.1 and earlier, anyway).
     *      This could put the sum of domain scores and whole seq score even
     *      further in disagreement.
     *
     * Note that you can't move the Forward calculation into
     * PostprocessSignificantHit(). The Forward score will exceed the
     * Viterbi score, so you can't apply thresholds until you
     * know the Forward score. Also, since PostprocessSignificantHit()
     * is wrapped by a mutex in the threaded implementation,
     * you'd destroy all useful parallelism if PostprocessSignificantHit()
     * did anything compute intensive.
     */
    if (do_forward) {

      sc = P7Forward(dsq, sqinfo->len, hmm, NULL);
      if (do_null2) {
        /* We need the trace - recalculate it if we didn't already do it */
#ifdef ALTIVEC
        if(tr == NULL) {
          if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx)) {
            /* Slower altivec version */
            P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
          } else {
            /* Low-memory C version */
            P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
          }
        }
#endif
        sc -= TraceScoreCorrection(hmm, tr, dsq);
      }
    }
    /* Store scores/pvalue for each HMM aligned to this sequence, overall
     */
    double  pvalue;    /* pvalue of an HMM score                  */
    double  evalue;    /* evalue of an HMM score                  */
    pvalue = PValue(hmm, sc);
    evalue = thresh->Z ? (double) thresh->Z * pvalue : (double) nhmm * pvalue;
    if (sc >= thresh->globT && evalue <= thresh->globE) {
      /* Recalculate trace if we used altivec */
#ifdef ALTIVEC
      if(tr == NULL) {
        if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx)) {
          /* Slower altivec version */
          P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
        } else {
          /* Low-memory C version */
          P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
        }
      }
#endif
      sc = PostprocessSignificantHit(ghit, dhit,
                                     tr, hmm, dsq, sqinfo->len,
                                     sqinfo->name, NULL, NULL, /* won't need acc or desc even if we have 'em */
                                     do_forward, sc,
                                     do_null2,
                                     thresh,
                                     true);
      /* true -> hmmpfam mode */
    }
    if(tr != NULL)
      P7FreeTrace(tr);

    tr = NULL;

    FreePlan7(hmm);
    nhmm++;
  }

  FreePlan7Matrix(mx);
  free(dsq);
  *ret_nhmm = nhmm;
  return;
}


/*****************************************************************
 * POSIX threads implementation.
 *
 * API:
 *      workpool_start()   (makes a workpool_s structure. Starts calculations.)
 *      workpool_stop()    (waits for threads to finish.)
 *      workpool_free()    (destroys the structure)
 *
 * Threads:
 *      worker_thread()    (the actual parallelized worker thread).
 *****************************************************************/
/* Function: main_loop_threaded()
 * Purpose:  Search a sequence against an HMM database;
 *           main loop for the threaded version.
 *
 *           On return, ghit and dhit contain info for all hits
 *           that satisfy the set thresholds. If an evalue
 *           cutoff is used at all, the lists will be overestimated --
 *           because the evalue will be underestimated until
 *           we know the final Z. (Thus the main program must recheck
 *           thresholds before printing any results.) If only
 *           score cutoffs are used, then the lists are correct,
 *           and may be printed exactly as they come (after
 *           appropriate sorting, anyway). This is especially
 *           important for dynamic thresholding using Pfam
 *           score cutoffs -- the main caller cannot afford to
 *           rescan the HMM file just to get the GA/TC/NC cutoffs
 *           back out for each HMM, and neither do I want to
 *           burn the space to store them as I make a pass thru
 *           Pfam.
 *
 * Args:     hmmfile - name of HMM file
 *           hmmfp   - open HMM file (and at start of file)
 *           dsq     - digitized sequence
 *           sqinfo  - ptr to SQINFO optional info for dsq
 *           thresh  - score/evalue threshold information
 *           do_xnu     - true to apply XNU filter to sequence
 *           do_forward - true to use Forward() scores
 *           do_null2   - true to adjust scores w/ ad hoc null2 model
 *           num_threads- number of threads, >= 1
 *           ghit    - global hits list
 *           dhit    - domain hits list
 *           ret_nhmm    - number of HMMs searched.
 *
 * Returns:  (void)
 */
static void
main_loop_threaded(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo,
                   struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
                   int num_threads,
                   struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm) {
  unsigned char     *dsq;       /* digitized sequence      */
  int                nhmm;  /* number of HMMs searched */
  struct workpool_s *wpool;     /* pool of worker threads  */

  /* Prepare sequence.
   */
  dsq = DigitizeSequence(seq, sqinfo->len);
  if (do_xnu && Alphabet_type == hmmAMINO) XNU(dsq, sqinfo->len);

  wpool = workpool_start(hmmfile, hmmfp, dsq, sqinfo->name, sqinfo->len,
                         do_forward, do_null2, thresh,
                         ghit, dhit, num_threads);
  workpool_stop(wpool);
  nhmm = wpool->nhmm;
  workpool_free(wpool);

  free(dsq);
  *ret_nhmm = nhmm;
  return;
}
/* Function: workpool_start()
 * Purpose:  Initialize a workpool_s structure, and return it.
 *
 * Args:     hmmfile    - name of HMM file
 *           hmmfp      - open HMM file, at start
 *           dsq        - ptr to sequence to search
 *           seqname    - ptr to name of dsq
 *           L          - length of dsq
 *           do_forward - true to score using Forward
 *           do_null2   - true to apply null2 ad hoc correction
 *           threshold  - evalue/score threshold settings
 *           ghit       - per-seq hit list
 *           dhit       - per-domain hit list
 *           num_threads- number of worker threads to run.
 *
 * Returns:  ptr to struct workpool_s.
 *           Caller must wait for threads to finish with workpool_stop(),
 *           then free the structure with workpool_free().
 */
static struct workpool_s *
workpool_start(char *hmmfile, HMMFILE *hmmfp, unsigned char *dsq, char *seqname, int L,
               int do_forward, int do_null2, struct threshold_s *thresh,
               struct tophit_s *ghit, struct tophit_s *dhit,
               int num_threads) {
  struct workpool_s *wpool;
  pthread_attr_t    attr;
  int i;
  int rtn;

  wpool             = MallocOrDie(sizeof(struct workpool_s));
  wpool->thread     = MallocOrDie(num_threads * sizeof(pthread_t));
  wpool->hmmfile    = hmmfile;
  wpool->dsq        = dsq;
  wpool->L          = L;
  wpool->seqname    = seqname;
  wpool->do_forward = do_forward;
  wpool->do_null2   = do_null2;
  wpool->thresh     = thresh;

  wpool->hmmfp      = hmmfp;
  wpool->nhmm       = 0;
  if ((rtn = pthread_mutex_init(&(wpool->input_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->ghit       = ghit;
  wpool->dhit       = dhit;
  if ((rtn = pthread_mutex_init(&(wpool->output_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->num_threads= num_threads;

  /* Create slave threads. See comments in hmmcalibrate.c at
   * this step regarding concurrency and system scope.
   */
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  for (i = 0; i < num_threads; i++)
    if ((rtn = pthread_create(&(wpool->thread[i]), &attr,
                              worker_thread, (void *) wpool)) != 0)
      Die("Failed to create thread %d; return code %d\n", i, rtn);

  pthread_attr_destroy(&attr);
  return wpool;
}

/* Function: workpool_stop()
 * Purpose:  Waits for threads in a workpool to finish.
 *
 * Args:     wpool -- ptr to the workpool structure
 *
 * Returns:  (void)
 */
static void
workpool_stop(struct workpool_s *wpool) {
  int i;
  /* wait for threads to stop */
  for (i = 0; i < wpool->num_threads; i++)
    if (pthread_join(wpool->thread[i],NULL) != 0)
      Die("pthread_join failed");
  return;
}

/* Function: workpool_free()
 * Purpose:  Free a workpool_s structure, after the threads
 *           have finished.
 *
 * Args:     wpool -- ptr to the workpool.
 *
 * Returns:  (void)
 */
static void
workpool_free(struct workpool_s *wpool) {
  free(wpool->thread);
  free(wpool);
  return;
}


/* Function: worker_thread()
 * Purpose:  The procedure executed by the worker threads.
 *
 * Args:     ptr  - (void *) that is recast to a pointer to
 *                  the workpool.
 *
 * Returns:  (void *)
 */
void *
worker_thread(void *ptr) {
  struct workpool_s *wpool;     /* our working threads structure   */
  struct plan7_s    *hmm;       /* an HMM to search with           */
  struct p7trace_s  *tr;        /* traceback from an alignment     */
  float  sc;      /* score of an alignment           */
  struct threshold_s thresh;  /* a local copy of thresholds      */
  struct dpmatrix_s *mx;        /* growable DP matrix              */

  tr = NULL;

  wpool = (struct workpool_s *) ptr;
  /* Because we might dynamically change the thresholds using
   * Pfam GA/NC/TC cutoffs, we make a local copy of the threshold
   * structure in this thread.
   */
  thresh.globT   = wpool->thresh->globT;
  thresh.globE   = wpool->thresh->globE;
  thresh.domT    = wpool->thresh->domT;
  thresh.domE    = wpool->thresh->domE;
  thresh.autocut = wpool->thresh->autocut;
  thresh.Z       = wpool->thresh->Z;

  /*
   * We'll create for at least N=300xM=300, and thus consume at least 1 MB,
   * regardless of RAMLIMIT -- this helps us avoid reallocating some weird
   * asymmetric matrices.
   *
   * We're growable in both M and N, because inside of P7SmallViterbi,
   * we're going to be calling P7Viterbi on subpieces that vary in size,
   * and for different models.
   */
  mx = CreatePlan7Matrix(300, 300, 25, 25);

  for (;;) {

    /* 1. acquire lock on HMM input, and get
     *    the next HMM to work on.
     */
    /* acquire a lock */
    int    rtn;      /* a return code from pthreads lib */
    if ((rtn = pthread_mutex_lock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));

    if (! HMMFileRead(wpool->hmmfp, &hmm)) {
      /* we're done. release lock, exit thread */
      if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
        Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
      FreePlan7Matrix(mx);
      pthread_exit(NULL);
    }
    wpool->nhmm++;
    SQD_DPRINTF1(("a thread is working on %s\n", hmm->name));
    /* release the lock */
    if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    if (hmm == NULL)
      Die("HMM file %s may be corrupt or in incorrect format; parse failed", wpool->hmmfile);
    P7Logoddsify(hmm, !(wpool->do_forward));

    if (!SetAutocuts(&thresh, hmm))
      Die("HMM %s did not have the right GA, NC, or TC cutoffs", hmm->name);

    /* 2. We have an HMM in score form.
     *    Score the sequence.
     */

#ifdef ALTIVEC

    /* By default we call an Altivec routine that doesn't save the full
     * trace. This also means the memory usage is just proportional to the
     * model (hmm->M), so we don't need to check if space is OK unless
     * we need the trace.
     */

    if(wpool->do_forward && wpool->do_null2) {
      /* Need the trace - first check space */
      if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx)) {
        /* Slower altivec version */
        sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
      } else {
        /* Low-memory C version */
        sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
      }
    } else {
      /* Fastest altivec version */
      sc = P7ViterbiNoTrace(wpool->dsq, wpool->L, hmm, mx);
      tr = NULL;
    }

#else /* not altivec */

    if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx)) {
      sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
    } else {
      sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
    }

#endif

    /* The Forward score override (see comments in serial vers)
     */
    if (wpool->do_forward) {
      sc  = P7Forward(wpool->dsq, wpool->L, hmm, NULL);
      if (wpool->do_null2) {
#ifdef ALTIVEC
        if(tr == NULL) {
          if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx)) {
            /* Slower altivec version */
            sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
          } else {
            /* Low-memory C version */
            sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
          }
        }
#endif
        sc -= TraceScoreCorrection(hmm, tr, wpool->dsq);
      }
    }

    /* 3. Save the output in tophits structures, after acquiring a lock
     */
    if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    SQD_DPRINTF1(("model %s scores %f\n", hmm->name, sc));

    double pvalue;    /* P-value of score                */
    double evalue;    /* E-value of a score              */
    pvalue = PValue(hmm, sc);
    evalue = thresh.Z ? (double) thresh.Z * pvalue : (double) wpool->nhmm * pvalue;
    if (sc >= thresh.globT && evalue <= thresh.globE) {
#ifdef ALTIVEC
      if(tr == NULL) {
        if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx)) {
          /* Slower altivec version */
          sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
        } else {
          /* Low-memory C version */
          sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
        }
      }
#endif
      sc = PostprocessSignificantHit(wpool->ghit, wpool->dhit,
                                     tr, hmm, wpool->dsq, wpool->L,
                                     wpool->seqname,
                                     NULL, NULL, /* won't need seq's acc or desc */
                                     wpool->do_forward, sc,
                                     wpool->do_null2,
                                     &thresh,
                                     true);
      /* true -> hmmpfam mode */
    }
    if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    if(tr != NULL)
      P7FreeTrace(tr);

    tr = NULL;

    FreePlan7(hmm);

  } /* end 'infinite' loop over HMMs in this thread */
}

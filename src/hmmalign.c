/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* hmmalign.c
 *
 * main() for aligning a set of sequences to an HMM.
 */

#include "hmmalign.h"    /* compile-time configuration constants */
#include "selex.h"


char banner[] = "hmmalign - align sequences to an HMM profile";

char usage[]  = "\
Usage: hmmalign [-options] <hmm file> <sequence file>\n\
Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -m     : only print symbols aligned to match states\n\
   -o <f> : save alignment in file <f>\n\
   -q     : quiet - suppress verbose banner\n\
";

char experts[] = "\
   --informat <s>  : sequence file is in format <s>\n\
   --mapali <f>    : include alignment in file <f> using map in HMM\n\
   --oneline       : output Stockholm fmt with 1 line/seq, not interleaved\n\
   --outformat <s> : output alignment in format <s> [default: Stockholm]\n\
                       formats include: MSF, Clustal, Phylip, SELEX\n\
   --withali <f>   : include alignment to (fixed) alignment in file <f>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h", true, sqdARG_NONE   },
  { "-m", true, sqdARG_NONE   },
  { "-o", true, sqdARG_STRING },
  { "-q", true, sqdARG_NONE   },
  { "--informat",  false, sqdARG_STRING },
  { "--mapali",    false, sqdARG_STRING },
  { "--oneline",   false, sqdARG_NONE },
  { "--outformat", false, sqdARG_STRING },
  { "--withali",   false, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))




int
DealignedLength(
  char *aseq
){
  int rlen; 
  for (rlen = 0; *aseq; aseq++)
    if (! isgap(*aseq)) rlen++;
  return rlen;
}


/* Function: GCGchecksum()
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


/* Function: MSAToSqinfo()
 * Purpose:  Take an MSA and generate a SQINFO array suitable
 *           for use in annotating the unaligned sequences.
 *           Return the array.
 *           
 *           Permanent temporary code. sqinfo was poorly designed.
 *           it must eventually be replaced, but the odds
 *           of this happening soon are nil, so I have to deal.
 *
 * Args:     msa   - the alignment
 *
 * Returns:  ptr to allocated sqinfo array.
 *           Freeing is ghastly: free in each individual sqinfo[i] 
 *           with FreeSequence(NULL, &(sqinfo[i])), then
 *           free(sqinfo).
 */
SQINFO *
MSAToSqinfo(
  MSA *msa
){
  SQINFO *sqinfo;

  sqinfo = MallocOrDie(sizeof(SQINFO) * msa->nseq);

  for (size_t idx = 0; idx < msa->nseq; idx++){
    sqinfo[idx].flags = 0;
    SetSeqinfoString(&(sqinfo[idx]), msa->sqname[idx], SQINFO_NAME);
    SetSeqinfoString(&(sqinfo[idx]), MSAGetSeqAccession(msa, idx), SQINFO_ACC);
    SetSeqinfoString(&(sqinfo[idx]), MSAGetSeqDescription(msa, idx), SQINFO_DESC);

    if (msa->ss != NULL && msa->ss[idx] != NULL) {
      MakeDealignedString(msa->aseq[idx], msa->alen, msa->ss[idx], &(sqinfo[idx].ss));
      sqinfo[idx].flags |= SQINFO_SS;
    }

    if (msa->sa != NULL && msa->sa[idx] != NULL) {
      MakeDealignedString(msa->aseq[idx], msa->alen, msa->sa[idx], &(sqinfo[idx].sa));
      sqinfo[idx].flags |= SQINFO_SA;
    }

    sqinfo[idx].len    = DealignedLength(msa->aseq[idx]);
    sqinfo[idx].flags |= SQINFO_LEN;
  }
  return sqinfo;
}


void
SeqinfoCopy(
  SQINFO *sq1, 
  SQINFO *sq2
){
  sq1->flags = sq2->flags;
  if (sq2->flags & SQINFO_NAME)  strcpy(sq1->name, sq2->name);
  if (sq2->flags & SQINFO_ID)    strcpy(sq1->id,   sq2->id);
  if (sq2->flags & SQINFO_ACC)   strcpy(sq1->acc,  sq2->acc);
  if (sq2->flags & SQINFO_DESC)  strcpy(sq1->desc, sq2->desc);
  if (sq2->flags & SQINFO_LEN)   sq1->len    = sq2->len;
  if (sq2->flags & SQINFO_START) sq1->start  = sq2->start;
  if (sq2->flags & SQINFO_STOP)  sq1->stop   = sq2->stop;
  if (sq2->flags & SQINFO_OLEN)  sq1->olen   = sq2->olen;
  if (sq2->flags & SQINFO_TYPE)  sq1->type   = sq2->type;
  if (sq2->flags & SQINFO_SS)    sq1->ss     = strdup(sq2->ss);
  if (sq2->flags & SQINFO_SA)    sq1->sa     = strdup(sq2->sa);
}


int
main(int argc, char **argv) {
  char            *hmmfile;  /* file to read HMMs from                  */
  HMMFILE         *hmmfp;       /* opened hmmfile for reading              */
  struct plan7_s  *hmm;         /* HMM to align to                         */
  char            *seqfile;     /* file to read target sequence from       */
  int              format;  /* format of seqfile                       */
  char           **rseq;        /* raw, unaligned sequences                */
  SQINFO          *sqinfo;      /* info associated with sequences          */
  unsigned char  **dsq;         /* digitized raw sequences                 */
  int              nseq;        /* number of sequences                     */
  float           *wgt;    /* weights to assign to alignment          */
  MSA             *msa;         /* alignment that's created                */
  int              i;
  struct dpmatrix_s *mx;        /* growable DP matrix                      */
  struct p7trace_s **tr;        /* traces for aligned sequences            */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   be_quiet;    /* true to suppress verbose banner          */
  int   matchonly;    /* true to show only match state syms       */
  char *outfile;                /* optional alignment output file           */
  int   outfmt;      /* code for output alignment format         */
  int   do_oneline;             /* true to do one line/seq, no interleaving */
  FILE *ofp;                    /* handle on alignment output file          */
  char *withali;                /* name of additional alignment file to align */
  char *mapali;                 /* name of additional alignment file to map   */

  /***********************************************
   * Parse command line
   ***********************************************/

  format     = SQFILE_UNKNOWN;    /* default: autodetect input format     */
  outfmt     = MSAFILE_STOCKHOLM; /* default: output in Stockholm format  */
  do_oneline = false;      /* default: interleaved format          */
  matchonly  = false;
  outfile    = NULL;      /* default: output alignment to stdout  */
  be_quiet   = false;
  withali    = NULL;
  mapali     = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-m")        == 0) matchonly  = true;
    else if (strcmp(optname, "-o")        == 0) outfile    = optarg;
    else if (strcmp(optname, "-q")        == 0) be_quiet   = true;
    else if (strcmp(optname, "--mapali")  == 0) mapali     = optarg;
    else if (strcmp(optname, "--oneline") == 0) do_oneline = true;
    else if (strcmp(optname, "--withali") == 0) withali    = optarg;
    else if (strcmp(optname, "--informat") == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == SQFILE_UNKNOWN)
        Die("unrecognized sequence file format \"%s\"", optarg);
    } else if (strcmp(optname, "--outformat") == 0) {
      outfmt = String2SeqfileFormat(optarg);
      if (outfmt == MSAFILE_UNKNOWN)
        Die("unrecognized output alignment file format \"%s\"", optarg);
      if (! IsAlignmentFormat(outfmt))
        Die("\"%s\" is not a multiple alignment format", optarg);
    } else if (strcmp(optname, "-h") == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
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
   * Open HMM file (might be in HMMERDB or current directory).
   * Read a single HMM from it.
   *
   * Currently hmmalign disallows the J state and
   * only allows one domain per sequence. To preserve
   * the S/W entry information, the J state is explicitly
   * disallowed, rather than calling a Plan7*Config() function.
   * this is a workaround in 2.1 for the 2.0.x "yo!" bug.
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm))
    Die("Failed to read any HMMs from %s\n", hmmfile);
  HMMFileClose(hmmfp);
  if (hmm == NULL)
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  hmm->xt[XTE][MOVE] = 1.;        /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  P7Logoddsify(hmm, true);
  /* do we have the map we might need? */
  if (mapali != NULL && ! (hmm->flags & PLAN7_MAP))
    Die("HMMER: HMM file %s has no map; you can't use --mapali.", hmmfile);

  /***********************************************
   * Open sequence file in current directory.
   * Read all seqs from it.
   ***********************************************/

  if (! ReadMultipleRseqs(seqfile, format, &rseq, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);

  /***********************************************
   * Show the banner
   ***********************************************/

  if (! be_quiet) {
    HMMERBanner(stdout, banner);
    printf(   "HMM file:             %s\n", hmmfile);
    printf(   "Sequence file:        %s\n", seqfile);
    printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  }

  /***********************************************
   * Do the work
   ***********************************************/

  /* Allocations and initializations.
   */
  dsq = MallocOrDie(sizeof(unsigned char *)    * nseq);
  tr  = MallocOrDie(sizeof(struct p7trace_s *) * nseq);
  mx  = CreatePlan7Matrix(1, hmm->M, 25, 0);

  /* Align each sequence to the model, collect traces
   */
  for (i = 0; i < nseq; i++) {
    dsq[i] = DigitizeSequence(rseq[i], sqinfo[i].len);

    if (P7ViterbiSpaceOK(sqinfo[i].len, hmm->M, mx))
      (void) P7Viterbi(dsq[i], sqinfo[i].len, hmm, mx, &(tr[i]));
    else
      (void) P7SmallViterbi(dsq[i], sqinfo[i].len, hmm, mx, &(tr[i]));
  }
  FreePlan7Matrix(mx);

  /* Include an aligned alignment, if desired.
   */
  if (mapali != NULL)
    include_alignment(mapali, hmm, true, &rseq, &dsq, &sqinfo, &tr, &nseq);
  if (withali != NULL)
    include_alignment(withali, hmm, false, &rseq, &dsq, &sqinfo, &tr, &nseq);

  /* Turn traces into a multiple alignment
   */
  wgt = MallocOrDie(sizeof(float) * nseq);
  FSet(wgt, nseq, 1.0);
  msa = P7Traces2Alignment(dsq, sqinfo, wgt, nseq, hmm->M, tr, matchonly);

  /***********************************************
   * Output the alignment
   ***********************************************/

  if (outfile != NULL && (ofp = fopen(outfile, "w")) != NULL) {
    MSAFileWrite(msa, outfmt, do_oneline);
    printf("Alignment saved in file %s\n", outfile);
    fclose(ofp);
  } else
    MSAFileWrite(msa, outfmt, do_oneline);


  /***********************************************
   * Cleanup and exit
   ***********************************************/

  for (i = 0; i < nseq; i++) {
    P7FreeTrace(tr[i]);
    FreeSequence(rseq[i], &(sqinfo[i]));
    free(dsq[i]);
  }
  MSAFree(msa);
  FreePlan7(hmm);
  free(sqinfo);
  free(rseq);
  free(dsq);
  free(wgt);
  free(tr);

  SqdClean();
  return 0;
}


/* Function: include_alignment()
 * Purpose:  Given the name of a multiple alignment file,
 *           align that alignment to the HMM, and add traces
 *           to an existing array of traces. If do_mapped
 *           is true, we use the HMM's map file. If not,
 *           we use P7ViterbiAlignAlignment().
 *
 * Args:     seqfile  - name of alignment file
 *           hmm      - model to align to
 *           do_mapped- true if we're to use the HMM's alignment map
 *           rsq      - RETURN: array of rseqs to add to
 *           dsq      - RETURN: array of dsq to add to
 *           sqinfo   - RETURN: array of SQINFO to add to
 *           tr       - RETURN: array of traces to add to
 *           nseq     - RETURN: number of seqs
 *
 * Returns:  new, realloc'ed arrays for rsq, dsq, sqinfo, tr; nseq is
 *           increased to nseq+ainfo.nseq.
 */
void
include_alignment(
  char *seqfile, 
  struct plan7_s *hmm, 
  int do_mapped,
  char ***rsq, 
  unsigned char ***dsq, 
  SQINFO **sqinfo,
  struct p7trace_s ***tr, 
  int *nseq
){
  int format;      /* format of alignment file */
  MSA   *msa;      /* alignment to align to    */
  MSAFILE *afp;
  SQINFO  *newinfo;    /* sqinfo array from msa */
  unsigned char **newdsq;
  char **newrseq;
  // idx      /* counter over aseqs       */
  struct p7trace_s *master;     /* master trace             */
  struct p7trace_s **addtr;     /* individual traces for aseq */

  format = MSAFILE_UNKNOWN;  /* invoke Babelfish */
  if ((afp = MSAFileOpen(seqfile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", seqfile);
  if ((msa = MSAFileRead(afp)) == NULL)
    Die("Failed to read an alignment from %s\n", seqfile);
  MSAFileClose(afp);
  for (size_t idx = 0; idx < msa->nseq; idx++)
    s2upper(msa->aseq[idx]);
  newinfo = MSAToSqinfo(msa);

  /* Verify checksums before mapping */
  if (do_mapped && GCGMultchecksum(msa->aseq, msa->nseq) != hmm->checksum)
    Die("The checksums for alignment file %s and the HMM alignment map don't match.",
        seqfile);
  /* Get a master trace */
  if (do_mapped) master = MasterTraceFromMap(hmm->map, hmm->M, msa->alen);
  else           master = P7ViterbiAlignAlignment(msa, hmm);

  /* convert to individual traces */
  ImposeMasterTrace(msa->aseq, msa->nseq, master, &addtr);
  /* add those traces to existing ones */
  *tr = MergeTraceArrays(*tr, *nseq, addtr, msa->nseq);

  /* additional bookkeeping: add to dsq, sqinfo */
  *rsq = ReallocOrDie((*rsq), sizeof(char *) * (*nseq + msa->nseq));
  DealignAseqs(msa->aseq, msa->nseq, &newrseq);
  for (size_t idx = *nseq; idx < *nseq + msa->nseq; idx++)
    (*rsq)[idx] = newrseq[idx - (*nseq)];
  free(newrseq);

  *dsq = ReallocOrDie((*dsq), sizeof(unsigned char *) * (*nseq + msa->nseq));
  DigitizeAlignment(msa, &newdsq);
  for (size_t idx = *nseq; idx < *nseq + msa->nseq; idx++)
    (*dsq)[idx] = newdsq[idx - (*nseq)];
  free(newdsq);
  /* unnecessarily complex, but I can't be bothered... */
  *sqinfo = ReallocOrDie((*sqinfo), sizeof(SQINFO) * (*nseq + msa->nseq));
  for (size_t idx = *nseq; idx < *nseq + msa->nseq; idx++)
    SeqinfoCopy(&((*sqinfo)[idx]), &(newinfo[idx - (*nseq)]));
  free(newinfo);

  *nseq = *nseq + msa->nseq;

  /* Cleanup */
  P7FreeTrace(master);
  MSAFree(msa);
  /* Return */
  return;
}

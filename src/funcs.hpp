/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* funcs.h
 *
 * Declarations of external functions in HMMER.
 */

#pragma once

#include "config.h"
#include "structs.h"
#include "squid.h"
#include "msa.h"

/* alphabet.c
 * Configuration of global alphabet information
 */
void
DetermineAlphabet(
  char **rseqs, 
  int  nseq);


void
SetAlphabet(
  int type);


unsigned char  
SymbolIndex(
  char sym);


unsigned char* 
DigitizeSequence(
  char *seq, 
  int L);


char*          
DedigitizeSequence(
  unsigned char *dsq, 
  int L);


void           
DigitizeAlignment(
  MSA *msa, 
  unsigned char ***ret_dsqs);


void           
P7CountSymbol(
  float *counters, 
  unsigned char sym, 
  float wt);


void           
DefaultGeneticCode(
  int *aacode);


void           
DefaultCodonBias(
  float *codebias);


/* from core_algorithms.c
 * Clean research/demonstration versions of basic algorithms.
 */


/* Function: CreatePlan7Matrix()
 *
 * Purpose:  Create a dynamic programming matrix for standard Forward,
 *           Backward, or Viterbi, with scores kept as scaled log-odds
 *           integers. Keeps 2D arrays compact in RAM in an attempt
 *           to maximize cache hits.
 *
 *           The mx structure can be dynamically grown, if a new
 *           HMM or seq exceeds the currently allocated size. Dynamic
 *           growing is more efficient than an alloc/free of a whole
 *           matrix for every new target. The ResizePlan7Matrix()
 *           call does this reallocation, if needed. Here, in the
 *           creation step, we set up some pads - to inform the resizing
 *           call how much to overallocate when it realloc's. If a pad
 *           is zero, we will not resize in that dimension.
 *
 * Args:     N     - N+1 rows are allocated, for sequence.
 *           M     - size of model in nodes
 *           padN  - over-realloc in seq/row dimension, or 0
 *           padM  - over-realloc in HMM/column dimension, or 0
 *
 * Return:   mx
 *           mx is allocated here. Caller frees with FreePlan7Matrix(mx).
 */
struct dpmatrix_s*
CreatePlan7Matrix(
  int N, 
  int M, 
  int padN, 
  int padM);


void   
ResizePlan7Matrix(
  struct dpmatrix_s *mx, 
  int N, 
  int M,
  int ***xmx, 
  int ***mmx, 
  int ***imx, 
  int ***dmx);


struct dpmatrix_s*
AllocPlan7Matrix(
  int rows, 
  int M,
  int ***xmx, 
  int ***mmx, 
  int ***imx, 
  int ***dmx);


struct dpshadow_s*
AllocShadowMatrix(
  int rows, 
  int M, 
  char ***xtb,
  char ***mtb,
  char ***itb, 
  char ***dtb);


void  
FreePlan7Matrix(
  struct dpmatrix_s *mx);


void  
FreeShadowMatrix(
  struct dpshadow_s *tb);


int   
P7ViterbiSpaceOK(
  int L, 
  int M, 
  struct dpmatrix_s *mx);


int   
P7ViterbiSize(
  int L, 
  int M);


int   
P7SmallViterbiSize(
  int L, 
  int M);


int   
P7WeeViterbiSize(
  int L, 
  int M);


float 
P7Forward(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s **ret_mx);


float 
P7Viterbi(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm, 
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);


float 
P7ViterbiNoTrace(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s *mx);


void  
P7ViterbiTrace(
  struct plan7_s *hmm, 
  unsigned char *dsq, 
  int L,
  struct dpmatrix_s *mx, 
  struct p7trace_s **ret_tr);


float 
P7SmallViterbi(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s *mx, 
  struct p7trace_s **ret_tr);


float 
P7ParsingViterbi(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr);


float 
P7WeeViterbi(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr);


float 
Plan7ESTViterbi(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s **ret_mx);


struct p7trace_s*
P7ViterbiAlignAlignment(
  MSA *msa, 
  struct plan7_s *hmm);


struct p7trace_s*
ShadowTrace(
  struct dpshadow_s *tb, 
  struct plan7_s *hmm, 
  int L);


float  
PostprocessSignificantHit(
  struct tophit_s *ghit, 
  struct tophit_s *dhit,
  struct p7trace_s   *tr, 
  struct plan7_s *hmm, 
  unsigned char *dsq,
  int L, 
  char *seqname, 
  char *seqacc, 
  char *seqdesc,
  int do_forward, 
  float sc_override, 
  int do_null2,
  struct threshold_s *thresh, 
  int hmmpfam_mode);


/* from debug.c
 * Debugging output of various sorts.
 */
char*
Statetype(
  char st);


char*
AlphabetType2String(
  int type);


void 
P7PrintTrace(
  FILE *fp, 
  struct p7trace_s *tr,
  struct plan7_s *hmm, 
  unsigned char *dsq);


void 
P7PrintPrior(
  FILE *fp, 
  struct p7prior_s *pri);


int  
TraceCompare(
  struct p7trace_s *t1, 
  struct p7trace_s *t2);


int  
TraceVerify(
  struct p7trace_s *tr, 
  int M, 
  int N);


/*
 * from display.c
 * Ian Holmes' functions for displaying HMMER2 data structures, especially
 * for posterior probabilities in alignments.
 */
void 
DisplayPlan7Matrix(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s *mx);


void 
DisplayPlan7Posteriors(
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s *forward, 
  struct dpmatrix_s *backward,
  struct p7trace_s *viterbi, 
  struct p7trace_s *optacc);


void 
DisplayPlan7PostAlign(
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s *forward, 
  struct dpmatrix_s *backward,
  struct p7trace_s **alignment, 
  int A);


/* from emit.c
 * Generation of sequences/traces from an HMM
 */
void 
EmitSequence(
  struct plan7_s *hmm, 
  unsigned char **ret_dsq, 
  int *ret_L,
  struct p7trace_s **ret_tr);


void 
EmitConsensusSequence(
  struct plan7_s *hmm, 
  char **ret_seq, 
  unsigned char **ret_dsq,
  int *ret_L, 
  struct p7trace_s **ret_tr);


void 
StateOccupancy(
  struct plan7_s *hmm, 
  float **ret_mp, 
  float **ret_ip, 
  float **ret_dp);


/* from emulation.c
 * Interfaces between HMMER and other software packages
 */
void 
WriteProfile(
  FILE *fp, 
  struct plan7_s *hmm, 
  int do_xsw);


/* from evolution.c
 * Phylogenetic extrapolation of profile HMMs.
 */
int 
EvolveOneTransitionVector(
  float *qs, 
  float ts, 
  int n, 
  float *q0, 
  float *qz,
  float t, 
  float *q);


/* from histogram.c
 * accumulation of scores
 */
struct histogram_s*
AllocHistogram(
  int min, 
  int max, 
  int lumpsize);


void 
FreeHistogram(
  struct histogram_s *h);


void 
UnfitHistogram(
  struct histogram_s *h);


void 
AddToHistogram(
  struct histogram_s *h, 
  float sc);


void 
PrintASCIIHistogram(
  FILE *fp, 
  struct histogram_s *h);


void 
PrintXMGRHistogram(
  FILE *fp, 
  struct histogram_s *h);


void 
PrintXMGRDistribution(
  FILE *fp, 
  struct histogram_s *h);


void 
PrintXMGRRegressionLine(
  FILE *fp, 
  struct histogram_s *h);


void 
EVDBasicFit(
  struct histogram_s *h);


int  
ExtremeValueFitHistogram(
  struct histogram_s *h, 
  int censor,
  float high_hint);


void 
ExtremeValueSetHistogram(
  struct histogram_s *h, 
  float mu, 
  float lambda,
  float low, 
  float high, 
  int ndegrees);


int  
GaussianFitHistogram(
  struct histogram_s *h);


void 
GaussianSetHistogram(
  struct histogram_s *h, 
  float mean, 
  float sd);


double 
EVDDensity(
  float x, 
  float mu, 
  float lambda);


double 
EVDDistribution(
  float x, 
  float mu, 
  float lambda);


double 
ExtremeValueP (
  float x, 
  float mu, 
  float lambda);


double 
ExtremeValueP2(
  float x, 
  float mu, 
  float lambda, 
  int N);


double 
ExtremeValueE (
  float x, 
  float mu, 
  float lambda, 
  int N);


float  
EVDrandom(
  float mu, 
  float lambda);


int    
EVDMaxLikelyFit(
  float *x, 
  int *y, 
  int n,
  float *ret_mu, 
  float *ret_lambda);


int    
EVDCensoredFit(
  float *x, 
  int *y, 
  int n, 
  int z, 
  float c,
  float *ret_mu, float *ret_lambda);


void   
Lawless416(
  float *x, 
  int *y, 
  int n, 
  float lambda,
  float *ret_f, 
  float *ret_df);


void   
Lawless422(
  float *x, 
  int *y, 
  int n, 
  int z, 
  float c,
  float lambda, 
  float *ret_f, 
  float *ret_df);


/* from hmmio.c
 * Input/output (saving/reading) of models
 */
HMMFILE*
HMMFileOpen(
  char *hmmfile, 
  char *env);


int      
HMMFileRead(
  HMMFILE *hmmfp, 
  struct plan7_s **ret_hmm);


void     
HMMFileClose(
  HMMFILE *hmmfp);


int      
HMMFileFormat(
  HMMFILE *hmmfp);


void     
HMMFileRewind(
  HMMFILE *hmmfp);


int      
HMMFilePositionByName(
  HMMFILE *hmmfp, 
  char *name);


int      
HMMFilePositionByIndex(
  HMMFILE *hmmfp, 
  int idx);


void     
WriteAscHMM(
  FILE *fp, 
  struct plan7_s *hmm);


void     
WriteBinHMM(
  FILE *fp, 
  struct plan7_s *hmm);


/* from infocontent.c
 * Evolving to specified information content
 */
void  
AdjustAveInfoContent (
  struct plan7_s *hmm, 
  float desired,
  char *matrixfile);


void  
EvolveEmits (
  double *temp_emits, 
  double *P, 
  int L);


float 
CalculateBackgroundEntropy ();


float 
CalculateEmitsEntropy (
  double *emits, 
  int L);


void  
NormalizeEmits (
  double *temp_emits, 
  int L);


void  
PrintAveInfoContent (
  struct plan7_s *hmm);


/* masks.c
 * Repetitive sequence masking.
 */
int   
XNU(
  unsigned char *dsq, 
  int len);


float 
TraceScoreCorrection(
  struct plan7_s *hmm, 
  struct p7trace_s *tr, 
  unsigned char *dsq);


/* mathsupport.c
 * Much of this code deals with Dirichlet prior mathematics.
 */
int   
Prob2Score(
  float p, 
  float null);


float 
Score2Prob(
  int sc, 
  float null);


float 
Scorify(
  int sc);


double 
PValue(
  struct plan7_s *hmm, 
  float sc);


float 
LogSum(
  float p1, 
  float p2);


int   
ILogsum(
  int p1, 
  int p2);


void  
LogNorm(
  float *vec, 
  int n);


float 
Logp_cvec(
  float *cvec, 
  int n, 
  float *alpha);


float 
P_PvecGivenDirichlet(
  float *p, 
  int n, 
  float *alpha);


/* from matrices.c
 * for matrix manipulation
 */
void   
ReadAAMatrices(
  double **ret_Sij, 
  double **ret_pi,
  char *matrixfile, 
  int environ, 
  int L);


void  
AssignWagMatrix(
  double **Sij, 
  double **pi);


void   
ReadMatrices (
  double **ret_Sij, 
  double **ret_pi,
  char *matrixfile, 
  int environ, 
  int L);


void   
PrintMatrices(
  double *prob, 
  int L, 
  int environ);


void   
UnlogAndPrintMatrices(
  double *prob, 
  int L, 
  int environ);


void   
SymToRateMatrices(
  double *Qij, 
  double *Sij, 
  double *pi, 
  int L,
  int environ);


void   
NormRateMatrices(
  double *Qij, 
  double *pi, 
  int L, 
  int environ);


void    
AssignMatrixNotLog (
  double *Qij, 
  int matrix_n, 
  double time,
  double *Pij);


double*
Cal_Id(
  int L);


double*
Cal_M_Exp(
  double *M, 
  int L, 
  double power);


void    
Comp_M_Exp(
  double *M, 
  int L, 
  double power);


void    
Comp_M_N_Prod(
  double *M, 
  double *N, 
  int L);


void    
CheckSingleProb(
  double *psingle, 
  int size);


void    
CopyMatrix (
  double *copyQ, 
  double *Q, 
  int N);


int     
Check_Accuracy(
  double *vec, 
  int L);


void  
LogifyMatrix(
  double *M, 
  int L);


/* from misc.c
 * Miscellaneous functions with no home
 */
void  
HMMERBanner(
  FILE *fp, 
  char *banner);


char*
Getword(
  FILE *fp, 
  int type);


char*
Getline(
  char *s, 
  int n, 
  FILE *fp);


int   
SetAutocuts(
  struct threshold_s *thresh, 
  struct plan7_s *hmm);


/* from modelmakers.c
 * Model construction algorithms
 */
void 
P7Handmodelmaker(
  MSA *msa, 
  unsigned char **dsq, 
  char *isfrag,
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr);


void 
P7Fastmodelmaker(
  MSA *msa, 
  unsigned char **dsq, 
  char *isfrag,
  float symfrac, 
  struct plan7_s **ret_hmm,
  struct p7trace_s ***ret_tr);


/* from plan7.c
 * Plan7 HMM structure support
 */
struct 
plan7_s*
AllocPlan7(
  int M);


struct 
plan7_s*
AllocPlan7Shell();


void 
AllocPlan7Body(
  struct plan7_s *hmm, 
  int M);


void 
FreePlan7(
  struct plan7_s *hmm);


void 
ZeroPlan7(
  struct plan7_s *hmm);


void 
Plan7SetName(
  struct plan7_s *hmm, 
  char *name);


void 
Plan7SetAccession(
  struct plan7_s *hmm, 
  char *acc);


void 
Plan7SetDescription(
  struct plan7_s *hmm, 
  char *desc);


void 
Plan7ComlogAppend(
  struct plan7_s *hmm, 
  int argc, 
  char **argv);


void 
Plan7SetCtime(
  struct plan7_s *hmm);


void 
Plan7SetNullModel(
  struct plan7_s *hmm, 
  float null[MAXABET], 
  float p1);


void 
P7Logoddsify(
  struct plan7_s *hmm, 
  int viterbi_mode);


void 
Plan7Rescale(
  struct plan7_s *hmm, 
  float scale);


void 
Plan7Renormalize(
  struct plan7_s *hmm);


void 
Plan7RenormalizeExits(
  struct plan7_s *hmm);


void 
Plan7NakedConfig(
  struct plan7_s *hmm);


void 
Plan7GlobalConfig(
  struct plan7_s *hmm);


void 
Plan7LSConfig(
  struct plan7_s *hmm);


void 
Plan7SWConfig(
  struct plan7_s *hmm, 
  float pentry, 
  float pexit);


void 
Plan7FSConfig(
  struct plan7_s *hmm, 
  float pentry, 
  float pexit);


void 
PrintPlan7Stats(
  FILE *fp, 
  struct plan7_s *hmm, 
  unsigned char **dsq,
  int nseq, 
  struct p7trace_s **tr);


int  
DegenerateSymbolScore(
  float *p, 
  float *null, 
  int ambig);


void 
Plan9toPlan7(
  struct plan9_s *hmm, 
  struct plan7_s **ret_plan7);


/*
 * from plan9.c
 * Backwards compatibility for the Plan 9 data structures of HMMER 1.x
 */
struct plan9_s*
P9AllocHMM(
  int M);


void 
P9ZeroHMM(
  struct plan9_s *hmm);


int  
P9FreeHMM(
  struct plan9_s *hmm);


void 
P9Renormalize(
  struct plan9_s *hmm);


void 
P9DefaultNullModel(
  float *null);


/*
 * from postprob.c
 * Functions for working with posterior probabilities within alignments
 */
float 
P7OptimalAccuracy(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm, 
  struct p7trace_s **ret_tr);


float 
P7Backward(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,  
  struct dpmatrix_s **ret_mx);


void  
P7EmitterPosterior(
  int L, 
  struct plan7_s *hmm, 
  struct dpmatrix_s *forward,
  struct dpmatrix_s *backward, 
  struct dpmatrix_s *mx);


float 
P7FillOptimalAccuracy(
  int L, 
  int M, 
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx, 
  struct p7trace_s **ret_tr);


void  
P7OptimalAccuracyTrace(
  int L, 
  int M, 
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx, 
  struct p7trace_s **ret_tr);


char*
PostalCode(
  int L, 
  struct dpmatrix_s *mx, 
  struct p7trace_s *tr);


/* from prior.c
 * Dirichlet priors
 */
struct p7prior_s*
P7AllocPrior();


struct p7prior_s*
P7LaplacePrior();


struct p7prior_s*
P7DefaultPrior();


struct p7prior_s*
P7ReadPrior(
  char *prifile);


void 
P7FreePrior(
  struct p7prior_s *pri);


void 
PAMPrior(
  char *pamfile, 
  struct p7prior_s *pri, 
  float pamwgt);


void 
P7DefaultNullModel(
  float *null, 
  float *ret_p1);


void 
P7ReadNullModel(
  char *rndfile, 
  float *null, 
  float *ret_p1);


void 
P7PriorifyHMM(
  struct plan7_s *hmm, 
  struct p7prior_s *pri);


void 
P7PriorifyTransitionVector(
  float *t, 
  struct p7prior_s *prior,
  float tq[MAXDCHLET]);


void 
P7PriorifyEmissionVector(
  float *vec, 
  struct p7prior_s *pri,
  int num, 
  float eq[MAXDCHLET],
  float e[MAXDCHLET][MAXABET],
  float *ret_mix);


/* from threads.c
 * POSIX threads implementation
 */
int   
ThreadNumber();


/* from tophits.c
 * Support for keeping/sorting top scoring hit/alignment lists
 */
struct tophit_s*
AllocTophits(
  int lumpsize);


void   
GrowTophits(
  struct tophit_s *h);


void   
FreeTophits(
  struct tophit_s *h);


struct fancyali_s*
AllocFancyAli();


void   
FreeFancyAli(
  struct fancyali_s *ali);


void
RegisterHit(
  struct tophit_s *h, 
  double sortkey,
  double pvalue, 
  float score,
  double motherp, 
  float mothersc,
  char *name, 
  char *acc, 
  char *desc,
  int sqfrom, 
  int sqto, 
  int sqlen,
  int hmmfrom, 
  int hmmto, 
  int hmmlen,
  int domidx, 
  int ndom,
  struct fancyali_s *ali);


void 
GetRankedHit(
  struct tophit_s *h, 
  int rank,
  double *r_pvalue, 
  float *r_score,
  double *r_motherp, 
  float *r_mothersc,
  char **r_name, 
  char **r_acc, 
  char **r_desc,
  int *r_sqfrom, 
  int *r_sqto, 
  int *r_sqlen,
  int *r_hmmfrom, 
  int *r_hmmto, 
  int *r_hmmlen,
  int *r_domidx, 
  int *r_ndom,
  struct fancyali_s **r_ali);


int    
TophitsMaxName(
  struct tophit_s *h);


void   
FullSortTophits(
  struct tophit_s *h);


void   
TophitsReport(
  struct tophit_s *h, 
  double E, 
  int nseq);


/* from trace.c
 * Support for traceback (state path) structure
 */
void  
P7AllocTrace(
  int tlen, 
  struct p7trace_s **ret_tr);


void  
P7ReallocTrace(
  struct p7trace_s *tr, 
  int tlen);


void  
P7FreeTrace(
  struct p7trace_s *tr);


void  
TraceSet(
  struct p7trace_s *tr, 
  int tpos, 
  char type, 
  int idx, 
  int pos);


struct p7trace_s**
MergeTraceArrays(
  struct p7trace_s **t1, 
  int n1, 
  struct p7trace_s **t2, 
  int n2);


void  
P7ReverseTrace(
  struct p7trace_s *tr);


void  
P7TraceCount(
  struct plan7_s *hmm, 
  unsigned char *dsq, 
  float wt,
  struct p7trace_s *tr);


float 
P7TraceScore(
  struct plan7_s *hmm, 
  unsigned char *dsq, 
  struct p7trace_s *tr);


MSA*
P7Traces2Alignment(
  unsigned char **dsq, 
  SQINFO *sqinfo, 
  float *wgt,
  int nseq, 
  int M,
  struct p7trace_s **tr, 
  int matchonly);


int  
TransitionScoreLookup(
  struct plan7_s *hmm, 
  char st1,
  int k1,
  char st2, 
  int k2);


struct fancyali_s*
CreateFancyAli(
  struct p7trace_s *tr, 
  struct plan7_s *hmm,
  unsigned char *dsq, 
  char *name);


void 
PrintFancyAli(
  FILE *fp, 
  struct fancyali_s *ali);


void 
TraceDecompose(
  struct p7trace_s *otr, 
  struct p7trace_s ***ret_tr,
  int *ret_ntr);


int  
TraceDomainNumber(
  struct p7trace_s *tr);


void 
TraceSimpleBounds(
  struct p7trace_s *tr, 
  int *ret_i1, 
  int *ret_i2,
  int *ret_k1,  
  int *ret_k2);


struct p7trace_s*
MasterTraceFromMap(
  int *map, 
  int M, 
  int alen);


void 
ImposeMasterTrace(
  char **aseq, 
  int nseq, 
  struct p7trace_s *mtr,
  struct p7trace_s ***ret_tr);

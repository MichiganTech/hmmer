#pragma once

/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* core_algorithms.c
 *
 * Simple and robust "research" implementations of Forward, Backward,
 * and Viterbi for Plan7. For optimized replacements for some of these functions,
 * see fast_algorithms.c.
 */


#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include "config.hpp"
#include "structs.hpp"
#include "vectorops.hpp"


float get_wee_midpt(
  struct plan7_s *hmm,
  unsigned char *dsq,
  int k1,
  char t1,
  int s1,
  int k3,
  char t3,
  int s3,
  int *ret_k2,
  char *ret_t2,
  int *ret_s2);


/* Function: ResizePlan7Matrix()
 *
 * Purpose:  Reallocate a dynamic programming matrix, if necessary,
 *           for a problem of NxM: sequence length N, model size M.
 *           (N=1 for small memory score-only variants; we allocate
 *           N+1 rows in the DP matrix.)
 *
 *           We know (because of the way hmmsearch and hmmpfam are coded)
 *           that only one of the two dimensions is going to change
 *           in size after the first call to ResizePlan7Matrix();
 *           that is, for hmmsearch, we have one HMM of fixed size M
 *           and our target sequences may grow in N; for hmmpfam,
 *           we have one sequence of fixed size N and our target models
 *           may grow in M. What we have to watch out for is P7SmallViterbi()
 *           working on a divide and conquer problem and passing us N < maxN,
 *           M > maxM; we should definitely *not* reallocate a smaller N.
 *           Since we know that only one dimension is going to grow,
 *           we aren't scared of reallocating to maxN,maxM. (If both
 *           M and N could grow, we would be more worried.)
 *
 *           Returns individual ptrs to the four matrix components
 *           as a convenience.
 *
 * Args:     mx    - an already allocated model to grow.
 *           N     - seq length to allocate for; N+1 rows
 *           M     - size of model
 *           xmx, mmx, imx, dmx
 *                 - RETURN: ptrs to four mx components as a convenience
 *
 * Return:   (void)
 *           mx is (re)allocated here.
 */
void
ResizePlan7Matrix(
  struct dpmatrix_s *mx,
  int N,
  int M,
  int ***xmx,
  int ***mmx,
  int ***imx,
  int ***dmx);


/* Function: AllocPlan7Matrix()
 * Purpose:  Used to be the main allocator for dp matrices; we used to
 *           allocate, calculate, free. But this spent a lot of time
 *           in malloc(). Replaced with Create..() and Resize..() to
 *           allow matrix reuse in P7Viterbi(), the main alignment
 *           engine. But matrices are alloc'ed by other alignment engines
 *           too, ones that are less frequently called and less
 *           important to optimization of cpu performance. Instead of
 *           tracking changes through them, for now, provide
 *           an Alloc...() call with the same API that's just a wrapper.
 *
 * Args:     rows  - generally L+1, or 2; # of DP rows in seq dimension to alloc
 *           M     - size of model, in nodes
 *           xmx, mmx, imx, dmx
 *                 - RETURN: ptrs to four mx components as a convenience
 *
 * Returns:  mx
 *           Caller free's w/ FreePlan7Matrix()
 */
struct dpmatrix_s*
AllocPlan7Matrix(
  int rows,
  int M,
  int ***xmx,
  int ***mmx,
  int ***imx,
  int ***dmx);


/* Function: FreePlan7Matrix()
 *
 * Purpose:  Free a dynamic programming matrix allocated by CreatePlan7Matrix().
 *
 * Return:   (void)
 */
void
FreePlan7Matrix(
  struct dpmatrix_s *mx);


/* Function: AllocShadowMatrix()
 *
 * Purpose:  Allocate a dynamic programming traceback pointer matrix for
 *           a Viterbi algorithm.
 *
 * Args:     rows  - number of rows to allocate; typically L+1
 *           M     - size of model
 *           xtb, mtb, itb, dtb
 *                 - RETURN: ptrs to four mx components as a convenience
 *
 * Return:   mx
 *           mx is allocated here. Caller frees with FreeDPMatrix(mx).
 */

struct dpshadow_s*
AllocShadowMatrix(
  int rows,
  int M,
  char ***xtb,
  char ***mtb,
  char ***itb,
  char ***dtb);


/* Function: FreeShadowMatrix()
 *
 * Purpose:  Free a dynamic programming matrix allocated by AllocShadowMatrix().
 *
 * Return:   (void)
 */
void
FreeShadowMatrix(
  struct dpshadow_s *tb);


/* Function:  P7ViterbiSpaceOK()
 *
 * Purpose:   Returns true if the existing matrix allocation
 *            is already sufficient to hold the requested MxN, or if
 *            the matrix can be expanded in M and/or N without
 *            exceeding RAMLIMIT megabytes.
 *
 *            This gets called anytime we're about to do P7Viterbi().
 *            If it returns false, we switch into the appropriate
 *            small-memory alternative: P7SmallViterbi() or
 *            P7WeeViterbi().
 *
 *            Checking the DP problem size against P7ViterbiSize()
 *            is not enough, because the DP matrix may be already
 *            allocated in MxN. For example, if we're already
 *            allocated to maxM,maxN of 1447x979, and we try to
 *            do a problem of MxN=12x125000, P7ViterbiSize() may
 *            think that's fine - but when we resize, we can only
 *            grow, so we'll expand to 1447x125000, which is
 *            likely over the RAMLIMIT. [bug#h26; xref SLT7 p.122]
 *
 * Args:
 *
 * Returns:   true if we can run P7Viterbi(); false if we need
 *            to use a small memory variant.
 *
 * Xref:      STL7 p.122.
 */
int
P7ViterbiSpaceOK(
  int L,
  int M,
  struct dpmatrix_s *mx);


/* Function: P7ViterbiSize()
 *
 * Purpose:  Returns the ballpark predicted memory requirement for a
 *           P7Viterbi() alignment, in MB.
 *
 *           Currently L must fit in an int (< 2 GB), but we have
 *           to deal with LM > 2 GB - e.g. watch out for overflow, do
 *           the whole calculation in floating point. Bug here detected
 *           in 2.1.1 by David Harper, Sanger Centre.
 *
 * Args:     L  - length of sequence
 *           M  - length of HMM
 *
 * Returns:  # of MB
 */
int
P7ViterbiSize(
  int L,
  int M);


/* Function: P7Forward()
 *
 * Purpose:  The Forward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
P7Forward(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s **ret_mx);



/* Function: P7ViterbiTrace()
 *
 * Purpose:  Traceback of a Viterbi matrix: i.e. retrieval
 *           of optimum alignment.
 *
 * Args:     hmm    - hmm, log odds form, used to make mx
 *           dsq    - sequence aligned to (digital form) 1..N
 *           N      - length of seq
 *           mx     - the matrix to trace back in, N x hmm->M
 *           ret_tr - RETURN: traceback.
 *
 * Return:   (void)
 *           ret_tr is allocated here. Free using P7FreeTrace().
 */
void
P7ViterbiTrace(
  struct plan7_s *hmm,
  unsigned char *dsq,
  int N,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);


/* Function: P7SmallViterbi()
 *
 * Purpose:  Wrapper function, for linear memory alignment
 *           with same arguments as P7Viterbi().
 *
 *           Calls P7ParsingViterbi to break the sequence
 *           into fragments. Then, based on size of fragments,
 *           calls either P7Viterbi() or P7WeeViterbi() to
 *           get traces for them. Finally, assembles all these
 *           traces together to produce an overall optimal
 *           trace for the sequence.
 *
 *           If the trace isn't needed for some reason,
 *           all we do is call P7ParsingViterbi.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - DP matrix, growable
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *
 * Returns:  Score of optimal alignment in bits.
 */
float
P7SmallViterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);


/* Function: P7ParsingViterbi()
 *
 * Purpose:  The "hmmfs" linear-memory algorithm for finding
 *           the optimal alignment of a very long sequence to
 *           a looping, multihit (e.g. Plan7) model, parsing it into
 *           a series of nonoverlapping subsequences that match
 *           the model once. Other algorithms (e.g. P7Viterbi()
 *           or P7WeeViterbi()) are applied subsequently to
 *           these subsequences to recover complete alignments.
 *
 *           The hmmfs algorithm appears briefly in [Durbin98],
 *           but is otherwise unpublished.
 *
 *           The traceback structure returned is special: a
 *           "collapsed" trace S->B->E->...->B->E->T, where
 *           stateidx is unused and pos is used to indicate the
 *           position of B and E in the sequence. The matched
 *           subsequence is B_pos+1...E_pos. The number of
 *           matches in the trace is (tlen/2)-1.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model (log odds scores ready)
 *           ret_tr - RETURN: a collapsed traceback.
 *
 * Returns:  Score of the optimal Viterbi alignment, in bits.
 */
float
P7ParsingViterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr);


/* Function: P7WeeViterbi()
 * Purpose:  Hirschberg/Myers/Miller linear memory alignment.
 *           See [Hirschberg75,MyM-88a] for the idea of the algorithm.
 *           Adapted to HMM implementation.
 *
 *           Requires that you /know/ that there's only
 *           one hit to the model in the sequence: either
 *           because you're forcing single-hit, or you've
 *           previously called P7ParsingViterbi to parse
 *           the sequence into single-hit segments. The reason
 *           for this is that a cyclic model (a la Plan7)
 *           defeats the nice divide and conquer trick.
 *           (I think some trickery with propagated trace pointers
 *           could get around this but haven't explored it.)
 *           This is implemented by ignoring transitions
 *           to/from J state.
 *
 *           Will not work (in current implementation) for sequences
 *           of length 1. (xref bug #h30).
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq; must be L > 1!
 *           hmm    - the model
 *           ret_tr - RETURN: traceback.
 *
 * Returns:  Score of the optimal Viterbi alignment.
 */
float
P7WeeViterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr);


/* Function: get_wee_midpt()
 * Purpose:  The heart of the divide and conquer algorithm
 *           for P7WeeViterbi(). This function is called
 *           recursively to find successive optimal midpoints
 *           in the alignment matrix. See P7WeeViterbi() for
 *           further comments on the assumptions of this algorithm.
 *
 * Args:     hmm   - the model, set up for integer scores
 *           dsq   - the sequence, digitized
 *           L     - length of the sequence
 *           k1    - model node to start with, 1..M
 *           t1    - state type to start with, STM | STI | STN | STC; STS to start
 *           s1    - sequence position to start with, 1..L; 1 to start
 *           k3    - model node to end with, 1..M
 *           t3    - state type to end with, STM | STI | STN | STC; STT to start
 *           s3    - sequence position to end with, 1..L; L to start
 *          ret_k2 - RETURN: optimal midpoint, node position in model
 *          ret_t2 - RETURN: optimal midpoint, state type
 *          ret_s2 - RETURN: optimal midpoint, sequence position
 *
 * Returns: score of optimal alignment, in bits.
 */
// static float
// get_wee_midpt(struct plan7_s *hmm, unsigned char *dsq, int L,
//               int k1, char t1, int s1,
//               int k3, char t3, int s3,
//               int *ret_k2, char *ret_t2, int *ret_s2) {
float
get_wee_midpt(
  struct plan7_s *hmm,
  unsigned char *dsq,
  int k1,
  char t1,
  int s1,
  int k3,
  char t3,
  int s3,
  int *ret_k2,
  char *ret_t2,
  int *ret_s2);


/* Function: P7ViterbiAlignAlignment()
 *
 * Purpose:  Align a multiple alignment to an HMM without
 *           changing the multiple alignment itself.
 *           Adapted from P7Viterbi().
 *
 *           Heuristic; not a guaranteed optimal alignment.
 *           Guaranteeing an optimal alignment appears difficult.
 *           [cryptic note to myself:] In paths connecting to I* metastates,
 *           recursion breaks down; if there is a gap in the
 *           previous column for a given seq, we can't determine what state the
 *           I* metastate corresponds to for this sequence, unless we
 *           look back in the DP matrix. The lookback would either involve
 *           recursing back to the previous M* metastate (giving a
 *           O(MN^2) algorithm instead of O(MN)) or expanding the I*
 *           metastate into 3^nseq separate I* metastates to keep track
 *           of which of three states each seq is in. Since the second
 *           option blows up exponentially w/ nseq, it is not attractive.
 *           If the first option were used, the correct algorithm would be related to
 *           modelmakers.c:Maxmodelmaker(), but somewhat more difficult.
 *
 *           The heuristic approach here is to calculate a "consensus"
 *           sequence from the alignment, and align the consensus to the HMM.
 *           Some hackery is employed, weighting transitions and emissions
 *           to make things work (re: con and mocc arrays).
 *
 * Args:     aseq  - aligned sequences
 *           ainfo - info for aseqs (includes alen, nseq, wgt)
 *           hmm   - model to align to
 *
 * Returns:  Traceback. Caller must free with P7FreeTrace().
 *           pos[] contains alignment columns, indexed 1..alen.
 *           statetype[] contains metastates M*, etc. as STM, etc.
 */
struct p7trace_s *
P7ViterbiAlignAlignment(
  MSA *msa,
  struct plan7_s *hmm);


/* Function: ShadowTrace()
 * Purpose:  Given a shadow matrix, trace it back, and return
 *           the trace.
 *
 * Args:     tb  - shadow matrix of traceback pointers
 *           hmm - the model (needed for figuring out wing unfolding)
 *           L   - sequence length
 *
 * Returns:  traceback. Caller must free w/ P7FreeTrace().
 */
struct p7trace_s *
ShadowTrace(
  struct dpshadow_s *tb,
  struct plan7_s *hmm,
  int L);



/* Function: PostprocessSignificantHit()
 * Purpose:  Add a significant hit to per-seq and per-domain hit
 *           lists, after postprocessing the scores appropriately,
 *           and making sure per-domain scores add up to the per-seq
 *           score.
 *
 *          [doesn't really belong in core_algorithms.c, because
 *           it's more of a hack than an algorithm, but on the other
 *           hand it's now part of the core of how HMMER scores
 *           things. Maybe there should be a core_hacks.c.]
 *
 *           Given: active hit lists for per-seq and per-domain
 *           scores (e.g. hmmpfam and hmmsearch, collating their
 *           results), and a new hit that's significant enough
 *           that it may need to be reported in final output.
 *
 *           Breaks the traceback into individual domain traces;
 *           scores each one of them, then applies null2 correction
 *           for biased composition. Recalculates the per-seq score
 *           as the sum of the per-domain scores. Stores the hits
 *           in the lists, for eventual sorting and output by the
 *           caller.
 *
 * Notes:    In principle we've got the score, and a pvalue, and a traceback
 *           by doing the Viterbi algorithm, right? What else is left
 *           to do? Well, in practice, life is more complicated, because
 *           of the trace-dependent null2 score correction.
 *
 *           After a null2 score correction is carried out on
 *           each domain (the default) the number of detected domains
 *           with scores > 0 may have decreased. We want the
 *           global (per-seq) hit list to have the recalculated number of
 *           domains, not necessarily what Viterbi gave us.
 *
 *           Also, since we want the global score to be the sum of
 *           the individual domains, but the null2 correction is
 *           applied to each domain individually, we have to calculate
 *           an adjusted global score. (To do otherwise invites
 *           subtle inconsistencies; xref bug 2.)
 *
 *           We don't have final evalues, so we may put a few
 *           more hits into the hit lists than we end up reporting.
 *           The main output routine is responsible for final
 *           enforcement of the thresholds.
 *
 *           This routine is NOT THREADSAFE. When multithreaded,
 *           with using shared ghit/dhit output buffers, calls to
 *           PostprocessSignificantHit() need to be protected.
 *
 * Args:     ghit     - an active list of per-seq (global) hits
 *           dhit     - an active list of per-domain hits
 *           tr       - the significant HMM/seq traceback we'll report on
 *           hmm      - ptr to the HMM
 *           dsq      - digitized sequence (1..L)
 *           L        - length of dsq
 *           seqname  - name of sequence (same as targname, in hmmsearch)
 *           seqacc   - seq's accession (or NULL)
 *           seqdesc  - seq's description (or NULL)
 *           do_forward  - true if we've already calculated final per-seq score
 *           sc_override - per-seq score to use if do_forward is true
 *           do_null2    - true to apply the null2 scoring correction
 *           thresh      - contains the threshold/cutoff information.
 *           hmmpfam_mode - true if called by hmmpfam, else assumes hmmsearch;
 *                          affects how the lists' sort keys are set.
 *
 * Returns:  the recalculated per-seq score (or sc_override),
 *           as appropriate, for subsequent storage in the histogram
 */
float
PostprocessSignificantHit(
  struct tophit_s *ghit,
  struct tophit_s *dhit,
  struct p7trace_s *tr,
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


/* fast_algorithms.c
 *
 * Optimized routines to replace slower implementations in core_algorithms.c.
 *
 * The routines in core_algorithms.c are designed for clarity
 * and maintainability, not for speed. Implementations here
 * are designed for speed, not clarity. If you're trying to
 * understand the code, or optimize for a specific platform,
 * you are probably better off looking at core_algorithms.c.
 *
 * P7Viterbi() is the key function to target optimization to.
 * The implementation in core_algorithms.c is currently ifdef'ed
 * out of the code. The implementation that is used by default
 * is here, in fast_algorithms.c. A third implementation, from
 * Erik Lindahl at Stanford, is Mac/Altivec specific.
 *
 * Which implementation is used is controlled by ifdef's. The
 * default implementation uses a fast implementation of
 * P7Viterbi() from here. Other options (mutually exclusive):
 *
 * -DSLOW
 *   enable original core_algorithms.c code: slower than default,
 *   but might be easier to follow, for someone trying
 *   to understand the DP code.
 * -DALTIVEC
 *   enable Erik Lindahl's Altivec code for Macintosh OSX
 */



/* the DEFAULT P7Viterbi() is portably optimized; code follows:
 */
/* Function: P7Viterbi() - portably optimized version
 *
 * Purpose:  The Viterbi dynamic programming algorithm.
 *           Derived from core_algorithms.c:P7Viterbi().
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - re-used DP matrix
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
P7Viterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);



/* Echo a message about Altivec being used once:
 *   (1) Confirm that the Altivec kernel is indeed being used
 *   (2) Trace bugs to the Altivec kernel instead of the vanilla code
 */

void
AltivecMessage();


/*################################################################
 * The Altivec port, for Macintosh PowerPC.
 * Erik Lindahl, Stanford University, 2002.
 *
 * Replaces the following functions:
 *    AllocPlan7Body()      plan7.c                (data alignment on 16-byte boundaries)
 *    CreatePlan7Matrix()   core_algorithms.c      (data alignment on 16-byte boundaries)
 *    ResizePlan7Matrix()   core_algorithms.c      (data alignment on 16-byte boundaries)
 *    P7Viterbi()           core_algorithms.c      (total recode, w/ Altivec instructions)
 ################################################################*/
void
AllocPlan7Body(
  struct plan7_s *hmm,
  int M);


struct dpmatrix_s*
CreatePlan7Matrix(
  int N,
  int M,
  int padN,
  int padM
);


float
P7Viterbi(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);


/* Function: P7ViterbiNoTrace()
 *
 * Purpose:  The Viterbi dynamic programming algorithm, but a version
 *           that does not store the dynamic programming matrix, and thus
 *           doesn't return the trace - just the score.
 *           This boost the Altivec performance significantly, since
 *           the tuned code otherwise would be memory bandwidth limited.
 *           Altivec implementation by Erik Lindahl, Stanford University, 2004.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - DP matrix (may get grown here)
 *
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
P7ViterbiNoTrace(
  unsigned char *dsq,
  int L,
  struct plan7_s *hmm,
  struct dpmatrix_s *mx);

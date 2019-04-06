/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/


/* plan7.c
 *
 * Support for Plan 7 HMM data structure, plan7_s.
 */


#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "squid.hpp"
#include "structs.hpp"
#include "vectorops.hpp"


/* Functions: AllocPlan7(), AllocPlan7Shell(), AllocPlan7Body(), FreePlan7()
 *
 * Purpose:   Allocate or free a Plan7 HMM structure.
 *            Can either allocate all at one (AllocPlan7()) or
 *            in two steps (AllocPlan7Shell(), AllocPlan7Body()).
 *            The two step method is used in hmmio.c where we start
 *            parsing the header of an HMM file but don't
 *            see the size of the model 'til partway thru the header.
 */
struct plan7_s *
AllocPlan7(int M);


struct plan7_s *
AllocPlan7Shell();


void
FreePlan7(struct plan7_s *hmm);


/* Function: ZeroPlan7()
 *
 * Purpose:  Zeros the counts/probabilities fields in a model.
 *           Leaves null model untouched.
 */
void
ZeroPlan7(
  struct plan7_s *hmm);


/* Function: Plan7SetName()
 *
 * Purpose:  Change the name of a Plan7 HMM. Convenience function.
 *
 * Note:     Trailing whitespace and \n's are chopped.
 */
void
Plan7SetName(
  struct plan7_s *hmm,
  char *name);


/* Function: Plan7SetAccession()
 *
 * Purpose:  Change the accession number of a Plan7 HMM. Convenience function.
 *
 * Note:     Trailing whitespace and \n's are chopped.
 */
void
Plan7SetAccession(
  struct plan7_s *hmm,
  char *acc);


/* Function: Plan7SetDescription()
 *
 * Purpose:  Change the description line of a Plan7 HMM. Convenience function.
 *
 * Note:     Trailing whitespace and \n's are chopped.
 */
void
Plan7SetDescription(
  struct plan7_s *hmm,
  char *desc);


/* Function: Plan7ComlogAppend()
 *
 * Purpose:  Concatenate command line options and append to the
 *           command line log.
 */
void
Plan7ComlogAppend(
  struct plan7_s *hmm,
  int argc,
  char **argv);


/* Function: Plan7SetCtime()
 *
 * Purpose:  Set the ctime field in a new HMM to the current time.
 */
void
Plan7SetCtime(
  struct plan7_s *hmm);


/* Function: Plan7SetNullModel()
 *
 * Purpose:  Set the null model section of an HMM.
 *           Convenience function.
 */
void
Plan7SetNullModel(
  struct plan7_s *hmm,
  float null[MAXABET],
  float p1);


/* Function: P7Logoddsify()
 *
 * Purpose:  Take an HMM with valid probabilities, and
 *           fill in the integer log-odds score section of the model.
 *
 *    Notes on log-odds scores:
 *         type of parameter       probability        score
 *         -----------------       -----------        ------
 *         any emission             p_x           log_2 p_x/null_x
 *             N,J,C /assume/ p_x = null_x so /always/ score zero.
 *         transition to emitters   t_x           log_2 t_x/p1
 *            (M,I; N,C; J)
 *             NN and CC loops are often equal to p1, so usu. score zero.
 *         C->T transition          t_x            log_2 t_x/p2
 *            often zero, usu. C->T = p2.
 *         all other transitions    t_x           log_2 t_x
 *             (no null model counterpart, so null prob is 1)
 *
 *    Notes on entry/exit scores, B->M and M->E:
 *         The probability form model includes delete states 1 and M.
 *         these states are removed from a search form model to
 *         prevent B->D...D->E->J->B mute cycles, which would complicate
 *         dynamic programming algorithms. The data-independent
 *         S/W B->M and M->E transitions are folded together with
 *         data-dependent B->D...D->M and M->D...D->E paths.
 *
 *         This process is referred to in the code as "wing folding"
 *         or "wing retraction"... the analogy is to a swept-wing
 *         fighter in landing vs. high speed flight configuration.
 *
 *    Note on Viterbi vs. forward flag:
 *         Wing retraction must take forward vs. Viterbi
 *         into account. If forward, sum two paths; if Viterbi, take
 *         max. I tried to slide this by as a sum, without
 *         the flag, but Alex detected it as a bug, because you can
 *         then find cases where the Viterbi score doesn't match
 *         the P7TraceScore().
 *
 * Args:      hmm          - the hmm to calculate scores in.
 *            viterbi_mode - true to fold wings in Viterbi configuration.
 *
 * Return:    (void)
 *            hmm scores are filled in.
 */
void
P7Logoddsify(
  struct plan7_s *hmm,
  int viterbi_mode);


/* Function:  Plan7Rescale()
 *
 * Purpose:   Scale a counts-based HMM by some factor, for
 *            adjusting counts to a new effective sequence number.
 *
 * Args:      hmm        - counts based HMM.
 *            scale      - scaling factor (e.g. eff_nseq/nseq); 1.0= no scaling.
 *
 * Returns:   (void)
 */
void
Plan7Rescale(
  struct plan7_s *hmm,
  float scale);



/* Function: Plan7Renormalize()
 *
 * Purpose:  Take an HMM in counts form, and renormalize
 *           all of its probability vectors. Also enforces
 *           Plan7 restrictions on nonexistent transitions.
 *
 * Args:     hmm - the model to renormalize.
 *
 * Return:   (void)
 *           hmm is changed.
 */
void
Plan7Renormalize(
  struct plan7_s *hmm);


/* Function: Plan7RenormalizeExits()
 *
 * Purpose:  Renormalize just the match state transitions;
 *           for instance, after a Config() function has
 *           modified the exit distribution.
 *
 * Args:     hmm - hmm to renormalize
 *
 * Returns:  void
 */
void
Plan7RenormalizeExits(
  struct plan7_s *hmm);


/*****************************************************************
 * Plan7 configuration functions
 * The following few functions are the Plan7 equivalent of choosing
 * different alignment styles (fully local, fully global, global/local,
 * multihit, etc.)
 *
 * There is (at least) one constraint worth noting.
 * If you want per-domain scores to sum up to per-sequence scores,
 * then one of the following two sets of conditions must be met:
 *
 *   1) t(E->J) = 0
 *      e.g. no multidomain hits
 *
 *   2) t(N->N) = t(C->C) = t(J->J) = hmm->p1
 *      e.g. unmatching sequence scores zero, and
 *      N->B first-model score is equal to J->B another-model score.
 *
 * These constraints are obeyed in the default Config() functions below,
 * but in the future (when HMM editing may be allowed) we'll have
 * to remember this. Non-equality of the summed domain scores and
 * the total sequence score is a really easy "red flag" for people to
 * notice and report as a bug, even if it may make probabilistic
 * sense not to meet either constraint for certain modeling problems.
 *****************************************************************
 */

/* Function: Plan7NakedConfig()
 *
 * Purpose:  Set the alignment-independent, algorithm-dependent parameters
 *           of a Plan7 model so that no special states (N,C,J) emit anything:
 *           one simple, full global pass through the model.
 *
 * Args:     hmm - the plan7 model
 *
 * Return:   (void)
 *           The HMM is modified; algorithm dependent parameters are set.
 *           Previous scores are invalidated if they existed.
 */
void
Plan7NakedConfig(
  struct plan7_s *hmm);


/* Function: Plan7GlobalConfig()
 *
 * Purpose:  Set the alignment-independent, algorithm-dependent parameters
 *           of a Plan7 model to global (Needleman/Wunsch) configuration.
 *
 *           Like a non-looping hmmls, since we actually allow flanking
 *           N and C terminal sequence.
 *
 * Args:     hmm - the plan7 model
 *
 * Return:   (void)
 *           The HMM is modified; algorithm dependent parameters are set.
 *           Previous scores are invalidated if they existed.
 */
void
Plan7GlobalConfig(
  struct plan7_s *hmm);


/* Function: Plan7LSConfig()
 *
 * Purpose:  Set the alignment independent parameters of a Plan7 model
 *           to hmmls (global in HMM, local in sequence) configuration.
 *
 * Args:     hmm  - the plan7 model
 *
 * Return:   (void);
 *           the HMM probabilities are modified.
 */
void
Plan7LSConfig(
  struct plan7_s *hmm);


/* Function: Plan7SWConfig()
 *
 * Purpose:  Set the alignment independent parameters of
 *           a Plan7 model to hmmsw (Smith/Waterman) configuration.
 *
 * Notes:    entry/exit is asymmetric because of the left/right
 *           nature of the HMM/profile. Entry probability is distributed
 *           simply by assigning p_x = pentry / (M-1) to M-1
 *           internal match states. However, the same approach doesn't
 *           lead to a flat distribution over exit points. Exit p's
 *           must be corrected for the probability of a previous exit
 *           from the model. Requiring a flat distribution over exit
 *           points leads to an easily solved piece of algebra, giving:
 *                      p_1 = pexit / (M-1)
 *                      p_x = p_1 / (1 - (x-1) p_1)
 *
 * Args:     hmm    - the Plan7 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere;
 *                    will be distributed over M-1 match states.
 *
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
Plan7SWConfig(
  struct plan7_s *hmm,
  float pentry, float pexit);


/* Function: Plan7FSConfig()
 *
 * Purpose:  Set the alignment independent parameters of
 *           a Plan7 model to hmmfs (multihit Smith/Waterman) configuration.
 *
 *           See comments on Plan7SWConfig() for explanation of
 *           how pentry and pexit are used.
 *
 * Args:     hmm    - the Plan7 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere;
 *                    will be distributed over M-1 match states.
 *
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
Plan7FSConfig(
  struct plan7_s *hmm,
  float pentry,
  float pexit);


/* Function: DegenerateSymbolScore()
 *
 * Purpose:  Given a sequence character x and an hmm emission probability
 *           vector, calculate the log-odds (base 2) score of
 *           the symbol.
 *
 *           Easy if x is in the emission alphabet, but not so easy
 *           is x is a degenerate symbol. The "correct" Bayesian
 *           philosophy is to calculate score(X) by summing over
 *           p(x) for all x in the degenerate symbol X to get P(X),
 *           doing the same sum over the prior to get F(X), and
 *           doing log_2 (P(X)/F(X)). This gives an X a zero score,
 *           for instance.
 *
 *           Though this is correct in a formal Bayesian sense --
 *           we have no information on the sequence, so we can't
 *           say if it's random or model, so it scores zero --
 *           it sucks, big time, for scoring biological sequences.
 *           Sequences with lots of X's score near zero, while
 *           real sequences have average scores that are negative --
 *           so the X-laden sequences appear to be lifted out
 *           of the noise of a full histogram of a database search.
 *           Correct or not, this is highly undesirable.
 *
 *           So therefore we calculated the expected score of
 *           the degenerate symbol by summing over all x in X:
 *                 e_x log_2 (p(x)/f(x))
 *           where the expectation of x, e_x, is calculated from
 *           the random model.
 *
 *           Empirically, this works; it also has a wooly hand-waving
 *           probabilistic justification that I'm happy enough about.
 *
 * Args:     p      - probabilities of normal symbols
 *           null   - null emission model
 *           ambig  - index of the degenerate character in Alphabet[]
 *
 * Return:   the integer log odds score of x given the emission
 *           vector and the null model, scaled up by INTSCALE.
 */
int
DegenerateSymbolScore(
  float *p,
  float *null,
  int ambig);


/*****************************************************************
 *
 * Plan9/Plan7 interface
 *
 * Very important code during the evolutionary takeover by Plan7 --
 * convert between Krogh/Haussler and Plan7 models.
 *****************************************************************/

/* Function: Plan9toPlan7()
 *
 * Purpose:  Convert an old HMM into Plan7. Configures it in
 *           ls mode.
 *
 * Args:     hmm       - old ugly plan9 style HMM
 *           ret_plan7 - new wonderful Plan7 HMM
 *
 * Return:   (void)
 *           Plan7 HMM is allocated here. Free w/ FreePlan7().
 */
void
Plan9toPlan7(
  struct plan9_s *hmm,
  struct plan7_s **ret_plan7);

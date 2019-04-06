/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* prior.c
 *
 * Support for Dirichlet prior data structure, p7prior_s.
 */


#include <string.h>

#include "getopt.hpp"
#include "squid.hpp"
#include "structs.hpp"


/* Function: P7AllocPrior(), P7FreePrior()
 *
 * Purpose:  Allocation and free'ing of a prior structure.
 *           Very simple, but might get more complex someday.
 */
struct p7prior_s *
P7AllocPrior();


void
P7FreePrior(
  struct p7prior_s *pri);


/* Function: P7DefaultPrior()
 *
 * Purpose:  Set up a somewhat more realistic single component
 *           Dirichlet prior than Laplace.
 */
struct p7prior_s *
P7DefaultPrior();


/* Function: P7ReadPrior()
 *
 * Purpose:  Input a prior from disk file.
 */
struct p7prior_s *
P7ReadPrior(
  char *prifile);


/* Function: PAMPrior()
 *
 * Purpose:  Produces an ad hoc "Dirichlet mixture" prior for
 *           match emissions, using a PAM matrix.
 *
 *           Side effect notice: PAMPrior() replaces the match
 *           emission section of an existing Dirichlet prior,
 *           which is /expected/ to be a simple one-component
 *           kind of prior. The insert emissions /must/ be a
 *           one-component prior (because of details in how
 *           PriorifyEmissionVector() is done). However,
 *           the transitions /could/ be a mixture Dirichlet prior
 *           without causing problems. In other words, the
 *           -p and -P options of hmmb can coexist, but there
 *           may be conflicts. PAMPrior() checks for these,
 *           so there's no serious problem, except that the
 *           error message from PAMPrior() might be confusing to
 *           a user.
 */
void
PAMPrior(
  char *pamfile,
  struct p7prior_s *pri,
  float wt);


/* Function: P7DefaultNullModel()
 *
 * Purpose:  Set up a default random sequence model, using
 *           global aafq[]'s for protein or 1/Alphabet_size for anything
 *           else. randomseq is alloc'ed in caller. Alphabet information
 *           must already be known.
 */
void
P7DefaultNullModel(
  float *null,
  float *ret_p1);


void
P7ReadNullModel(
  char *rndfile,
  float *null,
  float *ret_p1);


/* Function: P7PriorifyHMM()
 *
 * Purpose:  Add pseudocounts to an HMM using Dirichlet priors,
 *           and renormalize the HMM.
 *
 * Args:     hmm -- the HMM to add counts to (counts form)
 *           pri -- the Dirichlet prior to use
 *
 * Return:   (void)
 *           HMM returns in probability form.
 */
void
P7PriorifyHMM(
  struct plan7_s *hmm,
  struct p7prior_s *pri);


/* Function: P7PriorifyEmissionVector()
 *
 * Purpose:  Add prior pseudocounts to an observed
 *           emission count vector and renormalize.
 *
 *           Can return the posterior mixture probabilities
 *           P(q | counts) if ret_mix[MAXDCHLET] is passed.
 *           Else, pass NULL.
 *
 * Args:     vec     - the 4 or 20-long vector of counts to modify
 *           pri     - prior data structure
 *           num     - pri->mnum or pri->inum; # of mixtures
 *           eq      - pri->mq or pri->iq; prior mixture probabilities
 *           e       - pri->i or pri->m; Dirichlet components
 *           ret_mix - filled with posterior mixture probabilities, or NULL
 *
 * Return:   (void)
 *           The counts in vec are changed and normalized to probabilities.
 */
void
P7PriorifyEmissionVector(
  float *vec,
  struct p7prior_s *pri,
  int num,
  float eq[MAXDCHLET],
  float e[MAXDCHLET][MAXABET],
  float *ret_mix);



/* Function: P7PriorifyTransitionVector()
 *
 * Purpose:  Add prior pseudocounts to transition vector,
 *           which contains three different probability vectors
 *           for m, d, and i.
 *
 * Args:     t     - state transitions, counts: 3 for M, 2 for I, 2 for D.
 *           prior - Dirichlet prior information
 *           tq    - prior distribution over Dirichlet components.
 *                   (overrides prior->iq[]; used for alternative
 *                   methods of conditioning prior on structural data)
 *
 * Return:   (void)
 *           t is changed, and renormalized -- comes back as
 *           probability vectors.
 */
void
P7PriorifyTransitionVector(
  float *t,
  struct p7prior_s *prior,
  float tq[MAXDCHLET]);


/* Function: default_amino_prior()
 *
 * Purpose:  Set the default protein prior.
 */
struct p7prior_s*
default_amino_prior();


/* Function: default_nucleic_prior()
 *
 * Purpose:  Set the default DNA prior. (for now, almost a Laplace)
 */
struct p7prior_s *
default_nucleic_prior();

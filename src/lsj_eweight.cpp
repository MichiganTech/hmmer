/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* lsjfuncs.c
 *
 * entropy targeting:
 * Code for setting effective sequence number (in hmmbuild) by
 * achieving a certain target entropy loss, relative to background
 * null distribution.
 */

#include <stdio.h>
#include <string.h>

#include "vectorops.hpp"
#include "lsj_eweight.hpp"


float
Eweight(
  struct plan7_s *hmm, 
  struct p7prior_s *pri,
  float numb_seqs,
  float targetent
){
  int i;
  int j;
  float eff_no;                  /* New effective sequence number */
  float current;                 /* Current mean match state entropy */
//  float prevent;                 /* Previous mean match state entropy */
  float scale;                   /* Current model counts scaling factor */
  float leftscale;               /* Bracket scaling value used in binary search. Lowest mean entropy value. */
  float rightscale;              /* Bracket scaling value used in binary search. Highest mean entropy value. */
  float *pmat;                   /* Temp array of match state counts */
  float *ent;                    /* Match state entropy values */
  int count;                     /* Counter for binary search */
//  int flag;                      /* Used to detect entropy adjustment failure */

  /**************
   * Allocations
   **************/
  pmat     = MallocOrDie (MAXABET * sizeof(float *));
  ent      = MallocOrDie ((hmm->M+1) * sizeof(float *));

  /*****************
   * Initializations
   *****************/
  current  = 0.;
  scale    = 1.;
  count    = 0;
//  flag     = 0;

  for(i = 0; i < Alphabet_size; i++) {
    pmat[i] = 0;
  }
  for(i = 0; i < hmm->M+1; i++) {
    ent[i] = 0;
  }

  /***************************************
   * Calculate the starting model entropy
   ***************************************/

  /* Copy model match state probabilities into our temporary pmat[] */
  for(i = 1; i < hmm->M+1; i++) {
    for(j = 0; j < Alphabet_size; j++) {
      pmat[j] = hmm->mat[i][j];
    }
    /* Add priors to the current match state prob dist. (prior.c) */
    P7PriorifyEmissionVector(pmat, pri, pri->mnum, pri->mq, pri->m, NULL);
    /* ent[] is assigned the current match state emmision entropy. */
    ent[i] = FEntropy(pmat, Alphabet_size);
  }
  /* Calculate the mean match state entropy. (squid/vectorops.c::FSum) */
  current = FSum(ent, hmm->M+1)/hmm->M;

  /****************************************
  * Initialize binary search bracket values
  *****************************************/

  /* The reason the values seem backwards is because I'm trying to
     bracket my target mean entropy with model count scaling
     factors. A higher scaling factor generally produces a lower
     entropy and a lower scaling factor produces a higher
     entropy. Thus, the leftscale produces the lowest mean entropy
     bracket and rightscale produces the highest mean entropy
     bracket */
  if(current < targetent) {
    leftscale  = 1;
    rightscale = 0;
  } else {
    /* Current model has a higher entropy than our target.
       Calculated effective seq numb <= Number of seqs. Design decision.
    */
    printf("[e=%.2f >= %.2f] ...", current, targetent);
    free(pmat);
    free(ent);
    return(numb_seqs);
  }
  /***************************************
   * Binary search for target mean entropy
   ***************************************/
  /* Check to see if the current model mean entropy is within 0.01 bits of our target */
  while((current < targetent - 0.01) || (current > targetent + 0.01)) {
    count++;

    /* Emergency brake in case there is a bug in our binary search */
    if(count > 50) {
      printf("\nBUG: Problem with adjusting the model entropy. Please report.\n");
      break;
    }

    /* Calculate current scaling factor based on bracket values */
    scale = (leftscale + rightscale)/2;

    //prevent = current;

    /*******************************************
     * Scale the counts and re-calc the entropy
     *******************************************/
    /* Re-copy match state probabilities into pmat[] */
    for(i = 1; i < hmm->M+1; i++) {
      for(j = 0; j < Alphabet_size; j++) {
        pmat[j] = hmm->mat[i][j];
      }
      /* Re-scale the current counts by the previously determined amount. (squid/vectorops.c) */
      FScale(pmat, Alphabet_size, scale);
      /* Re-add priors to these scaled counts. (prior.c) */
      P7PriorifyEmissionVector(pmat, pri, pri->mnum, pri->mq, pri->m, NULL);
      /* Again, ent[] is assigned the current match emission entropy */
      ent[i] = FEntropy(pmat, Alphabet_size);
    }
    current = FSum(ent, hmm->M+1)/hmm->M;

    /* Adjust the brackets according to the new mean entropy value */
    if(current < targetent) {
      leftscale = scale;
    } else {
      /* We overshot the target. Replace right bracket with the current scale */
      rightscale = scale;
    }
  }
  free(pmat);
  free(ent);
  /**********************************************************************************************
   * End of binary search
   *********************************************************************************************/
  eff_no = numb_seqs * scale;
  return(eff_no);
}

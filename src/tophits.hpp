/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* tophits.c
 *
 * Routines for storing, sorting, displaying high scoring hits
 * and alignments.
 *
 *****************************************************************************
 *
 * main API:
 *
 * AllocTophits()       - allocation
 * FreeTophits()        - free'ing
 * RegisterHit()        - put information about a hit in the list
 * GetRankedHit()       - recovers information about a hit
 * FullSortTophits()    - sorts the top H hits.
 *
 *****************************************************************************
 * Brief example of use:
 *
 *   struct tophit_s   *yourhits;   // list of hits
 *   struct fancyali_s *ali;        // (optional structure) alignment of a hit
 *
 *   yourhits = AllocTophits(200);
 *   (for every hit in a search) {
 *        if (do_alignments)
 *           ali = Trace2FancyAli();  // You provide a function/structure here
 *        if (score > threshold)
 *           RegisterHit(yourhits, ...)
 *     }
 *
 *   FullSortTophits(yourhits);      // Sort hits by evalue
 *   for (i = 0; i < 100; i++)       // Recover hits out in ranked order
 *     {
 *       GetRankedHit(yourhits, i, ...);
 *                                   // Presumably you'd print here...
 *     }
 *   FreeTophits(yourhits);
 ***************************************************************************
 *
 * Estimated storage per hit:
 *        coords:   16 bytes
 *        scores:    8 bytes
 * name/acc/desc:  192 bytes
 *     alignment: 1000 bytes   total = ~1200 bytes with alignment;
 *                                   = ~200 bytes without
 *     Designed for: 10^5 hits (20 MB) or 10^4 alignments (10 MB)
 */

//#include "squidconf.h"

#include <string.h>
#include <float.h>
#include <limits.h>

#include "config.hpp"
#include "structs.hpp"


/* Function: AllocTophits()
 *
 * Purpose:  Allocate a struct tophit_s, for maintaining
 *           a list of top-scoring hits in a database search.
 *
 * Args:     lumpsize - allocation lumpsize
 *
 * Return:   An allocated struct hit_s. Caller must free.
 */
struct tophit_s *
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


/* Function: RegisterHit()
 *
 * Purpose:  Add a new hit to a list of top hits.
 *
 *           "ali", if provided, is a pointer to allocated memory
 *           for an alignment output structure.
 *           Management is turned over to the top hits structure.
 *           Caller should not free them; they will be free'd by
 *           the FreeTophits() call.
 *
 *           In contrast, "name", "acc", and "desc" are copied, so caller
 *           is still responsible for these.
 *
 *           Number of args is unwieldy.
 *
 * Args:     h        - active top hit list
 *           key      - value to sort by: bigger is better
 *           pvalue   - P-value of this hit
 *           score    - score of this hit
 *           motherp  - P-value of parent whole sequence
 *           mothersc - score of parent whole sequence
 *           name     - name of target
 *           acc      - accession of target (may be NULL)
 *           desc     - description of target (may be NULL)
 *           sqfrom   - 1..L pos in target seq  of start
 *           sqto     - 1..L pos; sqfrom > sqto if rev comp
 *           sqlen    - length of sequence, L
 *           hmmfrom  - 0..M+1 pos in HMM of start
 *           hmmto    - 0..M+1 pos in HMM of end
 *           hmmlen   - length of HMM, M
 *           domidx   - number of this domain
 *           ndom     - total # of domains in sequence
 *           ali      - optional printable alignment info
 *
 * Return:   (void)
 *           hitlist is modified and possibly reallocated internally.
 */
void
RegisterHit(
  struct tophit_s *h,
  double key,
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


/* Function: GetRankedHit()
 *
 * Purpose:  Recover the data from the i'th ranked hit.
 *           Any of the data ptrs may be passed as NULL for fields
 *           you don't want. hitlist must have been sorted first.
 *
 *           name, acc, desc, and ali are returned as pointers, not copies;
 *           don't free them!
 */
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


/* Function: TophitsMaxName()
 *
 * Purpose:  Returns the maximum name length in a top hits list;
 *           doesn't need to be sorted yet.
 */
int
TophitsMaxName(
  struct tophit_s *h);


/* Function: FullSortTophits()
 *
 * Purpose:  Completely sort the top hits list. Calls
 *           qsort() to do the sorting, and uses
 *           hit_comparison() to do the comparison.
 *
 * Args:     h - top hits structure
 */
int
hit_comparison(
  const void *vh1,
  const void *vh2);


void
FullSortTophits(
  struct tophit_s *h);



/* Function: TophitsReport()
 *
 * Purpose:  Generate a printout summarizing how much
 *           memory is used by a tophits structure,
 *           how many hits are stored, and how much
 *           waste there is from not knowing nseqs.
 *
 * Args:     h    - the sorted tophits list
 *           E    - the cutoff in Evalue
 *           nseq - the final number of seqs used for Eval
 *
 * Return:   (void)
 *           Prints information on stdout
 */
void
TophitsReport(
  struct tophit_s *h,
  double E,
  int nseq);

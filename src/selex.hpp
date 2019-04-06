/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* selex.c
 *
 * SELEX  obsolete as the preferred HMMER/SQUID format
 * replaced by Stockholm format
 * selex support retained for backwards compatibility
 * kludged to use the MSA interface
 *
 * #=SA side chain % surface accessibility annotation supported
 *
 * major revision. #= special comments and aliinfo_s optional
 * alignment info support added. Support for #=CS (consensus
 * secondary structure), #=SS (individual secondary structure),
 * #=RF (reference coordinate system), #=SQ (per-sequence header info),
 * and #=AU ("author") added.
 *
 * Reading and writing aligned sequences to/from disk files.
 * Implements a new, broader specification of SELEX format
 * and supercedes alignio.c.
 *
 * SELEX format is documented in Docs/formats.tex.
 */

#pragma once

#include <ctype.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.hpp"
#include "msa.hpp"


/* Function: ReadSELEX()
 * Purpose:  Parse an alignment read from an open SELEX format
 *           alignment file. (SELEX is a single alignment format).
 *           Return the alignment, or NULL if we've already read the
 *           alignment or there's no alignment data in the file.
 *          
 * Limitations: SELEX is the only remaining multipass parser for
 *           alignment files. It cannot read from gzip or from stdin.
 *           It Die()'s here if you try. The reason for this
 *           that SELEX allows space characters as gaps, so we don't
 *           know the borders of an alignment block until we've seen
 *           the whole block. I could rewrite to allow single-pass
 *           parsing (by storing the whole block in memory) but
 *           since SELEX is now legacy, why bother.
 *          
 *           Note that the interface is totally kludged: fastest
 *           possible adaptation of old ReadSELEX() to the new
 *           MSA interface. 
 *
 * Args:     afp  - open alignment file
 *
 * Returns:  MSA *  - an alignment object
 *                    caller responsible for an MSAFree()
 *           NULL if no alignment data.         
 */
MSA*
ReadSELEX(
  MSAFILE *afp);


/* Function: WriteSELEX()
 * Purpose:  Write a SELEX file in multiblock format.
 *
 * Args:     fp  - file that's open for writing
 *           msa - multiple sequence alignment object 
 *
 * Returns:  (void)
 */
void
WriteSELEX(
  FILE *fp,
  MSA *msa);


/* Function: WriteSELEXOneBlock()
 * Purpose:  Write a SELEX alignment file in Pfam's single-block
 *           format style. A wrapper for actually_write_selex().
 *
 * Args:     fp - file that's open for writing
 *           msa- alignment to write
 *
 * Returns:  (void)
 */
void
WriteSELEXOneBlock(
  FILE *fp,
  MSA *msa);


/* Function: actually_write_selex()
 * Purpose:  Write an alignment in SELEX format to an open
 *           file. This is the function that actually does
 *           the work. The API's WriteSELEX() and
 *           WriteSELEXOneBlock() are wrappers.
 *
 * Args:     fp  - file that's open for writing
 *           msa - alignment to write
 *           cpl - characters to write per line in alignment block
 *
 * Returns:  (void)
 */
 void
 actually_write_selex(
  FILE *fp,
  MSA *msa,
  int cpl);


/* Function: copy_alignment_line()
 *
 * Purpose:  Given a line from an alignment file, and bounds lcol,rcol
 *           on what part of it may be sequence, save the alignment into
 *           aseq starting at position apos.
 *          
 *           name_rcol is set to the rightmost column this aseqs's name
 *           occupies; if name_rcol >= lcol, we have a special case in
 *           which the name intrudes into the sequence zone.
 */
int
copy_alignment_line(
  char *aseq,
  int apos,
  int name_rcol,
  char *buffer,
  int lcol,
  int rcol,
  char gapsym);


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
DealignAseqs(
  char **aseqs,
  int num,
  char ***ret_rseqs);


/* Function: IsSELEXFormat()
 *
 * Return TRUE if filename may be in SELEX format.
 *
 * Accuracy is sacrificed for speed; a TRUE return does
 * *not* guarantee that the file will pass the stricter
 * error-checking of ReadSELEX(). All it checks is that
 * the first 500 non-comment lines of a file are
 * blank, or if there's a second "word" on the line
 * it looks like sequence (i.e., it's not kOtherSeq).
 *
 * Returns TRUE or FALSE.
 */
bool
IsSELEXFormat(
  char *filename);
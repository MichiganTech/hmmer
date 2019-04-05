/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* File: sqio.c
 * From: ureadseq.c in Don Gilbert's sequence i/o package
 *
 * Reads and writes nucleic/protein sequence in various
 * formats. Data files may have multiple sequences.
 *
 * Heavily modified from READSEQ package
 * Copyright (C) 1990 by D.G. Gilbert
 * Biology Dept., Indiana University, Bloomington, IN 47405
 * email: gilbertd@bio.indiana.edu
 * Thanks Don!
 *
 * SRE: Modifications as noted. Fri Jul  3 09:44:54 1992
 *      Packaged for squid, Thu Oct  1 10:07:11 1992
 *      ANSI conversion in full swing, Mon Jul 12 12:22:21 1993
 *
 * CVS $Id: sqio.c,v 1.29 2002/08/26 23:10:52 eddy Exp $
 * 
 *****************************************************************
 * Basic API for single sequence reading:
 * 
 * SQFILE *sqfp;
 * char   *seqfile;
 * int     format;        - see squid.h for formats; example: SQFILE_FASTA
 * char   *seq;
 * SQINFO  sqinfo;
 * 
 * if ((sqfp = SeqfileOpen(seqfile, format, "BLASTDB")) == NULL)
 *   Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);
 * while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo)) {
 *   do_stuff;
 *   FreeSequence(seq, &sqinfo);
 * }
 * SeqfileClose(sqfp);  
 * 
 *****************************************************************  
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <unistd.h>	

#include "squid.h"
#include "msa.h"
#include "ssi.h"


void 
SeqfileGetLine(
  SQFILE *V);


SQFILE*
seqfile_open(
  char *filename, 
  int format, 
  char *env, 
  int ssimode);


/* Function: SeqfileOpen()
 * 
 * Purpose : Open a sequence database file and prepare for reading
 *           sequentially.
 *           
 * Args:     filename - name of file to open
 *           format   - format of file
 *           env      - environment variable for path (e.g. BLASTDB)   
 *           ssimode  - -1, SSI_OFFSET_I32, or SSI_OFFSET_I64
 *
 *           Returns opened SQFILE ptr, or NULL on failure.
 */
SQFILE *
SeqfileOpen(
  char *filename, 
  int format, 
  char *env);


SQFILE *
SeqfileOpenForIndexing(
  char *filename, 
  int format, 
  char *env, 
  int ssimode);


SQFILE *
seqfile_open(
  char *filename, 
  int format, 
  char *env, 
  int ssimode);

/* Function: SeqfilePosition()
 * 
 * Purpose:  Move to a particular offset in a seqfile.
 *           Will not work on alignment files.
 */
void
SeqfilePosition(
  SQFILE *sqfp, 
  SSIOFFSET *offset);


/* Function: SeqfileRewind()
 * 
 * Purpose:  Set a sequence file back to the first sequence.
 * 
 *           Won't work on alignment files. Although it would
 *           seem that it could (just set msa->lastidx back to 0),
 *           that'll fail on "multiple multiple" alignment file formats
 *           (e.g. Stockholm).
 */
void
SeqfileRewind(
  SQFILE *sqfp);


/* Function: SeqfileLineParameters()
 *
 * Purpose:  After all the sequences have been read from the file,
 *           but before closing it, retrieve overall bytes-per-line and
 *           residues-per-line info. If non-zero, these mean that
 *           the file contains homogeneous sequence line lengths (except
 *           the last line in each record).  
 *           
 *           If either of bpl or rpl is determined to be inhomogeneous,
 *           both are returned as 0.
 *
 * Args:     *sqfp   - an open but fully read sequence file
 *           ret_bpl - RETURN: bytes per line, or 0 if inhomogeneous
 *           ret_rpl - RETURN: residues per line, or 0 if inhomogenous.
 *
 * Returns:  void
 */
void
SeqfileLineParameters(
  SQFILE *V, 
  int *ret_bpl, 
  int *ret_rpl);


void
SeqfileClose(
  SQFILE *sqfp);


/* Function: SeqfileGetLine()
 * Date:     SRE, Tue Jun 22 09:15:49 1999 [Sanger Centre]
 *
 * Purpose:  read a line from a sequence file into V->buf
 *           If the fgets() is NULL, sets V->buf[0] to '\0'.
 *
 * Args:     V
 *
 * Returns:  void
 */
void 
SeqfileGetLine(
  SQFILE *V);


void
FreeSequence(
  char *seq, 
  SQINFO *sqinfo);


int
SetSeqinfoString(
  SQINFO *sqinfo, 
  char *sptr, 
  int flag);


void
SeqinfoCopy(
  SQINFO *sq1, 
  SQINFO *sq2);


/* Function: ToDNA()
 * 
 * Purpose:  Convert a sequence to DNA.
 *           U --> T
 */
void
ToDNA(
  char *seq);


/* Function: ToRNA()
 * 
 * Purpose:  Convert a sequence to RNA.
 *           T --> U
 */
void
ToRNA(
  char *seq);


/* Function: ToIUPAC()
 * 
 * Purpose:  Convert X's, o's, other junk in a nucleic acid sequence to N's,
 *           to comply with IUPAC code. If is_aseq is TRUE, will allow gap
 *           characters though, so we can call ToIUPAC() on aligned seqs.
 *      
 *           NUCLEOTIDES is defined in squid.h as:
 *               "ACGTUNRYMKSWHBVDacgtunrymkswhbvd"
 *           gap chars allowed by isgap() are defined in squid.h as:
 *               " ._-~"
 *
 *           WU-BLAST's pressdb will
 *           choke on X's, for instance, necessitating conversion
 *           of certain genome centers' data. 
 */
void
ToIUPAC(
  char *seq, 
  int is_aseq);


/* Function: addseq()
 * 
 * Purpose:  Add a line of sequence to the growing string in V.
 *
 *           In the seven supported unaligned formats, all sequence
 *           lines may contain whitespace that must be filtered out;
 *           four formats (PIR, EMBL, Genbank, GCG) include coordinates
 *           that must be filtered out. Thus an (!isdigit && !isspace)
 *           test on each character before we accept it.
 */
void 
addseq(
  char *s, 
  struct ReadSeqVars *V);


void 
readLoop(
  int addfirst, 
  int (*endTest)(char *,int *), 
  struct ReadSeqVars *V);


int
endPIR(
  char *s, 
  int  *addend);


void
readPIR(
  struct ReadSeqVars *V);


int 
endIG(
  char *s, 
  int  *addend);


void 
readIG(
  struct ReadSeqVars *V);


int 
endStrider(
  char *s, 
  int *addend);


void 
readStrider(
  struct ReadSeqVars *V);


int 
endGB(
  char *s, 
  int *addend);


void 
readGenBank(
  struct ReadSeqVars *V);


int
endGCGdata(
  char *s, 
  int *addend);


void
readGCGdata(
  struct ReadSeqVars *V);


int
endPearson(
  char *s, 
  int *addend);


void 
readPearson(
  struct ReadSeqVars *V);


int
endEMBL(
  char *s, 
  int *addend);


void 
readEMBL(
  struct ReadSeqVars *V);


int
endZuker(
  char *s, 
  int *addend);


void
readZuker(
  struct ReadSeqVars *V);


void 
readUWGCG(
  struct ReadSeqVars *V);

    
/* Function: ReadSeq()
 * 
 * Purpose:  Read next sequence from an open database file.
 *           Return the sequence and associated info.
 *           
 * Args:     fp      - open sequence database file pointer          
 *           format  - format of the file (previously determined
 *                      by call to SeqfileFormat()).
 *                      Currently unused, since we carry it in V. 
 *           ret_seq - RETURN: sequence
 *           sqinfo  - RETURN: filled in w/ other information  
 *                     
 * Limitations: uses squid_errno, so it's not threadsafe.                    
 *           
 * Return:   1 on success, 0 on failure.
 *           ret_seq and some field of sqinfo are allocated here,
 *           The preferred call mechanism to properly free the memory is:
 *           
 *           SQINFO sqinfo;
 *           char  *seq;
 *           
 *           ReadSeq(fp, format, &seq, &sqinfo);
 *           ... do something...
 *           FreeSequence(seq, &sqinfo);
 */
int
ReadSeq(
  SQFILE *V, 
  //int format, 
  char **ret_seq, 
  SQINFO *sqinfo);


/* Function: SeqfileFormat()
 *
 * Purpose:  Determine format of an open file.
 *           Returns format code.
 *           Rewinds the file.
 *           
 *           Autodetects the following unaligned formats:
 *              SQFILE_FASTA
 *              SQFILE_GENBANK
 *              SQFILE_EMBL
 *              SQFILE_GCG
 *              SQFILE_GCGDATA
 *              SQFILE_PIR
 *           Also autodetects the following alignment formats:
 *              MSAFILE_STOCKHOLM
 *              MSAFILE_MSF
 *              MSAFILE_CLUSTAL
 *              MSAFILE_SELEX
 *              MSAFILE_PHYLIP
 *
 *           Can't autodetect MSAFILE_A2M, calls it SQFILE_FASTA.
 *           MSAFileFormat() does the opposite.
 *
 * Args:     sfp -  open SQFILE
 *           
 * Return:   format code, or SQFILE_UNKNOWN if unrecognized
 */          
int
SeqfileFormat(
  FILE *fp);


/* Function: GCGBinaryToSequence()
 * 
 * Purpose:  Convert a GCG 2BIT binary string to DNA sequence.
 *           0 = C  1 = T  2 = A  3 = G
 *           4 nts/byte
 *           
 * Args:     seq  - binary sequence. Converted in place to DNA.
 *           len  - length of DNA. binary is (len+3)/4 bytes
 */
int
GCGBinaryToSequence(
  char *seq, 
  int len);


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
GCGchecksum(
  char *seq, 
  int len);


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
GCGMultchecksum(
  char **seqs, 
  int nseq);


/* Function: Seqtype()
 * 
 * Purpose:  Returns a (very good) guess about type of sequence:
 *           kDNA, kRNA, kAmino, or kOtherSeq.
 *           
 *           Modified from, and replaces, Gilbert getseqtype().
 */
int
Seqtype(
  char *seq);


/* Function: GuessAlignmentSeqtype()
 *
 * Purpose:  Try to guess whether an alignment is protein 
 *           or nucleic acid; return a code for the
 *           type (kRNA, kDNA, or kAmino).
 *
 * Args:     aseq  - array of aligned sequences. (Could also
 *                   be an rseq unaligned sequence array)
 *           nseq  - number of aseqs
 *
 * Returns:  kRNA, kDNA, kAmino;
 *           kOtherSeq if inconsistency is detected.
 */
int
GuessAlignmentSeqtype(
  char **aseq, 
  int nseq);


/* Function: WriteSimpleFASTA()
 *
 * Purpose:  Just write a FASTA format sequence to a file;
 *           minimal interface, mostly for quick and dirty programs.
 *
 * Args:     fp   - open file handle (stdout, possibly)
 *           seq  - sequence to output
 *           name - name for the sequence
 *           desc - optional description line, or NULL.
 *
 * Returns:  void
 */
void
WriteSimpleFASTA(
  FILE *fp, 
  char *seq, 
  char *name, 
  char *desc);


int
WriteSeq(
  FILE *outf, 
  int outform, 
  char *seq, 
  SQINFO *sqinfo);


/* Function: ReadMultipleRseqs()
 * 
 * Purpose:  Open a data file and
 *           parse it into an array of rseqs (raw, unaligned
 *           sequences).
 * 
 *           Caller is responsible for free'ing memory allocated
 *           to ret_rseqs, ret_weights, and ret_names.
 *           
 *           Weights are currently only supported for MSF format.
 *           Sequences read from all other formats will be assigned
 *           weights of 1.0. If the caller isn't interested in
 *           weights, it passes NULL as ret_weights.
 * 
 * Returns 1 on success. Returns 0 on failure and sets
 * squid_errno to indicate the cause.
 */
int
ReadMultipleRseqs(
  char *seqfile,
	int fformat,
	char ***ret_rseqs,
  SQINFO **ret_sqinfo,
  int *ret_num);




int
String2SeqfileFormat(
  char *s);


char *
SeqfileFormat2String(
  int code);


/* Function: MSAToSqinfo()
 *
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
  MSA *msa);

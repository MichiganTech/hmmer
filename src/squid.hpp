/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

#pragma once

/* squid.h
 * Header file for my library of sequence functions.
 */


#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "msa.hpp"    /* for multiple sequence alignment support   */

#include "ssi.h"

/* Library version info is made available as a global to
 * any interested program. These are defined in iupac.c
 * with the other globals.
 */
extern char* squid_version;	/* version number  */
extern char* squid_date;	/* date of release */
extern int  squid_errno;	/* error codes     */


#define kStartLength  500


/****************************************************
 * Error codes returned by squid library functions (squid_errno)
 ****************************************************/

#define SQERR_OK        0	/* no error                     */
#define SQERR_UNKNOWN   1       /* generic error, unidentified  */
#define SQERR_NODATA    2	/* unexpectedly NULL stream     */
#define SQERR_MEM       3	/* malloc or realloc failed     */
#define SQERR_NOFILE    4	/* file not found               */
#define SQERR_FORMAT    5	/* file format not recognized   */
#define SQERR_PARAMETER 6	/* bad parameter passed to func */
#define SQERR_DIVZERO   7	/* error in sre_math.c          */
#define SQERR_INCOMPAT  8	/* incompatible parameters      */
#define SQERR_EOD       9	/* end-of-data (often normal)   */

/****************************************************
 * Single sequence information
 ****************************************************/
#define SQINFO_NAMELEN 64
#define SQINFO_DESCLEN 128

struct seqinfo_s {
  int      flags;               /* what extra data are available         */
  char     name[SQINFO_NAMELEN];/* up to 63 characters of name           */
  char     id[SQINFO_NAMELEN];	/* up to 63 char of database identifier  */
  char     acc[SQINFO_NAMELEN]; /* up to 63 char of database accession # */
  char     desc[SQINFO_DESCLEN];/* up to 127 char of description         */
  int      len;                 /* length of this seq                    */
  int      start;		/* (1..len) start position on source seq */
  int      stop;                /* (1..len) end position on source seq   */
  int      olen;                /* original length of source seq         */
  int      type;                /* kRNA, kDNA, kAmino, or kOther         */
  char    *ss;                  /* 0..len-1 secondary structure string   */
  char    *sa;			/* 0..len-1 % side chain surface access. */
};
typedef struct seqinfo_s SQINFO;

#define SQINFO_NAME  (1 << 0)
#define SQINFO_ID    (1 << 1)
#define SQINFO_ACC   (1 << 2)
#define SQINFO_DESC  (1 << 3)
#define SQINFO_START (1 << 4)
#define SQINFO_STOP  (1 << 5)
#define SQINFO_LEN   (1 << 6)
#define SQINFO_TYPE  (1 << 7)
#define SQINFO_OLEN  (1 << 8)
#define SQINFO_SS    (1 << 9)
#define SQINFO_SA    (1 << 10)


/****************************************************
 * Sequence alphabet: see also iupac.c
 ****************************************************/
				/* IUPAC symbols defined globally in iupac.c */
struct iupactype {
  char       sym;		/* character representation */
  char       symcomp;           /* complement (regular char */
  char       code;		/* my binary rep */
  char       comp;              /* binary encoded complement */
};
extern struct iupactype *iupac;
#define IUPACSYMNUM 17

extern char**    stdcode1;	/* 1-letter amino acid translation code */
extern char**    stdcode3;	/* 3-letter amino acid translation code */
extern float*    dnafq;        /* nucleotide occurrence frequencies    */
extern float*    aafq;		/* amino acid occurrence frequencies    */
extern char*     aa_alphabet;  /* amino acid alphabet                  */
extern int*      aa_index;     /* convert 0..19 indices to 0..26       */

				/* valid symbols in IUPAC code */
#define NUCLEOTIDES    "ACGTUNRYMKSWHBVDacgtunrymkswhbvd"
#define AMINO_ALPHABET "ACDEFGHIKLMNPQRSTVWY"
#define DNA_ALPHABET   "ACGT"
#define RNA_ALPHABET   "ACGU"
#define WHITESPACE     " \t\n"

#define isgap(c) ((c) == ' ' || (c) == '.' || (c) == '_' || (c) == '-' || (c) == '~')


/****************************************************
 * Sequence i/o: originally from Don Gilbert's readseq
 ****************************************************/

	/* buffer size for reading in lines from sequence files*/
#define LINEBUFLEN  4096

/* sequence types parsed by Seqtype()                          */
/* note that these must match hmmAMINO and hmmNUCLEIC in HMMER */
#define kOtherSeq   0		/* hmmNOTSETYET */
#define kDNA        1
#define kRNA        2		/* hmmNUCLEIC   */
#define kAmino      3		/* hmmAMINO     */

/* Unaligned sequence file formats recognized
 * Coexists with definitions of multiple alignment formats in msa.h:
 *   >100 reserved for alignment formats
 *   <100 reserved for unaligned formats
 *   0 reserved for unknown
 *  
 * Some "legacy" formats are supported only when explicitly
 * requested; not autodetected by SeqfileFormat().
 *
 * DON'T REASSIGN THESE CODES. They're written into
 * GSI index files. You can use new ones, but reassigning
 * the sense of old ones will break GSI indices.
 * Alignment format codes were reassigned with the creation
 * of msa.c, but before Stockholm format, there were no
 * indexed alignment databases.
 */
#define SQFILE_UNKNOWN  0	/* unknown format                  */
#define SQFILE_IG       1	/* Intelligenetics (!)             */
#define SQFILE_GENBANK  2	/* GenBank flatfile                */
				/* 3 was A2M. Now an alignment format  */
#define SQFILE_EMBL     4	/* EMBL or Swissprot flatfile      */
#define SQFILE_GCG      5	/* GCG single sequence files       */
#define SQFILE_STRIDER  6	/* MacStrider (!!)                 */
#define SQFILE_FASTA    7	/* FASTA format: default           */
#define SQFILE_ZUKER    8	/* Zuker MFOLD format (legacy)     */
#define SQFILE_IDRAW    9	/* Idraw-style PostScript (legacy) */
				/* 10 was SELEX. Now alignment format  */
				/* 11 was MSF. Now alignment format    */
#define SQFILE_PIR      12	/* PIR format                      */
#define SQFILE_RAW      13	/* raw sequence                    */
#define SQFILE_SQUID    14	/* my obsolete squid format        */
				/* 15 was kXPearson, extended FASTA; withdrawn */
#define SQFILE_GCGDATA  16	/* GCG data library file           */
				/* 17 was Clustal. Now alignment format*/

#define IsUnalignedFormat(fmt)  ((fmt) && (fmt) < 100)



struct ReadSeqVars {
  FILE   *f;                    /* open file pointer                  */
  char   *fname;                /* name of file; used for diagnostics */
  size_t     linenumber;           /* what line are we on in the file    */

  char   *buf;                  /* dynamically allocated sre_fgets() buffer */
  size_t     buflen;               /* allocation length for buf                */
 
  int       ssimode;		/* SSI_OFFSET_I32 or SSI_OFFSET_I64        */
  SSIOFFSET ssioffset;		/* disk offset to last line read into buf  */
  SSIOFFSET r_off;		/* offset to start of record               */
  SSIOFFSET d_off;		/* offset to start of sequence data        */

  int     rpl;			/* residues per data line for this file; -1 if unset, 0 if invalid */
  int     lastrpl;		/* rpl on last line seen */
  int     maxrpl;		/* max rpl on any line of the file */
  int     bpl;			/* bytes per data line; -1 if unset, 0 if invalid */
  int     lastbpl;		/* bpl on last line seen */
  int     maxbpl;		/* max bpl on any line of the file */

  char   *seq;                  /* growing sequence during parse */
  SQINFO *sqinfo;	        /* name, id, etc, gathered during parse */
  char   *sp;
  size_t     seqlen;		/* current sequence length */
  size_t     maxseq;		/* current allocation length for seq */

  int     format;		/* format of seqfile we're reading. */
  bool     do_gzip;		/* TRUE if f is a pipe from gzip -dc */
  bool     do_stdin;		/* TRUE if f is stdin */

  /* An (important) hack for sequential access of multiple alignment files:
   * we read the whole alignment in,
   * and then copy it one sequence at a time into seq and sqinfo.
   * It is active if msa is non NULL.
   * msa->lastidx is reused/overloaded: used to keep track of what
   * seq we'll return next.
   * afp->format is the real format, while SQFILE->format is kMSA.
   * Because we keep it in the SQFILE structure,
   * ReadSeq() and friends are always reentrant for multiple seqfiles.
   */
  MSA      *msa;
  MSAFILE  *afp;
};
typedef struct ReadSeqVars SQFILE;


/****************************************************
 * Cluster analysis and phylogenetic tree support
 ****************************************************/

/* struct phylo_s - a phylogenetic tree
 *                    
 * For N sequences, there will generally be an array of 0..N-2
 * phylo_s structures representing the nodes of a tree.
 * [0] is the root. The indexes of left and
 * right children are somewhat confusing so be careful. The
 * indexes can have values of 0..2N-2. If they are 0..N-1, they
 * represent pointers to individual sequences. If they are
 * >= N, they represent pointers to a phylo_s structure
 * at (index - N).
 */
struct phylo_s {
  int    parent;    /* index of parent, N..2N-2, or -1 for root */
  int    left;			/* index of one of the branches, 0..2N-2 */
  int    right;			/* index of other branch, 0..2N-2        */
  float  diff;			/* difference score between seqs         */
  float  lblen;     /* left branch length                    */
  float  rblen;     /* right branch length                   */
  char  *is_in;     /* 0..N-1 flag array, 1 if seq included  */
  int    incnum;    /* number of seqs included at this node  */
};


/* Strategies for cluster analysis; cluster by mean distance,
 * minimum distance, or maximum distance.
 */
enum clust_strategy { CLUSTER_MEAN, CLUSTER_MAX, CLUSTER_MIN };

/****************************************************
 * Generic data structure support
 ****************************************************/

/* a struct intstack_s implements a pushdown stack for storing
 * single integers.
 */
struct intstack_s {
  int data;
  struct intstack_s *nxt;
};

/****************************************************
 * Binary nucleotide alphabet support
 ****************************************************/

/* Binary encoding of the IUPAC code for nucleotides
 *
 *    four-bit "word", permitting rapid degenerate matching
 *         A  C  G  T/U
 *         0  0  1  0
 */
#define NTA 8
#define NTC 4
#define NTG 2
#define NTT 1
#define NTU 1
#define NTN 15			/* A|C|G|T */
#define NTR 10			/* A|G */
#define NTY 5			/* C|T */
#define NTM 12			/* A|C */
#define NTK 3			/* G|T */
#define NTS 6			/* C|G */
#define NTW 9			/* A|T */
#define NTH 13			/* A|C|T */
#define NTB 7			/* C|G|T */
#define NTV 14			/* A|C|G */
#define NTD 11			/* A|G|T */
#define NTGAP 16		/* GAP */
#define NTEND 0			/* null string terminator */

/* ntmatch(): bitwise comparison of two nuc's
 * note that it's sensitive to the order;
 * probe may be degenerate but target should not be
 */
#define ntmatch(probe, target)  ((probe & target) == target)


/****************************************************
 * Support for convenient Perl-y regexp matching
 * See hsregexp.c for copyright notice: this code is derived
 * from Henry Spencer's freely distributed regexp library.
 ****************************************************/

#define NSUBEXP  10
typedef struct sqd_regexp {
	char *startp[NSUBEXP];
	char *endp[NSUBEXP];
	char regstart;		/* Internal use only. */
	char reganch;		/* Internal use only. */
	char *regmust;		/* Internal use only. */
	int regmlen;		/* Internal use only. */
	char program[1];	/* Unwarranted chumminess with compiler. */
} sqd_regexp;


/* Strparse() defines and manages these.
 * sqd_parse[0] contains the substring that matched the pattern.
 * sqd_parse[1-9] contain substrings matched with ()'s.
 */
extern char *sqd_parse[10];




/* Function: Warn()
 *
 * Purpose:  Print an error message and return. The arguments
 *           are formatted exactly like arguments to printf().
 *          
 * Return:   (void)
 */         
/* VARARGS0 */
void
Warn(
  char *format,
  ...);


/* Function: Panic()
 *
 * Purpose:  Die from a lethal error that's not my problem,
 *           but instead a failure of a StdC/POSIX call that
 *           shouldn't fail. Call perror() to get the
 *           errno flag, then die.
 *          
 *           Usually called by the PANIC macro which adds
 *           the __FILE__ and __LINE__ information; see
 *           structs.h.
 *          
 *           Inspired by code in Donald Lewine's book, _POSIX
 *           Programmer's Guide_.
 */
void
Panic(
  char *file,
  int line);



/* PANIC is called for failures of Std C/POSIX functions,
 * instead of my own functions. Panic() calls perror() and exits
 * abnormally.
 */
#define PANIC   Panic(__FILE__, __LINE__)

/* Malloc/realloc calls are wrapped
 */
#define MallocOrDie(x)     sre_malloc(__FILE__, __LINE__, (x))
#define ReallocOrDie(x,y)  sre_realloc(__FILE__, __LINE__, (x), (y))

/****************************************************
 * Miscellaneous macros and defines
 ****************************************************/

#define SQDCONST_E    2.71828182845904523536028747135
#define SQDCONST_PI   3.14159265358979323846264338328

				/* must declare swapfoo to use SWAP() */
#define SWAP(a,b) {swapfoo = b; b = a; a = swapfoo;}
#define ScalarsEqual(a,b) (fabs((a)-(b)) < 1e-7)

#ifndef MIN
#define MIN(a,b)         (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)         (((a)>(b))?(a):(b))
#endif



int  
String2SeqfileFormat(
  char *s);

sqd_regexp*
sqd_regcomp(
  const char *re);

int        
sqd_regexec(
  sqd_regexp *rp,
  const char *s);


void
sqd_regsub(const sqd_regexp *rp, const char *src, char *dst);

void*
sre_malloc(
  char *file,
  int line,
  size_t size);

void*
sre_realloc(
  char *file,
  int line,
  void *p,
  size_t size);

int  
ReadMultipleRseqs(
  char *seqfile,
  int fformat,
  char ***ret_rseqs,
  SQINFO **ret_sqinfo,
  int *ret_num);


/* Function: Warn()
 *
 * Purpose:  Print an error message and return. The arguments
 *           are formatted exactly like arguments to printf().
 *          
 * Return:   (void)
 */         
void
Warn(
  char *format,
  ...);


/* Function: Die()
 *
 * Purpose:  Print an error message and die. The arguments
 *           are formatted exactly like arguments to printf().
 *          
 * Return:   None. Exits the program.
 */         
void
Die(
  char *format,
  ...);


/* Function: Strparse()
 *
 * Purpose:  Match a regexp to a string. Returns 1 if pattern matches,
 *           else 0.
 *
 *           Much like Perl, Strparse() makes copies of the matching
 *           substrings available via globals, sqd_parse[].
 *           sqd_parse[0] contains a copy of the complete matched
 *           text. sqd_parse[1-9] contain copies of up to nine
 *           different substrings matched within parentheses.
 *           The memory for these strings is internally managed and
 *           volatile; the next call to Strparse() may destroy them.
 *           If the caller needs the matched substrings to persist
 *           beyond a new Strparse() call, it must make its own
 *           copies.
 *          
 *           A minor drawback of the memory management is that
 *           there will be a small amount of unfree'd memory being
 *           managed by Strparse() when a program exits; this may
 *           confuse memory debugging (Purify, dbmalloc). The
 *           general cleanup function SqdClean() is provided;
 *           you can call this before exiting.
 *          
 *           Uses an extended POSIX regular expression interface.
 *           A copylefted GNU implementation is included in the squid
 *           implementation (gnuregex.c) for use on non-POSIX compliant
 *           systems. POSIX 1003.2-compliant systems (all UNIX,
 *           some WinNT, I believe) can omit the GNU code if necessary.
 *          
 *           I built this for ease of use, not speed nor efficiency.
 *
 * Example:  Strparse("foo-...-baz", "foo-bar-baz")  returns 0
 *           Strparse("foo-(...)-baz", "foo-bar-baz")
 *              returns 0; sqd_parse[0] is "foo-bar-baz";
 *              sqd_parse[1] is "bar".
 *             
 *           A real example:  
 *            s   = ">gnl|ti|3 G10P69425RH2.T0 {SUB 81..737}  /len=657"
 *            pat = "SUB ([0-9]+)"
 *            Strparse(pat, s, 1)
 *               returns 1; sqd_parse[1] is "81".
 *             
 * Args:     rexp  - regular expression, extended POSIX form
 *           s     - string to match against
 *           ntok  - number of () substrings we will save (maximum NSUBEXP-1)
 *                  
 * Return:   1 on match, 0 if no match
 */
int
Strparse(
  char *rexp,
  char *s,
  int ntok);


void*
sre_malloc(
  char *file,
  int line,
  size_t size);


void*
sre_realloc(
  char *file,
  int line,
  void *p,
  size_t size);

/* Function: EnvFileOpen()         
 *
 * Purpose:  Open a file, given a file name and an environment
 *           variable that contains a directory path. Files
 *           are opened read-only. Does not look at current directory
 *           unless "." is explicitly in the path specified by env.
 *          
 *           For instance:
 *             fp = EnvFileOpen("BLOSUM45", "BLASTMAT", NULL);
 *           or:
 *             fp = EnvFileOpen("swiss", "BLASTDB", NULL); 
 *            
 *           Environment variables may contain a colon-delimited
 *           list of more than one path; e.g.
 *             setenv BLASTDB /nfs/databases/foo:/nfs/databases/bar
 *            
 *           Sometimes a group of files may be found in
 *           one directory; for instance, an index file with a
 *           database. The caller can EnvFileOpen() the main
 *           file, and ask to get the name of the
 *           directory back in ret_dir, so it can construct
 *           the other auxiliary file names and fopen() them. (If it called
 *           EnvFileOpen(), it might get confused by
 *           file name clashes and open files in different
 *           directories.
 *            
 * Args:     fname   - name of file to open
 *           env     - name of environment variable containing path
 *           ret_dir - if non-NULL, RETURN: name of dir that was used.
 * 
 * Return:   FILE * to open file, or NULL on failure -- same as fopen()
 *           Caller must free ret_dir if it passed a non-NULL address.
 */
FILE*
EnvFileOpen(
  char *fname,
  char *env,
  char **ret_dir);



/* Function: IsBlankline()
 *
 * Purpose:  Returns TRUE if string consists solely of whitespace.
 *
 * Args:     s   - string to check
 */
bool
IsBlankline(
  char *s);


/* Function: IsInt()
 *
 * Returns TRUE if s points to something that atoi() will parse
 * completely and convert to an integer.
 */
bool
IsInt(
  char *s);


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


/* Function: SSIGetFilePosition()
 *
 * Purpose:  Fills {ret_offset} with the current disk
 *           offset of {fp}, relative to the start of the file.
 *           {mode} is set to either SSI_OFFSET_I32 or
 *           SSI_OFFSET_I64. If {mode} is _I32 (32 bit), just wraps
 *           a call to ftell(); otherwise deals with system-dependent
 *           details of 64-bit file offsets.
 *
 * Args:     fp         - open stream
 *           mode       - SSI_OFFSET_I32 or SSI_OFFSET_I64
 *           ret_offset - RETURN: file position      
 *
 * Returns:  0 on success. nonzero on error.
 */
int
SSIGetFilePosition(
  FILE *fp,
  int mode,
  SSIOFFSET *ret_offset);


/* Function: SeqfileGetLine()
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
SQFILE*
seqfile_open(
  char *filename,
  int format,
  char *env,
  int ssimode);


SQFILE*
SeqfileOpen(
  char *filename,
  int format,
  char *env);


SQFILE*
SeqfileOpenForIndexing(
  char *filename,
  int format,
  char *env,
  int ssimode);


void
SeqfileClose(
  SQFILE *sqfp);


void
FreeSequence(
  char *seq,
  SQINFO *sqinfo);


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
addseq(char *s, struct ReadSeqVars *V);


int
SetSeqinfoString(
  SQINFO *sqinfo,
  char *sptr,
  int flag);


int
GCGBinaryToSequence(
  char *seq,
  int len);


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
  int  fformat,
  char ***ret_rseqs,
  SQINFO **ret_sqinfo,
  int  *ret_num);


double
Gammln(
  double xx);


/* Function: FileConcat()
 *
 * Purpose:  Concatenate a directory path and a file name,
 *           returning a pointer to a malloc'ed string with the
 *           full filename. This isn't just a string concat,
 *           because we're careful about the dir slash.
 */
char *
FileConcat(
  char *dir,
  char *file);


void
s2upper(
  char *s);


/* Function: StringChop()
 *
 * Purpose:  Chop trailing whitespace off of a string.
 */
void
StringChop(
  char *s);



/* Function: Free2DArray(), Free3DArray()
 *
 * Purpose:  Convenience functions for free'ing 2D
 *           and 3D pointer arrays. Tolerates any of the
 *           pointers being NULL, to allow "sparse"
 *           arrays.
 *
 * Args:     p     - array to be freed
 *           dim1  - n for first dimension
 *           dim2  - n for second dimension
 *
 *           e.g. a 2d array is indexed p[0..dim1-1][]
 *                a 3D array is indexed p[0..dim1-1][0..dim2-1][]
 *          
 * Returns:  void
 *
 * Diagnostics: (void)
 *              "never fails"
 */
void
Free2DArray(
  void **p,
  int dim1);


void
Free3DArray(
  void ***p,
  int dim1,
  int dim2);


/* Function: ParsePAMFile()
 *
 * Purpose:  Given a pointer to an open file containing a PAM matrix,
 *           parse the file and allocate and fill a 2D array of
 *           floats containing the matrix. The PAM file is
 *           assumed to be in the format that NCBI distributes
 *           with BLAST. BLOSUM matrices also work fine, as
 *           produced by Henikoff's program "MATBLAS".
 *         
 *           Parses both old format and new format BLAST matrices.
 *           Old format just had rows of integers.
 *           New format includes a leading character on each row.
 *
 *           The PAM matrix is a 27x27 matrix, 0=A..25=Z,26=*.
 *           Note that it's not a 20x20 matrix as you might expect;
 *           this is for speed of indexing as well as the ability
 *           to deal with ambiguous characters.
 *          
 * Args:     fp        - open PAM file
 *           ret_pam   - RETURN: pam matrix, integers                  
 *           ret_scale - RETURN: scale factor for converting
 *                       to real Sij. For instance, PAM120 is
 *                       given in units of ln(2)/2. This may
 *                       be passed as NULL if the caller
 *                       doesn't care.
 *
 * Returns:  1 on success; 0 on failure and sets squid_errno to
 *           indicate the cause. ret_pam is allocated here and
 *           must be freed by the caller (use FreePAM).
 */
int
ParsePAMFile(
  FILE *fp,
  int ***ret_pam,
  float *ret_scale);


void
sqd_regerror(
  char *s);



/* Function: SqdClean()
 *
 * Purpose:  Clean up any squid library allocations before exiting
 *           a program, so we don't leave unfree'd memory around
 *           and confuse a malloc debugger like Purify or dbmalloc.
 */
void
SqdClean();


/* Function: String2SeqfileFormat()
 *
 * Purpose:  Convert a string (e.g. from command line option arg)
 *           to a format code. Case insensitive. Return
 *           MSAFILE_UNKNOWN/SQFILE_UNKNOWN if string is bad. 
 *           Uses codes defined in squid.h (unaligned formats) and
 *           msa.h (aligned formats).
 *
 * Args:     s   - string to convert; e.g. "stockholm"
 *
 * Returns:  format code; e.g. MSAFILE_STOCKHOLM
 */
int
String2SeqfileFormat(
  char *s);

char *
SeqfileFormat2String(
  int code);

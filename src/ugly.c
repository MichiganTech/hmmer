/* Generated automatically from squid.h.in by configure. */
/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/


#include "squid.h"

 char *aminos      = "ABCDEFGHIKLMNPQRSTVWXYZ*";
 char *primenuc    = "ACGTUN";
 char *protonly    = "EFIPQZ";


void
Warn(char *format, ...)
{
  va_list  argp;
        /* format the error mesg */
  fprintf(stderr, "WARNING: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
}


void
Die(char *format, ...)
{
  va_list  argp;
        /* format the error mesg */
  fprintf(stderr, "\nFATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
        /* exit  */
  exit(1);
}


char *sqd_parse[10];


int
Strparse(char *rexp, char *s, int ntok)
{
  sqd_regexp *pat;
  int         code;
  int         len;
  int         i;
        /* sanity check */
  if (ntok >= NSUBEXP )  Die("Strparse(): ntok must be <= %d", NSUBEXP-1); 

  /* Free previous global substring buffers
   */
  for (i = 0; i <= ntok; i++)
    if (sqd_parse[i] != NULL) 
      { 
  free(sqd_parse[i]);
  sqd_parse[i] = NULL;
      }

  /* Compile and match the pattern, using our modified 
   * copy of Henry Spencer's regexp library
   */
  if ((pat = sqd_regcomp(rexp)) == NULL) 
    Die("regexp compilation failed.");
  code = sqd_regexec(pat, s);

  /* Fill the global substring buffers
   */
  if (code == 1) 
    for (i = 0; i <= ntok; i++)
      if (pat->startp[i] != NULL && pat->endp[i] != NULL)
  {
    len = pat->endp[i] - pat->startp[i];
    sqd_parse[i] = (char *) MallocOrDie(sizeof(char) * (len+1));
    strncpy(sqd_parse[i], pat->startp[i], len);
    sqd_parse[i][len] = '\0';
  }

  free(pat);
  return code;
}



void *
sre_malloc(char *file, int line, size_t size)
{
  void *ptr;

  if ((ptr = malloc (size)) == NULL)
    Die("malloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}


void *
sre_realloc(char *file, int line, void *p, size_t size)
{
  void *ptr;

  if ((ptr = realloc(p, size)) == NULL)
    Die("realloc of %ld bytes failed: file %s line %d", size, file, line);
  return ptr;
}


FILE *
EnvFileOpen(char *fname, char *env, char **ret_dir)
{
  FILE *fp;
  char *path;
  char *s;                      /* ptr to indiv element in env list */
  char  full[1024];             /* constructed file name */

  if (env == NULL) return NULL;
  if ((path = strdup(getenv(env))) == NULL) return NULL;
  
  fp = NULL;
  s  = strtok(path, ":");
  while (s != NULL)
    {
      if (((int) strlen(fname) + (int) strlen(s) + 2) > 1024) 
  { free(path); return NULL; }
      sprintf(full, "%s%c%s", s, '/', fname);
      if ((fp = fopen(full, "r")) != NULL) break;
      s = strtok(NULL, ":");
    }

  /* Return the path we used, if caller wants it
   */
  if (ret_dir != NULL) *ret_dir = strdup(s);
  free(path);
  
  return fp;
}



bool
IsBlankline(char *s)
{
  for (; *s != '\0'; s++)
    if (! isspace(*s)) return false;
  return true;
}


bool
IsInt(char *s)
{
  int hex = 0;

  if (s == NULL) {squid_errno = SQERR_PARAMETER; return false; }

        /* skip whitespace */
  while (isspace((int) (*s))) s++;      
        /* skip leading sign */
  if (*s == '-' || *s == '+') s++;
        /* skip leading conversion signals */
  if ((strncmp(s, "0x", 2) == 0 && (int) strlen(s) > 2) ||
      (strncmp(s, "0X", 2) == 0 && (int) strlen(s) > 2))
    {
      s += 2;
      hex = 1;
    }
  else if (*s == '0' && (int) strlen(s) > 1)
    s++;
        /* examine remainder for garbage chars */
  if (!hex)
    while (*s != '\0')
      {
  if (!isdigit((int) (*s))) return false;
  s++;
      }
  else
    while (*s != '\0')
      {
  if (!isxdigit((int) (*s))) return false;
  s++;
      }

  return true;
}


int
Seqtype(char *seq)
{
  int  saw;     /* how many non-gap characters I saw */
  char c;
  int  po = 0;      /* count of protein-only */
  int  nt = 0;      /* count of t's */
  int  nu = 0;      /* count of u's */
  int  na = 0;      /* count of nucleotides */
  int  aa = 0;      /* count of amino acids */
  int  no = 0;      /* count of others */
  
  /* Look at the first 300 non-gap characters
   */
  for (saw = 0; *seq != '\0' && saw < 300; seq++)
    {
      c = toupper((int) *seq);
      if (! isgap(c)) 
  {
    if (strchr(protonly, c)) po++;
    else if (strchr(primenuc,c)) {
      na++;
      if (c == 'T') nt++;
      else if (c == 'U') nu++;
    }
    else if (strchr(aminos,c)) aa++;
    else if (isalpha((int) c)) no++;
    saw++;
  }
    }

  if (no > 0) return kOtherSeq;
  else if (po > 0) return kAmino;
  else if (na > aa) {
    if (nu > nt) return kRNA;
    else return kDNA;
    }
  else return kAmino;   /* ooooh. risky. */
}

         
int
SeqfileFormat(FILE *fp)
{
  char *buf;
  size_t len;
  int   fmt = SQFILE_UNKNOWN;
  int   ndataline;
  char *bufcpy, *s, *s1, *s2;
  bool   has_junk;

  buf       = NULL;
  len       = 0;
  ndataline = 0;
  has_junk  = false;
  while (getline(&buf, &len, fp))
    {
      if (IsBlankline(buf)) continue;

      /* Well-behaved formats identify themselves in first nonblank line.
       */
      if (ndataline == 0)
  {
    if (strncmp(buf, ">>>>", 4) == 0 && strstr(buf, "Len: "))
      { fmt = SQFILE_GCGDATA; goto DONE; }

    if (buf[0] == '>')
      { fmt = SQFILE_FASTA; goto DONE; }

    if (strncmp(buf, "!!AA_SEQUENCE", 13) == 0 ||
        strncmp(buf, "!!NA_SEQUENCE", 13) == 0)
      { fmt = SQFILE_GCG; goto DONE; }

    if (strncmp(buf, "# STOCKHOLM 1.", 14) == 0)
      { fmt = MSAFILE_STOCKHOLM; goto DONE; }

    if (strncmp(buf, "CLUSTAL", 7) == 0 && 
        strstr(buf, "multiple sequence alignment") != NULL)
      { fmt = MSAFILE_CLUSTAL; goto DONE; }

    if (strncmp(buf, "!!AA_MULTIPLE_ALIGNMENT", 23) == 0 ||
        strncmp(buf, "!!NA_MULTIPLE_ALIGNMENT", 23) == 0)
      { fmt = MSAFILE_MSF; goto DONE; }

          /* PHYLIP id: also just a good bet */
    bufcpy = strdup(buf);
    s = bufcpy;
    if ((s1 = strtok_r(s, WHITESPACE, &s)) != NULL &&
        (s2 = strtok_r(s, WHITESPACE, &s)) != NULL &&
        IsInt(s1) && 
        IsInt(s2))
      { free(bufcpy); fmt = MSAFILE_PHYLIP; goto DONE; }
    free(bufcpy);
  }

      /* We trust that other formats identify themselves soon.
       */
            /* dead giveaways for extended SELEX */
      if (strncmp(buf, "#=AU", 4) == 0 ||
          strncmp(buf, "#=ID", 4) == 0 ||
    strncmp(buf, "#=AC", 4) == 0 ||
    strncmp(buf, "#=DE", 4) == 0 ||
    strncmp(buf, "#=GA", 4) == 0 ||
    strncmp(buf, "#=TC", 4) == 0 ||
    strncmp(buf, "#=NC", 4) == 0 ||
    strncmp(buf, "#=SQ", 4) == 0 ||
    strncmp(buf, "#=SS", 4) == 0 ||
    strncmp(buf, "#=CS", 4) == 0 ||
    strncmp(buf, "#=RF", 4) == 0)
  { fmt = MSAFILE_SELEX; goto DONE; }
  
      if (strncmp(buf, "///", 3) == 0 || strncmp(buf, "ENTRY ", 6) == 0)
  { fmt = SQFILE_PIR; goto DONE; }

        /* a ha, diagnostic of an (old) MSF file */
      if ((strstr(buf, "..")    != NULL) && 
    (strstr(buf, "MSF:")  != NULL) &&
    (strstr(buf, "Check:")!= NULL))
  { fmt = MSAFILE_MSF; goto DONE; }

        /* unaligned GCG (must follow MSF test!) */
      if (strstr(buf, " Check: ") != NULL && strstr(buf, "..") != NULL)
  { fmt = SQFILE_GCG; goto DONE; }

      if (strncmp(buf,"LOCUS ",6) == 0 || strncmp(buf,"ORIGIN ",6) == 0)
  { fmt = SQFILE_GENBANK; goto DONE; }

      if (strncmp(buf,"ID   ",5) == 0 || strncmp(buf,"SQ   ",5) == 0)
  { fmt = SQFILE_EMBL; goto DONE; }

      /* But past here, we're being desperate. A simple SELEX file is
       * very difficult to detect; we can only try to disprove it.
       */
      s = buf;
      if ((s1 = strtok_r(s, WHITESPACE, &s)) == NULL) continue; /* skip blank lines */
      if (strchr("#%", *s1) != NULL) continue;   /* skip comment lines */

      /* Disproof 1. Noncomment, nonblank lines in a SELEX file
       * must have at least two space-delimited fields (name/seq)
       */
      if ((s2 = strtok_r(s, WHITESPACE, &s)) == NULL) 
  has_junk = true;

      /* Disproof 2. 
       * The sequence field should look like a sequence.
       */
      if (s2 != NULL && Seqtype(s2) == kOtherSeq) 
  has_junk = true;

      ndataline++;
      if (ndataline == 300) break; /* only look at first 300 lines */
    }

  if (ndataline == 0)
    Die("Sequence file contains no data");

  /* If we've made it this far, we've run out of data, but there
   * was at least one line of it; check if we've
   * disproven SELEX. If not, cross our fingers, pray, and guess SELEX. 
   */
  if (has_junk == true) fmt = SQFILE_UNKNOWN;
  else                  fmt = MSAFILE_SELEX;

 DONE:
  if (buf != NULL) free(buf);
  rewind(fp);
  return fmt;
}


// int 
// SSIGetFilePosition(FILE *fp, int mode, SSIOFFSET *ret_offset)
// {
//   if (mode == SSI_OFFSET_I32) 
//     {
//       ret_offset->mode    = SSI_OFFSET_I32;
//       ret_offset->off.i32 = ftell(fp);
//       if (ret_offset->off.i32 == (uint32_t)-1) return SSI_ERR_TELL_FAILED;
//     }
//   else if (mode != SSI_OFFSET_I64) abort(); /* only happens on a coding error */
//   else {
//     ret_offset->mode    = SSI_OFFSET_I64;
//     if ((ret_offset->off.i64 = ftello(fp)) == (uint64_t)-1) return SSI_ERR_TELL_FAILED;
//   }
//   return 0;
// }


void 
SeqfileGetLine(SQFILE *V)
{
  if (V->ssimode >= 0) 
    if (0 != SSIGetFilePosition(V->f, V->ssimode, &(V->ssioffset)))
      Die("SSIGetFilePosition() failed");
  if (!getline(&(V->buf), &(V->buflen), V->f))
    *(V->buf) = '\0';
  V->linenumber++;
}


SQFILE*
seqfile_open(char *filename, int format, char *env, int ssimode)
{
  SQFILE *dbfp;

  dbfp = (SQFILE *) MallocOrDie (sizeof(SQFILE));

  dbfp->ssimode  = ssimode;
  dbfp->rpl      = -1;    /* flag meaning "unset" */
  dbfp->lastrpl  = 0;
  dbfp->maxrpl   = 0;
  dbfp->bpl      = -1;    /* flag meaning "unset" */
  dbfp->lastbpl  = 0;
  dbfp->maxbpl   = 0;

  /* Open our file handle.
   * Three possibilities:
   *    1. normal file open
   *    2. filename = "-";    read from stdin
   *    3. filename = "*.gz"; read thru pipe from gzip 
   * If we're reading from stdin or a pipe, we can't reliably
   * back up, so we can't do two-pass parsers like the interleaved alignment   
   * formats.
   */
  if (strcmp(filename, "-") == 0)
    {
      dbfp->f         = stdin;
      dbfp->do_stdin  = true; 
      dbfp->do_gzip   = false;
      dbfp->fname     = strdup("[STDIN]");
    }
  #ifndef SRE_STRICT_ANSI
  /* popen(), pclose() aren't portable to non-POSIX systems; disable */
  else if (Strparse("^.*\\.gz$", filename, 0))
    {
      char cmd[256];

      /* Note that popen() will return "successfully"
       * if file doesn't exist, because gzip works fine
       * and prints an error! So we have to check for
       * existence of file ourself.
       */
      if ( access(filename, F_OK) != -1)
  Die("%s: file does not exist", filename);

      if (strlen(filename) + strlen("gzip -dc ") >= 256)
  Die("filename > 255 char in SeqfileOpen()"); 
      sprintf(cmd, "gzip -dc %s", filename);
      if ((dbfp->f = popen(cmd, "r")) == NULL)
  return NULL;

      dbfp->do_stdin = false;
      dbfp->do_gzip  = true;
      dbfp->fname    = strdup(filename);
    }
  #endif /*SRE_STRICT_ANSI*/
  else
    {
      if ((dbfp->f = fopen(filename, "r")) == NULL &&
    (dbfp->f = EnvFileOpen(filename, env, NULL)) == NULL)
  return NULL;

      dbfp->do_stdin = false;
      dbfp->do_gzip  = false;
      dbfp->fname    = strdup(filename);
    }
  

  /* Invoke autodetection if we haven't already been told what
   * to expect.
   */
  if (format == SQFILE_UNKNOWN)
    {
      if (dbfp->do_stdin == true || dbfp->do_gzip)
  Die("Can't autodetect sequence file format from a stdin or gzip pipe");
      format = SeqfileFormat(dbfp->f);
      if (format == SQFILE_UNKNOWN)
  Die("Can't determine format of sequence file %s", dbfp->fname);
    }

  /* The hack for sequential access of an interleaved alignment file:
   * read the alignment in, we'll copy sequences out one at a time.
   */
  dbfp->msa        = NULL;
  dbfp->afp        = NULL;
  dbfp->format     = format;
  dbfp->linenumber = 0;
  dbfp->buf        = NULL;
  dbfp->buflen     = 0;
  if (IsAlignmentFormat(format))  
    {
      /* We'll be reading from the MSA interface. Copy our data
       * to the MSA afp's structure.
       */
      dbfp->afp           = MallocOrDie(sizeof(MSAFILE));
      dbfp->afp->f        = dbfp->f;            /* just a ptr, don't close */
      dbfp->afp->do_stdin = dbfp->do_stdin;
      dbfp->afp->do_gzip  = dbfp->do_gzip;
      dbfp->afp->fname    = dbfp->fname;        /* just a ptr, don't free */
      dbfp->afp->format   = dbfp->format;       /* e.g. format */
      dbfp->afp->linenumber = dbfp->linenumber; /* e.g. 0 */
      dbfp->afp->buf      = NULL;
      dbfp->afp->buflen   = 0;

      if ((dbfp->msa = MSAFileRead(dbfp->afp)) == NULL)
  Die("Failed to read any alignment data from file %s", dbfp->fname);
        /* hack: overload/reuse msa->lastidx; indicates
           next seq to return upon a ReadSeq() call */
      dbfp->msa->lastidx = 0;

      return dbfp;
    }

  /* Load the first line.
   */
  SeqfileGetLine(dbfp); 
  return dbfp;
}


SQFILE *
SeqfileOpen(char *filename, int format, char *env)
{
  return seqfile_open(filename, format, env, -1);
}


SQFILE *
SeqfileOpenForIndexing(char *filename, int format, char *env, int ssimode)
{
  return seqfile_open(filename, format, env, ssimode);
}


void
SeqfileClose(SQFILE *sqfp)
{
  /* note: don't test for sqfp->msa being NULL. Now that
   * we're holding afp open and allowing access to multi-MSA
   * databases (e.g. Stockholm format, Pfam), msa ends
   * up being NULL when we run out of alignments.
   */
  if (sqfp->afp != NULL) {
    if (sqfp->msa      != NULL) MSAFree(sqfp->msa);
    if (sqfp->afp->buf != NULL) free(sqfp->afp->buf);
    free(sqfp->afp);
  }
#ifndef SRE_STRICT_ANSI /* gunzip functionality only on POSIX systems */
  if (sqfp->do_gzip)         pclose(sqfp->f);
#endif  
  else if (! sqfp->do_stdin) fclose(sqfp->f);
  if (sqfp->buf   != NULL) free(sqfp->buf);
  if (sqfp->fname != NULL) free(sqfp->fname);
  free(sqfp);
}


void
FreeSequence(char *seq, SQINFO *sqinfo)
{
  if (seq != NULL) free(seq);
  if (sqinfo->flags & SQINFO_SS)   free(sqinfo->ss);
  if (sqinfo->flags & SQINFO_SA)   free(sqinfo->sa);
}


int
MakeDealignedString(char *aseq, size_t alen, char *ss, char **ret_s)
{
  char *newStr;
  size_t apos, rpos;

  newStr = (char *) MallocOrDie ((alen+1) * sizeof(char));
  for (apos = rpos = 0; apos < alen; apos++)
    if (! isgap(aseq[apos])){
      newStr[rpos] = ss[apos];
      rpos++;
    }

  newStr[rpos] = '\0';

  if (alen != strlen(ss)){
    squid_errno = SQERR_PARAMETER; 
    free(newStr); 
    return 0; 
  }

  *ret_s = newStr;
  return 1;
}


 void 
addseq(char *s, struct ReadSeqVars *V)
{
  char *s0;
  char *sq;
  int   rpl;      /* valid residues per line */
  int   bpl;      /* characters per line     */

  if (V->ssimode == -1)
    {       /* Normal mode: keeping the seq */
      /* Make sure we have enough room. We know that s is <= buflen,
       * so just make sure we've got room for a whole new buflen worth
       * of sequence.
       */
      if (V->seqlen + V->buflen > V->maxseq) {
  V->maxseq += MAX(V->buflen, kStartLength);
  V->seq = ReallocOrDie (V->seq, V->maxseq+1);
      }

      sq = V->seq + V->seqlen;
      while (*s != 0) {
  if (! isdigit((int) *s) && ! isspace((int) *s)) { 
    *sq = *s;
    sq++;
  }
  s++;
      }
      V->seqlen = sq - V->seq;
    }
  else        /* else: indexing mode, discard the seq */
    {
      s0 = s;
      rpl = 0;
      while (*s != 0) {
  if (! isdigit((int) *s) && ! isspace((int) *s)) {
    rpl++;
  }
  s++;
      }
      V->seqlen += rpl;
      bpl = s - s0;

      /* Keep track of the global rpl, bpl for the file.
       * This is overly complicated because we have to 
       * allow the last line of each record (e.g. the last addseq() call
       * on each sequence) to have a different length - and sometimes
       * we'll have one-line sequence records, too.  Thus we only
       * do something with the global V->rpl when we have *passed over*
       * a line - we keep the last line's rpl in last_rpl. And because
       * a file might consist entirely of single-line records, we keep
       * a third guy, maxrpl, that tells us the maximum rpl of any line
       * in the file. If we reach the end of file and rpl is still unset,
       * we'll set it to maxrpl. If we reach eof and rpl is set, but is
       * less than maxrpl, that's a weird case where a last line in some
       * record is longer than every other line.
       */
      if (V->rpl != 0) {    /* 0 means we already know rpl is invalid       */
  if (V->lastrpl > 0) { /* we're on something that's not the first line */
    if (V->rpl > 0 && V->lastrpl != V->rpl)  V->rpl = 0; 
    else if (V->rpl == -1)                   V->rpl = V->lastrpl;
  }      
  V->lastrpl = rpl;
  if (rpl > V->maxrpl) V->maxrpl = rpl; /* make sure we check max length of final lines */
      }
      if (V->bpl != 0) {    /* 0 means we already know bpl is invalid       */
  if (V->lastbpl > 0) { /* we're on something that's not the first line */
    if (V->bpl > 0 && V->lastbpl != V->bpl)  V->bpl = 0; 
    else if (V->bpl == -1)                   V->bpl = V->lastbpl;
  }      
  V->lastbpl = bpl;
  if (bpl > V->maxbpl) V->maxbpl = bpl; /* make sure we check max length of final lines */
      }
    } /* end of indexing mode of addseq(). */

}


int
SetSeqinfoString(SQINFO *sqinfo, char *sptr, int flag)
{
  int len;
  int pos;

        /* silently ignore NULL. */
  if (sptr == NULL) return 1;

  while (*sptr == ' ') sptr++; /* ignore leading whitespace */
  for (pos = strlen(sptr)-1; pos >= 0; pos--)
    if (! isspace((int) sptr[pos])) break;
  sptr[pos+1] = '\0';        /* ignore trailing whitespace */

  switch (flag) {
  case SQINFO_NAME:
    if (*sptr != '-')
      { 
  strncpy(sqinfo->name, sptr, SQINFO_NAMELEN-1);
  sqinfo->name[SQINFO_NAMELEN-1] = '\0';
  sqinfo->flags   |= SQINFO_NAME;
      }
    break;

  case SQINFO_ID:
    if (*sptr != '-')
      { 
  strncpy(sqinfo->id, sptr, SQINFO_NAMELEN-1);
  sqinfo->id[SQINFO_NAMELEN-1] = '\0';
  sqinfo->flags |= SQINFO_ID;
      }
    break;

  case SQINFO_ACC:
    if (*sptr != '-')
      { 
  strncpy(sqinfo->acc, sptr, SQINFO_NAMELEN-1);
  sqinfo->acc[SQINFO_NAMELEN-1] = '\0';
  sqinfo->flags   |= SQINFO_ACC;
      }
    break;

  case SQINFO_DESC:
    if (*sptr != '-')
      { 
  if (sqinfo->flags & SQINFO_DESC) /* append? */
    {
      len = strlen(sqinfo->desc);
      if (len < SQINFO_DESCLEN-2) /* is there room? */
        {
    strncat(sqinfo->desc, " ", SQINFO_DESCLEN-1-len); len++;
    strncat(sqinfo->desc, sptr, SQINFO_DESCLEN-1-len);
        }
    }
  else      /* else copy */
    strncpy(sqinfo->desc, sptr, SQINFO_DESCLEN-1);
  sqinfo->desc[SQINFO_DESCLEN-1] = '\0';
  sqinfo->flags   |= SQINFO_DESC;
      }
    break;

  case SQINFO_START:
    if (!IsInt(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
    sqinfo->start = atoi(sptr);
    if (sqinfo->start != 0) sqinfo->flags |= SQINFO_START;
    break;

  case SQINFO_STOP:
    if (!IsInt(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
    sqinfo->stop = atoi(sptr);
    if (sqinfo->stop != 0) sqinfo->flags |= SQINFO_STOP;
    break;

  case SQINFO_OLEN:
    if (!IsInt(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
    sqinfo->olen = atoi(sptr);
    if (sqinfo->olen != 0) sqinfo->flags |= SQINFO_OLEN;
    break;

  default:
    Die("Invalid flag %d to SetSeqinfoString()", flag);
  }
  return 1;
}


int
GCGBinaryToSequence(char *seq, int len)
{
  int   bpos;     /* position in binary   */
  int   spos;     /* position in sequence */
  char  twobit;
  int   i;

  for (bpos = (len-1)/4; bpos >= 0; bpos--) 
    {
      twobit = seq[bpos];
      spos   = bpos*4;

      for (i = 3; i >= 0; i--) 
  {
    switch (twobit & 0x3) {
    case 0: seq[spos+i] = 'C'; break;
    case 1: seq[spos+i] = 'T'; break;
    case 2: seq[spos+i] = 'A'; break;
    case 3: seq[spos+i] = 'G'; break;
    }
    twobit = twobit >> 2;
  }
    }
  seq[len] = '\0';
  return 1;
}



 void 
readLoop(int addfirst, int (*endTest)(char *,int *), struct ReadSeqVars *V)
{
  int addend = 0;
  int done   = 0;

  V->seqlen = 0;
  V->lastrpl = V->lastbpl = 0;
  if (addfirst) {
    if (V->ssimode >= 0) V->d_off = V->ssioffset;
    addseq(V->buf, V);
  } else if (V->ssimode >= 0)
    if (0 != SSIGetFilePosition(V->f, V->ssimode, &(V->d_off)))
      Die("SSIGetFilePosition() failed");

  do {
    SeqfileGetLine(V);
  /* feof() alone is a bug; files not necessarily \n terminated */
    if (*(V->buf) == '\0' && feof(V->f))
      done = true;
    done |= (*endTest)(V->buf, &addend);
    if (addend || !done)
      addseq(V->buf, V);
  } while (!done);
}


 int
endPIR(char *s, int  *addend)
{
  *addend = 0;
  if ((strncmp(s, "///", 3) == 0) || 
      (strncmp(s, "ENTRY", 5) == 0))
    return 1;
  else
    return 0;
}


 void
readPIR(struct ReadSeqVars *V)
{
  char *sptr;
        /* load first line of entry  */
  while (!feof(V->f) && strncmp(V->buf, "ENTRY", 5) != 0) {
    SeqfileGetLine(V);
  }
  if (feof(V->f)) return;
  if (V->ssimode >= 0) V->r_off = V->ssioffset;

  if ((sptr = strtok(V->buf + 15, "\n\t ")) != NULL)
    {
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
    }
  do {
    SeqfileGetLine(V);
    if (!feof(V->f) && strncmp(V->buf, "TITLE", 5) == 0)
      SetSeqinfoString(V->sqinfo, V->buf+15, SQINFO_DESC);
    else if (!feof(V->f) && strncmp(V->buf, "ACCESSION", 9) == 0)
      {
  if ((sptr = strtok(V->buf+15, " \t\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
      }
  } while (! feof(V->f) && (strncmp(V->buf,"SEQUENCE", 8) != 0));
  SeqfileGetLine(V);      /* skip next line, coords */

  readLoop(0, endPIR, V);

  /* reading a real PIR-CODATA database file, we keep the source coords
   */
  V->sqinfo->start = 1;
  V->sqinfo->stop  = V->seqlen;
  V->sqinfo->olen  = V->seqlen;
  V->sqinfo->flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;

  /* get next line
   */
  while (!feof(V->f) && strncmp(V->buf, "ENTRY", 5) != 0) {
    SeqfileGetLine(V);
  }
}


 int 
endIG(char *s, int  *addend)
{
  *addend = 1; /* 1 or 2 occur in line w/ bases */
  return((strchr(s,'1')!=NULL) || (strchr(s,'2')!=NULL));
}


 void 
readIG(struct ReadSeqVars *V)
{
  char *nm;
        /* position past ';' comments */
  do {
    SeqfileGetLine(V);
  } while (! (feof(V->f) || ((*V->buf != 0) && (*V->buf != ';')) ));

  if (!feof(V->f))
    {
      if ((nm = strtok(V->buf, "\n\t ")) != NULL)
  SetSeqinfoString(V->sqinfo, nm, SQINFO_NAME);

      readLoop(0, endIG, V);
    }
  
  while (!(feof(V->f) || ((*V->buf != '\0') && (*V->buf == ';'))))
    SeqfileGetLine(V);
}


 int 
endStrider(char *s, int *addend)
{
  *addend = 0;
  return (strstr( s, "//") != NULL);
}


 void 
readStrider(struct ReadSeqVars *V)
{ 
  char *nm;
  
  while ((!feof(V->f)) && (*V->buf == ';')) 
    {
      if (strncmp(V->buf,"; DNA sequence", 14) == 0)
  {
    if ((nm = strtok(V->buf+16, ",\n\t ")) != NULL)
      SetSeqinfoString(V->sqinfo, nm, SQINFO_NAME);
  }
      SeqfileGetLine(V);
    }

  if (! feof(V->f))
    readLoop(1, endStrider, V);

  /* load next line
   */
  while ((!feof(V->f)) && (*V->buf != ';')) 
    SeqfileGetLine(V);
}


 int 
endGB(char *s, int *addend)
{
  *addend = 0;
  return ((strstr(s,"//") != NULL) || (strstr(s,"LOCUS") == s));
}


 void 
readGenBank(struct ReadSeqVars *V)
{
  char *sptr;
  int   in_definition;

  /* We'll map three genbank identifiers onto names:
   *     LOCUS     -> sqinfo.name
   *     ACCESSION -> sqinfo.acc   [primary accession only]
   *     VERSION   -> sqinfo.id
   * We don't currently store the GI number, or secondary accessions.    
   */
  while (strncmp(V->buf, "LOCUS", 5) != 0) {
    SeqfileGetLine(V);
  }
  if (V->ssimode >= 0) V->r_off = V->ssioffset;

  if ((sptr = strtok(V->buf+12, "\n\t ")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);

  in_definition = false;
  while (! feof(V->f))
    {
      SeqfileGetLine(V);
      if (! feof(V->f) && strstr(V->buf, "DEFINITION") == V->buf)
  {
    if ((sptr = strtok(V->buf+12, "\n")) != NULL)
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);
    in_definition = true;
  }
      else if (! feof(V->f) && strstr(V->buf, "ACCESSION") == V->buf)
  {
    if ((sptr = strtok(V->buf+12, "\n\t ")) != NULL)
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
    in_definition = false;
  }
      else if (! feof(V->f) && strstr(V->buf, "VERSION") == V->buf)
  {
    if ((sptr = strtok(V->buf+12, "\n\t ")) != NULL)
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
    in_definition = false;
  }
      else if (strncmp(V->buf,"ORIGIN", 6) != 0)
  {
    if (in_definition)
      SetSeqinfoString(V->sqinfo, V->buf, SQINFO_DESC);
  }
      else
  break;
    }

  readLoop(0, endGB, V);

  /* reading a real GenBank database file, we keep the source coords
   */
  V->sqinfo->start = 1;
  V->sqinfo->stop  = V->seqlen;
  V->sqinfo->olen  = V->seqlen;
  V->sqinfo->flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;


  while (!(feof(V->f) || ((*V->buf!=0) && (strstr(V->buf,"LOCUS") == V->buf))))
    SeqfileGetLine(V);
        /* SRE: V->s now holds "//", so sequential
           reads are wedged: fixed Tue Jul 13 1993 */
  while (!feof(V->f) && strstr(V->buf, "LOCUS  ") != V->buf)
    SeqfileGetLine(V);
}


 int
endGCGdata(char *s, int *addend) 
{
  *addend = 0;
  return (*s == '>');
}


 void
readGCGdata(struct ReadSeqVars *V)
{
  bool   binary = false;   /* whether data are binary or not */
  size_t   blen = 0;   /* length of binary sequence */
  
        /* first line contains ">>>>" followed by name */
  if (Strparse(">>>>([^ ]+) .+2BIT +Len: ([0-9]+)", V->buf, 2))
    {
      binary = true;
      SetSeqinfoString(V->sqinfo, sqd_parse[1], SQINFO_NAME);
      blen = atoi(sqd_parse[2]);
    } 
  else if (Strparse(">>>>([^ ]+) .+ASCII +Len: [0-9]+", V->buf, 1))
    SetSeqinfoString(V->sqinfo, sqd_parse[1], SQINFO_NAME);
  else 
    Die("bogus GCGdata format? %s", V->buf);

        /* second line contains free text description */
  SeqfileGetLine(V);
  SetSeqinfoString(V->sqinfo, V->buf, SQINFO_DESC);

  if (binary) {
    /* allocate for blen characters +3... (allow for 3 bytes of slop) */
    if (blen >= V->maxseq) {
      V->maxseq = blen;
      if ((V->seq = (char *) realloc (V->seq, sizeof(char)*(V->maxseq+4)))==NULL)
  Die("malloc failed");
    }
        /* read (blen+3)/4 bytes from file */
    if (fread(V->seq, sizeof(char), (blen+3)/4, V->f) < (size_t) ((blen+3)/4))
      Die("fread failed");
    V->seqlen = blen;
        /* convert binary code to seq */
    GCGBinaryToSequence(V->seq, blen);
  }
  else readLoop(0, endGCGdata, V);
  
  while (!(feof(V->f) || ((*V->buf != 0) && (*V->buf == '>'))))
    SeqfileGetLine(V);
}


 int
endPearson(char *s, int *addend)
{
  *addend = 0;
  return(*s == '>');
}


 void 
readPearson(struct ReadSeqVars *V)
{
  char *sptr;

  if (V->ssimode >= 0) V->r_off = V->ssioffset;

  if (*V->buf != '>') 
    Die("\
File %s does not appear to be in FASTA format at line %d.\n\
You may want to specify the file format on the command line.\n\
Usually this is done with an option --informat <fmt>.\n", 
  V->fname, V->linenumber);

  if ((sptr = strtok(V->buf+1, "\n\t ")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
  if ((sptr = strtok(NULL, "\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);

  readLoop(0, endPearson, V);

  while (!(feof(V->f) || ((*V->buf != 0) && (*V->buf == '>')))) {
    SeqfileGetLine(V);
  }
}


 int
endEMBL(char *s, int *addend)
{
  *addend = 0;
  /* Some people (Berlin 5S rRNA database, f'r instance) use
   * an extended EMBL format that attaches extra data after
   * the sequence -- watch out for that. We use the fact that
   * real EMBL sequence lines begin with five spaces.
   * 
   * We can use this as the sole end test because readEMBL() will
   * advance to the next ID line before starting to read again.
   */
  return (strncmp(s,"     ",5) != 0);
/*  return ((strstr(s,"//") != NULL) || (strstr(s,"ID   ") == s)); */
}

 void 
readEMBL(struct ReadSeqVars *V)
{
  char *sptr;

        /* make sure we have first line */
  while (!feof(V->f) && strncmp(V->buf, "ID  ", 4) != 0) {
    SeqfileGetLine(V);
  }
  if (V->ssimode >= 0) V->r_off = V->ssioffset;

  if ((sptr = strtok(V->buf+5, "\n\t ")) != NULL)
    {
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
    }

  do {
    SeqfileGetLine(V);
    if (!feof(V->f) && strstr(V->buf, "AC  ") == V->buf)
      {
  if ((sptr = strtok(V->buf+5, ";  \t\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
      }
    else if (!feof(V->f) && strstr(V->buf, "DE  ") == V->buf)
      {
  if ((sptr = strtok(V->buf+5, "\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);
      }
  } while (! feof(V->f) && strncmp(V->buf,"SQ",2) != 0);
  
  readLoop(0, endEMBL, V);

  /* Hack for Staden experiment files: convert - to N
   */
  if (V->ssimode == -1)   /* if we're in ssi mode, we're not keeping the seq */
    for (sptr = V->seq; *sptr != '\0'; sptr++)
      if (*sptr == '-') *sptr = 'N';

  /* reading a real EMBL database file, we keep the source coords
   */
  V->sqinfo->start = 1;
  V->sqinfo->stop  = V->seqlen;
  V->sqinfo->olen  = V->seqlen;
  V->sqinfo->flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;

        /* load next record's ID line */
  while (!feof(V->f) && strncmp(V->buf, "ID  ", 4) != 0) {
    SeqfileGetLine(V);
  }    

}


 int
endZuker(char *s, int *addend)
{
  *addend = 0;
  return( *s == '(' );
}


 void
readZuker(struct ReadSeqVars *V)
{
  char *sptr;

  SeqfileGetLine(V);  /*s == "seqLen seqid string..."*/

  if ((sptr = strtok(V->buf+6, " \t\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);

  if ((sptr = strtok(NULL, "\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);

  readLoop(0, endZuker, V);

  while (!(feof(V->f) | ((*V->buf != '\0') & (*V->buf == '('))))
    SeqfileGetLine(V);
}

 void 
readUWGCG(struct ReadSeqVars *V)
{
  char  *si;
  char  *sptr;
  int    done;

  V->seqlen = 0;

  if ((si = strstr(V->buf,"  Length: ")) != NULL) *si = 0;
  else if ((si = strstr(V->buf,"..")) != NULL)    *si = 0;

  if ((sptr = strtok(V->buf, "\n\t ")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);

  do {
    done = feof(V->f);
    SeqfileGetLine(V);
    if (! done) addseq(V->buf, V);
  } while (!done);
}

int
ReadSeq(SQFILE *V, char **ret_seq, SQINFO *sqinfo)
{
  int    gotuw;

  squid_errno = SQERR_OK;

  /* Here's the hack for sequential access of sequences from
   * the multiple sequence alignment formats
   */
  if (IsAlignmentFormat(V->format))
    {
      if (V->msa->lastidx >= V->msa->nseq) 
  { /* out of data. try to read another alignment */        
    MSAFree(V->msa);
    if ((V->msa = MSAFileRead(V->afp)) == NULL)
      return 0;
    V->msa->lastidx = 0;
  }
        /* copy and dealign the appropriate aligned seq */
      MakeDealignedString(V->msa->aseq[V->msa->lastidx], V->msa->alen, 
        V->msa->aseq[V->msa->lastidx], &(V->seq));
      V->seqlen = strlen(V->seq);

      /* Extract sqinfo stuff for this sequence from the msa.
       * Tedious; code that should be cleaned.
       */
      sqinfo->flags = 0;
      if (V->msa->sqname[V->msa->lastidx] != NULL) 
  SetSeqinfoString(sqinfo, V->msa->sqname[V->msa->lastidx], SQINFO_NAME);
      if (V->msa->sqacc != NULL && V->msa->sqacc[V->msa->lastidx] != NULL) 
  SetSeqinfoString(sqinfo, V->msa->sqacc[V->msa->lastidx], SQINFO_ACC);
      if (V->msa->sqdesc != NULL && V->msa->sqdesc[V->msa->lastidx] != NULL) 
  SetSeqinfoString(sqinfo, V->msa->sqdesc[V->msa->lastidx], SQINFO_DESC);
      if (V->msa->ss != NULL && V->msa->ss[V->msa->lastidx] != NULL) {
  MakeDealignedString(V->msa->aseq[V->msa->lastidx], V->msa->alen, 
          V->msa->ss[V->msa->lastidx], &(sqinfo->ss));
  sqinfo->flags |= SQINFO_SS;
      }
      if (V->msa->sa != NULL && V->msa->sa[V->msa->lastidx] != NULL) {
  MakeDealignedString(V->msa->aseq[V->msa->lastidx], V->msa->alen, 
          V->msa->sa[V->msa->lastidx], &(sqinfo->sa));
  sqinfo->flags |= SQINFO_SA;
      }
      V->msa->lastidx++;
    } 
  else {
    if (feof(V->f)) return 0;

    if (V->ssimode == -1) { /* normal mode */
      V->seq           = (char*) calloc (kStartLength+1, sizeof(char));
      V->maxseq        = kStartLength;
    } else {      /* index mode: discarding seq */
      V->seq           = NULL;
      V->maxseq        = 0;
    }
    V->seqlen        = 0;
    V->sqinfo        = sqinfo;
    V->sqinfo->flags = 0;

    switch (V->format) {
    case SQFILE_IG      : readIG(V);      break;
    case SQFILE_STRIDER : readStrider(V); break;
    case SQFILE_GENBANK : readGenBank(V); break;
    case SQFILE_FASTA   : readPearson(V); break;
    case SQFILE_EMBL    : readEMBL(V);    break;
    case SQFILE_ZUKER   : readZuker(V);   break;
    case SQFILE_PIR     : readPIR(V);     break;
    case SQFILE_GCGDATA : readGCGdata(V); break; 
  
    case SQFILE_GCG :
      do {      /* skip leading comments on GCG file */
  gotuw = (strstr(V->buf,"..") != NULL);
  if (gotuw) readUWGCG(V);
  SeqfileGetLine(V);
      } while (! feof(V->f));
      break;

    case SQFILE_IDRAW:   /* SRE: no attempt to read idraw postscript */
    default:
      squid_errno = SQERR_FORMAT;
      free(V->seq);
      return 0;
    }
    if (V->seq != NULL)   /* (it can be NULL in indexing mode) */
      V->seq[V->seqlen] = 0; /* stick a string terminator on it */
  }

  /* Cleanup
   */
  sqinfo->len    = V->seqlen; 
  sqinfo->flags |= SQINFO_LEN;
  *ret_seq = V->seq;
  if (squid_errno == SQERR_OK) return 1; else return 0;
}  




int
ReadMultipleRseqs(char              *seqfile,
      int                fformat,
      char            ***ret_rseqs,
      SQINFO **ret_sqinfo,
      int               *ret_num)
{
  SQINFO *sqinfo;               /* array of sequence optional info         */
  SQFILE *dbfp;                 /* open ptr for sequential access of file  */
  char  **rseqs;                /* sequence array                          */
  int     numalloced;           /* num of seqs currently alloced for       */
  int     num;


  num        = 0;
  numalloced = 16;
  rseqs  = (char **) MallocOrDie (numalloced * sizeof(char *));
  sqinfo = (SQINFO *) MallocOrDie (numalloced * sizeof(SQINFO));
  if ((dbfp = SeqfileOpen(seqfile, fformat, NULL)) == NULL) return 0;      

  while (ReadSeq(dbfp, &rseqs[num], &(sqinfo[num])))
    {
      num++;
      if (num == numalloced) /* more seqs coming, alloc more room */
  {
    numalloced += 16;
    rseqs  = (char **) ReallocOrDie (rseqs, numalloced*sizeof(char *));
    sqinfo = (SQINFO *) ReallocOrDie (sqinfo, numalloced * sizeof(SQINFO));
  }
    }
  SeqfileClose(dbfp);

  *ret_rseqs  = rseqs;
  *ret_sqinfo = sqinfo;
  *ret_num    = num;
  return 1;
}




double
Gammln(
  double xx
){
  double x,tmp,ser;
  const double cof[6]={
     76.18009173,   -86.50532033,  24.01409822,
    -1.231739516, 0.120858003e-2, -0.536382e-5};

  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (int j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}




char *
FileConcat(char *dir, char *file)
{
  char *full;

  full = (char *) MallocOrDie (sizeof(char) * (strlen(dir)+strlen(file)+2));

  if (*file == '/') 
    strcpy(full, file); /* file = "/foo", ignore directory. */
  else
    printf(full, "%s%c%s", dir, '/', file);

  return full;
}




void
s2upper(char *s)
{
  for (; *s != '\0'; s++)
    *s = (char) toupper(*s);
}



void
StringChop(char *s)
{
  int   i;

  i = strlen(s) - 1;             /* set i at last char in string     */
  while (i >= 0 && isspace((int) s[i])) i--;   /* i now at last non-whitespace char, or -1 */
  s[i+1] = '\0';
}


int
MakeAlignedString(char *aseq, size_t alen, char *ss, char **ret_s)
{
  char *new; 
  size_t   apos, rpos;

  new = (char *) MallocOrDie ((alen+1) * sizeof(char));
  for (apos = rpos = 0; apos < alen; apos++)
    if (! isgap(aseq[apos]))
      {
  new[apos] = ss[rpos];
  rpos++;
      }
    else
      new[apos] = '.';
  new[apos] = '\0';

  if (rpos != strlen(ss))
    { squid_errno = SQERR_PARAMETER; free(new); return 0; }
  *ret_s = new;
  return 1;
}


int
ParsePAMFile(FILE *fp, int ***ret_pam, float *ret_scale)
{
  int    **pam;
  char     buffer[512];   /* input buffer from fp                  */
  int      order[27];   /* order of fields, obtained from header */
  int      nsymbols;    /* total number of symbols in matrix     */
  char    *sptr;
  int      idx;
  int      row, col;
  float    scale;
  int      gotscale = false;
  
  scale = 0.0;    /* just to silence gcc uninit warnings */
  if (fp == NULL) { squid_errno = SQERR_NODATA; return 0; }
  
  /* Look at the first non-blank, non-comment line in the file.
   * It gives single-letter codes in the order the PAM matrix
   * is arrayed in the file. 
   */
  do {
    if (fgets(buffer, 512, fp) == NULL) 
      { squid_errno = SQERR_NODATA; return 0; }

    /* Get the scale factor from the header.
     * For BLOSUM files, we assume the line looks like:
     *     BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
     * and we assume that the fraction is always 1/x;
     * 
     * For PAM files, we assume the line looks like:
     *     PAM 120 substitution matrix, scale = ln(2)/2 = 0.346574
     * and we assume that the number following the final '=' is our scale
     */
    if (strstr(buffer, "BLOSUM Clustered Scoring Matrix") != NULL &&
  (sptr = strchr(buffer, '/')) != NULL)
      {
  sptr++;
  if (! isdigit((int) (*sptr))) { squid_errno = SQERR_FORMAT; return 0; }
  scale = (float) (log(2.0) / atof(sptr));
  gotscale = true;
      }
    else if (strstr(buffer, "substitution matrix,") != NULL)
      {
  while ((sptr = strrchr(buffer, '=')) != NULL) {
    sptr += 2;
    float tmpFloat;
    if (1 == sscanf(sptr, "%f", &tmpFloat)) {
      scale = atof(sptr);
      gotscale = true;
      break;
    }
  }
      }
  } while ((sptr = strtok(buffer, " \t\n")) == NULL || *sptr == '#');

  idx = 0;
  do {
    order[idx] = (int) *sptr - (int) 'A';
    if (order[idx] < 0 || order[idx] > 25) order[idx] = 26;
    idx++;
  } while ((sptr = strtok(NULL, " \t\n")) != NULL);
  nsymbols = idx;
  
  /* Allocate a pam matrix. For speed of indexing, we use
   * a 27x27 matrix so we can do lookups using the ASCII codes
   * of amino acid single-letter representations, plus one
   * extra field to deal with the "*" (terminators).
   */
  if ((pam = (int **) calloc (27, sizeof(int *))) == NULL)
    Die("calloc failed");
  for (idx = 0; idx < 27; idx++)
    if ((pam[idx] = (int *) calloc (27, sizeof(int))) == NULL)
      Die("calloc failed");

  /* Parse the rest of the file.
   */
  for (row = 0; row < nsymbols; row++)
    {
      if (fgets(buffer, 512, fp) == NULL) 
  { squid_errno = SQERR_NODATA; return 0; }

      if ((sptr = strtok(buffer, " \t\n")) == NULL)
  { squid_errno = SQERR_NODATA; return 0; }
      for (col = 0; col < nsymbols; col++)
  {
    if (sptr == NULL) { squid_errno = SQERR_NODATA; return 0; }

    /* Watch out for new BLAST format, with leading characters
     */
    if (*sptr == '*' || isalpha((int) *sptr))
      col--;  /* hack hack */
    else
      pam [order[row]] [order[col]] = atoi(sptr);

    sptr = strtok(NULL, " \t\n");
  }
    }
  
  /* Return
   */
  if (ret_scale != NULL)
    {
      if (gotscale) *ret_scale = scale;
      else
  {
    Warn("Failed to parse PAM matrix scale factor. Defaulting to ln(2)/2!");
    *ret_scale = log(2.0) / 2.0;
  }
    }
  *ret_pam = pam;
  return 1;
}


void
SqdClean(
){
  int i;

  /* Free global substring buffers that Strparse() uses
   */
  for (i = 0; i <= 9; i++)
    if (sqd_parse[i] != NULL) {
      free(sqd_parse[i]);
      sqd_parse[i] = NULL;
    }
}
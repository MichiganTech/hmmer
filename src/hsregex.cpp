/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/*****************************************************************
 * This code is an altered version of Henry Spencer's
 * regex library. Alterations are limited to minor streamlining,
 * and some name changes to protect the SQUID namespace.
 * Henry's copyright notice appears below.
 * You can obtain the original from
 *    ftp://ftp.zoo.toronto.edu/pub/bookregex.tar.Z
 * Thanks, Henry!
 *
 * The magic word for compiling a testdriver: NBA_TEAM_IN_STL
 * gcc -o test -g -DNBA_TEAM_IN_STL -L. hsregex.c -lsquid -lm
 * 
 * Usage:
 *  test <pattern> <ntok> <string>
 * 
 *****************************************************************/   

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "hsregex.hpp"


/* all code below is:
 * Copyright (c) 1986, 1993, 1995 by University of Toronto.
 * Written by Henry Spencer.  Not derived from licensed software.
 *
 * Permission is granted to anyone to use this software for any
 * purpose on any computer system, and to redistribute it in any way,
 * subject to the following restrictions:
 *
 * 1. The author is not responsible for the consequences of use of
 *   this software, no matter how awful, even if they arise
 *   from defects in it.
 *
 * 2. The origin of this software must not be misrepresented, either
 *   by explicit claim or by omission.
 *
 * 3. Altered versions must be plainly marked as such, and must not
 *   be misrepresented (by explicit claim or omission) as being
 *  the original software.
 *
 * 4. This notice must not be removed or altered.
 */
// This file has been modified to fall in line with newer coding standards.

/*
 * sqd_regcomp and sqd_regexec -- sqd_regsub and sqd_regerror are elsewhere
 */

sqd_regexp*
sqd_regcomp(
  const char *exp
){
  sqd_regexp *r;
  char *scan;
  int flags;
  struct comp co;

  if (exp == NULL)
    FAIL("NULL argument to sqd_regcomp");

  /* First pass: determine size, legality. */
  co.regparse = (char *)exp;
  co.regnpar = 1;
  co.regsize = 0L;
  co.regdummy[0] = NOTHING;
  co.regdummy[1] = co.regdummy[2] = 0;
  co.regcode = co.regdummy;
  regc(&co, SQD_REGMAGIC);
  if (reg(&co, 0, &flags) == NULL)
    return(NULL);

  /* Small enough for pointer-storage convention? */
  if (co.regsize >= 0x7fffL)  /* Probably could be 0xffffL. */
    FAIL("regexp too big");

  /* Allocate space. */
  r = (sqd_regexp *)malloc(sizeof(sqd_regexp) + (size_t)co.regsize);
  if (r == NULL)
    FAIL("out of space");

  /* Second pass: emit code. */
  co.regparse = (char *)exp;
  co.regnpar = 1;
  co.regcode = r->program;
  regc(&co, SQD_REGMAGIC);
  if (reg(&co, 0, &flags) == NULL)
    return(NULL);

  /* Dig out information for optimizations. */
  r->regstart = '\0';    /* Worst-case defaults. */
  r->reganch = 0;
  r->regmust = NULL;
  r->regmlen = 0;
  scan = r->program+1;    /* First BRANCH. */
  if (OP(regnext(scan)) == END) {  /* Only one top-level choice. */
    scan = OPERAND(scan);

    /* Starting-point info. */
    if (OP(scan) == EXACTLY)
      r->regstart = *OPERAND(scan);
    else if (OP(scan) == BOL)
      r->reganch = 1;

    /*
     * If there's something expensive in the r.e., find the
     * longest literal string that must appear and make it the
     * regmust.  Resolve ties in favor of later strings, since
     * the regstart check works with the beginning of the r.e.
     * and avoiding duplication strengthens checking.  Not a
     * strong reason, but sufficient in the absence of others.
     */
    if (flags&SPSTART) {
      char *longest = NULL;
      size_t len = 0;

      for (; scan != NULL; scan = regnext(scan))
        if (OP(scan) == EXACTLY && strlen(OPERAND(scan)) >= len) {
          longest = OPERAND(scan);
          len = strlen(OPERAND(scan));
        }
      r->regmust = longest;
      r->regmlen = (int)len;
    }
  }

  return(r);
}


char*
reg(
  struct comp *cp,
  bool paren,
  int *flagp
){
  char *ret = NULL;   /* SRE: NULL init added to silence gcc */
  char *br;
  char *ender;
  int parno = 0;  /* SRE: init added to silence gcc */
  int flags;

  *flagp = HASWIDTH;  /* Tentatively. */

  if (paren) {
    /* Make an OPEN node. */
    if (cp->regnpar >= NSUBEXP)
      FAIL("too many ()");
    parno = cp->regnpar;
    cp->regnpar++;
    ret = regnode(cp, OPEN+parno);
  }

  /* Pick up the branches, linking them together. */
  br = regbranch(cp, &flags);
  if (br == NULL)
    return(NULL);
  if (paren)
    regtail(cp, ret, br);  /* OPEN -> first. */
  else
    ret = br;
  *flagp &= ~(~flags&HASWIDTH);  /* Clear bit if bit 0. */
  *flagp |= flags&SPSTART;
  while (*cp->regparse == '|') {
    cp->regparse++;
    br = regbranch(cp, &flags);
    if (br == NULL)
      return(NULL);
    regtail(cp, ret, br);  /* BRANCH -> BRANCH. */
    *flagp &= ~(~flags&HASWIDTH);
    *flagp |= flags&SPSTART;
  }

  /* Make a closing node, and hook it on the end. */
  ender = regnode(cp, (paren) ? CLOSE+parno : END);
  regtail(cp, ret, ender);

  /* Hook the tails of the branches to the closing node. */
  for (br = ret; br != NULL; br = regnext(br))
    regoptail(cp, br, ender);

  /* Check for proper termination. */
  if (paren && *cp->regparse++ != ')') {
    FAIL("unterminated ()");
  } else if (!paren && *cp->regparse != '\0') {
    if (*cp->regparse == ')') {
      FAIL("unmatched ()");
    } else
      FAIL("internal error: junk on end");
    /* NOTREACHED */
  }

  return(ret);
}


char*
regbranch(
  struct comp *cp,
  int *flagp
){
  char *ret;
  char *chain;
  char *latest;
  int flags = 0;
  int c;

  *flagp = WORST;        /* Tentatively. */

  ret = regnode(cp, BRANCH);
  chain = NULL;
  while ((c = *cp->regparse) != '\0' && c != '|' && c != ')') {
    latest = regpiece(cp, &flags);
    if (latest == NULL)
      return(NULL);
    *flagp |= flags&HASWIDTH;
    if (chain == NULL)    /* First piece. */
      *flagp |= flags&SPSTART;
    else
      regtail(cp, chain, latest);
    chain = latest;
  }
  if (chain == NULL)      /* Loop ran zero times. */
    (void) regnode(cp, NOTHING);

  return(ret);
}


char*
regpiece(
  struct comp *cp,
  int *flagp
){
  char *ret;
  char op;
  char *next;
  int flags;

  ret = regatom(cp, &flags);
  if (ret == NULL)
    return(NULL);

  op = *cp->regparse;
  if (!ISREPN(op)) {
    *flagp = flags;
    return(ret);
  }

  if (!(flags&HASWIDTH) && op != '?')
    FAIL("*+ operand could be empty");
  switch (op) {
  case '*':  *flagp = WORST|SPSTART;      break;
  case '+':  *flagp = WORST|SPSTART|HASWIDTH;  break;
  case '?':  *flagp = WORST;        break;
  }

  if (op == '*' && (flags&SIMPLE))
    reginsert(cp, STAR, ret);
  else if (op == '*') {
    /* Emit x* as (x&|), where & means "self". */
    reginsert(cp, BRANCH, ret);    /* Either x */
    regoptail(cp, ret, regnode(cp, BACK));  /* and loop */
    regoptail(cp, ret, ret);    /* back */
    regtail(cp, ret, regnode(cp, BRANCH));  /* or */
    regtail(cp, ret, regnode(cp, NOTHING));  /* null. */
  } else if (op == '+' && (flags&SIMPLE))
    reginsert(cp, PLUS, ret);
  else if (op == '+') {
    /* Emit x+ as x(&|), where & means "self". */
    next = regnode(cp, BRANCH);    /* Either */
    regtail(cp, ret, next);
    regtail(cp, regnode(cp, BACK), ret);  /* loop back */
    regtail(cp, next, regnode(cp, BRANCH));  /* or */
    regtail(cp, ret, regnode(cp, NOTHING));  /* null. */
  } else if (op == '?') {
    /* Emit x? as (x|) */
    reginsert(cp, BRANCH, ret);    /* Either x */
    regtail(cp, ret, regnode(cp, BRANCH));  /* or */
    next = regnode(cp, NOTHING);    /* null. */
    regtail(cp, ret, next);
    regoptail(cp, ret, next);
  }
  cp->regparse++;
  if (ISREPN(*cp->regparse))
    FAIL("nested *?+");

  return(ret);
}


char *
regatom(
  struct comp *cp,
  int *flagp
){
  char *ret;
  int flags;

  *flagp = WORST;    /* Tentatively. */

  switch (*cp->regparse++) {
  case '^':
    ret = regnode(cp, BOL);
    break;
  case '$':
    ret = regnode(cp, EOL);
    break;
  case '.':
    ret = regnode(cp, ANY);
    *flagp |= HASWIDTH|SIMPLE;
    break;
  case '[': {
    int range;
    int rangeend;
    int c;

    if (*cp->regparse == '^') {  /* Complement of range. */
      ret = regnode(cp, ANYBUT);
      cp->regparse++;
    } else
      ret = regnode(cp, ANYOF);
    if ((c = *cp->regparse) == ']' || c == '-') {
      regc(cp, c);
      cp->regparse++;
    }
    while ((c = *cp->regparse++) != '\0' && c != ']') {
      if (c != '-')
        regc(cp, c);
      else if ((c = *cp->regparse) == ']' || c == '\0')
        regc(cp, '-');
      else {
        range = (unsigned char)*(cp->regparse-2);
        rangeend = (unsigned char)c;
        if (range > rangeend)
          FAIL("invalid [] range");
        for (range++; range <= rangeend; range++)
          regc(cp, range);
        cp->regparse++;
      }
    }
    regc(cp, '\0');
    if (c != ']')
      FAIL("unmatched []");
    *flagp |= HASWIDTH|SIMPLE;
    break;
    }
  case '(':
    ret = reg(cp, 1, &flags);
    if (ret == NULL)
      return(NULL);
    *flagp |= flags&(HASWIDTH|SPSTART);
    break;
  case '\0':
  case '|':
  case ')':
    /* supposed to be caught earlier */
    FAIL("internal error: \\0|) unexpected");
    break;
  case '?':
  case '+':
  case '*':
    FAIL("?+* follows nothing");
    break;
  case '\\':
    if (*cp->regparse == '\0')
      FAIL("trailing \\");
    ret = regnode(cp, EXACTLY);
    regc(cp, *cp->regparse++);
    regc(cp, '\0');
    *flagp |= HASWIDTH|SIMPLE;
    break;
  default: {
    size_t len;
    char ender;

    cp->regparse--;
    len = strcspn(cp->regparse, META);
    if (len == 0)
      FAIL("internal error: strcspn 0");
    ender = *(cp->regparse+len);
    if (len > 1 && ISREPN(ender))
      len--;    /* Back off clear of ?+* operand. */
    *flagp |= HASWIDTH;
    if (len == 1)
      *flagp |= SIMPLE;
    ret = regnode(cp, EXACTLY);
    for (; len > 0; len--)
      regc(cp, *cp->regparse++);
    regc(cp, '\0');
    break;
    }
  }

  return(ret);
}


char*      /* Location. */
regnode(
  struct comp *cp,
  char op
){
  char *const ret = cp->regcode;
  char *ptr;

  if (!EMITTING(cp)) {
    cp->regsize += 3;
    return(ret);
  }

  ptr = ret;
  *ptr++ = op;
  *ptr++ = '\0';    /* Null next pointer. */
  *ptr++ = '\0';
  cp->regcode = ptr;

  return(ret);
}


void
regc(
  struct comp *cp,
  char b
){
  if (EMITTING(cp))
    *cp->regcode++ = b;
  else
    cp->regsize++;
}


void
reginsert(
  struct comp *cp,
  char op,
  char *opnd
){
  char *place;

  if (!EMITTING(cp)) {
    cp->regsize += 3;
    return;
  }

  (void) memmove(opnd+3, opnd, (size_t)(cp->regcode - opnd));
  cp->regcode += 3;

  place = opnd;    /* Op node, where operand used to be. */
  *place++ = op;
  *place++ = '\0';
  *place++ = '\0';
}


void
regtail(
  struct comp *cp,
  char *p,
  char *val
){
  char *scan;
  char *temp;
  int offset;

  if (!EMITTING(cp))
    return;

  /* Find last node. */
  for (scan = p; (temp = regnext(scan)) != NULL; scan = temp)
    continue;

  offset = (OP(scan) == BACK) ? scan - val : val - scan;
  *(scan+1) = (offset>>8)&0177;
  *(scan+2) = offset&0377;
}


void
regoptail(
  struct comp *cp,
  char *p,
  char *val
){
  /* "Operandless" and "op != BRANCH" are synonymous in practice. */
  if (!EMITTING(cp) || OP(p) != BRANCH)
    return;
  regtail(cp, OPERAND(p), val);
}


struct exec {
  char *reginput;    /* String-input pointer. */
  char *regbol;    /* Beginning of input, for ^ check. */
  char **regstartp;  /* Pointer to startp array. */
  char **regendp;    /* Ditto for endp. */
};


#ifdef DEBUG
int regnarrate = 0;
#endif


int
sqd_regexec(
  sqd_regexp *prog,
  const char *str
){
  char *string = (char *)str;  /* avert const poisoning */
  char *s;
  struct exec ex;

  /* Be paranoid. */
  if (prog == NULL || string == NULL) {
    sqd_regerror("NULL argument to sqd_regexec");
    return(0);
  }

  /* Check validity of program. */
  if ((unsigned char)*prog->program != SQD_REGMAGIC) {
    sqd_regerror("corrupted regexp");
    return(0);
  }

  /* If there is a "must appear" string, look for it. */
  if (prog->regmust != NULL && strstr(string, prog->regmust) == NULL)
    return(0);

  /* Mark beginning of line for ^ . */
  ex.regbol = string;
  ex.regstartp = prog->startp;
  ex.regendp = prog->endp;

  /* Simplest case:  anchored match need be tried only once. */
  if (prog->reganch)
    return(regtry(&ex, prog, string));

  /* Messy cases:  unanchored match. */
  if (prog->regstart != '\0') {
    /* We know what char it must start with. */
    for (s = string; s != NULL; s = strchr(s+1, prog->regstart))
      if (regtry(&ex, prog, s))
        return(1);
    return(0);
  } else {
    /* We don't -- general case. */
    for (s = string; !regtry(&ex, prog, s); s++)
      if (*s == '\0')
        return(0);
    return(1);
  }
  /* NOTREACHED */
}


bool      /* 0 failure, 1 success */
regtry(
  struct exec *ep,
  sqd_regexp *prog,
  char *string
){
  int i;
  char **stp;
  char **enp;

  ep->reginput = string;

  stp = prog->startp;
  enp = prog->endp;
  for (i = NSUBEXP; i > 0; i--) {
    *stp++ = NULL;
    *enp++ = NULL;
  }
  if (regmatch(ep, prog->program + 1)) {
    prog->startp[0] = string;
    prog->endp[0] = ep->reginput;
    return true;
  } else
    return false;
}


bool      /* 0 failure, 1 success */
regmatch(
  struct exec *ep,
  char *prog
){
  char *scan;  /* Current node. */
  char *next;    /* Next node. */

#ifdef DEBUG
  if (prog != NULL && regnarrate)
    fprintf(stderr, "%s(\n", regprop(prog));
#endif
  for (scan = prog; scan != NULL; scan = next) {
#ifdef DEBUG
    if (regnarrate)
      fprintf(stderr, "%s...\n", regprop(scan));
#endif
    next = regnext(scan);

    switch (OP(scan)) {
    case BOL:
      if (ep->reginput != ep->regbol)
        return false;
      break;
    case EOL:
      if (*ep->reginput != '\0')
        return false;
      break;
    case ANY:
      if (*ep->reginput == '\0')
        return false;
      ep->reginput++;
      break;
    case EXACTLY: {
      size_t len;
      char *const opnd = OPERAND(scan);

      /* Inline the first character, for speed. */
      if (*opnd != *ep->reginput)
        return false;
      len = strlen(opnd);
      if (len > 1 && strncmp(opnd, ep->reginput, len) != 0)
        return false;
      ep->reginput += len;
      break;
      }
    case ANYOF:
      if (*ep->reginput == '\0' ||
          strchr(OPERAND(scan), *ep->reginput) == NULL)
        return false;
      ep->reginput++;
      break;
    case ANYBUT:
      if (*ep->reginput == '\0' ||
          strchr(OPERAND(scan), *ep->reginput) != NULL)
        return false;
      ep->reginput++;
      break;
    case NOTHING:
      break;
    case BACK:
      break;
    case OPEN+1: case OPEN+2: case OPEN+3:
    case OPEN+4: case OPEN+5: case OPEN+6:
    case OPEN+7: case OPEN+8: case OPEN+9: {
      const int no = OP(scan) - OPEN;
      char *const input = ep->reginput;

      if (regmatch(ep, next)) {
        /*
         * Don't set startp if some later
         * invocation of the same parentheses
         * already has.
         */
        if (ep->regstartp[no] == NULL)
          ep->regstartp[no] = input;
        return true;
      } else
        return false;
      break;
      }
    case CLOSE+1: case CLOSE+2: case CLOSE+3:
    case CLOSE+4: case CLOSE+5: case CLOSE+6:
    case CLOSE+7: case CLOSE+8: case CLOSE+9: {
      const int no = OP(scan) - CLOSE;
      char *const input = ep->reginput;

      if (regmatch(ep, next)) {
        /*
         * Don't set endp if some later
         * invocation of the same parentheses
         * already has.
         */
        if (ep->regendp[no] == NULL)
          ep->regendp[no] = input;
        return true;
      } else
        return false;
      break;
      }
    case BRANCH: {
      char *const save = ep->reginput;

      if (OP(next) != BRANCH)    /* No choice. */
        next = OPERAND(scan);  /* Avoid recursion. */
      else {
        while (OP(scan) == BRANCH) {
          if (regmatch(ep, OPERAND(scan)))
            return true;
          ep->reginput = save;
          scan = regnext(scan);
        }
        return false;
        /* NOTREACHED */
      }
      break;
      }
    case STAR: case PLUS: {
      const char nextch =
        (OP(next) == EXACTLY) ? *OPERAND(next) : '\0';
      size_t no;
      char *const save = ep->reginput;
      const size_t min = (OP(scan) == STAR) ? 0 : 1;

      for (no = regrepeat(ep, OPERAND(scan)) + 1; no > min; no--) {
        ep->reginput = save + no - 1;
        /* If it could work, try it. */
        if (nextch == '\0' || *ep->reginput == nextch)
          if (regmatch(ep, next))
            return(1);
      }
      return false;
      break;
      }
    case END:
      return true;  /* Success! */
      break;
    default:
      sqd_regerror("regexp corruption");
      return false;
      break;
    }
  }

  /*
   * We get here only if there's trouble -- normally "case END" is
   * the terminating point.
   */
  sqd_regerror("corrupted pointers");
  return false;
}


size_t
regrepeat(
  struct exec *ep,
  char *node
){
  size_t count;
  char *scan;
  char ch;

  switch (OP(node)) {
  case ANY:
    return(strlen(ep->reginput));
    break;
  case EXACTLY:
    ch = *OPERAND(node);
    count = 0;
    for (scan = ep->reginput; *scan == ch; scan++)
      count++;
    return(count);
    break;
  case ANYOF:
    return(strspn(ep->reginput, OPERAND(node)));
    break;
  case ANYBUT:
    return(strcspn(ep->reginput, OPERAND(node)));
    break;
  default:    /* Oh dear.  Called inappropriately. */
    sqd_regerror("internal error: bad call of regrepeat");
    return(0);  /* Best compromise. */
    break;
  }
  /* NOTREACHED */
}


char*
regnext(
  char *p
){
  const int offset = NEXT(p);

  if (offset == 0)
    return(NULL);

  return((OP(p) == BACK) ? p-offset : p+offset);
}

#ifdef DEBUG

char *regprop();


void
regdump(
  sqd_regexp *r
){
  char *s;
  char op = EXACTLY;  /* Arbitrary non-END op. */
  char *next;


  s = r->program + 1;
  while (op != END) {  /* While that wasn't END last time... */
    op = OP(s);
    printf("%2d%s", s-r->program, regprop(s));  /* Where, what. */
    next = regnext(s);
    if (next == NULL)    /* Next ptr. */
      printf("(0)");
    else
      printf("(%d)", (s-r->program)+(next-s));
    s += 3;
    if (op == ANYOF || op == ANYBUT || op == EXACTLY) {
      /* Literal string, where present. */
      while (*s != '\0') {
        putchar(*s);
        s++;
      }
      s++;
    }
    putchar('\n');
  }

  /* Header fields of interest. */
  if (r->regstart != '\0')
    printf("start `%c' ", r->regstart);
  if (r->reganch)
    printf("anchored ");
  if (r->regmust != NULL)
    printf("must have \"%s\"", r->regmust);
  printf("\n");
}


char *
regprop(
  char *op
){
  char *p;
  static char buf[50];

  (void) strcpy(buf, ":");

  switch (OP(op)) {
  case BOL:
    p = "BOL";
    break;
  case EOL:
    p = "EOL";
    break;
  case ANY:
    p = "ANY";
    break;
  case ANYOF:
    p = "ANYOF";
    break;
  case ANYBUT:
    p = "ANYBUT";
    break;
  case BRANCH:
    p = "BRANCH";
    break;
  case EXACTLY:
    p = "EXACTLY";
    break;
  case NOTHING:
    p = "NOTHING";
    break;
  case BACK:
    p = "BACK";
    break;
  case END:
    p = "END";
    break;
  case OPEN+1:
  case OPEN+2:
  case OPEN+3:
  case OPEN+4:
  case OPEN+5:
  case OPEN+6:
  case OPEN+7:
  case OPEN+8:
  case OPEN+9:
    sprintf(buf+strlen(buf), "OPEN%d", OP(op)-OPEN);
    p = NULL;
    break;
  case CLOSE+1:
  case CLOSE+2:
  case CLOSE+3:
  case CLOSE+4:
  case CLOSE+5:
  case CLOSE+6:
  case CLOSE+7:
  case CLOSE+8:
  case CLOSE+9:
    sprintf(buf+strlen(buf), "CLOSE%d", OP(op)-CLOSE);
    p = NULL;
    break;
  case STAR:
    p = "STAR";
    break;
  case PLUS:
    p = "PLUS";
    break;
  default:
    sqd_regerror("corrupted opcode");
    break;
  }
  if (p != NULL)
    (void) strcat(buf, p);
  return(buf);
}
#endif


void
sqd_regsub(
  const sqd_regexp *rp,
  const char *source,
  char *dest
){
  sqd_regexp * const prog = (sqd_regexp *)rp;
  char *src = (char *)source;
  char *dst = dest;
  char c;
  int no;
  size_t len;

  if (prog == NULL || source == NULL || dest == NULL) {
    sqd_regerror("NULL parameter to sqd_regsub");
    return;
  }
  if ((unsigned char)*(prog->program) != SQD_REGMAGIC) {
    sqd_regerror("damaged regexp");
    return;
  }

  while ((c = *src++) != '\0') {
    if (c == '&')
      no = 0;
    else if (c == '\\' && isdigit((int) (*src)))
      no = *src++ - '0';
    else
      no = -1;

    if (no < 0) {  /* Ordinary character. */
      if (c == '\\' && (*src == '\\' || *src == '&'))
        c = *src++;
      *dst++ = c;
    } else if (prog->startp[no] != NULL && prog->endp[no] != NULL &&
          prog->endp[no] > prog->startp[no]) {
      len = prog->endp[no] - prog->startp[no];
      (void) strncpy(dst, prog->startp[no], len);
      dst += len;
      if (*(dst-1) == '\0') {  /* strncpy hit NUL. */
        sqd_regerror("damaged match string");
        return;
      }
    }
  }
  *dst++ = '\0';
}


void
sqd_regerror(
  char *s
){
  fprintf(stderr, "regexp(3): %s\n", s);
  exit(EXIT_FAILURE);
  /* NOTREACHED */
}

#ifdef NBA_TEAM_IN_STL
int
main(
  int argc,
  char **argv
){
  char *pat;
  int   ntok;
  char *s;
  int   status;

  pat  = argv[1];
  ntok = atoi(argv[2]);
  s    = argv[3];

  status = Strparse(pat, s, ntok);
  if (status == 0) {
    printf("no match\n");
  } else {
    int i;
    printf("MATCH.\n");
    for (i = 1; i <= ntok; i++)
      printf("matched token %1d:  %s\n", i, sqd_parse[i]);
  }
}
#endif /*NBA_TEAM_IN_STL*/

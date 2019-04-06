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

#include "squid.hpp"


/* all code below is:
 * Copyright (c) 1986, 1993, 1995 by University of Toronto.
 * Written by Henry Spencer.  Not derived from licensed software.
 *
 * Permission is granted to anyone to use this software for any
 * purpose on any computer system, and to redistribute it in any way,
 * subject to the following restrictions:
 *
 * 1. The author is not responsible for the consequences of use of
 * 	this software, no matter how awful, even if they arise
 * 	from defects in it.
 *
 * 2. The origin of this software must not be misrepresented, either
 * 	by explicit claim or by omission.
 *
 * 3. Altered versions must be plainly marked as such, and must not
 * 	be misrepresented (by explicit claim or omission) as being
 *	the original software.
 *
 * 4. This notice must not be removed or altered.
 */
// This has been altered to fall in line with newer coding standards.

/*
 * sqd_regcomp and sqd_regexec -- sqd_regsub and sqd_regerror are elsewhere
 */

/*
 * The first byte of the regexp internal "program" is actually this magic
 * number; the start node begins in the second byte.
 */
#define	SQD_REGMAGIC	0234

/*
 * The "internal use only" fields in regexp.h are present to pass info from
 * compile to execute that permits the execute phase to run lots faster on
 * simple cases.  They are:
 *
 * regstart	char that must begin a match; '\0' if none obvious
 * reganch	is the match anchored (at beginning-of-line only)?
 * regmust	string (pointer into program) that match must include, or NULL
 * regmlen	length of regmust string
 *
 * Regstart and reganch permit very fast decisions on suitable starting points
 * for a match, cutting down the work a lot.  Regmust permits fast rejection
 * of lines that cannot possibly match.  The regmust tests are costly enough
 * that sqd_regcomp() supplies a regmust only if the r.e. contains something
 * potentially expensive (at present, the only such thing detected is * or +
 * at the start of the r.e., which can involve a lot of backup).  Regmlen is
 * supplied because the test in sqd_regexec() needs it and sqd_regcomp() is computing
 * it anyway.
 */

/*
 * Structure for regexp "program".  This is essentially a linear encoding
 * of a nondeterministic finite-state machine (aka syntax charts or
 * "railroad normal form" in parsing technology).  Each node is an opcode
 * plus a "next" pointer, possibly plus an operand.  "Next" pointers of
 * all nodes except BRANCH implement concatenation; a "next" pointer with
 * a BRANCH on both ends of it is connecting two alternatives.  (Here we
 * have one of the subtle syntax dependencies:  an individual BRANCH (as
 * opposed to a collection of them) is never concatenated with anything
 * because of operator precedence.)  The operand of some types of node is
 * a literal string; for others, it is a node leading into a sub-FSM.  In
 * particular, the operand of a BRANCH node is the first node of the branch.
 * (NB this is *not* a tree structure:  the tail of the branch connects
 * to the thing following the set of BRANCHes.)  The opcodes are:
 */

/* definition	number	opnd?	meaning */
#define	END	0	/* no	End of program. */
#define	BOL	1	/* no	Match beginning of line. */
#define	EOL	2	/* no	Match end of line. */
#define	ANY	3	/* no	Match any character. */
#define	ANYOF	4	/* str	Match any of these. */
#define	ANYBUT	5	/* str	Match any but one of these. */
#define	BRANCH	6	/* node	Match this, or the next..\&. */
#define	BACK	7	/* no	"next" ptr points backward. */
#define	EXACTLY	8	/* str	Match this string. */
#define	NOTHING	9	/* no	Match empty string. */
#define	STAR	10	/* node	Match this 0 or more times. */
#define	PLUS	11	/* node	Match this 1 or more times. */
#define	OPEN	20	/* no	Sub-RE starts here. */
			/*	OPEN+1 is number 1, etc. */
#define	CLOSE	30	/* no	Analogous to OPEN. */

/*
 * Opcode notes:
 *
 * BRANCH	The set of branches constituting a single choice are hooked
 *		together with their "next" pointers, since precedence prevents
 *		anything being concatenated to any individual branch.  The
 *		"next" pointer of the last BRANCH in a choice points to the
 *		thing following the whole choice.  This is also where the
 *		final "next" pointer of each individual branch points; each
 *		branch starts with the operand node of a BRANCH node.
 *
 * BACK		Normal "next" pointers all implicitly point forward; BACK
 *		exists to make loop structures possible.
 *
 * STAR,PLUS	'?', and complex '*' and '+', are implemented as circular
 *		BRANCH structures using BACK.  Simple cases (one character
 *		per match) are implemented with STAR and PLUS for speed
 *		and to minimize recursive plunges.
 *
 * OPEN,CLOSE	...are numbered at compile time.
 */

/*
 * A node is one char of opcode followed by two chars of "next" pointer.
 * "Next" pointers are stored as two 8-bit pieces, high order first.  The
 * value is a positive offset from the opcode of the node containing it.
 * An operand, if any, simply follows the node.  (Note that much of the
 * code generation knows about this implicit relationship.)
 *
 * Using two bytes for the "next" pointer is vast overkill for most things,
 * but allows patterns to get big without disasters.
 */
#define	OP(p)		(*(p))
#define	NEXT(p)		(((*((p)+1)&0177)<<8) + (*((p)+2)&0377))
#define	OPERAND(p)	((p) + 3)

/*
 * Utility definitions.
 */
#define	FAIL(m)		{ sqd_regerror(m); return(NULL); }
#define	ISREPN(c)	((c) == '*' || (c) == '+' || (c) == '?')
#define	META		"^$.[()|?+*\\"

/*
 * Flags to be passed up and down.
 */
#define	HASWIDTH	01	/* Known never to match null string. */
#define	SIMPLE		02	/* Simple enough to be STAR/PLUS operand. */
#define	SPSTART		04	/* Starts with * or +. */
#define	WORST		0	/* Worst case. */

/*
 * Work-variable struct for sqd_regcomp().
 */
struct comp {
	char *regparse;		/* Input-scan pointer. */
	int regnpar;		/* () count. */
	char *regcode;		/* Code-emit pointer; &regdummy = don't. */
	char regdummy[3];	/* NOTHING, 0 next ptr */
	long regsize;		/* Code size. */
};
#define	EMITTING(cp)	((cp)->regcode != (cp)->regdummy)



/*
 - sqd_regcomp - compile a regular expression into internal code
 *
 * We can't allocate space until we know how big the compiled form will be,
 * but we can't compile it (and thus know how big it is) until we've got a
 * place to put the code.  So we cheat:  we compile it twice, once with code
 * generation turned off and size counting turned on, and once "for real".
 * This also means that we don't allocate space until we are sure that the
 * thing really will compile successfully, and we never have to move the
 * code and thus invalidate pointers into it.  (Note that it has to be in
 * one piece because free() must be able to free it all.)
 *
 * Beware that the optimization-preparation code in here knows about some
 * of the structure of the compiled regexp.
 */
sqd_regexp*
sqd_regcomp(
	const char *exp
);


/*
 - reg - regular expression, i.e. main body or parenthesized thing
 *
 * Caller must absorb opening parenthesis.
 *
 * Combining parenthesis handling with the base level of regular expression
 * is a trifle forced, but the need to tie the tails of the branches to what
 * follows makes it hard to avoid.
 */
char*
reg(struct comp *cp, bool paren, int *flagp)


/*
 - regbranch - one alternative of an | operator
 *
 * Implements the concatenation operator.
 */
char*
regbranch(
  struct comp *cp,
  int *flagp);


/*
 - regpiece - something followed by possible [*+?]
 *
 * Note that the branching code sequences used for ? and the general cases
 * of * and + are somewhat optimized:  they use the same NOTHING node as
 * both the endmarker for their branch list and the body of the last branch.
 * It might seem that this node could be dispensed with entirely, but the
 * endmarker role is not redundant.
 */
char*
regpiece(
  struct comp *cp,
  int *flagp);


/*
 - regatom - the lowest level
 *
 * Optimization:  gobbles an entire sequence of ordinary characters so that
 * it can turn them into a single node, which is smaller to store and
 * faster to run.  Backslashed characters are exceptions, each becoming a
 * separate node; the code is simpler that way and it's not worth fixing.
 */
char*
regatom(
  struct comp *cp,
  int *flagp);


// regnode - emit a node
char*			/* Location. */
regnode(
  struct comp *cp,
  char op);


// regc - emit (if appropriate) a byte of code
void
regc(
  struct comp *cp,
  char b);


/*
 - reginsert - insert an operator in front of already-emitted operand
 *
 * Means relocating the operand.
 */
void
reginsert(
  struct comp *cp,
  char op,
  char *opnd);


// regtail - set the next-pointer at the end of a node chain
void
regtail(
  struct comp *cp,
  char *p,
  char *val);


// regoptail - regtail on operand of first argument; nop if operandless
void
regoptail(
  struct comp *cp,
  char *p,
  char *val);


//  sqd_regexec and friends

// Work-variable struct for sqd_regexec().
struct exec;


/*
 - sqd_regexec - match a regexp against a string
 */
int
sqd_regexec(
  sqd_regexp *prog,
  const char *str);


// regtry - try match at specific point
bool			/* 0 failure, 1 success */
regtry(
  struct exec *ep,
  sqd_regexp *prog,
  char *string);


/*
 - regmatch - main matching routine
 *
 * Conceptually the strategy is simple:  check to see whether the current
 * node matches, call self recursively to see whether the rest matches,
 * and then act accordingly.  In practice we make some effort to avoid
 * recursion, in particular by going through "ordinary" nodes (that don't
 * need to know whether the rest of the match failed) by a loop instead of
 * by recursion.
 */
bool			/* 0 failure, 1 success */
regmatch(
  struct exec *ep,
  char *prog);


// regrepeat - report how many times something simple would match
size_t
regrepeat(
  struct exec *ep,
  char *node);

// regnext - dig the "next" pointer out of a node
char*
regnext(
  char *p);

#ifdef DEBUG

// regdump - dump a regexp onto stdout in vaguely comprehensible form
void
regdump(
  sqd_regexp *r)


// regprop - printable representation of opcode
char*
regprop(
  char *op);

#endif


// sqd_regsub - perform substitutions after a regexp match
void
sqd_regsub(
  const sqd_regexp *rp,
  const char *source,
  char *dest);


void
sqd_regerror(
	char *s);

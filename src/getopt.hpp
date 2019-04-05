#include "squid.h"

#pragma once

/****************************************************
 * Support for a portable, flexible Getopt()
 ****************************************************/

enum arg_type{
  sqdARG_NONE,  // no argument
  sqdARG_INT,   // something that atoi() can grok
  sqdARG_FLOAT, // something that atof() can grok
  sqdARG_CHAR,  // require single character or digit
  sqdARG_STRING // anything goes
};


/* Structure: opt_s
 * 
 * Structure for declaring options to a main().
 */
struct opt_s {
  char *name;     /* name of option, e.g. "--option1" or "-o" */
  bool single;      /* TRUE if a single letter option           */
  enum arg_type argtype;    /* for typechecking, e.g. sqdARG_INT        */
};


int
Getopt(
  int argc, 
  char **argv, 
  struct opt_s *opt, 
  int nopts, 
  char *usage, 
  int *ret_optind, 
  char **ret_optname, 
  char **ret_optarg);



#ifdef GETOPT_TESTDRIVER 
/* cc -DGETOPT_TESTDRIVER -L ~/lib/squid.linux/ getopt.c -lsquid
 */
struct opt_s OPTIONS[] = {
  { "--test1", false, sqdARG_INT    },
  { "--test2", false, sqdARG_FLOAT  },
  { "--test3", false, sqdARG_STRING },
  { "--test4", false, sqdARG_CHAR   },
  { "-a",      true,  sqdARG_NONE   },
  { "-b",      true,  sqdARG_INT    },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

#endif /*GETOPT_TESTDRIVER*/

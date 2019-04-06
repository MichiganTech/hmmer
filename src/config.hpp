
/* Configurable compile-time parameters and options in HMMER.
 *
 * Because this header may configure the behavior of system headers
 * (for example, LFS support), it must be included before any other
 * header file.
 */

#pragma once

/*****************************************************************
 * 1. This first section consists of compile-time constants that control
 *    HMMER's computational behavior (memory and processor use).
 *    It can be edited and configured manually before compilation.
 *****************************************************************/

/* RAMLIMIT determines the point at which we switch from fast,
 * full dynamic programming to slow, linear-memory divide and conquer
 * dynamic programming algorithms. It is the minimum amount of available
 * RAM on the systems the package will run on. It can be overridden
 * from the Makefile.
 * By default, we assume we have 32 Mb RAM available (per thread).
 */
#ifndef RAMLIMIT
#define RAMLIMIT 1024
#endif


/*****************************************************************
 * 3. The next section probably shouldn't be edited at all, unless
 *    you really know what you're doing. It controls some fundamental
 *    parameters in HMMER that occasionally get reconfigured in
 *    experimental versions, or for variants of HMMER that work on
 *    non-biological alphabets.
 *****************************************************************/

#define INTSCALE    1000.0      /* scaling constant for floats to integer scores   */
#define MAXABET     20          /* maximum size of alphabet (4 or 20)              */
#define MAXCODE     24          /* maximum degenerate alphabet size (17 or 24)     */
#define MAXDCHLET   200          /* maximum # Dirichlet components in mixture prior */
#define NINPUTS     4          /* number of inputs into structural prior          */
#define INFTY       987654321   /* infinity for purposes of integer DP cells       */
#define NXRAY       4           /* number of structural inputs                */
#define LOGSUM_TBL  20000       /* controls precision of ILogsum()            */


/*****************************************************************
 * 4. The final section isn't human editable at all.
 *    It is configured automatically by the ./configure script.
 *    DO NOT EDIT BELOW THIS LINE.
 *****************************************************************/

/* Version info - set once for whole package in configure.ac
 */
#define PACKAGE_NAME "HMMER"
#define PACKAGE_VERSION "2.5j"
#define PACKAGE_DATE "December 2006"
#define PACKAGE_COPYRIGHT "Copyright (C) 1992-2006 HHMI Janelia Farm"
#define PACKAGE_LICENSE "Freely distributed under the GNU General Public License (GPL)"
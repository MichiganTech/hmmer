/* stopwatch.h
 * SRE, Fri Nov 26 14:54:21 1999 [St. Louis] [HMMER]
 * SRE, Thu Aug  3 08:00:35 2000 [St. Louis] [moved to SQUID]
 * CVS $Id: stopwatch.h,v 1.2 2000/08/03 22:24:38 eddy Exp $
 *
 * Header file for stopwatch.c module:
 * reporting of cpu/system/elapsed time used by a process.
 * See stopwatch.c comments for documentation of compile-time
 * configuration options and API.
 *
 *****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************
 */

#pragma once

#include <stdio.h>
#include <time.h>


struct stopwatch_s {
  time_t t0;			/* Wall clock time, ANSI time()  */
  struct tms cpu0;		/* CPU/system time, POSIX times()*/

  double elapsed;		/* elapsed time, seconds */
  double user;			/* CPU time, seconds */
  double sys;			/* system time, seconds */
};
typedef struct stopwatch_s Stopwatch_t;


/* Function: StopwatchStart()
 *
 * Purpose:  Start a stopwatch.
 *
 * Args:     w - the watch
 */
void
StopwatchStart(
  Stopwatch_t *w);


/* Function: StopwatchStop()
 *
 * Purpose:  Stop a stopwatch.
 *
 *           The implementation allows "split times":
 *           you can stop a watch multiple times, reporting
 *           times at multiple points during program
 *           execution.
 *
 * Args:     w - the watch
 */
void
StopwatchStop(
  Stopwatch_t *w);


/* Function: StopwatchInclude()
 *
 * Purpose:  Merge the cpu and system times from a slave into
 *           a master stopwatch. Both watches must be
 *           stopped, and should not be stopped again unless
 *           You Know What You're Doing.
 *          
 *           Elapsed time is *not* merged; master is assumed
 *           to be keeping track of the wall clock time,
 *           and the slave/worker watch is ignored.
 *          
 *           Used in threads, for broken pthreads/times()
 *           implementations that lose track of cpu times used
 *           by spawned threads.
 *             
 * Args:     w1 - the master stopwatch
 *           w2 - the slave/worker watch
 *
 */
void
StopwatchInclude(
  Stopwatch_t *w1,
  Stopwatch_t *w2);


/* Function: StopwatchAlloc(), StopwatchZero(), StopwatchCopy(),
 *           StopwatchFree()
 *
 * Purpose:  The usual creation/manipulation/destruction routines
 *           for a stopwatch object.
 */
Stopwatch_t*
StopwatchCreate();


void StopwatchZero(
  Stopwatch_t *w);


void
StopwatchCopy(
  Stopwatch_t *w1,
  Stopwatch_t *w2);


void
StopwatchFree(
  Stopwatch_t *w);


/* Function: StopwatchDisplay()
 *
 * Purpose:  Output a usage summary line from a *stopped*
 *           stopwatch (the times will reflect the last
 *           time StopwatchStop() was called.)
 *          
 *           For s = "CPU Time: " an example output line is:
 *           CPU Time: 142.55u 7.17s 149.72 Elapsed: 00:02:35.00
 *
 * Args:     fp - open file for writing (stdout, possibly)
 *           s  - prefix for the report line
 *           w  - a (recently stopped) stopwatch    
 *
 */
void
StopwatchDisplay(
  FILE *fp,
  char *s,
  Stopwatch_t *w);


/* Function: format_time_string()
 *
 * Purpose:  Given a number of seconds, format into
 *           hh:mm:ss.xx in a provided buffer.
 *
 * Args:     buf     - allocated space (128 is plenty!)
 *           sec     - number of seconds
 *           do_frac - TRUE (1) to include hundredths of a sec
 */
void
format_time_string(
  char *buf,
  double sec,
  int do_frac
);

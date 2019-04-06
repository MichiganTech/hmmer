/************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 ************************************************************/

/* stopwatch.c
 *
 * Reporting of cpu/system/elapsed time used by a process.
 * thanks to Warren Gish for assistance.
 *
 * Basic API:
 *
 *   Stopwatch_t *w;
 *   w = StopwatchCreate();
 *  
 *   StopwatchStart(w);
 *   do_lots_of_stuff;
 *   StopwatchStop(w);
 *   StopwatchDisplay(stdout, "CPU time: ", w);
 *  
 *   StopwatchFree(w);
 *  
 * Some behavior can be controlled at compile time by #define's:
 *
 *   SRE_STRICT_ANSI:  By default, stopwatch module assumes that a
 *         machine is POSIX-compliant (e.g. has struct tms, sys/times.h,
 *         and times()). If compiled with -DSRE_STRICT_ANSI, reverts to
 *         pure ANSI C conformant implementation. This simpler system
 *         won't report system times, only user and elapsed times.
 *        
 * One additional compile-time configuration note:       
 *   PTHREAD_TIMES_HACK: Linux pthreads, as of RH6.0/glibc-devel-2.1.1-6,
 *         appears to interact poorly with times() -- usage times in all
 *         but the master thread are lost. A workaround for this bug is
 *         to run stopwatches in each worker thread, and accumulate those
 *         times back into the master stopwatch using StopwatchInclude().
 *          In HMMER, this behavior is compiled in with -DPTHREAD_TIMES_HACK.
 *         No changes are made in stopwatch functions themselves, though;
 *         all the extra code is HMMER code. See hmmcalibrate.c for
 *         an example.
 *
 * See hmmcalibrate.c for examples of more complex usage
 * in dealing with pthreads.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "stopwatch.hpp"


void
format_time_string(
  char *buf,
  double sec,
  int do_frac
){
  int h, m, s, hs;
 
  h  = (int) (sec / 3600.);
  m  = (int) (sec / 60.) - h * 60;
  s  = (int) (sec) - h * 3600 - m * 60;
  if (do_frac) {
    hs = (int) (sec * 100.) - h * 360000 - m * 6000 - s * 100;
    sprintf(buf, "%02d:%02d:%02d.%02d", h,m,s,hs);
  } else {
    sprintf(buf, "%02d:%02d:%02d", h,m,s);
  }
}


void
StopwatchStart(
  Stopwatch_t *w
){
  w->t0 = time(NULL);
#ifdef SRE_STRICT_ANSI
  w->cpu0 = clock();
#else
  (void) times(&(w->cpu0));
#endif

  w->elapsed = 0.;
  w->user    = 0.;
  w->sys     = 0.;
}


void
StopwatchStop(
  Stopwatch_t *w
){
  time_t t1;
  struct tms cpu1;
  long       clk_tck;

  t1 = time(NULL);
  w->elapsed = difftime(t1, w->t0);

  (void) times(&cpu1);
 
  clk_tck = sysconf(_SC_CLK_TCK);
  w->user = (double) (cpu1.tms_utime + cpu1.tms_cutime -
		      w->cpu0.tms_utime - w->cpu0.tms_cutime) /
            (double) clk_tck;

  w->sys  = (double) (cpu1.tms_stime + cpu1.tms_cstime -
		      w->cpu0.tms_stime - w->cpu0.tms_cstime) /
            (double) clk_tck;
}


void
StopwatchInclude(Stopwatch_t *w1, Stopwatch_t *w2)
{
  w1->user    += w2->user;
  w1->sys     += w2->sys;
}


Stopwatch_t *
StopwatchCreate(
){
  Stopwatch_t *w;
  w = malloc(sizeof(Stopwatch_t));
  return w;
}


void
StopwatchZero(
  Stopwatch_t *w
){
  w->elapsed = 0.;
  w->user    = 0.;
  w->sys     = 0.;
}


void
StopwatchCopy(
  Stopwatch_t *w1,
  Stopwatch_t *w2
){
  w1->t0   = w2->t0;
  w1->cpu0.tms_utime = w2->cpu0.tms_utime;
  w1->cpu0.tms_stime = w2->cpu0.tms_stime;
  w1->cpu0.tms_cutime = w2->cpu0.tms_cutime;
  w1->cpu0.tms_cstime = w2->cpu0.tms_cstime;
  w1->elapsed = w2->elapsed;
  w1->user    = w2->user;
  w1->sys     = w2->sys;
}


void
StopwatchFree(
  Stopwatch_t *w
){
  free(w);
}



void
StopwatchDisplay(
  FILE *fp,
  char *s,
  Stopwatch_t *w
){
  char buf[128];	/* (safely holds up to 10^14 years) */
 
  if (s == NULL)
    fputs("CPU Time: ", fp);
  else
    fputs(s, fp);

  format_time_string(buf, w->user+w->sys, 1);
  fprintf(fp, "%.2fu %.2fs %s ", w->user, w->sys, buf);

  format_time_string(buf, w->elapsed, 0);
  fprintf(fp, "Elapsed: %s\n", buf);
}


#ifdef TESTDRIVER
int
main(
  int argc,
  char **argv
){
  Stopwatch_t stopwatch;

  StopwatchStart(&stopwatch);

  sleep(5);

  StopwatchStop(&stopwatch);
  StopwatchDisplay(stdout, "CPU Time: ", &stopwatch);
}
#endif

/*
 @copyright Copyright (C) 2017 Deutsches Klimarechenzentrum GmbH (DKRZ)

 @author Jörg Behrens <behrens@dkrz.de>
         Hendryk Bockelmann <bockelmann@dkrz.de>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the DKRZ GmbH nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * \file sct_reporter.c
 * all functions used to report the measurements
 */
#include "sct_reporter.h"

#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_STRING_H
#include <string.h>
#else
#error missing string.h header
#endif

#ifdef HAVE_STRINGS_H
#include <strings.h>
#else
#error missing strings.h header
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error missing unistd.h header
#endif

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#else
#error missing sys/stat.h header
#endif

#ifdef HAVE_LIBHDF5
#include "hdf5.h"
#endif

#include "sct_collector.h"
#include "sct_reduce.h"

#define NESTEDINDENT 2
#define EVENTINDENT  3
// minimum header length
#define MHL 223

static void report_reduction(sct_reduction_type *red, int timer_choice);
static void report_hdf5(int timer_choice);
static void time_sec_str(const double t, char *s, const int sn);

static FILE *outstream = NULL;
static char outfilename[256];

#ifdef HAVE_LIBHDF5
static hid_t   file_id;
static hsize_t timer_dims[3]; // timer matrix per context, nprocs * nthreads * ntimer
static hsize_t timername_dims[2]; // timername matrix per context, nprocs * ntimer
#endif

// implementation:
#ifdef NESTED_TIMER
static void print_report_hierarchical(const sct_reduction_type *red, const int it,
                                      const int tid, const int iproc, const int nd,
                                      const int mode, const sct_stats_type *x,
                                      const int timer_offset) {

  const int simple_mode = 0;
  const int callstats_mode = 1;
#ifndef HAVE_LIBPAPI
  const int mnl = sct_internal_get_max_name_length();
#else
  const int mnl = (sct_get_eventcounters() == 2) ? sct_internal_get_max_name_length()+7 : sct_internal_get_max_name_length();
#endif
  const double delta = 1.e-6; // used as save margin in plausibility tests

  int rthread_num = red->r_thread_num;

#ifdef HAVE_MPI
  char *cname = sct_get_global_timer_cname(red, iproc, it);
#else
  char *cname = sct_get_timer_cname(it);
#endif

  char sum_str[32];
#ifdef HAVE_LIBPAPI
  char ename[mnl];
  double val;
#endif

  if (mode == callstats_mode) {
    char min_str[32];
    char avg_str[32];
    char max_str[32];
    double a = x->tsum/x->cnum;

    if (x->tmax*(1.0+delta)+delta < a) sct_abort("internal error", __FILE__, __LINE__);
    time_sec_str(x->tmin, min_str, sizeof(min_str));
    time_sec_str(x->tmax, max_str, sizeof(max_str));
    time_sec_str(x->tsum, sum_str, sizeof(sum_str));
    time_sec_str(a, avg_str, sizeof(avg_str));

    // print timer on actual level
    if (nd == 0) {
#ifdef HAVE_MPI
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %-*s | %8d   serial |%7d %s%s%s%s\n", mnl, cname, iproc,
                x->cnum, min_str, avg_str, max_str, sum_str);
      else
        fprintf(outstream," %-*s | %8d %8d |%7d %s%s%s%s\n", mnl, cname, iproc, tid,
                x->cnum, min_str, avg_str, max_str, sum_str);
#else
      fprintf(outstream," %-*s | %8d |%7d %s%s%s%s\n", mnl, cname, iproc, x->cnum,
              min_str, avg_str, max_str, sum_str);
#endif
#else
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %-*s |   serial |%7d %s%s%s%s\n", mnl, cname, x->cnum,
                min_str, avg_str, max_str, sum_str);
      else
        fprintf(outstream," %-*s | %8d |%7d %s%s%s%s\n", mnl, cname, tid, x->cnum,
                min_str, avg_str, max_str, sum_str);
#else
      fprintf(outstream," %-*s |%7d %s%s%s%s\n", mnl, cname, x->cnum, min_str, avg_str,
              max_str, sum_str);
#endif
#endif

#ifdef HAVE_LIBPAPI
      for (int ie=0; ie<sct_get_event_num(); ie++) {
        strcpy(ename,sct_get_event_cname(ie));
        if (sct_get_eventcounters() == 1) {
          val = x->esum[ie]/x->cnum;
#ifdef HAVE_MPI
#ifdef _OPENMP
          fprintf(outstream," %*s%-*s |                   |        %10.2e%10.2e%10.2e%10.2e\n",
                  nd+EVENTINDENT, ".", mnl-nd-EVENTINDENT, ename,
                  (double)x->emin[ie], val, (double)x->emax[ie], (double)x->esum[ie]);
#else
          fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e%10.2e\n",
                  nd+EVENTINDENT, ".", mnl-nd-EVENTINDENT, ename,
                  (double)x->emin[ie], val, (double)x->emax[ie], (double)x->esum[ie]);
#endif
#else
#ifdef _OPENMP
            fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e%10.2e\n",
                    nd+EVENTINDENT, ".", mnl-nd-EVENTINDENT, ename, (double)x->emin[ie],
                    val, (double)x->emax[ie], (double)x->esum[ie]);
#else
          fprintf(outstream," %*s%-*s |        %10.2e%10.2e%10.2e%10.2e\n",
                  nd+EVENTINDENT, ".", mnl-nd-EVENTINDENT, ename, (double)x->emin[ie],
                  val, (double)x->emax[ie], (double)x->esum[ie]);
#endif
#endif
	} else {
          strcat(ename," (rate)");
          val = x->esum[ie]/x->tsum;
#ifdef HAVE_MPI
#ifdef _OPENMP
          fprintf(outstream," %*s%-*s |                   |        %10.2e%10.2e%10.2e\n",
                  nd+EVENTINDENT, ".", mnl-nd-EVENTINDENT, ename, x->rmin[ie], val, x->rmax[ie]);
#else
          fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e\n", nd+EVENTINDENT,
                  ".", mnl-nd-EVENTINDENT, ename, x->rmin[ie], val, x->rmax[ie]);
#endif
#else
#ifdef _OPENMP
	  fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e\n", nd+EVENTINDENT,
                  ".", mnl-nd-EVENTINDENT, ename, x->rmin[ie], val, x->rmax[ie]);
#else
          fprintf(outstream," %*s%-*s |        %10.2e%10.2e%10.2e\n", nd+EVENTINDENT, ".",
                  mnl-nd-EVENTINDENT, ename, x->rmin[ie], val, x->rmax[ie]);
#endif
#endif
	}
      }
#endif

    } else { // if (nd != 0)
#ifdef HAVE_MPI
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %*s%-*s | %8d   serial |%7d %s%s%s%s\n", nd+NESTEDINDENT,
                "L ", mnl-nd-NESTEDINDENT, cname, iproc, x->cnum, min_str, avg_str,
                max_str, sum_str);
      else
        fprintf(outstream," %*s%-*s | %8d %8d |%7d %s%s%s%s\n", nd+NESTEDINDENT, "L ",
                mnl-nd-NESTEDINDENT, cname, iproc, tid, x->cnum, min_str, avg_str, max_str,
                sum_str);
#else
      fprintf(outstream," %*s%-*s | %8d |%7d %s%s%s%s\n", nd+NESTEDINDENT, "L ",
              mnl-nd-NESTEDINDENT, cname, iproc, x->cnum, min_str, avg_str, max_str,
              sum_str);
#endif
#else
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %*s%-*s |   serial |%7d %s%s%s%s\n", nd+NESTEDINDENT, "L ",
                mnl-nd-NESTEDINDENT, cname, x->cnum, min_str, avg_str, max_str, sum_str);
      else
        fprintf(outstream," %*s%-*s | %8d |%7d %s%s%s%s\n", nd+NESTEDINDENT, "L ",
                mnl-nd-NESTEDINDENT, cname, tid, x->cnum, min_str, avg_str, max_str,
                sum_str);
#else
      fprintf(outstream," %*s%-*s |%7d %s%s%s%s\n", nd+NESTEDINDENT, "L ",
              mnl-nd-NESTEDINDENT, cname, x->cnum, min_str, avg_str, max_str, sum_str);
#endif
#endif

#ifdef HAVE_LIBPAPI
      for (int ie=0; ie<sct_get_event_num(); ie++) {
        strcpy(ename,sct_get_event_cname(ie));
        if (sct_get_eventcounters() == 1) {
          val = (double)x->esum[ie]/x->cnum;
#ifdef HAVE_MPI
#ifdef _OPENMP
          fprintf(outstream," %*s%-*s |                   |        %10.2e%10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  (double)x->emin[ie], val, (double)x->emax[ie], (double)x->esum[ie]);
#else
          fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  (double)x->emin[ie], val, (double)x->emax[ie], (double)x->esum[ie]);
#endif
#else
#ifdef _OPENMP
          fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  (double)x->emin[ie], val, (double)x->emax[ie], (double)x->esum[ie]);
#else
          fprintf(outstream," %*s%-*s |        %10.2e%10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  (double)x->emin[ie], val, (double)x->emax[ie], (double)x->esum[ie]);
#endif
#endif
	}
	else {
          strcat(ename," (rate)");
          val = (double)x->esum[ie]/x->tsum;
#ifdef HAVE_MPI
#ifdef _OPENMP
          fprintf(outstream," %*s%-*s |                   |        %10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  x->rmin[ie], val, x->rmax[ie]);
#else
          fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  x->rmin[ie], val, x->rmax[ie]);
#endif
#else
#ifdef _OPENMP
          fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  x->rmin[ie], val, x->rmax[ie]);
#else
          fprintf(outstream," %*s%-*s |        %10.2e%10.2e%10.2e\n",
                  nd+NESTEDINDENT+EVENTINDENT, ".", mnl-nd-NESTEDINDENT-EVENTINDENT, ename,
                  x->rmin[ie], val, x->rmax[ie]);
#endif
#endif
	}
      }
#endif
    } // endif (nd == 0)
  } // if (mode == callstats_mode)
  else if (mode == simple_mode) {
    time_sec_str(x->tsum, sum_str, sizeof(sum_str));
    // print timer on actual level
    if (nd == 0) {
#ifdef HAVE_MPI
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %-*s | %8d   serial |%s\n", mnl, cname, iproc, sum_str);
      else
        fprintf(outstream," %-*s | %8d %8d |%s\n", mnl, cname, iproc, tid, sum_str);
#else
      fprintf(outstream," %-*s | %8d |%s\n", mnl, cname, iproc, sum_str);
#endif
#else
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %-*s |   serial |%s\n", mnl, cname, sum_str);
      else
        fprintf(outstream," %-*s | %8d |%s\n", mnl, cname, tid, sum_str);
#else
      fprintf(outstream," %-*s | %s\n", mnl, cname, sum_str);
#endif
#endif

#ifdef HAVE_LIBPAPI
      for (int ie=0; ie<sct_get_event_num(); ie++) {
        strcpy(ename,sct_get_event_cname(ie));
        if (sct_get_eventcounters() == 1) val = (double)x->esum[ie];
	else{
          strcat(ename," (rate)");
	  val = (double)x->esum[ie]/x->tsum;
        }
#ifdef HAVE_MPI
#ifdef _OPENMP
        fprintf(outstream," %*s%-*s |                   | %10.2e\n", nd+EVENTINDENT, ".",
                mnl-nd-EVENTINDENT, ename, val);
#else
        fprintf(outstream," %*s%-*s |          | %10.2e\n", nd+EVENTINDENT, ".",
                mnl-nd-EVENTINDENT, ename, val);
#endif
#else
#ifdef _OPENMP
        fprintf(outstream," %*s%-*s |          | %10.2e\n", nd+EVENTINDENT, ".", mnl-nd-EVENTINDENT,
                ename, val);
#else
        fprintf(outstream," %*s%-*s | %10.2e\n", nd+EVENTINDENT, ".", mnl-nd-EVENTINDENT,
                ename, val);
#endif
#endif
      }
#endif
    } else { // if (nd != 0)
#ifdef HAVE_MPI
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %*s%-*s | %8d   serial |%s\n", nd+NESTEDINDENT, "L ",
                mnl-nd-NESTEDINDENT, cname, iproc, sum_str);
      else
        fprintf(outstream," %*s%-*s | %8d %8d |%s\n", nd+NESTEDINDENT, "L ",
                mnl-nd-NESTEDINDENT, cname, iproc, tid, sum_str);
#else
      fprintf(outstream," %*s%-*s | %8d |%s\n", nd+NESTEDINDENT, "L ", mnl-nd-NESTEDINDENT,
              cname, iproc, sum_str);
#endif
#else
#ifdef _OPENMP
      if ( red->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
        fprintf(outstream," %*s%-*s |   serial |%s\n", nd+NESTEDINDENT, "L ",
                mnl-nd-NESTEDINDENT, cname, sum_str);
      else
        fprintf(outstream," %*s%-*s | %8d |%s\n", nd+NESTEDINDENT, "L ",
                mnl-nd-NESTEDINDENT, cname, tid, sum_str);
#else
      fprintf(outstream," %*s%-*s | %s\n", nd+NESTEDINDENT, "L ", mnl-nd-NESTEDINDENT,
              cname, sum_str);
#endif
#endif

#ifdef HAVE_LIBPAPI
      for (int ie=0; ie<sct_get_event_num(); ie++) {
        strcpy(ename,sct_get_event_cname(ie));
        if (sct_get_eventcounters() == 1) val = (double)x->esum[ie];
	else{
          strcat(ename," (rate)");
	  val = (double)x->esum[ie]/x->tsum;
        }
#ifdef HAVE_MPI
#ifdef _OPENMP
        fprintf(outstream," %*s%-*s |                   | %10.2e\n", nd+NESTEDINDENT+EVENTINDENT, ".",
                mnl-nd-NESTEDINDENT-EVENTINDENT, ename, val);
#else
        fprintf(outstream," %*s%-*s |          | %10.2e\n", nd+NESTEDINDENT+EVENTINDENT, ".",
                mnl-nd-NESTEDINDENT-EVENTINDENT, ename, val);
#endif
#else
#ifdef _OPENMP
        fprintf(outstream," %*s%-*s |          | %10.2e\n", nd+NESTEDINDENT+EVENTINDENT, ".",
                mnl-nd-NESTEDINDENT-EVENTINDENT, ename, val);
#else
        fprintf(outstream," %*s%-*s | %10.2e\n", nd+NESTEDINDENT+EVENTINDENT, ".",
                mnl-nd-NESTEDINDENT-EVENTINDENT, ename, val);
#endif
#endif
      }
#endif
    } // endif (nd == 0)
  }

  // print all subtimer of actual one
  int offset = timer_offset;
  for (int i = 0; i<red->timer_num_per_rank[iproc]; i++) {
    sct_stats_type *y = &red->stats[offset];
    if ( y->active_under == it)
      print_report_hierarchical(red, i, tid, iproc, nd+1, mode, y, timer_offset);

    offset += rthread_num;
  }
}
#endif // ifdef NESTED_TIMER


static void set_outstream(){
  char *s = getenv("SCT_OUT");
  if (s) {
    if (!strcasecmp(s,"stdout") ) {
      outstream = stdout;
    } else if (!strcasecmp(s,"stderr") ) {
      outstream = stderr;
    } else {
      if (getenv("SCT_FILENAME"))
        strcpy(outfilename, getenv("SCT_FILENAME"));
      else
#ifdef HAVE_LIBHDF5
	strcpy(outfilename, "sct-timings.h5");
#else
	strcpy(outfilename, "sct-timings.txt");
#endif

#ifdef HAVE_LIBHDF5
      if (!strcasecmp(s,"hdf5") ) {
	outstream = NULL;
      }
      else {
#endif
	if (!(outstream=fopen(outfilename, "a"))) {
	  fprintf(stderr,"# SCT-Warning: Cannot open %s. Fallback to stdout output.\n", outfilename);
	  outstream = stdout;
	}
#ifdef HAVE_LIBHDF5
      }
#endif
    }
  } else {
    outstream = stdout;
  }
}


static void unset_outstream(){
  if ( (outstream == stdin) || (outstream == stdout) ) {
    // nothing to do;
  } else {
    if (outstream) {
      fclose(outstream);
    }
  }
  outstream = NULL;
  strcpy(outfilename, "");
}


void sct_single_report(int timer_choice, int proc_choice, int thread_choice, int sp_merging) {
  if (!outstream) set_outstream();

  if (timer_choice > sct_get_timer_num()) {
    sct_warn("output for requested timer exceeds number of registered timers - falling back to output of all timers", __FILE__, __LINE__);

#ifdef HAVE_LIBHDF5
    if (!outstream)
      report_hdf5(SCT_SELECT_ALL);
    else
#endif
      sct_stream_report(NULL, SCT_SELECT_ALL, proc_choice, thread_choice, sp_merging);
  }
  else {
#ifdef HAVE_LIBHDF5
    if (!outstream)
      report_hdf5(SCT_SELECT_ALL);
    else
#endif
      sct_stream_report(NULL, timer_choice, proc_choice, thread_choice, sp_merging);
  }

  if (outstream) unset_outstream();
}


void sct_report(int proc_choice, int thread_choice, int sp_merging) {
  if (!outstream) set_outstream();

#ifdef HAVE_LIBHDF5
  if (!outstream)
    report_hdf5(SCT_SELECT_ALL);
  else
#endif
    sct_stream_report(NULL, SCT_SELECT_ALL, proc_choice, thread_choice, sp_merging);

  if (outstream) unset_outstream();
}


static int get_selection_from_env(const char *key) {
  char *s = getenv(key);
  if (s) {
    if (!strcasecmp(s,"SCT_SELECT_ALL") ) return SCT_SELECT_ALL;
    if (!strcasecmp(s,"SCT_REDUCE_ALL") ) return SCT_REDUCE_ALL;
    int i = atoi(s);
    if (i<0) goto err;
    return i;
  }

 err:
  //fprintf(stderr,"# SCT-Warning: Use SCT_SELECT_ALL as fallback value for %s.\n", key);
  return SCT_SELECT_ALL;
}


static int get_sp_merging_from_env() {
  //fprintf(stdout,"# SCT-Info: get envVar SCT_SP_MERGING.\n");
  char *s = getenv ("SCT_SP_MERGING");
  if ( !s ) goto err;
  if ( !strcasecmp(s,"SCT_SP_MERGE_SIMPLE") ) return SCT_SP_MERGE_SIMPLE;
  if ( !strcasecmp(s,"SCT_SP_SERIAL_ONLY") ) return SCT_SP_SERIAL_ONLY;
  if ( !strcasecmp(s,"SCT_SP_PARALLEL_ONLY") ) return SCT_SP_PARALLEL_ONLY;
  if ( !strcasecmp(s,"SCT_SP_SELECT_ALL") ) return SCT_SP_SELECT_ALL;

 err:
  //fprintf(stderr,"# SCT-Warning: Use SCT_SP_MERGE_SIMPLE as fallback value for SCT_SP_MERGING.\n");
  return SCT_SP_SELECT_ALL;
}


static void report_hdf5(int timer_choice){
  if (timer_choice > sct_get_timer_num())
    sct_abort("output for requested timer exceeds number of registered timers", __FILE__, __LINE__);

#ifdef HAVE_LIBHDF5
  herr_t  status;
  int i, num_timer, num_events;
  const int root_pid = 0;
  int root = 0;
  int my_debug = 0;

  // stop timers if not done already
  if (timer_choice < 0)
    sct_stop_all();
  else
    sct_try_stop(timer_choice);

  // ... for the time being just write out everything ...
  int proc_choice = SCT_SELECT_ALL;
  int thread_choice = SCT_SELECT_ALL;
  int callstatsmode = sct_get_callstats();

  // loop over all contexts and write to file
  sct_context_type *con;
  int jcon = 0;
  for (int icon = 0; icon<sct_get_context_num(); icon++) {
    if ( !(con=sct_get_context(icon)) ) continue; // skip deleted contexts
    if (con->pid == root_pid)
      root = 1;
    else
      root = 0;

#ifdef HAVE_MPI
    sct_create_global_timer_map(con);
#endif
    sct_reduction_type *red = sct_reduction_new(icon,              // int context_choice
                                                proc_choice,       // int proc_choice,
                                                thread_choice,     // int thread_choice
                                                SCT_SP_SELECT_ALL  // int sp_choice
                                                );

    if (con->pid == root_pid) {
      if (red) {
	hid_t attr_type = H5Tcopy(H5T_C_S1);
	hsize_t proc_num = red->proc_num;

	if (jcon == 0) {
	  /* open HDF5 file and set generic info */
	  if ( (file_id = H5Fcreate(outfilename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT)) < 0) {

	    fprintf(stderr, "# SCT-Warning: HDF5 file named '%s' already present\n",
		    outfilename);

	    if( access( outfilename, W_OK ) != -1 ) {
	      struct stat attr;
	      if (stat(outfilename, &attr) == -1)
		sct_abort("not able to access stat of outputfile", __FILE__, __LINE__);

	      /* file already present, try to modify old file by datestring suffix and retry */
	      static char newoutfilename[256], date[20];
	      strcpy(newoutfilename, outfilename);
	      strcat(newoutfilename, "_");
	      strftime(date, 20, "%F_%T", localtime(&(attr.st_mtime)));
	      strcat(newoutfilename, date);
	      if ( rename(outfilename, newoutfilename) != 0 ) {
		sct_abort("not able to rename existing HDF5 outputfile", __FILE__, __LINE__);
	      }
	      fprintf(stderr, "               renaming old file to '%s'\n", newoutfilename);
	     
	    } else {
	      /* need to use a different name for current file since no modification on old file is allowed */
	      strcat(outfilename, "_");
	      strcat(outfilename, strrchr(tmpnam(NULL),'/')+1);
	      fprintf(stderr, "               using '%s' instead\n", outfilename);
	    }
	    /* now retry opening file */
	    if ( (file_id = H5Fcreate(outfilename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
	      sct_abort("not able to create HDF5 output", __FILE__, __LINE__);
	    }
	  }
	  
	  timer_dims[0] = proc_num;
	  timername_dims[0] = proc_num;
	  int threadnum = sct_get_pr_thread_count();
	  timer_dims[1] = 1;
#ifdef _OPENMP
	  timer_dims[1] += threadnum; // one for each thread plus one for pure seq part
#endif
	  // number of timers might vary between tasks !!!
	  num_timer = 0;
	  for (int iproc = 0; iproc < proc_num; iproc++)
            num_timer = (red->timer_num_per_rank[iproc] > num_timer) ? red->timer_num_per_rank[iproc]:num_timer;
	  timername_dims[1] = num_timer;
	  timer_dims[2] = num_timer;
#ifdef HAVE_LIBPAPI
	  num_events = sct_get_event_num();
#else
	  num_events = 0;
#endif
	  
#if (defined(HAVE_MPI) && defined(_OPENMP))
	  char prog_mode[] = "hybrid MPI-OpenMP";
#elif defined(HAVE_MPI)
	  char prog_mode[] = "MPI";
#elif defined(_OPENMP)
	  char prog_mode[] = "OpenMP";
#else
	  char prog_mode[] = "serial";
#endif
	  hid_t attr_id = H5Screate(H5S_SCALAR);

	  /* separate sct internal attributes into group */
	  hid_t group_id = H5Gcreate(file_id, "sct_internal_attributes", H5P_DEFAULT,
	                           H5P_DEFAULT, H5P_DEFAULT);

	  status = H5Tset_size(attr_type, strlen(prog_mode)+1);
	  status = H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
	  hid_t attr = H5Acreate(group_id, "sct program mode", attr_type, attr_id,
				 H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Awrite(attr, attr_type, prog_mode);
	  status = H5Aclose(attr);
	  status = H5Tset_size(attr_type, strlen(print_prg_start_time())+1);
	  status = H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
	  attr = H5Acreate(group_id, "sct start time", attr_type, attr_id,
			   H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Awrite(attr, attr_type, print_prg_start_time());
	  status = H5Aclose(attr);

          char prg_stop_time_string[40];
#if (HAVE_TIME_H || HAVE_SYS_TIME_H)
          struct timeval prg_stop_time;
          struct tm * ptm;
          gettimeofday(&prg_stop_time, NULL);
          ptm = localtime(&prg_stop_time.tv_sec);
          strftime(prg_stop_time_string, sizeof(prg_stop_time_string), "%Y-%m-%d %H:%M:%S", ptm);
#else
          strcpy(prg_stop_time_string,"not available");
#endif
	  status = H5Tset_size(attr_type, strlen(prg_stop_time_string)+1);
	  status = H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
	  attr = H5Acreate(group_id, "sct stop time", attr_type, attr_id,
			   H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Awrite(attr, attr_type, prg_stop_time_string);
	  status = H5Aclose(attr);
          
	  status = H5Gclose(group_id);

	  /* add various attributes if attribute_table is not NULL */
	  const sct_attribute_table *attr_table = sct_internal_get_attribute_table();
	  if ( attr_table != NULL ) {

	    group_id = H5Gcreate(file_id, "report_attributes", H5P_DEFAULT,
				 H5P_DEFAULT, H5P_DEFAULT);

	    sct_attribute_type *a;
	    for (a = attr_table->first; a != NULL; a = a->next) {
	      switch (a->type) {
	      case SCT_INT:
		attr = H5Acreate(group_id, a->key, H5T_NATIVE_INT, attr_id,
				 H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attr, H5T_NATIVE_INT, a->value);
		break;
	      case SCT_LONG:
		attr = H5Acreate(group_id, a->key, H5T_NATIVE_LONG, attr_id,
				 H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attr, H5T_NATIVE_LONG, a->value);
		break;
	      case SCT_FLOAT:
		attr = H5Acreate(group_id, a->key, H5T_NATIVE_FLOAT, attr_id,
				 H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attr, H5T_NATIVE_FLOAT, a->value);
		break;
	      case SCT_DOUBLE:
		attr = H5Acreate(group_id, a->key, H5T_NATIVE_DOUBLE, attr_id,
				 H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attr, H5T_NATIVE_DOUBLE, a->value);
		break;
	      case SCT_STRING:
		status = H5Tset_size(attr_type, strlen(a->value)+1);
		status = H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
		attr = H5Acreate(group_id, a->key, attr_type, attr_id,
				 H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attr, attr_type, a->value);
	      	break;
	      default:
		sct_abort("internal error", __FILE__, __LINE__);
	      }
	      status = H5Aclose(attr);
	    }

	    status = H5Gclose(group_id);
	  }
	  
#ifdef HAVE_MPI
	  attr = H5Acreate(file_id, "number of MPI tasks", H5T_NATIVE_INT, attr_id,
			   H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Awrite(attr, H5T_NATIVE_INT, &timer_dims[0]);
	  status = H5Aclose(attr);
#endif
#ifdef _OPENMP
	  attr = H5Acreate(file_id, "number of OpenMP threads", H5T_NATIVE_INT, attr_id,
			   H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Awrite(attr, H5T_NATIVE_INT, &threadnum);
	  status = H5Aclose(attr);
#endif
	  status = H5Sclose(attr_id);

#ifdef HAVE_MPI
	  int max_pname_mem = red->max_proc_name_len + 1;
	  char (*gpname)[max_pname_mem] = (char (*)[max_pname_mem]) red->proc_names;

	  attr_id = H5Screate_simple(1, &proc_num, NULL);
	  status = H5Tset_size(attr_type, max_pname_mem);
	  status = H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
	  attr = H5Acreate(file_id, "hosts", attr_type, attr_id, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Awrite(attr, attr_type, gpname);
	  status = H5Aclose(attr);
	  status = H5Sclose(attr_id);
#endif
	}

	int r_proc_num = red->red_proc_num;
	int r_thread_num = red->r_thread_num;
        double timer_tmin[r_proc_num][r_thread_num][num_timer];
        double timer_tmax[r_proc_num][r_thread_num][num_timer];
        double timer_tsum[r_proc_num][r_thread_num][num_timer];
        int timer_cnum[r_proc_num][r_thread_num][num_timer];
	char timer_names[r_proc_num][num_timer][sct_internal_get_max_name_length()];
	memset(timer_tmin, 0, r_proc_num*r_thread_num*num_timer*sizeof(double));
	memset(timer_tmax, 0, r_proc_num*r_thread_num*num_timer*sizeof(double));
	memset(timer_tsum, 0, r_proc_num*r_thread_num*num_timer*sizeof(double));
	memset(timer_cnum, 0, r_proc_num*r_thread_num*num_timer*sizeof(int));
	memset(timer_names, 0, r_proc_num*num_timer*sct_internal_get_max_name_length()*sizeof(char));

#ifdef HAVE_LIBPAPI
	char counter_name[sct_internal_get_max_name_length()];
        double counter_val[num_events][r_proc_num][r_thread_num][num_timer];
        double counter_min[num_events][r_proc_num][r_thread_num][num_timer];
        double counter_max[num_events][r_proc_num][r_thread_num][num_timer];
	memset(counter_name, 0, sct_internal_get_max_name_length()*sizeof(char));
	memset(counter_val, 0, r_proc_num*r_thread_num*num_timer*num_events*sizeof(double));
	memset(counter_min, 0, r_proc_num*r_thread_num*num_timer*num_events*sizeof(double));
	memset(counter_max, 0, r_proc_num*r_thread_num*num_timer*num_events*sizeof(double));
#endif

	int timer_offset = 0;

	for (int iproc = 0; iproc < r_proc_num; iproc++) {
	  for (int it = 0; it<red->timer_num_per_rank[iproc]; it++) {
	    if (timer_choice>=0 && timer_choice!=it) {
	      timer_offset += r_thread_num; continue;
	    }
#ifdef HAVE_MPI
	    char *cname = sct_get_global_timer_cname(red, iproc, it);
#else
	    char *cname = sct_get_timer_cname(it);
#endif
	    strcpy(timer_names[iproc][it], cname);

	    for (int tid = 0; tid < r_thread_num; tid++) {

	      sct_stats_type *x = &red->stats[timer_offset];
	      
	      // check whether timer was called at all
	      if (x->cnum < 1) {
		timer_offset++;

		if (my_debug) {
		  fprintf(stderr,"[report_hdf5] at %d %d %d is timer %s not used\n",
			  iproc, tid, it, cname);
		  fflush(stderr);
		}

		continue;
	      }

	      if (callstatsmode) {
		timer_tmin[iproc][tid][it] = x->tmin;
		timer_tmax[iproc][tid][it] = x->tmax;
	      }	
	      timer_tsum[iproc][tid][it] = x->tsum;
	      timer_cnum[iproc][tid][it] = x->cnum;

	      if (my_debug) {
	        fprintf(stderr,"[report_hdf5] at %d %d %d is timer %s with (%d,%f,%f,%f)\n",
			iproc, tid, it, cname, x->cnum, x->tmin, x->tmax, x->tsum);
		fflush(stderr);
	      }

	      // now report eventcounters ... TODO: counter_val should NOT have innermost loop over first index !!!! but actually needed to pass correct memory region to hdf5
#ifdef HAVE_LIBPAPI
	      for (int ie = 0; ie < num_events; ie++) {
		if (sct_get_eventcounters() == 1) {
		  // use sum of raw counters in this case
		  counter_val[ie][iproc][tid][it] = (double)x->esum[ie];
		  if (callstatsmode) {
		    counter_min[ie][iproc][tid][it] = (double)x->emin[ie];
		    counter_max[ie][iproc][tid][it] = (double)x->emax[ie];
		  }
		} else {
		  // use avg of counter rates in this case
		  counter_val[ie][iproc][tid][it] = (double)x->esum[ie]/x->tsum;
		  if (callstatsmode) {
		    counter_min[ie][iproc][tid][it] = x->rmin[ie];
		    counter_max[ie][iproc][tid][it] = x->rmax[ie];
		  }
		}
	      }
#endif
              timer_offset++;
	    }
          }
        }

        /* group for the context which is accessible via "red->context_name.cs" */
	hid_t group_id = H5Gcreate(file_id, red->context_name.cs, H5P_DEFAULT,
				   H5P_DEFAULT, H5P_DEFAULT);

	/* timer names are valid for all following groups */
        hid_t dataspace_id = H5Screate_simple(2, timername_dims, NULL);

	status = H5Tset_size(attr_type, sct_internal_get_max_name_length());
	status = H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
        hid_t dataset_id = H5Dcreate(group_id, "timer_names", attr_type, dataspace_id,
				     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, attr_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			  timer_names);
        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);

	/* subgroup for timer values */
	hid_t timergroup_id = H5Gcreate(group_id, "time", H5P_DEFAULT,
					H5P_DEFAULT, H5P_DEFAULT);

        dataspace_id = H5Screate_simple(3, timer_dims, NULL);

	if (callstatsmode) {
	  dataset_id = H5Dcreate(timergroup_id, "timer_tmin", H5T_NATIVE_DOUBLE, dataspace_id,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			    timer_tmin);
	  status = H5Dclose(dataset_id);
	  
	  dataset_id = H5Dcreate(timergroup_id, "timer_tmax", H5T_NATIVE_DOUBLE, dataspace_id,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			    timer_tmax);
	  status = H5Dclose(dataset_id);
	}

        dataset_id = H5Dcreate(timergroup_id, "timer_tsum", H5T_NATIVE_DOUBLE, dataspace_id,
	                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			  timer_tsum);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(timergroup_id, "timer_cnum", H5T_NATIVE_INT, dataspace_id,
	                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			  timer_cnum);
        status = H5Dclose(dataset_id);
	status = H5Gclose(timergroup_id);

	/* subgroup for eventcounter values */
#ifdef HAVE_LIBPAPI
	hid_t countergroup_id = H5Gcreate(group_id, "eventcounter", H5P_DEFAULT,
	                                  H5P_DEFAULT, H5P_DEFAULT);
	
	for (int ie = 0; ie < num_events; ie++) {
	  strcpy(counter_name, sct_get_event_cname(ie));
	  if (sct_get_eventcounters() != 1) {
	    strcat(counter_name," (rate)");
	  }

	  dataset_id = H5Dcreate(countergroup_id, counter_name, H5T_NATIVE_DOUBLE, dataspace_id,
	                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			    counter_val[ie]);
	  status = H5Dclose(dataset_id);

	  if (callstatsmode) {
	    strcpy(counter_name, sct_get_event_cname(ie));
	    if (sct_get_eventcounters() != 1) {
	      strcat(counter_name," (rate)");
	    }
	    strcat(counter_name," min");

	    dataset_id = H5Dcreate(countergroup_id, counter_name, H5T_NATIVE_DOUBLE, dataspace_id,
				   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			      counter_min[ie]);
	    status = H5Dclose(dataset_id);

	    strcpy(counter_name, sct_get_event_cname(ie));
	    if (sct_get_eventcounters() != 1) {
	      strcat(counter_name," (rate)");
	    }
	    strcat(counter_name," max");

	    dataset_id = H5Dcreate(countergroup_id, counter_name, H5T_NATIVE_DOUBLE, dataspace_id,
				   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			      counter_max[ie]);
	    status = H5Dclose(dataset_id);
	  }
	}

	status = H5Gclose(countergroup_id);
#endif
	status = H5Gclose(group_id);

        status = H5Sclose(dataspace_id);
      } /* if (red) */
      else
        fprintf(stderr, "# SCT-Warning: empty reduction returned.\n");
    } /* if (root) */

    sct_reduction_delete(red);
    jcon++;
  }

  if (root) H5Fclose(file_id);
#endif
}


void sct_stream_report(FILE *outstream_arg,
                       int timer_choice   ,
                       int proc_choice    ,
                       int thread_choice  ,
                       int sp_merging     ) {

  if (outstream_arg) outstream = outstream_arg;
  if (!outstream) sct_abort("missing outstream", __FILE__, __LINE__);
  if (timer_choice > sct_get_timer_num())
    sct_abort("output for requested timer exceeds number of registered timers", __FILE__, __LINE__);

  // try to get a clean canvas:
  if ( (outstream == stdout) || (outstream == stderr) ) fflush(outstream);

  const int root_pid = 0;
  int root = 0;
  int my_debug = 0;

#ifdef HAVE_MPI
  if (proc_choice == SCT_GETENV) proc_choice = get_selection_from_env("SCT_PROC_CHOICE");
#else
  proc_choice = SCT_SELECT_ALL;
#endif

#ifdef _OPENMP
  if (thread_choice == SCT_GETENV) thread_choice = get_selection_from_env("SCT_THREAD_CHOICE");
  if (sp_merging == SCT_GETENV) sp_merging = get_sp_merging_from_env();
#else
  thread_choice = 0;
  sp_merging = SCT_SP_SELECT_ALL;
#endif

  if (my_debug) {
    fprintf(stderr,"[-] sct_stream_report with %d %d %d %d\n",
            timer_choice, proc_choice, thread_choice, sp_merging);
    fflush(stderr);
  }

  // stop timers if not done already
  if (timer_choice < 0)
    sct_stop_all();
  else
    sct_try_stop(timer_choice);

  sct_context_type *con;
  int jcon = 0;
  for (int icon = 0; icon<sct_get_context_num(); icon++) {
    if ( !(con=sct_get_context(icon)) ) continue; // skip deleted contexts
    if (con->pid == root_pid)
      root = 1;
    else
      root = 0;

    if (jcon == 0 && root) {
      if (timer_choice<0)
        fprintf(outstream, "\n --> sct_report for all timers:\n");
      else
        fprintf(outstream, "\n --> sct_report for timer \"%s\":\n",
                sct_get_timer_cname(timer_choice));
#if (defined(HAVE_MPI) && defined(_OPENMP))
      fprintf(outstream, "program_mode  = MPI-OpenMP\n");
#elif defined(HAVE_MPI)
      fprintf(outstream, "program_mode  = MPI\n");
#elif defined(_OPENMP)
      fprintf(outstream, "program_mode  = OpenMP\n");
#else
      fprintf(outstream, "program_mode  = serial\n");
#endif
      fprintf(outstream, "program start = %s\n", print_prg_start_time());

#ifdef HAVE_MPI
      if (proc_choice >= 0) {
        fprintf(outstream, "proc_choice   = %d\n", proc_choice);
      } else {
        fprintf(outstream, "proc_choice   = %s\n", sct_internal_reduce_doc(proc_choice));
      }
#endif
#ifdef _OPENMP
      if (thread_choice >= 0) {
        fprintf(outstream, "thread_choice = %d\n", thread_choice);
      } else {
        fprintf(outstream, "thread_choice = %s\n", sct_internal_reduce_doc(thread_choice));
      }
      fprintf(outstream,   "sp_merging    = %s\n", sct_internal_sp_doc(sp_merging));
#endif
    }

    /// \todo in case of report for specific timer, reduction of all timers might be to much
#ifdef HAVE_MPI
    sct_create_global_timer_map(con); // \todo also needed in seq/openmp case ???
#endif
    sct_reduction_type *red = sct_reduction_new(icon,          // int context_choice
                                                proc_choice,   // int proc_choice,
                                                thread_choice, // int thread_choice
                                                sp_merging     // int sp_choice
                                                );
    if (root) {
      if (red)
        report_reduction(red, timer_choice);
      else
        fprintf(outstream,   "# sct: empty report\n");
    }

    sct_reduction_delete(red);
    jcon++;
  }

  if ( (outstream == stdout) || (outstream == stderr) ) fflush(outstream);
}


static void print_ranks2host(int a, int b, char *name) {
  if (a<0 || b<0) return;
  const int cbuf_size = 64;
  char cbuf[cbuf_size];
  if (a==b)
    snprintf(cbuf, cbuf_size, "host(%i)",a);
  else
    snprintf(cbuf, cbuf_size, "host(%i:%i)",a,b);
  cbuf[cbuf_size-1]=0;
  fprintf(outstream, "%-13s = %s\n",cbuf,name);
}


static void report_ranks(sct_reduction_type *r) {
#ifdef HAVE_MPI
  if (!r->proc_names) return;
  int proc_num = r->proc_num;
  int max_pname_mem =  r->max_proc_name_len + 1;
  char (*gpname)[max_pname_mem] =  (char (*)[max_pname_mem]) r->proc_names;
  int a = -1;
  int b = -1;

  for (int iproc = 0; iproc < proc_num; iproc++) {
    if (a<0) {
      // start section:
      a = iproc;
      b = iproc;
      continue;
    }
    if (strcmp(&gpname[iproc][0], &gpname[a][0]) == 0) {
      // continue section:
      b++;
      continue;
    } else {
      // report section and start new one:
      print_ranks2host(a, b, &gpname[a][0]);
      a=iproc;
      b=iproc;
    }
  }
  if (a >= 0) {
    // report dangling section:
    print_ranks2host(a, b, &gpname[a][0]);
  }
#endif
}


static void report_reduction(sct_reduction_type *r, int timer_choice) {

  int show_threads;
  int show_procs;
  int thread_stats;
  int proc_stats;
  int mpi_reduction_factor;
  int omp_reduction_factor;
  int proc_num;
  if (!r) return;
#ifdef HAVE_MPI
  const int mpi = 1;
  proc_num = r->proc_num;
  if (r->proc_choice == SCT_REDUCE_ALL) {
    if (r->red_proc_num != 1) sct_abort("unexpected case", __FILE__, __LINE__);
    proc_stats = 1;
    show_procs = 0;
    mpi_reduction_factor = r->proc_num;
  } else {
    proc_stats = 0;
    show_procs = 1;
    mpi_reduction_factor = 1;
  }
  fprintf(outstream, "nprocs        = %i\n", r->proc_num);
#else
  const int mpi = 0;
  proc_num = 1;
  if (r->proc_num != 1 || r->red_proc_num != 1) sct_abort("internal error", __FILE__, __LINE__);
  proc_stats = 0;
  show_procs = 0;
  mpi_reduction_factor = 1;
#endif

#ifdef _OPENMP
  const int omp = 1;
  if (r->thread_choice == SCT_REDUCE_ALL) {
    thread_stats = 1;
    show_threads = 0;
    if (r->r_thread_num != 1) sct_abort("unexpected case", __FILE__, __LINE__);
    omp_reduction_factor = r->m_thread_num;
  } else {
    thread_stats = 0;
    show_threads = 1;
    omp_reduction_factor = 1;
  }
  fprintf(outstream, "nthreads      = %i\n", r->p_thread_num);
#else
  const int omp = 0;
  if (r->m_thread_num != 1 || r->r_thread_num != 1) sct_abort("internal error", __FILE__, __LINE__);
  thread_stats = 0;
  show_threads = 0;
  omp_reduction_factor = 1;
#endif

  const double eps = 1.e-12; // near zero; eps should be below 1.e-9 for nano second resolution
  const double delta = 1.e-6; // used as save margin in plausibility tests
  char min_str[32];
  char avg_str[32];
  char max_str[32];
  char sum_str[32];

  fprintf(outstream, "context       = %s\n", r->context_name.cs);

  int reduce_factor = mpi_reduction_factor * omp_reduction_factor;
  if (reduce_factor < 1) sct_abort("unexpected case", __FILE__, __LINE__);

  char *hsep = NULL;
  char *hsep_proc = NULL;

  const int simple_mode = 0;
  const int callstats_mode = 1;
  const int reduction_mode = 2;
#ifndef HAVE_LIBPAPI
  const int mnl = sct_internal_get_max_name_length();
#else
  const int mnl = (sct_get_eventcounters() == 2) ? sct_internal_get_max_name_length()+7 : sct_internal_get_max_name_length();
#endif

  int mode;
  char head[MHL];
  head[0] = 0;
  if ( (!mpi || !proc_stats) && (!omp || !thread_stats) ) {
    // no reduction at all - we can use callstats if present
    if (sct_get_callstats()) {
      mode = callstats_mode;
      fprintf(outstream, "statistics    = call-statistics\n");
      if (mpi && omp) {
        report_ranks(r);
        snprintf(head,MHL," %-*s | %8s %8s |%7s %10s%10s%10s%10s\n", mnl,"name", "proc", "thread", "#calls", "min", "avg", "max", "sum");
      } else if (mpi) {
        report_ranks(r);
        snprintf(head,MHL," %-*s | %8s |%7s %10s%10s%10s%10s\n",mnl,"name","proc", "#calls", "min", "avg", "max", "sum");
      } else if (omp) {
        snprintf(head,MHL," %-*s | %8s |%7s %10s%10s%10s%10s\n", mnl,"name", "thread", "#calls", "min", "avg", "max", "sum");
      } else {
        snprintf(head,MHL," %-*s |%7s %10s%10s%10s%10s\n", mnl,"name", "#calls", "min", "avg", "max", "sum");
      }
    } else {
      mode = simple_mode;
      int thread_num = r->m_thread_num;
      if (omp && r->sp_merging == SCT_SP_SELECT_ALL && r->r_thread_num > 1)
        thread_num--; // we have an extra portion for the serial part allocated
      if (mpi && omp) {
        fprintf(outstream, "statistics    = none(%d x %d)\n", r->proc_num, thread_num);
        report_ranks(r);
        snprintf(head,MHL," %-*s | %8s %8s |%10s\n", mnl,"name", "proc", "thread", "sum");
      } else if (mpi) {
        fprintf(outstream, "statistics    = none(%d)\n",  r->proc_num);
        report_ranks(r);
        snprintf(head,MHL," %-*s | %8s |%10s\n", mnl,"name", "proc", "sum");
      } else if (omp) {
        fprintf(outstream, "statistics    = none(%d)\n", thread_num);
        snprintf(head,MHL," %-*s | %8s |%10s\n", mnl,"name", "thread", "sum");
      } else {
        fprintf(outstream, "statistics    = none\n");
        snprintf(head,MHL," %-*s |%10s\n", mnl,"name", "sum");
      }
    }
  } else {
    mode = reduction_mode; // means parallel statistics
    //fprintf(outstream,"show_procs=%i, thread_stats=%i\n",show_procs, thread_stats);
    if (show_procs && show_threads) {
      sct_abort("unexpected case", __FILE__, __LINE__);
    } else if (show_procs && thread_stats) {
      fprintf(outstream, "statistics    = thread(%d)-statistics for each process(%d)\n", r->m_thread_num,  r->proc_num);
      snprintf(head,MHL," %-*s | %8s |%10s%10s%10s%10s  %7s\n", mnl,"name", "proc", "min", "avg", "max", "sum", "lbe[%]");
    } else if (proc_stats && show_threads) {
      fprintf(outstream, "statistics    = process(%d)-statistics for each thread(%d)\n", r->proc_num,  r->m_thread_num);
      snprintf(head,MHL," %-*s | %8s |%10s%10s%10s%10s  %7s\n", mnl,"name", "thread", "min", "avg", "max", "sum", "lbe[%]");
    } else if (proc_stats && thread_stats) {
      fprintf(outstream, "statistics    = hybrid(%d x %d)-statistics\n",r->proc_num,  r->m_thread_num);
      snprintf(head,MHL," %-*s |%10s%10s%10s%10s  %7s\n", mnl,"name", "min", "avg", "max", "sum", "lbe[%]");
    } else if (!mpi && thread_stats) {
      fprintf(outstream, "statistics    = thread(%d)-statistics\n", r->m_thread_num);
      snprintf(head,MHL," %-*s |%10s%10s%10s%10s  %7s\n", mnl,"name", "min", "avg", "max", "sum", "lbe[%]");
    } else if (!omp && proc_stats) {
      fprintf(outstream, "statistics    = process-statistics\n");
      snprintf(head,MHL," %-*s |%10s%10s%10s%10s  %7s\n", mnl,"name", "min", "avg", "max", "sum", "lbe[%]");
    } else {
      // all other cases should have beeen handled in callstats_mode or simple_mode
      sct_abort("unexpected case", __FILE__, __LINE__);
    }
  }
  if (head[0]) {
    int hsep_len = strlen(head);
    hsep      = calloc(hsep_len+2, sizeof(*hsep));
    hsep_proc = calloc(3*hsep_len+3, sizeof(*hsep_proc));
    memset(hsep, '-', hsep_len);
    memset(hsep_proc, '-', hsep_len);
    hsep[hsep_len] = '\n';
    hsep_proc[hsep_len] = '\n';
    memcpy(hsep_proc+hsep_len+1, head, hsep_len);
    memset(hsep_proc+2*hsep_len+1, '-', hsep_len);
    hsep_proc[3*hsep_len+1] = '\n';
    fprintf(outstream, "%s", hsep_proc);
  } else {
    hsep = NULL;
  }

  int rproc_num = r->red_proc_num;
  int rthread_num = r->r_thread_num;
  int report;

#ifdef NESTED_TIMER
  // report hierarchical in case of no reduction
  if ( (mode == callstats_mode || mode == simple_mode) &&
       (r->sp_merging != SCT_SP_MERGE_SIMPLE) && sct_get_nestedtimers() ){

    int timer_offset = 0;
    int timer_offset_proc_start;
    int timer_offset_thread_start;

    for (int iproc = 0; iproc < rproc_num; iproc++) {
      timer_offset_proc_start = timer_offset;

#ifdef _OPENMP
      if ( r->sp_merging == SCT_SP_SELECT_ALL && rthread_num > 1) {
        // show pure serial part first
        timer_offset_thread_start = timer_offset;

        for (int it = 0; it<r->timer_num_per_rank[iproc]; it++) {
          if (timer_choice>=0 && timer_choice!=it) {timer_offset += rthread_num; continue;}
          sct_stats_type *x = &r->stats[timer_offset + rthread_num-1];
	  //printf("... [%d][%d][%d] : %d , %.2e : %d\n",iproc,it,rthread_num-1,x->cnum,x->tsum,timer_offset + rthread_num-1);
	  if ( mode == simple_mode && x->tsum < eps ) {timer_offset += rthread_num; continue;}
          if ( mode == callstats_mode && x->cnum<1 ) {timer_offset += rthread_num; continue;}
          if ( x->active_under == -1 ) print_report_hierarchical(r, it, rthread_num-1, iproc, 0,
                                                                 mode, x, timer_offset_thread_start + rthread_num-1);
          timer_offset += rthread_num;
        }

        timer_offset = timer_offset_thread_start;
        if (hsep) fprintf(outstream, "%s", hsep);

        // show parallel part for each thread
        for (int ithread = 0; ithread < rthread_num-1; ithread++) {
          timer_offset_thread_start = timer_offset;
          report = 0;

          for (int it = 0; it<r->timer_num_per_rank[iproc]; it++) {
            if (timer_choice>=0 && timer_choice!=it) {timer_offset += rthread_num; continue;}
            sct_stats_type *x = &r->stats[timer_offset];
	    //printf("... [%d][%d][%d] : %d , %.2e : %d\n",iproc,it,ithread,x->cnum,x->tsum, timer_offset);
     	    if ((mode == simple_mode && x->tsum < eps )
            	|| (mode == callstats_mode && x->cnum<1 ) ) {
              timer_offset += rthread_num;
	      continue;
            }
            if (x->active_under == -1) {
              report = 1;
              print_report_hierarchical(r, it, ithread, iproc, 0, mode, x, timer_offset_thread_start); // top-level timer
	    }
            timer_offset += rthread_num;
          }

          timer_offset = timer_offset_thread_start + 1;
          if (hsep && report && ithread<rthread_num-2) fprintf(outstream, "%s", hsep);
        }
      }
      else {
#endif
        for (int ithread = 0; ithread < rthread_num; ithread++) {
          timer_offset_thread_start = timer_offset;

          for (int it = 0; it<r->timer_num_per_rank[iproc]; it++) {
            if (timer_choice>=0 && timer_choice!=it) {timer_offset += rthread_num; continue;}
            sct_stats_type *x = &r->stats[timer_offset];
	    //printf("... [%d][%d][%d] : %d , %.2e : %d\n",iproc,it,ithread,x->cnum,x->tsum, timer_offset);
            // r->stats is organized as [proc][timer][thread] so next timer is timer_offset += rthread_num

     	    if ( mode == simple_mode && x->tsum < eps ) {timer_offset += rthread_num; continue;}
            if ( mode == callstats_mode && x->cnum<1 ) {timer_offset += rthread_num; continue;}

            if ( x->active_under == -1 ) print_report_hierarchical(r, it, ithread, iproc, 0,
                                                                   mode, x, timer_offset_thread_start); // top-level timer
            timer_offset += rthread_num; // increase by number of threads == one timer
          }

          timer_offset = timer_offset_thread_start + 1; // increase by one thread
          if (hsep && ithread<rthread_num-1) fprintf(outstream, "%s", hsep);
        }
#ifdef _OPENMP
      }
#endif

      timer_offset = timer_offset_proc_start + r->timer_num_per_rank[iproc]*rthread_num; // increase by one proc
      if (hsep_proc && iproc<rproc_num-1) fprintf(outstream, "%s",hsep_proc);
    }
  }

  // report in standard way else
  else
#endif
// NESTED_TIMER
  {

#ifdef HAVE_LIBPAPI
  char ename[mnl];
  double min,avg,max,sum;
#endif

  if (mode == simple_mode) {

    int timer_offset = 0;

    for (int iproc = 0; iproc < rproc_num; iproc++) {

      for (int it = 0; it<r->timer_num_per_rank[iproc]; it++) {
        if (timer_choice>=0 && timer_choice!=it) {timer_offset += rthread_num; continue;}

#ifdef HAVE_MPI
        char *cname = sct_get_global_timer_cname(r,iproc, it);
#else
        char *cname = sct_get_timer_cname(it);
#endif

        for (int tid = 0; tid < rthread_num; tid++) {

          sct_stats_type *x = &r->stats[timer_offset];

          if (x->tsum < eps) {timer_offset++; continue;}

          time_sec_str(x->tsum, sum_str, sizeof(sum_str));
          if (mpi && omp) {
            if ( r->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
              fprintf(outstream," %-*s | %8d   serial |%s\n", mnl,cname, iproc, sum_str);
	    else
              fprintf(outstream," %-*s | %8d %8d |%s\n", mnl,cname, iproc, tid, sum_str);
          } else if (mpi) {
            fprintf(outstream," %-*s | %8d |%s\n", mnl,cname, iproc, sum_str);
          } else if (omp) {
            if ( r->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
              fprintf(outstream," %-*s |   serial |%s\n", mnl,cname, sum_str);
            else
              fprintf(outstream," %-*s | %8d |%s\n", mnl,cname, tid, sum_str);
          } else {
            fprintf(outstream," %-*s |%s\n", mnl,cname, sum_str);
          }

          // now report eventcounters
#ifdef HAVE_LIBPAPI
          for (int ie=0; ie<sct_get_event_num(); ie++) {
            strcpy(ename,sct_get_event_cname(ie));
            if (sct_get_eventcounters() == 1) {
              sum = (double)x->esum[ie];
	    } else {
              strcat(ename," (rate)");
              sum = (double)x->esum[ie]/x->tsum;
	    }

            if (mpi && omp) {
              fprintf(outstream," %*s%-*s |                   |%10.2e\n",
                EVENTINDENT, ".", mnl-EVENTINDENT, ename, sum);
            } else if (mpi) {
              fprintf(outstream," %*s%-*s |          |%10.2e\n",
                EVENTINDENT, ".", mnl-EVENTINDENT, ename, sum);
            } else if (omp) {
              fprintf(outstream," %*s%-*s |          |%10.2e\n",
	        EVENTINDENT, ".", mnl-EVENTINDENT, ename, sum);
            } else {
              fprintf(outstream," %*s%-*s |%10.2e\n",
                EVENTINDENT, ".", mnl-EVENTINDENT, ename, sum);
	    }
          }
#endif
          timer_offset++;
        }
      }
      if ( iproc < rproc_num - 1 && hsep_proc ) fprintf(outstream, "%s", hsep_proc);
    }

  } else if (mode == callstats_mode) {

    int timer_offset = 0;

    for (int iproc = 0; iproc < rproc_num; iproc++) {
      for (int it = 0; it<r->timer_num_per_rank[iproc]; it++) {
        if (timer_choice>=0 && timer_choice!=it) {timer_offset += rthread_num; continue;}

#ifdef HAVE_MPI
        char *cname = sct_get_global_timer_cname(r, iproc, it);
#else
        char *cname = sct_get_timer_cname(it);
#endif

        for (int tid = 0; tid < rthread_num; tid++) {
          sct_stats_type *x = &r->stats[timer_offset];
          if (x->cnum<1) {timer_offset++; continue;}
          double a = x->tsum/x->cnum;
          if (x->tmax*(1.0+delta)+delta < a) sct_abort("internal error", __FILE__, __LINE__);
          time_sec_str(x->tmin, min_str, sizeof(min_str));
          time_sec_str(x->tmax, max_str, sizeof(max_str));
          time_sec_str(x->tsum, sum_str, sizeof(sum_str));
          time_sec_str(a, avg_str, sizeof(avg_str));
          if (mpi && omp) {
            if ( r->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
              fprintf(outstream," %-*s | %8d   serial |%7d %s%s%s%s\n", mnl,cname, iproc,
                      x->cnum, min_str, avg_str, max_str, sum_str);
	    else
              fprintf(outstream," %-*s | %8d %8d |%7d %s%s%s%s\n", mnl,cname, iproc, tid,
                      x->cnum, min_str, avg_str, max_str, sum_str);
          } else if (mpi) {
            fprintf(outstream," %-*s | %8d |%7d %s%s%s%s\n", mnl,cname, iproc, x->cnum,
                    min_str, avg_str, max_str, sum_str);
          } else if (omp) {
            if ( r->sp_merging == SCT_SP_SELECT_ALL && tid == rthread_num-1)
              fprintf(outstream," %-*s |   serial |%7d %s%s%s%s\n", mnl, cname, x->cnum,
                      min_str, avg_str, max_str, sum_str);
	    else
              fprintf(outstream," %-*s | %8d |%7d %s%s%s%s\n", mnl, cname, tid, x->cnum,
                      min_str, avg_str, max_str, sum_str);
          } else {
	    fprintf(outstream," %-*s |%7d %s%s%s%s\n", mnl,cname, x->cnum, min_str,
                    avg_str, max_str, sum_str);
          }
          // now report eventcounters
#ifdef HAVE_LIBPAPI
          for (int ie=0; ie<sct_get_event_num(); ie++) {
            strcpy(ename,sct_get_event_cname(ie));
            if (sct_get_eventcounters() == 1) {
              avg = (double)x->esum[ie]/x->cnum;
              min = (double)x->emin[ie];
              max = (double)x->emax[ie];
	      sum = (double)x->esum[ie];
	    } else {
              strcat(ename," (rate)");
              avg = (double)x->esum[ie]/x->tsum;
              min = x->rmin[ie];
              max = x->rmax[ie];
              sum = 0.; /// \todo Sum over rates does not make any sense
            }

            if (mpi && omp) {
              fprintf(outstream," %*s%-*s |                   |        %10.2e%10.2e%10.2e%10.2e\n",
                      EVENTINDENT, ".", mnl-EVENTINDENT, ename, min, avg, max, sum);
            } else if (mpi) {
              fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e%10.2e\n",
                      EVENTINDENT, ".", mnl-EVENTINDENT, ename, min, avg, max, sum);
            } else if (omp) {
              fprintf(outstream," %*s%-*s |          |        %10.2e%10.2e%10.2e%10.2e\n",
                      EVENTINDENT, ".", mnl-EVENTINDENT, ename, min, avg, max, sum);
            } else {
              fprintf(outstream," %*s%-*s |        %10.2e%10.2e%10.2e%10.2e\n",
                      EVENTINDENT, ".", mnl-EVENTINDENT, ename, min, avg, max, sum);
	    }
          }
#endif
          timer_offset++;
        }
      }
      if ( iproc < rproc_num - 1 && hsep_proc ) fprintf(outstream, "%s", hsep_proc);
    }

  } else if (mode == reduction_mode) {

    int timer_num = r->global_timer_num;

    sct_stats_type (*rstats)[timer_num][rthread_num] = 
                   (sct_stats_type (*)[timer_num][rthread_num]) r->stats;

    double q_reduce = 1.0/reduce_factor;
    int entity;

    for (int it = 0; it<timer_num; it++) {
      if (timer_choice>=0 && timer_choice!=it) continue;

      char *cname = sct_get_global_timer_cname(r, 0, it);
      for (int iproc = 0; iproc < rproc_num; iproc++) {
        for (int tid = 0; tid < rthread_num; tid++) {
          sct_stats_type *x = &rstats[iproc][it][tid];
	  //printf("DEBUG [%d]: %7.3f %7.3f %7.3f\n",it,x->tmin,x->tmax,x->tsum);
          if (x->tsum < eps) continue;
          time_sec_str(x->tmin, min_str, sizeof(min_str));
          time_sec_str(x->tmax, max_str, sizeof(max_str));
          time_sec_str(x->tsum, sum_str, sizeof(sum_str));
          double a = x->tsum * q_reduce;
          if (x->tmax*(1.0+delta)+delta < a) sct_abort("internal error", __FILE__, __LINE__);

          time_sec_str(a, avg_str, sizeof(avg_str));
          double e;
          if (x->tmax > eps)
            e   = a/x->tmax;
          else
            e   = 1.0;
          e = e*100; // percent
          if (show_procs && thread_stats)      entity = iproc;
          else if (proc_stats && show_threads) entity = tid;
          else if (   (proc_stats && thread_stats)
                   || (!mpi && thread_stats)
                   || (!omp && proc_stats))    entity = -1;
          else sct_abort("unexpected case", __FILE__, __LINE__);

	  if (entity == -1)
            fprintf(outstream," %-*s |%s%s%s%s  %7.3f\n", mnl, cname,
                    min_str, avg_str, max_str, sum_str, e);
          else
            fprintf(outstream," %-*s | %8d |%s%s%s%s  %7.3f\n", mnl, cname, entity,
                    min_str, avg_str, max_str, sum_str, e);

#ifdef HAVE_LIBPAPI
          for (int ie=0; ie<sct_get_event_num(); ie++) {
            strcpy(ename, sct_get_event_cname(ie));
            if (sct_get_eventcounters() == 1) {
              a = x->esum[ie] * q_reduce;
              e = (x->emax[ie] > eps) ? a/x->emax[ie] : 1.0;
              e = e*100; // percent
	      if ( entity == -1)
                fprintf(outstream," %*s%-*s |%10.2e%10.2e%10.2e%10.2e  %7.3f\n",
                        EVENTINDENT,".", mnl-EVENTINDENT, ename, (double)x->emin[ie], a,
                        (double)x->emax[ie], (double)x->esum[ie], e);
	      else
                fprintf(outstream," %*s%-*s |          |%10.2e%10.2e%10.2e%10.2e  %7.3f\n",
                        EVENTINDENT,".", mnl-EVENTINDENT, ename, (double)x->emin[ie], a,
                        (double)x->emax[ie], (double)x->esum[ie], e);
            } else {
              a = x->esum[ie] / x->tsum;
              e = (x->rmax[ie] > eps) ? a/x->rmax[ie] : 1.0;
              e = e*100; // percent
              strcat(ename," (rate)");
	      if ( entity == -1)
                fprintf(outstream," %*s%-*s |%10.2e%10.2e%10.2e%10s  %7.3f\n",
                        EVENTINDENT,".", mnl-EVENTINDENT, ename, (double)x->rmin[ie], a,
                        (double)x->rmax[ie], " ", e);
	      else
                fprintf(outstream," %*s%-*s |          |%10.2e%10.2e%10.2e%10s  %7.3f\n",
                        EVENTINDENT,".", mnl-EVENTINDENT, ename, (double)x->rmin[ie], a,
                        (double)x->rmax[ie], " ", e);
            }
	  }
#endif
        }
      }
    }

  }
  } /* if no reduction ... else ...*/

  if (hsep) free(hsep);
  if (hsep_proc) free(hsep_proc);

}


static void time_sec_str(const double t, char *s, const int sn) {

  if (sn<12) sct_abort("internal error - string size too small", __FILE__, __LINE__);

  if (t < 0.0) {
    sprintf(s, "%10s","    ??????");
  } else if (t < 5.0e-5) {
    sprintf(s, "    0.0000");
  } else if (t < 1.e1) {
    sprintf(s, "%10.4f", t);
    //printf("time_sec_str: exact_t=%.13e t=%10.4f, sn=%i, [s]=[%s]\n",t, t,sn,s);
  } else if (t < 1.e2) {
    sprintf(s, "%10.3f", t);
  } else if (t < 1.e3) {
    sprintf(s, "%10.2f", t);
  } else if (t < 1.e4) {
    sprintf(s, "%10.1f", t);
  } else {
    sprintf(s, "%10.0f", t);
  }
}


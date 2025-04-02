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
 * \file sct_mach.c
 * abstraction of machine dependencies, including special libraries
 */
#include "sct_mach.h"

#include <stdio.h>
#include <stdlib.h>
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
#include <unistd.h>

#ifdef _AIX
#if HAVE_TIME_H
#  include <time.h>
#endif
#endif

#if (defined (__linux) && (defined(__x86_64__) || defined(__i686__) || defined(__i386__)))
#if HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif
#endif

#ifdef _AIX
/* default pmapi counter group 127 for POWER6 ... FPU flop events */
static char default_event_list[] =
  "PAPI_FP_OPS,PAPI_FMA_INS,PAPI_TOT_CYC";
//  "PM_FPU_FLOP,PM_RUN_INST_CMPL,PM_RUN_CYC";
#elif (defined (__linux) && (defined(__GNUC__) || defined(__ICC) || defined(__SUNPRO_C) || defined(__PGIC__)) && (defined(__x86_64__)||defined(__i686__)||defined(__i386__)))
/* default PAPI FPU flop events */
static char default_event_list[] =
  "PAPI_LD_INS,PAPI_TOT_INS,PAPI_TOT_CYC";
#elif (defined(__APPLE__) && defined(__GNUC__) && defined(__x86_64__))
static char default_event_list[] = "";
#else
#  error unknown system
#endif

// estimated frequency of the system
static double mhz_inv = 1.;

static int on_err(const char *funstr, const int errcode, const char *fname, const int line);

#ifdef HAVE_LIBPAPI
#ifdef _OPENMP
static unsigned long wrap_omp_get_thread_num(void);
#endif
static int get_event_names(char ***event_name, int *event_num);
static int gen_eset_id(eset_type *eset);
#endif

// reference point for absolute time
#if (HAVE_TIME_H || HAVE_SYS_TIME_H)
static struct timeval prg_start_time;
#endif
static char prg_start_time_string[40];

// implementation:
void sct_mach_init() {
#if ( SCT_RTM & SCT_RTM_RDTSCP )
  // estimate mhz by 1 sec sleep
  int mhz = 0;
  tmark_type tm1;

  tm1 = rdtscp();
  sleep(1);
  mhz = (int) ((rdtscp() - tm1)/1000000);
  mhz_inv = 1. / ((double)mhz);

  fprintf(stderr,"(sct_mach_init) SCT_RTM = %i using mhz = %d\n", SCT_RTM, mhz);
#endif
}

void set_prg_start_time() {
#if (HAVE_TIME_H || HAVE_SYS_TIME_H)
  struct tm * ptm;
  gettimeofday(&prg_start_time, NULL);
  ptm = localtime(&prg_start_time.tv_sec);
  strftime(prg_start_time_string, sizeof(prg_start_time_string), "%Y-%m-%d %H:%M:%S", ptm);
#else
strcpy(prg_start_time_string,"not available");
#endif
}

const char* print_prg_start_time() {
  return prg_start_time_string;
}

inline double get_tdiff(tmark_type *ptm1, tmark_type *ptm2) {
  double dt;

  // time diff in secs
#if (SCT_RTM == SCT_RTM_READ_REAL_TIME)

  int secs, n_secs;
  // man page recommends to call the conversion routines unconditionally, so we do
  time_base_to_time(ptm1, TIMEBASE_SZ);
  time_base_to_time(ptm2, TIMEBASE_SZ);

  secs   = ptm2->tb_high - ptm1->tb_high;
  n_secs = ptm2->tb_low  - ptm1->tb_low;

  dt = (double)secs + n_secs * 1.e-9;

#elif ( SCT_RTM & (SCT_RTM_OMP_GET_WTIME | SCT_RTM_MPI_WTIME) )

  dt = *ptm2 - *ptm1;

#elif ( SCT_RTM & (SCT_RTM_CLOCK_GETTIME_MONOTONIC | SCT_RTM_CLOCK_GETTIME_REALTIME) )

  dt = (double) (ptm2->tv_sec - ptm1->tv_sec) + 1.e-9 * (double) (ptm2->tv_nsec - ptm1->tv_nsec);

#elif ( SCT_RTM & SCT_RTM_GETTIMEOFDAY )

  dt = (double) (ptm2->tv_sec - ptm1->tv_sec) + 1.e-6 * (double) (ptm2->tv_usec - ptm1->tv_usec);

#elif ( SCT_RTM & SCT_RTM_RDTSCP )

  dt = (double)(*ptm2 - *ptm1) * mhz_inv * 1.e-6;

#else
#error unknown system
#endif
  return dt;
}


void sct_abort(const char *reason, const char *fname, const int line) {

  if (reason||fname) fprintf(stderr, "[sct] ERROR in file %s, line %d: %s\n",fname, line, reason);

  fflush(stdout);
  fflush(stderr);

#ifdef HAVE_MPI
  int errorcode;
  errorcode = 1;
  MPI_Abort(MPI_COMM_WORLD, errorcode);
#endif

  abort();
}


void sct_warn(const char *reason, const char *fname, const int line) {

  if (reason||fname) fprintf(stderr, "[sct] WARNING in file %s, line %d: %s\n",fname, line, reason);
}


#ifdef HAVE_LIBPAPI
eset_type *sct_new_eventset() {
  eset_type *eset;

  if ( !(eset = calloc(1,sizeof(eset_type))) ) {
    fprintf(stderr,"sct_new_eventdef: calloc failed\n");
    return NULL;
  }

  get_event_names(&eset->name, &eset->en);

  if (eset->en < 1) {
    free(eset);
    fprintf(stderr,"sct_new_eventdef: eset->en < 1\n");
    return NULL;
  }

  eset->valid = 0;
  eset->tn = 0;
  eset->id = NULL;

  if (!gen_eset_id(eset)) {
    sct_del_eventset(eset);
    fprintf(stderr,"sct_new_eventdef: !gen_eset_id(eset)\n");
    return NULL;
  }

  return eset;
}
#endif


#ifdef HAVE_LIBPAPI
void sct_del_eventset(eset_type *eset) {

  if (!eset) return;
  
  for (int ie=eset->en-1; ie>=0; ie--) {
    if (eset->name[ie]) free(eset->name[ie]);
  }

  if (eset->id) free(eset->id);

  free(eset);
}
#endif


#ifdef HAVE_LIBPAPI
int sct_start_eventset(eset_type *eset) {

  if (!eset) return 1;

  int ierr;
  int istate = 0;
  int err = 0;

  if (omp_in_parallel()) sct_abort("sct_start_eventset called in parallel region", __FILE__, __LINE__);

  // start counting; returns true if successful

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int tid = omp_get_thread_num();

#ifdef _OPENMP
#pragma omp critical
#endif
    {
#ifdef DEBUG
      printf("start_eventset: tid=%d, id[tid] = %d, en = %d\n",tid, eset->id[tid], eset->en);
#endif
      if ( (ierr = PAPI_start(eset->id[tid])) != PAPI_OK ) err++;
    }
  }

  if (err) return on_err("PAPI_start", ierr, __FILE__, __LINE__);
  
  return 1;

}
#endif


#ifdef HAVE_LIBPAPI
int sct_stop_eventset(eset_type *eset) {

  if (!eset) return 1;

  int ierr;
  eval_type eval[eset->en];

  if (omp_in_parallel()) sct_abort("sct_stop_eventset called in parallel region", __FILE__, __LINE__);

  // stop counting; returns true if successful
  // we throw the end values away

  for (int tid=0; tid< eset->tn; tid++) {
    if ( (ierr = PAPI_stop(eset->id[tid], eval)) != PAPI_OK )
      return on_err("PAPI_stop", ierr, __FILE__, __LINE__);
  }

  return 1;
}
#endif


#ifdef HAVE_LIBPAPI
#ifdef _OPENMP
  static unsigned long wrap_omp_get_thread_num(void) {
    return (unsigned long)omp_get_thread_num();
  }
#endif
#endif


#ifdef HAVE_LIBPAPI
static int get_event_names(char ***event_name, int *event_num) {

  // returns number of successfully parsed tokens
  char *evlist;
  char *t, *w;
  int tn, tlen;
  char *evnames[10];

  *event_num = 0;

  if ( !(evlist = getenv ("SCT_EVENT_LIST")) ) evlist = default_event_list;
  tn = 0;
  t = strtok(evlist," ,:;\t\r\n");
  while(t) {
    if (tn>=10) {
      fprintf(stderr, "ERROR: event_list has too many elements (resticted to 10)\n");
      return 0;
    }
    tlen = strlen(t);
    w = (char*)malloc(tlen+1);
    strncpy(w, t, tlen);
    w[tlen]='\0';
    evnames[tn] = w;
    tn++;
    t = strtok(NULL," ,:;");
  }

  if (!tn) {
    fprintf(stderr, "ERROR: empty event_list\n");
    return 0;
  }

  // now update passed event_name and event_num
  *event_num = tn;

  if ( !(*event_name = calloc(tn, sizeof(char*))) )
    sct_abort("memory allocation for event_names failed", __FILE__, __LINE__);

  for (int i=0; i<tn; i++) (*event_name)[i] = evnames[i];

  return *event_num;
}
#endif



static int on_err(const char *funstr, const int errcode, const char *fname, const int line) {

  // error handler; always returns false
#ifdef HAVE_LIBPAPI
  fprintf(stderr, "PAPI ERROR in file %s line %d: %s failed. %d %s\n",
	  fname, line, funstr, errcode, PAPI_strerror(errcode));
#else
  fprintf(stderr, "PAPI ERROR in file %s line %d: %s failed. %d\n",
	  fname, line, funstr, errcode);

#endif
  return 0;
}


#ifdef HAVE_LIBPAPI
static int gen_eset_id(eset_type *eset) {

  if (!eset) return 1;

  const PAPI_hw_info_t *hw_info = NULL;

  // initialization of counters; returns true if successful
  int ierr, code, retcode;

  if ( (ierr = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT)
    return on_err("PAPI_library_init", ierr, __FILE__, __LINE__);

#ifdef _OPENMP
  if ( sizeof(unsigned long) == sizeof(int) ) { 
    ierr = PAPI_thread_init((unsigned long (*)(void))(omp_get_thread_num));
  } else {
    fprintf(stderr, "sct - warning: installing wrapper for omp_get_thread_num\n");
    ierr = PAPI_thread_init(wrap_omp_get_thread_num);
  }
  if ( ierr != PAPI_OK ) return on_err("PAPI_thread_init", ierr, __FILE__, __LINE__);
#endif

  // debug setting
#ifdef DEBUG
  if ( (ierr = PAPI_set_debug(PAPI_VERB_ESTOP)) != PAPI_OK )
    return on_err("PAPI_library_init", ierr, __FILE__, __LINE__);
#else
  if ( (ierr = PAPI_set_debug(PAPI_QUIET)) != PAPI_OK )
    return on_err("PAPI_library_init", ierr, __FILE__, __LINE__);
#endif

#ifdef DEBUG
  if ( !(hw_info = PAPI_get_hardware_info()) ) return on_err("PAPI_get_hardware_info", ierr, __FILE__, __LINE__);
  printf("PAPI_get_hardware_info: %d ncpu, %d nnodes, %d totalcpus, %s vendor, %f MHz\n",
         hw_info->ncpu, hw_info->nnodes, hw_info->totalcpus, hw_info->vendor_string, hw_info->mhz);
#endif

  eset->tn = omp_get_max_threads();

  if ( !(eset->id  = calloc(eset->tn, sizeof(int))) ) sct_abort("calloc failed", __FILE__, __LINE__);;

#ifdef _OPENMP
#pragma omp parallel private(ierr) reduction(+:retcode)
#endif
  {
    int tid = omp_get_thread_num();

    eset->id[tid] = PAPI_NULL;
    retcode = 1; // 1==OK, 0==failure
    if ( (ierr = PAPI_create_eventset(&eset->id[tid])) != PAPI_OK ) {
#ifdef DEBUG
      printf("PAPI_create_eventset failed for tid [%d]\n",tid);
#endif
      retcode = on_err("PAPI_create_eventset", ierr, __FILE__, __LINE__);
    }
    else {
      for (int ie=0; ie<eset->en; ie++) {
	if ( (ierr = PAPI_event_name_to_code(eset->name[ie], &code)) != PAPI_OK ) {
#ifdef DEBUG
          printf("PAPI_event_name_to_code failed for event [%s]\n",eset->name[ie]);
#endif
          retcode = on_err("PAPI_event_name_to_code", ierr, __FILE__, __LINE__);
          break;
	}
	if ( (ierr = PAPI_add_event(eset->id[tid], code)) != PAPI_OK ) {
#ifdef DEBUG
          printf("PAPI_add_event failed for code [%d]\n",code);
#endif
          retcode = on_err("PAPI_add_event", ierr, __FILE__, __LINE__);
	  break;
	}
      }
    }
  }
#ifdef DEBUG
  printf("PAPI event init finished with [%d]\n",retcode);
#endif
  if (retcode == eset->tn) eset->valid = 1;
  return retcode;
}
#endif


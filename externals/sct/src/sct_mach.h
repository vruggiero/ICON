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

/*
 * sct_mach.h: abstraction of machine dependencies, including special libraries
 */

#ifndef _H_SCT_MACH
#define _H_SCT_MACH

#if HAVE_CONFIG_H
#  ifndef _H_CONFIG
#    define _H_CONFIG
#    include <config.h>
#  endif
#endif

#if HAVE_STDIO_H
#  include <stdio.h>
#endif

// system includes:
#if HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif
#if HAVE_TIME_H
#  include <time.h>
#endif

#if (defined(_AIX) && defined(HAVE_SYS_SYSTEMCFG_H))
#  include <sys/systemcfg.h>
#elif (defined (__linux) && (defined(__x86_64__) || defined(__i686__) || defined(__i386__)))
#  if HAVE_STDINT_H
#    include <stdint.h>
#  endif
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  endif
#elif (defined(__APPLE__) && defined(__x86_64__))
#  if HAVE_STDINT_H
#    include <stdint.h>
#  endif
#  if HAVE_INTTYPES_H
#    include <inttypes.h>
#  endif
#else
#  error unknown system
#endif

//config:
#include "sct_config.h"

void sct_mach_init();

//public+

// OMP:
#ifdef _OPENMP
#  include <omp.h>
#else
/*! \brief fall back definition if we don't have OpenMP:
           we are never inside parallel region */
static int inline omp_in_parallel() { return 0; }
/*! \brief fall back definition if we don't have OpenMP:
           only thread id 0 is present */
static int inline omp_get_thread_num() { return 0; }
/*! \brief fall back definition if we don't have OpenMP:
           only 1 thread is defined */
static int inline omp_get_max_threads() { return 1; }
#endif

// MPI:
#ifdef HAVE_MPI
#  include <mpi.h>
#endif
//public-

// PAPI:
#ifdef HAVE_LIBPAPI
#  include "papi.h"
typedef long long eval_type;
#  define SCT_MAX_EVAL LLONG_MAX
#endif

// time mark type
#if   (SCT_RTM == SCT_RTM_READ_REAL_TIME)

typedef  timebasestruct_t tmark_type;

#elif ( SCT_RTM & (SCT_RTM_OMP_GET_WTIME | SCT_RTM_MPI_WTIME) )

typedef double tmark_type;

#elif ( SCT_RTM & (SCT_RTM_CLOCK_GETTIME_MONOTONIC | SCT_RTM_CLOCK_GETTIME_REALTIME) )

typedef struct timespec tmark_type;

#elif ( SCT_RTM & SCT_RTM_GETTIMEOFDAY )

typedef struct timeval tmark_type;

#elif ( SCT_RTM & SCT_RTM_RDTSCP )

typedef uint64_t tmark_type;

#else
#error unknown read time method
#endif

static inline int equal_tmarks(tmark_type *ptm1, tmark_type *ptm2) {
#if   (SCT_RTM == SCT_RTM_READ_REAL_TIME)

  return ( (ptm1->tb_high == ptm2->tb_high) && (ptm1->tb_low == ptm2->tb_low) );

#elif ( SCT_RTM & (SCT_RTM_OMP_GET_WTIME | SCT_RTM_MPI_WTIME) )

  return ( *ptm1 == *ptm2 );

#elif ( SCT_RTM & (SCT_RTM_CLOCK_GETTIME_MONOTONIC | SCT_RTM_CLOCK_GETTIME_REALTIME ) )

  return ( (ptm1->tv_sec == ptm2->tv_sec) && (ptm1->tv_nsec == ptm2->tv_nsec) );

#elif ( SCT_RTM & SCT_RTM_GETTIMEOFDAY )

  return ( (ptm1->tv_sec == ptm2->tv_sec) && (ptm1->tv_usec == ptm2->tv_usec) );

#elif ( SCT_RTM & SCT_RTM_RDTSCP )

  return ( *ptm1 == *ptm2 );

#else
#error unknown read time method
#endif
}

#if ( SCT_RTM & SCT_RTM_RDTSCP )

static inline uint64_t rdtscp(void) {
  uint32_t low, high;

  /* rdtscp changes RAX,RCX and RDX */
  uint64_t rax, rcx, rdx;
  __asm__ __volatile__( ""
  			: "=a" (rax), "=c" (rcx), "=d" (rdx)/* output*/
  			: /* input */
  			:);

#ifdef DEBUG
  fprintf(stderr,"pre: rax 0x%016"PRIx64"\n",myrax);
  fprintf(stderr,"pre: rcx 0x%016"PRIx64"\n",myrcx);
  fprintf(stderr,"pre: rdx 0x%016"PRIx64"\n",myrdx);
#endif

  __asm__ __volatile__ ("rdtscp"
                        : "=a" (low), "=d" (high));

  /* rewrite RAX,RCX and RDX */
  __asm__ __volatile__ (""
                        :
  			: "a" (rax), "c" (rcx), "d" (rdx)
  			:);

#ifdef DEBUG
  fprintf(stderr,"post: rax 0x%016"PRIx64"\n",myrax);
  fprintf(stderr,"post: rcx 0x%016"PRIx64"\n",myrcx);
  fprintf(stderr,"post: rdx 0x%016"PRIx64"\n",myrdx);
#endif

  return ((low) | ((uint64_t)(high) << 32));
}
#endif


// event definition
typedef struct {
  int valid;
  int *id;      // PAPI-id; this is an array with one id per thread
  int tn;       // number of threads
  int en;       // number of events
  char **name;  // array of event-names (size == en)
} eset_type;

// attributes:
#if (defined(__linux) && (defined(__GNUC__) || defined(__PGI) || defined(__ICC) || defined(__SUNPRO_C) || defined(__PGIC__)) && (defined(__x86_64__)||defined(__i686__)))
#   define ATTRIBUTE_NORETURN __attribute__ ((noreturn))
#else
#   define ATTRIBUTE_NORETURN
#endif

// not time critical; implemented in sct_mach.c:
#ifdef HAVE_LIBPAPI
eset_type *sct_new_eventset();
void sct_del_eventset(eset_type *eset);
int sct_start_eventset(eset_type *eset);
int sct_stop_eventset(eset_type *eset);
#endif
void sct_abort(const char *reason, const char *fname, const int line) ATTRIBUTE_NORETURN;
void sct_warn(const char *reason, const char *fname, const int line);

static inline void read_time(tmark_type *tm) {
#if (SCT_RTM == SCT_RTM_READ_REAL_TIME)
  read_real_time(tm, TIMEBASE_SZ);
#elif (SCT_RTM == SCT_RTM_OMP_GET_WTIME)
  *tm = omp_get_wtime();
#elif (SCT_RTM == SCT_RTM_MPI_WTIME)
  *tm = MPI_Wtime();
#elif (SCT_RTM == SCT_RTM_CLOCK_GETTIME_MONOTONIC)
  clock_gettime(CLOCK_MONOTONIC, tm);
#elif (SCT_RTM == SCT_RTM_CLOCK_GETTIME_REALTIME)
  clock_gettime(CLOCK_REALTIME, tm);
#elif (SCT_RTM == SCT_RTM_GETTIMEOFDAY)
  gettimeofday(tm, NULL);
#elif (SCT_RTM == SCT_RTM_RDTSCP)
  *tm = rdtscp();
#else
#error unknown read time method
#endif
}

void set_prg_start_time();
const char* print_prg_start_time();

#ifdef HAVE_LIBPAPI
static inline int sct_read_events(const int tid, const eset_type *eset, eval_type *eval) {
  int err;
  // returns true on success
  if ( (err = PAPI_read(eset->id[tid], eval)) != PAPI_OK) {
#ifdef DEBUG
    fprintf(stderr,"sct_read_events failed, err = %d\n", err);
#endif
    return 0;
  } else {
#ifdef DEBUG
    fprintf(stderr,"sct_read_events: eset->id = %i, eval[0] = %lld\n", *eset->id, eval[0]);
#endif
  return 1;
  }
}
#endif

// declarations:
double get_tdiff(tmark_type *ptm1, tmark_type *ptm2);

#ifdef HAVE_LIBPAPI
static inline void get_ediff(const int evn, const eval_type *evec1, const eval_type *evec2, eval_type *ediff) {
  int ie;

  for (ie = 0; ie < evn; ie++) {
    ediff[ie] = evec2[ie] - evec1[ie];
  }
}
#endif

static inline void update_tsum(tmark_type *ptm1, double *tsum) {
  tmark_type tm2;

  read_time(&tm2);
  *tsum += get_tdiff(ptm1, &tm2);
}

#endif

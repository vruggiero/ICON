#if HAVE_CONFIG_H
#  ifndef _H_CONFIG
#    define _H_CONFIG
#    include <config.h>
#  endif
#endif

#if HAVE_STDLIB_H
#  include <stdlib.h>
#endif
#if HAVE_STRING_H
#  include <string.h>
#endif
#if HAVE_STRINGS_H
#  include <strings.h>
#endif
#if HAVE_MATH_H
#  include <math.h>
#endif
#if _OPENMP
#  include <omp.h>
#endif
#if HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif
#if HAVE_TIME_H
#  include <time.h>
#endif

#include "sct.h"

#define INDEX 100

#define MY_TIMER_MAX 4
static const int timer_max = MY_TIMER_MAX;
static int timer[MY_TIMER_MAX];

static int pid, tid, nbt, istart, iend;
static double blocksize;
double **matrixa;
double **matrixb;
static double **mresult;

#pragma omp threadprivate(istart, iend, blocksize, tid, nbt, mresult)

void
dummy( void *array )
{
/* Confuse the compiler so as not to optimize
   away the flops in the calling routine    */
/* Cast the array as a void to eliminate unused argument warning */
        ( void ) array;
}


/* trying to implement rdtsc for each architecture :-)
   to determine actual cpu-frequency */

#if defined(__i386__)

static __inline__ unsigned long long int rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}

#elif defined(__x86_64__)

static __inline__ unsigned long long int rdtsc(void)
{
  unsigned int x = 0, y = 0;
  __asm__ __volatile__ ("rdtsc" : "=a"(x), "=d"(y));
  return ( ((unsigned long)x) | (((unsigned long)y)<<32) );
}

#elif defined(__powerpc__)

static __inline__ unsigned long long int rdtsc(void)
{
  timebasestruct_t t;
  unsigned long long int retval;
  int type;

  type = read_real_time(&t, TIMEBASE_SZ);
  time_base_to_time(&t, TIMEBASE_SZ);
  retval = t.tb_high;
  retval = retval<<32;
  retval = retval|t.tb_low;

  return (retval);
}

/* static __inline__ unsigned long long int rdtsc(void) */
/* { */
/*   unsigned long long int result=0; */
/*   unsigned long int upper, lower,tmp; */
/*   __asm__ volatile( */
/*                 "loop:                  \n" */
/*                 "\tmftbu   %0           \n" */
/*                 "\tmftb    %1           \n" */
/*                 "\tmftbu   %2           \n" */
/*                 "\tcmpw    %2,%0        \n" */
/*                 "\tbne     loop         \n" */
/*                 : "=r"(upper),"=r"(lower),"=r"(tmp) */
/*                 ); */
/*   result = upper; */
/*   result = result<<32; */
/*   result = result|lower; */

/*   return(result); */
/* } */

#else

#error "No tick counter is available!"

#endif


static __inline__ double getTime(void) {
  struct timeval tbuf;

  gettimeofday(&tbuf, NULL);

  return ( (double) tbuf.tv_sec + (double) tbuf.tv_usec*1.e-6 );
}


double determineFrequency(void) {
  double start, elapsed, accum=0.0, y;
  int i, flipper=1;
  unsigned long long int x;

  x = rdtsc();
  start = getTime();
  
  while ( (elapsed=getTime()-start) < 1 ) {
    if (flipper == 1) flipper = -1; else flipper = 1;
    for (i=0; i<1000000; i++) {
      accum = accum + (i * (double) flipper);
    }
  }

  x = rdtsc() - x;
  elapsed = getTime() - start;

  y = (double) x/elapsed;

  printf("outputcheck: cpu frequency estimation = %1.3f GHz\n", 1e-9*y);
  return y;
}


void initMatrix(int mysize) {
  int i,j,k;

  matrixa = (double **) malloc(mysize*sizeof(double *));
  matrixb = (double **) malloc(mysize*sizeof(double *));
  for (i=0; i<mysize; i++) {
    matrixa[i] = (double *) malloc(mysize*sizeof(double));
    matrixb[i] = (double *) malloc(mysize*sizeof(double));
  }

  for ( i = 0; i < mysize; i++ ) {
    for ( j = 0; j < mysize; j++) {
      matrixa[i][j] = ( float ) rand(  ) * ( float ) 1.1;
      matrixb[i][j] = ( float ) rand(  ) * ( float ) 1.1;
    }
  }

#pragma omp parallel private(i,j)
  {
#if _OPENMP
    tid = omp_get_thread_num();
    nbt = omp_get_num_threads();
#else
    tid = 0;
    nbt = 1;
#endif
    iend = 0;
    blocksize = 0.0;
    /* find some blocksize which produces inbalance in threaded region */
    for (i=0; i<nbt; i++) blocksize += 1./pow(1.5,i);
    blocksize = mysize / blocksize;
    for (i=0; i<=tid; i++) {
      istart = iend;
      iend += (int)(blocksize/pow(1.5,i));
    }
    if (tid == nbt-1) iend = mysize;

    mresult = (double **) malloc((iend-istart)*sizeof(double *));
    for (i=istart; i<iend; i++)
      mresult[i-istart] = (double *) malloc(mysize*sizeof(double));

    for (int i=istart; i<iend; i++) {
      for (int j=0; j<mysize; j++) {
        mresult[i-istart][j] = 0.0;
      }
    }
  }
}


void freeMatrix(int mysize) {
  int i;

  for (i=0; i<mysize; i++) {
    free(matrixa[i]);
    free(matrixb[i]);
  }
  free(matrixa);
  free(matrixb);
#pragma omp parallel private(i)
  {
    for (i=istart; i<iend; i++) free(mresult[i-istart]);
    free(mresult);
  }
}


int main( int argc, char **argv ) {
  int i, j, k;
  int myINDEX;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#else
  pid = 0;
#endif

  /* init SCT*/
  sct_init(timer_max, "test_papi", SCT_COMM_WORLD);
  for (i=0; i < timer_max; i++) timer[i] = 0;
  char label0[] = "sct-timer 0";
  char label1[] = "sct-timer 1";
  char label2[] = "seq-part";
  timer[0] = sct_new_timer(label0);
  timer[1] = sct_new_timer(label1);
  timer[2] = sct_new_timer(label2);

  /* Initialize the Matrix arrays for round 0*/
  myINDEX = (pid+1)*INDEX;
  initMatrix(myINDEX);

#if _OPENMP
  sct_start(timer[2]);
#endif

#pragma omp parallel private(i,j,k)
  {
    sct_start(timer[0]);

    /* Matrix-Matrix multiply */
    for (i=istart; i<iend; i++)
      for (j=0; j<myINDEX; j++)
        for (k=0; k<myINDEX; k++)
          mresult[i-istart][j] = mresult[i-istart][j] + matrixa[i][k] * matrixb[k][j];

    dummy( ( void * ) mresult );

    sct_stop(timer[0]);
  }

#if _OPENMP
  sct_stop(timer[2]);
#endif

  /* free memory */
  freeMatrix(myINDEX);

  /* Initialize the Matrix arrays for round 1*/
  myINDEX = (pid+1)*INDEX*2;
  initMatrix(myINDEX);

#if _OPENMP
  sct_start(timer[2]);
#endif

#pragma omp parallel private(i,j,k)
  {
    sct_start(timer[1]);

    /* Matrix-Matrix multiply */
    for (i=istart; i<iend; i++)
      for (j=0; j<myINDEX; j++)
        for (k=0; k<myINDEX; k++)
          mresult[i-istart][j] = mresult[i-istart][j] + matrixa[i][k] * matrixb[k][j];

    dummy( ( void * ) mresult );

    sct_stop(timer[1]);
  }

#if _OPENMP
  sct_stop(timer[2]);
#endif

  /* free memory */
  freeMatrix(myINDEX);

  sct_report(SCT_GETENV, SCT_GETENV, SCT_GETENV);
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
#endif

  int evc;
  char *s = getenv ("SCT_EVENTCOUNTERS");
#ifdef HAVE_MPI
  if (!s) evc = 1;
#else
  if (!s) evc = 2;
#endif
  else {
    if ( !strcmp(s,"0") || !strcasecmp(s,"SCT_WITHOUT_EVENTCOUNTERS") ) evc = 0;
    else {
      if ( !strcmp(s,"1") || !strcasecmp(s,"SCT_WITH_EVENTCOUNTERS") ) evc = 1;
      else evc = 2;
    }
  }

  // check for tests that might have failed
  int err = 0;
  double act0, act1, ratio;
  if (evc == 1) {
#if _OPENMP
#pragma omp parallel private(act0, act1, ratio)
    {
#pragma omp critical
      {
	/* if FLOPS are counted correct */
        /* ref = (double)2*(iend-istart)*myINDEX*myINDEX; */
        /* act = sct_event(timer[0], "PAPI_FP_OPS"); */
        /* diff = fabs(ref-act)/ref*100.; */
        /* printf("outputcheck: PAPI_FP_OPS for proc %i, thread %i: ref %10.2e act %10.2e diff %5.2f %%\n", */
        /*        pid, tid, ref, act, diff); */
        /* printf("             criterion (diff < 1.) ... is allowed to fail if fma is used !\n"); */
	/* if (diff > 1.) err = 1; */

	/* else compare 'Load instructions' for both runs of matmult*/
        act0 = sct_event(timer[0], "PAPI_LD_INS");
        act1 = sct_event(timer[1], "PAPI_LD_INS");
        ratio = act1/act0;
        printf("outputcheck: PAPI_LD_INS ratio for proc %i, thread %i: act0 %10.2e act1 %10.2e ratio %5.2f\n",
               pid, tid, act0, act1, ratio);
	printf("             criterion (7 <= ratio <= 8)\n");
	if ((ratio > 8.) || (ratio < 7.)) err = 1;
      }
    }
#else
    /* if FLOPS are counted correct */
    /* ref = (double)2*(iend-istart)*myINDEX*myINDEX; */
    /* act = sct_event(timer[0], "PAPI_FP_OPS"); */
    /* diff = fabs(ref-act)/ref*100.; */
    /* printf("outputcheck: PAPI_FP_OPS for proc %i: ref %10.2e act %10.2e diff %5.2f %%\n", pid, ref, act, diff); */
    /* printf("             criterion (diff < 1.) ... is allowed to fail if fma is used !\n"); */
    /* if ( diff > 1. ) err = 1; */

    /* else compare 'Load instructions' for both runs of matmult*/
    act0 = sct_event(timer[0], "PAPI_LD_INS");
    act1 = sct_event(timer[1], "PAPI_LD_INS");
    ratio = act1/act0;
    printf("outputcheck: PAPI_LD_INS ratio for proc %i: act0 %10.2e act1 %10.2e ratio %5.2f\n",
	   pid, act0, act1, ratio);
    printf("             criterion (7 <= ratio <= 8)\n");
    if ((ratio > 8.) || (ratio < 7.)) err = 1;
#endif
  }
  /* skip check for clock rate ... not accurate enough on HSW and BDW ...
  else if (evc == 2 && pid == 0) {
    ref = determineFrequency();
#if _OPENMP
    act = sct_event(timer[1], "PAPI_TOT_CYC");
#else
    act = sct_event(timer[0], "PAPI_TOT_CYC");
#endif
    diff = fabs(ref-act)/ref*100.;
    printf("outputcheck: PAPI_TOT_CYC (rate) for all procs and threads: ref %10.2e act %10.2e diff %5.2f %%\n", ref, act, diff);
    printf("             criterion (diff < 25.)\n");
    if ( diff > 25. ) err = 1;
  }
  */

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return err;
}

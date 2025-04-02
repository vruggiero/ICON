/**
 * \file test_hdf5.c
 * A simple test whether hdf5 file has correct values.
 *
 * \copyright Copyright  (C)  2016 Hendryk Bockelmann <bockelmann@dkrz.de>
 *
 * \author Hendryk Bockelmann <bockelmann@dkrz.de>
 */

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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "hdf5.h"
#include "sct.h"

#define INDEX 100

#define MY_TIMER_MAX 3
static const int timer_max = MY_TIMER_MAX;
static int timer[MY_TIMER_MAX];

void
dummy( void *array )
{
/* Confuse the compiler so as not to optimize
   away the flops in the calling routine    */
/* Cast the array as a void to eliminate unused argument warning */
        ( void ) array;
}


int main( int argc, char **argv ) {
  int i, j, k, l;
  static int pid, tid, nprocs, nbthreads, istart, iend;
  static double blocksize, freq;
  int myINDEX;
  double **matrixa;
  double **matrixb;
  static double **mresult;

#pragma omp threadprivate(istart, iend, blocksize, tid, nbthreads, mresult)

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
  pid = 0;
  nprocs = 1;
#endif

  /* init SCT */
  sct_init(timer_max, "test_hdf5", SCT_COMM_WORLD);
  for (i=0; i < timer_max; i++) timer[i] = 0;
  char label0[] = "matmult";
  char label1[] = "seq-part";
  char label2[] = "init";
  timer[0] = sct_new_timer(label0);
  timer[1] = sct_new_timer(label1);
  timer[2] = sct_new_timer(label2);

  /* Initialize the Matrix arrays */
  myINDEX = (pid+1)*INDEX;

  sct_start(timer[2]);
  matrixa = (double **) calloc(myINDEX,sizeof(double *));
  matrixb = (double **) calloc(myINDEX,sizeof(double *));
  for (i=0; i<myINDEX; i++) {
    matrixa[i] = (double *) calloc(myINDEX,sizeof(double));
    matrixb[i] = (double *) calloc(myINDEX,sizeof(double));
  }

  for ( i = 0; i < myINDEX; i++ ) {
    for ( j = 0; j < myINDEX; j++) {
      matrixa[i][j] = ( float ) rand(  ) * ( float ) 1.1;
      matrixb[i][j] = ( float ) rand(  ) * ( float ) 1.1;
    }
  }
  sct_stop(timer[2]);

#pragma omp parallel private(i,j)
  {
#if _OPENMP
    tid = omp_get_thread_num();
    nbthreads = omp_get_num_threads();
#else
    tid = 0;
    nbthreads = 1;
#endif
    iend = 0;
    blocksize = 0.0;

    sct_start(timer[2]);

    /* find some blocksize which produces inbalance in threaded region */
    for (i=0; i<nbthreads; i++) blocksize += 1./pow(1.5,i);
    blocksize = myINDEX / blocksize;
    for (i=0; i<=tid; i++) {
      istart = iend;
      iend += (int)(blocksize/pow(1.5,i));
    }
    if (tid == nbthreads-1) iend = myINDEX;

    mresult = (double **) calloc((iend-istart),sizeof(double *));
    for (i=istart; i<iend; i++)
      mresult[i-istart] = (double *) calloc(myINDEX,sizeof(double));

    for (int i=istart; i<iend; i++) {
      for (int j=0; j<myINDEX; j++) {
        mresult[i-istart][j] = 0.0;
      }
    }

    sct_stop(timer[2]);
  }

#if _OPENMP
  sct_start(timer[1]);
#endif

#pragma omp parallel private(i,j,k,l)
  {
    for (l=0; l<=tid; l++) {
      sct_start(timer[0]);
      /* Matrix-Matrix multiply */
      for (i=istart; i<iend; i++)
	for (j=0; j<myINDEX; j++)
	  for (k=0; k<myINDEX; k++)
	    mresult[i-istart][j] = mresult[i-istart][j] + matrixa[i][k] * matrixb[k][j];
      
      dummy( ( void * ) mresult );
      
      sct_stop(timer[0]);
    }
  }

#if _OPENMP
  sct_stop(timer[1]);
#endif

  /* free memory */
  for (i=0; i<myINDEX; i++) {
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

  /* specify some dummy attributes (key-value pairs) for timer report */
#ifdef HAVE_C__GENERIC
  sct_add_report_attribute("fooINT", 42);
  sct_add_report_attribute("fooLONG", (long)42);
  sct_add_report_attribute("fooFLOAT", (float)15.7);
  sct_add_report_attribute("fooDOUBLE", 15.7);
  sct_add_report_attribute("bar", "test");
#else
  sct_add_report_attribute_int("fooINT", 42);
  sct_add_report_attribute_long("fooLONG", (long)42);
  sct_add_report_attribute_float("fooFLOAT", (float)15.7);
  sct_add_report_attribute_double("fooDOUBLE", 15.7);
  sct_add_report_attribute_string("bar", "test");
#endif

  sct_report(SCT_GETENV, SCT_GETENV, SCT_GETENV);
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
#endif

  /* check for tests that might have failed */
  int err = 0;

  /* I/ get timer info of each thread parallel part for timer 0 */
  double* timer_val;
  int cnt;
#if _OPENMP
  cnt = nbthreads+1;
#else
  cnt = nbthreads;
#endif
  timer_val = (double *) calloc(cnt,sizeof(double));

#pragma omp parallel
  {
    timer_val[tid] = sct_val(timer[0]);
  }

#if _OPENMP
  /* II/ get timer info of seq part for timer 1 */
  timer_val[cnt-1] = sct_val(timer[1]);
#endif

  /* gather timer info of each task */
#ifdef HAVE_MPI
  double* timer_val_glob;

  if (pid == 0) {
    timer_val_glob = (double *) calloc(cnt*nprocs,sizeof(double));
  }

  MPI_Gather(timer_val, cnt, MPI_DOUBLE, timer_val_glob, cnt, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* only the master performs the check */
  if (pid == 0) {
#endif

    /* open hdf5 outfile */
    hid_t file_id, dataset_id;
    herr_t status;
    char *fn = getenv("SCT_FILENAME");
    if (! fn) {
      fn = "sct-timings.h5";
    }
    if ( (file_id = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
      err = 1;
    }
    
    /* check timer_val against data from file */

    double timer_tsum[nprocs][cnt][timer_max];
    memset(timer_tsum, 0, nprocs*cnt*timer_max*sizeof(double));

    dataset_id = H5Dopen(file_id, "/test_hdf5/time/timer_tsum", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		     timer_tsum);
    status = H5Dclose(dataset_id);

    for (j=0; j<nprocs; j++) {
      for (i=0; i<cnt; i++) {
#ifdef HAVE_MPI
#ifdef _OPENMP
	if (i == cnt-1) {
          /* check II/ for seq part of timer[1] */
          if (timer_val_glob[j*cnt+i] != timer_tsum[j][i][1]) {err = 1;}
  	  printf("val[%d][%d] = %f , hdf5: %f\n",j,i,timer_val_glob[j*cnt+i],timer_tsum[j][i][1]);
        } else {
#endif
          /* check I/ for thread parallel part of timer[0] */
          if (timer_val_glob[j*cnt+i] != timer_tsum[j][i][0]) {err = 1;}
  	  printf("val[%d][%d] = %f , hdf5: %f\n",j,i,timer_val_glob[j*cnt+i],timer_tsum[j][i][0]);
#ifdef _OPENMP
        }
#endif
#else
#ifdef _OPENMP
	if (i == cnt-1) {
          /* check II/ for seq part of timer[1] */
          if (timer_val[i] != timer_tsum[j][i][1]) {err = 1;}
	  printf("val[%d][%d] = %f , hdf5: %f\n",j,i,timer_val[i],timer_tsum[j][i][1]);
        } else {
#endif
          /* check I/ for thread parallel part of timer[0] */
          if (timer_val[i] != timer_tsum[j][i][0]) {err = 1;}
	  printf("val[%d][%d] = %f , hdf5: %f\n",j,i,timer_val[i],timer_tsum[j][i][0]);
#ifdef _OPENMP
        }
#endif
#endif
      }
    }
    
    status = H5Fclose(file_id);

    /* free memory */
#ifdef HAVE_MPI
    free(timer_val_glob);
  }
#endif
  free(timer_val);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return err;
}

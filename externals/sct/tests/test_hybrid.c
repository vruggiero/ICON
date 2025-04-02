/**
 * \file test_hybrid.c
 * A simple code showing how to use SCT in hybrid MPI-OpenMP mode
 *
 * \copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
 *                                 Hendryk Bockelmann <bockelmann@dkrz.de>
 *
 * \author Jörg Behrens <behrens@dkrz.de>
 * \author Hendryk Bockelmann <bockelmann@dkrz.de>
 */

#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#if _OPENMP
#  include <omp.h>
#endif

#include "sct.h"

static struct timeval timeval0;
static const int timer_max = 4;

static const double pi = 3.14159265358979323846264338327950288;
static int nprocs = 1;
static int pid = 0;
static int debug = 1;
static int nthreads = 1;


void test_abort(const char *reason, const char *fname, const int line) {

  if (reason||fname) fprintf(stderr, "ERROR in file %s, line %d: %s\n",fname, line, reason);

  int errorcode;
  errorcode = 1;
  MPI_Abort(MPI_COMM_WORLD, errorcode);

  abort();
}


static inline int get_tid() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}


static double simple_walltime() {
  struct timeval tbuf;

  gettimeofday(&tbuf, NULL);

  double time_in_secs = (double) ( (tbuf.tv_sec  - timeval0.tv_sec) +
                                   (tbuf.tv_usec - timeval0.tv_usec)*1.0e-6
                                   );
  return time_in_secs;
}


static void init0(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  if ( MPI_Comm_size(MPI_COMM_WORLD, &nprocs) != MPI_SUCCESS) test_abort("MPI_Comm_size failed", __FILE__, __LINE__);
  if ( MPI_Comm_rank(MPI_COMM_WORLD, &pid) != MPI_SUCCESS) test_abort("MPI_Comm_rank failed", __FILE__, __LINE__);

#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 0;
#endif
  if (debug) fprintf(stderr,"# (test_hybrid) nprocs=%i, pid=%i, nthreads=%i\n",nprocs, pid, nthreads);
  gettimeofday(&timeval0, NULL);
  sct_init(timer_max, "test_hybrid", MPI_COMM_WORLD);
}


double work(int m) {
  int imax = 1000 * 100 * m;
  double s = 0.0;
  for (int i = 0; i < imax; i++) {
    double x = 1.0*i;
    s += sin(0.5*pi*exp(-sqrt(x)));
  }
  if (s<0.0) fprintf(stderr,"# (test_hybrid) s=%e\n",s); // fake interest in result
  return s;
}


void ref_mpi_dt(char *label, int pid, double dt) {
  // thread serial print
  fprintf(stdout, "ref_mpi_dt: label=%s, pid=%i, dt=%12.6f\n",label, pid,  dt);
}


void ref_mpi_omp_dt(char *label, int pid, int tid, double dt) {
  // thread parallel print
  fprintf(stdout, "ref_mpi_omp_dt: label=%s, pid=%i, tid=%i, dt=%12.6f\n",label, pid, tid,  dt);
}


int main (int argc, char *argv[] ) {
  init0(argc, argv);

  int ltimer=1;
  int lstimer=1;
  int lptimer=1;
  int lsptimer=1;

  char label[] = "total_timer";
  int timer = 0;
  double t0, t1;
  if (ltimer) timer = sct_new_timer(label);

  char s_label[] = "thread_serial";
  int s_timer = 0;
  double s_t0, s_t1;
  if (lstimer) s_timer = sct_new_timer(s_label);

  char p_label[] = "thread_parallel";
  int p_timer = 0;
  double p_t0, p_t1;
  if (lptimer) p_timer = sct_new_timer(p_label);

  char sp_label[] = "mixed";
  int sp_timer = 0;
  double sp_t0, sp_t1;
  if (lsptimer) sp_timer = sct_new_timer(sp_label);

  if (ltimer) {
    t0 = simple_walltime();
    sct_start(timer);
  }

  for (int irun = 0; irun < 2; irun++) {

    if (lstimer) {
      s_t0 = simple_walltime();
      sct_start(s_timer);
    }

    if (lsptimer) {
      sp_t0 = simple_walltime();
      sct_start(sp_timer);
    }

    work(10*(irun+1)*(pid+1));

    if (lsptimer) {
      sct_stop(sp_timer);
      sp_t1 = simple_walltime();
      ref_mpi_dt(sp_label, pid, sp_t1-sp_t0);
    }

#pragma omp parallel default(shared)
    {
      double sp_omp_t0, sp_omp_t1;
      double p_omp_t0, p_omp_t1;

      if (lptimer) {
        p_omp_t0 = simple_walltime();
        sct_start(p_timer);
      }

      if (lsptimer) {
        sp_omp_t0 = simple_walltime();
        sct_start(sp_timer);
      }

      int tid = get_tid();
      work((irun+1)*(tid+1));

      if (lsptimer) {
        sct_stop(sp_timer);
        sp_omp_t1 = simple_walltime();
        ref_mpi_omp_dt(sp_label, pid, tid, sp_omp_t1-sp_omp_t0);
      }

      if (lptimer) {
        sct_stop(p_timer);
        p_omp_t1 = simple_walltime();
        ref_mpi_omp_dt(p_label, pid, tid, p_omp_t1-p_omp_t0);
      }
    }

    if (lstimer) {
      sct_stop(s_timer);
      s_t1 = simple_walltime();
      ref_mpi_dt(s_label, pid, s_t1-s_t0);
    }

  }

  if (ltimer) {
    sct_stop(timer);
    t1 = simple_walltime();
    ref_mpi_dt(label, pid, t1-t0);
  }

  fflush(stdout);

  sct_report(SCT_GETENV, SCT_GETENV, SCT_GETENV);

  MPI_Finalize();

  return 0;
}

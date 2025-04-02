/**
 * \file test_omp.c
 * A simple code showing how to use SCT in OpenMP mode.
 *
 * \copyright Copyright  (C)  2013 Joerg Behrens <behrens@dkrz.de>
 *                                 Hendryk Bockelmann <bockelmann@dkrz.de>
 *
 * \author Joerg Behrens <behrens@dkrz.de>
 * \author Hendryk Bockelmann <bockelmann@dkrz.de>
 */

#include <math.h>
#include <sys/time.h>
#include <stdlib.h>

#if _OPENMP
#  include <omp.h>
#endif

#include "sct.h"

#define MY_TIMER_MAX 3
static const int timer_max = MY_TIMER_MAX;
static int timer[MY_TIMER_MAX];

static int thread_id;
#pragma omp threadprivate (thread_id)

static struct timeval timeval0;

static double simple_walltime() {
  struct timeval tbuf;

  gettimeofday(&tbuf, NULL);

  double time_in_secs = (double) ( (tbuf.tv_sec  - timeval0.tv_sec) +
                                   (tbuf.tv_usec - timeval0.tv_usec)*1.0e-6);
  return time_in_secs;
}

void error(const char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

static void init0() {
  gettimeofday(&timeval0, NULL);

  sct_init(timer_max, "test_omp", 0);
  for (int i = 0; i < timer_max; i++) timer[i] = 0;

#pragma omp parallel
  {
    thread_id = omp_get_thread_num();
  }
}

void ref_serial_dt(char *label, double dt) {
  // thread serial print
  fprintf(stdout, "ref_serial_dt: label=%s, dt=%12.6f\n",label,  dt);
}

void ref_omp_dt(char *label, int tid, double dt) {
  // thread serial print
  fprintf(stdout, "ref_omp_dt: label=%s, tid=%i, dt=%12.6f\n",label, tid,  dt);
}

void check_tid(int tid) {
  if (tid != thread_id) {
    error("unexpected change of thread pool identifier\n");
  }
}


double work(int i1, int i2) {
  double s = 0.0;
  for (int i=i1; i<=i2; i++) {
    double x = 1.0*i;
    s += floor(10*sin(sqrt(x)));
  }
  return s;
}


int main ( int argc, char *argv[] ) {

  init0();

  // setup sct timer
  char label0[] = "pure-seq-timer-0";
  timer[0] = sct_new_timer(label0);

  char label1[] = "mixed-called-timer-1";
  timer[1] = sct_new_timer(label1);

  char label2[] = "pure-par-timer-2";
  timer[2] = sct_new_timer(label2);

  double ts[2], te[2];
  int icount = 10*1000*1000;

  // setup omp
  int max_threads= omp_get_max_threads();
  printf("# test_omp: max_threads=%i\n",max_threads);
  if (max_threads<1) {
    fprintf(stderr,"unexpected zero count of threads\n");
    return 1;
  }

  // start total timer
  sct_start(timer[0]);
  ts[0] = simple_walltime();

  // perform pure serial workload
  work(0, icount);

  // reference timer resolution:
  double ref_res;
  {
    double t0 = simple_walltime();
    double t1;
    while ( (t1=simple_walltime()) <= t0 ) {};
    ref_res = t1-t0;
  }
  printf("# ref_res = %.12e\n",ref_res);

  // sct timer resolution (not using internal method):
#pragma omp parallel default(shared)
  {
    int tid = omp_get_thread_num();
    check_tid(tid);
    double ts, te, t0, t1;

    sct_start(timer[2]);
    ts = simple_walltime();

    t0 = sct_val(timer[2]);
    while ( (t1=sct_val(timer[2])) <= t0 ) {};
    double res = t1-t0;

    te = simple_walltime();
    sct_stop(timer[2]);

    printf("# tid=%i, sct_res = %.12e\n",tid, res);
  }

  // serial measurement of mixed timer
  sct_start(timer[1]);
  ts[1] = simple_walltime();

  //define work load
  int iwork[max_threads];
  iwork[0] = icount/2;
  int wsum = iwork[0];
  for(int it = 1; it < max_threads; it++) {
    iwork[it] = iwork[it-1]/2;
    wsum += iwork[it];
  }
  fprintf(stderr,"# test_omp: icount=%i, wsum=%i\n", icount, wsum);
  iwork[0] += (icount-wsum);

  int istart[max_threads];
  int iend[max_threads];
  istart[0] = 0;
  iend[0] = iwork[0] -1;
  for(int it = 1; it < max_threads; it++) {
    istart[it] = iend[it-1] + 1;
    iend[it] = istart[it] + iwork[it] - 1;
  }
  int it = max_threads-1;
  if (iend[max_threads-1] != icount-1) {
    error("# test_omp: internal error\n");
  }

  // thread parallel measurement:
  const int ref_test_sum = -4984107;
  int test_sum = 0;
  double thread_dt[max_threads];
#pragma omp parallel default(shared) reduction(+:test_sum)
  {
    int tid = omp_get_thread_num();
    check_tid(tid);
    double ts[2], te[2], t0, t1;

    sct_start(timer[1]);
    ts[0] = simple_walltime();

    sct_start(timer[2]);
    ts[1] = simple_walltime();

    t0 = simple_walltime();

    test_sum += work(istart[tid], iend[tid]);

    t1 = simple_walltime();
    double res = t1-t0;

    sct_stop(timer[2]);
    te[1] = simple_walltime();
    ref_omp_dt(label2, tid, te[1]-ts[1]);

    sct_stop(timer[1]);
    te[0] = simple_walltime();
    ref_omp_dt(label1, tid, te[0]-ts[0]);

    // critical not strictly required here:
#pragma omp critical
    {
      thread_dt[tid] = res;
    }
  }

  sct_stop(timer[1]);
  te[1] = simple_walltime();
  ref_serial_dt(label1, te[1]-ts[1]);

  if (test_sum != ref_test_sum) {
    printf("# unexpected test-sum=%i\n",test_sum);
    error("internal error\n");
  }

  sct_stop(timer[0]);
  te[0] = simple_walltime();
  ref_serial_dt(label0, te[0]-ts[0]);

  sct_report(SCT_GETENV, SCT_GETENV, SCT_GETENV);

  return 0;
}

/**
 * \file test_serial.c
 * A simple code showing how to use SCT in serial mode.
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

#include "sct.h"

#define MY_TIMER_MAX 4
static const int timer_max = MY_TIMER_MAX;
static int timer[MY_TIMER_MAX];
static int ltimer[MY_TIMER_MAX];

static struct timeval timeval0;

static const int debug = 0;

void error(const char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

static double simple_walltime() {
  struct timeval tbuf;

  gettimeofday(&tbuf, NULL);

  double time_in_secs = (double) ( (tbuf.tv_sec  - timeval0.tv_sec) +
                                   (tbuf.tv_usec - timeval0.tv_usec)*1.0e-6
                                   );
  return time_in_secs;
}

static void init0() {
  gettimeofday(&timeval0, NULL);

  // enforce a smaller number of timer_max to check realloc of resources
  sct_init(timer_max/2, "test_serial", 0);

  for (int i = 0; i < timer_max; i++) {
    timer[i] = 0;
    ltimer[i] = 1;
  }
}

void ref_serial_dt(char *label, double dt) {
  // thread serial print
  fprintf(stdout, "ref_serial_dt: label=%s, dt=%12.6f\n",label,  dt);
}


int main ( int argc, char *argv[] ) {

  init0();

  double ts[4], te[4];

  // sct timer:
  char label0[] = "long-label-timer-0";
  if (ltimer[0]) timer[0] = sct_new_timer(label0);

  char label1[] = "even-longer-label-timer-1";
  if (ltimer[1]) timer[1] = sct_new_timer(label1);

  char label2[] = "timer-2";
  if (ltimer[2]) timer[2] = sct_new_timer(label2);

  char label3[] = "timer-3";
  if (ltimer[3]) timer[3] = sct_new_timer(label3);

  int icount0 = 5*1000*1000;

  // total
  if (ltimer[0]) {
    sct_start(timer[0]);
    ts[0] = simple_walltime();
  }

  // reference timer resolution:
  double ref_res;
  {
    double t0 = simple_walltime();
    double t1;
    while ( (t1=simple_walltime()) == t0 ) {};
    ref_res = t1-t0;
  }
  printf("# ref_res = %.12e\n",ref_res);

  for (int irun = 0; irun<4; irun++) {
    int icount = icount0 * (1+irun);
    //define work load
    int iwork = icount;
    //printf("# iwork = %i\n",iwork);
    int wsum = iwork;
    int istart;
    int iend;
    istart = 0;
    iend = iwork -1;

    if (ltimer[1]) {
      sct_start(timer[1]);
      ts[1] = simple_walltime();
    }
    if (ltimer[2] && irun/2) {
      sct_start(timer[2]);
      ts[2] = simple_walltime();
    }
    if (ltimer[3]) {
      sct_start(timer[3]);
      ts[3] = simple_walltime();
    }

    int sum = 0;
    for (int i=istart; i<=iend; i++) {
      double x = 1.0*i;
      sum += floor(10*sin(sqrt(x)));
    }
    if (sum == 0.12345) printf("nop");

    if (ltimer[3]) {
      sct_stop(timer[3]);
      te[3] = simple_walltime();
      ref_serial_dt(label3, te[3]-ts[3]);
    }
    if (ltimer[2] && irun/2) {
      sct_stop(timer[2]);
      te[2] = simple_walltime();
      ref_serial_dt(label2, te[2]-ts[2]);
    }
    if (ltimer[1]) {
      sct_stop(timer[1]);
      te[1] = simple_walltime();
      ref_serial_dt(label1, te[1]-ts[1]);
    }

    // test for intermediate sct_report
    //if (irun == 2) sct_single_report(timer[1], SCT_GETENV, SCT_GETENV, SCT_GETENV);
  }

  if (ltimer[0]) {
    sct_stop(timer[0]);
    double last_dt = sct_last_dt(timer[0]);
    if (last_dt != sct_val(timer[0])) error("last_dt != sct_val(it)");
    te[0] = simple_walltime();
    ref_serial_dt(label0, te[0]-ts[0]);
  }

  sct_report(SCT_GETENV, SCT_GETENV, SCT_GETENV);

  sct_finalize();

  return 0;
}

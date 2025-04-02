#include <math.h>
#include <sys/time.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "sct.h"

#define MY_TIMER_MAX 3
static const int timer_max = MY_TIMER_MAX;
static int timer[MY_TIMER_MAX];

static struct timeval timeval0;

static int err = 0;

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
  sct_init(timer_max, "test_precision", SCT_COMM_WORLD);
  for (int i = 0; i < timer_max; i++) timer[i] = 0;
}

int main ( int argc, char *argv[] ) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  init0();

  // sct timer:
  timer[0] = sct_new_timer("timer 0");
  timer[1] = sct_new_timer("timer 1");
  timer[2] = sct_new_timer("timer 2");

  double d = 0.0;
  int icount = 10*1000*1000;

  // reference timer:
  double ref_start[timer_max];
  double ref_end[timer_max];

  // reference timer resolution:
  double ref_res;
  {
    double t0 = simple_walltime();
    double t1;
    while ( (t1=simple_walltime()) <= t0 ) {};
    ref_res = t1-t0;
  }
  printf("ref_res = %.12e\n",ref_res);


  // sct timer resolution (not using internal method):
  double sct_res;
  {
    double t0,t1;
    sct_start(timer[0]);
    t0 = sct_val(timer[0]);
    while ( (t1 = sct_val(timer[0])) <= t0 ) {};
    sct_res = t1-t0;
    sct_stop(timer[0]);
  }
  printf("sct_res = %.12e\n",sct_res);

  // sct timer internal resolution:
  double sct_res_internal = sct_resolution();
  printf("sct_res_internal = %.12e\n",sct_res_internal);

  sct_start(timer[0]);
  sct_stop(timer[0]);

  ref_start[0] = simple_walltime();
  ref_end[0] = simple_walltime();

  // nested calls:
  sct_start(timer[1]);
  ref_start[1] = simple_walltime();

  sct_start(timer[2]);
  ref_start[2] = simple_walltime();

  for (int i =0; i<icount; i++) {
    double x = 1.0*i;
    d = d + exp(-sqrt(x));
  }

  ref_end[2] = simple_walltime();
  sct_stop(timer[2]);

  ref_end[1] = simple_walltime();
  sct_stop(timer[1]);

  if (d>1.e+3) printf("result=%.12e\n",d); //fake interest in the result to the compiler

  double ref_dt[3];
  double sct_dt[3];
  for (int i = 0; i < timer_max; i++) {
    sct_dt[i] = sct_val(timer[i]);
    ref_dt[i] = ref_end[i] - ref_start[i];
  }

  // reference call over head (timer 0)
  if (ref_dt[0] <= ref_res)
    printf("ref call overhead <= ref_res\n");
  else
    printf("ref call overhead  = %.12e\n",ref_dt[0]);

  // sct call over head (timer 0)
  if (sct_dt[0] <= sct_res)
    printf("sct call overhead  <= sct_res\n");
  else
    printf("sct call overhead  = %.12e\n",sct_dt[0]);

  // check nested characteristics
  double ref_dt1dt2 = ref_dt[1] - ref_dt[2];
  double sct_dt1dt2 = sct_dt[1] - sct_dt[2];

  if (ref_dt1dt2<0.0) {
    printf("unexpected nesting behaviour: ref_dt1dt2=%.12e\n",ref_dt1dt2);
    err++;
  };
  if (sct_dt1dt2<0.0) { printf("unexpected nesting behaviour: sct_dt1dt2=%.12e\n",sct_dt1dt2); err++; }

  double tol =  sct_res + ref_res;
  double d1 = sct_dt[1]-ref_dt[1];
  double d2 = sct_dt[2]-ref_dt[2];
  if (d2<-tol) { printf("unexpected nesting behaviour: sct_dt[2]=%.12e, ref_dt[2]=%.12e, diff=%e\n",sct_dt[2], ref_dt[2], d2); err++; };
  if (d1<-tol) { printf("unexpected nesting behaviour: sct_dt[1]=%.12e, ref_dt[1]=%.12e, diff=%e\n",sct_dt[1], ref_dt[1], d1); err++; };

  sct_report(SCT_GETENV, SCT_GETENV, SCT_GETENV);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  printf("err=%i\n",err);
  return err;
}

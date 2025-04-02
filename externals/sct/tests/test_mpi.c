/**
 * \file test_mpi.c
 * A simple code showing how to use SCT in MPI mode
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
#include <string.h>
#include <mpi.h>

#include "sct.h"

#define DIFFERENT_NUMBER_TIMERS 1
#define DIFFERENT_NAMES_TIMERS 0
#define DIFFERENT_ORDER_TIMERS 0
#define MIN_NUMBER_TIMERS 4
static int *timer;

static struct timeval timeval0;

static int err = 0;
static int nprocs = 0;
static int pid = 0;

static const int debug = 0;

void test_abort(const char *reason, const char *fname, const int line) {

  if (reason||fname) fprintf(stderr, "ERROR in file %s, line %d: %s\n",fname, line, reason);

  int errorcode;
  errorcode = 1;
  MPI_Abort(MPI_COMM_WORLD, errorcode);

  abort();
}

static double simple_walltime() {
  struct timeval tbuf;

  gettimeofday(&tbuf, NULL);

  double time_in_secs = (double) ( (tbuf.tv_sec  - timeval0.tv_sec) +
                                   (tbuf.tv_usec - timeval0.tv_usec)*1.0e-6
                                   );
  return time_in_secs;
}

static void init0(int timer_max) {
  gettimeofday(&timeval0, NULL);

  sct_init(timer_max, "test_mpi", SCT_COMM_WORLD);

  timer = calloc(timer_max, sizeof(int));
  for (int i = 0; i < timer_max; i++) timer[i] = 0;
}

void ref_mpi_dt(char *label, int pid, double dt) {
  // thread serial print
  fprintf(stdout, "ref_mpi_dt: label=%s, pid=%i, dt=%12.6f\n",label, pid,  dt);
}


int main ( int argc, char *argv[] ) {

  MPI_Init(&argc, &argv);

  if ( MPI_Comm_size(MPI_COMM_WORLD, &nprocs) != MPI_SUCCESS) test_abort("MPI_Comm_size failed", __FILE__, __LINE__);
  if ( MPI_Comm_rank(MPI_COMM_WORLD, &pid) != MPI_SUCCESS) test_abort("MPI_Comm_rank failed", __FILE__, __LINE__);

  if (DIFFERENT_NUMBER_TIMERS && (pid%2 == 1))
    init0(MIN_NUMBER_TIMERS+1);
  else
    init0(MIN_NUMBER_TIMERS);

  double ts[MIN_NUMBER_TIMERS+1], te[MIN_NUMBER_TIMERS+1];

  // sct timer:
  char label0[] = "long-label-timer-0";
  timer[0] = sct_new_timer(label0);

  char label1[] = "even-longer-label-timer-1";
  char label2[] = "timer-2";
  if (DIFFERENT_ORDER_TIMERS && (pid%2 == 0)) {
    timer[2] = sct_new_timer(label2);
    timer[1] = sct_new_timer(label1);
  }
  else {
    timer[1] = sct_new_timer(label1);
    timer[2] = sct_new_timer(label2);
  }

  char label3[11];
  if (DIFFERENT_NAMES_TIMERS && (pid%2==0))
    strncpy(label3, "timer-even\0", sizeof(label3));
  else
    strncpy(label3, "timer-oddd\0", sizeof(label3));
  timer[3] = sct_new_timer(label3);

  // fake inbalance of used timers between MPI-tasks
  char label4[] = "timer-4";
  if (DIFFERENT_NUMBER_TIMERS && (pid%2 == 1)) {
    fprintf(stdout,"pid = %d sets extra timer '%s' ...\n", pid, label4);
    timer[4] = sct_new_timer(label4);
  }

  int icount0 = 10*1000*1000;

  // total
  sct_start(timer[0]);
  ts[0] = simple_walltime();

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
  double sct_res;
  //printf("test_mpi: nprocs=%i\n",nprocs);
  if (nprocs<1) {
    fprintf(stderr,"unexpected zero count of procs\n");
    return 1;
  }
  for (int irun = 0; irun<=3; irun++) {
    int icount = icount0 * (1+irun);
    //define work load
    int iwork[nprocs];
    iwork[0] = icount/2;
    int wsum = iwork[0];
    for(int it = 1; it < nprocs; it++) {
      iwork[it] = iwork[it-1]/2;
      wsum += iwork[it];
    }
    //fprintf(stderr,"test_mpi[%d]: icount=%i, wsum=%i\n", pid, icount, wsum);
    iwork[0] += (icount-wsum);

    int istart[nprocs];
    int iend[nprocs];
    istart[0] = 0;
    iend[0] = iwork[0] -1;
    for(int it = 1; it < nprocs; it++) {
      istart[it] = iend[it-1] + 1;
      iend[it] = istart[it] + iwork[it] - 1;
    }
    int it = nprocs-1;
    if (iend[nprocs-1] != icount-1) {
      fprintf(stderr,"# test_mpi: internal error\n");
      return 1;
    }

    // proc parallel measurement:

    int part_sum = 0;
    double proc_dt[nprocs];
    sct_start(timer[1]);
    ts[1] = simple_walltime();

    if (irun/2) {
      sct_start(timer[2]);
      ts[2] = simple_walltime();
    }

    sct_start(timer[3]);
    ts[3] = simple_walltime();

    double t0 = simple_walltime();

    for (int i=istart[pid]; i<=iend[pid]; i++) {
      double x = 1.0*i;
      part_sum += floor(10*sin(sqrt(x)));
      //fprintf(stderr,"pid=%i, it=%i, part_sum=%i\n",pid,i,part_sum);
    }
    //fprintf(stderr,"sub-final:pid=%i, part_sum=%i\n",pid,part_sum);
    double t1 = simple_walltime();
    double dt = t1 - t0;

    // fake inbalance of used timers between MPI-tasks
    if (DIFFERENT_NUMBER_TIMERS && (pid%2 == 1)) {
      sct_start(timer[4]);
      ts[4] = simple_walltime();

      double foo = 0.;
      for (int i=istart[pid]; i<=iend[pid]; i++) {
        double x = 1.0*i;
        foo += 10*sin(sqrt(x))*sin(sqrt(x));
      }

      sct_stop(timer[4]);
      te[4] = simple_walltime();
      ref_mpi_dt(label4, pid, te[4]-ts[4]);
      if (foo <= 0.) fprintf(stderr,"PANIC: foo = %f\n",foo);
    }

    sct_stop(timer[3]);
    te[3] = simple_walltime();
    ref_mpi_dt(label3, pid, te[3]-ts[3]);

    if (irun/2) {
      sct_stop(timer[2]);
      te[2] = simple_walltime();
      ref_mpi_dt(label2, pid, te[2]-ts[2]);
    }

    sct_stop(timer[1]);
    te[1] = simple_walltime();
    ref_mpi_dt(label1, pid, te[1]-ts[1]);


    // gather partial results and check sum
    {
      int full_sum;
      int recvbuf[nprocs];
      for (int i = 0; i < nprocs; i++) { recvbuf[i] = 1234567; };
      //fprintf(stderr,"test_mpi[%d] part_sum=%d\n",pid,part_sum);
      if (
          MPI_Allgather(&part_sum, // sendbuf
                        1, // sendcount,
                        MPI_INT, //MPI_Datatype sendtype,
                        recvbuf,//void* recvbuf,
                        1,//int recvcount,
                        MPI_INT,//MPI_Datatype recvtype,
                        MPI_COMM_WORLD) //MPI_Comm comm)
          != MPI_SUCCESS) test_abort("MPI_Allgather failed", __FILE__, __LINE__);
      full_sum = 0;
      for (int i = 0; i < nprocs; i++) {
        full_sum += recvbuf[i];
        //fprintf(stderr,"test_mpi[%d] recvbuf[%d]=%d\n",pid,i,recvbuf[i]);
      }
      //fprintf(stderr,"test_mpi[%d]: full_sum=%d\n",pid,full_sum);
      if (icount == 10*1000*1000) {
        const int ref_test_sum = -4984107;
        if (full_sum == ref_test_sum) {
          fprintf(stderr,"# test_mpi: test_sum okay\n");
        } else {
          fprintf(stderr,"# test_mpi: internal error: unexpected sum=%i\n",full_sum);
          err++;
          return err;
        }
      }
    }

  // gather time measurement:
    double all_dt[nprocs];
    {
      for (int i = 0; i < nprocs; i++) { all_dt[i] = -1.0; };
      //fprintf(stderr,"test_mpi[%d] dt=%e\n",pid,dt);
      if (
          MPI_Allgather(&dt, // sendbuf
                        1, // sendcount,
                        MPI_DOUBLE, //MPI_Datatype sendtype,
                        all_dt,//void* recvbuf,
                        1,//int recvcount,
                        MPI_DOUBLE,//MPI_Datatype recvtype,
                        MPI_COMM_WORLD) //MPI_Comm comm)
          != MPI_SUCCESS) test_abort("MPI_Allgather failed", __FILE__, __LINE__);
      /*
      for (int i = 0; i < nprocs; i++) {
        fprintf(stderr,"test_mpi[%d] all_dt[%d]=%e\n",pid,i,all_dt[i]);
      }
      */
    }

    {
      double tmin, tsum, tavg, tmax, eff;
      tmin = all_dt[0];
      tmax = all_dt[0];
      tsum = all_dt[0];
      for (int pid = 1; pid < nprocs; pid++) {
        double t = all_dt[pid];
        if (t<tmin) tmin = t;
        if (t>tmax) tmax = t;
        tsum += t;
      }
      tavg = tsum/nprocs;
      eff = tavg/tmax;
      if (debug) printf("# test_mpi[%d]: reference stats for timer3: tmin=%.12e, tavg=%.12e, tmax=%.12e, eff=%.12e\n",
                        pid, tmin, tavg, tmax, eff);
    }
  }

  sct_stop(timer[0]);
  te[0] = simple_walltime();
  ref_mpi_dt(label0, pid, te[0]-ts[0]);

  fflush(stdout);

  sct_report(SCT_GETENV, SCT_GETENV, SCT_GETENV);
  //sct_single_report(3, SCT_GETENV, SCT_GETENV, SCT_GETENV);

  MPI_Finalize();

  //printf("err=%i\n",err);
  return err;
}

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
 * \file sct_reduce.c
 * all functions used to reduce various measurement datasets from tasks/threads
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sct_collector.h"
#include "sct_reduce.h"
#include "sct_reporter.h"
#include "sct_mergesort.h"
#ifndef SCT_GETENV
#  error undef problem
#endif
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static const char doc_sp_merge_simple[] = "SCT_SP_MERGE_SIMPLE";
static const char doc_sp_serial_only[] = "SCT_SP_SERIAL_ONLY";
static const char doc_sp_parallel_only[] = "SCT_SP_PARALLEL_ONLY";
static const char doc_sp_select_all[] = "SCT_SP_SELECT_ALL";

static const char doc_select_all[] = "SCT_SELECT_ALL";
static const char doc_reduce_all[] = "SCT_REDUCE_ALL";
static const char doc_select_one[] = "";

// used timers over all MPI ranks:
static int global_timer_num = 0;
static char (*timer_list)[SCT_LABEL_SIZE] = NULL;
static int *timer_map = NULL; // timer_map(local_ind) = global_ind in timer_list


const char *sct_internal_reduce_doc(int choice) {

  if (choice == SCT_SELECT_ALL) {
    return doc_select_all;
  } else if (choice == SCT_REDUCE_ALL) {
    return doc_reduce_all;
  } else if (choice >= 0) {
    return doc_select_one;
  } else {
    sct_abort("sct_reduce_doc: unexpected value for choice var", __FILE__, __LINE__);
  }
  return NULL;
}


const char *sct_internal_sp_doc(int sp_merging) {

  switch(sp_merging) {
  case SCT_SP_MERGE_SIMPLE:
    return doc_sp_merge_simple;
    break;
  case SCT_SP_SERIAL_ONLY:
    return doc_sp_serial_only;
    break;
  case SCT_SP_PARALLEL_ONLY:
    return doc_sp_parallel_only;
    break;
  case SCT_SP_SELECT_ALL:
    return doc_sp_select_all;
    break;
  default:
    sct_abort("sct_sp_doc: unexpected value for sp_merging", __FILE__, __LINE__);
  }
  return NULL;
}


/*
  get timer of timer over all ranks (duplicates removed)
 */
int sct_get_global_timer_num() {
#ifdef HAVE_MPI
  return global_timer_num;
#else
  return sct_get_timer_num();
#endif
}

/*
  map process local timer index it to index in global timer_list
 */
int sct_get_global_idx(int it) {
  if (it<0 || it>=sct_get_timer_num())
    sct_abort("sct_get_global_idx: timer out of bounds", __FILE__, __LINE__);
#ifdef HAVE_MPI
  return timer_map[it];
#else
  return it;
#endif
}


/*
  get name of process pid local timer it
 */
char *sct_get_global_timer_cname(const sct_reduction_type *r, int pid, int it) {
#ifdef HAVE_MPI
  if (r->proc_choice == SCT_REDUCE_ALL) {
    if (!r || it<0 || it>=r->global_timer_num)
      sct_abort("sct_get_global_timer_cname: timer out of bounds", __FILE__, __LINE__);

    return timer_list[it];
  }
  else {
    if (!r || it<0 || it>=r->timer_num_per_rank[pid])
      sct_abort("sct_get_global_timer_cname: timer out of bounds", __FILE__, __LINE__);

    // determine offset in r->timer_map_per_rank
    int offset = 0;
    for (int ip=0; ip<pid; ip++) offset += r->timer_num_per_rank[ip];

    return timer_list[r->timer_map_per_rank[offset+it]];
  }
#else
  if (!r || pid != 0 || it<0 || it>=r->timer_num_per_rank[pid])
    sct_abort("sct_get_global_timer_cname: timer out of bounds", __FILE__, __LINE__);

  return sct_get_timer_cname(it);
#endif
}


#ifdef HAVE_MPI
void sct_create_global_timer_map(sct_context_type *con) {
  const int my_debug = 0;

  typedef char sct_string[SCT_LABEL_SIZE];
  MPI_Datatype mpi_string;

  int *displs, *rcounts;
  int timer_num = sct_get_timer_num();
  const int pid = con->pid;

  sct_string *timer_table;

  if (my_debug) {
    fprintf(stderr,"[%d] local timer_num %d\n", pid, timer_num);
    fflush(stderr);
  }

  // determine global_timer_num
  displs = calloc(con->procnum, sizeof(int));
  rcounts = calloc(con->procnum, sizeof(int));
  if ( (MPI_Allgather(&timer_num, 1, MPI_INT,
		      rcounts, 1, MPI_INT, con->comm) != MPI_SUCCESS) )
    sct_abort("MPI_Allgather failed", __FILE__, __LINE__);

  displs[0] = 0;
  global_timer_num = rcounts[0];
  for (int ip=1; ip<con->procnum; ip++) {
    global_timer_num += rcounts[ip];
    displs[ip] = displs[ip-1] + rcounts[ip-1];
  }

  if (my_debug) {
    fprintf(stderr,"[%d] global_timer_num %d ...\n", pid, global_timer_num);
    for (int ip=0; ip<con->procnum; ip++)
      fprintf(stderr,"[%d] rank %d rcounts %d displs %d\n", pid, ip, rcounts[ip], displs[ip]);
  }

  // fill global timer_table with local part
  timer_table = calloc(global_timer_num, sizeof(sct_string));
  for (int it=0; it<timer_num; it++) {
    strcpy(timer_table[displs[pid] + it], sct_get_timer_cname(it));
  }

  // allgatherv timer from any proc
  MPI_Type_vector(1,SCT_LABEL_SIZE,SCT_LABEL_SIZE,MPI_CHAR,&mpi_string);
  MPI_Type_commit(&mpi_string);
  if ( (MPI_Allgatherv(MPI_IN_PLACE, timer_num, mpi_string,
                       &timer_table[0], rcounts, displs, mpi_string, con->comm) != MPI_SUCCESS) )
    sct_abort("MPI_Allgatherv failed", __FILE__, __LINE__);

  if (my_debug) {
    for (int it=0; it<global_timer_num; it++)
      fprintf(stderr,"[%d] pre %d -> %s\n",pid,it,timer_table[it]);
  }

  // sort global timer_table and remove duplicates
  sct_mergesort_index(timer_table, global_timer_num, NULL, 0);

  char *act_timer = timer_table[0];
  int end = global_timer_num;
  global_timer_num = 1;
  for (int it=1; it<end; it++) {
    if (strcmp(timer_table[it], act_timer) != 0) {
      global_timer_num++;
      act_timer = timer_table[it];
    }
    else {
      strcpy(timer_table[it],"");
    }
  }

  timer_list = calloc(global_timer_num, sizeof(sct_string));
  timer_map = calloc(timer_num, sizeof(int));
  int c = 0;
  int matched;
  for (int it=0; it<end; it++) {
    if (strcmp(timer_table[it], "") != 0) {
      strcpy(timer_list[c],timer_table[it]);
      matched = 0;
      for (int ind=0; ind<timer_num; ind++) {
        if (strcmp(sct_get_timer_cname(ind),timer_list[c]) == 0) {
          if (!matched) {
            timer_map[ind] = c;
            matched = 1;
	  }
	  else {
            timer_map[ind] = c;
            char warn[256];
            sprintf(warn,"two local timer with same name [%s] declared - wasting results\n",
                    sct_get_timer_cname(ind));
            sct_warn(warn, __FILE__, __LINE__);
	  }
	}
      }
      c++;
    }
  }
  if (c != global_timer_num) sct_abort("mergesort failed", __FILE__, __LINE__);

  if (my_debug) {
    for (int it=0; it<global_timer_num; it++)
      fprintf(stderr,"[%d] post %d -> %s\n",pid,it,timer_list[it]);
    for (int it=0; it<timer_num; it++)
      fprintf(stderr,"[%d] map %d -> %d (%s)\n", pid,
	      it, timer_map[it], timer_list[timer_map[it]]);
  }

  MPI_Type_free(&mpi_string);
  if(displs) free(displs);
  if(rcounts) free(rcounts);
  if(timer_table) free(timer_table);
}
#endif


/*
static int get_context_choice_from_env() {
  char *s = getenv ("SCT_CONTEXT_CHOICE");
  if (!s) goto err;
  int i = atoi(s);
  if (i<0) goto err;
  return i;

 err:
  fprintf(stderr,"SCT-Warning: use context 0 as fallback for context_choice\n");
  return 0;
}
*/


/*
  add operation for sct_stats_type: s = a+b
*/
static void add_stats(const sct_stats_type *a, const sct_stats_type *b, sct_stats_type *s, int en) {
  s->tsum = a->tsum + b->tsum;
  s->cnum = a->cnum + b->cnum;
  s->tmin = MIN(a->tmin, b->tmin);
  s->tmax = MAX(a->tmax, b->tmax);
#ifdef HAVE_LIBPAPI
  for (int ie=0; ie<en; ie++) {
    s->esum[ie] = a->esum[ie] + b->esum[ie];
    if (sct_get_callstats()) {
      s->emin[ie] = MIN(a->emin[ie], b->emin[ie]);
      s->emax[ie] = MAX(a->emax[ie], b->emax[ie]);
      s->rmin[ie] = MIN(a->rmin[ie], b->rmin[ie]);
      s->rmax[ie] = MAX(a->rmax[ie], b->rmax[ie]);
    }
  }
#endif
}


/*
  accumulate operation for sct_stats_type: acc = acc + x
*/
static void accumulate_stats(sct_stats_type *x, sct_stats_type *acc, int en) {
  acc->tsum = acc->tsum + x->tsum;
  acc->cnum = acc->cnum + x->cnum;
  acc->tmin = MIN(acc->tmin, x->tmin);
  acc->tmax = MAX(acc->tmax, x->tmax);
#ifdef HAVE_LIBPAPI
  for (int ie=0; ie<en; ie++) {
    acc->esum[ie] = acc->esum[ie] + x->esum[ie];
    if (sct_get_callstats()) {
      acc->emin[ie] = MIN(acc->emin[ie], x->emin[ie]);
      acc->emax[ie] = MAX(acc->emax[ie], x->emax[ie]);
      acc->rmin[ie] = MIN(acc->rmin[ie], x->rmin[ie]);
      acc->rmax[ie] = MAX(acc->rmax[ie], x->rmax[ie]);
    }
  }
#endif
}


static void get_sp(sct_context_type *con, int proc_choice, int sp_num, sct_stats_type *sp_stats_mem) {

  // copy the timer-stats from (possibly) thread-private memory to shared memory

  sct_stats_type(*sp_stats)[sp_num] = (sct_stats_type(*)[sp_num]) sp_stats_mem;

  const int en = sct_get_event_num();
  const int timer_num = sct_get_timer_num();
  int timer_idx;

  // serial phase will be stored at the very end of sp_stats such that
  // sp_stats[0...nb_timer-1][0...nb_threads-1] for parallel executed time
  // sp_stats[0...nb_timer-1][nb_threads]       for serial executed time
  if (sp_num<1) sct_abort("get_sp: (sp_num<1)", __FILE__, __LINE__);
  for (int it = 0; it < timer_num; it++) {
    if (proc_choice == SCT_REDUCE_ALL)
      timer_idx = sct_get_global_idx(it);
    else
      timer_idx = it;

    sct_internal_copy_stats(&con->s[it], &sp_stats[timer_idx][sp_num-1], en);
  }

#ifdef _OPENMP
  // transpose thread-private timer data into shared mem
  const int p_num = sct_get_pr_thread_count();
  if (p_num+1 != sp_num) sct_abort("get_sp: (p_num > sp_num-1)", __FILE__, __LINE__);
  if (omp_in_parallel()) sct_abort("get_sp: require serial phase", __FILE__, __LINE__);
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      int tid = sct_get_pr_thread_id();
      // critical region not required, but avoids some cache line trashing
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        for (int it = 0; it < timer_num; it++) {
          if (proc_choice == SCT_REDUCE_ALL)
            timer_idx = sct_get_global_idx(it);
          else
            timer_idx = it;
          sct_internal_copy_stats(&con->p[tid][it], &sp_stats[timer_idx][tid], en);
        }
      }
    }
#endif

}


static void cleanup_stats(int th_num, sct_stats_type *x_stats_mem) {

  // after the measurement we still might have huge init values for tmin, emin and rmin

  sct_stats_type(*x_stats)[th_num] = (sct_stats_type(*)[th_num]) x_stats_mem;

  const int en = sct_get_event_num();
  const int timer_num = sct_get_timer_num();

  for (int it = 0; it < timer_num; it++) {
    for (int tid = 0; tid < th_num; tid++) {
      sct_stats_type *x = &x_stats[it][tid];
      if (x->tmin > x->tmax) {
        x->tmin = 0.0;
#ifdef HAVE_LIBPAPI
        for (int ie=0; ie<en; ie++) {
          x->emin[ie] = 0;
          x->rmin[ie] = 0;
        }
#endif
      }
    }
  }
}


void show_stats(char *label, sct_stats_type *stats) {

  // for debuging only

  fprintf(stderr, "show_stats (%s):\ntsum=%e, tmin=%e, tmax=%e, cnum=%i\n",
          label, stats->tsum, stats->tmin, stats->tmax, stats->cnum);
#ifdef HAVE_LIBPAPI
  const int en = sct_get_event_num();
  for (int ie=0; ie<en; ie++) {
    fprintf(stderr, "ie=%i (%s): esum=%lld, emin=%lld, emax=%lld\n",
            ie, sct_get_event_cname(ie), stats->esum[ie], stats->emin[ie], stats->emax[ie]);
  }
#endif

}


#ifdef _OPENMP
void merge_sp(int sp_num, int m_num, int sp_merging, sct_stats_type* sp_stats_mem,
              sct_stats_type *m_stats_mem) {

  sct_stats_type(*sp_stats)[sp_num] = (sct_stats_type(*)[sp_num]) sp_stats_mem;
  sct_stats_type(*m_stats)[m_num] = (sct_stats_type(*)[m_num]) m_stats_mem;

  const int timer_num = sct_get_timer_num();
  const int en = sct_get_event_num();

  if ( sp_merging == SCT_SP_SERIAL_ONLY)  {
    if (m_num != 1) sct_abort("merge_sp unexpected number of slots for merged threads", __FILE__, __LINE__);
    for (int it = 0; it < timer_num; it++) {
      sct_internal_copy_stats(&sp_stats[it][sp_num-1], &m_stats[it][0], en);
    }
    return;
  }

  const int p_num = sct_get_pr_thread_count();
  if ( sp_merging == SCT_SP_SELECT_ALL ) {
    if (m_num != p_num + 1)  sct_abort("merge_sp unexpected number of slots for merged threads", __FILE__, __LINE__);
  } else {
    if (m_num != p_num) sct_abort("merge_sp unexpected number of slots for merged threads", __FILE__, __LINE__);
  }

  if ( sp_merging == SCT_SP_PARALLEL_ONLY) {
    for (int it = 0; it < timer_num; it++) {
      for (int tid = 0; tid < p_num; tid++) {
        sct_internal_copy_stats(&sp_stats[it][tid], &m_stats[it][tid], en);
      }
    }
  } else if ( sp_merging == SCT_SP_MERGE_SIMPLE)  {
    if (sp_num != p_num+1) sct_abort("merge_sp: (sp_num != p_num+1)", __FILE__, __LINE__);
    for (int it = 0; it < timer_num; it++) {
      add_stats(&sp_stats[it][0], &sp_stats[it][sp_num-1], &m_stats[it][0], en);
      for (int tid = 1; tid < p_num; tid++) {
        sct_internal_copy_stats(&sp_stats[it][tid], &m_stats[it][tid], en);
      }
    }
  } else if ( sp_merging == SCT_SP_SELECT_ALL) {
    if (sp_num != p_num+1) sct_abort("merge_sp: (sp_num != p_num+1)", __FILE__, __LINE__);
    for (int it = 0; it < timer_num; it++) {
      for (int tid = 0; tid <= p_num; tid++) {
        sct_internal_copy_stats(&sp_stats[it][tid], &m_stats[it][tid], en);
      }
    }
  } else {
    sct_abort("merge_sp: unsupported sp_merging value", __FILE__, __LINE__);
  }

}
#endif


void sct_reduction_delete(sct_reduction_type *r) {
  if (!r) return;

  int timer_num = (r->proc_choice == SCT_REDUCE_ALL) ? global_timer_num : sct_get_timer_num();

  sct_string_delete(&r->context_name);
  for (int i = 0; i < (timer_num * r->r_thread_num); i++)
    sct_internal_free_event_arrays(&r->stats[i]);
  free(r->stats);
  free(r);
}


static void reduce_threads(int thread_choice, int m_num, int r_num, int en, sct_stats_type *m_stats, sct_stats_type *r_stats) {

  if (thread_choice >= 0) {
    if (m_num-1<thread_choice || r_num != 1) sct_abort("reduce_threads: internal error", __FILE__, __LINE__);
    sct_internal_copy_stats(&m_stats[thread_choice], &r_stats[0], en);
  } else if (thread_choice == SCT_SELECT_ALL) {
    if (r_num != m_num) sct_abort("reduce_threads: internal error", __FILE__, __LINE__);
    for (int tid = 0; tid < m_num; tid++) {
      sct_internal_copy_stats(&m_stats[tid], &r_stats[tid], en);
    }
  } else if (thread_choice == SCT_REDUCE_ALL) {
    if (m_num<1 || r_num != 1) sct_abort("reduce_threads: internal error", __FILE__, __LINE__);
    {
      double t = m_stats[0].tsum; /* initialise stats with sum data from thread 0 */
      r_stats->tsum = t;
      r_stats->tmin = t; /* sum data will be reduced over thread-space, not callstats !!! */
      r_stats->tmax = t;
    }
#ifdef HAVE_LIBPAPI
    for (int ie = 0; ie < en; ie++) {
      eval_type e = m_stats[0].esum[ie];
      double r = m_stats[0].esum[ie]/m_stats[0].tsum;
      r_stats->esum[ie] = e;
      r_stats->emin[ie] = e;
      r_stats->emax[ie] = e;
      r_stats->rmin[ie] = r;
      r_stats->rmax[ie] = r;
    }
#endif
    for (int tid = 1; tid < m_num; tid++) {
      double t = m_stats[tid].tsum;
      r_stats->tsum += t;
      r_stats->tmin = MIN(r_stats->tmin, t);
      r_stats->tmax = MAX(r_stats->tmax, t);
#ifdef HAVE_LIBPAPI
      for (int ie = 0; ie < en; ie++) {
        eval_type e = m_stats[tid].esum[ie];
        double r = m_stats[tid].esum[ie]/m_stats[tid].tsum;
        r_stats->esum[ie] += e;
        r_stats->emin[ie] = MIN(r_stats->emin[ie], e);
        r_stats->emax[ie] = MAX(r_stats->emax[ie], e);
        r_stats->rmin[ie] = MIN(r_stats->rmin[ie], r);
        r_stats->rmax[ie] = MAX(r_stats->rmax[ie], r);
      }
#endif
    }
  }

}


#ifdef HAVE_MPI
static void reduce_mpi(sct_reduction_type *r) {

  const int my_debug = 0;
  int red_proc_num = r->red_proc_num;
  int proc_num = r->proc_num;
  int thread_num = r->r_thread_num;

  if (red_proc_num != 1) sct_abort("reduce_all: unexpected reduction situation", __FILE__, __LINE__);
  if (proc_num < 1) sct_abort("reduce_all: red_proc_num out of range", __FILE__, __LINE__);

  // we don't collect names
  r->proc_names = NULL;
  r->max_proc_name_len = 0;

  const int root_pid = 0;
  int is_root;
  if (r->pid == root_pid)
    is_root=1;
  else
    is_root = 0;

  int timer_num = r->global_timer_num;
  sct_stats_type (*f)[thread_num] =  (sct_stats_type (*)[thread_num]) r->stats;

#ifdef HAVE_LIBPAPI
  const int en = sct_get_event_num();

  if (en != 0) {
    long long eval[timer_num][en], red_eval[timer_num][en];
    double rval[timer_num][en], red_rval[timer_num][en];

    for (int tid = 0; tid < thread_num; tid++) {

      for (int it = 0; it < timer_num; it++) {
        for (int ie=0; ie<en; ie++) {
          eval[it][ie] = f[it][tid].esum[ie];
          rval[it][ie] = f[it][tid].esum[ie]/f[it][tid].tsum;
          /* fprintf(stderr,"reduce_mpi: it=%i, ie=%i, eval=%lld\n", it, ie, eval[it][ie]); */
          /* fprintf(stderr,"reduce_mpi: it=%i, ie=%i, rval=%e\n", it, ie, rval[it][ie]); */
        }
      }

      // reduce event values

      // sum:
      if ( (MPI_Reduce(eval, red_eval, timer_num*en, MPI_LONG_LONG, MPI_SUM, root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (e1)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          for (int ie=0; ie<en; ie++) {
            f[it][tid].esum[ie] = red_eval[it][ie];
          }
        }
      }

      // min:
      if ( (MPI_Reduce(eval, red_eval, timer_num*en, MPI_LONG_LONG, MPI_MIN, root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (e1)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          for (int ie=0; ie<en; ie++) {
            f[it][tid].emin[ie] = red_eval[it][ie];
          }
        }
      }

      // max:
      if ( (MPI_Reduce(eval, red_eval, timer_num*en, MPI_LONG_LONG, MPI_MAX, root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (e1)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          for (int ie=0; ie<en; ie++) {
            f[it][tid].emax[ie] = red_eval[it][ie];
          }
        }
      }

      // reduce event rates

      // min:
      if ( (MPI_Reduce(rval, red_rval, timer_num*en, MPI_DOUBLE, MPI_MIN, root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (e1)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          for (int ie=0; ie<en; ie++) {
            f[it][tid].rmin[ie] = red_rval[it][ie];
          }
        }
      }

      // max:
      if ( (MPI_Reduce(rval, red_rval, timer_num*en, MPI_DOUBLE, MPI_MAX, root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (e1)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          for (int ie=0; ie<en; ie++) {
            f[it][tid].rmax[ie] = red_rval[it][ie];
          }
        }
      }

    } /* tid-loop */

  } /* if en != 0 */
#endif

  {
    double tval[timer_num], red_tval[timer_num];

    for (int tid = 0; tid < thread_num; tid++) {

      // time-reduction:
      for (int it = 0; it < timer_num; it++) {
        tval[it] = f[it][tid].tsum;
        if (my_debug)
          fprintf(stderr,"reduce_mpi [%d]: it=%i, tval=%e\n", r->pid, it, tval[it]);
      }

      // sum:
      if ( (MPI_Reduce(tval, red_tval, timer_num, MPI_DOUBLE, MPI_SUM,
                       root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (t1)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          f[it][tid].tsum = red_tval[it];
          if (my_debug && r->pid == root_pid)
            fprintf(stderr,"reduce_mpi: it=%i, tid=%i, sum=%e\n", it, tid, red_tval[it]);
        }
      }

      // min:
      if ( (MPI_Reduce(tval, red_tval, timer_num, MPI_DOUBLE, MPI_MIN,
                       root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (t2)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          f[it][tid].tmin = red_tval[it];
          if (my_debug && r->pid == root_pid)
            fprintf(stderr,"reduce_mpi: it=%i, min=%e\n", it, red_tval[it]);
        }
      }

      // max:
      if ( (MPI_Reduce(tval, red_tval, timer_num, MPI_DOUBLE, MPI_MAX,
                       root_pid, r->comm) != MPI_SUCCESS) )
        sct_abort("MPI_Reduce failed (t3)", __FILE__, __LINE__);
      if (is_root) {
        for (int it = 0; it < timer_num; it++) {
          f[it][tid].tmax = red_tval[it];
          if (my_debug && r->pid == root_pid)
            fprintf(stderr,"reduce_mpi: it=%i, max=%e\n", it, red_tval[it]);
        }
      }

    } //tid-loop

  }

}
#endif


#ifdef HAVE_MPI
static void gather_mpi(sct_reduction_type *r) {
#ifdef DEBUG
  const int i_r_thread_num = 0;
  const int i_stats_size = 1;
  const int i_pid = 2;
  const int meta_rec_size = 3;
#endif

  const int my_debug = 0;

  const int max_evn = sct_get_event_num();

  int red_proc_num = r->red_proc_num;
  int proc_num = r->proc_num;
  int timer_num = sct_get_timer_num(); // this is the process local number of used timers

  if (red_proc_num != proc_num) sct_abort("gather_all: called in a mpi reduction situation", __FILE__, __LINE__);
  if (proc_num < 1) sct_abort("gather_all: red_proc_num out of range", __FILE__, __LINE__);

  // verify proc_num - we cannot afford to be wrong with that
  {
    int pnum = 0;
    if (MPI_Comm_size(r->comm, &pnum) != MPI_SUCCESS) sct_abort("MPI_Comm_size failed", __FILE__, __LINE__);
    if (pnum != proc_num)  sct_abort("internal error", __FILE__, __LINE__);
  }

  // my process name:
  const int root_pid = 0;
  int is_root;
  if (r->pid == root_pid)
    is_root = 1;
  else
    is_root = 0;
  const int r_thread_num = r->r_thread_num;

#ifdef DEBUG
  int meta_rec[meta_rec_size];
  meta_rec[i_r_thread_num] = r->r_thread_num;
  meta_rec[i_stats_size] = sizeof(sct_stats_type);
  meta_rec[i_pid] = r->pid;

  int *meta_buf;
  if (is_root) {
    if ( !(meta_buf = (int *) calloc(proc_num *  meta_rec_size, sizeof(*meta_buf))
	   ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  } else {
    meta_buf = NULL;
  }

  if (MPI_Gather(meta_rec, //void* sendbuf,
                 meta_rec_size, //int sendcount,
                 MPI_INT, //MPI_Datatype sendtype,
                 meta_buf, //void* recvbuf,
                 meta_rec_size, //int recvcount,
                 MPI_INT, //MPI_Datatype recvtype,
                 root_pid, //int root,
                 r->comm //MPI_Comm comm
                 )  != MPI_SUCCESS)
    sct_abort("MPI_Gather failed", __FILE__, __LINE__);

  // check meta data:
  if (meta_buf) {
    int (*meta)[meta_rec_size] =  (int (*)[meta_rec_size]) meta_buf;
    for (int i = 0; i < proc_num; i++) {
      if (my_debug)
        fprintf(stderr, "meta[%d][*]=[%d,%d,%d]\n", i, meta[i][i_r_thread_num],
                meta[i][i_stats_size], meta[i][i_pid]);
      if ((meta[i][i_r_thread_num] != meta_rec[i_r_thread_num]) ||
          (meta[i][i_pid] != i) ||
          (meta[i][i_stats_size] != meta_rec[i_stats_size])
          ) sct_abort("gather_all: meta data error", __FILE__, __LINE__);
    }
    if (my_debug) fprintf(stderr,"gather_all: meta test passed\n");
    free(meta_buf);
  }
#endif

  sct_stats_type *global_stats;
  if (is_root) {
    // we need to use timer_num_per_rank to alloc global mem in case of different number of timers per rank
    int timer_sum = r->timer_num_per_rank[0];
    for (int i=1; i<r->proc_num; i++) timer_sum += r->timer_num_per_rank[i];

    if ( !(global_stats = (sct_stats_type *) calloc(timer_sum *  r_thread_num, sizeof(sct_stats_type))
	   ) ) sct_abort("calloc failed", __FILE__, __LINE__);
    for (int it = 0; it < timer_sum *  r_thread_num; it++)
      sct_internal_alloc_event_arrays(&(global_stats[it]), sct_get_event_num());
  }
  else
    global_stats = NULL;

  /* --- gather hostnames --- */

  char pname[MPI_MAX_PROCESSOR_NAME];
  int pname_len = 0;
  int max_pname_len = 0;
  if (MPI_Get_processor_name(pname, &pname_len) != MPI_SUCCESS)
    sct_abort("MPI_Get_processor_name", __FILE__, __LINE__);
  if (pname_len>=MPI_MAX_PROCESSOR_NAME) sct_abort("bad process name", __FILE__, __LINE__);
  //fprintf(stderr,"[pname]=[%s], len=%i, max_len=%i\n",pname, pname_len,MPI_MAX_PROCESSOR_NAME-1);
  if (MPI_Allreduce(&pname_len, &max_pname_len, 1, MPI_INT, MPI_MAX, r->comm) != MPI_SUCCESS)
    sct_abort("MPI_Allreduce failed", __FILE__, __LINE__);
  const int max_pname_mem = max_pname_len + 1;
  //fprintf(stderr,"[max_pname_len]=%i\n",max_pname_len);
  char *global_pnames;
  if (is_root) {
    if ( !(global_pnames = (char *) calloc(proc_num * max_pname_mem, sizeof(*global_pnames))
	   ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  } else {
    global_pnames = NULL;
  }
  if (MPI_Gather(pname, //void* sendbuf,
                 max_pname_mem, //int sendcount,
                 MPI_CHAR, //MPI_Datatype sendtype,
                 global_pnames, //void* recvbuf,
                 max_pname_mem, //int recvcount,
                 MPI_CHAR, //MPI_Datatype recvtype,
                 root_pid, //int root,
                 r->comm //MPI_Comm comm
                 )  != MPI_SUCCESS)
    sct_abort("MPI_Gather failed", __FILE__, __LINE__);

  if (my_debug) {
    if (is_root) {
      char (*gpname)[max_pname_mem] =  (char (*)[max_pname_mem]) global_pnames;
      for (int iproc = 0; iproc < proc_num; iproc++) {
        fprintf(stderr,"iproc=%d, gpname=%s\n",iproc, &gpname[iproc][0]);
      }
    }
  }
  r->proc_names = global_pnames; // we give up ownership of global_names memory
  r->max_proc_name_len = max_pname_len; // we give up ownership of global_names memory

  /* ---  gather stats --- */

  /* pack sct_stats_type for MPI comm */
  int sct_stats_type_size = sct_internal_get_stats_datatype_size();
  int position;
  void *sendbuf;
  void *recvbuf;
  int sendbufsize = sct_stats_type_size * timer_num * r_thread_num;
  int recvbufsize;

  int sendbufsize_per_rank[proc_num];
  int displs[proc_num];

  if ( !(sendbuf = malloc(sendbufsize)) ) sct_abort("malloc failed", __FILE__, __LINE__);
  if (my_debug) fprintf(stderr,"sendbufsize = %d\n", sendbufsize);

  if ( MPI_Gather(&sendbufsize, 1, MPI_INT, sendbufsize_per_rank, 1, MPI_INT,
                  root_pid, r->comm) != MPI_SUCCESS )
    sct_abort("MPI_Gather failed (sendbufsize_per_rank)", __FILE__, __LINE__);

  if (is_root) {
    recvbufsize = sendbufsize_per_rank[0];
    displs[0] = 0;
    for (int i=1; i<proc_num; i++) {
      recvbufsize += sendbufsize_per_rank[i];
      displs[i] = displs[i-1]+sendbufsize_per_rank[i-1];
    }
    if (my_debug) fprintf(stderr,"recvbufsize = %d\n", recvbufsize);
    if ( !(recvbuf = malloc(recvbufsize)) ) sct_abort("malloc failed", __FILE__, __LINE__);
  }
  else {
    recvbuf = NULL;
  }

  sct_stats_type (*my_stats)[r_thread_num] = (sct_stats_type (*)[r_thread_num]) r->stats;

  if (my_debug) {
    for (int pid=0; pid<proc_num; pid++) {
      if (pid == r->pid) {
        fprintf(stderr,"timer_num=%d, r_thread_num=%d\n", timer_num, r_thread_num);
        for (int it = 0; it < timer_num; it++) {
          for (int tid =  0; tid < r_thread_num; tid++) {
            sct_stats_type *s = &my_stats[it][tid];
            fprintf(stderr,"[%d]: my_stats [%d -> %d][%d] = %d %.2e %.2e %.2e %d\n",
                    pid, it, timer_map[it], tid, s->cnum, s->tsum, s->tmin, s->tmax, s->active_under);
#ifdef HAVE_LIBPAPI
            if (max_evn != 0) {
              for (int ie = 0; ie < max_evn; ie++) {
                fprintf(stderr,"             event[%d]: esum=%lld, emin=%lld, emax=%lld, rmin=%.2e, rmax=%.2e\n",
                        ie, s->esum[ie], s->emin[ie], s->emax[ie], s->rmin[ie], s->rmax[ie]);
              }
	    }
#endif
          }
        }
      }
      fflush(stderr);
      MPI_Barrier(r->comm);
    }
  }

  /* pack MPI with tsum, last_dt, tmin, tmax, ... */
  position = 0;
  for (int it = 0; it < timer_num; it++) {
    for (int tid =  0; tid < r_thread_num; tid++) {
      sct_stats_type *s = &my_stats[it][tid];
      MPI_Pack(&(s->tsum), 4, MPI_DOUBLE, sendbuf, sendbufsize, &position, r->comm);
      MPI_Pack(&(s->cnum), 2, MPI_INT, sendbuf, sendbufsize, &position, r->comm);
#ifdef HAVE_LIBPAPI
      if (max_evn != 0) {
        MPI_Pack(&(s->esum[0]), max_evn, MPI_LONG_LONG, sendbuf, sendbufsize, &position, r->comm);
        MPI_Pack(&(s->emin[0]), max_evn, MPI_LONG_LONG, sendbuf, sendbufsize, &position, r->comm);
        MPI_Pack(&(s->emax[0]), max_evn, MPI_LONG_LONG, sendbuf, sendbufsize, &position, r->comm);
        MPI_Pack(&(s->rmin[0]), max_evn, MPI_DOUBLE, sendbuf, sendbufsize, &position, r->comm);
        MPI_Pack(&(s->rmax[0]), max_evn, MPI_DOUBLE, sendbuf, sendbufsize, &position, r->comm);
      }
#endif
    }
  }

  if (MPI_Gatherv(sendbuf, sendbufsize, MPI_PACKED,
                  recvbuf, sendbufsize_per_rank, displs, MPI_PACKED,
                  root_pid, r->comm
                 )  != MPI_SUCCESS)
    sct_abort("MPI_Gather failed", __FILE__, __LINE__);


  /* inflate root stats: */
  if (is_root) {

    /* free any previously defined reduction memory */
    if (r->stats) {
      for (int it = 0; it < timer_num * r->r_thread_num; it++)
        sct_internal_free_event_arrays(&(r->stats[it]));

      free(r->stats);
    }

    /* unpack */
    position = 0;
    int timer_offset = 0;

    for (int iproc = 0; iproc < proc_num; iproc++) {
      int timer_num_iproc = r->timer_num_per_rank[iproc];

      for (int it = 0; it < timer_num_iproc; it++) {
        for (int tid =  0; tid < r_thread_num; tid++) {

          sct_stats_type *s = &global_stats[timer_offset];

          MPI_Unpack(recvbuf, recvbufsize, &position, &(s->tsum), 4, MPI_DOUBLE, r->comm);
          MPI_Unpack(recvbuf, recvbufsize, &position, &(s->cnum), 2, MPI_INT, r->comm);
#ifdef HAVE_LIBPAPI
          if (max_evn != 0) {
            MPI_Unpack(recvbuf, recvbufsize, &position, &(s->esum[0]), max_evn, MPI_LONG_LONG, r->comm);
            MPI_Unpack(recvbuf, recvbufsize, &position, &(s->emin[0]), max_evn, MPI_LONG_LONG, r->comm);
            MPI_Unpack(recvbuf, recvbufsize, &position, &(s->emax[0]), max_evn, MPI_LONG_LONG, r->comm);
            MPI_Unpack(recvbuf, recvbufsize, &position, &(s->rmin[0]), max_evn, MPI_DOUBLE, r->comm);
            MPI_Unpack(recvbuf, recvbufsize, &position, &(s->rmax[0]), max_evn, MPI_DOUBLE, r->comm);
	  }
#endif
          timer_offset++;
	}
      }
    }

    if (my_debug) {
      int timer_sum = r->timer_num_per_rank[0];
      for (int i=1; i<r->proc_num; i++) timer_sum += r->timer_num_per_rank[i];
      if ( timer_offset != (timer_sum * r_thread_num) )
        sct_abort("unpack into global_stats failed", __FILE__, __LINE__);
    }

    /* save global stats in reduce type stats */
    r->stats = global_stats;

    if (my_debug) {
      timer_offset = 0;
      for (int iproc = 0; iproc < proc_num; iproc++) {
        int timer_num_iproc = r->timer_num_per_rank[iproc];

        for (int it = 0; it < timer_num_iproc; it++) {
          for (int tid =  0; tid < r_thread_num; tid++) {

            sct_stats_type *s = &global_stats[timer_offset];
            fprintf(stderr,"glob_stats [%d]: [%d][%d] %d %.2e %.2e %.2e %d\n",
		    iproc, it, tid, s->cnum, s->tsum, s->tmin, s->tmax, s->active_under);
#ifdef HAVE_LIBPAPI
            if (max_evn != 0) {
              for (int ie = 0; ie < max_evn; ie++) {
                fprintf(stderr,"             event[%d]: esum=%lld, emin=%lld, emax=%lld, rmin=%.2e, rmax=%.2e\n",
                        ie, s->esum[ie], s->emin[ie], s->emax[ie], s->rmin[ie], s->rmax[ie]);
	      }
	    }
#endif
            timer_offset++;
          }
        }
      }
    }

    if (recvbuf) free(recvbuf);
  }
  if (my_debug) {
    fflush(stderr);
    MPI_Barrier(r->comm);
  }
  if (sendbuf) free(sendbuf);

}
#endif


sct_reduction_type *sct_reduction_new(int context_choice, int proc_choice, int thread_choice, int sp_merging) {

  const int my_debug = 0;

  // check context:
  const int max_context_num = sct_get_context_num();
  if (context_choice < 0 || context_choice >= max_context_num)
    sct_abort("sct_reduction_new: invalid context_choice", __FILE__, __LINE__);

  sct_context_type *con = sct_get_context(context_choice);
  if (!con) sct_abort("sct_reduction_new: invalid context", __FILE__, __LINE__);

  // check timer space:
  int timer_num = sct_get_timer_num();
  if (timer_num < 1) return NULL;
  const int global_timer_num = sct_get_global_timer_num();

  // prepare return value and set some meta info
  sct_reduction_type *res;
  if ( !(res = (sct_reduction_type *) calloc(1, sizeof(sct_reduction_type))
	 ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  res->context_choice = context_choice;
  res->proc_choice = proc_choice;
  res->thread_choice = thread_choice;
  res->sp_merging = sp_merging;
  res->global_timer_num = global_timer_num;
  res->context_name = sct_string_zero;
  sct_string_copy(&con->name, &res->context_name);
  const int callstats = sct_get_callstats();

  // processes:
  int proc_num;
  int red_proc_num;
  const int root_pid = 0;
  int is_root;
#ifdef HAVE_MPI
  proc_num = con->procnum;
  if (proc_choice >= 0) {
    if (proc_choice >= proc_num) sct_abort("sct_reduction_new: invalid proc_choice", __FILE__, __LINE__);
    red_proc_num = 1;
  } else if (proc_choice == SCT_REDUCE_ALL) {
    red_proc_num = 1;
  } else if (proc_choice == SCT_SELECT_ALL) {
    red_proc_num = proc_num;
  }
  res->pid = con->pid;
  res->comm = con->comm;

  if (res->pid == root_pid)
    is_root = 1;
  else
    is_root = 0;
#else
  proc_num = 1;
  red_proc_num = 1;
  res->pid = 0;
  is_root = 1;
#endif
  res->proc_num = proc_num;
  res->red_proc_num = red_proc_num;

  // timer:
  if (is_root) {
    if ( !(res->timer_num_per_rank = (int *) calloc(proc_num, sizeof(int))
          ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  }
  else {
    res->timer_num_per_rank = NULL;
    res->timer_map_per_rank = NULL;
  }
#ifdef HAVE_MPI
  // gather number of timers per rank
  if (MPI_Gather(&timer_num, 1, MPI_INT, res->timer_num_per_rank, 1, MPI_INT,
                 root_pid, res->comm) != MPI_SUCCESS)
    sct_abort("MPI_Gather failed (timer_num_per_rank)", __FILE__, __LINE__);
  int displs[proc_num];
  if (is_root) {
    int timer_sum = res->timer_num_per_rank[0];
    displs[0] = 0;
    for (int ip=1; ip<proc_num; ip++) {
      timer_sum += res->timer_num_per_rank[ip];
      displs[ip] = displs[ip-1] + res->timer_num_per_rank[ip-1];
    }
    if ( !(res->timer_map_per_rank = (int *) calloc(timer_sum, sizeof(int))
          ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  }
  // gather timer_map from each rank
  if (MPI_Gatherv(timer_map, timer_num, MPI_INT,
                  res->timer_map_per_rank, res->timer_num_per_rank, displs, MPI_INT,
                  root_pid, res->comm) != MPI_SUCCESS)
    sct_abort("MPI_Gatherv failed (timer_map_per_rank)", __FILE__, __LINE__);
  if (my_debug && is_root) {
    int c=0;
    for (int ip=0; ip<proc_num; ip++) {
      for (int it=0; it<res->timer_num_per_rank[ip]; it++) {
	printf("DEBUG [%d] %d -> %d\n",ip,it,res->timer_map_per_rank[c]);
	c++;
      }
    }
  }
#else
  res->timer_num_per_rank[0] = timer_num;
  //  if ( !(res->timer_map_per_rank = (int *) calloc(timer_num, sizeof(int))
	 //        ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  //  timer_map = calloc(timer_num, sizeof(int));
  //  memcpy(res->timer_map_per_rank, timer_map, timer_num*sizeof(int));
  res->timer_map_per_rank = NULL;
#endif

  // threads:
  int p_num; // parallel region thread number
  int m_num; // thread number after reduction over sp-dim
  int red_thread_num; // number of threads after reduction along thread dim
#ifdef _OPENMP
  if (omp_in_parallel()) sct_abort("sct_reduction_new: call within parallel region not allowed", __FILE__, __LINE__);
  if (proc_choice == SCT_REDUCE_ALL && thread_choice == SCT_REDUCE_ALL)
    sct_abort("sct_reduction_new: Combined reduction over threads and procs is not supported.", __FILE__, __LINE__);
  p_num = sct_get_pr_thread_count(); // number of threads in parallel region
  if (sp_merging == 0)  {
    sct_abort("sct_reduction_new: (sp_merging == 0) is not supported", __FILE__, __LINE__);
  } else if ( sp_merging == SCT_SP_SERIAL_ONLY)  {
    m_num = 1; // just the serial part
    if (thread_choice == 0 || thread_choice == SCT_REDUCE_ALL || thread_choice == SCT_SELECT_ALL) {
      red_thread_num = 1; // serial
    } else {
      sct_abort("sct_reduction_new: cannot select thread>0 within serial phase", __FILE__, __LINE__);
    }
  } else if ( sp_merging == SCT_SP_PARALLEL_ONLY ||  sp_merging == SCT_SP_MERGE_SIMPLE)  {
    m_num = p_num; // parallel region part, possibly combined with serial part
    if (thread_choice >= 0) {
      if (thread_choice >= p_num) sct_abort("sct_reduction_new: invalid thread_choice", __FILE__, __LINE__);
      red_thread_num = 1;
    } else if (thread_choice >= 0 || thread_choice == SCT_REDUCE_ALL) {
      red_thread_num = 1; // reduce merged set
    } else if (thread_choice == SCT_SELECT_ALL) {
      red_thread_num = m_num;
    } else {
      sct_abort("sct_reduction_new: cannot select thread>0 within serial phase", __FILE__, __LINE__);
    }
  }
  else if ( sp_merging == SCT_SP_SELECT_ALL ) {
    m_num = p_num + 1; // parallel region part and serial part separately
    if (thread_choice == SCT_SELECT_ALL) {
      red_thread_num = m_num;
    } else {
      sct_abort("sct_reduction_new: cannot select thread or reduce over threads if sp_merging is SCT_SP_SELECT_ALL",
                __FILE__, __LINE__);
    }
  }
  else {
    sct_abort("sct_reduction_new: invalid sp_merging choice", __FILE__, __LINE__);
  }
#else
  p_num = 0;
  m_num = 1;
  red_thread_num = 1;
#endif
  res->p_thread_num = p_num;
  res->m_thread_num = m_num;
  res->r_thread_num = red_thread_num;


  if (my_debug)
    printf("sct_reduction_new: proc_num=%d, red_proc_num=%d, p_num=%d, m_num=%d, red_thread_num=%d\n",
	   proc_num, red_proc_num, p_num, m_num, red_thread_num );

  const int sp_num = p_num + 1;
  int en = sct_get_event_num();
  res->event_num = en;

#ifdef _OPENMP
  //sct_stats_type sp_stats[timer_num][sp_num];
  sct_stats_type *sp_stats_mem;
  if ( !(sp_stats_mem = (sct_stats_type *) calloc(timer_num * sp_num, sizeof(sct_stats_type))
	 ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  for (int it = 0; it < timer_num *  sp_num; it++)
    sct_internal_alloc_event_arrays(&(sp_stats_mem[it]), en);
  sct_stats_type (*sp_stats)[sp_num] =  (sct_stats_type (*)[sp_num]) sp_stats_mem;

  //sct_stats_type m_stats[timer_num][m_num];
  sct_stats_type *m_stats_mem;
  if ( !(m_stats_mem = (sct_stats_type *) calloc(timer_num * m_num, sizeof(sct_stats_type))
	 ) ) sct_abort("calloc failed", __FILE__, __LINE__);
  for (int it = 0; it < timer_num *  m_num; it++)
    sct_internal_alloc_event_arrays(&(m_stats_mem[it]), en);
  sct_stats_type (*m_stats)[m_num] =  (sct_stats_type (*)[m_num]) m_stats_mem;
#endif

  //sct_stats_type red_local[timer_num][red_thread_num];
  sct_stats_type *red_local_mem;
  if (red_proc_num == 1) { // SCT_REDUCE_ALL
    if ( !(red_local_mem = (sct_stats_type *) calloc(global_timer_num * red_thread_num,
                                                     sizeof(sct_stats_type))
          ) ) sct_abort("calloc failed", __FILE__, __LINE__);
    for (int it = 0; it < global_timer_num *  red_thread_num; it++)
      sct_internal_alloc_event_arrays(&(red_local_mem[it]), en);
  }
  else { // SCT_SELECT_ALL ... uses proc local number of timers
    if ( !(red_local_mem = (sct_stats_type *) calloc(timer_num * red_thread_num,
                                                     sizeof(sct_stats_type))
          ) ) sct_abort("calloc failed", __FILE__, __LINE__);
    for (int it = 0; it < timer_num *  red_thread_num; it++)
      sct_internal_alloc_event_arrays(&(red_local_mem[it]), en);
  }
  sct_stats_type (*red_local)[red_thread_num] = (sct_stats_type (*)[red_thread_num]) red_local_mem;

#ifdef _OPENMP
  get_sp(con, proc_choice, sp_num, sp_stats_mem);

  if (my_debug) {
    for (int it = 0; it < timer_num; it++)
      for (int tid = 0; tid < sp_num; tid++)
        printf("debug I: sp_stats[%d][%d] = %d , %.2e , %d\n", it, tid,
	       sp_stats[it][tid].cnum, sp_stats[it][tid].tsum,
	       sp_stats[it][tid].active_under);
  }

  merge_sp(sp_num, m_num, sp_merging, sp_stats_mem, m_stats_mem);

  if (my_debug) {
    for (int it = 0; it < timer_num; it++)
      for (int tid = 0; tid < m_num; tid++)
        printf("debug II: m_stats[%d][%d] = %d , %.2e , %d\n", it, tid,
	       m_stats[it][tid].cnum, m_stats[it][tid].tsum,
	       m_stats[it][tid].active_under);
  }

  for (int it = 0; it<timer_num; it++) {
    if (&m_stats[it][0] != &m_stats_mem[it*m_num+0])
      sct_abort("sct_reduction_new: internal error", __FILE__, __LINE__);
    if (&sp_stats[it][0] != &sp_stats_mem[it*sp_num+0])
      sct_abort("sct_reduction_new: internal error", __FILE__, __LINE__);
    reduce_threads(thread_choice, m_num, red_thread_num, en, &m_stats[it][0], &red_local[it][0]);
  }
  cleanup_stats(red_thread_num, red_local_mem);

  if (my_debug) {
    for (int it = 0; it < timer_num; it++)
      for (int tid = 0; tid < red_thread_num; tid++)
        printf("debug III: red_local[%d][%d] = %d , %.2e , %d\n", it, tid,
	       red_local[it][tid].cnum, red_local[it][tid].tsum,
	       red_local[it][tid].active_under);
  }

#else
  if (sp_num != 1 || red_thread_num != 1) sct_abort("sct_reduction_new: internal error", __FILE__, __LINE__);
  get_sp(con, proc_choice, red_thread_num, red_local_mem);

  if (my_debug) {
    int tnum;
    if (red_proc_num == 1)
      tnum = global_timer_num;
    else
      tnum = timer_num;

    for (int it = 0; it < tnum; it++)
      for (int tid = 0; tid < red_thread_num; tid++)
        printf("[%d]: red_local[%d][%d] = %d %.2e %.2e %.2e %d\n", con->pid, it, tid,
               red_local[it][tid].cnum, red_local[it][tid].tsum,
               red_local[it][tid].tmin, red_local[it][tid].tmax,
               red_local[it][tid].active_under);
  }

  cleanup_stats(red_thread_num, red_local_mem);

#endif

#ifdef HAVE_MPI
  res->stats = red_local_mem;

  if (proc_choice >= 0) {
    sct_abort("sct_reduction_new: single mpi task selection not implemented yet", __FILE__, __LINE__);
  } else if (proc_choice == SCT_SELECT_ALL) {
    gather_mpi(res); // this enlarges res->stats at root_pid
  } else if (proc_choice == SCT_REDUCE_ALL) {
    reduce_mpi(res); // this overwrites res->stats by reduces statistics
  } else {
    sct_abort("sct_reduction_new: unexpected case", __FILE__, __LINE__);
  }
#else
  res->stats = (sct_stats_type*)red_local;
#endif

#ifdef _OPENMP
  for (int it = 0; it < timer_num *  sp_num; it++)
    sct_internal_free_event_arrays(&sp_stats_mem[it]);
  free(sp_stats_mem);

  for (int it = 0; it < timer_num *  m_num; it++)
    sct_internal_free_event_arrays(&m_stats_mem[it]);
  free(m_stats_mem);
#endif

  return res;
}

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
 * \file sct_collector.c
 * measure and collect time and events
 */

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <limits.h>

#include "sct_mach.h"
#include "sct_collector.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// internal function declarations:

static void time_sec_str(const double t, char *s, const int sn);
static void abort_timer(const int it, const char *reason, const char *fname, const int line);
static inline double dmax(double a, double b);

#ifdef NESTED_TIMER
static void update_timer_stack_atstart(sct_stack_type *s, sct_stats_type *stats, const int it);
static void update_timer_stack_atstop(sct_stack_type *s, const int it);
#endif

// real timer states:
static const int sct_undef_state = 0;
static const int sct_on_state    = 1;
static const int sct_off_state   = 2;

static eset_type *eset;

// number of used timers:
static int timer_num = 0;

// depth of timer nesting
static int timer_nest_depth = 0;

// MPI-parallelization data:
static const int root = 0;

// meta data:
static sct_meta_type *meta_pool;

// attribute table:
static sct_attribute_table *attribute_table = NULL;

// context free thread private marks:
/// \todo Threadprivate pragma can be replaced by double allocation, access might be faster for some compilers
static sct_mark_type *mark_pool = NULL; // mark_pool[itimer]
#ifdef _OPENMP
#pragma omp threadprivate (mark_pool)
#endif

static int static_pr_thread_count = 0;
static int static_pr_thread_id = -1;
#ifdef _OPENMP
#pragma omp threadprivate (static_pr_thread_id)
#endif

// contexts:
static sct_context_type *context_pool = NULL;
static int context_poolsize = 0;
static int context_num = 0;

static int default_icontext = -1;
static int icontext = -1;
static int icontext_stack[SCT_MAX_CONTEXT_DEPTH+1];
static int icontext_stack_pos = 0;

// index range for pool data: [0 ... timer_size-1]
static int timer_size = 0;

// usage of callstats
static int callstats = -1;

// usage of eventcounters
static int eventcounters = -1;

// usage of nestedtimers
static int nestedtimers = -1;

#ifdef HAVE_MPI
static MPI_Datatype stats_mpi_datatype = MPI_DATATYPE_NULL;
static MPI_Comm sct_comm_self = MPI_COMM_NULL;
static int stats_mpi_datatype_size = 0;
#else
static sct_datatype stats_mpi_datatype = 0;
static sct_comm_type sct_comm_self = 0;
static int stats_mpi_datatype_size = 0;
#endif

static void init_datatypes();
static void abort_timer(const int it, const char *reason, const char *fname, const int line) ATTRIBUTE_NORETURN;

// --- implementation ---

sct_datatype sct_internal_get_stats_datatype() {
#ifdef HAVE_MPI
  if (stats_mpi_datatype == MPI_DATATYPE_NULL) init_datatypes();
  return stats_mpi_datatype;
#else
  sct_abort("sct_internal_get_stats_datatype: not defined without HAVE_MPI", __FILE__, __LINE__);
#endif
}

int sct_internal_get_stats_datatype_size() {
#ifdef HAVE_MPI
  if (stats_mpi_datatype == MPI_DATATYPE_NULL) init_datatypes();
  return stats_mpi_datatype_size;
#else
  sct_abort("sct_internal_get_stats_datatype_size: not defined without HAVE_MPI", __FILE__, __LINE__);
#endif
}


void sct_internal_alloc_event_arrays(sct_stats_type *s, int en) {
#ifdef HAVE_LIBPAPI
  if (!s) return;

  if (en == 0) {
    s->esum = NULL;
    s->emin = NULL;
    s->emax = NULL;
    s->rmin = NULL;
    s->rmax = NULL;
  }
  else {
    if (   !(s->esum = calloc(en, sizeof(eval_type)))
        || !(s->emin = calloc(en, sizeof(eval_type)))
        || !(s->emax = calloc(en, sizeof(eval_type)))
        || !(s->rmin = calloc(en, sizeof(eval_type)))
        || !(s->rmax = calloc(en, sizeof(eval_type)))
       ) sct_abort("memory allocation failed", __FILE__, __LINE__);
  }
#endif
}


void sct_internal_free_event_arrays(sct_stats_type *s) {
#ifdef HAVE_LIBPAPI
  if ( s->esum ) free(s->esum);
  if ( s->emin ) free(s->emin);
  if ( s->emax ) free(s->emax);
  if ( s->rmin ) free(s->rmin);
  if ( s->rmax ) free(s->rmax);
#endif
}


int sct_string_recap(sct_string_type *v, int mincap) {

  if (mincap<0) mincap = 0;
  int new_cap = 8*( ( (mincap+1)/8) + 1) - 1; //n*8 -1 with n>0
  if (new_cap < mincap) sct_abort("sct_string_recap: internal error", __FILE__, __LINE__);
  if (new_cap < v->cap) return v->cap;

  if (!(v->cs=realloc(v->cs, new_cap+1))) sct_abort("sct_string_recap: realloc failed", __FILE__, __LINE__);
  v->cap = new_cap;

  return v->cap;
}

sct_string_type sct_string_new(const char *cstring) {
  sct_string_type v = sct_string_zero;
  if (cstring) {
    v.n = strlen(cstring);
    if (!v.n) return v;
    if (v.n>SCT_LABEL_SIZE) sct_abort("sct_string_new: unexpected very long string", __FILE__, __LINE__);
    sct_string_recap(&v, v.n);
    memcpy(v.cs, cstring, v.n);
  } else {
    v.n = 0;
    sct_string_recap(&v, v.n);
  }

  v.cs[v.n] = 0;
  return v;
}

void sct_string_copy(const sct_string_type *src, sct_string_type *dst) {
  if (!src || !src->n) {
    dst->n = 0;
    sct_string_recap(dst, 0);
    dst->cs[0]=0;
    return;
  }

  sct_string_recap(dst, src->n);
  memcpy(dst->cs, src->cs, src->n);
  dst->n = src->n;
  dst->cs[dst->n] = 0;
  //printf("sct_string_copy: [src]=[%s], [dst]=[%s]\n",src->cs, dst->cs);
}

void sct_string_delete(sct_string_type *s) {
  if (s) {
    if (s->cs) free(s->cs);
    s->cap = 0;
    s->n = 0;
    s->cs = NULL;
  }
}

int sct_get_callstats() {
  if (callstats>=0) return callstats;
  char *s = getenv ("SCT_CALLSTATS");
  if (!s) return -1;
  if ( !strcmp(s,"0") || !strcasecmp(s,"SCT_WITHOUT_CALLSTATS") ) return 0;
  if ( !strcmp(s,"1") || !strcasecmp(s,"SCT_WITH_CALLSTATS") ) return 1;
  return -1;
}

void sct_set_callstats(const int val) {
  callstats = 1 ? (val!=0) : 0;
}

int sct_get_eventcounters() {
  if (eventcounters>=0) return eventcounters;
  char *s = getenv ("SCT_EVENTCOUNTERS");
  if (!s) return -1;
  if ( !strcmp(s,"0") || !strcasecmp(s,"SCT_WITHOUT_EVENTCOUNTERS") ) return 0;
  if ( !strcmp(s,"1") || !strcasecmp(s,"SCT_WITH_EVENTCOUNTERS") ) return 1;
  if ( !strcmp(s,"2") || !strcasecmp(s,"SCT_WITH_EVENTRATES") ) return 2;
  return -1;
}

void sct_set_eventcounters(const int val) {
  eventcounters = 1 ? (val!=0) : 0;
}

int sct_get_nestedtimers() {
  if (nestedtimers>=0) return nestedtimers;
  char *s = getenv ("SCT_NESTEDTIMERS");
  if (!s) return -1;
  if ( !strcmp(s,"0") || !strcasecmp(s,"SCT_WITHOUT_NESTEDTIMERS") ) return 0;
  if ( !strcmp(s,"1") || !strcasecmp(s,"SCT_WITH_NESTEDTIMERS") ) return 1;
  return -1;
}

void sct_set_nestedtimers(const int val) {
  nestedtimers = 1 ? (val!=0) : 0;
}

int sct_get_timer_num() {
  return timer_num;
}

int sct_get_pr_thread_count() {
  // returns in the number of threads used in a parallel region
  // notice that a serial run of an OMP-enabled program is treated differently than a serial program without OMP
#ifdef _OPENMP
  return static_pr_thread_count;
#else
  return 0; // we count the {parallel region threads} and the {single serial phase thread}  separately
#endif
}

int sct_get_pr_thread_id() {
#ifdef _OPENMP
  // check if we step on uninitialized thread private data or if we are beyond our initial thread count;
  if (static_pr_thread_id<0 || static_pr_thread_id>=static_pr_thread_count) sct_abort("sct_get_pr_thread_id: static_pr_thread_id exceeds static range", __FILE__, __LINE__);
  return static_pr_thread_id;
#else
  sct_abort("sct_get_pr_thread_id: call outside parallel region is forbidden", __FILE__, __LINE__);
#endif
}

int sct_get_event_num() {
#ifdef HAVE_LIBPAPI
  return (eventcounters) ? eset->en : 0;
#else
  return 0;
#endif
}

int sct_internal_get_max_name_length() {
  int mnl = 0;

  for (int it = 0; it < timer_num; it++) {
    if (meta_pool[it].used) {
      mnl = MAX(mnl, meta_pool[it].name.n + timer_nest_depth );
    }
  }
  return mnl;
}

int sct_internal_get_name_length(const int it) {
  if (it<0 || it>=timer_num) sct_abort("timer out of bounds", __FILE__, __LINE__);
  if (!meta_pool[it].used) sct_abort("access to meta data of unused timer", __FILE__, __LINE__);

  return meta_pool[it].name.n;
}

char *sct_get_timer_cname(const int it) {
  if (it<0 || it>=timer_num) sct_abort("timer out of bounds", __FILE__, __LINE__);
  if (!meta_pool[it].used) sct_abort("access to meta data of unused timer", __FILE__, __LINE__);
  return meta_pool[it].name.cs;
}

char *sct_get_event_cname(const int ie) {
#ifdef HAVE_LIBPAPI
  if (ie<0 || ie>=eset->en) sct_abort("event out of bounds", __FILE__, __LINE__);
  return eset->name[ie];
#else
  return NULL;
#endif
}


// attributes:

sct_attribute_table *sct_internal_create_attribute_table(void) {
  sct_attribute_table *t;

  if ( !(t = (sct_attribute_table*) malloc(sizeof(sct_attribute_table))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);
  t->first = NULL;
  return t;
}

const sct_attribute_table *sct_internal_get_attribute_table(void) {
  return attribute_table;
}

void sct_add_report_attribute_int(const char *key, int val) {
  if (attribute_table == NULL) attribute_table = sct_internal_create_attribute_table();

  sct_attribute_type *a;
  if ( !(a = (sct_attribute_type*) malloc(sizeof(sct_attribute_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);

  a->type = SCT_INT;

  char *k = (char*) malloc(strlen(key)+1);
  memcpy(k, key, strlen(key)+1);
  a->key = k;

  int *v = (int*) malloc(sizeof(int));
  *v = val;
  a->value = v;

  a->next = attribute_table->first;
  attribute_table->first = a;
}

void sct_add_report_attribute_long(const char *key, long val) {
  if (attribute_table == NULL) attribute_table = sct_internal_create_attribute_table();

  sct_attribute_type *a;
  if ( !(a = (sct_attribute_type*) malloc(sizeof(sct_attribute_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);

  a->type = SCT_LONG;

  char *k = (char*) malloc(strlen(key)+1);
  memcpy(k, key, strlen(key)+1);
  a->key = k;

  long *v = (long*) malloc(sizeof(long));
  *v = val;
  a->value = v;

  a->next = attribute_table->first;
  attribute_table->first = a;
}

void sct_add_report_attribute_float(const char *key, float val) {
  if (attribute_table == NULL) attribute_table = sct_internal_create_attribute_table();

  sct_attribute_type *a;
  if ( !(a = (sct_attribute_type*) malloc(sizeof(sct_attribute_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);

  a->type = SCT_FLOAT;

  char *k = (char*) malloc(strlen(key)+1);
  memcpy(k, key, strlen(key)+1);
  a->key = k;

  float *v = (float*) malloc(sizeof(float));
  *v = val;
  a->value = v;

  a->next = attribute_table->first;
  attribute_table->first = a;
}

void sct_add_report_attribute_double(const char *key, double val) {
  if (attribute_table == NULL) attribute_table = sct_internal_create_attribute_table();

  sct_attribute_type *a;
  if ( !(a = (sct_attribute_type*) malloc(sizeof(sct_attribute_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);

  a->type = SCT_DOUBLE;

  char *k = (char*) malloc(strlen(key)+1);
  memcpy(k, key, strlen(key)+1);
  a->key = k;

  double *v = (double*) malloc(sizeof(double));
  *v = val;
  a->value = v;

  a->next = attribute_table->first;
  attribute_table->first = a;
}

void sct_add_report_attribute_string(const char *key, char *val) {
  if (attribute_table == NULL) attribute_table = sct_internal_create_attribute_table();

  sct_attribute_type *a;
  if ( !(a = (sct_attribute_type*) malloc(sizeof(sct_attribute_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);

  a->type = SCT_STRING;

  char *k = (char*) malloc(strlen(key)+1);
  memcpy(k, key, strlen(key)+1);
  a->key = k;

  char *v = (char*) malloc(strlen(val)+1);
  memcpy(v, val, strlen(val)+1);
  a->value = v;

  a->next = attribute_table->first;
  attribute_table->first = a;
}

void sct_internal_free_attribute_table(void) {
  if (attribute_table != NULL) {
    sct_attribute_type *a, *nexta;
    for (a = attribute_table->first; nexta != NULL; a = nexta) {
      nexta = a->next;
      free(a->key);
      free(a->value);
      free(a);
    }
    free(attribute_table);
  }
}


// contexts:

int sct_get_context_num() {
  return context_num;
}


sct_context_type *sct_get_context(const int icon) {
  if (icon<0 || icon>=context_num) return NULL;
  if (!context_pool[icon].valid) return NULL;
  return &context_pool[icon];
}


void sct_context_start(const int icon) {
#ifdef CHECK_TIMER
  if (icon<0 || icon>=context_num) sct_abort("sct_context_start: context out of bounds", __FILE__, __LINE__);
  if (icontext_stack_pos >= SCT_MAX_CONTEXT_DEPTH-1) sct_abort("sct_context_start: context too deep", __FILE__, __LINE__);
  //  if (context_pool[icon].active) sct_abort("sct_context_start: context already active", __FILE__, __LINE__);
#endif
  context_pool[icon].active = 1;

#ifdef NESTED_TIMER
#ifdef _OPENMP
  const int p_num = sct_get_pr_thread_count();
#pragma omp parallel
  {
    int tid;
    if ( (tid = sct_get_pr_thread_id()) >= p_num) sct_abort("internal error", __FILE__, __LINE__);
    context_pool[icon].active_timer_p[tid].stack[0] = -1;
    context_pool[icon].active_timer_p[tid].pos = 0;
  }
#endif
  context_pool[icon].active_timer.stack[0] = -1;
  context_pool[icon].active_timer.pos = 0;
#endif

  icontext = icon;
  icontext_stack_pos++;
  icontext_stack[icontext_stack_pos]=icontext;
}


void sct_context_stop(const int icon) {
#ifdef CHECK_TIMER
  if (icon<0 || icon>=context_num) sct_abort("sct_context_stop: context out of bounds", __FILE__, __LINE__);
  if (!context_pool[icon].active) sct_abort("sct_context_stop: context not active", __FILE__, __LINE__);
  if (icon != icontext) sct_abort("sct_context_stop: contexts badly nested", __FILE__, __LINE__);
  if (icontext_stack_pos<1) sct_abort("sct_context_stop: internal error", __FILE__, __LINE__);
  if (icontext_stack[icontext_stack_pos] != icontext) sct_abort("sct_context_stop: internal error", __FILE__, __LINE__);
#endif
  context_pool[icon].active = 0;
  icontext_stack_pos--;
  icontext = icontext_stack[icontext_stack_pos];
}


int sct_new_global_context(const char* name) {
#ifdef HAVE_MPI
  return sct_new_context(name, MPI_COMM_WORLD);
#else
  return sct_new_context(name, sct_comm_self);
#endif
}

int sct_new_context(const char* name, sct_comm_type comm) {

  // returns handle for a new context;
  // handle is always > 0 since 0 is used for the default_context

#ifdef _OPENMP
  if (omp_in_parallel())
    sct_abort("ERROR: Cannot create new context in a parallel region",__FILE__, __LINE__);
#endif

  if (!timer_size) sct_abort("ERROR: missing sct_init call", __FILE__, __LINE__);

  if (context_num > context_poolsize) sct_abort("internal error", __FILE__, __LINE__);

  if (context_num == context_poolsize) {
    // we run out out of pre-allocated resources and need to realloc our pool
    if (!context_poolsize)
      context_poolsize = 1;
    else
      context_poolsize *= 2;

    if ( !(context_pool = realloc(context_pool, sizeof(sct_context_type)*context_poolsize)) )
      sct_abort("memory allocation failed", __FILE__, __LINE__);
    for (int ic = context_num; ic < context_poolsize; ic++) {
      // zero new memory since we do not have recalloc
      sct_context_type *c = &context_pool[ic];
      c->valid = 0;
      c->name = (sct_string_type) {.n=0, .cs=NULL};
#ifdef NESTED_TIMER
      c->active_timer.stack = NULL;
      c->active_timer.pos = 0;
      c->active_timer.n = 0;
#endif
#ifdef HAVE_MPI
      c->comm = MPI_COMM_NULL;
#endif
      c->procnum = 0;
      c->pid = 0;
#ifdef _OPENMP
#ifdef NESTED_TIMER
      c->active_timer_p = NULL;
#endif
      c->p = NULL;
#endif
      c->s = NULL;
    }
  }

  // new context:
  context_num++;
  int icon = context_num-1;
  sct_context_type *con = &context_pool[icon];

  con->name = sct_string_new(name);
  con->valid = 1;
#ifdef HAVE_MPI
  //con->comm = comm;
  if (MPI_Comm_dup(comm, &con->comm) != MPI_SUCCESS) sct_abort("MPI_Comm_dup failed", __FILE__, __LINE__);
  //printf("sct_new_context: con->comm = %d\n",con->comm);
  if (MPI_Comm_size(con->comm, &con->procnum) != MPI_SUCCESS) sct_abort("MPI_Comm_size failed", __FILE__, __LINE__);
  if (MPI_Comm_rank(con->comm, &con->pid) != MPI_SUCCESS) sct_abort("MPI_Comm_rank failed", __FILE__, __LINE__);
  if (con->procnum < 1 || con->pid<0) sct_abort("internal error", __FILE__, __LINE__);
#else
  con->procnum = 1;
  con->pid = 0;
#endif

  // serial measurement data:
  if ( !(con->s = calloc(timer_size, sizeof(sct_stats_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);

  for (int it=0; it<timer_size; it++) {
    con->s[it].tmin = SCT_MAX_TVAL;
    con->s[it].active_under = -2;
#ifdef HAVE_LIBPAPI
    // allocate event arrays
    sct_internal_alloc_event_arrays(&(con->s[it]), sct_get_event_num());
    // init arrays of min values
    for (int i=0; i<sct_get_event_num(); i++) {
      con->s[it].emin[i] = SCT_MAX_EVAL;
      con->s[it].rmin[i] = SCT_MAX_TVAL;
    }
#endif
  }

#ifdef NESTED_TIMER
  if ( !(con->active_timer.stack = calloc(SCT_MAX_NEST_DEPTH, sizeof(int))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);
  con->active_timer.stack[0] = -1;
  con->active_timer.pos = 0;
  con->active_timer.n = SCT_MAX_NEST_DEPTH;
  timer_nest_depth = MAX(timer_nest_depth, con->active_timer.n);
#endif


#ifdef _OPENMP
  const int p_num = sct_get_pr_thread_count();
  // shared pointers to threadprivate space, major dimension:
  if ( !(con->p = calloc(p_num, sizeof(sct_stats_type*))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);
#ifdef NESTED_TIMER
  if ( !(con->active_timer_p = calloc(p_num, sizeof(sct_stack_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);
#endif

  // threadprivate space, minor dimension:
#pragma omp parallel
  {
    int tid;
    if ( (tid = sct_get_pr_thread_id()) >= p_num) sct_abort("internal error", __FILE__, __LINE__);

#pragma omp critical
    {

#ifdef NESTED_TIMER
      if ( !(con->active_timer_p[tid].stack = calloc(SCT_MAX_NEST_DEPTH, sizeof(int))) )
        sct_abort("memory allocation failed", __FILE__, __LINE__);
      con->active_timer_p[tid].stack[0] = -1;
      con->active_timer_p[tid].pos = 0;
      con->active_timer_p[tid].n = SCT_MAX_NEST_DEPTH;
      timer_nest_depth = MAX(timer_nest_depth, con->active_timer_p[tid].n);
#endif

      if ( !(con->p[tid]  = calloc(timer_size, sizeof(sct_stats_type))) )
        sct_abort("memory allocation failed", __FILE__, __LINE__);

      for (int it=0; it<timer_size; it++) {
        con->p[tid][it].active_under = -2;
        con->p[tid][it].tmin = SCT_MAX_TVAL;

#ifdef HAVE_LIBPAPI
	// allocate event arrays
	sct_internal_alloc_event_arrays(&(con->p[tid][it]), sct_get_event_num());
	// init arrays of min values
	for (int i=0; i<sct_get_event_num(); i++) {
	  con->p[tid][it].emin[i] = SCT_MAX_EVAL;
	  con->p[tid][it].rmin[i] = SCT_MAX_TVAL;
	}
#endif
      }
    }
  }
#endif

  return icon;
}

static int stats_are_equal(sct_stats_type *a, sct_stats_type *b, int en) {
  if (
      (a->tsum != b->tsum) ||
      (a->last_dt != b->last_dt) ||
      (a->tmin != b->tmin) ||
      (a->tmax != b->tmax) ||
      (a->cnum != b->cnum) ||
      (a->active_under != b->active_under)
      ) return 0;

#ifdef HAVE_LIBPAPI
  for(int ie = 0; ie < en; ie++) {
    if (
        (a->esum[ie] != b->esum[ie]) ||
        (a->emin[ie] != b->emin[ie]) ||
        (a->emax[ie] != b->emax[ie]) ||
        (a->rmin[ie] != b->rmin[ie]) ||
        (a->rmax[ie] != b->rmax[ie])
        ) return 0;
  }
#endif
  return 1;
}

void sct_internal_copy_stats(sct_stats_type *src, sct_stats_type *dst, int en) {
  dst->tsum = src->tsum;
  dst->last_dt = src->last_dt;
  dst->tmin = src->tmin;
  dst->tmax = src->tmax;
  dst->cnum = src->cnum;
  dst->active_under = src->active_under;
#ifdef HAVE_LIBPAPI
  for (int ie=0; ie<en; ie++) {
    dst->esum[ie] = src->esum[ie];
    dst->emin[ie] = src->emin[ie];
    dst->emax[ie] = src->emax[ie];
    dst->rmin[ie] = src->rmin[ie];
    dst->rmax[ie] = src->rmax[ie];
  }
#endif
}

static void sct_internal_zero_stats(sct_stats_type *s, int en) {
  s->tsum = 0.0;
  s->last_dt = 0.0;
  s->tmin = 0.0;
  s->tmax = 0.0;
  s->cnum = 0;
  s->active_under = -2;
#ifdef HAVE_LIBPAPI
  for (int ie=0; ie<en; ie++) {
    s->esum[ie] = 0;
    s->emin[ie] = 0;
    s->emax[ie] = 0;
    s->rmin[ie] = 0.0;
    s->rmax[ie] = 0.0;
  }
#endif
}

#ifdef HAVE_MPI
static void init_datatypes() {

  if (stats_mpi_datatype != MPI_DATATYPE_NULL) sct_abort("init_datatypes called twice", __FILE__, __LINE__);

#ifndef HAVE_LIBPAPI
  int count = 6;
#else
  int count = 11;
#endif
  int array_of_blocklengths[count];
  MPI_Aint array_of_displacements[count];
  MPI_Datatype array_of_types[count];

  sct_stats_type s; //scalar test variable
  sct_internal_alloc_event_arrays(&s, sct_get_event_num());
  MPI_Aint a0; // start address of structure
  MPI_Aint a; // component address

  int ib;// block pos

  // base address:
  if(MPI_Get_address(&s, &a0) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);

  //.tsum:
  ib = 0;
  if(MPI_Get_address(&s.tsum, &a) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
  array_of_blocklengths[ib] = 1;
  array_of_displacements[ib] = a-a0;
  array_of_types[ib] = MPI_DOUBLE;

  //.last_dt:
  ib = 1;
  if(MPI_Get_address(&s.last_dt, &a) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
  array_of_blocklengths[ib] = 1;
  array_of_displacements[ib] = a-a0;
  array_of_types[ib] = MPI_DOUBLE;

  //.tmin:
  ib = 2;
  if(MPI_Get_address(&s.tmin, &a) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
  array_of_blocklengths[ib] = 1;
  array_of_displacements[ib] = a-a0;
  array_of_types[ib] = MPI_DOUBLE;

  //.tmax:
  ib = 3;
  if(MPI_Get_address(&s.tmax, &a) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
  array_of_blocklengths[ib] = 1;
  array_of_displacements[ib] = a-a0;
  array_of_types[ib] = MPI_DOUBLE;

  //.cnum:
  ib = 4;
  if(MPI_Get_address(&s.cnum, &a) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
  array_of_blocklengths[ib] = 1;
  array_of_displacements[ib] = a-a0;
  array_of_types[ib] = MPI_INT;

  //.active_under:
  ib = 5;
  if(MPI_Get_address(&s.active_under, &a) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
  array_of_blocklengths[ib] = 1;
  array_of_displacements[ib] = a-a0;
  array_of_types[ib] = MPI_INT;

#ifdef HAVE_LIBPAPI
  const int max_evn = sct_get_event_num();
  MPI_Aint b0; // start address of one event vector
  MPI_Aint b; // component address of event vector entry
  MPI_Aint eventarray_bytelength;

  if (max_evn != 0) { /* no events are presents -> esum, emin, ... in sct_stats_type will be NULL */

    /* determine bytelength of eventarray containing MPI_LONG_LONG */
    if(MPI_Get_address(&(s.esum[0]), &b0) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
    if(MPI_Get_address(&(s.esum[max_evn]), &b) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
    eventarray_bytelength = b-b0;

    //.esum:
    ib = 6;
    if(MPI_Get_address(&s.esum, &a) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
    array_of_blocklengths[ib] = max_evn;
    /* starting address of first event-vector is valid since address was taken from s.esum pointer */
    array_of_displacements[ib] = a-a0;
    array_of_types[ib] = MPI_LONG_LONG;

    //.emin:
    ib = 7;
    array_of_blocklengths[ib] = max_evn;
    /* starting address of second event-vector cannot be determined by address of s.emin since this is only a pointer to array */
    array_of_displacements[ib] = array_of_displacements[ib-1] + eventarray_bytelength;
    array_of_types[ib] = MPI_LONG_LONG;

    //.emax:
    ib = 8;
    array_of_blocklengths[ib] = max_evn;
    array_of_displacements[ib] = array_of_displacements[ib-1] + eventarray_bytelength;
    array_of_types[ib] = MPI_LONG_LONG;

    //.rmin:
    ib = 9;
    array_of_blocklengths[ib] = max_evn;
    array_of_displacements[ib] = array_of_displacements[ib-1] + eventarray_bytelength;
    array_of_types[ib] = MPI_DOUBLE;

    /* determine bytelength of eventarray containing MPI_DOUBLE */
    if(MPI_Get_address(&(s.rmin[0]), &b0) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
    if(MPI_Get_address(&(s.rmin[max_evn]), &b) != MPI_SUCCESS) sct_abort("MPI_Get_address failed", __FILE__, __LINE__);
    eventarray_bytelength = b-b0;

    //.rmax:
    ib = 10;
    array_of_blocklengths[ib] = max_evn;
    array_of_displacements[ib] = array_of_displacements[ib-1] + eventarray_bytelength;
    array_of_types[ib] = MPI_DOUBLE;
  }
  else
    count = 6; /* exclude event-arrays */
#else
  const int max_evn = 0;
#endif

  /*for (int foo=0; foo<count; foo++)
    printf("array_of_displacements[%d] = %d\n",foo,array_of_displacements[foo]);
  */

  if (MPI_Type_create_struct(count,
                             array_of_blocklengths,
                             array_of_displacements,
                             array_of_types,
                             &stats_mpi_datatype)
      != MPI_SUCCESS) sct_abort("MPI_Type_create_struct failed", __FILE__, __LINE__);

  if (MPI_Type_commit(&stats_mpi_datatype) != MPI_SUCCESS) sct_abort("MPI_Type_commit failed", __FILE__, __LINE__);

  // test:

  sct_stats_type z; // zero stats
  sct_internal_alloc_event_arrays(&z, max_evn);
  sct_internal_zero_stats(&z, max_evn);
  s.tsum = 1.0;
  s.last_dt = 1.5;
  s.tmin = 2.0;
  s.tmax = 3.0;
  s.cnum = 4;
  s.active_under = -2;
#ifdef HAVE_LIBPAPI
  for(int ie = 0; ie < max_evn; ie++) {
    s.esum[ie] = 10000 + 10*ie + 1;
    s.emin[ie] = 10000 + 10*ie + 2;
    s.emax[ie] = 10000 + 10*ie + 3;
    s.rmin[ie] = 10000 + 10*ie + 4;
    s.rmax[ie] = 10000 + 10*ie + 5;
  }
#endif

  sct_stats_type v[3];
  for (int j = 0; j < 3; j++) sct_internal_alloc_event_arrays(&v[j], max_evn);
  sct_stats_type ref_v1;
  sct_internal_alloc_event_arrays(&ref_v1, max_evn);

  sct_internal_zero_stats(&ref_v1,max_evn);
  sct_internal_copy_stats(&s, &ref_v1, max_evn);

  for (int j = 0; j < 3; j++) {
    sct_internal_zero_stats(&v[j],max_evn);
  }

  MPI_Request request;

  /* pack sct_stats_type for MPI comm */
  int position;
  void *sendbuf;
  void *recvbuf;
  MPI_Type_size(stats_mpi_datatype, &stats_mpi_datatype_size);
  sendbuf = malloc(stats_mpi_datatype_size);
  recvbuf = malloc(stats_mpi_datatype_size);

  /* tsum, last_dt, tmin and tmax */
  position = 0;
  MPI_Pack(&(s.tsum), 4, MPI_DOUBLE, sendbuf, stats_mpi_datatype_size, &position, sct_comm_self);

  /* cnum and active_under */
  MPI_Pack(&(s.cnum), 2, MPI_INT, sendbuf, stats_mpi_datatype_size, &position, sct_comm_self);

#ifdef HAVE_LIBPAPI
  if (max_evn != 0) {
    MPI_Pack(&(s.esum[0]), max_evn, MPI_LONG_LONG, sendbuf, stats_mpi_datatype_size, &position, sct_comm_self);
    MPI_Pack(&(s.emin[0]), max_evn, MPI_LONG_LONG, sendbuf, stats_mpi_datatype_size, &position, sct_comm_self);
    MPI_Pack(&(s.emax[0]), max_evn, MPI_LONG_LONG, sendbuf, stats_mpi_datatype_size, &position, sct_comm_self);
    MPI_Pack(&(s.rmin[0]), max_evn, MPI_DOUBLE, sendbuf, stats_mpi_datatype_size, &position, sct_comm_self);
    MPI_Pack(&(s.rmax[0]), max_evn, MPI_DOUBLE, sendbuf, stats_mpi_datatype_size, &position, sct_comm_self);
  }
#endif

  if (MPI_Isend(sendbuf, position, MPI_PACKED, 0, 0, sct_comm_self, &request) != MPI_SUCCESS)
    sct_abort("MPI_Isend failed", __FILE__, __LINE__);

  if (MPI_Recv(recvbuf, stats_mpi_datatype_size, MPI_PACKED, 0, 0, sct_comm_self, MPI_STATUS_IGNORE) != MPI_SUCCESS)
    sct_abort("MPI_Recv failed", __FILE__, __LINE__);

  if (MPI_Request_free(&request) != MPI_SUCCESS)
    sct_abort("MPI_Request_free", __FILE__, __LINE__);

  /* unpack */
  position = 0;
  MPI_Unpack(recvbuf, stats_mpi_datatype_size, &position, &(v[1].tsum), 4, MPI_DOUBLE, sct_comm_self);
  MPI_Unpack(recvbuf, stats_mpi_datatype_size, &position, &(v[1].cnum), 2, MPI_INT, sct_comm_self);
#ifdef HAVE_LIBPAPI
  if (max_evn != 0) {
    MPI_Unpack(recvbuf, stats_mpi_datatype_size, &position, &(v[1].esum[0]), max_evn, MPI_LONG_LONG, sct_comm_self);
    MPI_Unpack(recvbuf, stats_mpi_datatype_size, &position, &(v[1].emin[0]), max_evn, MPI_LONG_LONG, sct_comm_self);
    MPI_Unpack(recvbuf, stats_mpi_datatype_size, &position, &(v[1].emax[0]), max_evn, MPI_LONG_LONG, sct_comm_self);
    MPI_Unpack(recvbuf, stats_mpi_datatype_size, &position, &(v[1].rmin[0]), max_evn, MPI_DOUBLE, sct_comm_self);
    MPI_Unpack(recvbuf, stats_mpi_datatype_size, &position, &(v[1].rmax[0]), max_evn, MPI_DOUBLE, sct_comm_self);
  }
#endif

  //check:
  if ( ! stats_are_equal(&v[0], &z, max_evn) )
    sct_abort("stats datatype does not work as expected", __FILE__, __LINE__);
  if ( ! stats_are_equal(&v[1], &ref_v1, max_evn) ) {
    printf("ref_v1.tsum = %f\n",ref_v1.tsum);
    printf("ref_v1.last_dt = %f\n",ref_v1.last_dt);
    printf("ref_v1.tmin = %f\n",ref_v1.tmin);
    printf("ref_v1.tmax = %f\n",ref_v1.tmax);
    printf("ref_v1.tcnum = %d\n",ref_v1.cnum);
    printf("ref_v1.active_under = %d\n",ref_v1.active_under);
#ifdef HAVE_LIBPAPI
    for(int ie = 0; ie < max_evn; ie++) {
      printf("ref_v1.esum[i] = %lld\n", ref_v1.esum[ie]);
      printf("ref_v1.emin[i] = %lld\n", ref_v1.emin[ie]);
      printf("ref_v1.emax[i] = %lld\n", ref_v1.emax[ie]);
      printf("ref_v1.rmin[i] = %f\n", ref_v1.rmin[ie]);
      printf("ref_v1.rmax[i] = %f\n", ref_v1.rmax[ie]);
    }
#endif

    printf("\nv1.tsum = %f\n",v[1].tsum);
    printf("v1.last_dt = %f\n",v[1].last_dt);
    printf("v1.tmin = %f\n", v[1].tmin);
    printf("v1.tmax = %f\n", v[1].tmax);
    printf("v1.tcnum = %d\n",v[1].cnum);
    printf("v1.active_under = %d\n",v[1].active_under);
#ifdef HAVE_LIBPAPI
    for(int ie = 0; ie < max_evn; ie++) {
      printf("v1.esum[i] = %lld\n", v[1].esum[ie]);
      printf("v1.emin[i] = %lld\n", v[1].emin[ie]);
      printf("v1.emax[i] = %lld\n", v[1].emax[ie]);
      printf("v1.rmin[i] = %f\n", v[1].rmin[ie]);
      printf("v1.rmax[i] = %f\n", v[1].rmax[ie]);
    }
#endif

    sct_abort("stats datatype does not work as expected for MPI-comm", __FILE__, __LINE__);
  }
  if ( ! stats_are_equal(&v[2], &z, max_evn) )
    sct_abort("stats datatype does not work as expected", __FILE__, __LINE__);

  sct_internal_free_event_arrays(&s);
  sct_internal_free_event_arrays(&z);
  for (int j = 0; j < 3; j++) sct_internal_free_event_arrays(&v[j]);
  sct_internal_free_event_arrays(&ref_v1);

  free(sendbuf);
  free(recvbuf);

#ifdef DEBUG
  fprintf(stderr, "(debug) stats datatype passed all tests\n");
#endif
}
#endif


int sct_init(const int tsize, const char* default_context_name,
             const sct_comm_type default_comm) {

  // initialization of timer memory
  // returns the number of allowed timers
  if (timer_size) {
    // second call:
    // for now, we only allow compatible second calls
    fprintf(stderr,"# SCT-Warning: sct_init called again: call will be ignored\n");
    return timer_size;
  }

  /* check if callstats should be recorded */
  callstats = sct_get_callstats();

  /* set defaults if given values are nonsense */
  if (callstats < 0  || callstats > 1) {
#if (defined(HAVE_MPI) || defined(_OPENMP))
    callstats = 0;
#else
    callstats = 1;
#endif
#ifdef DEBUG
    fprintf(stderr,"# SCT-Warning: Set callstats to fallback value (%i).\n",callstats);
#endif
  }

  /* check if eventcounter should be recorded */
  eventcounters = sct_get_eventcounters();

  /* set defaults if given values are nonsense */
  if (eventcounters < 0  || eventcounters > 2) {
#ifdef HAVE_LIBPAPI
    eventcounters = 2; /* show eventrates */
#else
    eventcounters = 0;
#endif
#ifdef DEBUG
    fprintf(stderr,"# SCT-Warning: Set eventcounters to fallback value (%i).\n",eventcounters);
#endif
  }

#ifdef HAVE_LIBPAPI
  /* define PAPI eventset if SCT_EVENTCOUNTERS are needed */
  if (eventcounters) {
    if (!(eset = sct_new_eventset())) sct_abort("sct_new_eventset failed", __FILE__, __LINE__);
    if (!sct_start_eventset(eset)) sct_abort("sct_start_eventset failed", __FILE__, __LINE__);
  }
  else
    eset = NULL;
#else
  eset = NULL;
#endif

#ifdef HAVE_MPI
  if (MPI_Comm_dup(MPI_COMM_SELF, &sct_comm_self) != MPI_SUCCESS)
    sct_abort("MPI_Comm_dup failed", __FILE__, __LINE__);
  init_datatypes();
#endif

  /* check if nestedtimers should be used */
  nestedtimers = sct_get_nestedtimers();

  /* set defaults if given values are nonsense */
  if (nestedtimers < 0  || nestedtimers > 1) {
    nestedtimers = 0;
#ifdef DEBUG
    fprintf(stderr,"# SCT-Warning: Set nestedtimers to fallback value (%i).\n", nestedtimers);
#endif
  }
#ifndef NESTED_TIMER
  if (nestedtimers == 1) sct_abort("sct_init: SCT_NESTEDTIMERS cannot be set if compiled without NESTED_TIMER macro set", __FILE__, __LINE__);
#endif

  // first call:
  sct_mach_init();

  // set current timer_size:
  timer_size = tsize;
#ifdef _OPENMP
    // enforce some thread-data separation:
    if (timer_size<16) timer_size = 16;
#endif

  set_prg_start_time();

  // new meta data:
  if ( !(meta_pool = calloc(timer_size, sizeof(sct_meta_type))) )
    sct_abort("memory allocation failed", __FILE__, __LINE__);

  // thread_num:
#ifdef _OPENMP
  static_pr_thread_count = omp_get_max_threads();
#else
  static_pr_thread_count = 0;  // no parallel region
#endif

  // thread private data:
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int tid = omp_get_thread_num();
    static_pr_thread_id = tid;
    if ( !(mark_pool = calloc(timer_size, sizeof(sct_mark_type))) )
      sct_abort("memory allocation failed", __FILE__, __LINE__);
#ifdef HAVE_LIBPAPI
    if (eventcounters) {
      for (int i=0; i<timer_size; i++) {
        if ( !(mark_pool[i].em = calloc(eset->en,sizeof(eval_type))) )
          sct_abort("memory allocation failed", __FILE__, __LINE__);
      }
    }
#endif
  }

  if (context_num == 0) {
    // default context:
    if (default_context_name)
      default_icontext = sct_new_context(default_context_name, default_comm);
    else
      default_icontext = sct_new_context("default context", default_comm);
    if (default_icontext != 0) sct_abort("unexpected handle for default context", __FILE__, __LINE__);
    icontext = default_icontext;
  }

  return timer_size;
}


void sct_finalize() {
#ifdef HAVE_LIBPAPI

  sct_stop_eventset(eset);

  if (eventcounters) {

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int tid = omp_get_thread_num();
    int ierr;

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      if ( (ierr=PAPI_cleanup_eventset(eset->id[tid])) != PAPI_OK)
        fprintf(stderr, "PAPI ERROR in file %s line %d: %s failed (%s)\n",
	        __FILE__, __LINE__, "PAPI_cleanup_eventset", PAPI_strerror(ierr));
      if ( (ierr=PAPI_destroy_eventset(&eset->id[tid])) != PAPI_OK)
        fprintf(stderr, "PAPI ERROR in file %s line %d: %s failed (%s)\n",
	        __FILE__, __LINE__, "PAPI_destroy_eventset", PAPI_strerror(ierr));

      /// \todo Fix error "EventSet is currently counting"

    }
  }

    PAPI_shutdown();
  }
#endif

  // cleanup any attributes
  sct_internal_free_attribute_table();
}


void sct_reset_timer(const int it) {

#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer out of bounds", __FILE__, __LINE__);
#endif

#ifdef _OPENMP
  if (omp_in_parallel())
    sct_abort("cannot reset timer inside parallel region", __FILE__, __LINE__);
#endif

  for(int icon = 0; icon < context_num; icon++) {
    sct_context_type *con = &context_pool[icon];

    con->s[it].tsum  = 0.0;
    con->s[it].active_under = -2;
    con->s[it].cnum = 0;
    con->s[it].tmin = SCT_MAX_TVAL;
    con->s[it].tmax = 0.0;
#ifdef HAVE_LIBPAPI
    for (int ie=0; ie<sct_get_event_num(); ie++) {
      con->s[it].esum[ie] = 0;
      con->s[it].emin[ie] = SCT_MAX_EVAL;
      con->s[it].emax[ie] = 0;
      con->s[it].rmin[ie] = SCT_MAX_EVAL;
      con->s[it].rmax[ie] = 0.0;
    }
#endif


#ifdef _OPENMP

#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      con->p[tid][it].tsum  = 0.0;
      con->p[tid][it].active_under = -2;
      con->p[tid][it].cnum = 0;
      con->p[tid][it].tmin = 1.0e+20;
      con->p[tid][it].tmax = 0.0;
#ifdef HAVE_LIBPAPI
      for (int ie=0; ie<sct_get_event_num(); ie++) {
        con->p[tid][it].esum[ie] = 0;
        con->p[tid][it].emin[ie] = SCT_MAX_EVAL;
        con->p[tid][it].emax[ie] = 0;
        con->p[tid][it].rmin[ie] = SCT_MAX_EVAL;
        con->p[tid][it].rmax[ie] = 0.0;
      }
#endif
    }
#endif

    //fprintf(stderr, "(reset) icon=%i, it=%i, cnum=%i, .tmin=%e\n",icon, it,  context_pool[icon].p[0][it].cnum, context_pool[icon].p[0][it].tmin);
  }
}


void sct_reset_all() {
  for (int it=0; it<timer_num; it++) {
    sct_reset_timer(it);
  }
}


int sct_new_timer(const char* name) {
  int it;

#ifdef _OPENMP
  if (omp_in_parallel())
    sct_abort("ERROR: Cannot create new timer in a parallel region", __FILE__, __LINE__);
#endif

  if (!context_num) sct_init(SCT_DEFAULT_TIMER_SIZE, NULL, SCT_COMM_WORLD);

#ifdef DEBUG
  fprintf(stderr,"sct_new_timer: name=%s, it=%i\n",name,timer_num);
#endif

  if (timer_num == timer_size) {
    // we run out of pre-allocated resources and need to realloc our pools
    sct_warn("inititial timer_size too small - need to realloc resources", __FILE__, __LINE__);

    // enlarge meta data pool
    if ( !(meta_pool = realloc(meta_pool, 2*timer_size*sizeof(sct_meta_type))) )
      sct_abort("memory reallocation failed", __FILE__, __LINE__);

    // enlarge mark pool
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      int tid = omp_get_thread_num();
      static_pr_thread_id = tid;
      if ( !(mark_pool = realloc(mark_pool, 2*timer_size*sizeof(sct_mark_type))) )
	sct_abort("memory reallocation failed", __FILE__, __LINE__);
#ifdef HAVE_LIBPAPI
      for (int i=timer_size; i<2*timer_size; i++) {
	if ( !(mark_pool[i].em = calloc(sct_get_event_num(), sizeof(eval_type))) )
	  sct_abort("memory allocation failed", __FILE__, __LINE__);
      }
#endif
    }

    // enlarge contexts
    for(int icon = 0; icon < context_num; icon++) {
      sct_context_type *con = &context_pool[icon];

      // serial measurement data:
      if ( !(con->s = realloc(con->s, 2*timer_size*sizeof(sct_stats_type))) )
	sct_abort("memory reallocation failed", __FILE__, __LINE__);

      for (int i=timer_size; i<2*timer_size; i++) {
	con->s[i].tmin = SCT_MAX_TVAL;
	con->s[i].active_under = -2;
#ifdef HAVE_LIBPAPI
	// allocate event arrays
	sct_internal_alloc_event_arrays(&(con->s[i]), sct_get_event_num());
	// init arrays of min values
	for (int j=0; j<sct_get_event_num(); j++) {
	  con->s[i].emin[j] = SCT_MAX_EVAL;
	  con->s[i].rmin[j] = SCT_MAX_TVAL;
	}
#endif
      }

#ifdef _OPENMP
      // threadprivate space, minor dimension:
#pragma omp parallel
      {
	int tid;
	if ( (tid = sct_get_pr_thread_id()) >= sct_get_pr_thread_count())
	  sct_abort("internal error", __FILE__, __LINE__);

#pragma omp critical
	{
	  if ( !(con->p[tid] = realloc(con->p[tid], 2*timer_size*sizeof(sct_stats_type))) )
	    sct_abort("memory reallocation failed", __FILE__, __LINE__);

	  for (int i=timer_size; i<2*timer_size; i++) {
	    con->p[tid][i].tmin = SCT_MAX_TVAL;
	    con->p[tid][i].active_under = -2;
#ifdef HAVE_LIBPAPI
	    // allocate event arrays
	    sct_internal_alloc_event_arrays(&(con->p[tid][i]), sct_get_event_num());
	    // init arrays of min values
	    for (int j=0; j<sct_get_event_num(); j++) {
	      con->p[tid][i].emin[j] = SCT_MAX_EVAL;
	      con->p[tid][i].rmin[j] = SCT_MAX_TVAL;
	    }
#endif
	  }
	}
      }
#endif
    }

    // set new timer_size
    timer_size *= 2;
  }

  // new timer
  it = timer_num;
  timer_num++;

  // meta
  meta_pool[it].name = sct_string_new(name);
  meta_pool[it].used = 1;

  // mark
  mark_pool[it].state = sct_off_state;

  sct_reset_timer(it);

  return it;
}

void sct_del_timer(const int it) {

#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer out of bounds", __FILE__, __LINE__);
#endif

  meta_pool[it].used = 0;
  sct_string_delete(&meta_pool[it].name);

  while(timer_num>0) {
    if (meta_pool[timer_num-1].used) break;
    timer_num--;
  }
}

int sct_active(const int it) {
#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer out of bounds", __FILE__, __LINE__);
  if (icontext<0) abort_timer(it, "missing context", __FILE__, __LINE__);
#endif
  if (mark_pool[it].state == sct_on_state)
    return 1;
  else
    return 0;
}

void sct_start(const int it) {

  sct_context_type *con = &context_pool[icontext];

#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer out of bounds", __FILE__, __LINE__);
  if (mark_pool[it].state == sct_on_state) abort_timer(it, "sct_stop call missing", __FILE__, __LINE__);
  if (icontext<0) abort_timer(it, "missing context", __FILE__, __LINE__);
  mark_pool[it].icon = icontext;
#endif

#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  const int tid = 0;
#endif

#ifdef _OPENMP
  if (omp_in_parallel()) {
#ifdef DEBUG
    fprintf(stderr,"sct_start: tid=%d, stat = p\n",tid);
#endif
#ifdef NESTED_TIMER
    if (nestedtimers) update_timer_stack_atstart(&con->active_timer_p[tid], con->p[tid], it);
#endif
  }
  else {
#ifdef DEBUG
    fprintf(stderr,"sct_start: tid=%d, stat = s\n",tid);
#endif
#ifdef NESTED_TIMER
    if (nestedtimers) update_timer_stack_atstart(&con->active_timer, con->s, it);
#endif
  }
#else
#ifdef NESTED_TIMER
  if (nestedtimers) update_timer_stack_atstart(&con->active_timer, con->s, it);
#endif
#endif

  /* read timer first, then event counter */

  read_time(&mark_pool[it].tm);

#ifdef HAVE_LIBPAPI
  if (eventcounters) {
    if (!sct_read_events(tid, eset, mark_pool[it].em)) sct_abort("read_events failed", __FILE__, __LINE__);
  }
#endif

  mark_pool[it].state = sct_on_state;
}

#ifdef HAVE_LIBPAPI

static void inline update_stats(sct_stats_type *restrict stats, double dt, eval_type *restrict de, int en) {

  stats->tsum += dt;
  stats->last_dt = dt;
  stats->cnum++;
  if (callstats) {
    double r;
    if (dt < stats->tmin) stats->tmin = dt;
    if (dt > stats->tmax) stats->tmax = dt;
    for (int ie=0; ie<en; ie++) {
      r = de[ie]/dt;
      stats->esum[ie] += de[ie];
      if (de[ie] < stats->emin[ie]) stats->emin[ie] = de[ie];
      if (de[ie] > stats->emax[ie]) stats->emax[ie] = de[ie];
      if (r < stats->rmin[ie]) stats->rmin[ie] = r;
      if (r > stats->rmax[ie]) stats->rmax[ie] = r;
    }
  } else {
    for (int ie=0; ie<en; ie++) {
      stats->esum[ie] += de[ie];
    }
  }
}

#else

static void inline update_stats(sct_stats_type *restrict stats, double dt) {

  stats->tsum += dt;
  stats->last_dt = dt;
  stats->cnum++;
  if (callstats) {
    if (dt < stats->tmin) stats->tmin = dt;
    if (dt > stats->tmax) stats->tmax = dt;
  }
}

#endif


void sct_stop(const int it) {

  sct_context_type *con = &context_pool[icontext];

#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("sct_stop: timer out of bounds", __FILE__, __LINE__);
  if (mark_pool[it].state != sct_on_state) {
    if (mark_pool[it].state == sct_off_state) {
      abort_timer(it, "sct_stop: sct_start call missing", __FILE__, __LINE__);
    } else {
      abort_timer(it, "sct_stop: undefined timer", __FILE__, __LINE__);
    }
  }
  if (icontext<0 || icontext>=context_num) abort_timer(it, "invalid context", __FILE__, __LINE__);
  if (icontext != mark_pool[it].icon) abort_timer(it, "mixed contexts", __FILE__, __LINE__);
#endif

#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  const int tid = 0;
#endif

  /* read event counter first, then timer */

  tmark_type tm;
#ifdef HAVE_LIBPAPI
  const int en = (eventcounters) ? eset->en : 0;
  eval_type em[en];
  if (eventcounters) {
    if (!sct_read_events(tid, eset, em)) sct_abort("read_events failed", __FILE__, __LINE__);
  }
#endif
  read_time(&tm);
  const double dt = get_tdiff(&mark_pool[it].tm, &tm);

#ifdef HAVE_LIBPAPI
  eval_type de[en];
  if (eventcounters) get_ediff(en, mark_pool[it].em, em, de);
#endif

#ifdef _OPENMP
  if (omp_in_parallel()) {
    // inside parallel region: add differences to thread private counter
#ifdef HAVE_LIBPAPI
    update_stats(&con->p[tid][it], dt, de, en);
#else
    update_stats(&con->p[tid][it], dt);
#endif
#ifdef NESTED_TIMER
    if (nestedtimers) update_timer_stack_atstop(&con->active_timer_p[tid], it);
#endif
  } else {
    // outside parallel region: add differences to serial phase counter
#ifdef HAVE_LIBPAPI
    update_stats(&con->s[it], dt, de, en);
#else
    update_stats(&con->s[it], dt);
#endif
#ifdef NESTED_TIMER
    if (nestedtimers) update_timer_stack_atstop(&con->active_timer, it);
#endif
  }
  // end of _OPENMP branch

#else

  // no OpenMP:
#ifdef HAVE_LIBPAPI
  update_stats(&con->s[it], dt, de, en);
#else
  update_stats(&con->s[it], dt);
#endif
#ifdef NESTED_TIMER
  if (nestedtimers) update_timer_stack_atstop(&con->active_timer, it);
#endif
  // end of undef(_OPENMP)

#endif

  mark_pool[it].state = sct_off_state;
}


double sct_val(const int it) {
  tmark_type tm;
  double t0;

#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer_val: timer out of bounds", __FILE__, __LINE__);
  if (icontext<0) abort_timer(it, "missing context", __FILE__, __LINE__);
#endif

  sct_context_type *con = &context_pool[icontext];

#ifdef _OPENMP
  if (omp_in_parallel()) {
    int tid = omp_get_thread_num();
    t0 = con->p[tid][it].tsum;
  } else {
    t0 = con->s[it].tsum;
  }
#else
  t0 = con->s[it].tsum;
#endif

  // add current delta time if timer is running:
  if (mark_pool[it].state == sct_on_state) {
    read_time(&tm);
    return t0 + get_tdiff(&mark_pool[it].tm, &tm);
  } else {
    return t0;
  }

}


double sct_last_dt(const int it) {
  tmark_type tm;
  double t0;

#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer_val: timer out of bounds", __FILE__, __LINE__);
  if (icontext<0) abort_timer(it, "missing context", __FILE__, __LINE__);
#endif

  sct_context_type *con = &context_pool[icontext];

#ifdef _OPENMP
  if (omp_in_parallel()) {
    int tid = omp_get_thread_num();
    return con->p[tid][it].last_dt;
  } else {
    return con->s[it].last_dt;
  }
#else
  return con->s[it].last_dt;
#endif

}


double sct_event(const int it, const char *ev) {
#ifdef HAVE_LIBPAPI

  int ie, err = 1;
#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer_val: timer out of bounds", __FILE__, __LINE__);
  if (icontext<0) abort_timer(it, "missing context", __FILE__, __LINE__);
#endif

  sct_context_type *con = &context_pool[icontext];

  for (ie=0; ie<eset->en; ie++) {
    if (!strcmp(eset->name[ie], ev)) {
      err = 0;
      break;
    }
  }
  if (err) {
    sct_warn("requesting value of unknown event - return 0 instead as fallback", __FILE__, __LINE__);
    return 0.0;
  }

  // return value based on type of event counting
  if (sct_get_eventcounters() == 1) {
#ifdef _OPENMP
    if (omp_in_parallel()) {
      int tid = omp_get_thread_num();
      return (double) con->p[tid][it].esum[ie];
    } else {
      return (double) con->s[it].esum[ie];
    }
#else
    return (double) con->s[it].esum[ie];
#endif
  }
  else if (sct_get_eventcounters() == 2) {
#ifdef _OPENMP
    if (omp_in_parallel()) {
      int tid = omp_get_thread_num();
      return (double) con->p[tid][it].esum[ie] / con->p[tid][it].tsum;
    } else {
      return (double) con->s[it].esum[ie] / con->s[it].tsum;
    }
#else
    return (double) con->s[it].esum[ie] / con->s[it].tsum;
#endif
  }
  else
    return 0.0;

#else
  return 0.0;
#endif
}


void sct_stop_all() {
  for (int it = 0; it<timer_num; it++) {
    if (meta_pool[it].used) {
      if (mark_pool[it].state == sct_on_state) {
        fprintf(stderr,"# SCT-Warning: closing timer: id=%i, [label]= [%s]\n",
                it, meta_pool[it].name.cs);
        sct_stop(it);
      }
    }
  }
}

void sct_try_stop(const int it) {
  if (meta_pool[it].used) {
    if (mark_pool[it].state == sct_on_state) {
      fprintf(stderr,"# SCT-Warning: closing timer: id=%i, [label]= [%s]\n",
              it, meta_pool[it].name.cs);
      sct_stop(it);
    }
  }
}

double sct_resolution() {
  int n;
  double dt1, dt2;
  tmark_type tmark1, tmark2;

#ifdef _OPENMP
  if (omp_in_parallel())
    sct_abort("cannot measure resolution inside parallel region", __FILE__, __LINE__);
#endif

  // estimate upper bound for internal resolution:
  n=0;
  read_time(&tmark1);
  read_time(&tmark2);
  while( equal_tmarks(&tmark1, &tmark2) ) {
    read_time(&tmark2);
    n++;
    if (n>1000000) sct_abort("sct_resolution: I give up.", __FILE__, __LINE__);
  }
  dt1=get_tdiff(&tmark1, &tmark2);

  /*

  // we don't need to check the call overhead here

  int it;

  if (!context_num) sct_init(SCT_DEFAULT_TIMER_SIZE);

  it = sct_new_timer("test");

#ifdef CHECK_TIMER
  if (it<0 || it>=timer_num) sct_abort("timer_start: timer out of bounds", __FILE__, __LINE__);
#endif

  // estimate call overhead:
  for (int i=0; i<10; i++) {
    sct_start(it);
    sct_stop(it);
  }

  sct_context_type *con = &context_pool[icontext];
x
  dt2 = con->p[0][it].tsum * 0.1;
  if (dt2 <= dt1) fprintf(stderr,"sct_resolution: call overhead <= resolution (%e <= %e)\n",dt2,dt1);

  sct_del_timer(it);
  return dmax(dt1,dt2);

  */

  return dt1;
}

static void abort_timer(const int it, const char *reason, const char *fname, const int line) {
  if (it < 0 || it>=timer_num) {
    fprintf(stderr, "ERROR in file %s, line %d: %s\n",fname, line, reason);
  } else {
    fprintf(stderr, "ERROR in [file]=[%s], [line]=[%d], [label]=[%s]: [reason]=[%s]\n",fname, line, meta_pool[it].name.cs, reason);
  }
  fflush(stdout);
  fflush(stderr);
  sct_abort(NULL, NULL, 0);
}

// internal low level inline functions:

static inline double dmax(double a, double b) {
  return a > b ? a : b;
}


// internal function to update the simple calltree of timers
#ifdef NESTED_TIMER
static void update_timer_stack_atstart(sct_stack_type *s, sct_stats_type *stats, const int it) {

  sct_stats_type *stat = &stats[it];

#ifdef DEBUG
  printf("DEBUG start timer_stack:\n");
  for (int i=0; i<s->n; i++)
    printf("%d ", s->stack[i]);
  printf("pos -> %d\n",s->pos);
#endif

  /*
   * Find the "nearest" common node of all callchain paths
   * from the current stat to the root.
   */
  if (stat->active_under == -2) { /* No ancestor set yet */
    stat->active_under = s->stack[s->pos];
  } else if(stat->active_under != -1) {
    /*
     * Determine a new common ancestor by simultaneously searching the new
     * chain as given by the current stack and the old chain of common
     * ancestors, the linked list rooted at stat->active_under, for a new
     * common ancestor.
     *
     * Note that all elements on the stack have underwent that procedure
     * already and thus, the following invariants hold:
     * i)  For every element on the stack, every list rooted at that element's
     *     ->active_under is contained in the stack (below that element of
     *     course), because it is only reachable through its ->active_under.
     * ii) Let i and j be two positions on the stack with i below j.
     *     Then either the position of j->active_under is above or equal to that
     *     of i (i is some common ancestor of j) or the position of
     *     j->active_under is below or equal to that of i->active_under.
     * These remarks allow us to not to examine the whole stack in our search
     * for common ancestors, but only the linked ->active_under list rooted at
     * its top.
     */
    for (int it_new_sup = s->stack[s->pos]; it_new_sup >= 0;
      it_new_sup = stats[it_new_sup].active_under) {
      for (int it_old_sup = stat->active_under; it_old_sup >= 0;
        it_old_sup = stats[it_old_sup].active_under) {
        if (it_old_sup == it_new_sup) {
          stat->active_under = it_new_sup;
          goto found_common_ancestor;
        }
      }
    }
    /* no common ancestor found */
    stat->active_under = -1;
  }
found_common_ancestor:

  /* increase stack pointer */
  (s->pos)++;

  /* realloc stack if needed */
  if (s->pos >= s->n) {
    sct_warn("number of simultaneously active timers caused realloc of stack", __FILE__, __LINE__);

    if ( !(s->stack = realloc(s->stack, 2 * s->n * sizeof(int))) )
      sct_abort("memory reallocation failed", __FILE__, __LINE__);
    s->n *= 2;

#ifdef _OPENMP
    if (omp_in_parallel()) {
#pragma omp critical
      {
        timer_nest_depth = MAX(timer_nest_depth, s->n);
      }
    }
    else
#endif
      timer_nest_depth = MAX(timer_nest_depth, s->n);
  }

  /* push actual timer on stack */
  s->stack[s->pos] = it;
}


static void update_timer_stack_atstop(sct_stack_type *s, const int it) {

#ifdef DEBUG
  printf("DEBUG stop timer_stack\n");
  for (int i=0; i<s->n; i++)
    printf("%d ", s->stack[i]);
  printf("pos -> %d\n",s->pos);
#endif

  if (s->stack[s->pos] != it)
    abort_timer(it,"timer_stop: a subsidary timer is still active, stop that first", __FILE__, __LINE__);

  /* pop stack, ie. decrease stack pointer */
  (s->pos)--;
}
#endif

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
 * sct_collector.h: measure and collect time and events
 */

#ifndef _H_SCT_COLLECTOR
#define _H_SCT_COLLECTOR

#include "sct_mach.h"

//public+
#ifdef HAVE_MPI
/*! \brief null communicator in case that MPI is used */
#  define SCT_COMM_NULL MPI_COMM_NULL
/*! \brief comm world in case that MPI is used */
#  define SCT_COMM_WORLD MPI_COMM_WORLD
typedef MPI_Comm sct_comm_type;
typedef MPI_Datatype sct_datatype;
#else
/*! \brief null communicator that is also valid without MPI */
#  define SCT_COMM_NULL 0
/*! \brief comm world that is also valid without MPI */
#  define SCT_COMM_WORLD 0
/*! \brief use integer as comm_type if no MPI is used */
typedef int sct_comm_type;
/*! \brief use integer as datatype if no MPI is used */
typedef int sct_datatype;
#endif

#define SCT_GETENV -255
#define SCT_WITHOUT_CALLSTATS 0
#define SCT_WITH_CALLSTATS 1
#define SCT_DEFAULT_CALLSTATS 2
//public-
#define SCT_MAX_TVAL HUGE_VAL


/*! \brief special string implementation */
typedef struct {
  int n;    //!< used length without terminating zero (fortran style)
  int cap;  //!< character capacity
  char *cs; //!< null-terminated cstring of size cap+1
} sct_string_type;

#define sct_string_zero (sct_string_type){.n=0, .cap=0, .cs=NULL}


/*! \brief shared meta data type per timer */
typedef struct {
  int used;             //!< 0: unused, 1: used
  sct_string_type name; //!< label for report
} sct_meta_type;


/*! \brief type for adding time and events */
typedef struct {
  double tsum;     //!< sum over all measured timeframes
  double last_dt;  //!< duration of last timeframes
  double tmin;     //!< minimum timeframe duration
  double tmax;     //!< maximum timeframe duration
  int cnum;        //!< number of timer calls
  int active_under;//!< timerID of superordinate timer; -2: "undefined", -1: "none, root-timer", >=0: id
#ifdef HAVE_LIBPAPI
  eval_type *esum; //!< event sum over all measured timeframes
  eval_type *emin; //!< minimum event measure ever measured
  eval_type *emax; //!< maximum event measure ever measured
  double *rmin;    //!< minimum event-rate measure ever measured
  double *rmax;    //!< maximum event-rate measure ever measured
#endif
} sct_stats_type;


#ifdef NESTED_TIMER
/*! \brief dynamic stack to remind timer nesting */
typedef struct {
  int n;       //!< length of dynamic stack
  int *stack;  //!< data stack
  int pos;     //!< stack top pointer
} sct_stack_type;
#endif


/*! \brief holds all members needed to specify a context */
typedef struct {
  int valid;                    //!< 0: invalid; 1: valid
  int active;                   //!< 1/0 = yes/no
  sct_string_type name;         //!< context name
#ifdef NESTED_TIMER
  sct_stack_type active_timer;    //!< serial stack
#ifdef _OPENMP
  sct_stack_type *active_timer_p; //!< threadprivate stack
#endif
#endif
#ifdef HAVE_MPI
  MPI_Comm comm;                //!< communicator
#endif
  int procnum;                  //!< number of mpi tasks
  int pid;                      //!< MPI process id
#ifdef _OPENMP
  /*! data for thread-parallel (private) part of timer data
      each thread gets its own instance
      used for measurements within a parallel region */
  sct_stats_type **p;             //!< p[ithread][itimer]
#endif
  /*! serial-phase part of timer data (shared) outside of parallel regions */
  sct_stats_type *s;              //!< s[itimer]
} sct_context_type;


/*! \brief mark type to hold all data of one single context free measurement */
typedef struct {
  tmark_type tm;                  //!< time mark
#ifdef HAVE_LIBPAPI
  eval_type *em;                  //!< event mark
#endif
  int state;                      //!< e.g., on, off, undef
#ifdef CHECK_TIMER
  /*! allow us to verify that start and stop belong to the same context */
  int icon;                       //!< context handle
#endif
} sct_mark_type;


typedef enum {
  SCT_INT = 0,
  SCT_LONG,
  SCT_FLOAT,
  SCT_DOUBLE,
  SCT_STRING,
  SCT_VALUE_UNDEFINED
} sct_attribute_value_type;

/*! \brief attribute type to hold a generic key-value pair to be used for reporting */
typedef struct sct_attribute {
  char* key;
  void* value;
  sct_attribute_value_type type;
  struct sct_attribute *next;
} sct_attribute_type;

/*! \brief attribute table implemented as linked list of attributes*/
typedef struct {
  sct_attribute_type *first;
} sct_attribute_table;


// declarations:

//public+
/*! \brief initialisation function of sct
    \param[in] tsize number of timers to be used
    \param[in] default_context_name name of the context to be used
    \param[in] default_comm communicator used for this context
    \return timer_size; number of timers that can be used
    \details This function needs to be called in the very beginning before any timer can be used. It sets the whole needed environment.
 */
int sct_init(const int tsize, const char* default_context_name, const sct_comm_type default_comm);

/*! \brief finalising function of sct
    \details This function should be called in the very end of your program or whenever you do not want
             to use timers anymore. It mainly frees dynamic memory, such that no call to sct is avail afterwards!
 */
void sct_finalize();

/*! \brief defines a new local context
    \param[in] name name of the new context
    \param[in] comm communicator to be ised for reduction within this context
    \return handle for a new context; handle is always > 0 since 0 is used for the default_context
    \details a new local context is created using the same number and naming for timers as the default context
*/
int sct_new_context(const char* name, sct_comm_type comm);

/*! \brief defines a new global context, ie. sct_new_context using MPI_COMM_WORLD or sct_comm_self
    \param[in] name name of the new context
    \return handle for a new context; handle is always > 0 since 0 is used for the default_context
    \details a new global context is created using the same number and naming for timers as the default context
*/
int sct_new_global_context(const char* name);

/*! \brief initialise a context switch
    \param[in] icon handle of context to be used
    \details from default context to context icon
*/
void sct_context_start(const int icon);

/*! \brief finalise a context switch
    \param[in] icon handle of context to stop
    \details from context icon back to default context
*/
void sct_context_stop(const int icon);

/*! \brief initialise a new timer
    \param[in] name name of timer to be used for report
    \return handle of new timer
*/
int sct_new_timer(const char* name);

/*! \brief start timer if not active
    \param[in] it handle of timer to be started
*/
void sct_start(const int it);

/*! \brief stop timer if active - report error else
    \param[in] it handle of timer to be stopped
*/
void sct_stop(const int it);

/*! \brief stop all active timer
*/
void sct_stop_all();

/*! \brief get actual value of timer (accumulated time of all calls)
    \param[in] it handle of timer to be evaluated
    \return actual timer value
*/
double sct_val(const int it);

/*! \brief get value of last timed interval
    \param[in] it handle of timer to be evaluated
    \return time of last timed interval
*/
double sct_last_dt(const int it);

/*! \brief get actual value of event (accumulated over all calls so far)
    \param[in] it handle of timer to be evaluated
    \param[in] ev name of event to be evaluated
    \return actual event value
*/
double sct_event(const int it, const char* ev);

/*! \brief get timer resolution of underlying system
    \return timer resolution
*/
double sct_resolution();

/*! \brief reset all statistics of timer
    \param[in] it handle of timer to be reset
*/
void sct_reset_timer(const int it);

/*! \brief reset all statistics of all timer
*/
void sct_reset_all();

/*! \brief check whether timer is in use
    \param[in] it handle of timer to be checked
    \return 0: timer not active, 1: timer already started
*/
int sct_active(const int it);

/*! \brief delete timer if not active
    \param[in] it handle of timer to be deleted
*/
void sct_del_timer(const int it);

/*! \brief set output mode for callstatistics
    \param[in] val 0 = disable output, 1 = enable output
*/
void sct_set_callstats(const int val);

/*! \brief set output mode for eventcounters
    \param[in] val 0 = disable output, 1 = enable output
*/
void sct_set_eventcounters(const int val);

/*! \brief set output mode for nested timers
    \param[in] val 0 = disable output, 1 = enable output
*/
void sct_set_nestedtimers(const int val);

/*! \brief add attribute for timer report in key-value pair style
    \param[in] key
    \param[in] val
 */
void sct_add_report_attribute_int(const char *key, int val);
void sct_add_report_attribute_long(const char *key, long val);
void sct_add_report_attribute_float(const char *key, float val);
void sct_add_report_attribute_double(const char *key, double val);
void sct_add_report_attribute_string(const char *key, char *val);
#ifdef HAVE_C__GENERIC
#define sct_add_report_attribute(_1, _2) _Generic((_2),                              \
                                           int: sct_add_report_attribute_int,        \
                                           long: sct_add_report_attribute_long,      \
                                           float: sct_add_report_attribute_float,    \
                                           double: sct_add_report_attribute_double,  \
                                           char*: sct_add_report_attribute_string    \
					     )(_1,_2)
#endif

//public-

/*! \brief stop timer if active - do nothing else
    \param[in] it handle of timer to be stopped
*/
void sct_try_stop(const int it);


char *sct_get_timer_cname(const int it);
char *sct_get_event_cname(const int ie);
int sct_get_context_num();
int sct_get_timer_num();
int sct_get_event_num();
int sct_get_pr_thread_count();
int sct_get_pr_thread_id();
sct_context_type *sct_get_context(const int icon);
void sct_del_timer(const int it);
int sct_get_callstats();
int sct_get_eventcounters();
int sct_get_nestedtimers();
int sct_string_recap(sct_string_type *v, int mincap);
sct_string_type sct_string_new(const char *cstring);
void sct_string_copy(const sct_string_type *src, sct_string_type *dst);
void sct_string_delete(sct_string_type *s);

sct_datatype sct_internal_get_stats_datatype();
int sct_internal_get_stats_datatype_size();
void sct_internal_copy_stats(sct_stats_type *src, sct_stats_type *dst, int en);
void sct_internal_alloc_event_arrays(sct_stats_type *s, int en);
void sct_internal_free_event_arrays(sct_stats_type *s);
int sct_internal_get_max_name_length();
int sct_internal_get_name_length(const int it);

sct_attribute_table *sct_internal_create_attribute_table(void);
const sct_attribute_table *sct_internal_get_attribute_table(void);
void sct_internal_free_attribute_table(void);

#endif

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
 * sct_reduce.h: reduce measurment-space after the measurement
 */

#ifndef _H_SCT_REDUCE
#define _H_SCT_REDUCE

#if HAVE_CONFIG_H
#  ifndef _H_CONFIG
#    define _H_CONFIG
#    include <config.h>
#  endif
#endif

//public+
#define SCT_SELECT_ALL -1
#define SCT_REDUCE_ALL -2
//public-

/*! \brief all needed data to do reduction operations */
typedef struct {
#ifdef HAVE_MPI
  // in case we want to distribute the reduction variable:
  MPI_Comm comm;
  int max_proc_name_len;   //!< maximum length of precess name
  char *proc_names;        //!< (char (*)[max_proc_name_len+1]) ; only valid at root
#endif
  int context_choice;      //!< context id
  int proc_choice;         //!< SCT_SELECT_ALL, SCT_REDUCE_ALL, etc
  int thread_choice;       //!< SCT_SELECT_ALL, SCT_REDUCE_ALL, etc
  int sp_merging;          //!< SCT_SP_SELECT_ALL, SCT_SP_SERIAL_ONLY, etc
  int global_timer_num;    //!< number of timers used (accumulated over all MPI-tasks)
  int *timer_num_per_rank; //!< number of timers used per MPI-task
  int *timer_map_per_rank; //!< timer_map to get locally used name from timer_list
  int pid;                 //!< local process id within comm
  int proc_num;            //!< comm_size
  int red_proc_num;        //!< comm_size
  int p_thread_num;        //!< number of threads in parallel region or zero
  int m_thread_num;        //!< number of threads after merging: 1, or p_num
  int r_thread_num;        //!< number of threads after reduction along thread_dim
  int event_num;           //!< number of events used (assumed equal for all MPI-tasks)
  sct_string_type context_name; //!< context name
  sct_stats_type *stats;        //!< p_stats[red_proc_num][timer_num][r_thread_num]
} sct_reduction_type;

typedef struct {
  int psize;
  sct_stats_type *s;
  sct_stats_type **p;
} sct_thread_reduction_type;

typedef struct {
  int proc_size;
  sct_thread_reduction_type *thread_data;
} sct_proc_reduction_type;

//public+
// In the following we define different ways to combine the serial phase measurement
// with the parallel phase measurement - this only matters for OpenMP:

//! The serial phase measurement is added to the master thread of parallel regions.
//! This is the default:
#define SCT_SP_MERGE_SIMPLE 1

//! Only the serial phase measurement is used:
#define SCT_SP_SERIAL_ONLY 2

//! Only the thread-parallel measurement is used:
#define SCT_SP_PARALLEL_ONLY 4

//! The serial phase measurement and the thread-parallel measurement are reported
#define SCT_SP_SELECT_ALL 8
//public-

sct_reduction_type *sct_reduction_new(int context_choice, int proc_choice, int thread_choice, int sp_choice);
void sct_reduction_delete(sct_reduction_type *r);
const char *sct_internal_reduce_doc(int choice);
const char *sct_internal_sp_doc(int sp_merging);
void show_stats(char *label, sct_stats_type *stats);

int sct_get_global_timer_num();
int sct_get_global_idx(int local_ind);
#ifdef HAVE_MPI
void sct_create_global_timer_map(sct_context_type *con);
#endif
char *sct_get_global_timer_cname(const sct_reduction_type *r, int pid, int it);

#endif

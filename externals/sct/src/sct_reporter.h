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

#ifndef _H_SCT_REPORTER
#define _H_SCT_REPORTER

#if HAVE_CONFIG_H
#  ifndef _H_CONFIG
#    define _H_CONFIG
#    include <config.h>
#  endif
#endif

//public+
#include <stdio.h>

/*! \brief print specific timer results
    \param[in] timer_choice SCT_SELECT_ALL or timer_id for output
    \param[in] proc_choice SCT_SELECT_ALL or SCT_REDUCE_ALL MPI-tasks for output
    \param[in] thread_choice SCT_SELECT_ALL or SCT_REDUCE_ALL OpenMP-threads for output
    \param[in] sp_merging merging operation for serial part of OpenMP timer
*/
void sct_single_report(int timer_choice, int proc_choice, int thread_choice, int sp_merging);

/*! \brief print all timer results
    \param[in] proc_choice SCT_SELECT_ALL or SCT_REDUCE_ALL MPI-tasks for output
    \param[in] thread_choice SCT_SELECT_ALL or SCT_REDUCE_ALL OpenMP-threads for output
    \param[in] sp_merging merging operation for serial part of OpenMP timer
*/
void sct_report(int proc_choice, int thread_choice, int sp_merging);

/*! \brief print all or specific  timer results to specified output stream
    \param[in] outstream_arg output stream to write results to
    \param[in] timer_choice SCT_SELECT_ALL or timer_id for output
    \param[in] proc_choice SCT_SELECT_ALL or SCT_REDUCE_ALL MPI-tasks for output
    \param[in] thread_choice SCT_SELECT_ALL or SCT_REDUCE_ALL OpenMP-threads for output
    \param[in] sp_merging merging operation for serial part of OpenMP timer
*/
void sct_stream_report(FILE *outstream_arg, int timer_choice, int proc_choice, int thread_choice, int sp_merging);

//public-

#endif

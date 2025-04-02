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
 * sct_config.h: configuration of sct package
 */

#ifndef _H_SCT_CONFIG
#define _H_SCT_CONFIG


//! number of timers if not set by user
//! can be set to any positive value
#define SCT_DEFAULT_TIMER_SIZE 300


//! maximum length of timer label
#define SCT_LABEL_SIZE 256


//! how deep can timer contexts go?
//! costs: 4 Bytes per context level
#define SCT_MAX_CONTEXT_DEPTH 15


//! how many nested timer levels are allowed without reallocation of datastructures
#define SCT_MAX_NEST_DEPTH 10


//! known real time methods
#define SCT_RTM_UNDEF                    0
#define SCT_RTM_READ_REAL_TIME           1
#define SCT_RTM_OMP_GET_WTIME            2
#define SCT_RTM_MPI_WTIME                4
#define SCT_RTM_CLOCK_GETTIME_MONOTONIC  8
#define SCT_RTM_CLOCK_GETTIME_REALTIME   16
#define SCT_RTM_GETTIMEOFDAY             32
#define SCT_RTM_RDTSCP                   64

#endif

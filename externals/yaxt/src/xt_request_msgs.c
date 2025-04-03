/**
 * @file xt_request_msgs.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_config.h"
#include "xt/xt_request_msgs.h"
#include "xt_request_msgs_internal.h"
#include "xt_mpi_internal.h"
#include "xt_config_internal.h"
#include "xt_request_internal.h"

static void xt_request_msgs_wait(Xt_request request);
static int xt_request_msgs_test(Xt_request request);

#ifdef _OPENMP
static void xt_request_msgs_wait_omp(Xt_request request);
static int xt_request_msgs_test_omp(Xt_request request);
#endif

static const struct Xt_request_vtable request_msgs_vtable = {

  .wait = xt_request_msgs_wait,
  .test = xt_request_msgs_test,
}
#ifdef _OPENMP
  , request_msgs_auto_omp_vtable = {
  .wait = xt_request_msgs_wait_omp,
  .test = xt_request_msgs_test_omp,
};
#endif
  ;

typedef struct Xt_request_msgs_ *Xt_request_msgs;

struct Xt_request_msgs_ {

  const struct Xt_request_vtable *vtable;
  int n;
  MPI_Comm comm;
  MPI_Request requests[];
};

Xt_request
xt_request_msgs_alloc(int n, MPI_Comm comm, Xt_config config)
{
  Xt_request_msgs request
    = xmalloc(sizeof(*request) + (size_t)n * sizeof(MPI_Request) +
              (size_t)n * sizeof(int));
#ifndef _OPENMP
  (void)config;
#else
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  if (mthread_mode == XT_MT_OPENMP)
    request->vtable = &request_msgs_auto_omp_vtable;
  else
#endif
    request->vtable = &request_msgs_vtable;
  request->n = n;
  request->comm = comm;
  return (Xt_request)request;
}


Xt_request
xt_request_msgs_new(int n,
                    const MPI_Request requests[n],
                    MPI_Comm comm)
{
  return xt_request_msgs_custom_new(n, requests, comm, &xt_default_config);
}

Xt_request
xt_request_msgs_custom_new(int n,
                           const MPI_Request requests[n],
                           MPI_Comm comm,
                           Xt_config config)
{
  assert(n >= 0);
  Xt_request_msgs request
    = (Xt_request_msgs)xt_request_msgs_alloc(n, comm, config);

  memcpy(request->requests, requests, (size_t)n * sizeof(*request->requests));

  return (Xt_request)request;
}

MPI_Request *
xt_request_msgs_get_req_ptr(Xt_request request)
{
  Xt_request_msgs request_msgs = (Xt_request_msgs)request;
#ifdef _OPENMP
  assert(request->vtable == &request_msgs_vtable
         || request->vtable == &request_msgs_auto_omp_vtable);
#else
  assert(request->vtable == &request_msgs_vtable);
#endif
  return request_msgs->requests;
}


static void xt_request_msgs_wait(Xt_request request) {

  Xt_request_msgs request_msgs = (Xt_request_msgs)request;

#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  xt_mpi_call(MPI_Waitall(request_msgs->n, request_msgs->requests,
                          MPI_STATUSES_IGNORE), request_msgs->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif

  free(request_msgs);
}

static int xt_request_msgs_test(Xt_request request) {

  Xt_request_msgs request_msgs = (Xt_request_msgs)request;

  size_t n = (size_t)request_msgs->n;
  int *ops_completed_buffer = (int *)(request_msgs->requests + n);
  int flag = xt_mpi_test_some(&(request_msgs->n), request_msgs->requests,
                              ops_completed_buffer, request_msgs->comm);

  if (flag) free(request_msgs);

  return flag;
}

#ifdef _OPENMP
static void
xt_request_msgs_wait_omp(Xt_request request)
{
#pragma omp parallel
  {
    Xt_request_msgs request_msgs = (Xt_request_msgs)request;

#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
    int num_threads = omp_get_num_threads(),
      tid = omp_get_thread_num();
    int nreq = request_msgs->n,
      start_req = (nreq * tid) / num_threads,
      nreq_ = (nreq * (tid+1)) / num_threads - start_req;
    xt_mpi_call(MPI_Waitall(nreq_, request_msgs->requests+start_req,
                            MPI_STATUSES_IGNORE), request_msgs->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif
  }

  free(request);
}

static int
xt_request_msgs_test_omp(Xt_request request) {

  bool flag = true;
#pragma omp parallel reduction(&: flag)
  {
    Xt_request_msgs request_msgs = (Xt_request_msgs)request;
    size_t n = (size_t)request_msgs->n;
    int *ops_completed_buffer = (int *)(request_msgs->requests + n);

    flag &= xt_mpi_test_some_mt(&request_msgs->n, request_msgs->requests,
                                ops_completed_buffer, request_msgs->comm);
  }
  if (flag) free(request);

  return flag;
}
#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

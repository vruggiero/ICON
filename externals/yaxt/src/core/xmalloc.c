/**
 * @file xmalloc.c
 * @brief fail-safe [cm]alloc wrappers
 *
 * @copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
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
#  include <config.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "core/ppm_xfuncs.h"
#include "core/core.h"

#include "core/symprefix.h"

void *
SymPrefix(xcalloc)(size_t nmemb, size_t size, const char *source, int line)
{
  void *p = calloc(nmemb, size);
  if (p || !nmemb || !size)
    ;
  else
    SymPrefix(abort)(SymPrefix(default_comm), strerror(errno), source, line);
  return p;
}


void *
SymPrefix(xmalloc)(size_t size, const char *source, int line)
{
  void *p = malloc(size);
  if (p || !size)
    ;
  else
    SymPrefix(abort)(SymPrefix(default_comm), strerror(errno), source, line);
  return p;
}


void *
SymPrefix(xrealloc)(void *ptr, size_t size, const char *source, int line)
{
  void *p = realloc(ptr, size);
  if (p || !size)
    ;
  else
    SymPrefix(abort)(SymPrefix(default_comm), strerror(errno), source, line);
  return p;
}

/*
 * Local Variables:
 * license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */

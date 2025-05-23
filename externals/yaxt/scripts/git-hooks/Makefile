#
# scripts/release/Makefile --- some automation for development of
#			       release scripts
#
# Copyright  (C)  2021  Thomas Jahns <jahns@dkrz.de>
#
# Author: Thomas Jahns <jahns@dkrz.de>
# Maintainer: Thomas Jahns <jahns@dkrz.de>
# URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#

all: checklicense.exclusions.pcre shellcheck

shellcheck:
	for f in * ; do \
           if { head -n 1 "$$f" | grep -E -q '^#! *(/bin/bash|/bin/sh|/usr/bin/env bash)$$' ; } ; then \
		shellcheck -x "$$f" ; \
	   fi ; done

CC = gcc
CFLAGS = -O2 -march=native -Wall -Wextra

bre2pcre: bre2pcre.c
	$(CC) $(CFLAGS) -o $@ $^ $$(pcre2-config --cflags --libs8)

checklicense.exclusions.pcre: checklicense.exclusions bre2pcre.c
	$(MAKE) bre2pcre && ./bre2pcre <$^ >$@

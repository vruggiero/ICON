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
 * @file mergesort.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sct_mergesort.h"

static void insertionsort_idxpos(strpos_type *w, int n);
static void insertionsort_index(char (*val)[SCT_LABEL_SIZE], int n, int *pos, int reset_pos);


static void insertionsort_idxpos(strpos_type *w, int n) {

  int i, j;
  for(i = 1; i < n; i++) {

    strpos_type tw = w[i];

    for(j = i;  j > 0 && (strcmp(w[j-1].str, tw.str)>0); j--) {
      w[j] = w[j-1];
    }
    w[j] = tw;
  }
}

static void insertionsort_index(char (*val)[SCT_LABEL_SIZE], int n, int *pos, int reset_pos) {

  int i, j;

  if (pos != NULL) {

    if (reset_pos) {
      for(i = 0; i < n; i++) {
        pos[i]=i;
      }
    }

    for(i = 1; i < n; i++) {

      char tv[SCT_LABEL_SIZE];
      memcpy(tv, val[i], SCT_LABEL_SIZE);
      int tp = pos[i];

      for(j = i;  j > 0 && (strcmp(val[j-1],tv)>0); j--) {
        memcpy(val[j], val[j-1], SCT_LABEL_SIZE);
        pos[j] = pos[j-1];
      }

      memcpy(val[j], tv, SCT_LABEL_SIZE);
      pos[j] = tp;
    }

  } else {

    for(i = 1; i < n; i++) {

      char tv[SCT_LABEL_SIZE];
      memcpy(tv, val[i], SCT_LABEL_SIZE);

      for(j = i;  j > 0 && (strcmp(val[j-1],tv)>0); j--) {
        memcpy(val[j], val[j-1], SCT_LABEL_SIZE);
      }

      memcpy(val[j], tv, SCT_LABEL_SIZE);
    }

  }

}

static inline void
my_merge_idxpos(strpos_type * restrict v, strpos_type * restrict w, int lb1, int ub1, int lb2, int ub2) {

  int p  = lb1;
  int i1 = lb1;
  int i2 = lb2;

  while(i1<=ub1 && i2<=ub2) {
    if (strcmp(v[i1].str, v[i2].str) <= 0) {
      w[p]=v[i1];
      i1++;
    } else {
      w[p]=v[i2];
      i2++;
    }
    p++;
  }
  while(i1<=ub1) {
    w[p]=v[i1];
    i1++;
    p++;
  }
  while(i2<=ub2) {
    w[p]=v[i2];
    i2++;
    p++;
  }
}

static void
my_mergesort_idxpos(strpos_type * restrict v, strpos_type * restrict w, int lb, int ub)
{

  int n = ub - lb + 1;
  if (n<2) return;

  if (n<9) {
    insertionsort_idxpos(&v[lb],n);
    return;
  }

  int m   = n/4;

  // 4-way subsort:
  int lb1 = lb;
  int ub1 = lb1 + m - 1;

  int lb2 = ub1+1;
  int ub2 = lb2 + m - 1;

  int lb3 = ub2+1;
  int ub3 = lb3 + m - 1;

  int lb4 = ub3+1;
  int ub4 = ub;

  my_mergesort_idxpos(v, w, lb1, ub1);
  my_mergesort_idxpos(v, w, lb2, ub2);
  my_mergesort_idxpos(v, w, lb3, ub3);
  my_mergesort_idxpos(v, w, lb4, ub4);

  // 2 x 2-way merge v -> w
  my_merge_idxpos(v, w, lb1, ub1, lb2, ub2);
  my_merge_idxpos(v, w, lb3, ub3, lb4, ub4);
  // final merge: w -> v
  my_merge_idxpos(w, v, lb1, ub2, lb3, ub4);

}

void sct_mergesort_idxpos(strpos_type * restrict v, int n)
{

  strpos_type *w;

  w = malloc((size_t)n * sizeof(w[0]));
  my_mergesort_idxpos(v, w, 0, n-1);
  free(w);

}

void sct_mergesort_index(char (*val)[SCT_LABEL_SIZE], int n, int *pos, int reset_pos) {

  // insertion sort for few elements:
  if (n<9) {
    insertionsort_index(val,n,pos,reset_pos);
    return;
  }

  strpos_type *v;
  strpos_type *w;

  v = malloc((size_t)n * sizeof(v[0]));
  w = malloc((size_t)n * sizeof(w[0]));

  if (reset_pos || pos == NULL) {
    for(int i=0; i<n; i++) {
      memcpy(v[i].str, val[i], SCT_LABEL_SIZE);
      v[i].pos = i;
    }
  } else {
    for(int i=0; i<n; i++) {
      memcpy(v[i].str, val[i], SCT_LABEL_SIZE);
      v[i].pos = pos[i];
    }
  }

  my_mergesort_idxpos(v, w, 0, n-1);
  if (pos) {
    for(int i=0; i<n; i++) {
      memcpy(val[i], v[i].str, SCT_LABEL_SIZE);
      pos[i] = v[i].pos;
    }
  } else {
    for(int i=0; i<n; i++) {
      memcpy(val[i], v[i].str, SCT_LABEL_SIZE);
    }
  }

  free(w);
  free(v);
}


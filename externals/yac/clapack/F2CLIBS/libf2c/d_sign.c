// Copyright (c) 1992-2008 The University of Tennessee.  All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause

#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
double d_sign(a,b) doublereal *a, *b;
#else
double d_sign(doublereal *a, doublereal *b)
#endif
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}
#ifdef __cplusplus
}
#endif


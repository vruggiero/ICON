/*                                                                           */
/*+ dace_md5_ifc.c - DACE interface to external MD5 message digest algorithm */
/*                                                                           */
/* Description                                                               */
/*   Provide Fortran-callable interface routines to an RFC-1321 conforming   */
/*   implementation of the MD5 message digest algorithm.                     */
/*   Includes self-test.                                                     */
/*                                                                           */
/* Current Code Owner: DWD, Harald Anlauf                                    */
/*    phone: +49 69 8062 4941                                                */
/*    fax:   +49 69 8062 3721                                                */
/*    email: harald.anlauf@dwd.de                                            */
/*                                                                           */
/* History:                                                                  */
/* Version      Date       Name                                              */
/* ------------ ---------- ----                                              */
/* V2_11        2022/04/14 Harald Anlauf                                     */
/*  initial version                                                          */
/*                                                                           */
/* Code Description:                                                         */
/* Language: C                                                               */
/* Software Standards:                                                       */
/*                                                                           */
/* Authors:                                                                  */
/* Copyright (C) 2022  Harald Anlauf                                         */
/*****************************************************************************/

#include "md5.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

static bool
selftest_ok = false;

static void
md5_get_digest (md5_byte_t digest[16], const char *text, int len)
{
  md5_state_t state;

  /* Use full (C) string length if len is not given.  */
  if (len < 0)
    len = strlen (text);
  md5_init (&state);
  md5_append (&state, (const md5_byte_t *)text, len);
  md5_finish (&state, digest);
}

/* Run the self-test. */
static int
selftest (void)
{
  static const char *const test[] = {
    "", "d41d8cd98f00b204e9800998ecf8427e",
    "a", "0cc175b9c0f1b6a831c399e269772661",
    "abc", "900150983cd24fb0d6963f7d28e17f72",
    "abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b",
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
    "d174ab98d277d9f5a5611c2c9f419d9f",
    NULL, NULL
  };
  int status = 0;
  md5_byte_t digest[16];
  char hex_output[16*2 + 1];

  for (int i = 0; test[i] != NULL; i += 2)
    {
      md5_get_digest (digest, test[i], -1);
      for (int di = 0; di < 16; ++di)
	sprintf (hex_output + di * 2, "%02x", digest[di]);
      if (strcmp (hex_output, test[i+1])) {
	printf ("MD5 (\"%s\") = %s\n", test[i], hex_output);
	printf ("**** ERROR, should be: %s\n", test[i+1]);
	status = 1;
      }
    }
  if (status == 0)
    selftest_ok = true;
  return status;
}

static int
to_base32 (char *dest, const char *src, int len)
{
  static const char base32tab[32] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
  int i, j, k;
  int olen = (len * 8 + 4) / 5;
  for (i = 0, j = 0; i < len; )
    {
      /* Process up to 5 octets = 40 bits with zero-padding.  */
      uint64_t work = 0;
      for (k = 0; k < 5; k++)
	work = (work << 8) + (i < len ? (unsigned char) src[i++] : 0);

      for (k = 7; k >= 0 && olen > 0; k--, olen--)
	dest[j++] = base32tab[ (work >> (5*k)) & 0x1f ];
    }

  /* Pad output to a multiple of 8.  */
  k = j % 8;
  if (k > 0)
    {
      k = 8 - k;
      for (; k > 0; k--)
	dest[j++] = '=';
    }
  return j;
}

static int
to_base64 (char *dest, const char *src, int len)
{
  static const char base64tab[64] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  int i, j, k;
  int olen = (len * 4 + 2) / 3;
  for (i = 0, j = 0; i < len; )
    {
      /* Process up to 3 octets = 24 bits with zero-padding.  */
      uint32_t work = 0;
      for (k = 0; k < 3; k++)
	work = (work << 8) + (i < len ? (unsigned char) src[i++] : 0);

      for (k = 3; k >= 0 && olen > 0; k--, olen--)
	dest[j++] = base64tab[ (work >> (6*k)) & 0x3f ];
    }

  /* Pad if input length was not a multiple of 3.  */
  k = len % 3;
  if (k > 0)
    {
      i = (k == 1) ? 2 : 1;
      for (; i > 0; i--)
	dest[j++] = '=';
    }
  return j;
}

static int
test_to_base ()
{
  /* Test vectors and results from RFC4648 for base64 and base32.  */
  static const char *const test[] = {
    "", "", "",
    "f", "Zg==", "MY======",
    "fo", "Zm8=", "MZXQ====",
    "foo", "Zm9v", "MZXW6===",
    "foob", "Zm9vYg==", "MZXW6YQ=",
    "fooba", "Zm9vYmE=", "MZXW6YTB",
    "foobar", "Zm9vYmFy", "MZXW6YTBOI======",
    NULL, NULL, NULL
  };
  int status = 0;

  for (int i = 0; test[i] != NULL; i += 3)
    {
      char d[65];
      int j, l = strlen (test[i]);
      j = to_base64 (d, test[i], l);
      d[j] = '\0';
      if (strcmp (d, test[i+1])) {
	printf ("base64 (\"%s\") = %s != %s\n", test[i], d, test[i+1]);
	status = 1;
      }
      j = to_base32 (d, test[i], l);
      d[j] = '\0';
      if (strcmp (d, test[i+2])) {
	printf ("base32 (\"%s\") = %s != %s\n", test[i], d, test[i+2]);
	status = 1;
      }
    }
  return status;
}

/*
  interface
     subroutine dace_md5 (text, len, base,      &
                          digest, string, status) bind(c,name="dace_md5")
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       character(kind=c_char), intent(in)  :: text(*)
       integer(c_int), value,  intent(in)  :: len
       integer(c_int), value,  intent(in)  :: base
       character(kind=c_char), intent(out) :: digest(16)
       character(kind=c_char), intent(out) :: string(32)
       integer(c_int),         intent(out) :: status
     end subroutine dace_md5
  end interface
 */

void
dace_md5 (const char *text, int len, int base,
	  char *digest, char *string, int *status)
{
  md5_byte_t md5_digest[16];

  if (!selftest_ok && (selftest () != 0 || test_to_base () != 0))
    {
      *status = 2;
      return;
    }

  md5_get_digest (md5_digest, text, len);
  memcpy (digest, md5_digest, sizeof (md5_digest));
  *status = 0;

  if (string)
    {
      if (base == 16)
	{
	  for (int i = 0; i < 16; ++i)
	    sprintf (string + i * 2, "%02x", md5_digest[i]);
	}
      else if (base == 32)
	{
	  memset (string, ' ', 32);
	  to_base32 (string, digest, sizeof (md5_digest));
	}
      else if (base == 64)
	{
	  memset (string, ' ', 32);
	  to_base64 (string, digest, sizeof (md5_digest));
	}
      else
	{
	  memset (string, '*', 32);
	  *status = 1;
	}
    }
  return;
}

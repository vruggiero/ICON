# taken from
# http://git.savannah.gnu.org/gitweb/?p=autoconf.git;a=blob_plain;f=lib/autoconf/c.m4;hb=refs/heads/master
# to define search for c11 compatible compiler

# AC_PROG_CC([COMPILER ...])
# --------------------------
# COMPILER ... is a space separated list of C compilers to search for.
# This just gives the user an opportunity to specify an alternative
# search list for the C compiler.
AN_MAKEVAR([CC],  [AC_PROG_CC])
AN_PROGRAM([cc],  [AC_PROG_CC])
AN_PROGRAM([gcc], [AC_PROG_CC])
AC_DEFUN_ONCE([AX_PROG_CC11],
[AC_LANG_PUSH(C)dnl
AC_ARG_VAR([CC],     [C compiler command])dnl
AC_ARG_VAR([CFLAGS], [C compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
_AC_ARG_VAR_LIBS()dnl
_AC_ARG_VAR_CPPFLAGS()dnl
m4_ifval([$1],
      [AC_CHECK_TOOLS(CC, [$1])],
[AC_CHECK_TOOL(CC, gcc)
if test -z "$CC"; then
  dnl Here we want:
  dnl	AC_CHECK_TOOL(CC, cc)
  dnl but without the check for a tool without the prefix.
  dnl Until the check is removed from there, copy the code:
  if test -n "$ac_tool_prefix"; then
    AC_CHECK_PROG(CC, [${ac_tool_prefix}cc], [${ac_tool_prefix}cc])
  fi
fi
if test -z "$CC"; then
  AC_CHECK_PROG(CC, cc, cc, , , /usr/ucb/cc)
fi
if test -z "$CC"; then
  AC_CHECK_TOOLS(CC, cl.exe)
fi
if test -z "$CC"; then
  AC_CHECK_TOOL(CC, clang)
fi
])

test -z "$CC" && AC_MSG_FAILURE([no acceptable C compiler found in \$PATH])

# Provide some information about the compiler.
_AS_ECHO_LOG([checking for _AC_LANG compiler version])
set X $ac_compile
ac_compiler=$[2]
for ac_option in --version -v -V -qversion -version; do
  _AC_DO_LIMIT([$ac_compiler $ac_option >&AS_MESSAGE_LOG_FD])
done

m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
_AC_LANG_COMPILER_GNU
if test $ac_compiler_gnu = yes; then
  GCC=yes
else
  GCC=
fi
_AC_PROG_CC_G
dnl
dnl Set ac_prog_cc_stdc to the supported C version.
dnl Also set the documented variable ac_cv_prog_cc_stdc;
dnl its name was chosen when it was cached, but it is no longer cached.
_AC_PROG_CC_C11([ac_prog_cc_stdc=c11
		 ac_cv_prog_cc_stdc=$ac_cv_prog_cc_c11],
  [_AC_PROG_CC_C99([ac_prog_cc_stdc=c99
		    ac_cv_prog_cc_stdc=$ac_cv_prog_cc_c99],
     [_AC_PROG_CC_C89([ac_prog_cc_stdc=c89
		       ac_cv_prog_cc_stdc=$ac_cv_prog_cc_c89],
		      [ac_prog_cc_stdc=no
		       ac_cv_prog_cc_stdc=no])])])
dnl
AC_LANG_POP(C)dnl
])# AX_PROG_CC11


# _AC_PROG_CC_C89 ([ACTION-IF-AVAILABLE], [ACTION-IF-UNAVAILABLE])
# ----------------------------------------------------------------
# If the C compiler is not in ANSI C89 (ISO C90) mode by default, try
# to add an option to output variable CC to make it so.  This macro
# tries various options that select ANSI C89 on some system or
# another.  It considers the compiler to be in ANSI C89 mode if it
# handles function prototypes correctly.
AC_DEFUN([_AC_PROG_CC_C89],
[_AC_C_STD_TRY([c89],
[[#include <stdarg.h>
#include <stdio.h>
struct stat;
/* Most of the following tests are stolen from RCS 5.7's src/conf.sh.  */
struct buf { int x; };
FILE * (*rcsopen) (struct buf *, struct stat *, int);
static char *e (p, i)
     char **p;
     int i;
{
  return p[i];
}
static char *f (char * (*g) (char **, int), char **p, ...)
{
  char *s;
  va_list v;
  va_start (v,p);
  s = g (p, va_arg (v,int));
  va_end (v);
  return s;
}

/* OSF 4.0 Compaq cc is some sort of almost-ANSI by default.  It has
   function prototypes and stuff, but not '\xHH' hex character constants.
   These don't provoke an error unfortunately, instead are silently treated
   as 'x'.  The following induces an error, until -std is added to get
   proper ANSI mode.  Curiously '\x00'!='x' always comes out true, for an
   array size at least.  It's necessary to write '\x00'==0 to get something
   that's true only with -std.  */
int osf4_cc_array ['\x00' == 0 ? 1 : -1];

/* IBM C 6 for AIX is almost-ANSI by default, but it replaces macro parameters
   inside strings and character constants.  */
#define FOO(x) 'x'
int xlc6_cc_array[FOO(a) == 'x' ? 1 : -1];

int test (int i, double x);
struct s1 {int (*f) (int a);};
struct s2 {int (*f) (double a);};
int pairnames (int, char **, FILE *(*)(struct buf *, struct stat *, int), int, int);
int argc;
char **argv;]],
[[return f (e, argv, 0) != argv[0]  ||  f (e, argv, 1) != argv[1];]],
dnl Don't try gcc -ansi; that turns off useful extensions and
dnl breaks some systems' header files.
dnl AIX circa 2003	-qlanglvl=extc89
dnl old AIX		-qlanglvl=ansi
dnl Ultrix, OSF/1, Tru64	-std
dnl HP-UX 10.20 and later	-Ae
dnl HP-UX older versions	-Aa -D_HPUX_SOURCE
dnl SVR4			-Xc -D__EXTENSIONS__
[-qlanglvl=extc89 -qlanglvl=ansi -std \
	-Ae "-Aa -D_HPUX_SOURCE" "-Xc -D__EXTENSIONS__"], [$1], [$2])[]dnl
])# _AC_PROG_CC_C89


# _AC_C_C99_TEST_HEADER
# ---------------------
# A C header suitable for testing for C99.
AC_DEFUN([_AC_C_C99_TEST_HEADER],
[[#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <wchar.h>
#include <stdio.h>

// Check varargs macros.  These examples are taken from C99 6.10.3.5.
#define debug(...) fprintf (stderr, __VA_ARGS__)
#define showlist(...) puts (#__VA_ARGS__)
#define report(test,...) ((test) ? puts (#test) : printf (__VA_ARGS__))
static void
test_varargs_macros (void)
{
  int x = 1234;
  int y = 5678;
  debug ("Flag");
  debug ("X = %d\n", x);
  showlist (The first, second, and third items.);
  report (x>y, "x is %d but y is %d", x, y);
}

// Check long long types.
#define BIG64 18446744073709551615ull
#define BIG32 4294967295ul
#define BIG_OK (BIG64 / BIG32 == 4294967297ull && BIG64 % BIG32 == 0)
#if !BIG_OK
  your preprocessor is broken;
#endif
#if BIG_OK
#else
  your preprocessor is broken;
#endif
static long long int bignum = -9223372036854775807LL;
static unsigned long long int ubignum = BIG64;

struct incomplete_array
{
  int datasize;
  double data[];
};

struct named_init {
  int number;
  const wchar_t *name;
  double average;
};

typedef const char *ccp;

static inline int
test_restrict (ccp restrict text)
{
  // See if C++-style comments work.
  // Iterate through items via the restricted pointer.
  // Also check for declarations in for loops.
  for (unsigned int i = 0; *(text+i) != '\0'; ++i)
    continue;
  return 0;
}

// Check varargs and va_copy.
static bool
test_varargs (const char *format, ...)
{
  va_list args;
  va_start (args, format);
  va_list args_copy;
  va_copy (args_copy, args);

  const char *str = "";
  int number = 0;
  float fnumber = 0;

  while (*format)
    {
      switch (*format++)
	{
	case 's': // string
	  str = va_arg (args_copy, const char *);
	  break;
	case 'd': // int
	  number = va_arg (args_copy, int);
	  break;
	case 'f': // float
	  fnumber = va_arg (args_copy, double);
	  break;
	default:
	  break;
	}
    }
  va_end (args_copy);
  va_end (args);

  return *str && number && fnumber;
}]])# _AC_C_C99_TEST_HEADER

# _AC_C_C99_TEST_BODY
# -------------------
# A C body suitable for testing for C99, assuming the corresponding header.
AC_DEFUN([_AC_C_C99_TEST_BODY],
[[
  // Check bool.
  _Bool success = false;

  // Check restrict.
  if (test_restrict ("String literal") == 0)
    success = true;
  char *restrict newvar = "Another string";

  // Check varargs.
  success &= test_varargs ("s, d' f .", "string", 65, 34.234);
  test_varargs_macros ();

  // Check flexible array members.
  struct incomplete_array *ia =
    malloc (sizeof (struct incomplete_array) + (sizeof (double) * 10));
  ia->datasize = 10;
  for (int i = 0; i < ia->datasize; ++i)
    ia->data[i] = i * 1.234;

  // Check named initializers.
  struct named_init ni = {
    .number = 34,
    .name = L"Test wide string",
    .average = 543.34343,
  };

  ni.number = 58;

  int dynamic_array[ni.number];
  dynamic_array[ni.number - 1] = 543;

  // work around unused variable warnings
  return (!success || bignum == 0LL || ubignum == 0uLL || newvar[0] == 'x'
	  || dynamic_array[ni.number - 1] != 543);
]])


# _AC_PROG_CC_C99 ([ACTION-IF-AVAILABLE], [ACTION-IF-UNAVAILABLE])
# ----------------------------------------------------------------
# If the C compiler is not in ISO C99 mode by default, try to add an
# option to output variable CC to make it so.  This macro tries
# various options that select ISO C99 on some system or another.  It
# considers the compiler to be in ISO C99 mode if it handles _Bool,
# // comments, flexible array members, inline, long long int, mixed
# code and declarations, named initialization of structs, restrict,
# va_copy, varargs macros, variable declarations in for loops and
# variable length arrays.
AC_DEFUN([_AC_PROG_CC_C99],
[_AC_C_STD_TRY([c99],
[_AC_C_C99_TEST_HEADER],
[_AC_C_C99_TEST_BODY],
dnl Try
dnl GCC		-std=gnu99 (unused restrictive modes: -std=c99 -std=iso9899:1999)
dnl IBM XL C	-qlanglvl=extc1x (V12.1; does not pass C11 test)
dnl IBM XL C	-qlanglvl=extc99
dnl		(pre-V12.1; unused restrictive mode: -qlanglvl=stdc99)
dnl HP cc	-AC99
dnl Intel ICC	-std=c99, -c99 (deprecated)
dnl IRIX	-c99
dnl Solaris	-D_STDC_C99=
dnl		cc's -xc99 option uses linker magic to define the external
dnl		symbol __xpg4 as if by "int __xpg4 = 1;", which enables C99
dnl		behavior for C library functions.  This is not wanted here,
dnl		because it means that a single module compiled with -xc99
dnl		alters C runtime behavior for the entire program, not for
dnl		just the module.  Instead, define the (private) symbol
dnl		_STDC_C99, which suppresses a bogus failure in <stdbool.h>.
dnl		The resulting compiler passes the test case here, and that's
dnl		good enough.  For more, please see the thread starting at:
dnl            http://lists.gnu.org/archive/html/autoconf/2010-12/msg00059.html
dnl Tru64	-c99
dnl with extended modes being tried first.
[[-std=gnu99 -std=c99 -c99 -AC99 -D_STDC_C99= -qlanglvl=extc1x -qlanglvl=extc99]], [$1], [$2])[]dnl
])# _AC_PROG_CC_C99


# _AC_PROG_CC_C11 ([ACTION-IF-AVAILABLE], [ACTION-IF-UNAVAILABLE])
# ----------------------------------------------------------------
# If the C compiler is not in ISO C11 mode by default, try to add an
# option to output variable CC to make it so.  This macro tries
# various options that select ISO C11 on some system or another.  It
# considers the compiler to be in ISO C11 mode if it handles _Alignas,
# _Alignof, _Noreturn, _Static_assert, UTF-8 string literals,
# duplicate typedefs, and anonymous structures and unions.
AC_DEFUN([_AC_PROG_CC_C11],
[_AC_C_STD_TRY([c11],
[_AC_C_C99_TEST_HEADER[
// Check _Alignas.
char _Alignas (double) aligned_as_double;
char _Alignas (0) no_special_alignment;
extern char aligned_as_int;
char _Alignas (0) _Alignas (int) aligned_as_int;

// Check _Alignof.
enum
{
  int_alignment = _Alignof (int),
  int_array_alignment = _Alignof (int[100]),
  char_alignment = _Alignof (char)
};
_Static_assert (0 < -_Alignof (int), "_Alignof is signed");

// Check _Noreturn.
int _Noreturn does_not_return (void) { for (;;) continue; }

// Check _Static_assert.
struct test_static_assert
{
  int x;
  _Static_assert (sizeof (int) <= sizeof (long int),
                  "_Static_assert does not work in struct");
  long int y;
};

// Check UTF-8 literals.
#define u8 syntax error!
char const utf8_literal[] = u8"happens to be ASCII" "another string";

// Check duplicate typedefs.
typedef long *long_ptr;
typedef long int *long_ptr;
typedef long_ptr long_ptr;

// Anonymous structures and unions -- taken from C11 6.7.2.1 Example 1.
struct anonymous
{
  union {
    struct { int i; int j; };
    struct { int k; long int l; } w;
  };
  int m;
} v1;
]],
[_AC_C_C99_TEST_BODY[
  v1.i = 2;
  v1.w.k = 5;
  _Static_assert ((offsetof (struct anonymous, i)
		   == offsetof (struct anonymous, w.k)),
		  "Anonymous union alignment botch");
]],
dnl Try
dnl GCC		-std=gnu11 (unused restrictive mode: -std=c11)
dnl with extended modes being tried first.
dnl
dnl Do not try -qlanglvl=extc1x, because IBM XL C V12.1 (the latest version as
dnl of September 2012) does not pass the C11 test.  For now, try extc1x when
dnl compiling the C99 test instead, since it enables _Static_assert and
dnl _Noreturn, which is a win.  If -qlanglvl=extc11 or -qlanglvl=extc1x passes
dnl the C11 test in some future version of IBM XL C, we'll add it here,
dnl preferably extc11.
[[-std=gnu11]], [$1], [$2])[]dnl
])# _AC_PROG_CC_C11


# AC_C__GENERIC
# -------------
# Define HAVE_C__GENERIC if _Generic works, a la C11.
AN_IDENTIFIER([_Generic], [AC_C__GENERIC])
AC_DEFUN([AC_C__GENERIC],
[AC_CACHE_CHECK([for _Generic], ac_cv_c__Generic,
[AC_COMPILE_IFELSE(
   [AC_LANG_SOURCE(
      [[int
	 main (int argc, char **argv)
	 {
	   int a = _Generic (argc, int: argc = 1);
	   int *b = &_Generic (argc, default: argc);
	   char ***c = _Generic (argv, int: argc, default: argv ? &argv : 0);
	   _Generic (1 ? 0 : b, int: a, default: b) = &argc;
	   _Generic (a = 1, default: a) = 3;
	   return a + !b + !c;
	 }
      ]])],
   [ac_cv_c__Generic=yes],
   [ac_cv_c__Generic=no])])
if test $ac_cv_c__Generic = yes; then
  AC_DEFINE([HAVE_C__GENERIC], 1,
	    [Define to 1 if C11-style _Generic works.])
fi
])# AC_C__GENERIC

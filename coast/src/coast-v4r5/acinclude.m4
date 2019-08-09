################################################
#
# COAST configuration macros
#
# Author: Ralf Ulrich 
# Date: Mo 30. Jul 14:01:57 CEST 2007
# Version: $Id: acinclude.m4,v 1.3 2007-08-23 12:41:19 rulrich Exp $
#
################################################

# ----------------------------------------------
#                 DOXYGEN
# ----------------------------------------------


# COAST_WITH_DOXYGEN
#   Check for presence of the DOXYGEN package, either from command line
#   or by searching.
#
#   If found, provide substitution variables:
#      DOXYGEN
#
AC_DEFUN([COAST_WITH_DOXYGEN],
[
  AM_CONDITIONAL([HAVE_DOXYGEN], test xyes = xyes)
])



# ----------------------------------------------
#                 COAST_USERLIB
# ----------------------------------------------


# COAST_WITH_USERLIB
#   Check for presence of a COAST USER LIB, and include it 
#   in build
#
#   If found, provide substitution variables:
#      COAST_USER_LIB
#
AC_DEFUN([COAST_WITH_USERLIB],
[
#  if test "x$COAST_USER_LIB" = "x"; then
#    AC_MSG_ERROR([\$COAST_DIR environment variable needs to be defined! ABORTING ...])
#  fi
#  AM_CONDITIONAL([HAVE_USERLIB], test xyes = xyes)
])



# ----------------------------------------------
#                 ROOT
# ----------------------------------------------


# COAST_WITH_ROOT
#   Check for presence of the ROOT package, either from command line
#   or by searching.
#
#   If found, provide substitution variables:
#      ROOTSYS
#      ROOT_LIBS
#      ROOT_GLIBS
#      ROOT_CXXFLAGS
#      ROOTCINT
#
AC_DEFUN([COAST_WITH_ROOT],
[_COAST_WITH_PACKAGE([root],
  [AC_HELP_STRING([--with-root],
                  [use ROOT (default ROOTSYS=ARG or use ROOTSYS or search)])],
  [/usr /usr/local /opt/local /usr/local/root /usr/local/root/pro/root $ROOTSYS],
  [bin/root-config])
if test x$coast_cv_root_root != xno ; then
  ROOTSYS=$coast_cv_root_root
  export ROOTSYS
  AC_SUBST([ROOTSYS])
  ROOT_LIBS=`$ROOTSYS/bin/root-config --libs`
  AC_SUBST(ROOT_LIBS)
  ROOT_GLIBS=`$ROOTSYS/bin/root-config --glibs`
  AC_SUBST(ROOT_GLIBS)
  ROOT_LIBDIR=`$ROOTSYS/bin/root-config --libdir 2>/dev/null` || \
  ROOT_LIBDIR=`$ROOTSYS/bin/root-config --libs | $AWK '{print substr($[]1, 3)}'`
  AC_SUBST(ROOT_LIBDIR)
  ROOT_CXXFLAGS=`$ROOTSYS/bin/root-config --cflags`
  AC_SUBST(ROOT_CXXFLAGS)
  ROOTCINT=$ROOTSYS/bin/rootcint
  AC_SUBST(ROOTCINT)
  
  AC_LANG_PUSH([C++])
  _COAST_SAVE_COMPILE_FLAGS
  CXXFLAGS="$CXXFLAGS $ROOT_CXXFLAGS"
  LDFLAGS="-Wl,--no-as-needed $LDFLAGS $ROOT_LIBS"
  AC_MSG_CHECKING([whether root can be used])
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <TCanvas.h>],[TCanvas* c = new TCanvas("test");])],
      [
        AC_MSG_RESULT([yes])
        AC_DEFINE([HAVE_ROOT], 1)
      ],[
        AC_MSG_RESULT([no - not binary compatible to selection])
        coast_cv_root_root="no"
      ])
  _COAST_RESTORE_COMPILE_FLAGS
  AC_LANG_POP([C++])
  
else
  COAST_MISSING([ROOT])
fi
AH_TEMPLATE([HAVE_ROOT], [Define to 1 if ROOT is available on the system])
AM_CONDITIONAL([HAVE_ROOT], test x$coast_cv_root_root != xno)
])



# COAST_REQUIRE_ROOT([version])
#   Require that ROOT is installed, if not, configuration fails
#   If version argument present, require that ROOT as at least
#   this version.
#
AC_DEFUN([COAST_REQUIRE_ROOT],
[AC_REQUIRE([COAST_WITH_ROOT])
if test x$coast_cv_root_root = xno ; then
  COAST_MISSING_REQUIRED([ROOT])
fi
m4_if($1,[],[dnl],
  [AC_LANG_PUSH([C++])
   _COAST_SAVE_COMPILE_FLAGS
  CPPFLAGS="$CPPFLAGS $ROOT_CPPFLAGS"
  LDFLAGS="$LDFLAGS $ROOT_LIBS"
  LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ROOT_LIBDIR"
  export LD_LIBRARY_PATH
  AC_RUN_IFELSE([#include <RVersion.h>
                 int main ()
                 {
                   if (ROOT_VERSION_CODE >= 
                       ROOT_VERSION(_COAST_NORMALISE_VERSION_3($1))) return 0;
                   else return 1;
                 }], [], [AC_MSG_NOTICE([ROOT > $1 required], 1)
                          COAST_MISSING_REQUIRED([ROOT])])
  _COAST_RESTORE_COMPILE_FLAGS
  AC_LANG_POP([C++])
  ])
])



# COAST_INIT
#   Initial code for extra scripts to generate
#
AC_DEFUN([COAST_INIT],
[coast_missing=
coast_missing_required=
])


# COAST_MISSING
#   Deal with missing pacakges
#
AC_DEFUN([COAST_MISSING],
[coast_missing="$coast_missing $1"
])


# COAST_MISSING_REQUIRED
#   Deal with missing pacakges
#
AC_DEFUN([COAST_MISSING_REQUIRED],
[coast_missing_required="$coast_missing_required $1"
])


# COAST_CHECK_MISSING_REQUIRED_PACKAGES
#   Code to test if we are missing any required packages.
#   Calls COAST_FINISH if this is the case.
#
AC_DEFUN([COAST_CHECK_MISSING_REQUIRED_PACKAGES],
[
if test "x$coast_missing_required" != x
then
  COAST_FINISH
fi  
])

# COAST_FINISH
#   Final code after all tests for external packages
#   Aborts execution if we are missing any required packages.
#
AC_DEFUN([COAST_FINISH],
[
if test "x$coast_missing" != x
then
  AC_MSG_NOTICE([Optional packages not found or deactivated on this system:])
  AC_MSG_NOTICE([    $coast_missing])
fi
if test "x$coast_missing_required" != x
then
  AC_MSG_NOTICE([Missing required packages: $coast_missing_required])
  AC_MSG_ERROR([Aborting configuration])
fi
])

# _COAST_WITH_PACKAGE(name, help, path, file)
#   Wrapper for AC_ARG_WITH and AC_CACHE_CHECK
#   Provides a --with-name option with help.
#   If package location not given, try to find `file' in `path' and
#   keep the *last* match
#
AC_DEFUN([_COAST_WITH_PACKAGE], 
[AC_ARG_WITH([$1], [$2],
             [_coast_with_$1=$withval])
AC_CACHE_CHECK([for $1],
               [coast_cv_root_$1],
               [if test x$_coast_with_$1 = xno
                then
                  coast_cv_root_$1=no
                fi
                if test x$_coast_with_$1 = x || test x$_coast_with_$1 = xyes
                then
                  for p in $3
                  do
                    test -f $p/$4 && coast_cv_root_$1=$p
                  done
                  test x$coast_cv_root_$1 = x && coast_cv_root_$1=no
                else
                  coast_cv_root_$1=$_coast_with_$1
                fi
               ])
])


##########################################################
#
#     M4 macro section
#


# _COAST_BASENAME
#  M4 analogue to the basename shell utility
#
m4_define([_COAST_BASENAME], [m4_bregexp($1, [\(.*/\|\)\(.+\)], [\2])])

# _COAST_NORMALISE_VERSION_3
#   Convert dot-separated version number into a 
#   comma-separated list for individual processing
#   Up to 3 components accepted, missing componented
#   replaced by 0.
#
m4_define([_COAST_NORMALISE_VERSION_3], 
          [_COAST_NORMALISE_VERSION_AUX_3(m4_bpatsubst($1, [\.], [,]))])

m4_define([_COAST_NORMALISE_VERSION_AUX_3], 
  [ifelse($#, 3, [$1, $2, $3], $#, 2, [$1, $2, 0], $#, 1, [$1, 0, 0], [AC_FATAL([Maximal 3 components in version expected - got $*])])])

# _COAST_SAVE_COMPILE_FLAGS
#   Save flags for compilation so we can modify them for test compilation
#   NB: This macro cannot be nested
#
m4_define([_COAST_SAVE_COMPILE_FLAGS],
[_coast_save_CXXFLAGS="$CXXFLAGS"
_coast_save_CPPFLAGS="$CPPFLAGS"
_coast_save_FFLAGS="$FFLAGS"
_coast_save_CFLAGS="$CFLAGS"
_coast_save_LIBS="$LIBS"
_coast_save_LDFLAGS="$LDFLAGS"
_coast_save_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
])

# _COAST_SAVE_RESTORE_FLAGS
#   Restore flags after compilation. so we can modify them for test compilation
#
m4_define([_COAST_RESTORE_COMPILE_FLAGS],
[CXXFLAGS="$_coast_save_CXXFLAGS"
CPPFLAGS="$_coast_save_CPPFLAGS"
FFLAGS="$_coast_save_FFLAGS"
CFLAGS="$_coast_save_CFLAGS"
LIBS="$_coast_save_LIBS"
LDFLAGS="$_coast_save_LDFLAGS"
LD_LIBRARY_PATH="$_coast_save_LD_LIBRARY_PATH"
])

AC_INIT([spade], 1.0.0)
AC_PREREQ([2.59])

m4_include([m4/ac_openmp.m4])

# Find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "could not determine R_HOME"
    exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP


	
AC_LANG_PUSH([C++])

AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS,${OPENMP_CXXFLAGS})

# Check if parallel STL is available (only in newer g++ versions)
AC_CHECK_HEADERS([parallel/algorithm],
		[AC_DEFINE([HAVE_PARALLEL_STL],[1],[Define if parallel STL available.])],
		[AC_MSG_RESULT("Parallel STL not found using default")])

AC_LANG_POP([C++])

AC_SUBST(FLOWSPD_CXXFLAGS,${DEFS})
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

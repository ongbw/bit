# Process this file with autoconf to produce a configure script.

AC_INIT([libtree3d], [0.0], [ongbw@msu.edu])
AC_CONFIG_SRCDIR([src/tree3d.cpp])
AM_INIT_AUTOMAKE([1.9])
AC_PROG_LIBTOOL

# use the C++ compiler for the following checks
AC_LANG([C++])

# Checks for programs.
AC_PROG_CXX
AC_OPENMP

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([tests/Makefile])

AC_OUTPUT

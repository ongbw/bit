# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = --pedantic -Wall $(OPENMP_CXXFLAGS)

lib_LTLIBRARIES = libtree3d.la
libtree3d_la_SOURCES = tree3d.cpp
include_HEADERS = problem.h constants.h solver.h tree3d.h
libtree3d_la_LDFLAGS = -version-info 0:0:0

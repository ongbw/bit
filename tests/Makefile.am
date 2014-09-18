AM_CXXFLAGS = --pedantic -Wall $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = -I$(top_srcdir)/src

TESTS = directsum_mp treecode_mp treecode
check_PROGRAMS = directsum_mp treecode_mp treecode
directsum_mp_SOURCES = directsum_mp.cpp  
directsum_mp_LDADD = $(top_builddir)/src/libtree3d.la 
treecode_mp_SOURCES = treecode_mp.cpp  
treecode_mp_LDADD = $(top_builddir)/src/libtree3d.la 
treecode_SOURCES = treecode.cpp  
treecode_LDADD = $(top_builddir)/src/libtree3d.la 

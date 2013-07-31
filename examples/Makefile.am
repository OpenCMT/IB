
include $(top_srcdir)/Common.am

#AM_DEFAULT_SOURCE_EXT = .cpp

# if HAVE_CHECK
TESTS =  LNL_SijGuess \
	 LNL_SijGuess_Trim

# else
# TEST =
# endif

LDADD = ../libmutomIB-0.3.5.la

check_PROGRAMS = $(TESTS)


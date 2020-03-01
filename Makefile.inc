CC := gcc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef HAVE_LAPACK
	CPPFLAGS += -DHAVE_LAPACK
	LDLIBS += -lopenblas -lgfortran -lpthread	
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
else
  CFLAGS += -g -Wall -Wextra -ansi -Wpedantic #-Werror
endif

CF := $(CFLAGS) $(CPPFLAGS) -fPIC



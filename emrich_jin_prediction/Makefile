# Builds with OpenMP when the compiler supports it (gcc, Homebrew gcc, clang+libgomp).  On macOS the default `cc`
# is Apple clang, which has no OpenMP runtime: it uses Homebrew libomp if the header is actually installed
# (`brew install libomp`); otherwise it builds serially and says so.  Force a compiler with `make CC=gcc-14`.
CC ?= gcc
CFLAGS = -O2
OMP_OK := $(shell echo 'int main(void){return 0;}' | $(CC) -fopenmp -x c - -o /dev/null 2>/dev/null && echo yes)
ifeq ($(OMP_OK),yes)
  OMPFLAGS = -fopenmp
  OMPLIBS =
else
  LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(wildcard $(LIBOMP)/include/omp.h),)
    OMPFLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
    OMPLIBS = -L$(LIBOMP)/lib -lomp
  else
    OMPFLAGS =
    OMPLIBS =
    $(warning OpenMP not found: serial build.  macOS: `brew install libomp`, or `make CC=gcc-14` with Homebrew gcc)
  endif
endif
all: emrich_jin_sde
emrich_jin_sde: emrich_jin_sde.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $< $(OMPLIBS) -lm
clean:
	rm -f emrich_jin_sde

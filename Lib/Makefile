#####################################################################
## MAKEFILE to compile the Kalman Libraby 						   ##
#####################################################################

# To just compile the library, do:
#       make
# resp.  
#       make VERBOSE=1
#       make CC=gcc OPENMP=1
#       make CC=/opt/local/bin/clang-mp-3.7 OPT=1 OPENMP=1
#
# To run a demo using the library 
#       cd ../Demo ; make
#

LIBRARY = libkal

include ../Makefile.inc

RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f
KalIncDir = ../Include
LahIncDir = ../../linear-algebra-helpers/Include
I = -I$(KalIncDir) -I$(LahIncDir)

AR_TARGET = $(LIBRARY).a
SO_TARGET = $(LIBRARY).so

all: install

KAL_OBJ = \
 kal_util.o \
 kal_predict.o \
 kal_correct.o \
 ukal_util.o \
 ukal_predict.o \
 ukal_correct.o \
 ukal_sigmaPoints.o

$(KAL_OBJ): $(KalIncDir)/kalman.h $(LahIncDir)/lah.h Makefile

%.o: ../Source/%.c  $(KalIncDir)/kalman.h $(LahIncDir)/lah.h 
	$(CC) $(CF) $(I) -c $< 

static: $(AR_TARGET)

$(SO_TARGET): $(KAL_OBJ)
	$(CC) $(CF) $(KAL_OBJ) -o $@ $(LDFLAGS) -shared

$(AR_TARGET): $(KAL_OBJ)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install archive in this directory
install: $(AR_TARGET) $(SO_TARGET)

.PHONY: clean purge

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*


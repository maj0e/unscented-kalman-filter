include ../Makefile.inc

LDLIBS += -lm
LahDir = ../../linear-algebra-helpers
CS := $(LDFLAGS) ../Lib/libkal.a $(LahDir)/Lib/liblah.a $(LDLIBS)

KalIncDir = ../Include
LahIncDir = $(LahDir)/Include
I = -I$(KalIncDir) -I$(LahIncDir)

INC_DIR := ../../Include
BUILD_DIR := ./build

all : test_UnscentedKalman

test_Kalman: test_Kalman.o
		$(CC) $(CF) $(I) -o $(BUILD_DIR)/$@ $(BUILD_DIR)/$< $(CS)
		
%.o : %.c 
		$(CC) $(CF) $(I) -o $(BUILD_DIR)/$@ -c $<
	
.PHONY: clean purge

clean:
	- $(RM) .build/*.o

purge: clean
	- $(RM) -r ./build
	- mkdir ./build

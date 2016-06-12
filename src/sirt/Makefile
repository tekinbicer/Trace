# Compiler
CC = clang++

# Directories
ROOTDIR = ../..
# Bin directory
BINDIR = ${ROOTDIR}/bin/


# ------ Flags ------
# TIMERON: 				Collects time information and prints this at the end ot exec.
# PREFETCHON:			Enables hw level prefetching during the reconstruction. May
# 								not work on all CPUs.
# BGPM_PROFILER: 	Collects BGPM profiler information and prints at the end of 
# 								exec. Currently only L2 cache and FLOPS information are 
# 								collected. Available only on BG/Q machines.
CFLAGS = -std=c++11 -O2 -Wall -Wextra 
CFLAGS += -g -ffast-math -fPIC -DTIMERON #-DBGPM_PROFILER #-DPREFETCHON #-mprfchw 
#CFLAGS += -Rpass-analysis=loop-vectorize -S -mllvm --x86-asm-syntax=intel

# Common data structures and utilities classes
COMMONDIR = ${ROOTDIR}/src/common
# DISP Framework src dir
DISPDIR = ${COMMONDIR}/DISP/core

HDF5LIB = ${HDF5_HOME}/lib
HDF5INC = ${HDF5_HOME}/include

MPILIB = ${MPI_HOME}/lib
MPIINC = ${MPI_HOME}/inc

BGPMLIB = ${BGPM_HOME}/lib
BGPMINC = ${BGPM_HOME}/inc

LIBDIRS = -L$(HDF5LIB) -L$(MPILIB) #-L${BGPMLIB}
INCLUDES = -I$(DISPDIR) -I${COMMONDIR} -I$(HDF5INC) -I$(MPIINC) #-I${BGPMINC}

LIBS = -lm -lhdf5 -lz -lmpi -lrt 
#LIBS += -lbgpm

# Executable/reconstruction objects
SIRT_OBJS = sirt.o sirt_main.o #sirt_main_model.o #sirt_main.o
COMMON_OBJS = trace_h5io.o trace_utils.o

# Executables
PROGS = sirt #sirt_model 

.PHONY: default clean

all: $(PROGS)

# SIRT
sirt_model: $(COMMON_OBJS) $(SIRT_OBJS)
	$(CC) $(CFLAGS) -o $@ $(SIRT_OBJS) $(COMMON_OBJS) $(LIBDIRS) $(LIBS) 

sirt: $(COMMON_OBJS) $(SIRT_OBJS)
	$(CC) $(CFLAGS) -o $@ $(SIRT_OBJS) $(COMMON_OBJS) $(LIBDIRS) $(LIBS) 

sirt.o: sirt.cc sirt.h
	$(CC) $(CFLAGS) -c sirt.cc $(INCLUDES)

sirt_main_model.o: sirt_main_model.cc
	$(CC) $(CFLAGS) -c sirt_main_model.cc $(INCLUDES)

sirt_main.o: sirt_main.cc
	$(CC) $(CFLAGS) -c sirt_main.cc $(INCLUDES)

trace_h5io.o: $(COMMONDIR)/trace_h5io.cc $(COMMONDIR)/trace_h5io.h
	$(CC) $(CFLAGS) -c $(COMMONDIR)/trace_h5io.cc $(INCLUDES)

trace_utils.o: $(COMMONDIR)/trace_utils.cc $(COMMONDIR)/trace_utils.h
	$(CC) $(CFLAGS) -c $(COMMONDIR)/trace_utils.cc $(INCLUDES)

install:
	mv $(PROGS) $(BINDIR)

clean:
	rm -f $(PROGS) *.o *.a *~ *.lst *.tmp .pure *.bak *.log

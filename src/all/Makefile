# Compiler
CC=mpiicpc

# Directories
ROOTDIR = ../..
# Bin directory
BINDIR = $(ROOTDIR)/bin/

# Flags
#CFLAGS = -O3 -Wall -Wextra -std=c++11 
# KNL Flags
CFLAGS = -O3 -xMIC-AVX512 -mcmodel medium -shared-intel -std=c++11 
#CFLAGS += -qopt-report=4 -ipo -qopt-report-phase=loop,vec,cg,ipo -qopt-report-file=make-report.txt -qopt-report-embed  # -qopt-report-phase=loop,vec,cg 
#CFLAGS += -guide -guide-file=guide-file.txt
CFLAGS += -g -S -fsource-asm
CFLAGS += -restrict -qopt-prefetch=3 # -qopt-prefetch=5 # Test the prefetcher
CFLAGS += -Wall -Wextra -DTIMERON # -fast # Less precision
CFLAGS += -DT_KNL_OPTIMIZED # Enable KNL optimizations
# -qopt-report=[0-5] : report what optimizations are performed
#  From ALCF presentation:
#   -c -S -xmic-avx512 -O3 -qopt-report=4 -qopt-report-file=stderr
#   -qopt-report-phase=loop,vec,cg -qopt-report-embed

# Common data structures and utilities classes
COMMONDIR = $(ROOTDIR)/src/common
# DISP Framework src dir
DISPDIR = $(COMMONDIR)/disp

HDF5_HOME=${KNL_LOCAL_HOME}
HDF5INC = ${HDF5_HOME}/include
HDF5LIB = ${HDF5_HOME}/lib

LIBS = -lhdf5 -lz # -lpthread
LIBDIRS = -L$(HDF5LIB) 
INCLUDES = -I$(DISPDIR) -I$(COMMONDIR) -I$(HDF5INC) 

# Executable/reconstruction objects
TRACE_OBJS = trace.o
RECON_OBJS =  recon_space.o sirt.o mlem.o pml.o # Reconstruction codes
COMMON_OBJS = trace_utils.o trace_h5io.o trace_runtime_config.o trace_engine.o trace_ph5io.o

# Executables
PROGS = trace

.PHONY: default clean

all: $(PROGS)

# Trace
trace: $(TRACE_OBJS) $(COMMON_OBJS) $(RECON_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TRACE_OBJS) $(RECON_OBJS) $(COMMON_OBJS) $(LIBDIRS) $(LIBS)

trace.o: trace.cc
	$(CC) $(CFLAGS) -c trace.cc $(INCLUDES)

trace_h5io.o: $(COMMONDIR)/trace_h5io.cc $(COMMONDIR)/trace_h5io.h
	$(CC) $(CFLAGS) -c $(COMMONDIR)/trace_h5io.cc $(INCLUDES)

trace_utils.o: $(COMMONDIR)/trace_utils.cc $(COMMONDIR)/trace_utils.h
	$(CC) $(CFLAGS) -c $(COMMONDIR)/trace_utils.cc $(INCLUDES)
trace_runtime_config.o: trace_runtime_config.cc trace_runtime_config.h
	$(CC) $(CFLAGS) -c trace_runtime_config.cc $(INCLUDES)
trace_engine.o: trace_engine.cc trace_engine.h
	$(CC) $(CFLAGS) -c trace_engine.cc $(INCLUDES)
trace_ph5io.o: trace_ph5io.cc trace_ph5io.h
	$(CC) $(CFLAGS) -c trace_ph5io.cc $(INCLUDES)

recon_space.o: recon_space.cc recon_space.h
	$(CC) $(CFLAGS) -c recon_space.cc $(INCLUDES)
sirt.o: sirt.cc sirt.h
	$(CC) $(CFLAGS) -c sirt.cc $(INCLUDES)
mlem.o: mlem.cc mlem.h
	$(CC) $(CFLAGS) -c mlem.cc $(INCLUDES)
pml.o: pml.cc pml.h
	$(CC) $(CFLAGS) -c pml.cc $(INCLUDES)

install:
	mv $(PROGS) $(BINDIR)

clean:
	rm -f $(PROGS) *.s *.o *.a *~ *.lst *.tmp .pure *.bak *.log

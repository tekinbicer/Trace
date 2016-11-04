# Compiler
CC = clang++

# Directories
ROOTDIR = ../..
# Bin directory
BINDIR = $(ROOTDIR)/bin/

# Flags
CFLAGS = -O3 -Wall -Wextra -std=c++11 
CFLAGS += -g -fPIC -DTIMERON

# Common data structures and utilities classes
COMMONDIR = $(ROOTDIR)/src/common
# DISP Framework src dir
DISPDIR = $(COMMONDIR)/DISP/core

HDF5LIB = /usr/lib/x86_64-linux-gnu/hdf5/mpich
HDF5INC = /usr/include/hdf5/mpich
MPIINC = /usr/include/mpich
MPILIB =

LIBS = -lm -lmpi -lhdf5 -lpthread
LIBDIRS = -L$(HDF5LIB)
INCLUDES = -I$(DISPDIR) -I$(COMMONDIR) -I$(HDF5INC) -I$(MPIINC)

# Executable/reconstruction objects
TRACE_OBJS = trace.o
RECON_OBJS =  # Reconstruction codes
COMMON_OBJS = trace_h5io.o trace_utils.o

# Executables
PROGS = trace

.PHONY: default clean

all: $(PROGS)

# Trace
trace: $(TRACE_OBJS) $(COMMON_OBJS) $(RECON_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TRACE_OBJS) $(COMMON_OBJS) $(RECON_OBJS) $(LIBDIRS) $(LIBS)

trace.o: trace.cc
	$(CC) $(CFLAGS) -c trace.cc $(INCLUDES)


trace_h5io.o: $(COMMONDIR)/trace_h5io.cc $(COMMONDIR)/trace_h5io.h
	$(CC) $(CFLAGS) -c $(COMMONDIR)/trace_h5io.cc $(INCLUDES)

trace_utils.o: $(COMMONDIR)/trace_utils.cc $(COMMONDIR)/trace_utils.h
	$(CC) $(CFLAGS) -c $(COMMONDIR)/trace_utils.cc $(INCLUDES)

install:
	mv $(PROGS) $(BINDIR)

clean:
	rm -f $(PROGS) *.o *.a *~ *.lst *.tmp .pure *.bak *.log

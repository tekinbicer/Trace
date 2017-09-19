# Compiler
CC = mpiicpc

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
DISPDIR = $(COMMONDIR)/disp

HDF5_HOME = "/home/bicer/projects/phdf5"
HDF5INC = ${HDF5_HOME}/include
HDF5LIB = ${HDF5_HOME}/lib

LIBS = -lhdf5 #-lmpich # -lpthread
LIBDIRS = -L$(HDF5LIB) #-L/usr/lib64/mpich-3.2/lib
INCLUDES = -I$(DISPDIR) -I$(COMMONDIR) -I$(HDF5INC) #-I/usr/include/mpich-3.2-x86_64

# Executable/reconstruction objects
TRACE_OBJS = trace.o
RECON_OBJS =  recon_space.o sirt.o mlem.o pml.o apmlr.o # Reconstruction codes
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
apmlr.o: apmlr.cc apmlr.h
	$(CC) $(CFLAGS) -c apmlr.cc $(INCLUDES)

install:
	mv $(PROGS) $(BINDIR)

clean:
	rm -f $(PROGS) *.o *.a *~ *.lst *.tmp .pure *.bak *.log

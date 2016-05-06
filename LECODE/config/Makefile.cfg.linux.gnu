# FC MPI compiler to use
FC= mpif90 
# C MPI compiler to use
CC=mpicc 
# C MPI compiler to use
CPP=mpicxx 

# Executable name
EXE=lecode

# Shared library name
SHAREFLAGS=-shared -Wl,-soname,libLECODE.so.4
PLUGLIB= libLECODE.so.4.0.1
PLUGSLN= libLECODE.so.4
SHARELIB= libLECODE.so

# Rules to make library
AR = ar -rcs

# C/C++ optimisation flags
CFLAGS= -O2 -fno-common -fPIC
CPPFLAGS= -O2 -fno-common -fPIC

# Fortran optimisation flags
FCFLAGS=-O3 -funroll-loops --param max-unroll-times=2 -cpp -ftree-vectorize -ffast-math -lstdc++
#FCFLAGS=-O2 -cpp -lstdc++
#FCFLAGS= -O0 -g -fcheck=bounds -Wall -fbacktrace -finit-real=zero\
	-fimplicit-none -ffpe-trap=zero,overflow,invalid -ffree-form -fno-common\
	-Wtabs -Wunused-parameter -Wuninitialized  -ffree-line-length-none -cpp \
	-fdump-fortran-optimized -fdump-tree-original -lstdc++

FOX=/usr/local/gnu/fox/bin

# Linked libraries

TRIANGLELIB=/usr/local/gnu/triangle

HDF5=/usr/local/gnu/hdf5
H5LDFLAGS = -L${HDF5}/lib 
H5FLAGS = -I${HDF5}/include
H5LIBS =  -lhdf5_fortran -lhdf5  -lhdf5hl_fortran -lhdf5_hl -L/usr/local/gnu/zlib/lib -lz -ldl -lstdc++


KDTREE=/usr/local/gnu/kdtree2 
KDTREELDFLAGS = -L${KDTREE}
KDTREEFLAGS = -I${KDTREE}
KDTREELIBS = -lkdtree
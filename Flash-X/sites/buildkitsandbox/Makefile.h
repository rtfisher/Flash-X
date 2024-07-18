# FLASH makefile definitions for x86-64 Linux (GNU compilers)
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

#HDF5_PATH   = /usr/local/hdf5-1.10.10
#HDF5_PATH=/usr/local/hdf5-1.14.4.2
HDF5_PATH_INC=/usr/include/hdf5/openmpi
HDF5_PATH_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial

ZLIB_PATH   =

PAPI_PATH   =
PAPI_FLAGS  =

LIB_NCMPI   = 
MPE_PATH    =
#MPI_PATH    = /usr/local/openmpi-4.1.5/bin/
MPI_PATH=/usr/bin

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------
FCOMP   = $(MPI_PATH)/mpif90
CCOMP   = $(MPI_PATH)/mpicc
CPPCOMP = $(MPI_PATH)/mpiCC
LINK    = $(MPI_PATH)/mpif90

# pre-processor flag
PP      = -D

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for
#  flash_test, and is set for quick code generation, and (sometimes)
#  profiling.  The Makefile generated by setup will assign the generic token
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

FFLAGS_OPT = -ggdb -c -O2 -fdefault-real-8 -fdefault-double-8 \
-Wuninitialized -Wno-argument-mismatch -I${HDF5_PATH_INC} -fcommon \

FFLAGS_DEBUG = -ggdb -c -O0 -fdefault-real-8 -fdefault-double-8 \
-pedantic -Wall -Wextra -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbounds-check \
-fimplicit-none -fstack-protector-all \
-Wno-argument-mismatch -I${HDF5_PATH_INC}

FFLAGS_TEST = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -Wno-argument-mismatch -I${HDF5_PATH_INC}

FFLAGS_HYPRE = -I${HYPRE_PATH}/include


F90FLAGS =

#The macro _FORTIFY_SOURCE adds some lightweight checks for buffer
#overflows at both compile time and run time (only active at -O1 or higher)
#http://gcc.gnu.org/ml/gcc-patches/2004-09/msg02055.html
CFLAGS_OPT = -c -O2 -Wno-implicit-function-declaration

CFLAGS_DEBUG = -c -O0 -Wno-div-by-zero -Wundef \
-Wconversion -Wstrict-prototypes -Wunreachable-code \
-pedantic -Wall -Wextra -Winit-self -ftree-vrp -Wfloat-equal \
-Wunsafe-loop-optimizations -Wpadded -fstack-protector-all \
-Wno-implicit-function-declaration

CFLAGS_TEST = -c

# Platform symbol
CDEFINES += -DDarwin

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I${HDF5_PATH_INC}

CFLAGS_NCMPI = -I$(LIB_NCMPI)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -ggdb -fcommon -Wl,-ld_classic -o 
LFLAGS_DEBUG = -ggdb -O0 -o
LFLAGS_TEST  = -ggdb -o


#----------------------------------------------------------------------------
# Library specific linking
#
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#----------------------------------------------------------------------------

LIB_OPT   =
LIB_DEBUG =
LIB_TEST  =

#LIB_HDF5  = -L${HDF5_PATH_LIB} -lhdf5 /usr/lib64/libz.a
LIB_HDF5  = -L${HDF5_PATH_LIB} -lhdf5 

LIB_PAPI  =
LIB_MATH  =

LIB_MPI   = 
#LIB_NCMPI = -L $(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE   =

LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE

# Uncomment the following line to use electic fence memory debugger.
# Need the following environmental variable (see env.sh):
# export EF_ALLOW_MALLOC_0=1
#CONFIG_LIB = -L/usr/lib64 -lefence

#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
#----------------------------------------------------------------------------

MACHOBJ =

#----------------------------------------------------------------------------
# Additional commands
#----------------------------------------------------------------------------

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo

#----------------------------------------------------------------------------
# Fake existence of iso_c_bindings module to prevent unnecessary recompilations.
#---------------------------------------------------------------------------- 
ifeq ($(FLASHBINARY),true)
iso_c_binding.mod :
	touch $@
endif

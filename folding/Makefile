# -*- Makefile -*-
# contact: Andrea Pozzer MPIC (andrea.pozzer@mpic.de)
# ----------------------------------------------
export
SHELL    = sh

# hydra intel fortran
# F90       = ifort
F90  = gfortran -ffree-line-length-512

#debugging
#F90FLAGS  = -autodouble -cpp -g -debug full -traceback -O3
#production
# F90FLAGS  = -autodouble -cpp -O3
F90FLAGS  = -fdefault-real-8 -cpp -O2 -ffree-line-length-512 -fcheck=all

FC        = $(F90)
FFLAGS    = $(F90FLAGS)

# NETCDF_ROOT = /ptmp/mpcdata/software/x86_64-suse-linux/netcdf/v3.6.3_i
NETCDF_LIB = /usr/lib64
NETCDF_INC = /usr/include
NETCDF_MOD = /usr/lib64/gfortran/modules

# LIBS      = -L$(NETCDF_ROOT)/lib -lnetcdf
# INCLUDES  = -I$(NETCDF_ROOT)/include
LIBS      = -L$(NETCDF_ROOT)/lib -lnetcdff -lnetcdf
INCLUDES  = -I$(NETCDF_INC) -I$(NETCDF_MOD)

INSTALLDIR = .
# ----------------------------------------------
include main.mk
# ----------------------------------------------

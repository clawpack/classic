# ==============================================================================
#  Program:     Advection Test 1D
#  Created:     Date
# ==============================================================================
#  Available rules:
#  Available rules:"
#    all - Build execubtable"
#    advection1d - Build main executable
#    clean - Deletes local objects
#    clobber - Deletes local and library objects
#    help - Prints this help message
# ==============================================================================
# Basic MAKE settings
SHELL = /bin/sh
INSTALL_PROGRAM ?= $(INSTALL)
INSTALL_DATA ?= $(INSTALL) -m 644

CLAW_PKG = classic                  # Clawpack package to use
EXE = advection1d                 	# Executable to create
SETRUN_FILE = setrun.py        		# File containing function to make data
OUTDIR = _output               		# Directory for output
SETPLOT_FILE = setplot.py      		# File containing function to set plots
PLOTDIR = _plots               		# Directory for plots

# ==============================================================================
# Compiler settings and flags
FC ?= gfortran
FFLAGS = -fbounds-check -ffpe-trap=invalid,overflow,zero  -fopenmp -g

# ==============================================================================
# Source lists
LIB_PATH = ../../src

MODULES = \
	$(LIB_PATH)/precision_module.f90 \
	$(LIB_PATH)/clawdata_module.f90 \
	$(LIB_PATH)/1d/solution_module.f90 \
	$(LIB_PATH)/1d/solver_module.f90 \
	$(LIB_PATH)/1d/rp1_interface.f90 \
	./rp1_advection.f90

SOURCES = \
	./main.f90 \
	$(LIB_PATH)/1d/evolve_to_time.f90 \
	$(LIB_PATH)/1d/hyperbolic_step.f90 \
	$(LIB_PATH)/1d/before_step.f90 \
	$(LIB_PATH)/1d/source_term.f90 \
	$(LIB_PATH)/1d/boundary_conditions.f90 \
	$(LIB_PATH)/open_data_file.f90


# ==============================================================================
# Custom Targets

# Generic Include File
include $(CLAWUTIL)/src/Makefile.common

### DO NOT remove this line - make depends on it ###
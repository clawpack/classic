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

CLAW_PKG = clawpack                 # Clawpack package to use
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
LIB_PATH = ../../src/1d

MODULES = \
	$(LIB_PATH)/solution_module.f90 \
	$(LIB_PATH)/solver_module.f90

SOURCES = \
	./main.f90 \
	./rp1.f90 \
	$(LIB_PATH)/step1.f90


# ==============================================================================
# Custom Targets

# Generic Include File
include $(CLAWUTIL)/src/Makefile.common

### DO NOT remove this line - make depends on it ###
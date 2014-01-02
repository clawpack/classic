
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = classic                  # Clawpack package to use
EXE = xclaw                         # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=  

# ---------------------------------
# List of sources for this program:
# ---------------------------------

MODULES = \

SOURCES = \
  qinit.f \
  setprob.f \
  bc3.f \
  setaux.f \
  $(CLAW)/riemann/src/rpn3_vc_acoustics.f90 \
  $(CLAW)/riemann/src/rpt3_vc_acoustics.f90 \
  $(CLAW)/riemann/src/rptt3_vc_acoustics.f90 \
  $(CLAW)/classic/src/3d/driver.f90 \
  $(CLAW)/classic/src/3d/claw3ez.f \
  $(CLAW)/classic/src/3d/b4step3.f90 \
  $(CLAW)/classic/src/3d/chkmth.f \
  $(CLAW)/classic/src/3d/claw3.f \
  $(CLAW)/classic/src/3d/step3.f90 \
  $(CLAW)/classic/src/3d/step3ds.f90 \
  $(CLAW)/classic/src/3d/dimsp3.f \
  $(CLAW)/classic/src/3d/flux3.f90 \
  $(CLAW)/classic/src/3d/copyq3.f \
  $(CLAW)/classic/src/3d/limiter.f90 \
  $(CLAW)/classic/src/3d/philim.f90 \
  $(CLAW)/classic/src/3d/src3.f90 \
  $(CLAW)/classic/src/3d/out3.f \
  $(CLAW)/classic/src/3d/restart3.f \
  $(CLAW)/classic/src/3d/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


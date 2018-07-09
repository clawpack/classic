
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

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=

# ---------------------------------
# package sources for this program:
# ---------------------------------

include $(CLAW)/classic/src/2d/Makefile.classic_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \


# ---------------------------------
# List of custom sources for this program:
# ---------------------------------

MODULES = \

SOURCES = \
  qinit.f \
  setprob.f \
  $(CLAW)/riemann/src/rpn2_acoustics.f90 \
  $(CLAW)/riemann/src/rpt2_acoustics.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

#-------------------------------------------------------------------
# Add make options that use absorbing boundary layer setrun
.data_abl: setrun_abl.py $(MAKEFILE_LIST) ;
	$(MAKE) data_abl -f $(MAKEFILE_LIST)

data_abl: $(MAKEFILE_LIST);
	-rm -f .data_abl
	$(CLAW_PYTHON) setrun_abl.py $(CLAW_PKG)
	touch .data_abl

.output_abl: $(EXE) .data_abl $(MAKEFILE_LIST);
	$(MAKE) output_abl -f $(MAKEFILE_LIST)

output_abl: $(MAKEFILE_LIST);
	-rm -f .output_abl
	$(CLAW_PYTHON) $(CLAW)/clawutil/src/python/clawutil/runclaw.py $(EXE) $(OUTDIR) \
	$(OVERWRITE) $(RESTART) . $(GIT_STATUS) $(NOHUP) $(NICE)
	@echo $(OUTDIR) > .output_abl

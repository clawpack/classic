
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

CLAW_PKG = classic
EXE = xclaw
FFLAGS = -g -fbounds-check -ffpe-trap=zero,overflow,invalid
LFLAGS = -g

SOURCES = \
  qinit.f \
  setprob.f \
  setaux.f \
  $(CLAW)/riemann/src/rp1_acoustics_variable.f90 \
  $(CLAW)/classic/src/1d/b4step1.f \
  $(CLAW)/classic/src/1d/bc1.f \
  $(CLAW)/classic/src/1d/claw1ez.f \
  $(CLAW)/classic/src/1d/claw1.f \
  $(CLAW)/classic/src/1d/copyq1.f \
  $(CLAW)/classic/src/1d/driver.f90 \
  $(CLAW)/classic/src/1d/inlinelimiter.f90 \
  $(CLAW)/classic/src/1d/opendatafile.f \
  $(CLAW)/classic/src/1d/out1.f \
  $(CLAW)/classic/src/1d/src1.f \
  $(CLAW)/classic/src/1d/step1.f90

include $(CLAWMAKE)

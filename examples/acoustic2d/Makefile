
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

CLAW_PKG = classic
EXE = xclaw
FFLAGS = -g -fbounds-check -ffpe-trap=zero,overflow,invalid
LFLAGS = -g

SOURCES = \
  qinit.f \
  setprob.f \
  $(CLAW)/riemann/src/rpn2_acoustics.f90 \
  $(CLAW)/riemann/src/rpt2_acoustics.f90 \
  $(CLAW)/classic/src/2d/setaux.f90 \
  $(CLAW)/classic/src/2d/driver.f90 \
  $(CLAW)/classic/src/2d/claw2ez.f \
  $(CLAW)/classic/src/2d/claw2.f \
  $(CLAW)/classic/src/2d/bc2.f \
  $(CLAW)/classic/src/2d/b4step2.f \
  $(CLAW)/classic/src/2d/step2.f90 \
  $(CLAW)/classic/src/2d/step2ds.f90 \
  $(CLAW)/classic/src/2d/dimsp2.f \
  $(CLAW)/classic/src/2d/flux2.f90 \
  $(CLAW)/classic/src/2d/copyq2.f \
  $(CLAW)/classic/src/2d/inlinelimiter.f90 \
  $(CLAW)/classic/src/2d/src2.f \
  $(CLAW)/classic/src/2d/out2.f \
  $(CLAW)/classic/src/2d/restart2.f \
  $(CLAW)/classic/src/2d/opendatafile.f

include $(CLAWMAKE)


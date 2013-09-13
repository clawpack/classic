
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

CLAW_PKG = classic
EXE = xclaw
FFLAGS = -g -fbounds-check -ffpe-trap=zero,overflow,invalid
LFLAGS = -g

# No modules in 3D Classic (yet)

SOURCES = \
  qinit.f \
  $(CLAW)/riemann/src/rpn3_vc_acoustics.f90 \
  $(CLAW)/riemann/src/rpt3_vc_acoustics.f90 \
  $(CLAW)/riemann/src/rptt3_vc_acoustics.f90 \
  setprob.f \
  bc3.f \
  setaux.f \
  $(CLAW)/classic/src/3d/driver.f90 \
  $(CLAW)/classic/src/3d/claw3ez.f \
  $(CLAW)/classic/src/3d/b4step3.f \
  $(CLAW)/classic/src/3d/chkmth.f \
  $(CLAW)/classic/src/3d/claw3.f \
  $(CLAW)/classic/src/3d/step3.f90 \
  $(CLAW)/classic/src/3d/step3ds.f90 \
  $(CLAW)/classic/src/3d/dimsp3.f \
  $(CLAW)/classic/src/3d/flux3.f90 \
  $(CLAW)/classic/src/3d/copyq3.f \
  $(CLAW)/classic/src/3d/limiter.f90 \
  $(CLAW)/classic/src/3d/philim.f90 \
  $(CLAW)/classic/src/3d/src3.f \
  $(CLAW)/classic/src/3d/out3.f \
  $(CLAW)/classic/src/3d/restart3.f \
  $(CLAW)/classic/src/3d/opendatafile.f

include $(CLAWMAKE)

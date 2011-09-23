EXECUTABLE = xclaw
LIB_PATHS = ./lib ./lib/solution ./lib/solver ./lib/rp

FC = gfortran
LINK = $(FC)
FFLAGS = -g
LFLAGS = $(FFLAGS)
MODULE_PATH = ./

MODULE_SRC = 
SRC = main.f90
LIB_OBJ = \
	./lib/utility_module.o \
	./lib/solution/grid_module.o \
	./lib/solution/state_1d_module.o \
	./lib/solution/solution_1d_module.o \
	./lib/solver/solver_types_module.o \
	./lib/rp/rp_euler_module.o

OBJ = $(subst .f90,.o,$(SRC))
MODULES = $(addprefix $(MODULE_PATH),$(subst .f90,.mod,$(MODULE_SRC)))
INCLUDE = $(sort $(addprefix -I,$(MODULE_PATH)) $(addprefix -I,$(LIB_PATHS)))

%.mod : %.f90
	$(FC) -c $(INCLUDE) -J$(MODULE_PATH) $(FFLAGS) $< -o $(subst .mod,.o,$@)
%.o : %.f90
	$(FC) -c $(INCLUDE) $(FFLAGS) $< -o $@
%.o : %.f
	$(FC) -c $(INCLUDE) $(FFLAGS) $< -o $@
	
$(EXECUTABLE): lib $(MODULES) $(MODULE_OBJ) $(OBJ)
	$(LINK) $(LIB_OBJ) $(MODULE_OBJ) $(OBJ) $(INCLUDE) $(LFLAGS) -o $(EXECUTABLE)
	
lib: 
	$(MAKE) -C./lib lib
	
clean:
	-rm $(MODULES) $(MODULE_OBJ) $(OBJ) $(EXECUTABLE)
	
clobber:
	$(MAKE) clean
	$(MAKE) -C./lib clean
	

EXECUTABLE = xclaw
LIB_PATH = ./lib

FC = gfortran
LINK = $(FC)
FFLAGS = -g
LFLAGS = $(FFLAGS)

SRC = main.f90
MOD_SRC = \
	$(LIB_PATH)/utility_module.f90 \
	$(LIB_PATH)/solution/grid_module.f90 \
	$(LIB_PATH)/solution/state_1d_module.f90 \
	$(LIB_PATH)/solution/solution_1d_module.f90 \
	$(LIB_PATH)/solver/solver_types_module.f90 \
	$(LIB_PATH)/rp/rp_euler_module.f90 \
	$(LIB_PATH)/rp/rp_swe_module.f90   

OBJ = $(subst .f90,.o,$(SRC))
MOD_FILES = $(subst .f90,.mod,$(MOD_SRC))
MOD_OBJ = $(subst .f90,.o,$(MOD_SRC))
INCLUDE = $(sort $(addprefix -I,$(dir $(MOD_SRC))))

%.mod : %.f90
	$(FC) -c $(INCLUDE) -J$(dir $<) $(FFLAGS) $< -o $(subst .mod,.o,$@)
%.o : %.f90
	$(FC) -c $(INCLUDE) $(FFLAGS) $< -o $@
%.o : %.f
	$(FC) -c $(INCLUDE) $(FFLAGS) $< -o $@

$(EXECUTABLE): $(MOD_FILES) $(MOD_OBJ) $(OBJ)
	$(LINK) $(MOD_OBJ) $(OBJ) $(INCLUDE) $(LFLAGS) -o $(EXECUTABLE)
	
clean:
	-rm $(OBJ) $(EXECUTABLE)
	
clobber: clean
	-rm $(MOD_OBJ) $(MOD_FILES)


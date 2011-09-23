FC ?= gfortran
FFLAGS = -g
LFLAGS = $(FFLAGS)

MODULE_SRC = \
	./utility_module.f90 \
	./solution/grid_module.f90 \
	./solution/state_1d_module.f90 \
	./solution/solution_1d_module.f90 \
	./solver/solver_types_module.f90 \
	./rp/rp_euler_module.f90 \
	./rp/rp_swe_module.f90
	
SRC = 
	
OBJ = $(subst .f90,.o,$(SRC))
MODULE_OBJ = $(subst .f90,.o,$(MODULE_SRC))
MODULES = $(subst .f90,.mod,$(MODULE_SRC))
INCLUDE = $(sort $(addprefix -I,$(dir $(MODULES))))

%.mod : %.f90
	$(FC) -c $(INCLUDE) -J$(dir $<) $(FFLAGS) $< -o $(subst .mod,.o,$@)
%.o : %.f90
	$(FC) -c $(INCLUDE) $(FFLAGS) $< -o $@
%.o : %.f
	$(FC) -c $(INCLUDE) $(FFLAGS) $< -o $@
	
lib: $(MODULES) $(OBJ) 
	
clean:
	-rm $(MODULES) $(OBJ)
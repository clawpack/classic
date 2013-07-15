set(SRC
    clawdata_module.f90  open_data_file.f90  precision_module.f90
    utils.f90

    1d/before_step.f90          1d/hyperbolic_step.f90  1d/solution_module.f90
    1d/boundary_conditions.f90  1d/limiter.f90          1d/solver_module.f90
    1d/controller_module.f90    1d/qinit.f90            1d/source_term.f90
    1d/evolve_to_time.f90       1d/rp1_interface.f90
    1d/geometry_module.f90      1d/set_aux.f90
)

add_library(classic ${SRC})

install(TARGETS classic
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
    )

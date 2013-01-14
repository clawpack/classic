program advection1d

    use solver_module
    use solution_module

    use clawdata_module

    implicit none
    
    type(solution_type) :: solution
    type(solver_type) :: solver
    type(clawdata_type) :: clawdata

    ! Read input data
    clawdata = read_data_file('./claw.data')

    solution = new_solution(clawdata)
    solver = new_solver(clawdata)

    call evolve_to_time(clawdata%t_final,solution,solver)
 
end program advection1d

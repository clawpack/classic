program advection1d

    use solver_module
    use solution_module

    implicit none
    
    type(solution_type) :: solution
    type(solver_type) :: solver

    solution = new_solution(100,2,0,2)

    solution%q = 10.d2

    print "('q = ',2f9.1)", solution%q(1:2,12)

    solver = new_solver()

    call step1(solution,solver)
 
end program advection1d

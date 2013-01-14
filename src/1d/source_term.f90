subroutine source_term(t, dt, solution, solver)
    
    use solution_module, only: solution_type
    use solver_module, only: solver_type

    implicit none

    real(kind=8), intent(in) :: t, dt
    type(solution_type), intent(in out) :: solution
    type(solver_type), intent(in out) :: solver
    
end subroutine source_term
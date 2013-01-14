subroutine before_step(solution,solver)

    use solution_module, only: solution_type
    use solver_module, only: solver_type

    implicit none

    type(solution_type), intent(in out) :: solution
    type(solver_type), intent(in out) :: solver

end subroutine before_step
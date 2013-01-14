! Default qinit function
subroutine qinit(solution)

    use solution_module, only: solution_type

    implicit none

    type(solution_type), intent(in out) :: solution

end subroutine qinit
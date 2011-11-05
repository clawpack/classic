module problem_module

    use solution_module_1d
    
    integer, parameter :: gamma = 1.4d0
    
contains

    subroutine setup(grid,state,rp_parms,rp_work)

        implicit none
        type(grid_type), pointer, intent(in) :: grid,state,rp_parms,rp_work
        type(state_type), pointer, intent(inout) :: state
        type(rp_params), pointer, intent(inout) :: rp_params
        type(rp_work), pointer, intent(inout) :: rp_work

    end subroutine setup

    subroutine qinit(grid,state)

        implicit none
        type(grid_type), intent(in) :: grid
        type(state_type), intent(inout) :: state

        state%q(1,:) = 1.d0
        state%q(2,:) = 2.d0

    end subroutine qinit
    
    subroutine set_aux(grid,state)

        implicit none
        type(grid_type), intent(in) :: grid
        type(state_type), intent(inout) :: state

        state%aux(1,:) = 1.d0 - gamma
        state%aux(2,:) = gamma

    end subroutine set_aux
    
end module problem module
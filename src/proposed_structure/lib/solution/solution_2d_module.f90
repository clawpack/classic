module solution_1d_module

    use grid_module
    use state_1d_module

    implicit none
    
    integer, parameter, private :: MAX_GRIDS = 1000
    integer, parameter, private :: MAX_STATES = 15
    
    type solution_type
        type(grid_type_pointer), private :: grids(MAX_GRIDS)
        type(state_type_pointer), private :: states(MAX_GRIDS*MAX_STATES) 
        logical, private :: available_grids(MAX_GRIDS)
        logical, private :: available_states(MAX_GRIDS*MAX_STATES)
    end type solution_type
    
    interface add
        module procedure add_grid
        module procedure add_state
    end interface add
    
!     interface delete
!         module procedure delete_grid
!         module procedure delete_state
!     end interface delete
    
contains

    subroutine add_grid(self,grid)

        implicit none
        
        ! Inputs
        type(solution_type), intent(inout) :: self
        type(grid_type), intent(in) :: grid
        
        ! Locals
        integer :: i
        do i=1,MAX_GRIDS
            if (self%available_grids(i)) then
                self%grids(i)%grid => grid
                self%available_grids(i) = .false.
                exit
            endif
        enddo
        if (i > MAX_GRIDS) then
            print *,"ERROR:  Reached maximum allowed number of grids "
            print *,"        increase MAX_GRIDS parameter in "
            print *,"        solution_2d_mod.f90 to fix this problem!"
            stop
        endif
        
    end subroutine add_grid
    
    subroutine add_state(self,state)

        implicit none
        type(solution_type), intent(inout) :: self
        type(state_1d_type_pointer), intent(in) :: state
        
        integer :: i
        do i=1,MAX_GRIDS*MAX_STATES
            if (self%available_states(i)) then
                self%states(i)%state => state
                self%available_states(i) = .false.
                exit
            endif
        enddo
        if (i > MAX_GRIDS*MAX_STATES) then
            print *,"ERROR:  Reached maximum allowed number of states "
            print *,"        increase MAX_STATES parameter in "
            print *,"        solution_2d_mod.f90 to fix this problem!"
            stop
        endif
        
    end subroutine add_state
    
end module solution_2d_module
module solution_1d_module

    use grid_module
    use state_1d_module

    implicit none
    
    integer, parameter, private :: MAX_GRIDS = 1000
    integer, parameter, private :: MAX_STATES = 15
    
    type solution_type
        type(grid_type_pointer) :: grids(MAX_GRIDS)
        type(state_type_pointer) :: states(MAX_GRIDS*MAX_STATES)
        logical, private :: available_grids(MAX_GRIDS) = .true.
        logical, private :: available_states(MAX_GRIDS*MAX_STATES) = .true.
    end type solution_type
    
    interface add
        module procedure add_grid_pointer
        module procedure add_grid
        module procedure add_state_pointer
        module procedure add_state
    end interface add
    
contains

    subroutine add_grid(self,grid)

        implicit none
        
        ! Inputs
        type(solution_type), intent(inout) :: self
        type(grid_type), pointer, intent(in) :: grid
        
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
            print *,"        solution_1d_mod.f90 to fix this problem!"
            stop
        endif
        
    end subroutine add_grid
    
    subroutine add_grid_pointer(self,grid_pointer)

        implicit none
        type(solution_type), intent(inout) :: self
        type(grid_type_pointer), intent(in) :: grid_pointer
        
        ! Locals
        integer :: i
        
        do i=1,MAX_GRIDS
            if (self%available_grids(i)) then
                if (.not.associated(self%grids(i)%grid)) then
                    deallocate(self%grids(i)%grid)
                endif
                self%grids(i) = grid_pointer
                self%available_grids(i) = .false.
                exit
            endif
        enddo
        if (i > MAX_GRIDS) then
            print *,"ERROR:  Reached maximum allowed number of grids "
            print *,"        increase MAX_GRIDS parameter in "
            print *,"        solution_1d_mod.f90 to fix this problem!"
            stop
        endif
        
    end subroutine add_grid_pointer
    
    subroutine add_state(self,state)

        implicit none
        type(solution_type), intent(inout) :: self
        type(state_type), pointer, intent(in) :: state
        
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
            print *,"        solution_1d_mod.f90 to fix this problem!"
            stop
        endif
        
    end subroutine add_state
    
    subroutine add_state_pointer(self,state_pointer)

        implicit none
        type(solution_type), intent(inout) :: self
        type(state_type_pointer), intent(in) :: state_pointer
        
        integer :: i
        do i=1,MAX_GRIDS*MAX_STATES
            if (self%available_states(i)) then
                if (.not.associated(self%states(i)%state)) then
                    deallocate(self%states(i)%state)
                endif
                self%states(i) = state_pointer
                self%available_states(i) = .false.
                exit
            endif
        enddo
        if (i > MAX_GRIDS*MAX_STATES) then
            print *,"ERROR:  Reached maximum allowed number of states "
            print *,"        increase MAX_STATES parameter in "
            print *,"        solution_1d_mod.f90 to fix this problem!"
            stop
        endif
        
    end subroutine add_state_pointer
    
    function get_state(self,state_num) result(state)

        implicit none
        type(solution_type), intent(in) :: self
        integer, intent(in) :: state_num
        type(state_type), pointer :: state
        
        state => self%states(state_num)%state

    end function get_state

end module solution_1d_module
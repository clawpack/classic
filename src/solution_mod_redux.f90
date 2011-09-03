! Module containing pertinent data structures for Classic Clawpack
module solution_module

    implicit none
    
    ! Probably could make this dynamic, later enhancement...
    integer, private, parameter :: MAX_GRIDS = 1000
    integer, private, parameter :: MAX_STAGES = 10
    
    type grid_type
        integer, allocatable :: grid_no
        integer, allocatable :: level
        integer, allocatable :: n(:)
        double precision, allocatable :: d(:),lower(:),upper(:)
    end type grid_type
    
    type state_type
        double precision, allocatable :: q(:,:)
        double precision, allocatable :: aux(:,:)
    end type state_type
    
    type, private :: state_container
        type(state_type), pointer :: state
    end type state_container
    
    type, private :: grid_container
        type(grid_type), pointer :: grid
    end type grid_container
    
    type solution_type
        private
        logical :: available_grids(MAX_GRIDS)
        logical :: available_states(MAX_STAGES*MAX_GRIDS)
        type(grid_container) :: grids(MAX_GRIDS)
        type(state_container) :: states(MAX_STAGES*MAX_GRIDS)
    end type solution_type
    
    interface new
        module procedure new_solution
        module procedure new_grid
        module procedure new_state
    end interface new
    
    interface delete
        module procedure delete_grid
        module procedure delete_state
    end interface delete
    
    interface operator(.eq.)
        module procedure grids_equal
        module procedure states_equal
    end interface operator(.eq.)
    
contains

    function new_solution()
        implicit none
        
        type(solution_type), pointer :: new_solution

        allocate(new_solution)
        new_solution%available_grids = .true.
        new_solution%available_states = .true.

    end function new_solution
    
    subroutine add_grid(self,grid)

        implicit none
        
        ! Input
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
        if (i == MAX_GRIDS) then
            print *,"ERROR:  Reached maximum allowed number of grids "
            print *,"        reached, increase MAX_GRIDS parameter in "
            print *,"        solution_mod.f90 to fix this problem!"
            stop
        endif

    end subroutine add_grid
    
    ! Assume grid_no is unique???
    subroutine delete_grid(self,grid)
        
        implicit none
        
        ! Input
        type(solution_type), intent(inout) :: self
        type(grid_type), intent(in) :: grid
        
        ! Locals
        integer :: i
        
        do i=1,MAX_GRIDS
            if (.not.self%available_grids(i)) then
                if (self%grids(i)%grid == grid) then
                    deallocate(self%grids(i)%grid)
                    self%available_grids(i) = .false.
                    exit
                endif
            endif
        enddo
        if (i == MAX_GRIDS) then
            print *,"WARNING:  Did not find specified grid for deletion!"
        endif
        
    end subroutine delete_grid
    
    subroutine add_state(self,state)

        implicit none
        type(solution_type), intent(inout) :: self
        type(state_type), pointer, intent(in) :: state
        
        ! Locals
        integer :: i
        
        do i=1,MAX_GRIDS*MAX_STAGES
            if (self%available_states(i)) then
                self%states(i)%state => state
                self%available_states(i) = .false.
                exit
            endif
        enddo
        if (i == MAX_GRIDS*MAX_STAGES) then
            print *,"ERROR:  Reached maximum allowed number of states, "
            print *,"        increase MAX_STAGES parameter in "
            print *,"        solution_mod.f90 to fix this problem!"
            stop
        endif
        

    end subroutine add_state

    function new_grid(ndim,n)
        
        implicit none
        
        integer, intent(in) :: ndim,n(:)
        
        type(grid_type), pointer :: new_grid
        
        allocate(new_grid)
        allocate(new_grid%n(ndim),new_grid%d(ndim))
        allocate(new_grid%lower(ndim),new_grid%upper(ndim))
        new_grid%n = n
        
    end function new_grid
    
    function new_state(grid,meqn,maux)

        implicit none
        type(grid_type), intent(in) :: grid
        integer, intent(in) :: meqn,maux
        
        type(state_type), pointer :: new_state
        
        allocate(new_state)
        allocate(new_state%q(meqn,1:grid%n(1)))
        allocate(new_state%aux(maux,1:grid%n(1)))

    end function new_state
    
    logical function grids_equal(grid_1,grid_2)

        implicit none
        type(grid_type), intent(in) :: grid_1,grid_2
        
        grids_equal = (grid_1%grid_no == grid_2%grid_no)
        
    end function grids_equal

end module solution_module


program alloc_test

    use solution_module

    implicit none
    
    type(solution_type), pointer :: solution
    type(grid_type), pointer :: grid
    type(state_type), pointer :: state
    
    solution => new_solution()
    grid => new_grid(1,[10])
    state => new_state(grid,2,0)
    call add_grid(solution,grid)
    call add_state(solution,state)
    
end program alloc_test

module solution_module

    implicit none
    save
    
    ! Global parameter values for solutions, these are kept private
    integer, private, parameter :: MAX_GRIDS = 1000
    integer, private, parameter :: MAX_STAGES = 50
    
    ! Function specifications
    interface init
        module procedure init_solution
        module procedure init_grid
        module procedure init_state
    end interface init
    interface add
        module procedure add_grid
        module procedure add_state
    end interface add
    interface delete
        module procedure delete_grid
        module procedure delete_state
    end interface
    interface get
        module procedure get_grid
        module procedure get_state
    end interface get
!     interface copy
!         module procedure copy_solution
!         module procedure copy_grid
!     end interface copy
    interface operator(.eq.)
        module procedure grids_equal
    end interface operator(.eq.)

    ! Type definitions
    type grid_type
        integer, allocatable :: grid_no
        integer, allocatable :: level
        integer, allocatable :: n(:)
        double precision, allocatable :: d(:),lower(:),upper(:)
    end type grid_type
    type, private :: grid_container
        type(grid_type), pointer :: grid
    end type
    
    type state_type
        type(grid_type), pointer :: grid
        double precision, allocatable :: t
        double precision, allocatable :: q(:,:)
        double precision, allocatable :: aux(:,:)
    end type state_type
    type, private :: state_container
        type(state_type), pointer :: state
    end type state_container
        
    type solution_type    
        private
        logical :: available_grids(MAX_GRIDS)
        logical :: available_states(MAX_STAGES*MAX_GRIDS)
        type(grid_container) :: grids(MAX_GRIDS)
        type(state_container) :: states(MAX_STAGES*MAX_GRIDS)
    end type solution_type
    
contains

    subroutine init_solution(self)

        implicit none
        type(solution_type), pointer, intent(out) :: self
        
        if (.not.associated(self)) then
            allocate(self)
        endif
        self%available_grids = .true.
        self%available_states = .true.
!         self%grids(:)%grid => null()
!         self%states(:)%state => null()
        
    end subroutine init_solution
    
    subroutine init_grid(self,ndim,n)

        implicit none
        type(grid_type), pointer, intent(out) :: self
        integer, intent(in) :: ndim, n(:)
        
        if (.not.associated(self)) then
            allocate(self)
        endif
        
        allocate(self%n(ndim),self%d(ndim))
        allocate(self%lower(ndim),self%upper(ndim))
        self%n = n

    end subroutine init_grid
    
    subroutine init_state(self,grid,meqn,maux)

        implicit none
        type(state_type), pointer, intent(out) :: self
        type(grid_type), pointer, intent(in) :: grid
        integer, intent(in) :: meqn,maux
        
        if (.not.associated(self)) then
            allocate(self)
        endif
        
        allocate(self%q(meqn,1:grid%n(1)))
        allocate(self%aux(maux,1:grid%n(1)))

        self%grid => grid

    end subroutine init_state
    
    subroutine add_grid(self,grid)

        implicit none
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
        ! Note that i is incremented to MAX_GRIDS+1 before exiting loop!
        if (i > MAX_GRIDS) then
            print *,"ERROR:  Reached maximum allowed number of grids "
            print *,"        increase MAX_GRIDS parameter in "
            print *,"        solution_mod.f90 to fix this problem!"
            stop
        endif

    end subroutine add_grid
    
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
        if (i > MAX_GRIDS*MAX_STAGES) then
            print *,"ERROR:  Reached maximum allowed number of states, "
            print *,"        increase MAX_STAGES parameter in "
            print *,"        solution_mod.f90 to fix this problem!"
            stop
        endif
        

    end subroutine add_state
    
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
        if (i > MAX_GRIDS) then
            print *,"WARNING:  Did not find specified grid for deletion!"
        endif
        
    end subroutine delete_grid
    
    subroutine delete_state(self,state)
        
        implicit none
        
        ! Input
        type(solution_type), intent(inout) :: self
        type(state_type), intent(in) :: state
        
        ! Locals
!         integer :: i
        
!         do i=1,MAX_GRIDS*MAX_STAGES
!             if (.not.self%available_states(i)) then
!                 if (self%states(i)%state == state) then
!                     deallocate(self%states(i)%state)
!                     self%available_states(i) = .false.
!                     exit
!                 endif
!             endif
!         enddo
!         if (i == MAX_GRIDS*MAX_STAGES) then
!             print *,"WARNING:  Did not find specified state for deletion!"
!         endif
        print *,"WARNING:  Routine has not been implemented!"
        
    end subroutine delete_state
    
    logical function grids_equal(grid_1,grid_2)
        implicit none
        type(grid_type), intent(in) :: grid_1, grid_2
        
        grids_equal = (grid_1%grid_no == grid_2%grid_no)
    end function grids_equal
    
    logical function states_equal(state_1,state_2)
        implicit none
        type(state_type), intent(in) :: state_1, state_2
        
        states_equal = all(state_1%q == state_2%q)
    end function states_equal

    subroutine get_grid(self,grid_no,grid)
        implicit none
        integer, intent(in) :: grid_no
        type(solution_type), intent(in) :: self 
        type(grid_type), pointer, intent(out) :: grid

        grid = self%grids(grid_no)%grid
    end subroutine get_grid

    subroutine get_state(self,state_no,state)
        implicit none
        type(solution_type), intent(in) :: self
        integer, intent(in) :: state_no
        type(state_type), intent(out), pointer :: state

        state = self%states(state_no)%state
        
    end subroutine get_state

end module solution_module

program solution_test

    use solution_module

    implicit none
    
    type(solution_type), pointer :: solution
    type(grid_type), pointer :: grid
    type(state_type), pointer :: state
    
    call init(solution)
    call init(grid,1,[10])
    call init(state,grid,1,0)

    call add(solution,grid)
    call add(solution,state)

    call get(solution,1,grid)
    call get(solution,1,state)

    print *, grid%n(1)
    print *, state%grid%n

    pause 10

end program solution_test


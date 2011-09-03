module solution_module

    implicit none
    save 
    
    ! Function visibility
!     public new,delete,copy
!     private new_solution,delete_solution,copy_solution
!     private new_grid,delete_grid,copy_grid
!     private new_state,delete_state,copy_state
    
    ! Type definitions
    type grid_type
        
        integer, allocatable :: n(:)
        double precision, allocatable :: d(:),lower(:),upper(:)
        
    end type grid_type
    
    type state_type
        private
        type(grid_type), pointer :: grid
        double precision :: t
        
        double precision, allocatable :: q_work(:)
        double precision, allocatable :: aux_work(:)
    end type state_type
        
    type solution_type    
        private
        
        type(grid_type), pointer :: grids(:)
        type(state_type), pointer :: states(:)
        
        integer :: last_grid, last_state
        
    end type solution_type
    
    interface new
        module procedure new_solution
        module procedure new_grid
    end interface new
    
    interface delete
        module procedure delete_solution
        module procedure delete_grid
    end interface delete
    
    interface copy
        module procedure copy_solution
        module procedure copy_grid
    end interface copy
    
    ! Global parameter values for solutions, these are kept private
    integer, private, parameter :: MAX_GRIDS = 1000
    integer, private, parameter :: MAX_LEVELS = 10
    integer, private, parameter :: MAX_STATES = 1000
    
contains

    subroutine new_solution(self)

        implicit none
        type(solution_type), intent(out) :: self
        
        ! Initialize positional array values
        self%last_grid = 0
        self%last_state = 0
        
    end subroutine new_solution
    
    subroutine add_grid(self,ndim,n,d,lower,upper)
        implicit none
        
        ! This solution object
        type(solution_type), intent(in) :: self
        
        ! Input arguments for construction of grid
        integer, intent(in) :: ndim
        integer, optional, intent(in) :: n(:)
        double precision, optional, intent(in) :: d(:)
        double precision, optional, intent(in) :: lower(:)
        double precision, optional, intent(in) :: upper(:)
        
        ! Input checks
        ! Make sure we have not surpassed the end of our pointer list
        if (last_grid + 1 > MAX_GRIDS) then
            print *,"Maximum number of grids", MAX_GRIDS," reached."
            stop
        endif
        ! Check inputs for correct sizes
        if (present(n)) then
            if (len(n) /= ndim) then
                print *,"Incorrect length of input parameter n,"
                print *,"found len(n) == ",len(n),", expected ",ndim,"."
                stop
            endif
        endif
        if (present(d)) then
            if (len(d) /= ndim) then
                print *,"Incorrect length of input parameter d,"
                print *,"found len(d) == ",len(d),", expected ",ndim,"."
                stop
            endif
        endif
        if (present(lower)) then
            if (len(lower) /= ndim) then
                print *,"Incorrect length of input parameter lower,"
                print *,"found len(lower) == ",len(lower),", expected ",ndim,"."
                stop
            endif
        endif
        if (present(upper)) then
            if (len(upper) /= ndim) then
                print *,"Incorrect length of input parameter upper,"
                print *,"found len(upper) == ",len(upper),", expected ",ndim,"."
                stop
            endif
        endif

        ! Increment last_grid value to point to new grid
        last_grid = last_grid + 1
        
        ! Allocate dimensional arrays
        allocate(self%grids(last_grid)%n(ndim))
        allocate(self%grids(last_grid)%d(ndim))
        allocate(self%grids(last_grid)%lower(ndim))
        allocate(self%grids(last_grid)%upper(ndim))
        
        ! Depending on the input given, assign grid values
        if (present(n).and.present(lower).and.present(upper)) then
            self%grids(last_grid)%n = n
            self%grids(last_grid)%lower = lower
            self%grids(last_grid)%upper = upper
            forall (i=1:ndim)
                self%grids(last_grid)%d = (upper(i) - lower(i)) / n(i)
            end forall
        else
            print *, "Invalid set of parameters given for grid creation."
            stop
        endif
        
    end subroutine add_grid
    
    subroutine delete_grid(self)

        implicit none
        argument type, intent(inout) :: self
        

    end subroutine delete_grid
    
    subroutine add_state(self,grid)
        implicit none
        
        type(solution_type), intent(in) :: self
        type(grid_type), intent(in) :: grid
        
        print *,"Not implemented!"
        stop
        
    end subroutine add_state

end module solution_module

module solution_module

    implicit none
    
    integer, private, parameter :: MAX_GRIDS=1000
    integer, private, parameter :: MAX_STATES=20

    type grid_type
        integer, allocatable :: n(:)
        double precision, allocatable :: d(:), lower(:), upper(:)
    end type grid_type

    type state_type
        integer :: meqn,maux
        double precision :: t
        type(grid_type), pointer :: grid
        double precision, allocatable :: q(:,:)
        double precision, allocatable :: aux(:,:)
    end type state_type

    type solution_type
        type(grid_type) :: grids(MAX_GRIDS)
        type(state_type) :: states(MAX_STATES)
    end type solution_type

contains

    function new_grid(ndim,n)
        implicit none

        integer, intent(in) :: ndim, n(:)
        type(grid_type), pointer :: new_grid

        allocate(new_grid)

        ! Allocate grid
        allocate(new_grid%n(ndim), new_grid%d(ndim), new_grid%lower(ndim), &
            new_grid%upper(ndim))
        new_grid%n = n

    end function new_grid

    function new_state(meqn,maux,grid)
        
        implicit none
        
        integer, intent(in) :: meqn,maux
        type(grid_type), pointer, intent(in) :: grid
        type(state_type), pointer :: new_state

        integer :: err

        allocate(new_state)

        ! Set dimension variables
        new_state%meqn = meqn
        new_state%maux = maux

        ! Allocate state arrays
        allocate(new_state%q(meqn,grid.n(1)), stat=err)
        if (err /= 0) print *, "q: Allocation request denied"
        allocate(new_state%aux(maux,grid.n(1)), stat=err)
        if (err /= 0) print *, "aux: Allocation request denied"
        
        new_state%grid => grid

    end function new_state

end module solution_module

program solution_test
    
    integer, parameter :: ndim = 1

    type(solution_type) :: my_solution
    type(grid_type), pointer :: my_grid
    type(state_type), pointer :: my_state_1
    type(state_type), pointer :: my_state_2


    print *, "Starting Test:"

    print *,"  Setup grid..."
    ! Setup grid
    my_grid => new_grid(1,[10])

    print *,"  Setup states..."
    ! Setup states
    my_state_1 => new_state(2,1,my_grid)
    my_state_1%q(1,:) = 1.d0
    my_state_1%q(2,:) = 2.d0
    my_state_1%aux(1,:) = 4.d0
    my_state_2 => new_state(1,0,my_grid)
    my_state_2%q = 3.d0

    print *,"  Setup solution..."
    ! Setup Solution
    my_solution%states(1) = my_state_1
    my_solution%states(2) = my_state_2
    my_solution%grids(1) = my_grid

    print *,"Output:"
    print *,"n=",my_solution%grids(1)%n,my_solution%states(1)%grid%n,my_solution%states(2)%grid%n
    print *,my_solution%states(1)%q(2,:)
    print *,"============================================================================"
    print *,"q(1)=",my_solution%states(1)%q(:,5)
    print *,"============================================================================"
    print *,"aux(1)=",my_solution%states(1)%aux
    print *,"============================================================================"
    print *,"q(2)=",my_solution%states(2)%q
    print *,"============================================================================"
    print *,"Aux(2)=",my_solution%states(2)%aux

end program solution_test
>>>>>>> b6cd642c481f67e5e1e4b8d1419fe2d4ac3f01a9

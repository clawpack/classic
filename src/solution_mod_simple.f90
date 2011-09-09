module solution_module

    implicit none
    
    ! Array size parameters
    integer, private, parameter :: MAX_GRIDS = 1
    integer, private, parameter :: MAX_STAGES = 15
    
    ! Grid data
    integer, allocatable :: n(:)
    double precision, allocatable :: d(:),lower(:),upper(:)
    
    ! State data (Assumes we only have one state)
    integer :: meqn,maux
    double precision :: t
    double precision, private, allocatable :: q_work(:,:,:,:)
    double precision, private, allocatable :: aux_work(:,:,:,:)
    
contains

    subroutine init_grid(num_dim,new_n,new_d,new_lower,new_upper)

        implicit none
        
        ! Input parameters
        integer, intent(in) :: num_dim
        integer, intent(in), optional :: new_n(:)
        double precision, intent(in), optional :: new_d(:),new_lower(:),new_upper(:)
        
        ! Allocate all the arrays base on number of dimensions
        allocate(n(num_dim),d(num_dim),lower(num_dim),upper(num_dim))
        
        if (present(new_n)) then
            n = new_n
        endif

    end subroutine init_grid
    
    subroutine init_state(meqn,maux)

        implicit none
        integer, intent(in) :: meqn,maux
        
        ! Allocate state arrays
        allocate(q_work(meqn*product(n)),aux_work(maux*product(n)))

    end subroutine init_state

    function get_q_1d()
    
        implicit none
        
        double precision, pointer :: get_q_1d(:,:)
        
        get_q_1d => q_work
        
    end function get_q_1d
    
end module solution_module

program solution_module_test

    use solution_module

    implicit none
    
    call init_grid(1,[10])
    call init_state(2,1)
    
    q = get_q(ndim)
    aux = get_aux(ndim)
    
end program solution_module_test

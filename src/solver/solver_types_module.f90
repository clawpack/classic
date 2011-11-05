module solver_types_module

    ! Type returned by evolve_to_time giving the status of the previous
    ! series of time steps
    type solver_status_type
        ! Number of time steps taken
        integer :: steps_taken
        ! Largest and smallest time step taken
        double precision :: dt_max,dt_min
        ! Last CFL observed, maximum CFL reached
        double precision :: cfl,cfl_max
    end type solver_status_type

    ! Passed to evolve_to_time containing all the solver parameters
    type solver_parameters_type
        ! Number of waves used in the solver, number of ghost cells
        integer :: mwaves, mbc
        
        ! Limiter type for each wave, boundary condtions requested
        integer, allocatable :: limiters(:), bc_lower(:), bc_upper(:)
        
        ! Maximum number of time steps allowed
        integer :: max_steps
        
        ! Initial time step, maximum time step allowed
        double precision :: dt_initial,dt_max
        
        ! Maximum allowed CFL, desired CFL
        double precision :: cfl_max,cfl_desired

        ! Whether to take variable time steps
        logical :: dt_variable
    
        ! Order of method
        integer :: order,trans_order

        ! Source term splitting order
        integer :: src_order

        ! Capacity array index in aux array
        integer :: capa

        ! Verbosity control
        integer :: verbosity
    end type solver_parameters_type

contains
    
    function new_solver_status() result(status)

        implicit none
        type(solver_status_type) :: status

        status%steps_taken = 0
        status%dt_max = tiny(1.d0)
        status%dt_min = huge(1.d0)
        status%cfl = 0.d0
        status%cfl_max = tiny(1.d0)

    end function new_solver_status

    function new_solver_params(file_path,ndim) result(params)
        implicit none
        
        character(len=*), intent(in) :: file_path
        integer, intent(in) :: ndim

        type(solver_parameters_type), pointer :: params
        
        allocate(params)
        
        ! Read in data from file at file_path
        open(unit=13, file=file_path, status="old", action="read")
        
        read(13,"(i1)") params%mwaves
        allocate(params%limiters(params%mwaves))
        read(13,*) params%limiters

        read(13,"(i1)") params%mbc
        allocate(params%bc_lower(ndim),params%bc_upper(ndim))
        read(13,*) params%bc_lower
        read(13,*) params%bc_upper

        close(13)

    end function new_solver_params

end module solver_types_module

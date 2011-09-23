module solver_types_module

    type solver_status_type
        double precision :: dt_max,dt_min
        integer :: max_steps
    end type solver_status_type

    type solver_parameters_type
        integer :: mwaves, mbc
        integer, allocatable :: limiters(:), bc_lower(:), bc_upper(:)
    end type solver_parameters_type

contains
    
    function new_solver_status() result(status)

        implicit none
        type(solver_status_type), pointer :: status

        allocate(status)
        status%dt_max = tiny(1.d0)
        status%dt_min = huge(1.d0)
        status%max_steps = 0

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
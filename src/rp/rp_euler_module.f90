module rp_euler_module

    implicit none
    
    type rp_params_type
        double precision :: gamma
    end type rp_params_type
    
    type rp_work_type
        double precision, allocatable :: roe_averages(:,:)
    end type rp_work_type
    
contains

    function new_rp_params(file_path) result(rp_params)

        implicit none
        character(len=*), intent(in) :: file_path
        type(rp_params_type), pointer :: rp_params

        integer :: ios

        allocate(rp_params)

        ! Read in data file
        open(unit=13, file=file_path, iostat=ios, status="old", action="read")
        if ( ios /= 0 ) stop "Error opening file rp.data"
        read(13,*) rp_params%gamma
        close(13)
        
    end function new_rp_params
    
    function new_rp_work(solution) result(rp_work)

        use state_1d_module
        use grid_module
        use solution_1d_module

        implicit none
        type(solution_type), intent(in) :: solution
        type(rp_work_type), pointer :: rp_work
        
        type(state_type), pointer :: state
        integer, parameter :: NUM_THREADS = 1
        integer :: ios
        
        allocate(rp_work)
        
        state => get_state(solution,1)
        allocate(rp_work%roe_averages(NUM_THREADS,state%grid%n(1)))

    end function new_rp_work
    
end module rp_euler_module

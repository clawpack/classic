module rp_swe_module

    implicit none

    type rp_params_type
        double precision :: g
        double precision, pointer :: bathymetry(:)
    end type rp_params_type

    type rp_work_type

    end type rp_work_type

contains

    function new_rp_params(file_path,grid) result(rp_params)
        
        use grid_module

        implicit none
        
        character(len=*), intent(in) :: file_path
        type(grid_type), pointer :: grid
        
        type(rp_params_type), pointer :: rp_params
        
        integer :: ios

        allocate(rp_params)

        ! Read in data file
        open(13,file=file_path,iostat=ios,status="old",action="read")
        if ( ios /= 0 ) then
            print *, "Error opening file ",file_path
            stop
        endif

        read(13,*) rp_params%g
        close(13)

        ! Set bathy
        allocate(rp_params%bathymetry(grid%n(1)))
        rp_params%bathymetry = -1.d0

    end function new_rp_params

    function new_rp_work() result(rp_work)
        implicit none

        type(rp_work_type), pointer :: rp_work

    end function new_rp_work
    
end module rp_swe_module

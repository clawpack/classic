module clawdata_module

    implicit none
    
    ! Representation of the clawpack intput data file
    type clawdata_type

        ! Dimensions of q and aux
        integer :: num_dim, num_cells, num_eqn, num_aux

        ! Domain specification
        real(kind=8) :: lower, upper, t0

        ! Output
        integer :: num_output, out_style, num_step_out, n_stop
        real(kind=8) :: t_final
        real(kind=8), allocatable :: t_out(:)

        ! Solver attributes
        integer :: max_num_steps, order, transverse_order, verbosity
        integer :: source_splitting, capa_index, num_waves
        logical :: dt_variable, use_fwave
        real(kind=8) :: cfl_max_allowed, cfl_desired, dt_max_allowed, dt_init
        real(kind=8), allocatable :: limiters(:)

        ! Boundary conditions
        integer :: num_ghost, bc_lower, bc_upper
    end type clawdata_type


contains

    type(clawdata_type) function read_data_file(path) result(clawdata)
        
        implicit none

        character(len=*), intent(in), optional :: path

        ! Locals
        integer, parameter :: iounit = 55
        integer :: stat, ios, mw, i

        ! Open data file
        if (present(path)) then
            open(unit=iounit, file=path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) then
                print *,"Error opening data file at ",path
                stop
            endif
        else
            open(unit=iounit, file='./claw.data', iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error opening file ./claw.data"
        end if

        ! Begin reading file
        read(iounit,'(i2)') clawdata%num_dim
        if (clawdata%num_dim /= 1) then
            stop "Dimension is not 1, wrong data source files used."
        end if

        ! Number of grid cells
        read(iounit,*) clawdata%num_cells
        if (clawdata%num_cells <= 0) stop "Need number of cells > 0."

        ! I/O Controls
        read(iounit,*) clawdata%num_output
        read(iounit,*) clawdata%out_style
        if (clawdata%out_style == 1) then
            read(iounit,*) clawdata%t_final
            clawdata%num_step_out = 1
        else if (clawdata%out_style == 2) then
            allocate(clawdata%t_out(clawdata%num_output),stat=stat)
            if (stat /= 0) stop "Allocation of t_out failed!"
            read(iounit,*) (clawdata%t_out(i), i=1,clawdata%num_output)
            clawdata%num_step_out = 1
        else if (clawdata%out_style == 3) then
            read(iounit,*) clawdata%num_step_out, clawdata%n_stop
            clawdata%num_output = clawdata%n_stop
        else
            stop "Output style > 3 unimplemented."
        end if

        ! Time Stepping controls
        read(iounit,*) clawdata%dt_init
        read(iounit,*) clawdata%dt_max_allowed
        read(iounit,*) clawdata%cfl_max_allowed
        read(iounit,*) clawdata%cfl_desired
        read(iounit,*) clawdata%max_num_steps

        ! Input parameters for clawpack algorithms
        read(iounit,*) clawdata%dt_variable
        read(iounit,*) clawdata%order
        read(iounit,*) clawdata%transverse_order
        read(iounit,*) clawdata%verbosity
        read(iounit,*) clawdata%source_splitting
        read(iounit,*) clawdata%capa_index
        read(iounit,*) clawdata%num_aux
!         read(iounit,*) clawdata%use_fwave
        clawdata%use_fwave = .false.

        ! Number of equations and waves
        read(iounit,*) clawdata%num_eqn
        read(iounit,*) clawdata%num_waves
        allocate(clawdata%limiters(clawdata%num_waves),stat=stat)
        if (stat /= 0) stop "Allocation of limiters array failed!"
        read(iounit,*) (clawdata%limiters(mw), mw=1,clawdata%num_waves)

        ! Domain
        read(iounit,*) clawdata%t0
        read(iounit,*) clawdata%lower
        read(iounit,*) clawdata%upper

        ! Boundary conditions
        read(iounit,*) clawdata%num_ghost
        read(iounit,*) clawdata%bc_lower
        read(iounit,*) clawdata%bc_upper
        if ((clawdata%bc_lower == 2 .and. clawdata%bc_upper /= 2) .or.  &
            (clawdata%bc_lower /= 2 .and. clawdata%bc_upper == 2)) then

            print *, "Periodic boundary conditions requested but bc_lower and"
            print *, "bc_upper were not BOTH set to 2."
            stop
        end if

        close(iounit)

    end function read_data_file

end module clawdata_module

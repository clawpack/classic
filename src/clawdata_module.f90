module clawdata_module

    implicit none
    
    ! Representation of the clawpack intput data file
    type clawdata_type

        ! Dimensions of q and aux
        integer :: num_dim, num_eqn, num_aux
        integer, allocatable :: num_cells(:)

        ! Domain specification
        real(kind=8) :: t0
        real(kind=8), allocatable :: lower(:), upper(:)

        ! Output
        integer :: output_style, num_output_times, output_step_interval
        integer :: total_steps, output_format
        integer, allocatable :: output_q_components(:), output_aux_components(:)
        logical :: output_t0, output_aux_onlyonce
        real(kind=8) :: t_final
        real(kind=8), allocatable :: t_out(:)

        ! Solver attributes
        integer :: steps_max, order, transverse_waves, dimensional_split
        integer :: source_splitting, capa_index, num_waves, verbosity
        logical :: dt_variable, use_fwaves
        real(kind=8) :: cfl_max_allowed, cfl_desired, dt_max_allowed, dt_init
        real(kind=8), allocatable :: limiters(:)

        ! Boundary conditions
        integer :: num_ghost
        integer, allocatable :: bc_lower(:), bc_upper(:)

        ! Restart settings
        logical :: restart
        character(len=128) :: restart_file

    end type clawdata_type


contains

    type(clawdata_type) function read_data_file(path) result(clawdata)
        
        implicit none

        character(len=*), intent(in), optional :: path

        ! Locals
        integer, parameter :: iounit = 55
        integer :: stat, ios, mw, i, err

        ! Open data file
        if (present(path)) then
            call open_data_file(iounit, path)
        else
            call open_data_file(iounit, 'claw.data')
        end if

        ! Begin reading file
        read(iounit,'(i2)') clawdata%num_dim
        if (clawdata%num_dim /= 1) then
            stop "Dimension is not 1, wrong data source files used."
        end if

        ! Allocate all variable dimension arrays
        allocate(clawdata%num_cells(clawdata%num_dim), stat=err)
        if (err /= 0) stop "*** ERROR *** - Allocation request denied."
        allocate(clawdata%lower(clawdata%num_dim), stat=err)
        if (err /= 0) stop "*** ERROR *** - Allocation request denied."
        allocate(clawdata%upper(clawdata%num_dim), stat=err)
        if (err /= 0) stop "*** ERROR *** - Allocation request denied."
        allocate(clawdata%bc_lower(clawdata%num_dim), stat=err)
        if (err /= 0) stop "*** ERROR *** - Allocation request denied."
        allocate(clawdata%bc_upper(clawdata%num_dim), stat=err)
        if (err /= 0) stop "*** ERROR *** - Allocation request denied."
        
        ! Domain
        read(iounit,*) clawdata%lower
        read(iounit,*) clawdata%upper
        read(iounit,*) clawdata%num_cells
        if (.not. all(clawdata%num_cells > 0, dim=1)) stop "Need number of cells > 0."
        read(iounit,*)

        ! Base array sizes
        read(iounit,*) clawdata%num_eqn
        read(iounit,*) clawdata%num_waves
        read(iounit,*) clawdata%num_aux
        read(iounit,*)

        read(iounit,*) clawdata%t0
        read(iounit,*)

        ! I/O Controls
        read(iounit,*) clawdata%output_style
        if (clawdata%output_style == 1) then
            read(iounit,*) clawdata%num_output_times
            read(iounit,*) clawdata%t_final
            read(iounit,*) clawdata%output_t0
        else if (clawdata%output_style == 2) then
            read(iounit,*) clawdata%num_output_times
            allocate(clawdata%t_out(clawdata%num_output_times),stat=stat)
            if (stat /= 0) stop "Allocation of t_out failed!"
            read(iounit,*) (clawdata%t_out(i), i=1,clawdata%num_output_times)
        else if (clawdata%output_style == 3) then
            read(iounit,*) clawdata%output_step_interval
            read(iounit,*) clawdata%total_steps
            read(iounit,*) clawdata%output_t0
        else
            stop "Output style > 3 unimplemented."
        end if
        read(iounit,*)

        ! Output details
        read(iounit,*) clawdata%output_format
        read(iounit,*) clawdata%output_q_components
        if (clawdata%num_aux > 0) then
            read(iounit,*) clawdata%output_aux_components
            read(iounit,*) clawdata%output_aux_onlyonce
        endif
        read(iounit,*)

        ! Time Stepping controls
        read(iounit,*) clawdata%dt_init
        read(iounit,*) clawdata%dt_max_allowed
        read(iounit,*) clawdata%cfl_max_allowed
        read(iounit,*) clawdata%cfl_desired
        read(iounit,*) clawdata%steps_max
        read(iounit,*)

        ! Input parameters for clawpack algorithms
        read(iounit,*) clawdata%dt_variable
        read(iounit,*) clawdata%order
        if (clawdata%num_dim > 1) then
            read(iounit,*) clawdata%transverse_waves
            read(iounit,*) clawdata%dimensional_split
        endif
        read(iounit,*) clawdata%verbosity
        read(iounit,*) clawdata%source_splitting
        read(iounit,*) clawdata%capa_index
        read(iounit,*) clawdata%use_fwaves
        read(iounit,*)

        ! Limiters to be used
        allocate(clawdata%limiters(clawdata%num_waves),stat=stat)
        if (stat /= 0) stop "Allocation of limiters array failed!"
        read(iounit,*) (clawdata%limiters(mw), mw=1,clawdata%num_waves)
        read(iounit,*)

        ! Boundary conditions
        read(iounit,*) clawdata%num_ghost
        read(iounit,*) clawdata%bc_lower
        read(iounit,*) clawdata%bc_upper
        if ((clawdata%bc_lower(1) == 2 .and. clawdata%bc_upper(1) /= 2) .or.  &
            (clawdata%bc_lower(1) /= 2 .and. clawdata%bc_upper(1) == 2)) then

            print *, "Periodic boundary conditions requested but bc_lower and"
            print *, "bc_upper were not BOTH set to 2."
            stop
        end if

        ! Restart settings
        read(iounit,*) clawdata%restart
        read(iounit,*) clawdata%restart_file

        close(iounit)

    end function read_data_file

end module clawdata_module

module clawdata_module

    use precision_module, only: dp
    use utils, only: stop_error, open_data_file
    implicit none
    
    ! Representation of the clawpack intput data file
    type clawdata_type

        ! Dimensions of q and aux
        integer :: num_dim, num_eqn, num_aux
        integer, allocatable :: num_cells(:)

        ! Domain specification
        integer :: capa_index
        real(dp) :: t0
        real(dp), allocatable :: lower(:), upper(:)

        ! Output
        integer :: output_style, num_output_times, output_step_interval
        integer :: total_steps, output_format
        integer, allocatable :: output_q_components(:), output_aux_components(:)
        logical :: output_t0, output_aux_onlyonce
        real(dp) :: t_final
        real(dp), allocatable :: t_out(:)

        ! Solver attributes
        integer :: steps_max, order, transverse_waves, dimensional_split
        integer :: source_splitting, num_waves, verbosity
        logical :: dt_variable, use_fwaves
        real(kind=8) :: cfl_max_allowed, cfl_desired, dt_max_allowed, dt_init
        integer, allocatable :: limiters(:)

        ! Boundary conditions
        integer :: num_ghost
        integer, allocatable :: bc_lower(:), bc_upper(:)

        ! Restart settings
        logical :: restart
        character(len=128) :: restart_file

    end type clawdata_type

    interface new
        module procedure new_clawdata
    end interface new

contains

    subroutine new_clawdata(self,path)
        
        implicit none

        ! Input
        character(len=*), intent(in), optional :: path

        ! Output
        type(clawdata_type), intent(out) :: self

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
        read(iounit,'(i2)') self%num_dim
        if (self%num_dim /= 1) then
            call stop_error("Dimension is not 1, wrong data source files used.")
        end if
        
        ! Domain
        allocate(self%lower(self%num_dim), stat=err)
        if (err /= 0) call stop_error("*** ERROR *** - Allocation request denied.")
        read(iounit,*) self%lower
        allocate(self%upper(self%num_dim), stat=err)
        if (err /= 0) call stop_error("*** ERROR *** - Allocation request denied.")
        read(iounit,*) self%upper
        allocate(self%num_cells(self%num_dim), stat=err)
        if (err /= 0) call stop_error("*** ERROR *** - Allocation request denied.")
        read(iounit,*) self%num_cells
        if (.not. all(self%num_cells > 0, dim=1)) call stop_error("Need number of cells > 0.")
        read(iounit,*)

        ! Base array sizes
        read(iounit,*) self%num_eqn
        read(iounit,*) self%num_waves
        read(iounit,*) self%num_aux
        read(iounit,*)

        read(iounit,*) self%t0
        read(iounit,*)

        ! I/O Controls
        read(iounit,*) self%output_style
        self%output_step_interval = 1
        if (self%output_style == 1) then
            read(iounit,*) self%num_output_times
            read(iounit,*) self%t_final
            read(iounit,*) self%output_t0

        else if (self%output_style == 2) then
            read(iounit,*) self%num_output_times
            allocate(self%t_out(self%num_output_times),stat=stat)
            if (stat /= 0) call stop_error("Allocation of t_out failed!")
            read(iounit,*) (self%t_out(i), i=1,self%num_output_times)

        else if (self%output_style == 3) then
            read(iounit,*) self%output_step_interval
            read(iounit,*) self%total_steps
            read(iounit,*) self%output_t0
            self%num_output_times = self%total_steps

        else
            call stop_error("Output style > 3 unimplemented.")
        end if
        read(iounit,*)

        ! Output details
        read(iounit,*) self%output_format
        allocate(self%output_q_components(self%num_eqn), stat=stat)
        if (stat /= 0) call stop_error("Allocation of output_q_components failed!" )
        allocate(self%output_aux_components(self%num_aux), stat=stat)
        if (stat /= 0) call stop_error("Allocation of output_aux_components failed!" )
        read(iounit,*) self%output_q_components
        if (self%num_aux > 0) then
            read(iounit,*) self%output_aux_components
            read(iounit,*) self%output_aux_onlyonce
        endif
        read(iounit,*)

        ! Time Stepping controls
        read(iounit,*) self%dt_init
        read(iounit,*) self%dt_max_allowed
        read(iounit,*) self%cfl_max_allowed
        read(iounit,*) self%cfl_desired
        read(iounit,*) self%steps_max
        read(iounit,*)

        ! Input parameters for clawpack algorithms
        read(iounit,*) self%dt_variable
        read(iounit,*) self%order
        if (self%num_dim > 1) then
            read(iounit,*) self%transverse_waves
            read(iounit,*) self%dimensional_split
        endif
        read(iounit,*) self%verbosity
        read(iounit,*) self%source_splitting
        read(iounit,*) self%capa_index
        read(iounit,*) self%use_fwaves
        read(iounit,*)

        ! Limiters to be used
        allocate(self%limiters(self%num_waves),stat=stat)
        if (stat /= 0) call stop_error("Allocation of limiters array failed!")
        read(iounit,*) (self%limiters(mw), mw=1,self%num_waves)
        read(iounit,*)

        ! Boundary conditions
        read(iounit,*) self%num_ghost
        allocate(self%bc_lower(self%num_dim), stat=stat)
        if (stat /= 0) call stop_error("Allocation of bc_lower failed!")
        read(iounit,*) self%bc_lower
        allocate(self%bc_upper(self%num_dim), stat=stat)
        if (stat /= 0) call stop_error("Allocation of bc_upper failed!")
        read(iounit,*) self%bc_upper
        if ((self%bc_lower(1) == 2 .and. self%bc_upper(1) /= 2) .or.  &
            (self%bc_lower(1) /= 2 .and. self%bc_upper(1) == 2)) then

            print *, "Periodic boundary conditions requested but bc_lower and"
            print *, "bc_upper were not BOTH set to 2."
            stop
        end if

        ! Restart settings
        read(iounit,*) self%restart
        read(iounit,*) self%restart_file

        close(iounit)

    end subroutine new_clawdata

end module clawdata_module

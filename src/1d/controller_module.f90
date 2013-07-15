module controller_module

    implicit none

contains

    logical function run_simulation(solution, solver, clawdata) result(success)
        
        use solution_module, only: solution_type, output
        use solver_module, only: solver_type, evolve_to_time
        use clawdata_module, only: clawdata_type

        implicit none

        ! Arguments
        type(solution_type), intent(in out) :: solution
        type(solver_type), intent(in out) :: solver
        type(clawdata_type), intent(in out) :: clawdata

        ! Locals
        integer :: frame, stat
        integer, allocatable :: no_aux_output(:)
        real(kind=8) :: dt_out, t_start

        ! Log file settings
        integer, parameter :: LOG_UNIT = 213
        character(len=*), parameter :: LOG_FILE_NAME = "claw.log"

        ! Format parameters
        character(len=*), parameter :: OUTPUT_FORMAT = &
            "('CLAW: Frame ',i4,' output files done at time t =',d12.4,/)"
        character(len=*), parameter :: LOG_FORMAT = "('  ',i4)"

        ! Initialize return value
        success = .false.

        ! No auxillary output array
        allocate(no_aux_output(solution%num_aux), stat=stat)
        if (stat /= 0) stop "Not able to allocate work array!"
        no_aux_output = 0

        ! Output initial time if requested
        if (clawdata%output_t0) then
            if (clawdata%output_aux_onlyonce) then
                call output(solution, clawdata%t0, 0,                         &
                                                clawdata%output_q_components, &
                                                clawdata%output_aux_components)
            else
                call output(solution, clawdata%t0, 0,                         &
                                                clawdata%output_q_components, &
                                                no_aux_output)
            endif
        endif

        ! Setup time output steps
        if (clawdata%output_style == 1) then
                dt_out = (clawdata%t_final - clawdata%t0) &
                                    / real(clawdata%num_output_times,kind=8)
        end if

        ! Open log file
        open(unit=LOG_UNIT, file=LOG_FILE_NAME, iostat=stat, status="unknown", action="write")
        if ( stat /= 0 ) then
            print *, "Error opening log file with error ",stat
            success = .false.
            stop
        endif

        ! Primary output loop
        do frame=1,clawdata%num_output_times
            ! Calculate next time to end at
            select case(clawdata%output_style)
                case(1)
                    t_start = solution%t 
                    success = evolve_to_time(t_start + frame * dt_out,solution,solver)
                case(2)
                    success = evolve_to_time(clawdata%t_out(frame),solution,solver)
                case(3)
                    ! Single step
                    success = evolve_to_time(huge(1.d0), solution, solver, .true.)
            end select

            ! Break out of stepping loop
            if (.not.success) then
                print *,"ERROR - An error occured!"
                exit
            endif

            if (frame / clawdata%output_step_interval == frame) then
                if (clawdata%output_aux_onlyonce) then
                    call output(solution, solution%t, frame,                  &
                                            clawdata%output_q_components,     &
                                            clawdata%output_aux_components)
                else
                    call output(solution, solution%t, frame,                  &
                                            clawdata%output_q_components,     &
                                            no_aux_output)
                endif

                print OUTPUT_FORMAT, frame, solution%t
                write(LOG_UNIT, OUTPUT_FORMAT) frame, solution%t
                write(LOG_UNIT, LOG_FORMAT) frame
            endif
        end do

        ! Close log file
        close(LOG_UNIT)

    end function run_simulation

end module controller_module
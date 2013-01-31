program advection1d

    use solver_module
    use solution_module
    use rp1_advection, only: u, rp1, rp_ptwise

    use clawdata_module

    implicit none
    
    type(solution_type) :: solution
    type(solver_type) :: solver
    type(clawdata_type) :: clawdata

    real(kind=8) :: beta

    ! Read in advective speed and initial condition data
    call open_data_file(13,'./setprob.data')
    read(13,*) u
    read(13,*) beta
    close(13)

    ! Read input data
    call new(clawdata,'./claw.data')

    ! Initialize solution
    call new(solution,clawdata)
    solution%q(1,:) = exp(-beta * (solution%centers - 0.3d0)**2)

    ! Initialize solver
    call new(solver,clawdata)
    solver%rp1 => rp_ptwise

    ! Start output loop
    call run_simulation(solution,solver,clawdata)

contains

    subroutine run_simulation(solution, solver, clawdata)
        
        use solution_module, only: solution_type
        use solver_module, only: solver_type
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

        ! Format parameters
        character(len=*), parameter :: OUTPUT_FORMAT = &
            "('CLAW: Frame ',i4,' output files done at time t =',d12.4,/)"

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

        ! Primary output loop
        do frame=1,clawdata%num_output_times
            ! Calculate next time to end at
            select case(clawdata%output_style)
                case(1)
                    t_start = solution%t 
                    call evolve_to_time(t_start + frame * dt_out,solution,solver)
                case(2)
                    call evolve_to_time(clawdata%t_out(frame),solution,solver)
                case(3)
                    ! Single step
                    call evolve_to_time(huge(1.d0), solution, solver, .true.)
            end select

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
                ! TODO: Log this output step
            endif
        end do

    end subroutine run_simulation

end program advection1d

logical function evolve_to_time(t_end, solution, solver, single_step) result success
    
    use solution_module, only: solution_type
    use solver_module, only: solver_type

    implicit none

    ! Arguments
    real(kind=8), intent(in) :: t_end
    type(solution_type), intent(in out) :: solution
    type(solver_type), intent(in out) :: solver
    logical, optional, intent(in) :: single_step

    ! Locals
    integer :: num_steps
    real(kind=8) :: t_old, t_new, t_start

    character(len=*), parameter :: stat_msg = "('CLAW1... Step',i4,"  // &
                         "'   Courant number =',f6.3,'  dt =',d12.4," // &
                         "'  t =',d12.4)" 

    ! Initialize return status
    success = .false.

    ! Initialize status storage
    num_steps = 0
    t_start = solution%t

    ! Support single stepping mode
    if (present(single_step)) then
        if (single_step) solver%steps_max = 1
    endif

    ! Primary Loop
    primary_loop: do num_steps=1,solver%steps_max

        t_old = solution%t
        t_new = solution%t + solver%dt

        ! Adjust dt to hit t_end exactly if we are near end of computation
        if (t_old + solver%dt > t_end .and. t_start < t_end) then
            solver%dt = t_end - t_old
        endif

        ! Save q in case we need to retake this step
        if (solver%dt_variable) then
            solver%q_old => solution%q
        endif

        ! Main steps in algorithm
        
        ! Extend data from actual domain to bordering boundary cells
        call set_boundary_conditions(solution,solver)

        ! Call user-supplied routine
        call before_step(solution,solver)

        ! If Strang splitting is used take a step on the source term
        if (solver%source_splitting == 2) then
            call source_term(t_old, solver%dt / 2.d0, solution, solver)
        endif

        ! Take a single-step on the homogeneous conservation law
        call hyperbolic_step(solution,solver)

        ! Take a step on the source term, if Strang splitting used take another
        ! half time step, if Godunov splitting take a full time step
        if (solver%source_splitting == 1) then
            call source_term(solution%t + solver%dt,solver%dt,solution,solver)
        else if (solver%source_splitting == 2) then
            call source_term(solution%t + solver%dt * 0.5d0,solver%dt * 0.5d0,solution,solver)
        end if

        ! Write out status message if verbosity /= 0
        if (solver%verbosity > 0) then
            print stat_msg, num_steps, solver%cfl, solver%dt, t_old + solver%dt
        end if

        ! Choose new time step if allowed
        call choose_new_dt(solver)

        ! Check CFL condition
        if (.not.check_CFL(solver)) then
            ! Reject this step and reset q
            print *, 'CLAW - Rejecting step... Courant number too large'
            solution%q = solver%q_old

            if (.not.solver%dt_variable) then
                ! Not taking variable time steps so exit due to CFL violation
                solver%cfl_max = max(solver%cfl,solver%cfl_max)
                success = .false.
                return
            endif

            ! Go back to beginning of loop, already should have choosen the 
            ! correct time step from above
            cycle primary_loop
        else
            ! Accept this step
            solver%cfl_max = max(solver%cfl,solver%cfl_max)
        end if

        ! Update solution time
        solution%t = t_new

        ! Check to see if we are done
        if (solution%t >= t_end) then
            exit primary_loop
        endif

    end do primary_loop

    success = .true.
    
end function evolve_to_time
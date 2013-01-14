subroutine evolve_to_time(t_end,solution,solver)
    
    use solution_module, only: solution_type
    use solver_module, only: solver_type, choose_new_dt, check_CFL

    implicit none

    real(kind=8), intent(in) :: t_end

    type(solution_type), intent(in out) :: solution
    type(solver_type), intent(in out) :: solver

    ! Locals
    integer :: num_steps
    real(kind=8) :: t_old, t_new, t_start

    character(len=*), parameter :: stat_msg = "('CLAW1... Step',i4,"  // &
                         "'   Courant number =',f6.3,'  dt =',d12.4," // &
                         "'  t =',d12.4)" 

    ! Initialize status storage
    num_steps = 0
    t_start = solution%t

    ! Primary Loop
    primary_loop: do num_steps=1,solver%max_num_steps

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
        if (solver%source_split == 2) then
            call source_term(t_old, solver%dt / 2.d0, solution, solver)
        endif

        ! Take a single-step on the homogeneous conservation law
        call step1()

        ! Take a step on the source term, if Strang splitting used take another
        ! half time step, if Godunov splitting take a full time step
        if (solver%source_split == 1) then
            call source_term(solution%t + solver%dt,solver%dt,solution,solver)
        else if (solver%source_split == 2) then
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
            ! Reject this step
            stop "Rejected CFL step not supported."
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
    
end subroutine evolve_to_time
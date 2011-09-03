module solver_module

    ! Status of solver
    integer :: num_steps
    double precision :: t,cfl,cfl_max,dt_min,dt_max
    
    ! Parameters for the solver, note that these are not actually solvers as
    ! they are input parameters
    integer :: MAX_STEPS
    logical :: dt_variable,fwave
    double precision :: cfl_desired
    integer :: limiters,order,src_split,verbosity
    
    external start_step
    
contains

    ! This needs to be filled in
    subroutine status_message(out)

        implicit none
        integer, intent(in) :: out
        
        print *,"here"
        
    end subroutine status_message

    subroutine evolve_to_time(solution,t)

        use solution_module

        implicit none
    
        ! Input
        double precision, optional, intent(in) :: t
        type(solution_type), intent(inout) :: solution
        
        ! Locals for time stepping
        logical :: retake_step, take_one_step
        double precision :: t_start,t_old,t_half,dt_2
        
        ! Parse input arguments
        if (present(t)) then
            take_one_step = .false.
        else
            take_one_step = .true.
        endif
        
        ! Reset status of solver
        t_start = t
        cfl_max = cfl_desired
        dt_min = dt
        dt_max = dt
        num_steps = 0
        
        ! Setup time step for evolution
        if (.not.dt_variable) then
            if (take_one_step) then
                max_steps = 1
            else
                max_steps = int((t_end - t_start + 1d-10) / dt)
                if (abs(max_steps*dt - (t_end - t_start)) > 1e-5 * t_end-t_start)
                    stop "dt does not divide (t_end-t_start) and dt is fixed!"
                endif
            endif
        endif
        if (dt_variable.and.(cfl_desired > cfl_max)) then
            stop "Variable time stepping and desired CFL > maximum CFL."
        endif
        if (t_end <= t_start) then
            print *, "Already at or beyond end tiem: no evolution required."
            max_steps = 0
        endif
    
        ! ========================================================================
        !  Main loop
        ! ========================================================================
        do num_steps=1,MAX_STEPS
        
            ! Strang splitting time step needs half time steps
            dt_2 = 0.5d0 * dt 
            t_half = t + dt_2
            t = t_old + dt  ! Time at end of step
        
            ! Extend state q to boundary ghost cells
            q = bc() ! *** 
            
            ! Adjust dt so that we hit t_end exactly if we are near t_end
            if ((t + dt > t_end).and.(t_start < t_end)    &
                                .and.(.not.take_one_step)) then 
                dt = t_end - t
            endif
            
            ! Keep a backup in case we need to retake a time step
            if (dt_variable) then ! ***
                q_old = q.copy()
                t_old = t
            endif
            retake_step = .false.
            
            ! Take a full time step dt on the equations
            call step(solution)
            
            ! Check to make sure that the Courant number was not too large
            if (cfl <= cfl_max) then
                ! Accept this step
                cfl_max = max(cfl,cfl_max)
                if (dt_variable) then
                    t = t + dt
                else
                    ! Avoid roundoff error if dt_variable == false
                    t = t_start + (n+1) * dt
                endif
                
                ! Verbose messaging
                if (verbose > 0) then
                    call status_message() !***
                endif
                
                ! See if we are finished yet
                if ((t >= t_end).or.take_one_step) then
                    exit
                endif
            else
                ! Reject this step
                print *,"Rejecting time step, CFL number too large..."
                if (dt_variable) then
                   ! Recopy back backup
                   t = t_old
                   retake_step = .true.
                else
                    ! Giveup, we cannot adapt
                    cfl_max = max(cfl,cfl_max)
                    print *,"CFL",cfl_max," too large, giving up!"
                    stop
                endif
            endif
            
            ! Choose a new time step
            if (dt_variable) then
                if (cfl > 0.d0) then
                    dt = min(dt_max,dt * cfl_desired / cfl)
                    dt_min = min(dt,dt_min)
                    dt_max = max(dt,dt_max)
                else
                    dt = dt_max
                endif
            endif
        enddo
        ! ====== End of Main Time Stepping Loop ==============================
        
        if (dt_variable.and.(t < t_end).and.(num_steps == max_steps)) then
            stop "Maximum number of time step have been taken."
        endif
        
    end subroutine evolve_to_time
    
    subroutine step(solution)

        implicit none
        
        ! Call user defined step setup routine
        call start_step()
        
        ! Call source term if Strang splitting is being used
        if (src_split == 2) then
            call src(solution,t,dt/2.d0)
        endif
        
        ! Take a time step on the homogeneous hyperbolic problem
        call homogeneous_step(solution)
        
        ! Check if we have violated the CFL condition, if we did, return
        ! immedieatly to evolve_to_time and let it deal with picking a new dt
        if (cfl >= cfl_max) then
            return
        endif
        
        ! Source term splitting steps
        if (src_split > 0) then
            ! Godunov splitting
            if (src_split == 1) then
                call src(solution,t,dt)
            ! Strang splitting
            else if (src_split == 2) then
               call src(solution,t + dt / 2.d0,dt/2.d0)
            endif
        endif
    end subroutine step

    ! This routine assumes there is one grid and one solution, be careful with
    ! it's application
    subroutine homogenous_step(solution)

        implicit none
        type(solution_type), intent(intout) :: solution
        
        ! Local
        type(grid_type), pointer :: grid
        type(state_type), pointer :: state
        integer :: mx
        double precision :: dtdx(-mbc:mx+mbc)
        
        ! Extract relevant quantities
        grid = solution%grids(1)%grid
        state = solution%states(1)%state
        mx = grid%n(1)
        dx = grid%d(1)

        ! Take into account grid mapping
        if (mcapa >= 0) then
            dtdx = dt / (dx * aux(mcapa,:))
        else
            dtdx = dt/dx
        endif
        
        ! Solve Riemann problem
        q_l = q_bc(:,:-1)
        q_r = q_bc(:,-1,:)
        aux_l = aux_bc(:,:-1)
        aux_r = aux_bc(:,0:)
        call rp()
        
        ! Update loops from fluxes
        do m=1,meqn
            do i=1,mx+1
                q(m,i) = q(m,i) - dtdx(i) * apdq(m,i)
            enddo
        enddo
        
        ! Compute maximum wave speed
        cfl = 0.d0
        do mw=1,mwaves
            do i=1,mx+1
                cfl = max(cfl,dtdx(i)*s(i,mw), -dtdx(i-1)*s(i,mw))
            enddo
        enddo
        
        ! Perform higher order limiting
        if (order == 2) then
            f = 0.d0
            
            ! Apply limiters
            if (any(limiter > 0)) then
                wave = limiters() !***
            endif
            
            ! Compute correction fluxes for second order q_{xx} terms
            dtdxave = 0.5d0 * (dtdx()) !***
            if (fwave) then
                
            else
                
            endif
            
            ! Update q by differencing correction fluxes
            do m=1,meqn
                q = 
            enddo
        endif
        
        ! Update q in state
        state%q = q_bc
        cfl = update_global_mx
        
    end subroutine homogenous_step
    
end module solver_module


! ==============================================================================
!  Take one time step, updating q.
! ------------------------------------------------------------------------------
subroutine hyperbolic_step(solution,solver)

    use solution_module, only: solution_type
    use precision_module, only: dp
    use geometry_module, only: geometry_type
    use solver_module, only: solver_type, limiter
    use utils, only: stop_error

    implicit none

    ! Function Arguments
    type(solution_type), intent(in out) :: solution
    type(solver_type), intent(in out) :: solver

    ! Local storage
    integer :: i, m, mw
    type(geometry_type) :: geometry
    real(dp) :: cfl

    ! Take apart solution and solver data types for ease of reference below
    associate(aux => solution%aux, f => solver%f,           &
              apdq => solver%apdq, amdq => solver%amdq, wave => solver%wave, &
              s => solver%s, dtdx => solver%dtdx)

        ! Calculate dt/dx ratio, handles non-uniform grids
        if (solution%capa_index == 0) then
            ! Uniform grid being used, set dtdx(:) = dt/dx
            dtdx = solver%dt / solution%dx(1)
        else
            ! Use capacity function to compute spatially variable dx
            dtdx = solver%dt / (solution%dx(1)*aux(:,solution%capa_index))
        endif

        ! Solve Riemann problem at each interface - Note that if both the
        ! vectorized solver function and the point-wise function pointers are
        ! set the vectorized form is used
        if (associated(solver%rp1_vectorized)) then
            ! Not sure if x should be included due to vector construction overhead
            geometry%x = 0.0_dp  !solution%lower(1) + (i-0.5_dp) * solution%dx(1)
            geometry%dx = solution%dx(1)
            geometry%t = solution%t
            geometry%dt = solver%dt
            call solver%rp1_vectorized(solution%num_eqn,             &
                                      solution%num_aux,              &
                                      solver%num_ghost,              &
                                      solution%num_cells(1),         &
                                      solver%num_waves,              &
                                      solver%rp_data,                &
                                      geometry,                      &
                                      solution%q, solution%q,                          &
                                      aux, aux,                      &
                                      solver%wave,                   &
                                      solver%s,                      &
                                      solver%amdq,                   &
                                      solver%apdq)
        else if (associated(solver%rp1_ptwise)) then
            do i=2-solver%num_ghost,solution%num_cells(1)+solver%num_ghost
                geometry%x = solution%lower(1) + (i-0.5_dp) * solution%dx(1)
                geometry%dx = solution%dx(1)
                geometry%t = solution%t
                geometry%dt = solver%dt
                call solver%rp1_ptwise(solution%num_eqn,              &
                                      solution%num_aux,               &
                                      solver%num_waves,               &
                                      solver%rp_data,                 &
                                      geometry,                       &
                                      solution%q(:,i-1), solution%q(:,i),               &
                                      aux(:,i-1), aux(:,i),           &
                                      wave(:,:,i),                    &
                                      s(:,i),                         &
                                      amdq(:,i),                      &
                                      apdq(:,i))
            end do
        else
            call stop_error("No Riemann solver function specified.")
        endif

        ! Modify q for Godunov update
        !  Note this may not correspond to a conservative flux-differencing for
        !  equations not in conservative form.  It is conservative if 
        !  amdq + apdq = f(q(i)) - f(q(i-1)).
        forall (i=1:solution%num_cells(1)+1, m=1:solution%num_eqn)
            solution%q(m,i) = solution%q(m,i) - dtdx(i) * apdq(m,i)
            solution%q(m,i-1) = solution%q(m,i-1) - dtdx(i-1) * amdq(m,i)
        end forall

        ! Compute maximum wave speed for CFL condition
        do mw=1,solver%num_waves
            cfl = max(maxval(dtdx * s(mw,:)),maxval(-dtdx * s(mw,:)))
        end do

        ! Set new CFL from this time step
        solver%cfl = cfl

        if (solver%order > 1) then
            ! Compute correction fluxes for second order q_{xx} terms
            solver%f = 0.d0

            ! Apply limiters
            if (any(solver%limiters /= 0)) then
                call limiter(solution%num_cells(1), solver%num_ghost,        &
                             solution%num_eqn, solver%num_waves, wave, s,    &
                             solver%limiters)
            end if

            ! Construct flux corrections
            do i = 1, solution%num_cells(1)+1
            do m = 1, solution%num_eqn
            do mw = 1, solver%num_waves
                f(m,i) = f(m,i) + 0.5d0 * abs(s(mw,i))      &
                    * (1.d0 - abs(s(mw,i)) * (0.5d0 * (dtdx(i-1) + dtdx(i)))) &
                    * wave(m,mw,i)
            end do
            end do
            end do

            ! Update q by differencing correction fluxes
            do m=1,solution%num_eqn
                solution%q(m,1:solution%num_cells(1)) = solution%q(m,1:solution%num_cells(1))  &
                                             - dtdx(1:solution%num_cells(1)) &
            * (f(m,2:solution%num_cells(1)+1) - f(m,1:solution%num_cells(1)))
            enddo
        endif

    end associate

end subroutine hyperbolic_step

! ==============================================================================
!  Take one time step, updating q.
! ------------------------------------------------------------------------------
subroutine hyperbolic_step(solution,solver)

    use solution_module, only: solution_type
    use solver_module, only: solver_type

    implicit none

    ! Function Arguments
    type(solution_type), intent(in out) :: solution
    type(solver_type), intent(in out) :: solver

    ! Local storage
    integer :: i, m, mw
    real(kind=8) :: cfl

    ! Take apart solution and solver data types for ease of reference below
    associate(q => solution%q, aux => solution%aux)

        ! Solve Riemann problem at each interface
        call solver%rp1(solution%num_eqn,              &
                        solution%num_aux,              &
                        solver%num_ghost,              &
                        solution%num_cells(1),         &
                        solver%num_waves,              &
                        q, q,                          &
                        aux, aux,                      &
                        solver%wave,                   &
                        solver%s,                      &
                        solver%amdq,                   &
                        solver%apdq)

        ! Calculate dt/dx ratio, handles non-uniform grids
        if (solution%capa_index == 0) then
            ! Uniform grid being used, set dtdx(:) = dt/dx
            solver%dtdx = solver%dt / solution%dx(1)
        else
            ! Use capacity function to compute spatially variable dx
            solver%dtdx = solver%dt / (solution%dx(1)*aux(:,solution%capa_index))
        endif

        ! Modify q for Godunov update
        !  Note this may not correspond to a conservative flux-differencing for
        !  equations not in conservative form.  It is conservative if 
        !  amdq + apdq = f(q(i)) - f(q(i-1)).
        forall (i=1:solution%num_cells(1)+1, m=1:solution%num_eqn)
            q(m,i) = q(m,i) - solver%dtdx(i) * solver%apdq(m,i)
            q(m,i-1) = q(m,i-1) - solver%dtdx(i-1) * solver%apdq(m,i)
        end forall

        ! Compute maximum wave speed for CFL condition
        do mw=1,solver%num_waves
            cfl = maxval(solver%dtdx * solver%s(mw,:))
            cfl = max(cfl,maxval(-solver%dtdx * solver%s(mw,:)))
        end do

        ! Set new CFL from this time step
        solver%cfl = cfl

        if (solver%order > 1) then
            ! Compute correction fluxes for second order q_{xx} terms
            stop "Second order not implemented yet"
        endif

    end associate

end subroutine hyperbolic_step
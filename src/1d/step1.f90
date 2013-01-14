! ==============================================================================
!
!     # Take one time step, updating q.
!
!     method(1) = 1   ==>  Godunov method
!     method(1) = 2   ==>  Slope limiter method
!     mthlim(p)  controls what limiter is used in the pth family
!
!
!     amdq, apdq, wave, s, and f are used locally:
!
!     amdq(1-mbc:maxmx+mbc, meqn) = left-going flux-differences
!     apdq(1-mbc:maxmx+mbc, meqn) = right-going flux-differences
!        e.g. amdq(i,m) = m'th component of A^- \Delta q from i'th Riemann
!                         problem (between cells i-1 and i).
!
!     wave(1-mbc:maxmx+mbc, meqn, mwaves) = waves from solution of
!                                           Riemann problems,
!            wave(i,m,mw) = mth component of jump in q across
!                           wave in family mw in Riemann problem between
!                           states i-1 and i.
!
!     s(1-mbc:maxmx+mbc, mwaves) = wave speeds,
!            s(i,mw) = speed of wave in family mw in Riemann problem between
!                      states i-1 and i.
!
!     f(1-mbc:maxmx+mbc, meqn) = correction fluxes for second order method
!            f(i,m) = mth component of flux at left edge of ith cell 
! ------------------------------------------------------------------------------
subroutine step1(solution,solver)

    use solution_module
    use solver_module

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
        call rp1(solution%num_eqn,              &
                 solution%num_aux,              &
                 solution%num_ghost,            &
                 solution%num_cells,            &
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
            solver%dtdx = solver%dt / solution%dx
        else
            ! Use capacity function to compute spatially variable dx
            solver%dtdx = solver%dt / (solution%dx * aux(:,solution%capa_index))
        endif

        ! Modify q for Godunov update
        !  Note this may not correspond to a conservative flux-differencing for
        !  equations not in conservative form.  It is conservative if 
        !  amdq + apdq = f(q(i)) - f(q(i-1)).
        forall (i=1:solution%num_cells+1, m=1:solution%num_eqn)
            q(m,i) = q(m,i) - solver%dtdx(i) * solver%apdq(m,i)
            q(m,i-1) = q(m,i-1) - solver%dtdx(i-1) * solver%apdq(m,i)
        end forall

        ! Compute maximum wave speed for CFL condition
        do mw=1,solver%num_waves
            cfl = maxval(solver%dtdx * solver%s(:,mw))
            cfl = max(cfl,maxval(-solver%dtdx * solver%s(:,mw)))
        end do
!         forall(i=1:solution%num_cells+1,mw=1:solver%num_waves)
!             cfl = max(cfl, solver%dtdx(i) * s(i,mw), -solver%dtdx(i-1) * s(i,mw))
!         end forall

        ! Set new CFL from this time step
        solver%cfl = cfl

        if (solver%order > 1) then
            ! Compute correction fluxes for second order q_{xx} terms
            stop "Second order not implemented yet"
        endif

    end associate

end subroutine step1
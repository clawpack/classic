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

end subroutine step1
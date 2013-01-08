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

    with(solution)
    

end subroutine step1

!       dimension    q(1-mbc:maxmx+mbc, meqn)
!       dimension  aux(1-mbc:maxmx+mbc, *)
!       dimension    f(1-mbc:maxmx+mbc, meqn)
!       dimension    s(1-mbc:maxmx+mbc, mwaves)
!       dimension wave(1-mbc:maxmx+mbc, meqn, mwaves)
!       dimension amdq(1-mbc:maxmx+mbc, meqn)
!       dimension apdq(1-mbc:maxmx+mbc, meqn)
!       dimension dtdx(1-mbc:maxmx+mbc)
!       dimension method(7),mthlim(mwaves)
!       logical limit
! c
! c     # check if any limiters are used:
!       limit = .false.
!       do 5 mw=1,mwaves
!          if (mthlim(mw) .gt. 0) limit = .true.
!    5     continue
! c
!       mcapa = method(6)
!       do 10 i=1-mbc,mx+mbc
!          if (mcapa.gt.0) then
!              if (aux(i,mcapa) .le. 0.d0) then
!                 write(6,*) 'Error -- capa must be positive'
!                 stop
!                 endif
!              dtdx(i) = dt / (dx*aux(i,mcapa))
!             else
!              dtdx(i) = dt/dx
!             endif
!    10    continue
! c
! c
! c
! c     # solve Riemann problem at each interface 
! c     -----------------------------------------
! c
!       call rp1(maxmx,meqn,mwaves,mbc,mx,q,q,aux,aux,wave,s,amdq,apdq)
! c
! c     # Modify q for Godunov update:
! c     # Note this may not correspond to a conservative flux-differencing
! c     # for equations not in conservation form.  It is conservative if
! c     # amdq + apdq = f(q(i)) - f(q(i-1)).
! c
!       do 40 i=1,mx+1
!          do 40 m=1,meqn
!             q(i,m) = q(i,m) - dtdx(i)*apdq(i,m)
!             q(i-1,m) = q(i-1,m) - dtdx(i-1)*amdq(i,m)
!    40       continue

! c
! c     # compute maximum wave speed:
!       cfl = 0.d0
!       do 50 mw=1,mwaves
!          do 45 i=1,mx+1
! c          # if s>0 use dtdx(i) to compute CFL,
! c          # if s<0 use dtdx(i-1) to compute CFL:
!            cfl = dmax1(cfl, dtdx(i)*s(i,mw), -dtdx(i-1)*s(i,mw))
!    45      continue
!    50    continue
! c
!       if (method(2) .eq. 1) go to 900
! c
! c     # compute correction fluxes for second order q_{xx} terms:
! c     ----------------------------------------------------------
! c
!       do 100 m = 1, meqn
!             do 100 i = 1-mbc, mx+mbc
!                f(i,m) = 0.d0
!   100          continue
! c
! c      # apply limiter to waves:
!       if (limit) call limiter(maxmx,meqn,mwaves,mbc,mx,wave,s,mthlim)
! c
!       do 120 i=1,mx+1
!          do 120 m=1,meqn
!             do 110 mw=1,mwaves
!                dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
!                f(i,m) = f(i,m) + 0.5d0 * dabs(s(i,mw))
!      &             * (1.d0 - dabs(s(i,mw))*dtdxave) * wave(i,m,mw)
!   110          continue
!   120       continue
! c
! c
!   140 continue
! c
! c     # update q by differencing correction fluxes 
! c     ============================================
! c
! c     # (Note:  Godunov update has already been performed above)
! c
!       do 150 m=1,meqn
!          do 150 i=1,mx
!             q(i,m) = q(i,m) - dtdx(i) * (f(i+1,m) - f(i,m))
!   150       continue
! c
!   900 continue
!       return
!       end

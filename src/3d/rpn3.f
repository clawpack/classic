c
c
c
c     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,
     &			auxl,auxr,wave,s,amdq,apdq)
c     ==================================================================
c
c     # Solve Riemann problems for the 3D hyperbolic problem.
c     # This is a dummy routine and is only intended
c     # to illustrate the format of this routine.  
c     # See various example directories for better examples.
c
c       -----------------------------------------------------------
c
c     # solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.  
c
c     # Note that the i'th Riemann problem has left state qr(:,i-1)
c     #                                    and right state ql(:,i)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
      implicit real*8(a-h,o-z)
c
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension amdq(meqn, 1-mbc:maxm+mbc)
      dimension apdq(meqn, 1-mbc:maxm+mbc)
      dimension auxl(maux, 1-mbc:maxm+mbc)
      dimension auxr(maux, 1-mbc:maxm+mbc)

      write(6,*) 'You must provide a valid Riemann solver'
      stop
c
      return
      end

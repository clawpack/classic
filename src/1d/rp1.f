c
c
c =========================================================
      subroutine rp1(meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &           wave,s,amdq,apdq,maux)
c =========================================================
c
c     # Solve Riemann problems for the 1D hyperbolic problem.
c     # This is a dummy routine and is only intended
c     # to illustrate the format of this routine.  
c     # See various example directories for better examples.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c
      implicit double precision (a-h,o-z)
      dimension   ql(meqn,1-mbc:mx+mbc)
      dimension   qr(meqn,1-mbc:mx+mbc)
      dimension    s(mwaves,1-mbc:mx+mbc)
      dimension wave(meqn,mwaves,1-mbc:mx+mbc)
      dimension amdq(meqn,1-mbc:mx+mbc)
      dimension apdq(meqn,1-mbc:mx+mbc)
      dimension auxl(maux,1-mbc:mx+mbc)
      dimension auxr(maux,1-mbc:mx+mbc)
c
c
      write(6,*) 'You must provide a valid Riemann solver'
      stop
c
      return
      end


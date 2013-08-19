c
c
c =========================================================
       subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      common /cqinit/ beta
c
c
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         q(1,i) = dexp(-beta * (xcell-0.3d0)**2)  
         q(2,i) = 0.d0
  150    continue
c
      return
      end


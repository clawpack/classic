
c
c
c =========================================================
      subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c =========================================================
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
c
c     do nothing... no source term
c
      return
      end


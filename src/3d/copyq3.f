c
c
c     ==================================================================
      subroutine copyq3(meqn,mbc,mx,my,mz,q1,q2)
c     ==================================================================
c
c     # copy the contents of q1 into q2
c
      implicit real*8(a-h,o-z)
      dimension q1(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 
     &             1-mbc:mz+mbc)
      dimension q2(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 
     &             1-mbc:mz+mbc)
c
c
      !$OMP parallel do collapse(4)
      do 10 k = 1-mbc, mz+mbc
         do 10 j = 1-mbc, my+mbc
            do 10 i = 1-mbc, mx+mbc
               do 10 m=1,meqn
                  q2(m,i,j,k) = q1(m,i,j,k)
 10            continue
c
      return
      end

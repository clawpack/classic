c
c
c =========================================================
      subroutine copyq2(meqn,mbc,mx,my,q1,q2)
c =========================================================
c
c     # copy the contents of q1 into q2
c
      implicit double precision (a-h,o-z)
      dimension q1(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension q2(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
c
c
      do 10 j = 1-mbc, my+mbc
        do 10 i = 1-mbc, mx+mbc
          do 10 m=1,meqn
             q2(m,i,j) = q1(m,i,j)
   10        continue
      return
      end


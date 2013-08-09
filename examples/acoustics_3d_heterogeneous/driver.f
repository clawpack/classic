      program driver
c
c  Generic driver routine for claw3
c
c
      implicit double precision (a-h,o-z)

c     # set parameters for maximum array sizes used in declarations
c     # these must be increased for larger problems.
c
c
      parameter (maxmx =    80)
      parameter (maxmy =    80)
      parameter (maxmz =    80)
      parameter (mwork = 2388876)


      parameter (mbc = 2)
      parameter (meqn = 4)
      parameter (mwaves = 2)
      parameter (maux = 2)

      dimension    q(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc)

      dimension  aux(maux, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc)

      dimension mthlim(mwaves)
      dimension work(mwork)
c
      call claw3ez(maxmx,maxmy,maxmz,meqn,mwaves,mbc,maux,mwork,mthlim,
     &           q,work,aux)

      stop
      end

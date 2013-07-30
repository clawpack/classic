
c
c
c     =====================================================
      subroutine bc3(meqn,mbc,mx,my,mz,xlower,
     &               ylower,zlower,dx,dy,dz,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw3
c
c     # At each boundary  k = 1 (xlower),  2 (xupper), 
c     #                       3 (ylower),  4 (yupper),
c     #                       5 (zlower),  6 (zupper):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  or 4'th (for k=5,6) component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my, 1:mz)
c     # to a layer of mbc ghost cells outside the region.
c
      implicit double precision (a-h,o-z)
      dimension    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
      dimension  aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 
     &               1-mbc:mz+mbc)
      dimension mthbc(6)

c
c
c-------------------------------------------------------
c     # left boundary (xlower):
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc3'
      stop
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 m=1,meqn
         do 115 ibc=1,mbc
            do 115 j = 1-mbc, my+mbc
               do 115 k = 1-mbc, mz+mbc
                  q(m,1-ibc,j,k) = q(m,1,j,k)
  115       continue
      go to 199

  120 continue
c     # periodic:  
      do 125 m=1,meqn
         do 125 ibc=1,mbc
            do 125 j = 1-mbc, my+mbc
               do 125 k = 1-mbc, mz+mbc
                  q(m,1-ibc,j,k) = q(m,mx+1-ibc,j,k)
  125       continue
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 ibc=1,mbc
            do 135 j = 1-mbc, my+mbc
               do 135 k = 1-mbc, mz+mbc
                  q(m,1-ibc,j,k) = q(m,ibc,j,k)
  135       continue
c     # negate the normal velocity:
      do 136 ibc=1,mbc
         do 136 j = 1-mbc, my+mbc
            do 136 k = 1-mbc, mz+mbc
               q(2,1-ibc,j,k) = -q(2,ibc,j,k)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary (xupper):
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc3'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 m=1,meqn
         do 215 ibc=1,mbc
            do 215 j = 1-mbc, my+mbc
               do 215 k = 1-mbc, mz+mbc
                  q(m,mx+ibc,j,k) = q(m,mx,j,k)
  215       continue
      go to 299

  220 continue
c     # periodic:  
      do 225 m=1,meqn
         do 225 ibc=1,mbc
            do 225 j = 1-mbc, my+mbc
               do 225 k = 1-mbc, mz+mbc
                  q(m,mx+ibc,j,k) = q(m,ibc,j,k)
  225       continue
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 ibc=1,mbc
            do 235 j = 1-mbc, my+mbc
               do 235 k = 1-mbc, mz+mbc
                  q(m,mx+ibc,j,k) = q(m,mx+1-ibc,j,k)
  235       continue
c     # negate the normal velocity:
      do 236 ibc=1,mbc
         do 236 j = 1-mbc, my+mbc
            do 236 k = 1-mbc, mz+mbc
               q(2,mx+ibc,j,k) = -q(2,mx+1-ibc,j,k)
  236    continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary (ylower):
c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc3'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 m=1,meqn
         do 315 jbc=1,mbc
            do 315 i = 1-mbc, mx+mbc
               do 315 k = 1-mbc, mz+mbc
                  q(m,i,1-jbc,k) = q(m,i,1,k)
  315       continue
      go to 399

  320 continue
c     # periodic:  
      do 325 m=1,meqn
         do 325 jbc=1,mbc
            do 325 i = 1-mbc, mx+mbc
               do 325 k = 1-mbc, mz+mbc
                  q(m,i,1-jbc,k) = q(m,i,my+1-jbc,k)
  325       continue
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 m=1,meqn
         do 335 jbc=1,mbc
            do 335 i = 1-mbc, mx+mbc
               do 335 k = 1-mbc, mz+mbc
                  q(m,i,1-jbc,k) = q(m,i,jbc,k)
  335       continue
c     # negate the normal velocity:
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, mx+mbc
            do 336 k = 1-mbc, mz+mbc
               q(3,i,1-jbc,k) = -q(3,i,jbc,k)
  336    continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary (yupper):
c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc3'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 m=1,meqn
         do 415 jbc=1,mbc
            do 415 i = 1-mbc, mx+mbc
               do 415 k = 1-mbc, mz+mbc
                  q(m,i,my+jbc,k) = q(m,i,my,k)
  415       continue
      go to 499

  420 continue
c     # periodic:  
      do 425 m=1,meqn
         do 425 jbc=1,mbc
            do 425 i = 1-mbc, mx+mbc
               do 425 k = 1-mbc, mz+mbc
                  q(m,i,my+jbc,k) = q(m,i,jbc,k)
  425       continue
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 m=1,meqn
         do 435 jbc=1,mbc
            do 435 i = 1-mbc, mx+mbc
               do 435 k = 1-mbc, mz+mbc
                  q(m,i,my+jbc,k) = q(m,i,my+1-jbc,k)
  435       continue
c     # negate the normal velocity:
      do 436 jbc=1,mbc
         do 436 i = 1-mbc, mx+mbc
            do 436 k = 1-mbc, mz+mbc
               q(3,i,my+jbc,k) = -q(3,i,my+1-jbc,k)
  436    continue
      go to 499

  499 continue

c
c-------------------------------------------------------
c     # boundary (zlower):
c-------------------------------------------------------
      go to (500,510,520,530) mthbc(5)+1
c
  500 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(5)=0 and no BCs specified in bc3'
      stop
      go to 599
c
  510 continue
c     # zero-order extrapolation:
      do 515 m=1,meqn
         do 515 kbc=1,mbc
            do 515 i = 1-mbc, mx+mbc
               do 515 j = 1-mbc, my+mbc
                  q(m,i,j,1-kbc) = q(m,i,j,1)
  515       continue
      go to 599

  520 continue
c     # periodic:  
      do 525 m=1,meqn
         do 525 kbc=1,mbc
            do 525 i = 1-mbc, mx+mbc
               do 525 j = 1-mbc, my+mbc
                  q(m,i,j,1-kbc) = q(m,i,j,mz+1-kbc)
  525       continue
      go to 599

  530 continue
c     # solid wall (assumes 4'rd component is velocity or momentum in y):
      do 535 m=1,meqn
         do 535 kbc=1,mbc
            do 535 i = 1-mbc, mx+mbc
               do 535 j = 1-mbc, my+mbc
                  q(m,i,j,1-kbc) = q(m,i,j,kbc)
  535       continue
c     # negate the normal velocity:
      do 536 kbc=1,mbc
         do 536 i = 1-mbc, mx+mbc
            do 536 j = 1-mbc, my+mbc
               q(4,i,j,1-kbc) = -q(4,i,j,kbc)
  536    continue
      go to 599

  599 continue
c
c-------------------------------------------------------
c     # boundary (zupper):
c-------------------------------------------------------
      go to (600,610,620,630) mthbc(6)+1
c
  600 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(6)=0 and no BCs specified in bc3'
      stop
      go to 699

  610 continue
c     # zero-order extrapolation:
      do 615 m=1,meqn
         do 615 kbc=1,mbc
            do 615 i = 1-mbc, mx+mbc
               do 615 j = 1-mbc, my+mbc
                  q(m,i,j,mz+kbc) = q(m,i,j,mz)
  615       continue
      go to 699

  620 continue
c     # periodic:  
      do 625 m=1,meqn
         do 625 kbc=1,mbc
            do 625 i = 1-mbc, mx+mbc
               do 625 j = 1-mbc, my+mbc
                  q(m,i,j,mz+kbc) = q(m,i,j,kbc)
  625       continue
      go to 699

  630 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 635 m=1,meqn
         do 635 kbc=1,mbc
            do 635 i = 1-mbc, mx+mbc
               do 635 j = 1-mbc, my+mbc
                  q(m,i,j,mz+kbc) = q(m,i,j,mz+1-kbc)
  635       continue
c     # negate the normal velocity:
      do 636 kbc=1,mbc
         do 636 i = 1-mbc, mx+mbc
            do 636 j = 1-mbc, my+mbc
               q(4,i,j,mz+kbc) = -q(4,i,j,mz+1-kbc)
  636    continue
      go to 699

  699 continue

      return
      end


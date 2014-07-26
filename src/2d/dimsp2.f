c
c
c     ==========================================================
      subroutine dimsp2(maxm,meqn,mwaves,maux,mbc,mx,my,
     &                  qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &                  cflv,qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                  aux1,aux2,aux3,work,mwork,use_fwaves,rpn2,rpt2)
c     ==========================================================
c
c     # Take one time step, updating q, using dimensional
c     # splitting. Two choices are available:
c     #
c     # method(3) = -1   gives Godunov splitting:
c     #    time step dt in x-direction
c     #    time step dt in y-direction
c
c     # method(3) = -2   gives Strang splitting
c     #    time step dt/2 in x-direction
c     #    time step dt   in y-direction
c     #    time step dt/2 in x-direction
c
c     # Godunov splitting is recommended over Strang splitting normally
c     # since it typically works as well, is faster, and boundary
c     # conditions are handled properly.
c
      implicit double precision (a-h,o-z)
      external rpn2,rpt2
      dimension qold(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension qnew(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension  q1d(meqn,1-mbc:maxm+mbc)
      dimension cflv(4)
      dimension qadd(meqn,1-mbc:maxm+mbc)
      dimension fadd(meqn,1-mbc:maxm+mbc)
      dimension gadd(meqn,1-mbc:maxm+mbc,2)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux1(maux,1-mbc:maxm+mbc)
      dimension aux2(maux,1-mbc:maxm+mbc)
      dimension aux3(maux,1-mbc:maxm+mbc)

      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)

      logical :: use_fwaves

c
c     # If method(3) = -1, take a full time step in x.
c     # If method(3) = -2, take a half time step in x.
c
      dt2 = dt/2.d0
c
      if( method(3) .eq. -2 )then
          call step2ds(maxm,meqn,mwaves,maux,mbc,mx,my,
     &                 qold,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,1,use_fwaves,rpn2,rpt2)
      else
          call step2ds(maxm,meqn,mwaves,maux,mbc,mx,my,
     &                 qold,qnew,aux,dx,dy,dt,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,1,use_fwaves,rpn2,rpt2)
      endif
c
      if (cflx .gt. cflv(1)) then
c        # Abort if the Courant number was too large in x-sweep
         cfl = cflx
         return
         endif
c
c     # Take full step in y-direction
c
      call step2ds(maxm,meqn,mwaves,maux,mbc,mx,my,
     &             qnew,qnew,aux,dx,dy,dt,method,mthlim,cfly,
     &             qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &             aux1,aux2,aux3,work,mwork,2,use_fwaves,rpn2,rpt2)
c
      cfl = dmax1(cflx,cfly)
c
c     # Finally, take a half time step in the x-direction
c     # if Strang splitting is used.  NOTE: boundary conditions may
c     # not be set properly for this sweep.
c
      if( method(3) .eq. -2 )then
         if (cfly .gt. cflv(1)) then
c           # Abort if the Courant number was too large in y-sweep
            cfl = cfly
            return
            endif
          call step2ds(maxm,meqn,mwaves,maux,mbc,mx,my,
     &                 qnew,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,1,use_fwaves,rpn2,rpt2)
          cfl = dmax1(cfl,cflx)
      endif
c
      return
      end


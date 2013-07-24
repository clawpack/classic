c
c
c
c     =================================================================
      subroutine claw3ez(mx_in,my_in,mz_in,meqn,mbc,maux,mwork,
     &                   q,work,aux)
c     =================================================================
c
c     An easy-to-use clawpack driver routine for simple applications
c     Documentation is available at
c                 http://www.amath.washington.edu/~claw/doc.html
c
c     Author: Randall J. LeVeque
c     Version of August, 2002 --  CLAWPACK Version 4.1
c
c
      implicit double precision (a-h,o-z)
      external bc3,rpn3,rpt3,rptt3,src3,b4step3

      dimension    q(meqn, 1-mbc:mx_in+mbc, 1-mbc:my_in+mbc,
     &               1-mbc:mz_in+mbc)
      dimension  aux(maux, 1-mbc:mx_in+mbc, 1-mbc:my_in+mbc,
     &               1-mbc:mz_in+mbc)
      dimension work(mwork)
c
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(6)
      dimension tout(100)
      logical rest
c
      common /restrt_block/ tinitial, iframe
c
      integer, dimension(:), allocatable :: mthlim

      open(55,file='claw3ez.data',status='old',form='formatted')
      open(10,file='fort.info',status='unknown',form='formatted')
c
c
c     # Read the input in standard form from claw2ez.data:

c     domain variables
      read(55,*) mx
      read(55,*) my
      read(55,*) mz

c     i/o variables
      read(55,*) nout
      read(55,*) outstyle
      if (outstyle.eq.1) then
          read(55,*) tfinal
          nstepout = 1
        elseif (outstyle.eq.2) then
          read(55,*) (tout(i), i=1,nout)
          nstepout = 1
        elseif (outstyle.eq.3) then
          read(55,*) nstepout, nstop
          nout = nstop
        endif


c     timestepping variables
      read(55,*) dtv(1)
      read(55,*) dtv(2)
      read(55,*) cflv(1)
      read(55,*) cflv(2)
      read(55,*) nv(1)
c


c     # input parameters for clawpack routines
      read(55,*) method(1)
      read(55,*) method(2)
      read(55,*) method(3)
      read(55,*) method(4)
      read(55,*) method(5)
      read(55,*) method(6)
      read(55,*) method(7)

      read(55,*) meqn1
      read(55,*) mwaves
      allocate(mthlim(mwaves))
      read(55,*) (mthlim(mw), mw=1,mwaves)

      read(55,*) t0
      read(55,*) xlower
      read(55,*) xupper
      read(55,*) ylower
      read(55,*) yupper
      read(55,*) zlower
      read(55,*) zupper
c
      read(55,*) mbc1
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)
      read(55,*) mthbc(3)
      read(55,*) mthbc(4)
      read(55,*) mthbc(5)
      read(55,*) mthbc(6)

c     # check to see if we are restarting:
      rest = .false.
c     # The next two lines may not exist in old versions of claw3ez.data.
c     # Jump over the second read statement if the 1st finds an EOF:
      read(55,*,end=199,err=199) rest
      read(55,*) iframe   !# restart from data in fort.qN file, N=iframe
 199  continue


      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(1) and mthbc(2) BOTH be set to 2'
         stop
         endif

      if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &    (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(3) and mthbc(4) BOTH be set to 2'
         stop
         endif

      if ((mthbc(5).eq.2 .and. mthbc(6).ne.2) .or.
     &    (mthbc(6).eq.2 .and. mthbc(5).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(5) and mthbc(6) BOTH be set to 2'
         stop
         endif

c     # These values were passed in, but check for consistency:
c
      if (method(7) .ne. maux) then
         write(6,*) '*** ERROR ***  method(7) should equal maux'
         stop
         endif
      if (meqn1 .ne. meqn) then
         write(6,*) '*** ERROR ***  meqn set wrong in input or driver'
         stop
         endif
      if (mbc1 .ne. mbc) then
         write(6,*) '*** ERROR ***  mbc set wrong in input or driver'
         stop
         endif
c
c     # check that enough storage has been allocated:
c
      if (method(5).lt.2) then
          narray = 1   !# only need one qwork array
        else
          narray = 2   !# need two qwork arrays for Strang splitting
        endif

      maxm = max0(mx_in, my_in, mz_in)
      mwork1 = (maxm+2*mbc)*(46*meqn + mwaves + meqn*mwaves
     &                      + 9*maux + 3)
     &          + narray * (mx_in + 2*mbc) * (my_in + 2*mbc)
     &                   * (mz_in + 2*mbc) * meqn
c
c
      if (mx.gt.mx_in .or. my.gt.my_in .or. mz.gt.mz_in .or.
     &    mwork.lt.mwork1) then
c        # insufficient storage
         maxmx1 = max0(mx,mx_in)
         maxmy1 = max0(my,my_in)
         maxmz1 = max0(mz,mz_in)
         maxm1 = max0(maxmx1,maxmy1,maxmz1)

         mwork1 = (maxm1+2*mbc)*(46*meqn + mwaves + meqn*mwaves
     &                      + 9*maux + 3)
     &          + narray * (mx_in + 2*mbc) * (my_in + 2*mbc)
     &                   * (mz_in + 2*mbc) * meqn

         ! This error should never happen with the new dynamic
         ! allocation setup
         write(6,*) ' '
         write(6,*) '*** ERROR *** Insufficient storage allocated'
         write(6,*) 'Recompile after increasing values in driver.f:'
         write(6,611) maxmx1
         write(6,612) maxmy1
         write(6,613) maxmz1
         write(6,614) mwork1
 611     format(/,'parameter (mx_in = ',i5,')')
 612     format('parameter (my_in = ',i5,')')
 613     format('parameter (maxmz = ',i5,')')
 614     format('parameter (mwork = ',i9,')',/)
         stop
         endif

      call chkmth(method,info)
      if( info .eq. 6) stop
c
c
      write(6,*) 'running...'
      write(6,*) ' '
c
c     # grid spacing
      dx = (xupper - xlower) / float(mx)
      dy = (yupper - ylower) / float(my)
      dz = (zupper - zlower) / float(mz)
c


c     # time increments between outputing solution:
      if (outstyle .eq. 1) then
         dtout = (tfinal - t0)/float(nout)
      endif
c
c
c     # call user's routine setprob to set any specific parameters
c     # or other initialization required.
c
      call setprob
c
c     # set aux array:
c
      if (maux .gt. 0)  then
         call setaux(mbc,mx,my,mz,xlower,ylower,
     &               zlower,dx,dy,dz,maux,aux,t0)
         endif
c
c     # set initial conditions:

      if (rest) then
          call restart(meqn,mbc,mx,my,mz,
     &          xlower,ylower,zlower,dx,dy,dz,q)
          t0 = tinitial
        else
          call qinit(meqn,mbc,mx,my,mz,xlower,
     &           ylower,zlower,dx,dy,dz,q,maux,aux)
          iframe = 0
        endif
c
c
c
      if (.not. rest) then
c        # output initial data
         call out3(meqn,mbc,mx,my,mz,xlower,ylower,
     &          zlower,dx,dy,dz,q,t0,iframe,aux,maux)
         write(6,601) iframe, t0
         endif
c
c     ----------
c     Main loop:
c     ----------
c
      tend = t0
      n0   = iframe*nstepout + 1
      do 100 n=n0,nout
         tstart = tend
         if (outstyle .eq. 1)  tend = tstart + dtout
         if (outstyle .eq. 2)  tend = tout(n)
         if (outstyle .eq. 3)  tend = tstart - 1.d0  !# single-step mode
c
         call claw3(meqn,mwaves,mbc,mx,my,mz,maux,
     &           q,aux,xlower,ylower,zlower,dx,dy,dz,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,info,bc3,rpn3,rpt3,rptt3,src3,b4step3)
c
c        # check to see if an error occured:
         if (info .ne. 0) then
            write(6,*) 'claw3ez aborting: Error return from claw3',
     &                 ' with info =',info
            go to 999
            endif
c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c        # output solution at this time
c        ------------------------------
c
c        # if outstyle=1 or 2, then nstepout=1 and we output every time
c        # we reach this point, since claw1 was called for the entire time
c        # increment between outputs.
c
c        # if outstyle=3 then we only output if we have taken nstepout
c        # time steps since the last output.

c        # iframe is the frame number used to form file names in out3
         if (mod(n,nstepout) .eq. 0) then
            iframe = iframe + 1
            call out3(meqn,mbc,mx,my,mz,xlower,ylower,
     &            zlower,dx,dy,dz,q,tend,iframe,aux,maux)
            write(6,601) iframe,tend
            write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
     &           cflv(3),cflv(4),nv(2)
            endif

c
c        # formats for writing out information about this call to claw:

  601    format('CLAW3EZ: Frame ',i4,
     &         ' matlab plot files done at time t =',
     &         d12.4,/)
c
 1010    format('tend =',d15.4,/,
     &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
     &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
     &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
c
  100    continue
c
  999 continue
c
      deallocate(mthlim)
      return
      end


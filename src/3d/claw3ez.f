c
c
c
c     =================================================================
      subroutine claw3ez    ! No arguments
c     =================================================================
c
c     An easy-to-use clawpack driver routine for simple applications
c     Documentation is available at
c                 http://www.amath.washington.edu/~claw/doc.html
c
c     Author: Randall J. LeVeque, Grady I. Lemoine
c
c     Modified July 2013 to use dynamic allocation, new input file
c     format.  Allocation of all the large arrays (aux, q, work) is
c     delayed until the last possible moment in case any of the user
c     code called before them is a memory hog.
c

      implicit double precision (a-h,o-z)

c     F77 way to get OpenMP function declarations
c$    include 'omp_lib.h'

      external bc3,rpn3,rpt3,rptt3,src3,b4step3

      double precision, dimension(:,:,:,:), allocatable :: q, aux
      double precision, dimension(:), allocatable :: work

c
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(6)
      double precision, dimension(:), allocatable :: tout
      logical rest
c
      common /restrt_block/ tinitial, iframe
c
      integer, dimension(:), allocatable :: mthlim

      logical :: outaux_init_only, output_t0, use_fwaves
      logical :: outaux_always, dt_variable
      integer :: output_format, dimensional_split, outstyle
      integer, dimension(:), allocatable :: iout_q, iout_aux
      integer :: allocate_status

      call opendatafile(55, 'claw.data')
      open(10,file='fort.info',status='unknown',form='formatted')
c
c
c     # Read the input in standard form from claw2ez.data:

      read(55,*) ndim    ! Not actually used, but the modern Python setup writes it
      read(55,*) xlower, ylower, zlower
      read(55,*) xupper, yupper, zupper
      read(55,*) mx, my, mz
      read(55,*) meqn
      read(55,*) mwaves
      read(55,*) maux
      read(55,*) t0

      read(55,*) outstyle
      if (outstyle.eq.1) then
         read(55,*) nout
         read(55,*) tfinal
         read(55,*) output_t0    ! Not currently used
         nstepout = 1
      else if (outstyle.eq.2) then
         read(55,*) nout
         allocate(tout(nout), stat=allocate_status)
         if (allocate_status .ne. 0) then
            print *, '*** Error allocating tout array; exiting claw3ez'
            go to 900
         end if
         read(55,*) (tout(i), i=1,nout)
         nstepout = 1
      else if (outstyle.eq.3) then
         read(55,*) nstepout
         read(55,*) nstop
         read(55,*) output_t0
         nout = nstop
      else
         print *, '*** Unrecognized output style ', outstyle
         print *, '*** Exiting claw3ez'
         go to 900
      end if

      read(55,*) output_format    ! Not used yet
      ! These iout variables are not currently used, but hang onto them
      ! anyway in case somebody wants to use them at a future date.  The
      ! same goes for outaux_init_only.
      allocate(iout_q(meqn), stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating iout_q array; exiting claw3ez'
         go to 900    ! Exception handling, old school style
      end if
      read(55,*) (iout_q(i), i = 1, meqn)
      if (maux > 0) then
         allocate(iout_aux(maux), stat=allocate_status)
         if (allocate_status .ne. 0) then
            print *, '*** Error allocating iout_aux array;',
     &               ' exiting claw3ez'
            go to 900
         end if
         read(55,*) (iout_aux(i), i = 1, maux)
         read(55,*) outaux_init_only
         ! Not implementing selective output of aux fields yet
         if (any(iout_aux .ne. 0)) then
            outaux_always = .not. outaux_init_only
         else
            outaux_always = .false.
            outaux_init_only = .false.
         end if
      else
         outaux_always = .false.
         outaux_init_only = .false.    ! Just for the sake of initialization
      end if

      read(55,*) dtv(1)     ! Initial dt
      read(55,*) dtv(2)     ! Max dt
      read(55,*) cflv(1)    ! Max CFL number
      read(55,*) cflv(2)    ! Desired CFL number
      read(55,*) nv(1)      ! Maximum number of steps

      read(55,*) dt_variable  ! Variable or fixed dt
      if (dt_variable) then
         method(1) = 1
      else
         method(1) = 0
      end if
      read(55,*) method(2)    ! Order
      read(55,*) method(3)    ! Transverse propagation style
      read(55,*) dimensional_split    ! Whether to use dimensional splitting

      ! Translate new-style Python specification of transverse
      ! order/dimensional splitting into something the older code
      ! understands.
      if (dimensional_split > 0) method(3) = -1

      read(55,*) method(4)    ! Verbosity
      read(55,*) method(5)    ! Source term splitting style
      read(55,*) method(6)    ! Index into aux for capacity function
      method(7) = maux    ! Number of aux variables

      read(55,*) use_fwaves

      allocate(mthlim(mwaves), stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating mthlim array; exiting claw3ez'
         go to 900
      end if
      read(55,*) (mthlim(i), i = 1, mwaves)

      read(55,*) mbc
      read(55,*) mthbc(1), mthbc(3), mthbc(5)
      read(55,*) mthbc(2), mthbc(4), mthbc(6)

      rest = .false.
      read(55, *, err=199, end=199) rest
      if (rest) then
         print *, 'Doing a restart run'
         print *, 'Attempting to read restart frame number'
         print *, 'You may need to hand-edit claw.data'
         read(55,*) iframe      ! restart from data in fort.qN file, N=iframe
         print *, 'Planning restart from frame ', iframe
      end if
 199  continue

      close(unit=55)


      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(1) and mthbc(2) BOTH be set to 2'
         go to 900
         endif

      if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &    (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(3) and mthbc(4) BOTH be set to 2'
         go to 900
         endif

      if ((mthbc(5).eq.2 .and. mthbc(6).ne.2) .or.
     &    (mthbc(6).eq.2 .and. mthbc(5).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(5) and mthbc(6) BOTH be set to 2'
         go to 900
         endif

c
c     # Figure out how much work space is necessary
c
      if (method(5).lt.2) then
          narray = 1   !# only need one qwork array
        else
          narray = 2   !# need two qwork arrays for Strang splitting
        endif

      maxm = max0(mx, my, mz)
      nthreads=1    ! Serial
c$omp parallel
c$omp single
c$    nthreads = omp_get_num_threads()
c$    maxthreads = omp_get_max_threads()
c$omp end single
c$omp end parallel
      mwork = (maxm+2*mbc)*(46*meqn + mwaves + meqn*mwaves
     &                     + 9*maux + 3)*nthreads
     &         + narray * (mx + 2*mbc) * (my + 2*mbc)
     &                  * (mz + 2*mbc) * meqn

c
c

      call chkmth(method,info)
      if( info .eq. 6) then
         print *, '*** Error in method data; exiting claw3ez'
         go to 900
      end if
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

c     # Allocate aux
      if (maux > 0) then
         allocate(aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc),
     &            stat=allocate_status)
      else
         ! Allocate a dummy array so that dereferencing it doesn't cause
         ! a segfault.  This may not be necessary.
         allocate(aux(1,1,1,1), stat=allocate_status)
      end if
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating aux array; exiting claw3ez'
         go to 900
      end if
c
c     # set aux array:
c
      if (maux .gt. 0)  then
         call setaux(mbc,mx,my,mz,xlower,ylower,
     &               zlower,dx,dy,dz,maux,aux)
         endif

c     # Allocate q
      allocate(q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc),
     &         stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error allocating q array; exiting claw3ez'
         go to 900
      end if

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
     &          zlower,dx,dy,dz,q,t0,iframe,aux,maux,
     &          outaux_init_only .or. outaux_always)
         write(6,601) iframe, t0
         endif

      ! Allocate work array
      allocate(work(mwork), stat=allocate_status)
      if (allocate_status .ne. 0) then
         print *, '*** Error (nthreads=',nthreads,
     &        ') allocating work array; exiting claw3ez'
         go to 900
      end if

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
     &           work,mwork,use_fwaves,info,bc3,rpn3,rpt3,rptt3,src3,
     &           b4step3)
c
c        # check to see if an error occured:
         if (info .ne. 0) then
            write(6,*) 'claw3ez aborting: Error return from claw3',
     &                 ' with info =',info
            go to 900
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
     &            zlower,dx,dy,dz,q,tend,iframe,aux,maux,
     &            outaux_always)
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
 900  continue
c
      ! Destructor, old school style
      if (allocated(work))     deallocate(work)
      if (allocated(q))        deallocate(q)
      if (allocated(aux))      deallocate(aux)
      if (allocated(mthlim))   deallocate(mthlim)
      if (allocated(iout_aux)) deallocate(iout_aux)
      if (allocated(iout_q))   deallocate(iout_q)
      if (allocated(tout))     deallocate(tout)

      return
      end


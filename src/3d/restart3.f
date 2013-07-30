c
c     =====================================================
      subroutine restart(meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
c
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &      1-mbc:mz+mbc)
      character*10 fname1, fname2
      logical outt0
      common /restrt_block/ t0, iframe

      iunit = 16

c     # first create the file name and open file
      fname1 = 'fort.q'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
      open(iunit,file=fname1)

c     # Read grid parameters.
      read(iunit,*) igrid
      read(iunit,*) level
      read(iunit,*) mx_in
      read(iunit,*) my_in
      read(iunit,*) mz_in
      read(iunit,*) xlow_in
      read(iunit,*) ylow_in
      read(iunit,*) zlow_in
      read(iunit,*) dx_in
      read(iunit,*) dy_in
      read(iunit,*) dz_in

c     # Test for compatibility of grid resolution.
      if (mx_in .ne. mx .or. my_in .ne. my .or. mz_in .ne. mz) then
         stop 'rstart.f : data not compatible'
      endif

c     # Read variables in from old fort.qXXXX file.
      do k = 1,mz
         do j = 1,my
            do i = 1,mx
               read(iunit,*) (q(m,i,j,k),m=1,meqn)
            enddo
         enddo
      enddo
      close(iunit)

c     # Read initial time in from fort.tXXXX file.      
      fname2 = 'fort.t'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
      open(iunit,file=fname2)
      read(iunit,*) t0
      close(iunit)

      write(*,*) 'Restarting from old output file ', fname1
      write(*,*) 'Simulation will be restarted at time t = ', t0
      write(*,*) 'Inital condition will not be output to a matlab ',
     &     'plot file'
      write(*,*) 

      return
      end


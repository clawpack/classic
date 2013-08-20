c
c     ==================================================================
      subroutine chkmth(method,info)
c     ==================================================================
c
c     # Checks whether the method parameters are correct or not.
c     # Note that method(3) < 0 yields dimensional splitting.
c
      implicit double precision (a-h,o-z)
c
      dimension method(7)
c
      info = 0
c
      if(method(1) .ne. 0 .and. method(1) .ne. 1)then
         write(6,*) ' '
         write(6,*) 'CLAW3 error .... method(1) should be:'
         write(6,*) '   0 - fixed time steps'
         write(6,*) '   1 - variable time steps'
         write(6,*) ' '
         info = 6 
      endif
c
      if(method(2) .ne. 1 .and. method(2) .ne. 2)then
         write(6,*) ' '
         write(6,*) 'CLAW3 error .... method(2) should be:'
         write(6,*) '   1 - first order method'
         write(6,*) '   2 - second order method'
         write(6,*) ' '
         info = 6 
      endif
c
      if(method(2) .eq. 1 .and. method(3) .ge. 0)then
         if(method(3) .ne. 0 .and. method(3) .ne. 10 .and.
     &      method(3) .ne. 11)then
            write(6,*) ' '
            write(6,*) 'CLAW3 error .... when method(2) = 1,'
            write(6,*) 'method(3) should be:'
            write(6,*) '   0  - donor cell'
            write(6,*) '   10 - 2D wave propagation of increment waves'
            write(6,*) '   11 - corner transport upwind'
            write(6,*) ' '
            info = 6
         endif
      endif       
c
      if(method(2) .eq. 2 .and. method(3) .ge. 0)then
         if(method(3) .ne. 0  .and. method(3) .ne. 10 .and.
     &      method(3) .ne. 11 .and. method(3) .ne. 20 .and.
     &      method(3) .ne. 21 .and. method(3) .ne. 22)then
            write(6,*) ' '
            write(6,*) 'CLAW3 error .... when method(2) = 2,'
            write(6,*) 'method(3) should be:'
            write(6,*) '   0  - donor cell'
            write(6,*) '   10 - 2D wave propagation of increment waves'
            write(6,*) '        and 1D propagation of correction waves'
            write(6,*) '        UNCONDITIONALLY UNSTABLE'
            write(6,*) '   11 - corner transport upwind with 2D'
            write(6,*) '        propagation of the correction waves'
            write(6,*) '        UNCONDITIONALLY UNSTABLE'
            write(6,*) '   20 - 2D propagation of both increment and'
            write(6,*) '        correction waves.'
            write(6,*) '   21 - corner transport upwind with 2D '
            write(6,*) '        propagation of the correction waves.'
            write(6,*) '   22 - Full 3D wave propagation method.'
            write(6,*) ' '
            info = 6
         endif
      endif       
c
      return
      end

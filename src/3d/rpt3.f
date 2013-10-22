c
c
c     ==================================================================
      subroutine rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,asdq,
     &                  bmasdq,bpasdq)
c     ==================================================================
c
c     # Riemann solver in the transverse direction.
c     # This is a dummy routine that returns zeros and is only intended
c     # to illustrate the format of this routine.  See various example
c     # directories for better examples.

c     # This dummy routine can be used if transverse solves are not being
c     # used, i.e. if method(3) = 0.
c     #
c     # On input,
c
c     #    ql,qr is the data along some one-dimensional slice, as in rpn3
c     #         This slice is
c     #             in the x-direction if ixyz=1,
c     #             in the y-direction if ixyz=2, or
c     #             in the z-direction if ixyz=3.
c     #    asdq is an array of flux differences (A^*\Dq).
c     #         asdq(:,i) is the flux difference propagating away from
c     #         the interface between cells i-1 and i.
c     #    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
c
c     #    ixyz indicates the direction of the original Riemann solve,
c     #         called the x-like direction in the table below:
c
c     #               x-like direction   y-like direction   z-like direction
c     #      ixyz=1:        x                  y                  z         
c     #      ixyz=2:        y                  z                  x         
c     #      ixyz=3:        z                  x                  y         
c
c     #    icoor indicates direction in which the transverse solve should 
c     #         be performed.
c     #      icoor=2: split in the y-like direction.
c     #      icoor=3: split in the z-like direction.
c
c     #    For example,
c     #      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
c     #                        bmasdq = B^-A^*\Dq,
c     #                        bpasdq = B^+A^*\Dq.
c     #
c     #      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into 
c     #                        bmasdq = C^-B^*\Dq,
c     #                        bpasdq = C^+B^*\Dq.
c
c     #    The parameter imp is generally needed only if aux
c     #    arrays are being used, in order to access the appropriate
c     #    variable coefficients.
c
c
c
      implicit real*8(a-h,o-z)
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
      dimension   aux1(maux, 1-mbc:maxm+mbc, 3)
      dimension   aux2(maux, 1-mbc:maxm+mbc, 3)
      dimension   aux3(maux, 1-mbc:maxm+mbc, 3)
c
      do i = 2-mbc, mx+mbc
         do m=1,meqn
            bmasdq(m,i) = 0.d0
            bpasdq(m,i) = 0.d0
         enddo
      enddo
c
c
      return
      end

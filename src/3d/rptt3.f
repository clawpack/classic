c
c
c     ==================================================================
      subroutine rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,bsasdq,
     &                  cmbsasdq,cpbsasdq)
c     ==================================================================
c
c     # Double transverse Riemann solver.
c     # This is a dummy routine that returns zeros and is only intended
c     # to illustrate the format of this routine.  See various example
c     # directories for better examples.

c     # This dummy routine can be used if double transverse solves are not being
c     # used, i.e. if method(3) = 0, 10, or 20.
c     #
c     # On input,
c
c     #    ql,qr is the data along some one-dimensional slice, as in rpn3
c     #         This slice is
c     #             in the x-direction if ixyz=1,
c     #             in the y-direction if ixyz=2, or
c     #             in the z-direction if ixyz=3.
c
c     #    bsasdq is an array of flux differences that result from a
c     #         transverse splitting (a previous call to rpt3).  
c     #         This stands for B^* A^* \Dq but could represent any of 
c     #         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
c     #         and icoor (see below).
c     #         Moreover, each * represents either + or -, as specified by
c     #         imp and impt.
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
c     #        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be 
c     #                        split in z into 
c     #                           cmbsasdq = C^-B^*A^*\Dq,
c     #                           cpbsasdq = C^+B^*A^*\Dq.
c     #
c     #        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
c     #                        split in x into 
c     #                           cmbsasdq = A^-C^*B^*\Dq,
c     #                           cpbsasdq = A^+C^*B^*\Dq.
c
c     #    The parameters imp and impt are generally needed only if aux
c     #    arrays are being used, in order to access the appropriate
c     #    variable coefficients.

c
      implicit real*8(a-h,o-z)
      dimension       ql(meqn, 1-mbc:maxm+mbc)
      dimension       qr(meqn, 1-mbc:maxm+mbc)
      dimension   bsasdq(meqn, 1-mbc:maxm+mbc)
      dimension cmbsasdq(meqn, 1-mbc:maxm+mbc)
      dimension cpbsasdq(meqn, 1-mbc:maxm+mbc)
      dimension     aux1(maux, 1-mbc:maxm+mbc, 3)
      dimension     aux2(maux, 1-mbc:maxm+mbc, 3)
      dimension     aux3(maux, 1-mbc:maxm+mbc, 3)

      do i = 2-mbc, mx+mbc
         do m=1,meqn
            cmbsasdq(m,i) = 0.d0
            cpbsasdq(m,i) = 0.d0
         enddo
      enddo
c
      return
      end

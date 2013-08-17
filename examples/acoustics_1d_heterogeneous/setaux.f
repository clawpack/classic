c     ============================================
      subroutine setaux(mbc,mx,xlower,dx,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # variable coefficient acoustics
c     #  aux(i,1) = density rho in i'th cell
c     #  aux(i,2) = sound speed c in i'th cell
c
c     # Piecewise constant medium with single interface at x=0
c     # Density and sound speed to left and right are set in setprob.f
c
c     
      implicit double precision (a-h,o-z)
      dimension aux(maux, 1-mbc:mx+mbc)
      common /comaux/ rhol,cl,rhor,cr

      open(unit=31,file='fort.aux',status='unknown',form='formatted')
c

       do i=1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
          if (xcell .lt. 0.0d0) then
              aux(1,i) = rhol
              aux(2,i) = cl
	    else
              aux(1,i) = rhor
              aux(2,i) = cr
            endif
	  enddo


	do i=1,mx
          write(31,701) aux(1,i), aux(2,i)
  701     format(2e16.6)
          enddo

       close(unit=31)
c
       return
       end

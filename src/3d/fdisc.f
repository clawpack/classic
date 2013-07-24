c
c
c
c     =================================================
      function fdisc(x,y,z)
c     =================================================
c
c     # For computing cell averages for initial data or coefficients that
c     # have a discontinuity along some surface in 3d.  

c     # fdisc should be negative to the "left" of the surface and 
c     # positive to the "right".

c     # The cellave routine can then be used to compute the fraction wl of
c     # a grid cell that lies to the "left" of the surface.

c     # Sample code for the case where the surface is the unit sphere:

      implicit double precision (a-h,o-z)

      fdisc = x**2 + y**2 + z**2 - 1.d0
c
      return
      end

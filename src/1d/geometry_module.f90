! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module geometry_module

    use precision_module

    implicit none
    
    type geometry_type
        real(kind=dp) :: x, dx, t, dt
    end type geometry_type

end module geometry_module

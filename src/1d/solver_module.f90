! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module solver_module

    implicit none
    
    ! Solver type parameters - defaults to 1.d0
    integer, parameter :: F_TYPE = kind(1.d0)

    ! Solver type declaration
    type solver_type
        
        ! Status of solver

        ! Memory storage for this solver
        real(kind=F_TYPE), pointer :: f(:,:)


    end type solver_type

contains

    function new_solver() result(solver)

        implicit none

        type(solver_type) :: solver
        
    end function new_solver

end module solver_module

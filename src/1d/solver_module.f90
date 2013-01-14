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
    integer, parameter :: ASDQ_TYPE = kind(1.d0)
    integer, parameter :: WAVE_TYPE = kind(1.d0)
    integer, parameter :: S_TYPE = kind(1.d0)

    ! Solver type declaration
    type solver_type
        
        ! Status of solver
        real(kind=8) :: cfl, dt

        ! Solver parameters
        integer :: order, transverse_waves, dimensional_split, source_split
        integer :: verbosity
        integer, allocatable :: limiter(:)
        real(kind=8) :: dt_max, cfl_max, cfl_desired
        logical :: use_fwaves, dt_variable

        ! Memory storage for this solver
        real(kind=F_TYPE), pointer :: f(:,:)
        real(kind=ASDQ_TYPE), pointer :: apdq(:,:), amdq(:,:)
        real(kind=WAVE_TYPE), pointer :: wave(:,:,:)
        real(kind=S_TYPE), pointer :: s(:,:)

    end type solver_type

contains

    function new_solver() result(solver)

        implicit none

        type(solver_type) :: solver
        
    end function new_solver

end module solver_module

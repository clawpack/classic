! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module solution_module

    implicit none
    
    ! Array types - Default is double precision as defined by 1.d0
    integer, parameter :: Q_TYPE = kind(1.d0)
    integer, parameter :: AUX_TYPE = kind(1.d0)

    type solution_type

        ! Solution extents
        integer :: num_cells, num_eqn, num_aux, num_ghost
        real(kind=8) :: t

        ! Aux array descriptors
        integer :: capa_index

        ! Solution arrays
        real(kind=Q_TYPE), pointer :: q(:,:)
        real(kind=AUX_TYPE), pointer :: aux(:,:)
        
    end type solution_type

contains

    function new_solution(num_cells, num_eqn, num_aux, num_ghost) result(solution)

        implicit none

        integer, intent(in) :: num_cells, num_eqn, num_aux, num_ghost

        type(solution_type) :: solution

        ! Set new solution object's array extents
        solution%num_cells = num_cells
        solution%num_eqn = num_eqn
        solution%num_aux = num_aux
        solution%num_ghost = num_ghost

        ! Various other parameters - This needs to be set by the user
        solution%capa_index = 0

        ! Allocate memory for the solution arrays
        ! TODO: Probably should use the new solution's data, also should check
        !       for unsuccessful allocations
        allocate(solution%q(num_eqn,1-num_ghost:num_cells+num_ghost))
        allocate(solution%aux(num_aux,1-num_ghost:num_cells+num_ghost))
        
    end function new_solution

end module solution_module

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
        integer :: num_cells, num_eqn, num_aux
        real(kind=8) :: t, dx, lower, upper

        ! Aux array descriptors
        integer :: capa_index

        ! Solution arrays
        real(kind=Q_TYPE), pointer :: q(:,:)
        real(kind=AUX_TYPE), pointer :: aux(:,:)
        
    end type solution_type

contains

    type(solution_type) function new_solution(clawdata) result(solution)

        use clawdata_module, only: clawdata_type

        implicit none

        type(clawdata_type), intent(in) :: clawdata

        integer :: stat

        ! Set new solution object's array extents
        solution%num_cells = clawdata%num_cells
        solution%num_eqn = clawdata%num_eqn
        solution%num_aux = clawdata%num_aux
        solution%lower = clawdata%lower
        solution%upper = clawdata%upper
        solution%dx = (solution%upper - solution%lower) / solution%num_cells

        ! Set initial time
        solution%t = clawdata%t0

        ! Various other parameters - This needs to be set by the user
        solution%capa_index = clawdata%capa_index

        ! Allocate memory for the solution arrays
        associate(num_eqn => clawdata%num_eqn, &
                  num_ghost => clawdata%num_ghost, &
                  num_cells => clawdata%num_cells, &
                  num_aux => clawdata%num_aux)

            allocate(solution%q(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
            if (stat /= 0) stop "Allocation of solutions's q array failed!"
            allocate(solution%aux(num_aux,1-num_ghost:num_cells+num_ghost),stat=stat)
            if (stat /= 0) stop "Allocation of solutions's aux array failed!"

        end associate
        
    end function new_solution

end module solution_module

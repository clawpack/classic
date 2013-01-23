! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module solver_module

    use solution_module, only: Q_TYPE

    implicit none
    
    ! Solver type parameters - defaults to 1.d0
    integer, parameter :: F_TYPE = kind(1.d0)
    integer, parameter :: ASDQ_TYPE = kind(1.d0)
    integer, parameter :: WAVE_TYPE = kind(1.d0)
    integer, parameter :: S_TYPE = kind(1.d0)
    integer, parameter :: DTDX_TYPE = kind(1.d0)

    ! Solver type declaration
    type solver_type
        
        ! Status of solver
        integer :: num_steps
        real(kind=8) :: dt_min, dt_max, cfl, cfl_max, dt

        ! Solver parameters
        integer :: order, transverse_waves, dimensional_split, source_split
        integer :: num_waves, verbosity, steps_max
        integer, allocatable :: limiter(:)
        real(kind=8) :: dt_max_allowed, cfl_max_allowed, cfl_desired
        logical :: use_fwaves, dt_variable

        ! Boundary conditions
        integer :: num_ghost, bc_lower(1), bc_upper(1)

        ! Memory storage for this solver
        real(kind=F_TYPE), pointer :: f(:,:)
        real(kind=ASDQ_TYPE), pointer :: apdq(:,:), amdq(:,:)
        real(kind=WAVE_TYPE), pointer :: wave(:,:,:)
        real(kind=S_TYPE), pointer :: s(:,:)
        real(kind=DTDX_TYPE), pointer :: dtdx(:)

        real(kind=Q_TYPE), pointer :: q_old(:,:)

    end type solver_type

contains

    function new_solver(clawdata) result(solver)

        use clawdata_module, only: clawdata_type

        implicit none

        ! Input
        type(clawdata_type) :: clawdata

        ! Output
        type(solver_type) :: solver

        ! Locals
        integer :: stat

        ! Solver parameters
        solver%steps_max = clawdata%steps_max
        solver%dt_max_allowed = clawdata%dt_max_allowed
        solver%cfl_desired = clawdata%cfl_desired
        solver%cfl_max_allowed = solver%cfl_max_allowed

        ! Boundary conditions
        solver%num_ghost = clawdata%num_ghost
        solver%bc_lower = clawdata%bc_lower
        solver%bc_upper = clawdata%bc_upper

        ! Initialize solver status
        solver%num_steps = 0
        solver%dt_min = huge(1.d0)
        solver%dt_max = tiny(1.d0)
        solver%cfl = clawdata%cfl_desired
        solver%cfl_max = tiny(1.d0)
        solver%dt = clawdata%dt_init

        ! Allocate solver work arrays
        solver%num_waves = clawdata%num_waves

        associate(num_eqn => clawdata%num_eqn, &
                  num_cells => clawdata%num_cells(1), &
                  num_ghost => clawdata%num_ghost, &
                  num_waves => clawdata%num_waves)

        ! Based on number of waves and solution parameters allocate work arrays
        allocate(solver%f(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of f failed."
        allocate(solver%apdq(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of apdq failed."
        allocate(solver%amdq(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of amdq failed."
        allocate(solver%wave(num_eqn,num_waves,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of wave failed."
        allocate(solver%s(num_waves,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of s failed."
        allocate(solver%dtdx(1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of dtdx failed."
        
        end associate

    end function new_solver

    subroutine choose_new_dt(solver)

        implicit none

        type(solver_type), intent(in out) :: solver

        ! Change dt only if allowed
        if (solver%dt_variable) then
            ! Check if CFL makes sense and take proportionally the largest 
            ! desired time step
            if (solver%cfl > 0.d0) then
                solver%dt = min(solver%dt_max_allowed, solver%dt * solver%cfl_desired / solver%cfl)
                solver%dt_min = min(solver%dt_min, solver%dt)
                solver%dt_max = max(solver%dt_max, solver%dt)
            else
                ! CFL is negative, take the largest time step allowed
                solver%dt = solver%dt_max_allowed
            end if
        end if

    end subroutine choose_new_dt

    logical function check_CFL(solver)

        implicit none

        type(solver_type) :: solver

        if (solver%cfl > solver%cfl_max_allowed) then
            check_CFL = .false.
        else
            check_CFL = .true.
        endif

    end function check_CFL

end module solver_module

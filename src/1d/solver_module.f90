! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module solver_module

    use precision_module

    implicit none

    ! Riemann Solver abstract interface
!     abstract interface   

!         subroutine rp(num_eqn,num_aux,num_ghost,num_cells,num_waves,ql,qr, &
!                        auxl,auxr,wave,s,amdq,apdq)

!             use precision_module

!             implicit none

!             ! Input Arguments
!             integer, intent(in) :: num_eqn, num_aux, num_ghost, num_cells, num_waves
!             real(kind=Q_TYPE), intent(in) :: ql(num_eqn,1-num_ghost:num_cells+num_ghost)
!             real(kind=Q_TYPE), intent(in) :: qr(num_eqn,1-num_ghost:num_cells+num_ghost)
!             real(kind=AUX_TYPE), intent(in) :: auxl(num_aux,1-num_ghost:num_cells+num_ghost)
!             real(kind=AUX_TYPE), intent(in) :: auxr(num_aux,1-num_ghost:num_cells+num_ghost)

!             ! Output Arguments
!             real(kind=WAVE_TYPE), intent(in out) :: wave(num_eqn,num_waves,1-num_ghost:num_cells+num_ghost)
!             real(kind=S_TYPE), intent(in out) :: s(num_waves,1-num_ghost:num_cells+num_ghost)
!             real(kind=ASDQ_TYPE), intent(in out) :: amdq(num_eqn,1-num_ghost:num_cells+num_ghost)
!             real(kind=ASDQ_TYPE), intent(in out) :: apdq(num_eqn,1-num_ghost:num_cells+num_ghost)

!         end subroutine rp
!     end interface

    abstract interface
        subroutine rp(num_eqn, num_aux, num_waves, q_l, q_r,  &
                             aux_l, aux_r, wave, s, amdq, apdq)

            use precision_module

            implicit none

            ! Input Arguments
            integer, intent(in) :: num_eqn, num_aux, num_waves
            real(kind=Q_TYPE), intent(in) :: q_l(num_eqn), q_r(num_eqn)
            real(kind=AUX_TYPE), intent(in) :: aux_l(num_aux), aux_r(num_aux)

            ! Output arguments
            real(kind=WAVE_TYPE), intent(out) :: wave(num_eqn, num_waves)
            real(kind=S_TYPE), intent(out) :: s(num_waves)
            real(kind=ASDQ_TYPE), intent(out) :: apdq(num_eqn), amdq(num_eqn)

        end subroutine rp
    end interface

    ! Solver type declaration
    type solver_type
        
        ! Status of solver
        integer :: num_steps
        real(kind=8) :: dt_min, dt_max, cfl, cfl_max, dt

        ! Solver parameters
        integer :: order, transverse_waves, dimensional_split, source_splitting
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
        
        ! Function pointer to Riemann solver
        procedure (rp), nopass, pointer :: rp1

    end type solver_type

    interface new
        module procedure new_solver
    end interface

contains

    subroutine new_solver(self, clawdata)

        use clawdata_module, only: clawdata_type

        implicit none

        ! Input
        type(solver_type), intent(out) :: self
        type(clawdata_type), intent(in) :: clawdata

        ! Locals
        integer :: stat

        ! Solver parameters
        self%order = clawdata%order
        self%source_splitting = clawdata%source_splitting
        self%use_fwaves = clawdata%use_fwaves

        self%steps_max = clawdata%steps_max
        self%dt_max_allowed = clawdata%dt_max_allowed
        self%cfl_desired = clawdata%cfl_desired
        self%cfl_max_allowed = clawdata%cfl_max_allowed

        ! Boundary conditions
        self%num_ghost = clawdata%num_ghost
        self%bc_lower = clawdata%bc_lower
        self%bc_upper = clawdata%bc_upper

        ! Initialize solver status
        self%num_steps = 0
        self%dt_min = huge(1.d0)
        self%dt_max = tiny(1.d0)
        self%cfl = clawdata%cfl_desired
        self%cfl_max = tiny(1.d0)
        self%dt = clawdata%dt_init

        ! Allocate solver work arrays
        self%num_waves = clawdata%num_waves

        associate(num_eqn => clawdata%num_eqn, &
                  num_cells => clawdata%num_cells(1), &
                  num_ghost => clawdata%num_ghost, &
                  num_waves => clawdata%num_waves)

        ! Based on number of waves and solution parameters allocate work arrays
        allocate(self%f(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of f failed."
        allocate(self%apdq(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of apdq failed."
        allocate(self%amdq(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of amdq failed."
        allocate(self%wave(num_eqn,num_waves,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of wave failed."
        allocate(self%s(num_waves,1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of s failed."
        allocate(self%dtdx(1-num_ghost:num_cells+num_ghost),stat=stat)
        if (stat /= 0) stop "Allocation of dtdx failed."
        
        end associate

        ! Set the Riemann solver function pointer to null to avoid confusion
        self%rp1 => null()

    end subroutine new_solver

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

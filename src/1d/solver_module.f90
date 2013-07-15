! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module solver_module

    use iso_c_binding, only: c_ptr, C_NULL_PTR
    use precision_module, only: dp

    implicit none

    ! Riemann Solver abstract interface
    abstract interface   

        subroutine rp_vectorized(num_eqn, num_aux, num_ghost, num_cells, num_waves, rp_data, geometry,&
                                    ql, qr, auxl, auxr, wave, s, amdq, apdq)

            use iso_c_binding, only: c_ptr
            use precision_module, only: DP
            use geometry_module, only: geometry_type

            implicit none

            ! Input Arguments
            integer, intent(in) :: num_eqn, num_aux, num_ghost, num_cells, num_waves
            type(geometry_type), intent(in) :: geometry
            type(c_ptr), intent(in) :: rp_data
            real(kind=DP), intent(in) :: ql(num_eqn,1-num_ghost:num_cells+num_ghost)
            real(kind=DP), intent(in) :: qr(num_eqn,1-num_ghost:num_cells+num_ghost)
            real(kind=DP), intent(in) :: auxl(num_aux,1-num_ghost:num_cells+num_ghost)
            real(kind=DP), intent(in) :: auxr(num_aux,1-num_ghost:num_cells+num_ghost)

            ! Output Arguments
            real(kind=DP), intent(in out) :: wave(num_eqn,num_waves,1-num_ghost:num_cells+num_ghost)
            real(kind=DP), intent(in out) :: s(num_waves,1-num_ghost:num_cells+num_ghost)
            real(kind=DP), intent(in out) :: amdq(num_eqn,1-num_ghost:num_cells+num_ghost)
            real(kind=DP), intent(in out) :: apdq(num_eqn,1-num_ghost:num_cells+num_ghost)

        end subroutine rp_vectorized
    end interface

    abstract interface
        subroutine rp_ptwise(num_eqn, num_aux, num_waves, rp_data, geometry,   &
                        q_l, q_r, aux_l, aux_r, wave, s, amdq, apdq)

            use iso_c_binding, only: c_ptr
            use precision_module, only: DP
            use geometry_module, only: geometry_type

            implicit none

            ! Input Arguments
            integer, intent(in) :: num_eqn, num_aux, num_waves
            type(geometry_type), intent(in) :: geometry
            type(c_ptr), intent(in) :: rp_data
            real(kind=DP), intent(in) :: q_l(num_eqn), q_r(num_eqn)
            real(kind=DP), intent(in) :: aux_l(num_aux), aux_r(num_aux)

            ! Output arguments
            real(kind=DP), intent(in out) :: wave(num_eqn, num_waves)
            real(kind=DP), intent(in out) :: s(num_waves)
            real(kind=DP), intent(in out) :: apdq(num_eqn), amdq(num_eqn)
        end subroutine rp_ptwise
    end interface

    ! Solver type declaration
    type solver_type
        
        ! Status of solver
        integer :: num_steps
        real(dp) :: dt_min, dt_max, cfl, cfl_max, dt

        ! Solver parameters
        integer :: order, transverse_waves, dimensional_split, source_splitting
        integer :: num_waves, verbosity, steps_max
        integer, allocatable :: limiters(:)
        real(dp) :: dt_max_allowed, cfl_max_allowed, cfl_desired
        logical :: use_fwaves, dt_variable

        ! Time stepper function controls
!         procedure (CFL_check), nopass, pointer :: cfl_check_func
!         procedure (dt_chooser), nopass, pointer :: dt_chooser_func

        ! Boundary conditions
        integer :: num_ghost, bc_lower(1), bc_upper(1)

        ! Memory storage for this solver
        real(kind=DP), pointer :: f(:,:)
        real(kind=DP), pointer :: apdq(:,:), amdq(:,:)
        real(kind=DP), pointer :: wave(:,:,:)
        real(kind=DP), pointer :: s(:,:)
        real(kind=DP), pointer :: dtdx(:)

        real(kind=DP), pointer :: q_old(:,:)
        
        ! Function pointer to Riemann solver and Riemann solver data
        procedure (rp_ptwise), nopass, pointer :: rp1_ptwise
        procedure (rp_vectorized), nopass, pointer :: rp1_vectorized
        type(c_ptr) :: rp_data

    end type solver_type

    ! Time stepper control function interface definitions
!     abstract interface
!         logical function CFL_check(solver)
!             implicit none
!             type(solver_type), intent(in) :: solver
!         end function CFL_check
!     end interface

!     abstract interface
!         subroutine dt_chooser(solver)
!             implicit none
!             type(solver_type), intent(in out) :: solver
!         end subroutine dt_chooser
!     end interface

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

        ! Time stepping control
        self%steps_max = clawdata%steps_max
        self%dt_max_allowed = clawdata%dt_max_allowed
        self%cfl_desired = clawdata%cfl_desired
        self%cfl_max_allowed = clawdata%cfl_max_allowed

        ! Allocate solver work arrays
        self%num_waves = clawdata%num_waves

        allocate(self%limiters(self%num_waves))
        self%limiters = clawdata%limiters

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
        self%rp1_ptwise => null()
        self%rp1_vectorized => null()
        self%rp_data = C_NULL_PTR

    end subroutine new_solver


    logical function evolve_to_time(t_end, solution, solver, single_step) result(success)
        
        use solution_module, only: solution_type

        implicit none

        ! Arguments
        real(dp), intent(in) :: t_end
        type(solution_type), intent(in out) :: solution
        type(solver_type), intent(in out) :: solver
        logical, optional, intent(in) :: single_step

        ! Locals
        integer :: num_steps
        real(dp) :: t_old, t_new, t_start

        character(len=*), parameter :: stat_msg = "('CLAW1... Step',i4,"  // &
                             "'   Courant number =',f6.3,'  dt =',d12.4," // &
                             "'  t =',d12.4)" 

        ! Initialize return status
        success = .false.

        ! Initialize status storage
        num_steps = 0
        t_start = solution%t

        ! Support single stepping mode
        if (present(single_step)) then
            if (single_step) solver%steps_max = 1
        endif

        ! Primary Loop
        primary_loop: do num_steps=1,solver%steps_max

            t_old = solution%t
            t_new = solution%t + solver%dt

            ! Adjust dt to hit t_end exactly if we are near end of computation
            if (t_old + solver%dt > t_end .and. t_start < t_end) then
                solver%dt = t_end - t_old
            endif

            ! Save q in case we need to retake this step
            if (solver%dt_variable) then
                solver%q_old => solution%q
            endif

            ! Main steps in algorithm
            
            ! Extend data from actual domain to bordering boundary cells
            call set_boundary_conditions(solution,solver)

            ! Call user-supplied routine
            call before_step(solution,solver)

            ! If Strang splitting is used take a step on the source term
            if (solver%source_splitting == 2) then
                call source_term(t_old, solver%dt / 2.d0, solution, solver)
            endif

            ! Take a single-step on the homogeneous conservation law
            call hyperbolic_step(solution,solver)

            ! Take a step on the source term, if Strang splitting used take another
            ! half time step, if Godunov splitting take a full time step
            if (solver%source_splitting == 1) then
                call source_term(solution%t + solver%dt,solver%dt,solution,solver)
            else if (solver%source_splitting == 2) then
                call source_term(solution%t + solver%dt * 0.5d0,solver%dt * 0.5d0,solution,solver)
            end if

            ! Write out status message if verbosity /= 0
            if (solver%verbosity > 0) then
                print stat_msg, num_steps, solver%cfl, solver%dt, t_old + solver%dt
            end if

            ! Choose new time step if allowed
            call choose_new_dt(solver)

            ! Check CFL condition
            if (.not.check_CFL(solver)) then
                ! Reject this step and reset q
                print *, 'CLAW - Rejecting step... Courant number too large'
                solution%q = solver%q_old

                if (.not.solver%dt_variable) then
                    ! Not taking variable time steps so exit due to CFL violation
                    solver%cfl_max = max(solver%cfl,solver%cfl_max)
                    success = .false.
                    return
                endif

                ! Go back to beginning of loop, already should have choosen the 
                ! correct time step from above
                cycle primary_loop
            else
                ! Accept this step
                solver%cfl_max = max(solver%cfl,solver%cfl_max)
            end if

            ! Update solution time
            solution%t = t_new

            ! Check to see if we are done
            if (solution%t >= t_end) then
                exit primary_loop
            endif

        end do primary_loop

        success = .true.
        
    end function evolve_to_time


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

        type(solver_type), intent(in) :: solver

        if (solver%cfl > solver%cfl_max_allowed) then
            check_CFL = .false.
        else
            check_CFL = .true.
        endif

    end function check_CFL

end module solver_module

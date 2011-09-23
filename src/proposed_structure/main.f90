subroutine solution_status(solution)

    use solution_1d_module
    use state_1d_module

    implicit none
    
    type(solution_type), intent(inout) :: solution
    type(state_type) :: state
    type(grid_type) :: grid
    
    print *,"Fetching state"
    state = get_state(solution,1)
    print *,"Fetching grid"
    grid = state%grid
    
    print *,"Start output"
    print *,"  q(1,:)="
    print *, state%q(1,:)
    print *,"  q(2,:)="
    print *, state%q(2,:)
    print *,"  aux="
    print *, state%aux
    
    print *,grid%n
    
end subroutine solution_status

program clawpack_structure_test
    
    ! Clawpack utility functions
    use utility_module
    
    ! State date structures
    use grid_module
    use state_1d_module
    use solution_1d_module
    
    ! Solver date structures
    use solver_types_module
    
    ! Riemann solver
    use rp_euler_module
    
    implicit none
    
    ! Data types
    type(solution_type), pointer :: solution
    type(grid_type), pointer :: grid
    type(state_type), pointer :: state

    ! Solver data
    type(solver_parameters_type), pointer :: solver_params
    type(solver_status_type), pointer :: status
    
    ! Riemann solver static data
    type(rp_params_type), pointer :: rp_params
    type(rp_work_type), pointer :: rp_work
    
    double precision, pointer :: q(:,:),aux(:,:)
    
    integer :: ios
    
    ! TEMPORARY
    double precision :: gamma = 1.4d0
    
    ! ========================================================================
    !  Solution setup
    print *, "Starting"
    grid => new_grid([10])
    state => new_state(2,1,grid)
    print *, "Created base objects"

    print *, "Setting conserved quantities"
    q => get_q(state)
    q(1,:) = 1.d0
    q(2,:) = 2.d0
    
    print *,"Setting aux array"
    aux => get_aux(state)
    aux(1,:) = gamma
    
    print *,"Creating solution"
    allocate(solution)
    
    print *,"Adding grid to solution"
    call add_grid(solution,grid)
    print *,"Adding state to solution"
    call add_state(solution,state)
    
    call solution_status(solution)
    
    ! ========================================================================
    ! Solver setup
    solver_params => new_solver_params("./claw.data",1)
    
    rp_params => new_rp_params("./rp.data")
    rp_work => new_rp_work(solution)
    
    print *,rp_params%gamma
    
end program clawpack_structure_test

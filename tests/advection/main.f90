program advection1d

    use iso_c_binding, only: c_ptr, c_loc, c_associated

    use solver_module
    use solution_module
    use controller_module
    use rp1_advection, only: rp1_advection_ptwise, rp1_advection_vectorized
    use rp1_advection, only: rp_type
    use precision_module, only: dp

    use clawdata_module
    use utils, only: open_data_file

    implicit none
    
    type(solution_type) :: solution
    type(solver_type) :: solver
    type(clawdata_type) :: clawdata

    type(rp_type), target :: rp_aux

    real(dp) :: beta
    integer :: i

    ! Read in advective speed and initial condition data
    call open_data_file(13,'./setprob.data')
    read(13,*) rp_aux%u
    read(13,*) beta
    close(13)

    ! Read input data
    call new(clawdata,'./claw.data')

    ! Initialize solution
    call new(solution,clawdata)
    solution%q(1,:) = exp(-beta * (solution%centers - 0.2d0)**2)
    forall(i=1:solution%num_cells(1), solution%centers(i) > 0.7d0 .and.   &
                                      solution%centers(i) < 0.8d0)
        solution%q(1,i) = 1.d0
    end forall

    ! Initialize solver
    call new(solver,clawdata)
    solver%rp1_ptwise => rp1_advection_ptwise
!     solver%rp1_vectorized => rp1_advection_vectorized
    solver%rp_data = c_loc(rp_aux)

    ! Start output loop
    if (.not.run_simulation(solution, solver, clawdata)) then
        print *,"uh oh"
    endif

end program advection1d

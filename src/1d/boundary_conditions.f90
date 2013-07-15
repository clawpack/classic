subroutine set_boundary_conditions(solution,solver)

    use solution_module, only: solution_type
    use solver_module, only: solver_type
    use utils, only: stop_error

    implicit none

    type(solution_type), intent(in out) :: solution
    type(solver_type), intent(in out) :: solver

    integer :: i, m

    associate(num_cells => solution%num_cells(1), &
              num_ghost => solver%num_ghost, &
              q => solution%q)

        ! Lower boundary
        select case (solver%bc_lower(1))
            case (0)
                call stop_error("User defined boundary condintion requested but not implemented.")
            case (1)
                ! Zero-order extrapolation
                forall (i=1:num_ghost)
                    q(:,1-i) = q(:,1)
                end forall
            case (2)
                ! Periodic
                forall (i=1:num_ghost)
                    q(:,1-i) = q(:,num_cells+1-i)
                end forall
            case (3)
                ! Solid wall, (assumes 2nd component of q is velocity or momentum)
                forall (i=1:num_ghost, m=1:solution%num_eqn, m /= 2)
                    q(m,1-i) = q(m,i)
                end forall
                forall (i=1:num_ghost)
                    q(2,1-i) = -q(2,i)
                end forall
        end select

        ! Upper boundary
        select case (solver%bc_upper(1))
            case (0)
                call stop_error("User defined boundary condintion requested but not implemented.")
            case (1)
                ! Zero-order extrapolation
                forall (i=1:num_ghost)
                    q(:,num_cells + i) = q(:,num_cells)
                end forall
            case (2)
                ! Periodic
                forall (i=1:num_ghost)
                    q(:,num_cells+i) = q(:,i)
                end forall
            case (3)
                ! Solid wall, (assumes 2nd component of q is velocity or momentum)
                forall (i=1:num_ghost, m=1:solution%num_eqn, m /= 2)
                    q(m,num_cells+i) = q(m,num_cells+1-i)
                end forall
                forall (i=1:num_ghost)
                    q(2,num_cells+i) = -q(2,num_cells+1-i)
                end forall
        end select

    end associate


end subroutine set_boundary_conditions

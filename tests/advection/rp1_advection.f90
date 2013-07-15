module rp1_advection

    use precision_module, only: dp
    implicit none
    
    type rp_type
        ! Advection speed
        real(dp) :: u
    end type rp_type

contains

    ! Simple sweeping 1d advection Riemann solver
    subroutine rp1_advection_vectorized(num_eqn, num_aux, num_ghost, num_cells, num_waves, rp_data, geometry,&
                   ql, qr, auxl, auxr, wave, s, amdq, apdq)

        use iso_c_binding, only: c_ptr, c_loc, c_f_pointer
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

        ! Locals
        integer :: i
        type(rp_type), pointer :: rp_aux

        ! Caste riemann solver data
        call c_f_pointer(rp_data, rp_aux)

        ! Initialize cumulative amdq and apdq
        amdq = 0.0_dp
        apdq = 0.0_dp

        ! Speeds are all the same
        s = rp_aux%u

        ! Wave is equal to the jump at the grid cell interfaces
        forall(i=2-num_ghost:num_cells+num_ghost)
            wave(1,1,i) = ql(1,i) - qr(1,i-1)
        end forall

        ! Set fluctuations
        if (rp_aux%u < 0.0_dp) then
            amdq(1,:) = rp_aux%u * wave(1,1,:)
        else
            apdq(1,:) = rp_aux%u * wave(1,1,:)
        endif

    end subroutine rp1_advection_vectorized

    ! Point-wise constant advection Riemann solver
    subroutine rp1_advection_ptwise(num_eqn, num_aux, num_waves, rp_data, geometry,   &
                         q_l, q_r, aux_l, aux_r, wave, s, amdq, apdq)

        use iso_c_binding, only: c_ptr, c_f_pointer
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

        ! Local data
        type(rp_type), pointer :: rp_aux
        
        call c_f_pointer(rp_data, rp_aux)

        ! Initialize cumulative amdq and apdq
        amdq = 0.0_dp
        apdq = 0.0_dp

        ! Speeds are all the same
        s = rp_aux%u

        ! Wave is equal to the jump at the grid cell interfaces
        wave(1,1) = q_r(1) - q_l(1)

        ! Set fluctuations
        amdq = min(rp_aux%u, 0.0_dp) * wave(1,1)
        apdq = max(rp_aux%u, 0.0_dp) * wave(1,1)

    end subroutine rp1_advection_ptwise

end module rp1_advection

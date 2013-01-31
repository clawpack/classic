module rp1_advection

    implicit none
    
    ! Special data
    real(kind=8) :: u

contains

    ! Simple sweeping 1d advection Riemann solver
    subroutine rp1(num_eqn,num_aux,num_ghost,num_cells,num_waves,ql,qr,auxl,auxr, &
                   wave,s,amdq,apdq)

        use precision_module, only: Q_TYPE, AUX_TYPE, WAVE_TYPE, S_TYPE, ASDQ_TYPE 

        implicit none

        ! Input Arguments
        integer, intent(in) :: num_eqn, num_aux, num_ghost, num_cells, num_waves
        real(kind=Q_TYPE), intent(in) :: ql(num_eqn,1-num_ghost:num_cells+num_ghost)
        real(kind=Q_TYPE), intent(in) :: qr(num_eqn,1-num_ghost:num_cells+num_ghost)
        real(kind=AUX_TYPE), intent(in) :: auxl(num_aux,1-num_ghost:num_cells+num_ghost)
        real(kind=AUX_TYPE), intent(in) :: auxr(num_aux,1-num_ghost:num_cells+num_ghost)

        ! Output Arguments
        real(kind=WAVE_TYPE), intent(in out) :: wave(num_eqn,num_waves,1-num_ghost:num_cells+num_ghost)
        real(kind=S_TYPE), intent(in out) :: s(num_waves,1-num_ghost:num_cells+num_ghost)
        real(kind=ASDQ_TYPE), intent(in out) :: amdq(num_eqn,1-num_ghost:num_cells+num_ghost)
        real(kind=ASDQ_TYPE), intent(in out) :: apdq(num_eqn,1-num_ghost:num_cells+num_ghost)

        ! Locals
        integer :: i

        ! Initialize cumulative amdq and apdq
        amdq = 0.d0
        apdq = 0.d0

        ! Speeds are all the same
        s = u

        ! Wave is equal to the jump at the grid cell interfaces
        forall(i=2-num_ghost:num_cells+num_ghost)
            wave(1,1,i) = ql(1,i) - qr(1,i-1)
        end forall

        ! Set fluctuations
        if (u < 0.d0) then
            amdq(1,:) = u * wave(1,1,:)
        else
            apdq(1,:) = u * wave(1,1,:)
        endif

    end subroutine rp1

    ! Point-wise constant advection Riemann solver
    subroutine rp_ptwise(num_eqn,num_aux,num_waves,q_l,q_r,aux_l,aux_r,wave,s,amdq,apdq)

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

        ! Initialize cumulative amdq and apdq
        amdq = 0.d0
        apdq = 0.d0

        ! Speeds are all the same
        s = u

        ! Wave is equal to the jump at the grid cell interfaces
        wave(1,1) = q_r(1) - q_l(1)

        ! Set fluctuations
        amdq = min(u, 0.d0) * wave(1,1)
        apdq = max(u, 0.d0) * wave(1,1)

    end subroutine rp_ptwise

end module rp1_advection
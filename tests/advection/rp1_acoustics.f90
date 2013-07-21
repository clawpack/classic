module rp1_acoustics

    use precision_module, only: dp
    implicit none
    
    type rp_type
        real(dp) :: rho, bulk, cc, zz
    end type rp_type

contains


!     =====================================================
    !subroutine rp1(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr, &
    !wave,s,amdq,apdq,num_aux)
    subroutine rp1(num_eqn, num_aux, num_ghost, num_cells, num_waves, rp_data, geometry,&
                   ql, qr, auxl, auxr, wave, s, amdq, apdq)
!     =====================================================

!     # Riemann solver for the acoustics equations in 1d,

!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell

!     # On output, wave contains the waves,
!     #            s the speeds,
!     #
!     #            amdq = A^- Delta q,
!     #            apdq = A^+ Delta q,
!     #                   the decomposition of the flux difference
!     #                       f(qr(i-1)) - f(ql(i))
!     #                   into leftgoing and rightgoing parts respectively.
!     #

!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routines, this routine is called with ql = qr
        use iso_c_binding, only: c_ptr, c_loc, c_f_pointer
        use geometry_module, only: geometry_type

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


!     local arrays
!     ------------
    dimension delta(2)
    real(dp) :: a1, a2, delta

!     # density, bulk modulus, and sound speed, and impedence of medium:
!     # (should be set in setprob.f)
    integer :: i, m

    type(rp_type), pointer :: rp_aux

    ! Caste riemann solver data
    call c_f_pointer(rp_data, rp_aux)


!     # split the jump in q at each interface into waves

!     # find a1 and a2, the coefficients of the 2 eigenvectors:
    do 20 i = 2-num_ghost, num_cells+num_ghost
        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = ql(2,i) - qr(2,i-1)
        a1 = (-delta(1) + rp_aux%zz*delta(2)) / (2.d0*rp_aux%zz)
        a2 =  (delta(1) + rp_aux%zz*delta(2)) / (2.d0*rp_aux%zz)
    
    !        # Compute the waves.
    
        wave(1,1,i) = -a1*rp_aux%zz
        wave(2,1,i) = a1
        s(1,i) = -rp_aux%cc
    
        wave(1,2,i) = a2*rp_aux%zz
        wave(2,2,i) = a2
        s(2,i) = rp_aux%cc
    
    20 END DO


!     # compute the leftgoing and rightgoing flux differences:
!     # Note s(1,i) < 0   and   s(2,i) > 0.

    do 220 m=1,num_eqn
        do 220 i = 2-num_ghost, num_cells+num_ghost
            amdq(m,i) = s(1,i)*wave(m,1,i)
            apdq(m,i) = s(2,i)*wave(m,2,i)
    220 END DO

    return
    end subroutine rp1

end

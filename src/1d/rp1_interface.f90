! ==============================================================================
!  subroutine rp1(num_eqn,num_aux,num_ghost,num_cells,num_waves,ql,qr,auxl,auxr,
!                 wave,s,amdq,apdq)
!  Solve Riemann problems for the 1D hyperbolic problem
!    This is a dummy routine and is only intended to illustrate the function
!    signature of this routine.  See various examples from the Riemann
!    repository.
!
!  On input, ql contains the state vector at the left edge of each cell
!            qr contains the state vector at the right edge of each cell
!  On output, wave contains the waves
!             s the speeds
!             amdq the left-going flux differences A^- \Delta q
!             apdq the right-going flux differences A^+ \Delta q
!
!   Note that the ith Riemann problem is defined by qr(:,i-1) as the left state
!   and ql(:,i) as the right state with resulting waves, s, amdq, and apdq with
!   indices wave(:,:,i), s(:,i), amdq(:,i), and apdq(:,i)
!
!               C_{i-1}             C_{i}                C_{i+1}  
!   +----------------------+----------------------+----------------------+
!   | ql(:,i-1)   qr(:,i-1)|ql(:,i)        qr(:,i)|ql(:,i+1)    qr(:,i+1)|
!   |                      |                      |                      |
!   |                 wave(:,:,i)            wave(:,:,i+1)               |
!   |                    s(:,i)                 s(:,i+1)                 |
!   |                amdq,apdq(:,i)         amdq,apdq(:,i+1)             |
! ==============================================================================
module rp1_interface

contains

    subroutine rp1(num_eqn,num_aux,num_ghost,num_cells,num_waves,ql,qr,auxl,auxr, &
                   wave,s,amdq,apdq)

        use solver_module, only: WAVE_TYPE, S_TYPE, ASDQ_TYPE
        use solution_module, only: Q_TYPE, AUX_TYPE

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

        stop "You must provide a valid Riemann solver!"

    end subroutine rp1

end module rp1_interface
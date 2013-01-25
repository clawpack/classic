module precision_module

    implicit none

    ! All of these default to double precision

    ! Solution type parameters
    integer, parameter :: Q_TYPE = kind(1.d0)
    integer, parameter :: AUX_TYPE = kind(1.d0)

    ! Solver type parameters
    integer, parameter :: F_TYPE = kind(1.d0)
    integer, parameter :: ASDQ_TYPE = kind(1.d0)
    integer, parameter :: WAVE_TYPE = kind(1.d0)
    integer, parameter :: S_TYPE = kind(1.d0)
    integer, parameter :: DTDX_TYPE = kind(1.d0)

end module precision_module
module precision_module
    THIS WILL FAIL

    implicit none

    ! All of these default to double precision
    integer, parameter :: DP = kind(1.d0)

end module precision_module

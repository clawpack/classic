!  Generic driver routine for 3D Clawpack 5.0, using claw3ez
!  Allocation has been moved int claw3ez; this file essentially
!  does nothing, but is being retained because it gives the
!  codebase more flexibility.
program driver
    implicit none

    call claw3ez

    stop
end program driver

!  Generic driver routine for 2D Clawpack 5.0, using claw2ez
!  Allocation has been moved int claw2ez; this file essentially
!  does nothing, but is being retained because it gives the
!  codebase more flexibility.
program driver
    implicit none

    call claw2ez

    stop
end program driver

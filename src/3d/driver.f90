!  Generic driver routine for claw3 poroelasticity code
!  Uses dynamic allocation
!
!  Author: Grady Lemoine, based on version of Donna Calhoun
!  Date : 3/30/2013, based on Calhoun code of 3/22/2002
program driver
    implicit none

    ! These will always be the same for this poroelasticity code.
    ! narray should be 2 for Strang splitting.  Not specifying meqn
    ! here to leave things flexible for high-frequency extension at
    ! some point.
    integer, parameter :: mbc = 2, mwaves = 8, maux = 20, narray = 2

    integer :: mx, my, mz, meqn, mwork, maxm
    double precision, dimension(:,:,:,:), allocatable :: q, aux
    double precision, dimension(:), allocatable :: work
    double precision, dimension(mwaves) :: mthlim

    ! Working variables needed for reading claw.data
    integer :: tmp, outstyle, i, nout
    double precision :: dtmp
    double precision, dimension(:), allocatable :: tout

    ! Get problem dimensions, so we can size workspace
    open(55,file='claw3ez.data',status='old',form='formatted')

    ! Need to read a fair way into the file to get meqn
    read (55,*) mx
    read (55,*) my
    read (55,*) mz
    read (55,*) nout
    read (55,*) outstyle
    if (outstyle == 1) then
        read (55,*) dtmp
    else if (outstyle == 2) then
        allocate (tout(nout))
        read (55,*) (tout(i), i = 1,nout)
        deallocate (tout)
    else if (outstyle == 3) then
        read (55,*) tmp, tmp
    end if
    read (55,*) dtmp    ! dtv(1)
    read (55,*) dtmp    ! dtv(2)
    read (55,*) dtmp    ! cflv(1)
    read (55,*) dtmp    ! cflv(2)
    read (55,*) tmp     ! nv(1)
    do i = 1,7
        read (55,*) tmp    ! method(1:7)
    end do
    read (55,*) meqn
    ! Have meqn, done with file
    close(55)

    allocate(  q(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, meqn))
    allocate(aux(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, maux))

    ! Copied from work estimate in claw3ez
    maxm = max0(mx, my, mz)
    mwork = (maxm+2*mbc)*(46*meqn + mwaves + meqn*mwaves + 9*maux + 3) &
          + narray*(mx + 2*mbc)*(my + 2*mbc)*(mz + 2*mbc)*meqn
    allocate(work(mwork))

    call claw3ez(mx,my,mz,meqn,mwaves,mbc,maux,mwork,mthlim,q,work,aux)

    deallocate(q)
    deallocate(aux)

    stop
end program driver

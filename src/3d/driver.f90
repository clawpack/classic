!  Generic driver routine for 3D Clawpack 5.0, using claw3ez
!  Uses dynamic allocation
!
!  Author: Grady Lemoine, based on version of Donna Calhoun
!  Date: 7/24/2013, based on Lemoine poroelasticity code of 3/30/2013,
!                   based on Calhoun code of 3/22/2002
program driver
    implicit none

    integer :: mx, my, mz, meqn, mwork, maxm, maux, narray, mwaves, mbc
    double precision, dimension(:,:,:,:), allocatable :: q, aux
    double precision, dimension(:), allocatable :: work

    ! Working variables needed for reading claw.data
    integer :: tmp, outstyle, i, nout, splitstyle
    double precision :: dtmp

    ! Get problem dimensions, so we can size workspace
    call opendatafile(55, 'claw.data')
    !open(55,file='claw.data',status='old',form='formatted')

    ! Need to read a fair way into the file to get meqn
    read (55,*) tmp
    read (55,*) mx
    read (55,*) my
    read (55,*) mz
    read (55,*) nout
    read (55,*) outstyle
    if (outstyle == 1) then
        read (55,*) dtmp
    else if (outstyle == 2) then
        read (55,*) (dtmp, i = 1,nout)
    else if (outstyle == 3) then
        read (55,*) tmp, tmp
    end if
    read (55,*) dtmp    ! dtv(1)
    read (55,*) dtmp    ! dtv(2)
    read (55,*) dtmp    ! cflv(1)
    read (55,*) dtmp    ! cflv(2)
    read (55,*) tmp     ! nv(1)
    do i = 1,4
        read (55,*) tmp    ! method(1:4)
    end do
    read (55,*) splitstyle ! method(5) == splitting method
    if (splitstyle < 2) then
        narray = 1
    else
        narray = 2
    end if
    read (55,*) tmp        ! method(6)
    read (55,*) maux       ! method(7) == maux
    read (55,*) meqn
    read (55,*) mwaves
    read (55,*) (tmp, i = 1,mwaves)    ! Limiter type for each wave
    do i = 1,7
        read (55,*) dtmp    ! t0, lower and upper bounds on computational domain
    end do
    read (55,*) mbc    ! Number of ghost cells
    ! Have all data needed to allocate arrays, done with file
    close(55)

    allocate(  q(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, meqn))
    allocate(aux(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, maux))

    ! Copied from work estimate in claw3ez
    maxm = max0(mx, my, mz)
    mwork = (maxm+2*mbc)*(46*meqn + mwaves + meqn*mwaves + 9*maux + 3) &
          + narray*(mx + 2*mbc)*(my + 2*mbc)*(mz + 2*mbc)*meqn
    allocate(work(mwork))

    call claw3ez(mx,my,mz,meqn,mbc,maux,mwork,q,work,aux)

    deallocate(q)
    deallocate(aux)
    deallocate(work)

    stop
end program driver

!     ==========================================================
    subroutine step2(maxm,num_eqn,num_waves,num_aux,num_ghost,mx,my, &
    qold,qnew,aux,xlower,ylower,dx,dy,dt,method,mthlim,cfl, &
    qadd,fadd,gadd,q1d,dtdx1d,dtdy1d, &
    aux1,aux2,aux3,work,mwork,use_fwave,rpn2,rpt2)
!     ==========================================================

!     # Take one time step, updating q.
!     # On entry, qold and qnew should be identical and give the
!     #    initial data for this step
!     # On exit, qnew returns values at the end of the time step.
!     #    qold is unchanged.

!     # qadd is used to return increments to q from flux2
!     # fadd and gadd are used to return flux increments from flux2.
!     # See the flux2 documentation for more information.

    use abl_module, only: abl_type, set_scaling_factor

    implicit double precision (a-h,o-z)
    external :: rpn2,rpt2
    double precision :: qold(num_eqn,1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: qnew(num_eqn,1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision ::  q1d(num_eqn,1-num_ghost:maxm+num_ghost)
    double precision :: qadd(num_eqn,1-num_ghost:maxm+num_ghost)
    double precision :: fadd(num_eqn,1-num_ghost:maxm+num_ghost)
    double precision :: gadd(num_eqn, 2, 1-num_ghost:maxm+num_ghost)
    double precision :: aux(num_aux, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: aux1(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux2(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux3(num_aux, 1-num_ghost:maxm+num_ghost)

    double precision :: dtdx1d(1-num_ghost:maxm+num_ghost)
    double precision :: dtdy1d(1-num_ghost:maxm+num_ghost)
    integer :: method(7),mthlim(num_waves)
    logical ::          use_fwave
    double precision :: work(mwork)

    ! added for abl implementation
    integer :: i_abl_lower, i_abl_upper, j_abl_lower, j_abl_upper
    double precision :: xlower, ylower, dtdxm, dtdxc, dtdxp, dtdyp, dtdyc, dtdym
    double precision :: ablX_M(1-num_ghost:mx+num_ghost)
    double precision :: ablX_J(1-num_ghost:mx+num_ghost)
    double precision :: ablX_R(1-num_ghost:mx+num_ghost)
    double precision :: ablY_M(1-num_ghost:my+num_ghost)
    double precision :: ablY_J(1-num_ghost:my+num_ghost)
    double precision :: ablY_R(1-num_ghost:my+num_ghost)

    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

!f2py intent(out) cfl
!f2py intent(in,out) qnew
!f2py optional q1d, qadd, fadd, gadd, dtdx1d, dtdy1d

! Dummy interfaces just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn2(x)
!f2py x=rpt2(x)



!     # partition work array into pieces needed for local storage in
!     # flux2 routine.  Find starting index of each piece:

    i0wave = 1
    i0s = i0wave + (maxm+2*num_ghost)*num_eqn*num_waves
    i0amdq = i0s + (maxm+2*num_ghost)*num_waves
    i0apdq = i0amdq + (maxm+2*num_ghost)*num_eqn
    i0cqxx = i0apdq + (maxm+2*num_ghost)*num_eqn
    i0bmadq = i0cqxx + (maxm+2*num_ghost)*num_eqn
    i0bpadq = i0bmadq + (maxm+2*num_ghost)*num_eqn
    iused = i0bpadq + (maxm+2*num_ghost)*num_eqn - 1

    if (iused > mwork) then
    !        # This shouldn't happen due to checks in claw2
        write(6,*) '*** not enough work space in step2'
        write(6,*) '*** iused = ', iused, '   mwork =',mwork
        stop
    endif


    index_capa = method(6)
    num_aux = method(7)
    cfl = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
    dtdxm = dtdx
    dtdxc = dtdx
    dtdxp = dtdx
    dtdym = dtdy
    dtdyc = dtdy
    dtdyp = dtdy



    if (index_capa == 0) then
    !        # no capa array:
        do i=1-num_ghost,maxm+num_ghost
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
        end do
    endif

!   compute abl scaling factors if needed
    if (abl_type /= 0) then
        call set_scaling_factor(xlower,ylower,dx,dy,mx,my,num_ghost, &
                                ablX_M,ablX_J,ablX_R,ablY_M,ablY_J,ablY_R, &
                                i_abl_lower,i_abl_upper,j_abl_lower,j_abl_upper)
    else
        i_abl_lower = 0
        j_abl_lower = 0
        i_abl_upper = mx+1
        j_abl_upper = my+1
    end if

!     # perform x-sweeps
!     ==================

    do j = 0,my+1

    !        # copy data along a slice into 1d arrays:
        do i = 1-num_ghost, mx+num_ghost
            do m=1,num_eqn
                q1d(m,i) = qold(m,i,j)
            end do
        end do

        if (index_capa > 0)  then
            do i = 1-num_ghost, mx+num_ghost
                dtdx1d(i) = dtdx/aux(index_capa,i,j)
            end do
        endif

        if (num_aux > 0)  then
            do i = 1-num_ghost, mx+num_ghost
                do ma=1,num_aux
                    aux1(ma,i) = aux(ma,i,j-1)
                    aux2(ma,i) = aux(ma,i,j  )
                    aux3(ma,i) = aux(ma,i,j+1)
                end do
            end do
        endif

    !     # Store the value of j along this slice in the common block
    !        # comxyt in case it is needed in the Riemann solver (for
    !        # variable coefficient problems)
        jcom = j

        ! adjust dtdy* if ABL is present
        if (abl_type /= 0) then
            if (j-1 <= j_abl_lower .or. j-1 >= j_abl_upper) then
                dtdym = dtdy * ablY_R(j-1) * ablY_M(j-1)
            else
                dtdym = dtdy
            end if
            if (j <= j_abl_lower .or. j >= j_abl_upper) then
                dtdyc = dtdy * ablY_R(j) * ablY_M(j)
                dtdxc = dtdx * ablY_R(j)
            else
                dtdyc = dtdy
                dtdxc = dtdx
            end if
            if (j+1 <= j_abl_lower .or. j+1 >= j_abl_upper) then
                dtdyp = dtdy * ablY_R(j+1) * ablY_M(j+1)
            else
                dtdyp = dtdy
            end if
        end if

    !        # compute modifications fadd and gadd to fluxes along this slice:
        call flux2(1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx, &
        q1d,dtdx1d,aux1,aux2,aux3,method,mthlim, &
        qadd,fadd,gadd,cfl1d, &
        work(i0wave),work(i0s),work(i0amdq),work(i0apdq), &
        work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
        use_fwave,ablX_J,i_abl_lower,i_abl_upper)
        cfl = dmax1(cfl,cfl1d)

    !        # update qnew by flux differencing.
    !        # (rather than maintaining arrays f and g for the total fluxes,
    !        # the modifications are used immediately to update qnew
    !        # in order to save storage.)

        if (index_capa == 0) then

        !            # no capa array.  Standard flux differencing:
            do i=i_abl_lower+1,i_abl_upper-1
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + dtdxc * qadd(m,i) &
                    - dtdxc * (fadd(m,i+1) - fadd(m,i)) &
                    - dtdyc * (gadd(m,2,i) - gadd(m,1,i))
                    qnew(m,i,j-1) = qnew(m,i,j-1) - dtdym * gadd(m,1,i)
                    qnew(m,i,j+1) = qnew(m,i,j+1) + dtdyp * gadd(m,2,i)
                end do
            end do
            ! these loops are skipped if no ABL present
            do i=1,i_abl_lower
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + ablX_M(i)*( dtdxc * qadd(m,i) &
                    - dtdxc * (fadd(m,i+1) - fadd(m,i)) &
                    - dtdyc * (gadd(m,2,i) - gadd(m,1,i)) )
                    qnew(m,i,j-1) = qnew(m,i,j-1) - dtdym*ablX_M(i) * gadd(m,1,i)
                    qnew(m,i,j+1) = qnew(m,i,j+1) + dtdyp*ablX_M(i) * gadd(m,2,i)
                end do
            end do
            do i=i_abl_upper,mx
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + ablX_M(i)*( dtdxc * qadd(m,i) &
                    - dtdxc * (fadd(m,i+1) - fadd(m,i)) &
                    - dtdyc * (gadd(m,2,i) - gadd(m,1,i)) )
                    qnew(m,i,j-1) = qnew(m,i,j-1) - dtdym*ablX_M(i) * gadd(m,1,i)
                    qnew(m,i,j+1) = qnew(m,i,j+1) + dtdyp*ablX_M(i) * gadd(m,2,i)
                end do
            end do

        else
        !            # with capa array.
            if (abl_type /= 0) then
                print *, "ABL not implemented yet for capa array"
                stop
            end if
            do i=1,mx
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + dtdx * qadd(m,i) &
                    - (dtdx * (fadd(m,i+1) - fadd(m,i)) &
                    +  dtdy * (gadd(m,2,i) - gadd(m,1,i))) &
                    / aux(index_capa,i,j)
                    qnew(m,i,j-1) = qnew(m,i,j-1) - dtdy * gadd(m,1,i) &
                    / aux(index_capa,i,j-1)
                    qnew(m,i,j+1) = qnew(m,i,j+1) + dtdy * gadd(m,2,i) &
                    / aux(index_capa,i,j+1)
                end do
            end do
        endif
    end do    ! End x sweeps



!     # perform y sweeps
!     ==================

    do i = 0, mx+1

    !        # copy data along a slice into 1d arrays:
        do j = 1-num_ghost, my+num_ghost
            do m=1,num_eqn
                q1d(m,j) = qold(m,i,j)
            end do
        end do

        if (index_capa > 0)  then
            do j = 1-num_ghost, my+num_ghost
                dtdy1d(j) = dtdy / aux(index_capa,i,j)
            end do
        endif

        if (num_aux > 0)  then
            do j = 1-num_ghost, my+num_ghost
                do ma=1,num_aux
                    aux1(ma,j) = aux(ma,i-1,j)
                    aux2(ma,j) = aux(ma,i,  j)
                    aux3(ma,j) = aux(ma,i+1,j)
                end do
            end do
        endif

    !     # Store the value of i along this slice in the common block
    !        # comxyt in case it is needed in the Riemann solver (for
    !        # variable coefficient problems)
        icom = i

        ! adjust dtdx* if ABL is present
        if (abl_type /= 0) then
            if (i-1 <= i_abl_lower .or. i-1 >= i_abl_upper) then
                dtdxm = dtdx * ablX_R(i-1) * ablX_M(i-1)
            else
                dtdxm = dtdx
            end if
            if (i <= i_abl_lower .or. i >= i_abl_upper) then
                dtdxc = dtdx * ablX_R(i) * ablX_M(i)
                dtdyc = dtdy * ablX_R(i)
            else
                dtdxc = dtdx
                dtdyc = dtdy
            end if
            if (i+1 <= i_abl_lower .or. i+1 >= i_abl_upper) then
                dtdxp = dtdx * ablX_R(i+1) * ablX_M(i+1)
            else
                dtdxp = dtdx
            end if
        end if


    !        # compute modifications fadd and gadd to fluxes along this slice:
        call flux2(2,maxm,num_eqn,num_waves,num_aux,num_ghost,my, &
        q1d,dtdy1d,aux1,aux2,aux3,method,mthlim, &
        qadd,fadd,gadd,cfl1d, &
        work(i0wave),work(i0s),work(i0amdq),work(i0apdq), &
        work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
        use_fwave,ablY_J,j_abl_lower,j_abl_upper)

        cfl = dmax1(cfl,cfl1d)

    !        # update qnew by flux differencing.
    !        # Note that the roles of fadd and gadd are reversed for
    !        # the y-sweeps -- fadd is the modification to g-fluxes and
    !        # gadd is the modification to f-fluxes to the left and right.

        if (index_capa == 0) then

        !            # no capa array.  Standard flux differencing:
            do j=j_abl_lower+1,j_abl_upper-1
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + dtdyc * qadd(m,j) &
                    - dtdyc * (fadd(m,j+1) - fadd(m,j)) &
                    - dtdxc * (gadd(m,2,j) - gadd(m,1,j))
                    qnew(m,i-1,j) = qnew(m,i-1,j) - dtdxm * gadd(m,1,j)
                    qnew(m,i+1,j) = qnew(m,i+1,j) + dtdxp * gadd(m,2,j)
                end do
            end do
            ! these loops are skipped if no ABL present
            do j=1,j_abl_lower
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + ablY_M(j)*( dtdyc * qadd(m,j) &
                    - dtdyc * (fadd(m,j+1) - fadd(m,j)) &
                    - dtdxc * (gadd(m,2,j) - gadd(m,1,j)) )
                    qnew(m,i-1,j) = qnew(m,i-1,j) - dtdxm*ablY_M(j) * gadd(m,1,j)
                    qnew(m,i+1,j) = qnew(m,i+1,j) + dtdxp*ablY_M(j) * gadd(m,2,j)
                end do
            end do
            do j=j_abl_upper,my
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + ablY_M(j)*( dtdyc * qadd(m,j) &
                    - dtdyc * (fadd(m,j+1) - fadd(m,j)) &
                    - dtdxc * (gadd(m,2,j) - gadd(m,1,j)) )
                    qnew(m,i-1,j) = qnew(m,i-1,j) - dtdxm*ablY_M(j) * gadd(m,1,j)
                    qnew(m,i+1,j) = qnew(m,i+1,j) + dtdxp*ablY_M(j) * gadd(m,2,j)
                end do
            end do

        else

        !            # with capa array.
            if (abl_type /= 0) then
                print *, "ABL not implemented yet for capa array"
                stop
            end if
            do j=1,my
                do m=1,num_eqn
                    qnew(m,i,j) = qnew(m,i,j) + dtdy * qadd(m,j) &
                    - (dtdy * (fadd(m,j+1) - fadd(m,j)) &
                    + dtdx * (gadd(m,2,j) - gadd(m,1,j))) &
                    / aux(index_capa,i,j)
                    qnew(m,i-1,j) = qnew(m,i-1,j) - dtdx * gadd(m,1,j) &
                    / aux(index_capa,i-1,j)
                    qnew(m,i+1,j) = qnew(m,i+1,j) + dtdx * gadd(m,2,j) &
                    / aux(index_capa,i+1,j)
                end do
            end do
        endif
    end do    ! End y sweeps

    return
    end subroutine step2

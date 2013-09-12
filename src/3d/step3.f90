
!     ==================================================================
    subroutine step3(maxm,num_eqn,num_waves,num_ghost,mx,my, &
    mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl, &
    qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d, &
    aux1,aux2,aux3,num_aux,work,mwork,use_fwave,rpn3,rpt3,rptt3)
!     ==================================================================

!     # Take one time step, updating q.
!     # On entry, qold and qnew should be identical and give the
!     #    initial data for this step
!     # On exit, qnew returns values at the end of the time step.
!     #    qold is unchanged.

!     # qadd is used to return increments to q from flux3.
!     # fadd, gadd and hadd are used to return flux increments from flux3.
!     # See the flux3 documentation for more information.

    !$ use omp_lib

    implicit real*8(a-h,o-z)
    external rpn3,rpt3,rptt3
    dimension qold(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    dimension qnew(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    dimension  q1d(num_eqn,1-num_ghost:maxm+num_ghost,*)
    dimension qadd(num_eqn,1-num_ghost:maxm+num_ghost,*)
    dimension fadd(num_eqn,1-num_ghost:maxm+num_ghost,*)
    dimension gadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost,*)
    dimension hadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost,*)
    dimension aux(num_aux, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    dimension aux1(num_aux,1-num_ghost:maxm+num_ghost,3,*)
    dimension aux2(num_aux,1-num_ghost:maxm+num_ghost,3,*)
    dimension aux3(num_aux,1-num_ghost:maxm+num_ghost,3,*)
    dimension dtdx1d(1-num_ghost:maxm+num_ghost,*)
    dimension dtdy1d(1-num_ghost:maxm+num_ghost,*)
    dimension dtdz1d(1-num_ghost:maxm+num_ghost,*)
    dimension method(7),mthlim(num_waves)
    dimension work(mwork)

    logical :: use_fwave

!f2py intent(out) cfl
!f2py intent(in,out) qnew
!f2py optional q1d, qadd, fadd, gadd, hadd, dtdx1d, dtdy1d, dtdz1d

! Dummy interfaces just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn3(x)
!f2py x=rpt3(x)
!f2py x=rptt3(x)

    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

!     # partition work array into pieces needed for local storage in
!     # flux3 routine.  Find starting index of each piece:

    nthreads=1    ! Serial
    !$omp parallel
    !$omp single
    !$ nthreads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    nsiz = (maxm+2*num_ghost)*num_eqn
    nsiz_w = nsiz * num_waves
    nsiz_s = (maxm+2*num_ghost)*num_waves

    i0wave     = 1 
    i0s        = i0wave     + nsiz_w*nthreads
    i0amdq     = i0s        + nsiz_s*nthreads
    i0apdq     = i0amdq     + nsiz*nthreads
    i0cqxx     = i0apdq     + nsiz*nthreads
    i0bmamdq   = i0cqxx     + nsiz*nthreads
    i0bmapdq   = i0bmamdq   + nsiz*nthreads
    i0bpamdq   = i0bmapdq   + nsiz*nthreads
    i0bpapdq   = i0bpamdq   + nsiz*nthreads
    i0cmamdq   = i0bpapdq   + nsiz*nthreads
    i0cmapdq   = i0cmamdq   + nsiz*nthreads
    i0cpamdq   = i0cmapdq   + nsiz*nthreads
    i0cpapdq   = i0cpamdq   + nsiz*nthreads
    i0cmamdq2  = i0cpapdq   + nsiz*nthreads
    i0cmapdq2  = i0cmamdq2  + nsiz*nthreads
    i0cpamdq2  = i0cmapdq2  + nsiz*nthreads
    i0cpapdq2  = i0cpamdq2  + nsiz*nthreads
    i0bmcqxxp  = i0cpapdq2  + nsiz*nthreads
    i0bmcqxxm  = i0bmcqxxp  + nsiz*nthreads
    i0bpcqxxp  = i0bmcqxxm  + nsiz*nthreads
    i0bpcqxxm  = i0bpcqxxp  + nsiz*nthreads
    i0cmcqxxp  = i0bpcqxxm  + nsiz*nthreads
    i0cmcqxxm  = i0cmcqxxp  + nsiz*nthreads
    i0cpcqxxp  = i0cmcqxxm  + nsiz*nthreads
    i0cpcqxxm  = i0cpcqxxp  + nsiz*nthreads
    i0bmcmamdq = i0cpcqxxm  + nsiz*nthreads
    i0bmcmapdq = i0bmcmamdq + nsiz*nthreads
    i0bpcmamdq = i0bmcmapdq + nsiz*nthreads
    i0bpcmapdq = i0bpcmamdq + nsiz*nthreads
    i0bmcpamdq = i0bpcmapdq + nsiz*nthreads
    i0bmcpapdq = i0bmcpamdq + nsiz*nthreads
    i0bpcpamdq = i0bmcpapdq + nsiz*nthreads
    i0bpcpapdq = i0bpcpamdq + nsiz*nthreads
    iused      = i0bpcpapdq + nsiz*nthreads - 1

    if (iused > mwork) then
    !        # This shouldn't happen due to checks in claw3
        write(6,*) '*** not enough work space in step3'
        write(6,*) '*** iused = ', iused, '   mwork =',mwork
        stop
    endif



    index_capa = method(6)
!      num_aux = method(7)
    cfl = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz

    !$omp parallel private(me,me1,m,i,j,k,ma,ia,ja,ka,cfl1d,joffset,koffset) reduction(max:cfl)

    me = 0
    !$ me = omp_get_thread_num()
    me1 = me + 1

    if (index_capa == 0) then
    !        # no capa array:
        do 5 i=1-num_ghost,maxm+num_ghost
            dtdx1d(i,me1) = dtdx
            dtdy1d(i,me1) = dtdy
            dtdz1d(i,me1) = dtdz
        5 END DO
    endif

!     # perform x-sweeps
!     ==================

    ! This offset/strided approach keeps race conditions from
    ! happening with OpenMP parallelization.  Because the results of
    ! the flux3 call affect the column flux3 was called on and all
    ! adjacent columns, an easy way to keep threads from interfering
    ! with each other is to make sure that every thread takes a column
    ! that is at least three cells away from all the other threads in
    ! each direction.  Thus, take the parallelized loop with a stride
    ! of 3 on k, synchronize, then repeat at an offset.
    do 51 koffset=0,2

    ! Guided or dynamic scheduling is necessary to keep all the
    ! threads busy if the work per column is very nonuniform.  Dynamic
    ! with a chunk size of 1 seems to work well here, probably because
    ! there's a lot of work per iteration.

    !$omp do schedule(dynamic,1)
    do 50 k = koffset,mz+1,3
        do 50 j = 0,my+1

            forall (m = 1:num_eqn, i = 1-num_ghost:mx+num_ghost)
        !                 # copy data along a slice into 1d array:
            q1d(m,i,me1) = qold(m,i,j,k)
            end forall
        
            if (index_capa > 0)  then
                do 23 i = 1-num_ghost, mx+num_ghost
                    dtdx1d(i,me1) = dtdx / aux(index_capa,i,j,k)
                23 END DO
            endif
        
            if (num_aux > 0)  then
            ! This odd construct may help improve cache locality.
            ! (The F95 standard says each statement in the FORALL
            ! must be executed for all indices before the next
            ! statement is started, so this is different semantically
            ! from putting all three indices in the same FORALL.)
                forall (ka = -1:1)
                forall (ma = 1:num_aux, i = 1-num_ghost:mx+num_ghost)
                aux1(ma, i, 2+ka, me1) = aux(ma, i, j-1, k+ka)
                aux2(ma, i, 2+ka, me1) = aux(ma, i,   j, k+ka)
                aux3(ma, i, 2+ka, me1) = aux(ma, i, j+1, k+ka)
                end forall
                end forall
            endif
        
        !           # compute modifications fadd, gadd and hadd to fluxes along
        !           # this slice:
        
            call flux3(1,maxm,num_eqn,num_waves,num_ghost,mx, &
            q1d(1,1-num_ghost,me1),dtdx1d(1-num_ghost,me1),dtdy,dtdz, &
            aux1(1,1-num_ghost,1,me1),aux2(1,1-num_ghost,1,me1), &
            aux3(1,1-num_ghost,1,me1),num_aux, &
            method,mthlim,qadd(1,1-num_ghost,me1),fadd(1,1-num_ghost,me1), &
            gadd(1,1,-1,1-num_ghost,me1),hadd(1,1,-1,1-num_ghost,me1),cfl1d, &
            work(i0wave+me*nsiz_w),work(i0s+me*nsiz_s),work(i0amdq+me*nsiz), &
            work(i0apdq+me*nsiz),work(i0cqxx+me*nsiz), &
            work(i0bmamdq+me*nsiz),work(i0bmapdq+me*nsiz), &
            work(i0bpamdq+me*nsiz),work(i0bpapdq+me*nsiz), &
            work(i0cmamdq+me*nsiz),work(i0cmapdq+me*nsiz), &
            work(i0cpamdq+me*nsiz),work(i0cpapdq+me*nsiz), &
            work(i0cmamdq2+me*nsiz),work(i0cmapdq2+me*nsiz), &
            work(i0cpamdq2+me*nsiz),work(i0cpapdq2+me*nsiz), &
            work(i0bmcqxxp+me*nsiz),work(i0bpcqxxp+me*nsiz), &
            work(i0bmcqxxm+me*nsiz),work(i0bpcqxxm+me*nsiz), &
            work(i0cmcqxxp+me*nsiz),work(i0cpcqxxp+me*nsiz), &
            work(i0cmcqxxm+me*nsiz),work(i0cpcqxxm+me*nsiz), &
            work(i0bmcmamdq+me*nsiz),work(i0bmcmapdq+me*nsiz), &
            work(i0bpcmamdq+me*nsiz),work(i0bpcmapdq+me*nsiz), &
            work(i0bmcpamdq+me*nsiz),work(i0bmcpapdq+me*nsiz), &
            work(i0bpcpamdq+me*nsiz),work(i0bpcpapdq+me*nsiz), &
            rpn3,rpt3,rptt3,use_fwave)
        
            cfl = dmax1(cfl,cfl1d)
        
        !           # update qnew by flux differencing.
        !           # (rather than maintaining arrays f, g and h for the total fluxes,
        !           # the modifications are used immediately to update qnew
        !           # in order to save storage.)

            if(index_capa == 0)then
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i,me1) &
                - dtdx * (fadd(m,i+1,me1) - fadd(m,i,me1)) &
                - dtdy * (gadd(m,2,0,i,me1) - gadd(m,1,0,i,me1)) &
                - dtdz * (hadd(m,2,0,i,me1) - hadd(m,1,0,i,me1))
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j-1,k)   = qnew(m,i,j-1,k) &
                - dtdy * gadd(m,1,0,i,me1) &
                - dtdz * ( hadd(m,2,-1,i,me1) &
                -   hadd(m,1,-1,i,me1) )
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j-1,k-1) = qnew(m,i,j-1,k-1) &
                - dtdy * gadd(m,1,-1,i,me1) &
                - dtdz * hadd(m,1,-1,i,me1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j,k-1)   = qnew(m,i,j,k-1) &
                - dtdy * ( gadd(m,2,-1,i,me1) &
                -   gadd(m,1,-1,i,me1) ) &
                - dtdz * hadd(m,1,0,i,me1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j+1,k-1) = qnew(m,i,j+1,k-1) &
                + dtdy * gadd(m,2,-1,i,me1) &
                - dtdz * hadd(m,1,1,i,me1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j+1,k)   = qnew(m,i,j+1,k) &
                + dtdy * gadd(m,2,0,i,me1) &
                - dtdz * ( hadd(m,2,1,i,me1) &
                -   hadd(m,1,1,i,me1) )
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j+1,k+1) = qnew(m,i,j+1,k+1) &
                + dtdy * gadd(m,2,1,i,me1) &
                + dtdz * hadd(m,2,1,i,me1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j,k+1)   = qnew(m,i,j,k+1) &
                - dtdy * ( gadd(m,2,1,i,me1) &
                -   gadd(m,1,1,i,me1) ) &
                + dtdz * hadd(m,2,0,i,me1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j-1,k+1) = qnew(m,i,j-1,k+1) &
                - dtdy * gadd(m,1,1,i,me1) &
                + dtdz * hadd(m,2,-1,i,me1)
                end forall
            else
            !              # with capa array
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i,me1) &
                - (dtdx * (fadd(m,i+1,me1) - fadd(m,i,me1)) &
                +  dtdy * (gadd(m,2,0,i,me1) - gadd(m,1,0,i,me1)) &
                +  dtdz * (hadd(m,2,0,i,me1) - hadd(m,1,0,i,me1))) &
                / aux(index_capa,i,j,k)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j-1,k)   = qnew(m,i,j-1,k) &
                - (dtdy * gadd(m,1,0,i,me1) &
                +  dtdz * ( hadd(m,2,-1,i,me1) &
                -   hadd(m,1,-1,i,me1) )) &
                / aux(index_capa,i,j-1,k)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j-1,k-1) = qnew(m,i,j-1,k-1) &
                - (dtdy * gadd(m,1,-1,i,me1) &
                +  dtdz * hadd(m,1,-1,i,me1)) &
                / aux(index_capa,i,j-1,k-1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j,k-1)   = qnew(m,i,j,k-1) &
                - (dtdy * ( gadd(m,2,-1,i,me1) &
                -   gadd(m,1,-1,i,me1) ) &
                +  dtdz * hadd(m,1,0,i,me1)) &
                / aux(index_capa,i,j,k-1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j+1,k-1) = qnew(m,i,j+1,k-1) &
                + (dtdy * gadd(m,2,-1,i,me1) &
                -  dtdz * hadd(m,1,1,i,me1)) &
                / aux(index_capa,i,j+1,k-1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j+1,k)   = qnew(m,i,j+1,k) &
                + (dtdy * gadd(m,2,0,i,me1) &
                - dtdz * ( hadd(m,2,1,i,me1) &
                -   hadd(m,1,1,i,me1) )) &
                / aux(index_capa,i,j+1,k)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j+1,k+1) = qnew(m,i,j+1,k+1) &
                + (dtdy * gadd(m,2,1,i,me1) &
                +  dtdz * hadd(m,2,1,i,me1)) &
                / aux(index_capa,i,j+1,k+1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j,k+1)   = qnew(m,i,j,k+1) &
                - (dtdy * ( gadd(m,2,1,i,me1) &
                -   gadd(m,1,1,i,me1) ) &
                -  dtdz * hadd(m,2,0,i,me1)) &
                / aux(index_capa,i,j,k+1)
                end forall
                forall (m = 1:num_eqn, i = 1:mx)
                qnew(m,i,j-1,k+1) = qnew(m,i,j-1,k+1) &
                - (dtdy * gadd(m,1,1,i,me1) &
                -  dtdz * hadd(m,2,-1,i,me1)) &
                / aux(index_capa,i,j-1,k+1)
                end forall
            endif

    50 END DO
    ! Flush memory (may not be necessary), then make sure everybody
    ! synchronizes before the next offset

    !$omp flush
    !$omp barrier
    51 end do


!     # perform y sweeps
!     ==================


    do 101 koffset=0,2

    !$omp do schedule(dynamic,1)
    do 100 k = koffset, mz+1, 3
        do 100 i = 0, mx+1
        
            forall (m = 1:num_eqn, j = 1-num_ghost:my+num_ghost)
        !                 # copy data along a slice into 1d array:
            q1d(m,j,me1) = qold(m,i,j,k)
            end forall
        
            if (index_capa > 0)  then
                do 71 j = 1-num_ghost, my+num_ghost
                    dtdy1d(j,me1) = dtdy / aux(index_capa,i,j,k)
                71 END DO
            endif
        
            if (num_aux > 0)  then
            ! aux1, aux2, aux3 probably fit in cache, so optimize
            ! access to aux
                forall (ma=1:num_aux,ia= -1:1, j = 1-num_ghost:my+num_ghost)
                aux1(ma, j, 2+ia, me1) = aux(ma, i+ia, j, k-1)
                aux2(ma, j, 2+ia, me1) = aux(ma, i+ia, j,   k)
                aux3(ma, j, 2+ia, me1) = aux(ma, i+ia, j, k+1)
                end forall
            endif
        
        !           # compute modifications fadd, gadd and hadd to fluxes along this
        !           # slice:
        
            call flux3(2,maxm,num_eqn,num_waves,num_ghost,my, &
            q1d(1,1-num_ghost,me1),dtdy1d(1-num_ghost,me1),dtdz,dtdx, &
            aux1(1,1-num_ghost,1,me1),aux2(1,1-num_ghost,1,me1), &
            aux3(1,1-num_ghost,1,me1),num_aux, &
            method,mthlim,qadd(1,1-num_ghost,me1),fadd(1,1-num_ghost,me1), &
            gadd(1,1,-1,1-num_ghost,me1),hadd(1,1,-1,1-num_ghost,me1),cfl1d, &
            work(i0wave+me*nsiz_w),work(i0s+me*nsiz_s),work(i0amdq+me*nsiz), &
            work(i0apdq+me*nsiz),work(i0cqxx+me*nsiz), &
            work(i0bmamdq+me*nsiz),work(i0bmapdq+me*nsiz), &
            work(i0bpamdq+me*nsiz),work(i0bpapdq+me*nsiz), &
            work(i0cmamdq+me*nsiz),work(i0cmapdq+me*nsiz), &
            work(i0cpamdq+me*nsiz),work(i0cpapdq+me*nsiz), &
            work(i0cmamdq2+me*nsiz),work(i0cmapdq2+me*nsiz), &
            work(i0cpamdq2+me*nsiz),work(i0cpapdq2+me*nsiz), &
            work(i0bmcqxxp+me*nsiz),work(i0bpcqxxp+me*nsiz), &
            work(i0bmcqxxm+me*nsiz),work(i0bpcqxxm+me*nsiz), &
            work(i0cmcqxxp+me*nsiz),work(i0cpcqxxp+me*nsiz), &
            work(i0cmcqxxm+me*nsiz),work(i0cpcqxxm+me*nsiz), &
            work(i0bmcmamdq+me*nsiz),work(i0bmcmapdq+me*nsiz), &
            work(i0bpcmamdq+me*nsiz),work(i0bpcmapdq+me*nsiz), &
            work(i0bmcpamdq+me*nsiz),work(i0bmcpapdq+me*nsiz), &
            work(i0bpcpamdq+me*nsiz),work(i0bpcpapdq+me*nsiz), &
            rpn3,rpt3,rptt3,use_fwave)
        
            cfl = dmax1(cfl,cfl1d)
        
        !           # update qnew by flux differencing.
        !           # Note that the roles of the flux updates are changed.
        !           # fadd - modifies the g-fluxes
        !           # gadd - modifies the h-fluxes
        !           # hadd - modifies the f-fluxes
        
            if( index_capa == 0)then
            !               # no capa array.  Standard flux differencing:
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j,me1) &
                - dtdy * (fadd(m,j+1,me1) - fadd(m,j,me1)) &
                - dtdz * (gadd(m,2,0,j,me1) - gadd(m,1,0,j,me1)) &
                - dtdx * (hadd(m,2,0,j,me1) - hadd(m,1,0,j,me1))
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i,j,k+1)   = qnew(m,i,j,k+1) &
                + dtdz * gadd(m,2,0,j,me1) &
                - dtdx * ( hadd(m,2,1,j,me1) &
                -   hadd(m,1,1,j,me1) )
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i+1,j,k+1) = qnew(m,i+1,j,k+1) &
                + dtdz * gadd(m,2,1,j,me1) &
                +  dtdx * hadd(m,2,1,j,me1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i+1,j,k)   = qnew(m,i+1,j,k) &
                - dtdz * ( gadd(m,2,1,j,me1) &
                -   gadd(m,1,1,j,me1) ) &
                + dtdx * hadd(m,2,0,j,me1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i+1,j,k-1) = qnew(m,i+1,j,k-1) &
                - dtdz * gadd(m,1,1,j,me1) &
                + dtdx * hadd(m,2,-1,j,me1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i,j,k-1)   = qnew(m,i,j,k-1) &
                - dtdz * gadd(m,1,0,j,me1) &
                - dtdx * ( hadd(m,2,-1,j,me1) &
                -   hadd(m,1,-1,j,me1) )
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i-1,j,k-1) = qnew(m,i-1,j,k-1) &
                - dtdz * gadd(m,1,-1,j,me1) &
                - dtdx * hadd(m,1,-1,j,me1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i-1,j,k)   = qnew(m,i-1,j,k) &
                - dtdz * ( gadd(m,2,-1,j,me1) &
                -   gadd(m,1,-1,j,me1) ) &
                - dtdx * hadd(m,1,0,j,me1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i-1,j,k+1) = qnew(m,i-1,j,k+1) &
                + dtdz * gadd(m,2,-1,j,me1) &
                -  dtdx*hadd(m,1,1,j,me1)
                end forall
            else
            
            !              #with capa array.
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j,me1) &
                - (dtdy * (fadd(m,j+1,me1) - fadd(m,j,me1)) &
                +  dtdz * (gadd(m,2,0,j,me1) - gadd(m,1,0,j,me1)) &
                +  dtdx * (hadd(m,2,0,j,me1) - hadd(m,1,0,j,me1))) &
                / aux(index_capa,i,j,k)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i,j,k+1)   = qnew(m,i,j,k+1) &
                + (dtdz * gadd(m,2,0,j,me1) &
                -  dtdx * ( hadd(m,2,1,j,me1) &
                -   hadd(m,1,1,j,me1) )) &
                / aux(index_capa,i,j,k+1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i+1,j,k+1) = qnew(m,i+1,j,k+1) &
                + (dtdz * gadd(m,2,1,j,me1) &
                +  dtdx * hadd(m,2,1,j,me1)) &
                / aux(index_capa,i+1,j,k+1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i+1,j,k)   = qnew(m,i+1,j,k) &
                - (dtdz * ( gadd(m,2,1,j,me1) &
                -   gadd(m,1,1,j,me1) ) &
                -  dtdx * hadd(m,2,0,j,me1) ) &
                / aux(index_capa,i+1,j,k)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i+1,j,k-1) = qnew(m,i+1,j,k-1) &
                - (dtdz * gadd(m,1,1,j,me1) &
                -  dtdx * hadd(m,2,-1,j,me1)) &
                / aux(index_capa,i+1,j,k-1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i,j,k-1)   = qnew(m,i,j,k-1) &
                - (dtdz * gadd(m,1,0,j,me1) &
                +  dtdx * ( hadd(m,2,-1,j,me1) &
                -   hadd(m,1,-1,j,me1) )) &
                / aux(index_capa,i,j,k-1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i-1,j,k-1) = qnew(m,i-1,j,k-1) &
                - (dtdz * gadd(m,1,-1,j,me1) &
                +  dtdx * hadd(m,1,-1,j,me1)) &
                / aux(index_capa,i-1,j,k-1)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i-1,j,k)   = qnew(m,i-1,j,k) &
                - (dtdz * ( gadd(m,2,-1,j,me1) &
                -   gadd(m,1,-1,j,me1) ) &
                +  dtdx * hadd(m,1,0,j,me1)) &
                / aux(index_capa,i-1,j,k)
                end forall
                forall (m = 1:num_eqn, j = 1:my)
                qnew(m,i-1,j,k+1) = qnew(m,i-1,j,k+1) &
                + (dtdz * gadd(m,2,-1,j,me1) &
                -  dtdx*hadd(m,1,1,j,me1)) &
                / aux(index_capa,i-1,j,k+1)
                end forall
            endif

        
    100 END DO

    !$omp flush
    !$omp barrier
    101 end do


!     # perform z sweeps
!     ==================


    do 151 joffset=0,2

    !$omp do schedule(dynamic,1)
    do 150 j = joffset, my+1, 3
        do 150 i = 0, mx+1
        
            forall (m = 1:num_eqn, k = 1-num_ghost:mz+num_ghost)
        !                 # copy data along a slice into 1d array:
            q1d(m,k,me1) = qold(m,i,j,k)
            end forall
        
            if (index_capa > 0)  then
                do 130 k = 1-num_ghost, mz+num_ghost
                    dtdz1d(k,me1) = dtdz / aux(index_capa,i,j,k)
                130 END DO
            endif
        
            if (num_aux > 0)  then
            ! See the comment on the X sweeps.  This is semantically
            ! slightly different than putting all the indices in the
            ! same forall, and hopefully better for optimizing access
            ! to aux.
                forall (ja = -1:1, k = 1-num_ghost:mz+num_ghost)
                forall (ma = 1:num_aux)
                aux1(ma, k, 2+ja, me1) = aux(ma, i-1, j+ja, k)
                aux2(ma, k, 2+ja, me1) = aux(ma,   i, j+ja, k)
                aux3(ma, k, 2+ja, me1) = aux(ma, i+1, j+ja, k)
                end forall
                end forall
            endif
        
        !           # compute modifications fadd, gadd and hadd to fluxes along this
        !           # slice:
        
            call flux3(3,maxm,num_eqn,num_waves,num_ghost,mz, &
            q1d(1,1-num_ghost,me1),dtdz1d(1-num_ghost,me1),dtdx,dtdy, &
            aux1(1,1-num_ghost,1,me1),aux2(1,1-num_ghost,1,me1), &
            aux3(1,1-num_ghost,1,me1),num_aux, &
            method,mthlim,qadd(1,1-num_ghost,me1),fadd(1,1-num_ghost,me1), &
            gadd(1,1,-1,1-num_ghost,me1),hadd(1,1,-1,1-num_ghost,me1),cfl1d, &
            work(i0wave+me*nsiz_w),work(i0s+me*nsiz_s),work(i0amdq+me*nsiz), &
            work(i0apdq+me*nsiz),work(i0cqxx+me*nsiz), &
            work(i0bmamdq+me*nsiz),work(i0bmapdq+me*nsiz), &
            work(i0bpamdq+me*nsiz),work(i0bpapdq+me*nsiz), &
            work(i0cmamdq+me*nsiz),work(i0cmapdq+me*nsiz), &
            work(i0cpamdq+me*nsiz),work(i0cpapdq+me*nsiz), &
            work(i0cmamdq2+me*nsiz),work(i0cmapdq2+me*nsiz), &
            work(i0cpamdq2+me*nsiz),work(i0cpapdq2+me*nsiz), &
            work(i0bmcqxxp+me*nsiz),work(i0bpcqxxp+me*nsiz), &
            work(i0bmcqxxm+me*nsiz),work(i0bpcqxxm+me*nsiz), &
            work(i0cmcqxxp+me*nsiz),work(i0cpcqxxp+me*nsiz), &
            work(i0cmcqxxm+me*nsiz),work(i0cpcqxxm+me*nsiz), &
            work(i0bmcmamdq+me*nsiz),work(i0bmcmapdq+me*nsiz), &
            work(i0bpcmamdq+me*nsiz),work(i0bpcmapdq+me*nsiz), &
            work(i0bmcpamdq+me*nsiz),work(i0bmcpapdq+me*nsiz), &
            work(i0bpcpamdq+me*nsiz),work(i0bpcpapdq+me*nsiz), &
            rpn3,rpt3,rptt3,use_fwave)
        
            cfl = dmax1(cfl,cfl1d)
        
        !           # update qnew by flux differencing.
        !           # Note that the roles of the flux updates are changed.
        !           # fadd - modifies the h-fluxes
        !           # gadd - modifies the f-fluxes
        !           # hadd - modifies the g-fluxes
        
            if(index_capa == 0)then
            
            !              #no capa array. Standard flux differencing:
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k,me1) &
                - dtdz * (fadd(m,k+1,me1) - fadd(m,k,me1)) &
                - dtdx * (gadd(m,2,0,k,me1) - gadd(m,1,0,k,me1)) &
                - dtdy * (hadd(m,2,0,k,me1) - hadd(m,1,0,k,me1))
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i,j+1,k)   = qnew(m,i,j+1,k) &
                - dtdx * ( gadd(m,2,1,k,me1) &
                -   gadd(m,1,1,k,me1) ) &
                + dtdy * hadd(m,2,0,k,me1)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i+1,j+1,k) = qnew(m,i+1,j+1,k) &
                + dtdx * gadd(m,2,1,k,me1) &
                + dtdy * hadd(m,2,1,k,me1)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i+1,j,k)   = qnew(m,i+1,j,k) &
                + dtdx * gadd(m,2,0,k,me1) &
                - dtdy * ( hadd(m,2,1,k,me1) &
                -   hadd(m,1,1,k,me1) )
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i+1,j-1,k) = qnew(m,i+1,j-1,k) &
                + dtdx * gadd(m,2,-1,k,me1) &
                - dtdy * hadd(m,1,1,k,me1)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i,j-1,k)   = qnew(m,i,j-1,k) &
                - dtdx * ( gadd(m,2,-1,k,me1) &
                -   gadd(m,1,-1,k,me1) ) &
                - dtdy * hadd(m,1,0,k,me1)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i-1,j-1,k) = qnew(m,i-1,j-1,k) &
                - dtdx * gadd(m,1,-1,k,me1) &
                - dtdy * hadd(m,1,-1,k,me1)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i-1,j,k)   = qnew(m,i-1,j,k) &
                - dtdx * gadd(m,1,0,k,me1) &
                - dtdy * ( hadd(m,2,-1,k,me1) &
                -   hadd(m,1,-1,k,me1) )
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i-1,j+1,k) = qnew(m,i-1,j+1,k) &
                - dtdx * gadd(m,1,1,k,me1) &
                + dtdy * hadd(m,2,-1,k,me1)
                end forall
            else
            
            !              # with capa array
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k,me1) &
                - (dtdz * (fadd(m,k+1,me1) - fadd(m,k,me1)) &
                +  dtdx * (gadd(m,2,0,k,me1) - gadd(m,1,0,k,me1)) &
                +  dtdy * (hadd(m,2,0,k,me1) - hadd(m,1,0,k,me1))) &
                / aux(index_capa,i,j,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i,j+1,k) = qnew(m,i,j+1,k) &
                - (dtdx * ( gadd(m,2,1,k,me1) &
                -   gadd(m,1,1,k,me1) ) &
                -  dtdy * hadd(m,2,0,k,me1)) &
                / aux(index_capa,i,j+1,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i+1,j+1,k) = qnew(m,i+1,j+1,k) &
                + (dtdx * gadd(m,2,1,k,me1) &
                +  dtdy * hadd(m,2,1,k,me1)) &
                / aux(index_capa,i+1,j+1,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i+1,j,k)   = qnew(m,i+1,j,k) &
                + (dtdx * gadd(m,2,0,k,me1) &
                - dtdy * ( hadd(m,2,1,k,me1) &
                -   hadd(m,1,1,k,me1) )) &
                / aux(index_capa,i+1,j,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i+1,j-1,k) = qnew(m,i+1,j-1,k) &
                + (dtdx * gadd(m,2,-1,k,me1) &
                -  dtdy * hadd(m,1,1,k,me1)) &
                / aux(index_capa,i+1,j-1,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i,j-1,k)   = qnew(m,i,j-1,k) &
                - (dtdx * ( gadd(m,2,-1,k,me1) &
                -   gadd(m,1,-1,k,me1) ) &
                +  dtdy * hadd(m,1,0,k,me1)) &
                / aux(index_capa,i,j-1,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i-1,j-1,k) = qnew(m,i-1,j-1,k) &
                - (dtdx * gadd(m,1,-1,k,me1) &
                +  dtdy * hadd(m,1,-1,k,me1)) &
                / aux(index_capa,i-1,j-1,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i-1,j,k)   = qnew(m,i-1,j,k) &
                - (dtdx * gadd(m,1,0,k,me1) &
                + dtdy * ( hadd(m,2,-1,k,me1) &
                -   hadd(m,1,-1,k,me1) )) &
                / aux(index_capa,i-1,j,k)
                end forall
                forall (m = 1:num_eqn, k = 1:mz)
                qnew(m,i-1,j+1,k) = qnew(m,i-1,j+1,k) &
                - (dtdx * gadd(m,1,1,k,me1) &
                -  dtdy * hadd(m,2,-1,k,me1)) &
                / aux(index_capa,i-1,j+1,k)
                end forall
            endif

    150 END DO

    !$omp flush
    !$omp barrier
    151 end do
    !$omp end parallel

    return

    end subroutine step3



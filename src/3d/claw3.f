c
c
c     ==================================================================
      subroutine claw3(meqn,mwaves,mbc,mx,my,mz,maux,
     &           q,aux,xlower,ylower,zlower,dx,dy,dz,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,use_fwave,info,bc3,rpn3,rpt3,rptt3,src3,
     &           b4step3)
c     ==================================================================
c
c  Solves a hyperbolic system of conservation laws in three space
c  dimensions of the general form
c
c     capa * q_t + A q_x + B q_y + C q_z = psi
c
c  The "capacity function" capa(x,y,z) and source term psi are optional
c  (see below).
c
c  For a more complete description see the documentation at
c      http://www.amath.washington.edu/~claw
c
c  --------------------------------------------------------
c
c  The user must supply the following subroutines:
c
c    bc3, rpn3, rpt3, rptt3  subroutines specifying the boundary conditions
c                            and Riemann solvers.
c
c    b4step3            The routine b4step3 is called each time step and
c                       can be supplied by the user in order to perform
c                       other operations that are necessary every time
c                       step.  For example, if the variables stored in
c                       the aux arrays are time-dependent then these
c                       values can be set.
c
c
c  In addition, if the equation contains source terms psi, then the user
c  must provide:
c
c    src3               subroutine that solves capa * q_t = psi
c                       over a single time step.
c
c  These routines must be declared EXTERNAL in the main program.
c  For description of the calling sequences, see below.
c
c  Dummy routines b4step3.f and src3.f are available in
c       claw/clawpack/3d/lib
c
c  A subroutine implementing many standard boundary conditions is
c  available in claw/clawpack/3d/lib/bc3.f.
c
c
c
c  Description of parameters...
c  ----------------------------
c
c    meqn is the number of equations in the system of
c         conservation laws.
c
c    mwaves is the number of waves that result from the
c           solution of each Riemann problem.  Often mwaves = meqn but
c           for some problems these may be different, e.g. for the Euler
c           equations meqn = 5 but mwaves = 3 since there are only 3
c           distinct wave speeds.
c
c    mbc is the number of "ghost cells" that must be added on to each
c       side of the domain to handle boundary conditions.  The cells
c       actually in the physical domain are labelled from 1 to mx in x,
c       from 1 to my in y, and from 1 to mz in z.
c       The arrays are dimensioned actually indexed
c       from 1-mbc to mx+mbc, from 1-mbc to my+mbc and from 1-mbc
c       to mz+mbc.
c       For the methods currently implemented, mbc = 2 should be used.
c       If the user implements another method that has a larger stencil and
c       hence requires more ghost cells, a larger value of mbc could be used.
c       q is extended from the physical domain to the ghost cells by the
c       user-supplied routine bc2.
c
c    mx is the number of grid cells in the x-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c
c    my is the number of grid cells in the y-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c
c    mz is the number of grid cells in the z-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c
c    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
c        On input:  initial data at time tstart.
c        On output: final solution at time tend.
c        q(m,i,j,k) = value of mth component in the (i,j,k) cell.
c        Values within the physical domain are in q(m,i,j,k)
c                for i = 1,2,...,mx , j = 1,2,...,my, and
c                    k = 1,2,...,mz,
c        mbc extra cells on each end are needed for boundary conditions
c        as specified in the routine bc3.
c
c    aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:my+mbc)
c        Array of auxiliary variables that are used in specifying the problem.
c        If method(7) = 0 then there are no auxiliary variables and aux
c                         can be a dummy variable.
c        If method(7) = maux > 0 then there are maux auxiliary variables
c                         and aux must be dimensioned as above.
c
c        Capacity functions are one particular form of auxiliary variable.
c        These arise in some applications, e.g. the
c        determinant of the Jacobian if a mapped grid is used, or a density
c        or porosity function in some advection problems.
c        See Clawpack Note # 5 for examples.
c
c        If method(6) = 0 then there is no capacity function.
c        If method(6) = mcapa > 0  then there is a capacity function and
c            capa(i,j,k), the "capacity" of the (i,j,k) cell, is assumed to be
c            stored in aux(mcapa,i,j,k).
c            In this case we require method(7).ge.mcapa.
c
c    dx = grid spacing in x.
c         (for a computation in ax <= x <= bx,  set dx = (bx-ax)/mx.)
c
c    dy = grid spacing in y.
c         (for a computation in ay <= y <= by,  set dy = (by-ay)/my.)
c
c    dz = grid spacing in z.
c         (for a computation in az <= z <= bz,  set dz = (bz-az)/mz.)
c
c    tstart = initial time.
c
c    tend = Desired final time (on input).
c              If tend<tstart, then claw3 returns after a single successful
c                 time step has been taken (single-step mode).
c              Otherwise, as many steps are taken as needed to reach tend,
c                 up to a maximum of nv(1).
c         = Actual time reached (on output).
c
c    dtv(1:5) = array of values related to the time step:
c               (Note: method(1)=1 indicates variable size time steps)
c         dtv(1) = value of dt to be used in all steps if method(1) = 0
c                = value of dt to use in first step if method(1) = 1 
c         dtv(2) = unused if method(1) = 0.
c                = maximum dt allowed if method(1) = 1.
c         dtv(3) = smallest dt used (on output)
c         dtv(4) = largest dt used (on output)
c         dtv(5) = dt used in last step (on output)
c
c    cflv(1:4) = array of values related to Courant number:
c         cflv(1) = maximum Courant number to be allowed.
c                   With variable time steps the step is retracted and a
c                   smaller step taken if the Courant
c                   number is larger than this value.
c                   With fixed time steps the routine aborts.
c                   Usually cflv(1)=1.0 should work
c         cflv(2) = unused if method(1) = 0.
c                 = desired Courant number if method(1) = 1.
c                   Should be somewhat less than cflv(1), e.g. 0.9
c         cflv(3) = largest Courant number observed (on output).
c         cflv(4) = Courant number in last step (on output).
c
c    nv(1:2) = array of values related to the number of time steps:
c         nv(1) = unused if method(1) = 0
c               = maximum number of time steps allowed if method(1) = 1
c         nv(2) = number of time steps taken (on output).
c
c    method(1:7) = array of values specifying the numerical method to use
c                  and also indicating whether source terms, capacity
c                  function, auxiliary variables are present in the equation.
c
c         method(1) = 0 if fixed size time steps are to be taken.
c                       In this case, dt = dtv(1) in all steps.
c                   = 1 if variable time steps are to be used.
c                       In this case, dt = dtv(1) in the first step and
c                       thereafter the value cflv(2) is used to choose the
c                       next time step based on the maximum wave speed seen
c                       in the previous step.  Note that since this value
c                       comes from the previous step, the Courant number will
c                       not in general be exactly equal to the desired value
c                       If the actual Courant number in the next step is
c                       greater than cflv(1), then this step is redone with a
c                       smaller dt.
c
c         method(2) = 1 if only first order increment waves are to be used.
c                   = 2 if second order correction terms are to be added, with
c                       a flux limiter as specified by mthlim.
c
c         method(3) <  0 Gives dimensional splitting using Godunov
c                        splitting, i.e. formally first order
c                        accurate.
c                      0 Gives the Donor cell method. No transverse
c                        propagation of neither the increment wave
c                        nor the correction wave.
c                   = 10 Transverse propagation of the increment wave
c                        as in 2D. Note that method (2,10) is
c                        unconditionally unstable.
c                   = 11 Corner transport upwind of the increment
c                        wave. Note that method (2,11) also is
c                        unconditionally unstable.
c                   = 20 Both the increment wave and the correction
c                        wave propagate as in the 2D case. Only to
c                        be used with method(2) = 2.
c                   = 21 Corner transport upwind of the increment wave,
c                        and the correction wave propagates as in 2D.
c                        Only to be used with method(2) = 2.
c                   = 22 3D propagation of both the increment wave and
c                        the correction wave. Only to be used with
c                        method(2) = 2.
c
c         method(4) = 0 to suppress printing
c                   = 1 to print dt and Courant number every time step
c
c         method(5) = 0 if there is no source term psi.  In this case
c                       the subroutine src2 is never called so a dummy
c                       parameter can be given.
c                   = 1 if there is a source term.  In this case
c                       the subroutine src2 must be provided.
c
c         method(6) = 0 if there is no capacity function capa.
c                   = mcapa > 0 if there is a capacity function.  In this case
c                       aux(mcapa,i,j,k) is the capacity of cell (i,j,k)
c                       and you must also specify method(7) .ge. mcapa
c                       and set aux.
c
c         method(7) = 0 if there is no aux array used.
c                   = maux > 0  if there are maux auxiliary variables.
c
c         The recommended choice of methods for most problems is
c            method(1) = 1,  method(2) = 2,  method(3) = 22.
c
c    mthlim(1:mwaves) = array of values specifying the flux limiter to be used
c                     in each wave family mw.  Often the same value will be used
c                     for each value of mw, but in some cases it may be
c                     desirable to use different limiters.  For example,
c                     for the Euler equations the superbee limiter might be
c                     used for the contact discontinuity (mw=2) while another
c                     limiter is used for the nonlinear waves.  Several limiters
c                     are built in and others can be added by modifying the
c                     subroutine philim.
c
c        mthlim(mw) = 0 for no limiter
c                   = 1 for minmod
c                   = 2 for superbee
c                   = 3 for van Leer
c                   = 4 for monotonized centered
c
c    mthbc(1:6) = array of values specifying what boundary conditions should
c                 be used at each edge of the domain, if the standard
c                 bc3.f routine is used.  Passed to bc3.
c
c
c    work(mwork) = double precision work array of length at least mwork
c
c    mwork = length of work array.  Must be at least
c               N * (mx + 2*mbc)*(my + 2*mbc)*(mz + 2*mbc)*meqn
c               + (max(mx,my) + 2*mbc) * (46*meqn + mwaves + meqn*mwaves
c                                          + 9*maux + 2)
c            where N = 1 if method(5)=0  (no source term)
c                  N = 2 if method(5)=1  (source term)
c            If mwork is too small then the program returns with info = 4
c            and also prints the required value of mwork to unit 6.
c
c
c    info = output value yielding error information:
c         = 0 if normal return.
c         = 1 if mbc is too small.
c         = 2 if method(1)=0 and dt doesn't divide (tend - tstart).
c         = 3 if method(1)=1 and cflv(2) > cflv(1).
c         = 4 if mwork is too small.
c         = 5 if method(6) > method(7)
c         = 6 if illegal parameters are used in specifying the method.
c         = 11 if the code attempted to take too many time steps, n > nv(1).
c              This could only happen if method(1) = 1 (variable time steps).
c         = 12 if the method(1)=0 and the Courant number is greater than 1
c              in some time step.
c
c           Note: if info.ne.0, then tend is reset to the value of t actually
c           reached and q contains the value of the solution at this time.
c
c    User-supplied subroutines
c    -------------------------
c
c    bc3 = subroutine that specifies the boundary conditions.
c          This subroutine should extend the values of q from cells
c          (1:mx, 1:my, 1:mz) to the mbc ghost cells along each edge
c          of the domain.  The routine claw/clawpack/3d/lib/bc3.f
c          implements many standard choices.
c
c
c    rpn3 = user-supplied subroutine that implements the Riemann solver
c           along a one-dimensional slice of data.
c
c          The form of this subroutine is
c  ---------------------------------------------------------------------
c     subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,
c    &                  auxl,auxr,maux,wave,s,amdq,apdq)
c    implicit double precision (a-h,o-z)
c     dimension wave(mwaves, meqn, 1-mbc:maxm+mbc)
c     dimension    s(mwaves, 1-mbc:maxm+mbc)
c     dimension   ql(meqn, 1-mbc:maxm+mbc)
c     dimension   qr(meqn, 1-mbc:maxm+mbc)
c     dimension amdq(meqn, 1-mbc:maxm+mbc)
c     dimension apdq(meqn, 1-mbc:maxm+mbc)
c     dimension auxl(maux, 1-mbc:maxm+mbc)
c     dimension auxr(maux, 1-mbc:maxm+mbc)
c  ---------------------------------------------------------------------
c         solve Riemann problems along one slice of data.
c         This data is along a slice in the x-direction if ixyz=1
c                                       the y-direction if ixyz=2.
c                                       the z-direction if ixyz=3.
c
c         On input, ql contains the state vector at the left edge of each cell
c                   qr contains the state vector at the right edge of each cell
c
c         auxl(ma,i) contains auxiliary data for cells along this slice,
c            where ma=1,maux in the case where maux=method(7) > 0.
c
c         Note that the i'th Riemann problem has left state qr(:,i-1)
c                                            and right state ql(:,i)
c         From the basic clawpack routines, this routine is called with ql = qr
c
c         More flexibility is allowed
c         in case the user wishes to implement another solution method
c         that requires left and rate states at each interface.
c
c
c          On output,
c             wave(m,mw,i) is the mth component of the jump across
c                              wave number mw in the ith Riemann problem.
c             s(mw,i) is the wave speed of wave number mw in the
c                              ith Riemann problem.
c             amdq(m,i) is the m'th component of the left-going flux
c                       difference.
c             apdq(m,i) is the m'th component of the right-going flux
c                       difference.
c           It is assumed that each wave consists of a jump discontinuity
c           propagating at a single speed, as results, for example, from a
c           Roe approximate Riemann solver.  An entropy fix can be included
c           into the specification of amdq and apdq.
c
c
c    rpt3 = user-supplied subroutine that implements the splitting of
c           a flux difference asdq into waves in the transverse direction.
c           The form of this subroutine is
c  ---------------------------------------------------------------------
c     subroutine rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx,
c    &                  ql,qr,aux1,aux2,aux3,maux,asdq,
c    &                  bmasdq,bpasdq)
c     implicit double precision (a-h,o-z)
c     dimension     ql(meqn, 1-mbc:maxm+mbc)
c     dimension     qr(meqn, 1-mbc:maxm+mbc)
c     dimension   asdq(meqn, 1-mbc:maxm+mbc)
c     dimension bmasdq(meqn, 1-mbc:maxm+mbc)
c     dimension bpasdq(meqn, 1-mbc:maxm+mbc)
c     dimension   aux1(maux, 1-mbc:maxm+mbc, 3)
c     dimension   aux2(maux, 1-mbc:maxm+mbc, 3)
c     dimension   aux3(maux, 1-mbc:maxm+mbc, 3)
c  -------------------------------------------------
c       On input,
c
c          ql,qr is the data along some one-dimensional slice, as in rpn3
c               This slice is
c                   in the x-direction if ixyz=1,
c                   in the y-direction if ixyz=2, or
c                   in the z-direction if ixyz=3.
c          asdq is an array of flux differences (A^* \Delta q).
c               asdq(i,:) is the flux difference propagating away from
c               the interface between cells i-1 and i.
c          imp = 1 if asdq = A^- \Delta q,  the left-going flux difference
c                2 if asdq = A^+ \Delta q, the right-going flux difference
c
c          aux2 is the auxiliary array (if method(7)=maux>0) along
c               the plane where this slice lies, say at j=J if ixyz=1.
c               aux2(:,:,1) contains data along j=J, k=k-1
c               aux2(:,:,2) contains data along j=J, k=k
c               aux2(:,:,3) contains data along j=J, k=k+1
c          aux1 is the auxiliary array along the plane with j=J-1
c          aux3 is the auxiliary array along the plane with j=J+1
c
c            if ixyz=2 then aux2 is in some plane k=K, and
c               aux2(:,:,1)  contains data along i=I-1, k=K, etc.
c
c            if ixyz=3 then aux2 is in some plane i=I, and
c               aux2(:,:,1)  contains data along j=j-1, i=I, etc.
c
c       On output,

c       If data is in x-direction (ixyz=1) then this routine does the
c       splitting of  asdq (= A^* \Delta q, where * = + or -) ...
c
c       into down-going flux difference bmasdq (= B^- A^* \Delta q)
c          and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c          when icoor = 2,
c
c       or
c
c       into down-going flux difference bmasdq (= C^- A^* \Delta q)
c          and up-going flux difference bpasdq (= C^+ A^* \Delta q)
c          when icoor = 3.
c
c
c       More generally, ixyz specifies what direction the slice of data is
c       in, and icoor tells which transverse direction to do the splitting in:
c
c       If ixyz = 1,  data is in x-direction and then
c             icoor = 2  =>  split in the y-direction  (iuvw=2)
c             icoor = 3  =>  split in the z-direction  (iuvw=3)
c
c       If ixyz = 2,  data is in y-direction and then
c             icoor = 2  =>  split in the z-direction  (iuvw=3)
c             icoor = 3  =>  split in the x-direction  (iuvw=1)
c
c       If ixyz = 3,  data is in z-direction and then
c             icoor = 2  =>  split in the x-direction  (iuvw=1)
c             icoor = 3  =>  split in the y-direction  (iuvw=2)
c
c           For example, for a linear system q_t + Aq_x + Bq_y + Cq_z = 0,
c                   asdq = A^+ dq  or  A^- dq
c                   and this is then split into
c                       bmasdq = B^- asdq   and   bpasdq = B^+ asdq
c                   when ixyz = 1, and icoor = 2.
c
c   rptt3 = user-supplied subroutine that implements the splitting of
c           a flux difference bsasdq into waves in the double transverse
c           direction.
c           This subroutine has the form :
c  ---------------------------------------------------------------------
c     subroutine rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,
c    &                  ql,qr,aux1,aux2,aux3,bsasdq,
c    &                  cmbsasdq,cpbsasdq)
c     implicit double precision (a-h,o-z)
c     dimension      ql(meqn, 1-mbc:maxm+mbc)
c     dimension      qr(meqn, 1-mbc:maxm+mbc)
c     dimension   bsasdq(meqn, 1-mbc:maxm+mbc)
c     dimension cmbsasdq(meqn, 1-mbc:maxm+mbc)
c     dimension cpbsasdq(meqn, 1-mbc:maxm+mbc)
c     dimension   aux1(maux, 1-mbc:maxm+mbc, 3)
c     dimension   aux2(maux, 1-mbc:maxm+mbc, 3)
c     dimension   aux3(maux, 1-mbc:maxm+mbc, 3)
c  -------------------------------------------------
c       On input,
c
c          ql,qr is the data along some one-dimensional slice, as in rpn3
c               This slice is
c                   in the x-direction if ixyz=1,
c                   in the y-direction if ixyz=2, or
c                   in the z-direction if ixyz=3.
c          bsasdq is an array of flux differences (B^* A^* \Delta q).
c               bsasdq(i,:) is the flux difference propagating into cells
c               who share only a corner with interface between cells i and i-1.
c          imp = 1 if asdq = A^- \Delta q,  the left-going flux difference
c                2 if asdq = A^+ \Delta q, the right-going flux difference
c
c        impt = 1 if bsasdq = B^-A^* \Delta q,
c               2 if bsasdq = B^+A^* \Delta q
c
c          aux2 is the auxiliary array (if method(7)=maux>0) along
c               the plane where this slice lies, say at j=J if ixyz=1.
c               aux2(:,:,1) contains data along j=J, k=k-1
c               aux2(:,:,2) contains data along j=J, k=k
c               aux2(:,:,3) contains data along j=J, k=k+1
c          aux1 is the auxiliary array along the plane with j=J-1
c          aux3 is the auxiliary array along the plane with j=J+1
c
c            if ixyz=2 then aux2 is in some plane k=K, and
c               aux2(:,:,1)  contains data along i=I-1, k=K, etc.
c
c            if ixyz=3 then aux2 is in some plane i=I, and
c               aux2(:,:,1)  contains data along j=j-1, i=I, etc.
c
c       On output,

c       If data is in x-direction (ixyz=1) then this routine does the
c       splitting of  bsasdq (= B^*A^* \Delta q, where * = + or -) ...
c
c       into down-going flux difference cmbsasdq (= C^- B^* A^* \Delta q)
c          and up-going flux difference cpbsasdq (= c^+ B^* A^* \Delta q)
c          when icoor = 2,
c
c       or
c
c       into front-going (in negative y direction) flux difference
c          bmcsasdq (= B^- C^*A^* \Delta q)
c          and back-going (in positive y direction) flux difference
c          bpcsasdq (= B^+ C^* A^* \Delta q) when icoor = 3.
c
c
c       More generally, ixyz specifies what direction the slice of data is
c       in, and icoor tells which transverse direction to do the splitting in:
c
c       If ixyz = 1,  data is in x-direction and then
c             icoor = 2  =>  split in the y-direction  (iuvw=2)
c             icoor = 3  =>  split in the z-direction  (iuvw=3)
c
c       If ixyz = 2,  data is in y-direction and then
c             icoor = 2  =>  split in the z-direction  (iuvw=3)
c             icoor = 3  =>  split in the x-direction  (iuvw=1)
c
c       If ixyz = 3,  data is in z-direction and then
c             icoor = 2  =>  split in the x-direction  (iuvw=1)
c             icoor = 3  =>  split in the y-direction  (iuvw=2)
c
c           For example, for a linear system q_t + Aq_x + Bq_y + Cq_z = 0,
c                   bsasdq = B^+ A^* dq  or  B^- A^* dq
c                   and this is then split into
c                       cmbsasdq = C^- bsasdq   and   cpbsasdq = C^+ bsasdq
c                   when ixyz = 1, and icoor = 3.
c
c           Or, if ixyz = 1, and icoor = 2, then
c                   bsasdq = C^+ A^* dq or C^- A^* dq
c                   and this is then split into
c                       cmbsasdq = B^- bsasdq and cpbsasdq = B^+ bsasdq
c
c
c
c    src3 = user-supplied subroutine that takes one time step on the
c           source terms alone, solving
c               capa * q_t = psi
c           over time dt.
c
c           If method(5)=0 then the equation does not contain a source
c           term and this routine is never called.  A dummy argument can
c           be used with many compilers, or provide a dummy subroutine that
c           does nothing.
c
c           The form of this subroutine is
c  -------------------------------------------------
c      subroutine src3(meqn,mbc,mx,my,mz,q,aux,t,dt)
c      implicit double precision (a-h,o-z)
c      dimension    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
c     &               1-mbc:mz+mbc)
c      dimension aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
c     &              1-mbc:mz+mbc)
c  -------------------------------------------------
c      If method(7)=0  or the auxiliary variables are not needed in
c      this solver,
c      then the latter dimension statement can be omitted, but aux should
c      still appear in the argument list.
c
c      On input, q(m,i,j,k) contains the data for solving the
c                source term equation.
c      On output, q(m,i,j,k) should have been replaced by the solution to
c                 the source term equation after a step of length dt.
c
c      b4step3 = subroutine that is called from claw3 before each call to
c                step3.  Use to set time-dependent aux arrays or perform
c                other tasks which must be done every time step.
c
c          The form of this subroutine is
c  -------------------------------------------------
c      subroutine b4step3(mbc,mx,my,mz,meqn,q,
c    &            xlower,ylower,zlower,dx,dy,dz,time,dt,maux,aux)
c      implicit double precision (a-h,o-z)
c      dimension     q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
c     &                1-mbc:mz+mbc)
c      dimension   aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
c     &                1-mbc:mz+mbc)
c  -------------------------------------------------

c
c
c =========================================================================
c
c  Copyright 1996 -- 2002 R. J. LeVeque and J. O. Langseth
c
c  This software is made available for research and instructional use only.
c  You may copy and use this software without charge for these non-commercial
c  purposes, provided that the copyright notice and associated text is
c  reproduced on all copies.  For all other uses (including distribution of
c  modified versions), please contact the author at the address given below.
c
c  *** This software is made available "as is" without any assurance that it
c  *** will work for your purposes.  The software may in fact have defects, so
c  *** use the software at your own risk.
c
c  --------------------------------------
c    CLAWPACK Version 4.1,  August, 2002
c    Webpage: http://www.amath.washington.edu/~claw
c  --------------------------------------
c
c    Authors: Randall J. LeVeque
c             Applied Mathematics
c             Box 352420
c             University of Washington,
c             Seattle, WA 98195-2420
c             rjl@amath.washington.edu
c
c      and
c             Jan Olav Langseth
c             FFI/VM
c             Box 25
c             N-2007 Kjeller
c             Norway
c             jol@ffi.no
c
c     Modified in 2002 by
c             Donna Calhoun
c             Applied Mathematics
c             Box 352420
c             University of Washington,
c             Seattle, WA 98195-2420
c             calhoun@amath.washington.edu
c
c    ======================================================================
c    Beginning of claw3 code
c    ======================================================================
c
      implicit real*8(a-h,o-z)

c     F77 way to get OpenMP function declarations
c$    include 'omp_lib.h'

      external bc3,rpn3,rpt3, rptt3
      dimension     q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
      dimension   aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
      dimension work(mwork)
      dimension mthlim(mwaves),method(7),dtv(5),cflv(4),nv(2)
c
      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
c
      maxm = max0(mx, my, mz)
      info = 0
      t = tstart
      maxn = nv(1)
      dt = dtv(1)   !# initial dt
      cflmax = 0.d0
      dtmin = dt
      dtmax = dt
      nv(2) = 0
c      maux = method(7)
c
c     # check for errors in data:
c     ---------------------------
c
      call chkmth(method,info)
      if( info .eq. 6) go to 900
c
      if (mbc.lt.2) then
         info = 1
         go to 900
         endif
c
      if (method(1) .eq. 0) then
c        # fixed size time steps.  Compute the number of steps:
         if (tend .lt. tstart) then
c             # single step mode
              maxn = 1
           else
              maxn = nint((tend - tstart)/dt)    ! Round to nearest int
              if (dabs(maxn*dt - (tend-tstart)) .gt.
     &                          1d-5*(tend-tstart)) then
c                # dt doesn't divide time interval integer number of times
                 info = 2
                 write(6,*) 
     &               'CLAW3 ERROR... dt does not divide (tend-tstart)'
                 go to 900
                 endif
           endif
         endif

c
      if (method(1).eq.1 .and. cflv(2).gt.cflv(1)) then
         info = 3
         write(6,*) 'CLAW3 ERROR...  cflv(2) > cflv(1)'
         go to 900
         endif
c
      if (method(6).gt.method(7)) then
         info = 5
         write(6,*) 'CLAW3 ERROR...  method(6) > method(7)'
         go to 900
         endif
c
      if (method(5).lt.2) then
          narray = 1   !# no source terms -- only need one qwork array
        else
          narray = 2           !# need two qwork arrays
       endif
c
       nthreads=1               ! Serial
c$omp  parallel
c$omp  single
c$     nthreads = omp_get_num_threads()
c$     maxthreads = omp_get_max_threads()
c$omp  end single
c$omp  end parallel
       mwork0 = (maxm+2*mbc)*(46*meqn + mwaves + meqn*mwaves
     &                      + 9*maux + 3)*nthreads
     &          + narray * (mx + 2*mbc) * (my + 2*mbc)
     &                   * (mz + 2*mbc) * meqn
c
      if (mwork .lt. mwork0) then
         info = 4
         write(6,*) 'CLAW3 ERROR... mwork should be increased to ',
     &               mwork0
         go to 900
         endif
c
c     # partition work array into pieces needed for local storage in
c     # step2 routine. Find starting index of each piece:
c
      i0qadd = 1
      i0fadd = i0qadd + (maxm+2*mbc)*meqn*nthreads
      i0gadd = i0fadd + (maxm+2*mbc)*meqn*nthreads
      i0hadd = i0gadd + 6*(maxm+2*mbc)*meqn*nthreads
      i0q1d = i0hadd + 6*(maxm+2*mbc)*meqn*nthreads
      i0dtdx1d = i0q1d + (maxm+2*mbc)*meqn*nthreads
      i0dtdy1d = i0dtdx1d + (maxm+2*mbc)*nthreads
      i0dtdz1d = i0dtdy1d + (maxm+2*mbc)*nthreads
      i0qwrk1 = i0dtdz1d + (maxm+2*mbc)*nthreads
c
      nqwork = (mx + 2*mbc) * (my + 2*mbc)
     &       * (mz + 2*mbc) * meqn  !# size of q array
      if (method(5).lt.2) then
          i0qwrk2 = i0qwrk1  !# qwrk2 points to same storage as qwrk1
        else
          i0qwrk2 = i0qwrk1 + nqwork  !# second qwork array with source term
        endif
c
      i0aux1 = i0qwrk2 + nqwork
      i0aux2 = i0aux1 + 3*(maxm+2*mbc)*maux*nthreads
      i0aux3 = i0aux2 + 3*(maxm+2*mbc)*maux*nthreads
c
      i0next = i0aux3 + 3*(maxm+2*mbc)*maux*nthreads  !# next free space
      mused = i0next - 1                     !# space already used
      mwork1 = mwork - mused              !# remaining space (passed to step2)
c
c
c     -----------
c     # main loop
c     -----------
c
      if (maxn.eq.0) go to 900
      do 100 n=1,maxn
         told = t   !# time at beginning of time step.

c        # adjust dt to hit tend exactly if we're near end of computation
c        #  (unless tend < tstart, which is a flag to take only a single step)
         if (told+dt.gt.tend .and. tstart.lt.tend) dt = tend - told

c
   40    continue
c
c        # store dt and t in the common block comxyt in case they are needed
c        # in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt
         dxcom = dx
         dycom = dy
         dzcom = dz
c
c
c
c        ================================================================
c
c        -------------------------
c        # main steps in algorithm
c        -------------------------
c
c        # extend data from grid to bordering boundary cells:
c
         call bc3(meqn,mbc,mx,my,mz,xlower,
     &               ylower,zlower,dx,dy,dz,q,maux,aux,told,dt,mthbc)
c
c        # call user-supplied routine which might set aux arrays
c        # for this time step, for example.
         call b4step3(mbc,mx,my,mz,meqn,q,
     &            xlower,ylower,zlower,dx,dy,dz,told,dt,maux,aux)
c
c
         if (method(5).eq.2) then
c            # with source term:   use Strang splitting
c            # First need to store solution before taking
c            # step on source terms in case we need to redo everything with
c            # a smaller time step if the Courant number is too large in
c            # subroutine step3.
c
              call copyq3(meqn,mbc,mx,my,mz,q,
     &                   work(i0qwrk2))
c
c            # source terms over a half time step:
             dt2 = dt / 2.d0
             thalf = told + dt2 !# midpoint in time.
             call src3(meqn,mbc,mx,my,mz,
     &                 xlower,ylower,zlower,dx,dy,dz,
     &                 q,maux,aux,told,dt2)
             endif
c
c        # copy q into qwork1.  q is updated in step2 and qwork1 is
c        # preserved to provide data for Riemann problems.
c
         call copyq3(meqn,mbc,mx,my,mz,q,
     &               work(i0qwrk1))
c
c        # take one step on the conservation law:
c
         if(method(3) .ge. 0)then
c           # unsplit version
c
            call step3(maxm,meqn,mwaves,mbc,mx,my,
     &                 mz,work(i0qwrk1),q,aux,dx,dy,dz,dt,method,
     &                 mthlim,cfl,
     &                 work(i0qadd),work(i0fadd),work(i0gadd),
     &                 work(i0hadd),work(i0q1d),work(i0dtdx1d),
     &                 work(i0dtdy1d),work(i0dtdz1d),
     &                 work(i0aux1),work(i0aux2),work(i0aux3),maux,
     &                 work(i0next),mwork1,use_fwave,rpn3,rpt3, rptt3)
c
         else
c           # dimensional splitting (fractional steps)
c
            call dimsp3(maxm,meqn,mwaves,mbc,mx,my,
     &                  mz,work(i0qwrk1),q,aux,dx,dy,dz,dt,method,
     &                  mthlim,cfl,cflv,
     &                  work(i0qadd),work(i0fadd),work(i0gadd),
     &                  work(i0hadd),work(i0q1d),work(i0dtdx1d),
     &                  work(i0dtdy1d),work(i0dtdz1d),
     &                  work(i0aux1),work(i0aux2),work(i0aux3),maux,
     &                  work(i0next),mwork1,use_fwave,rpn3,rpt3,rptt3)
c
         endif
c
         t = told + dt
c
         if (method(4).eq.1) then
c            # verbose mode
             write(6,601) n,cfl,dt,t
  601        format('CLAW3... Step',i4,
     &                   '   Courant number =',f6.3,'  dt =',d12.4,
     &                   '  t =',d12.4)
             endif
c
c
c        # check to see if the Courant number was too large:
         if (cfl .le. cflv(1)) then
c             # accept this step
              cflmax = dmax1(cfl,cflmax)
            else
c             # Reject this step.  Reset q to qwork from previous time:
c             # Note that i0qwrk2 points to work space where previous
c             # solution is stored in both cases method(5) = 0 or 1.
              t = told
              call copyq3(meqn,mbc,mx,my,mz,
     &                    work(i0qwrk2),q)
c
              if (method(4).eq.1) then
c                # verbose mode
                 write(6,602)
  602            format('CLAW3 rejecting step... ',
     &                  'Courant number too large')
                 endif
              if (method(1).eq.1) then
c                 # if variable dt, go back and take a smaller step.
                  dt = dmin1(dtv(2), dt * cflv(2)/cfl)
                  go to 40
                else
c                 # if fixed dt, give up and return
                  cflmax = dmax1(cfl,cflmax)
                  go to 900
                endif
            endif
c
c        # claw3 step is accepted
c        # now apply source terms:
c
         if (method(5).eq.2) then
c            # source terms over a second half time step for Strang splitting:
c            # Note it is not so clear what time t should be used here if
c            # the source terms are time-dependent!
             call src3(meqn,mbc,mx,my,mz,
     &                 xlower,ylower,zlower,dx,dy,dz,
     &                 q,maux,aux,t,dt2)
             endif
c
         if (method(5).eq.1) then
c            # source terms over a full time step:
             call src3(meqn,mbc,mx,my,mz,
     &                 xlower,ylower,zlower,dx,dy,dz,
     &                 q,maux,aux,t,dt)
             endif
c
c        ================================================================
c
c
         if (method(1) .eq. 1) then
c           # choose new time step if variable time step
            if (cfl.eq.0.d0) then
                dt = dtv(2)
              else
                dt = dmin1(dtv(2), dt * cflv(2)/cfl)
              endif
            dtmin = dmin1(dt,dtmin)
            dtmax = dmax1(dt,dtmax)
            endif
c
c        # see if we are done:
c
         nv(2) = nv(2) + 1
         if (t .ge. tend) go to 900
  100    continue
c
  900  continue
c
c      # return information
c
       if (method(1).eq.1 .and. t.lt.tend .and. nv(2) .eq. maxn) then
c         # too many timesteps
          write(6,*) 'CLAW3 ERROR...  too many timesteps'
          info = 11
          endif
       if (method(1).eq.0 .and. cflmax .gt. cflv(1)) then
c         # Courant number too large with fixed dt
          write(6,*) 'CLAW3 ERROR...  Courant number too large'
          info = 12
          endif
       tend = t
       cflv(3) = cflmax
       cflv(4) = cfl
       dtv(3) = dtmin
       dtv(4) = dtmax
       dtv(5) = dt
c
       return
       end
c
c


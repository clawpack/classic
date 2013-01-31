! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module solution_module

    use precision_module

    implicit none

    type solution_type

        ! Solution extents
        integer :: num_dim, num_eqn, num_aux
        integer, allocatable :: num_cells(:)
        real(kind=8) :: t
        real(kind=8), allocatable :: dx(:), lower(:), upper(:), centers(:)

        ! Aux array descriptors
        integer :: capa_index

        ! Solution arrays
        real(kind=Q_TYPE), pointer :: q(:,:)
        real(kind=AUX_TYPE), pointer :: aux(:,:)
        
    end type solution_type

    interface new
        module procedure new_solution
    end interface

    interface output
        module procedure output_solution
    end interface

contains


    subroutine new_solution(self,clawdata)

        use clawdata_module, only: clawdata_type

        implicit none

        ! Input
        type(solution_type), intent(out) :: self
        type(clawdata_type) :: clawdata

        integer :: stat, i

        ! Set new solution object's array extents
        self%num_dim = clawdata%num_dim
        allocate(self%num_cells(self%num_dim), stat=stat)
        if (stat /= 0) stop "Allocation request denied for num_cells"
        self%num_cells = clawdata%num_cells
        self%num_eqn = clawdata%num_eqn
        self%num_aux = clawdata%num_aux

        ! Domain
        allocate(self%lower(self%num_dim), stat=stat)
        if (stat /= 0) print *, "Allocation request denied!"
        self%lower = clawdata%lower
        allocate(self%upper(self%num_dim), stat=stat)
        if (stat /= 0) print *, "Allocation request denied!"
        self%upper = clawdata%upper

        ! Set initial time
        self%t = clawdata%t0

        ! Various other parameters
        self%capa_index = clawdata%capa_index

        ! Allocate memory for the solution arrays
        associate(num_eqn => clawdata%num_eqn, &
                  num_ghost => clawdata%num_ghost, &
                  num_cells => clawdata%num_cells(1), &
                  num_aux => clawdata%num_aux)

            ! Allocate primary arrays
            allocate(self%q(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
            if (stat /= 0) stop "Allocation of solutions's q array failed!"
            allocate(self%aux(num_aux,1-num_ghost:num_cells+num_ghost),stat=stat)
            if (stat /= 0) stop "Allocation of solutions's aux array failed!"

            ! Calculated values
            allocate(self%dx(self%num_dim), stat=stat)
            if (stat /= 0) stop "Allocation of solution's dx array failed!"
            self%dx(1) = (self%upper(1) - self%lower(1)) / num_cells
            allocate(self%centers(1-num_ghost:num_cells+num_ghost),stat=stat)
            if (stat /= 0) stop "Allocation of solutions's centers array failed!"
            forall(i=1-num_ghost:num_cells+num_ghost)
                self%centers(i) = self%lower(1) + (i-0.5d0) * self%dx(1)
            end forall

        end associate
        
    end subroutine new_solution


    subroutine output_solution(solution,t,frame,q_components,aux_components)

        implicit none

        ! Arguments
        type(solution_type), intent(in out) :: solution
        real(kind=8), intent(in) :: t
        integer, intent(in) :: frame
        integer, intent(in) :: q_components(solution%num_eqn)
        integer, intent(in) :: aux_components(solution%num_aux)

        ! Locals
        character(len=128) :: q_file_name, t_file_name, aux_file_name
        integer, parameter :: IOUNIT = 50
        integer :: i, m

        ! Format strings
        character(len=9) :: q_output_format
        character(len=*), parameter :: Q_HEADER_FORMAT = &
                                    "(i5,'                 grid_number',/," // &
                                     "i5,'                 AMR_level',/,"   // &
                                     "i5,'                 mx',/,"          // &
                                     "e18.8,'    xlow', /,"                 // &
                                     "e18.8,'    dx', /)"

        character(len=*), parameter :: T_HEADER_FORMAT = &
                                    "(e26.16,'    time', /,"                // &
                                         "i5,'                 meqn'/,"     // &
                                         "i5,'                 ngrids'/,"   // &
                                         "i5,'                 maux'/,"     // &
                                         "i5,'                 ndim'/,/)"

        write(q_output_format,"('(',i2,'e16.8)')") solution%num_eqn


        ! Construct paths to output files and open them
        q_file_name = construct_file_name(frame,'fort.q')
        open(IOUNIT,file=q_file_name,status='unknown',form='formatted')

        ! Write q file header
        write(IOUNIT,Q_HEADER_FORMAT)                  1,    &
                                                       1,    &
                                      solution%num_cells,    &
                                          solution%lower,    &
                                             solution%dx

        ! Loop through q array fixing any exponents with more than 2 digits in 
        ! exponent cause problems, reset tiny value to zero
        forall(i=1:solution%num_cells(1), m=1:solution%num_eqn, &
               abs(solution%q(m,i)) < 1d-99)
            solution%q(m,i) = 0.d0
        end forall

        ! Write out each line
        do i=1,solution%num_cells(1)
            write(IOUNIT,q_output_format) (solution%q(m,i), m=1,solution%num_eqn)
        enddo

        close(IOUNIT)

        ! Write out fort.t file                                         
        t_file_name = construct_file_name(frame,'fort.t')
        open(IOUNIT,file=t_file_name,status='unknown',form='formatted')

        write(IOUNIT,T_HEADER_FORMAT)                t,     &
                                      solution%num_eqn,     &
                                                     1,     &
                                      solution%num_aux,     &
                                                     1
        
        close(IOUNIT)

        ! Output the aux array if requested
        if (.not.all(aux_components == 0)) then
            stop "Outputting the aux array not supported yet!"
        end if
        
    end subroutine output_solution



    character(len=10) function construct_file_name(frame,prefix) &
                                            result(file_name)

        implicit none

        ! Arguments
        integer, intent(in) :: frame
        character(len=6), intent(in) :: prefix

        ! Locals
        integer :: position, nstop, digits

        file_name = prefix // "xxxx"
        nstop = frame
        do position = 10, 7, -1
            digits = mod(nstop,10)
            file_name(position:position) = char(ichar('0') + digits)
            nstop = nstop / 10
        enddo

    end function construct_file_name


end module solution_module

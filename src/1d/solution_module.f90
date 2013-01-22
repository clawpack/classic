! ==============================================================================
!          Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module solution_module

    implicit none
    
    ! Array types - Default is double precision as defined by 1.d0
    integer, parameter :: Q_TYPE = kind(1.d0)
    integer, parameter :: AUX_TYPE = kind(1.d0)

    type solution_type

        ! Solution extents
        integer :: num_cells, num_eqn, num_aux
        real(kind=8) :: t, dx, lower, upper

        ! Aux array descriptors
        integer :: capa_index

        ! Solution arrays
        real(kind=Q_TYPE), pointer :: q(:,:)
        real(kind=AUX_TYPE), pointer :: aux(:,:)
        
    end type solution_type

contains



    type(solution_type) function new_solution(clawdata) result(solution)

        use clawdata_module, only: clawdata_type

        implicit none

        type(clawdata_type), intent(in) :: clawdata

        integer :: stat

        ! Set new solution object's array extents
        solution%num_cells = clawdata%num_cells
        solution%num_eqn = clawdata%num_eqn
        solution%num_aux = clawdata%num_aux
        solution%lower = clawdata%lower
        solution%upper = clawdata%upper
        solution%dx = (solution%upper - solution%lower) / solution%num_cells

        ! Set initial time
        solution%t = clawdata%t0

        ! Various other parameters - This needs to be set by the user
        solution%capa_index = clawdata%capa_index

        ! Allocate memory for the solution arrays
        associate(num_eqn => clawdata%num_eqn, &
                  num_ghost => clawdata%num_ghost, &
                  num_cells => clawdata%num_cells, &
                  num_aux => clawdata%num_aux)

            allocate(solution%q(num_eqn,1-num_ghost:num_cells+num_ghost),stat=stat)
            if (stat /= 0) stop "Allocation of solutions's q array failed!"
            allocate(solution%aux(num_aux,1-num_ghost:num_cells+num_ghost),stat=stat)
            if (stat /= 0) stop "Allocation of solutions's aux array failed!"

        end associate
        
    end function new_solution




    subroutine output_solution(solution,t,frame,output_aux,path)

        implicit none

        ! Arguments
        type(solution_type), intent(in out) :: solution
        real(kind=8), intent(in) :: t
        integer, intent(in) :: frame
        character(len=*), optional, intent(in) :: path
        logical, intent(in) :: output_aux

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

        q_output_format = "(" // char(solution%num_eqn) // "e16.8)"

        ! Construct paths to output files and open them
        q_file_name = path // construct_file_name(frame,'fort.q')
        open(IOUNIT,file=q_file_name,status='unknown',form='formatted')

        ! Write q file header
        write(IOUNIT,Q_HEADER_FORMAT)                  1,    &
                                                       1,    &
                                      solution%num_cells,    &
                                          solution%lower,    &
                                             solution%dx

        ! Loop through q array fixing any exponents with more than 2 digits in 
        ! exponent cause problems, reset tiny value to zero
        forall(i=1:solution%num_cells, m=1:solution%num_eqn, &
               abs(solution%q(m,i)) < 1d-99)
            solution%q(m,i) = 0.d0
        end forall

        ! Write out each line
        do i=1,solution%num_cells
            write(IOUNIT,q_output_format) (solution%q(m,i), m=1,solution%num_eqn)
        enddo

        close(IOUNIT)

        ! Write out fort.t file                                         
        t_file_name = path // construct_file_name(frame,'fort.t')
        open(IOUNIT,file=t_file_name,status='unknown',form='formatted')

        write(IOUNIT,T_HEADER_FORMAT)                t,     &
                                      solution%num_eqn,     &
                                                     1,     &
                                      solution%num_aux,     &
                                                     1
        
        close(IOUNIT)

        ! Output the aux array if requested
        if (output_aux) then
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

module utils
implicit none
private
public open_data_file

contains

subroutine open_data_file(iounit, fname)
! Open the file fname and determine how many leading lines are
! comments.  Then rewind and skip over the comment lines so that the
! file is ready for reading data from.
!
! All comment lines must start with # in the first column.


    ! Arguments
    integer, intent(in) :: iounit

    ! This allows for a variable name length and obsolves (I think) the
    ! length problem that appeared for some compilers (notably ifort)
    character(len=*), intent(in) :: fname

    ! Local
    integer :: line, commentlines
    character(len=1) :: firstchar
    logical :: foundFile
    

    inquire(file=fname,exist=foundFile)
    if (.not. foundFile) then
        ! # for backward compatability: maybe fname is declared character*12
        ! # in calling routine...
        ! fname12 = fname(1:12)
        ! #write(6,*) 'truncated fname12 = XXX',fname12,'XXX'
        ! inquire(file=fname12,exist=foundFile)
        if (.not. foundFile) then
            print "(2a)",'*** in opendatafile, file not found:', fname 
            stop
        endif
        ! open(unit=iounit,file=fname12,status='old',form='formatted')
        ! write(6,*) 'Reading data file: ', fname12
    else
        open(unit=iounit,file=fname,status='old',form='formatted')
        print "(2a)",'Reading data file: ', fname
    endif

    ! # this version may not work in f77
    ! firstchar = '#'
    ! commentlines = -1
    ! do while (firstchar == '#')
    !    read(iounit,*) firstchar
    !    commentlines = commentlines + 1
    !    enddo

    firstchar = '#'
    commentlines = -1
    do commentlines=-1,1000
        if (firstchar .eq. '#') then
            read(iounit,"(a1)") firstchar
        else
            exit
        endif
    enddo
      
    print "('         first',i2,' lines are comments and will be skipped')", &
        commentlines
     
    ! Rewind file and go back to data line
    !   print *,firstchar
    !   rewind(iounit)
    !   do line=1,commentlines
    !      read(iounit,*) firstchar
    !      print *,firstchar
    !      enddo
    !   

end subroutine open_data_file

end module

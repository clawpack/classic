module utils
use precision_module, only: dp
implicit none
private
public stop_error, str, open_data_file

interface str
    module procedure str_int, str_real, str_real_n
end interface

! This subroutine is defined in open_data_file.f90 as a global subroutine, so
! we interface is here is that user code can import it from here:
interface

    subroutine open_data_file(iounit, fname)
        implicit none
        integer, intent(in) :: iounit
        character(len=*), intent(in) :: fname
    end subroutine

end interface

contains

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine

pure integer function str_int_len(i) result(sz)
! Returns the length of the string representation of 'i'
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, '(i0)') i
sz = len_trim(s)
end function

pure function str_int(i) result(s)
! Converts integer "i" to string
integer, intent(in) :: i
character(len=str_int_len(i)) :: s
write(s, '(i0)') i
end function

pure integer function str_real_len(r, fmt) result(sz)
! Returns the length of the string representation of 'i'
real(dp), intent(in) :: r
character(len=*), intent(in) :: fmt
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fmt) r
sz = len_trim(s)
end function

pure function str_real(r) result(s)
! Converts the real number "r" to string with 7 decimal digits.
real(dp), intent(in) :: r
character(len=*), parameter :: fmt="(f0.6)"
character(len=str_real_len(r, fmt)) :: s
write(s, fmt) r
end function

pure function str_real_n(r, n) result(s)
! Converts the real number "r" to string with 'n' decimal digits.
real(dp), intent(in) :: r
integer, intent(in) :: n
character(len=str_real_len(r, "(f0." // str_int(n) // ")")) :: s
write(s, "(f0." // str_int(n) // ")") r
end function

end module

program test_solution

!     use solution_module
    
    implicit none

    type grid_type
        public
        integer :: stuff
        private
        integer :: private_stuff
    end type
    
    type(grid_type) :: my_grid

    print *,"Starting tests of solution structure..."

    my_grid%stuff = 1
    my_grid%private_stuff = 2

end program test_solution

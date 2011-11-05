module grid_module

    implicit none
    
    type grid_type
        integer, pointer :: n(:) => null()
    end type grid_type
    type grid_type_pointer
        type(grid_type), pointer :: grid => null()
    end type grid_type_pointer
    
    interface new
        module procedure new_grid
    end interface new
    interface delete
        module procedure delete_grid
    end interface delete
    
contains
    
    function new_grid(n)

        implicit none
        
        integer, intent(in) :: n(:)
        type(grid_type), pointer :: new_grid
        
        allocate(new_grid)
        allocate(new_grid%n(size(n)))
        new_grid%n = n
        
    end function new_grid
    
    subroutine delete_grid(grid)

        implicit none
        type(grid_type), pointer, intent(inout) :: grid
        
        deallocate(grid%n)
        deallocate(grid)

    end subroutine delete_grid

end module grid_module

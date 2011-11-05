module state_2d_module

    use grid_module

    implicit none

    type state_type
        private
        type(grid_type), pointer :: grid
        double precision, pointer :: q(:,:,:)
        double precision, pointer :: aux(:,:,:)
    end type state_type
    type state_type_pointer
        type(state_type), pointer :: state
    end type state_type_pointer
    
contains
    
    function new_state(meqn,maux,grid)
    
        implicit none
        
        integer, intent(in) :: meqn,maux
        type(grid_type), pointer :: grid
        type(state_type), pointer :: new_state
        
        ! Local storage
        integer :: err
        
        allocate(new_state)
        allocate(new_state%q(meqn,grid%n(1),grid%n(2)))
        allocate(new_state%aux(maux,grid%n(1),grid%n(2)))
        
    end function new_state
    
    subroutine delete_state(state)

        implicit none
        type(state_type), intent(inout) :: state
        
        deallocate(state%q)
        deallocate(state%aux)
        deallocate(state)
        
    end subroutine delete_state
    
    function get_q(state) result(q)
        implicit none
        type(state_type_pointer), intent(in) :: state
        double precision, pointer, intent(out) :: q(:,:,:)
        
        q => state%q
        
    end function get_q
    
    function get_aux(state) result(aux)
        implicit none
        type(state_type_pointer), intent(in) :: state
        double precision, pointer, intent(out) :: aux(:,:,:)
        
        aux => state%aux
        
    end function get_aux
end module state_2d_module

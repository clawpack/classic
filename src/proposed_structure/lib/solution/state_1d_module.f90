module state_1d_module

    use grid_module

    implicit none

    type state_type_pointer
        type(state_type), pointer :: state => null()
    end type state_type_pointer

    type state_type
        type(grid_type), pointer :: grid => null()
        double precision, pointer :: q(:,:) => null()
        double precision, pointer :: aux(:,:) => null()
    end type state_type
    
    interface get_q
        module procedure get_state_pointer_q
        module procedure get_state_q
    end interface get_q
    
    interface get_aux
        module procedure get_state_pointer_aux
        module procedure get_state_aux
    end interface get_aux
    
contains
    
    function new_state(meqn,maux,grid)
    
        implicit none
        
        integer, intent(in) :: meqn,maux
        type(grid_type), pointer, intent(in) :: grid
        
        type(state_type), pointer :: new_state
        
        allocate(new_state)
        new_state%grid => grid
        allocate(new_state%q(meqn,grid%n(1)))
        allocate(new_state%aux(maux,grid%n(1)))
        
    end function new_state
    
    subroutine delete_state(state)

        implicit none
        type(state_type), pointer, intent(inout) :: state
        
        deallocate(state%q)
        deallocate(state%aux)
        deallocate(state)
        
    end subroutine delete_state
    
    function get_state_pointer_q(self) result(q)

        implicit none
        type(state_type_pointer), intent(in) :: self
        double precision, pointer :: q(:,:)
        
        q => self%state%q

    end function get_state_pointer_q
    
    function get_state_q(self) result(q)
    
        implicit none
        type(state_type), intent(in) :: self
        double precision, pointer :: q(:,:)
    
        q => self%q

    end function get_state_q
    
    function get_state_pointer_aux(self) result(aux)

        implicit none
        type(state_type_pointer), intent(in) :: self
        double precision, pointer :: aux(:,:)
        
        aux => self%state%aux

    end function get_state_pointer_aux

    function get_state_aux(self) result(aux)

        implicit none
        type(state_type), intent(in) :: self
        double precision, pointer :: aux(:,:)
        
        aux => self%aux

    end function get_state_aux

end module state_1d_module

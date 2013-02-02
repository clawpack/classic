subroutine limiter(num_cells, num_ghost, num_eqn, num_waves, wave, s, mthlim)
    
    implicit none

    ! Input
    integer :: num_cells, num_ghost, num_eqn, num_waves
    real(kind=8), intent(in out) :: wave(num_eqn, num_waves, 1-num_ghost:num_cells + num_ghost)
    real(kind=8), intent(in out) :: s(num_waves, 1-num_ghost:num_cells + num_ghost)
    integer, intent(in) :: mthlim(num_waves)

    ! Locals
    integer :: i, mw
    real(kind=8) :: dot_product_right, dot_product_left, wave_norm, philim, r

    do mw=1,num_waves
        if (mthlim(mw) /= 0) then
            do i=0,num_cells+1
                dot_product_left = dot_product_right
                dot_product_right = sum(wave(:,mw,i+1) * wave(:,mw,i))
                wave_norm = sum(wave(:,mw,i)**2)
                if (i == 0) cycle
                if (wave_norm == 0.d0) cycle

                if (s(mw,i) > 0.d0) then
                    r = dot_product_left / wave_norm
                else
                    r = dot_product_right / wave_norm
                endif

                select case(mthlim(mw))
                    ! Min-mod
                    case(1)
                        philim = max(0.d0, min(1.d0, r))

                    ! Superbee
                    case(2)
                        philim = max(0.d0, min(1.d0, 2.d0 * r), min(2.d0, r))

                    ! MC
                    case(3)
                        philim = max(0.d0, min((1.d0 + r) / 2.d0, 2.d0, 2.d0 * r))

                    ! Van Leer
                    case(4)
                        philim = (r + abs(r)) / (1.d0 + abs(r))

                    ! Beam-Warming
                    case(5)
                        philim = r

                    case default
                        stop "*** ERROR *** Invalid limiter requested."
                end select

                wave(:,mw,i) = philim * wave(:,mw,i)
            enddo
        endif
    enddo
    
end subroutine limiter
module field_update
    use setup_parameter
    use load_model
    implicit none

    double precision :: vin
    contains

    subroutine update_feed
        implicit none
        vin = sin(omega * t)
        ez(feedx, feedy, feedz) = vin / dz
    end subroutine update_feed


    subroutine update_inc_field
        implicit none
        
    end subroutine update_inc_field


    subroutine update_efield
        implicit none
        integer :: i, j, k
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1

                    ex(i, j, k) = cex(i, j, k) * ex(i, j, k) &
                                + dez(i, j, k) * (hz(i, j, k) - hz(i, j-1, k)) &
                                - dey(i, j, k) * (hy(i, j, k) - hy(i, j, k-1))

                    ey(i, j, k) = cey(i, j, k) * ey(i, j, k) &
                                + dex(i, j, k) * (hx(i, j, k) - hx(i, j, k-1)) &
                                - dez(i, j, k) * (hz(i, j, k) - hz(i-1, j, k))

                    ez(i, j, k) = cez(i, j, k) * ez(i, j, k) &
                                + dey(i, j, k) * (hy(i, j, k) - hy(i-1, j, k)) &
                                - dex(i, j, k) * (hx(i, j, k) - hx(i, j-1, k))

                    ex(i, j, k) = idpecx(i, j, k) * ex(i, j, k)
                    ey(i, j, k) = idpecy(i, j, k) * ey(i, j, k)
                    ez(i, j, k) = idpecz(i, j, k) * ez(i, j, k)

                end do
            end do
        end do
    end subroutine update_efield


    subroutine update_hfield
        implicit none
        integer :: i, j, k
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1

                    hx(i, j, k) = chx(i, j, k) * hx(i, j, k) &
                                + dhz(i, j, k) * (ez(i, j, k) - ez(i, j+1, k)) &
                                - dhy(i, j, k) * (ey(i, j, k) - ey(i, j, k+1))

                    hy(i, j, k) = chy(i, j, k) * hy(i, j, k) &
                                + dhx(i, j, k) * (ex(i, j, k) - ex(i, j, k+1)) &
                                - dhz(i, j, k) * (ez(i, j, k) - ez(i+1, j, k))

                    hz(i, j, k) = chz(i, j, k) * hz(i, j, k) &
                                + dhy(i, j, k) * (ey(i, j, k) - ey(i+1, j, k)) &
                                - dhx(i, j, k) * (ex(i, j, k) - ex(i, j+1, k))

                end do
            end do
        end do
    end subroutine update_hfield

end module field_update



module calc_sar
    implicit none
    contains

    subroutine calc_ave_sar_wb
        implicit none
    
    end subroutine calc_ave_sar_wb


    subroutine calc_peak_sar_xg
        implicit none
    
    end subroutine calc_peak_sar_xg

end module calc_sar



module calc_scatter
    implicit none
    contains

    subroutine scatter_efield
        implicit none
    
    end subroutine scatter_efield

    subroutine calc_return_voltage
        implicit none
    
    end subroutine calc_return_voltage
end module calc_scatter
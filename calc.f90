module field_update
    use setup_parameter
    use load_model
    implicit none

    contains

    subroutine update_Vfeed
        implicit none
        vfeed = sin(omega * t)
        ez(feedx, feedy, feedz) = vfeed / dz

        vfeed_sub = cos(omega * t)
        ez_sub(feedx, feedy, feedz) = vfeed_sub / dz
    end subroutine update_Vfeed


    subroutine update_Ifeed
        implicit none
        ifeed = dy * (hx(feedx, feedy - 1, feedz) - hx(feedx, feedy, feedz)) &
              + dx * (hy(feedx, feedy, feedz) - hy(feedx - 1, feedy, feedz))
    end subroutine update_Ifeed


    subroutine update_inc_field
        implicit none
        
    end subroutine update_inc_field


    subroutine update_efield
        implicit none
        integer :: i, j, k
        double precision :: cex_buf, cey_buf, cez_buf
        double precision :: dex_buf, dey_buf, dez_buf
        integer :: idpecx_buf, idpecy_buf, idpecz_buf

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k, cex_buf, cey_buf, cez_buf, dex_buf, dey_buf, dez_buf, idpecx_buf, idpecy_buf, idpecz_buf)
        !$omp do
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1

                    !配列アクセス減らすため...
                    cex_buf = cex(i, j, k)
                    cey_buf = cey(i, j, k)
                    cez_buf = cez(i, j, k)
                    dex_buf = dex(i, j, k)
                    dey_buf = dey(i, j, k)
                    dez_buf = dez(i, j, k)
                    idpecx_buf = idpecx(i, j, k)
                    idpecy_buf = idpecy(i, j, k)
                    idpecz_buf = idpecz(i, j, k)

                    ex(i, j, k) = cex_buf * ex(i, j, k) &
                                + dez_buf * (hz(i, j, k) - hz(i, j-1, k)) &
                                - dey_buf * (hy(i, j, k) - hy(i, j, k-1))

                    ey(i, j, k) = cey_buf * ey(i, j, k) &
                                + dex_buf * (hx(i, j, k) - hx(i, j, k-1)) &
                                - dez_buf * (hz(i, j, k) - hz(i-1, j, k))

                    ez(i, j, k) = cez_buf * ez(i, j, k) &
                                + dey_buf * (hy(i, j, k) - hy(i-1, j, k)) &
                                - dex_buf * (hx(i, j, k) - hx(i, j-1, k))

                    ex(i, j, k) = idpecx_buf * ex(i, j, k)
                    ey(i, j, k) = idpecy_buf * ey(i, j, k)
                    ez(i, j, k) = idpecz_buf * ez(i, j, k)

                    !#####

                    ex_sub(i, j, k) = cex_buf * ex_sub(i, j, k) &
                                    + dez_buf * (hz_sub(i, j, k) - hz_sub(i, j-1, k)) &
                                    - dey_buf * (hy_sub(i, j, k) - hy_sub(i, j, k-1))

                    ey_sub(i, j, k) = cey_buf * ey_sub(i, j, k) &
                                    + dex_buf * (hx_sub(i, j, k) - hx_sub(i, j, k-1)) &
                                    - dez_buf * (hz_sub(i, j, k) - hz_sub(i-1, j, k))

                    ez_sub(i, j, k) = cez_buf * ez_sub(i, j, k) &
                                    + dey_buf * (hy_sub(i, j, k) - hy_sub(i-1, j, k)) &
                                    - dex_buf * (hx_sub(i, j, k) - hx_sub(i, j-1, k))

                    ex_sub(i, j, k) = idpecx(i, j, k) * ex_sub(i, j, k)
                    ey_sub(i, j, k) = idpecy(i, j, k) * ey_sub(i, j, k)
                    ez_sub(i, j, k) = idpecz(i, j, k) * ez_sub(i, j, k)

                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

    end subroutine update_efield


    subroutine update_hfield
        implicit none
        integer :: i, j, k
        double precision :: chx_buf, chy_buf, chz_buf
        double precision :: dhx_buf, dhy_buf, dhz_buf

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k, chx_buf, chy_buf, chz_buf, dhx_buf, dhy_buf, dhz_buf)
        !$omp do
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    
                    !配列アクセス減らすため...
                    chx_buf = chx(i, j, k)
                    chy_buf = chy(i, j, k)
                    chz_buf = chz(i, j, k)
                    dhx_buf = dhx(i, j, k)
                    dhy_buf = dhy(i, j, k)
                    dhz_buf = dhz(i, j, k)

                    hx(i, j, k) = chx_buf * hx(i, j, k) &
                                + dhz_buf * (ez(i, j, k) - ez(i, j+1, k)) &
                                - dhy_buf * (ey(i, j, k) - ey(i, j, k+1))

                    hy(i, j, k) = chy_buf * hy(i, j, k) &
                                + dhx_buf * (ex(i, j, k) - ex(i, j, k+1)) &
                                - dhz_buf * (ez(i, j, k) - ez(i+1, j, k))

                    hz(i, j, k) = chz_buf * hz(i, j, k) &
                                + dhy_buf * (ey(i, j, k) - ey(i+1, j, k)) &
                                - dhx_buf * (ex(i, j, k) - ex(i, j+1, k))

                    !#####

                    hx_sub(i, j, k) = chx_buf * hx_sub(i, j, k) &
                                    + dhz_buf * (ez_sub(i, j, k) - ez_sub(i, j+1, k)) &
                                    - dhy_buf * (ey_sub(i, j, k) - ey_sub(i, j, k+1))

                    hy_sub(i, j, k) = chy_buf * hy_sub(i, j, k) &
                                    + dhx_buf * (ex_sub(i, j, k) - ex_sub(i, j, k+1)) &
                                    - dhz_buf * (ez_sub(i, j, k) - ez_sub(i+1, j, k))

                    hz_sub(i, j, k) = chz_buf * hz_sub(i, j, k) &
                                    + dhy_buf * (ey_sub(i, j, k) - ey_sub(i+1, j, k)) &
                                    - dhx_buf * (ex_sub(i, j, k) - ex_sub(i, j+1, k))

                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine update_hfield

end module field_update



module calc_amp
    use setup_parameter
    use load_model
    implicit none
    contains

    ! this subroutine is called at convergence check times
    subroutine calc_field_amp
        implicit none
        integer :: i, j, k

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k)
        !$omp do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    
                    examp(i, j, k) = sqrt(ex(i, j, k)*ex(i, j, k) + ex_sub(i, j, k)*ex_sub(i, j, k))
                    eyamp(i, j, k) = sqrt(ey(i, j, k)*ey(i, j, k) + ey_sub(i, j, k)*ey_sub(i, j, k))
                    ezamp(i, j, k) = sqrt(ez(i, j, k)*ez(i, j, k) + ez_sub(i, j, k)*ez_sub(i, j, k))

                    hxamp(i, j, k) = sqrt(hx(i, j, k)*hx(i, j, k) + hx_sub(i, j, k)*hx_sub(i, j, k))
                    hyamp(i, j, k) = sqrt(hy(i, j, k)*hy(i, j, k) + hy_sub(i, j, k)*hy_sub(i, j, k))
                    hzamp(i, j, k) = sqrt(hz(i, j, k)*hz(i, j, k) + hz_sub(i, j, k)*hz_sub(i, j, k))
                    
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine calc_field_amp

    subroutine calc_sar
        implicit none
        integer :: i, j, k
        integer :: idbuf
        double precision :: ex_buf, ey_buf, ez_buf

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k, idbuf, ex_buf, ey_buf, ez_buf)
        !$omp do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    idbuf = idper(i, j, k)
                    ex_buf = examp(i, j, k)
                    ey_buf = eyamp(i, j, k)
                    ez_buf = ezamp(i, j, k)

                    sar(i, j, k) = (ex_buf * ex_buf + ey_buf * ey_buf + ez_buf * ez_buf) &
                                 * sigma(idbuf) / rho(idbuf)
                    
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine calc_sar

    ! this subroutine is called at convergence check times
    subroutine calc_ave_sar_wb
        implicit none
    
    end subroutine calc_ave_sar_wb

    ! this subroutine is called at convergence check times
    subroutine calc_peak_sar_xg
        implicit none
    
    end subroutine calc_peak_sar_xg

    ! this subroutine is called at end of programs
    subroutine calc_field_phase
        implicit none
        integer :: i, j, k

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k)
        !$omp do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    
                    exphase(i, j, k) = atan2(ex(i, j, k), ex_sub(i, j, k))
                    eyphase(i, j, k) = atan2(ey(i, j, k), ey_sub(i, j, k))
                    ezphase(i, j, k) = atan2(ez(i, j, k), ez_sub(i, j, k))

                    hxphase(i, j, k) = atan2(hx(i, j, k), hx_sub(i, j, k))
                    hyphase(i, j, k) = atan2(hy(i, j, k), hy_sub(i, j, k))
                    hzphase(i, j, k) = atan2(hz(i, j, k), hz_sub(i, j, k))
                    
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine calc_field_phase

end module calc_amp



module calc_scatter
    implicit none
    contains
    
    ! this subroutine is called at end of programs
    subroutine scatter_efield
        implicit none
    
    end subroutine scatter_efield

    ! this subroutine is called at end of programs
    subroutine calc_return_voltage
        implicit none
    
    end subroutine calc_return_voltage
end module calc_scatter
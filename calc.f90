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
        integer :: i, j, k

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k)
        !$omp do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    einx(i, j, k) = einx(i, j, k) * cphase
                    einy(i, j, k) = einy(i, j, k) * cphase
                    einz(i, j, k) = einz(i, j, k) * cphase

                    einx_sub(i, j, k) = einx_sub(i, j, k) * cphase
                    einy_sub(i, j, k) = einy_sub(i, j, k) * cphase
                    einz_sub(i, j, k) = einz_sub(i, j, k) * cphase
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine update_inc_field


    subroutine update_scatter_efield
        implicit none
        integer :: i, j, k
        integer :: idperx_buf, idpery_buf, idperz_buf

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k, idperx_buf, idpery_buf, idperz_buf)
        !$omp do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    idperx_buf = idperx(i, j, k)
                    idpery_buf = idpery(i, j, k)
                    idperz_buf = idperz(i, j, k)
                    
                    ex(i, j, k) = ex(i, j, k) - (sigma(idperx_buf) * imag(einx(i, j, k)) &
                                + (eps(idperx_buf) - eps0) * imag(im * omega * einx(i, j, k))) &
                                * dex(i, j, k) * dx

                    ey(i, j, k) = ey(i, j, k) - (sigma(idpery_buf) * imag(einy(i, j, k)) &
                                + (eps(idpery_buf) - eps0) * imag(im * omega * einy(i, j, k))) &
                                * dey(i, j, k) * dy

                    ez(i, j, k) = ez(i, j, k) - (sigma(idperz_buf) * imag(einz(i, j, k)) &
                                + (eps(idperz_buf) - eps0) * imag(im * omega * einz(i, j, k))) &
                                * dez(i, j, k) * dz
                                        
                    
                    ex_sub(i, j, k) = ex_sub(i, j, k) - (sigma(idperx_buf) * imag(einx_sub(i, j, k)) &
                                    + (eps(idperx_buf) - eps0) * imag(im * omega * einx_sub(i, j, k))) &
                                    * dex(i, j, k) * dx

                    ey_sub(i, j, k) = ey_sub(i, j, k) - (sigma(idpery_buf) * imag(einy_sub(i, j, k)) &
                                    + (eps(idpery_buf) - eps0) * imag(im * omega * einy_sub(i, j, k))) &
                                    * dey(i, j, k) * dy

                    ez_sub(i, j, k) = ez_sub(i, j, k) - (sigma(idperz_buf) * imag(einz_sub(i, j, k)) &
                                    + (eps(idperz_buf) - eps0) * imag(im * omega * einz_sub(i, j, k))) &
                                    * dez(i, j, k) * dz
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

    end subroutine update_scatter_efield

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

                    ex_sub(i, j, k) = idpecx_buf * ex_sub(i, j, k)
                    ey_sub(i, j, k) = idpecy_buf * ey_sub(i, j, k)
                    ez_sub(i, j, k) = idpecz_buf * ez_sub(i, j, k)

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

    subroutine calc_total_field
        implicit none
        integer :: i, j, k

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k)
        !$omp do
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
                    etx(i, j, k) &
                        = (ex(i, j, k) + ex(i, j, k + 1) + ex(i, j + 1, k) + ex(i, j + 1, k + 1)) * 0.25d0 &
                        + imag(einx(i, j, k) + einx(i, j, k + 1) &
                        +      einx(i, j + 1, k) + einx(i, j + 1, k + 1)) * 0.25d0
                    
                    ety(i, j, k) &
                        = (ey(i, j, k) + ey(i, j, k + 1) + ey(i, j + 1, k) + ey(i, j + 1, k + 1)) * 0.25d0 &
                        + imag(einy(i, j, k) + einy(i, j, k + 1) & 
                        +      einy(i, j + 1, k) + einy(i, j + 1, k + 1)) * 0.25d0
                    
                    etz(i, j, k) &
                        = (ez(i, j, k) + ez(i, j, k + 1) + ez(i, j + 1, k) + ez(i, j + 1, k + 1)) * 0.25d0 &
                        + imag(einz(i, j, k) + einz(i, j, k + 1) &
                        +      einz(i, j + 1, k) + einz(i, j + 1, k + 1)) * 0.25d0
            

                    etx_sub(i, j, k) &
                        = (ex_sub(i, j, k) + ex_sub(i, j, k + 1) + ex_sub(i, j + 1, k) + ex_sub(i, j + 1, k + 1)) * 0.25d0 &
                        + imag(einx_sub(i, j, k) + einx_sub(i, j, k + 1) &
                        +      einx_sub(i, j + 1, k) + einx_sub(i, j + 1, k + 1)) * 0.25d0
                    
                    ety_sub(i, j, k) &
                        = (ey_sub(i, j, k) + ey_sub(i, j, k + 1) + ey_sub(i, j + 1, k) + ey_sub(i, j + 1, k + 1)) * 0.25d0 &
                        + imag(einy_sub(i, j, k) + einy_sub(i, j, k + 1) &
                        +      einy_sub(i, j + 1, k) + einy_sub(i, j + 1, k + 1)) * 0.25d0
                    
                    etz_sub(i, j, k) &
                        = (ez_sub(i, j, k) + ez_sub(i, j, k + 1) + ez_sub(i, j + 1, k) + ez_sub(i, j + 1, k + 1)) * 0.25d0 &
                        + imag(einz_sub(i, j, k) + einz_sub(i, j, k + 1) &
                        +      einz_sub(i, j + 1, k) + einz_sub(i, j + 1, k + 1)) * 0.25d0
            
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine calc_total_field
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
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
                    
                    examp(i, j, k) = sqrt(etx(i, j, k) * etx(i, j, k) + etx_sub(i, j, k) * etx_sub(i, j, k))
                    eyamp(i, j, k) = sqrt(ety(i, j, k) * ety(i, j, k) + ety_sub(i, j, k) * ety_sub(i, j, k))
                    ezamp(i, j, k) = sqrt(etz(i, j, k) * etz(i, j, k) + etz_sub(i, j, k) * etz_sub(i, j, k))

                    hxamp(i, j, k) = sqrt(hx(i, j, k) * hx(i, j, k) + hx_sub(i, j, k) * hx_sub(i, j, k))
                    hyamp(i, j, k) = sqrt(hy(i, j, k) * hy(i, j, k) + hy_sub(i, j, k) * hy_sub(i, j, k))
                    hzamp(i, j, k) = sqrt(hz(i, j, k) * hz(i, j, k) + hz_sub(i, j, k) * hz_sub(i, j, k))
                    
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
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
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
    subroutine calc_body_sar
        implicit none
        integer :: i, j, k
        integer :: idbuf

        sar_ave_wb = 0.0d0
        mass_weight = 0.0d0
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
                    idbuf = idper(i, j, k)
                    sar_ave_wb = sar_ave_wb + sar(i, j, k) * rho(idbuf)
                    mass_weight = mass_weight + dx * dy * dz * rho(idbuf)
                end do
            end do
        end do
        sar_ave_wb = sar_ave_wb / mass_weight
    end subroutine calc_body_sar

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
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
                    
                    exphase(i, j, k) = -omega * t * atan2(etx(i, j, k), etx_sub(i, j, k))
                    eyphase(i, j, k) = -omega * t * atan2(ety(i, j, k), ety_sub(i, j, k))
                    ezphase(i, j, k) = -omega * t * atan2(etz(i, j, k), etz_sub(i, j, k))

                    hxphase(i, j, k) = -omega * t * atan2(hx(i, j, k), hx_sub(i, j, k))
                    hyphase(i, j, k) = -omega * t * atan2(hy(i, j, k), hy_sub(i, j, k))
                    hzphase(i, j, k) = -omega * t * atan2(hz(i, j, k), hz_sub(i, j, k))
                    
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine calc_field_phase

end module calc_amp


module calc_scatter
    use setup_parameter
    use load_model
    implicit none
    contains

    ! this subroutine is called at end of programs
    subroutine calc_return_voltage
        implicit none
    
    end subroutine calc_return_voltage
end module calc_scatter
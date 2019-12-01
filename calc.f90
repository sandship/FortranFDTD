module field_update
    use setup_parameter
    use load_model
    implicit none

    contains
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

                    eamp(i, j, k) = sqrt(examp(i, j, k) * examp(i, j, k) &
                                       + eyamp(i, j, k) * eyamp(i, j, k) &
                                       + ezamp(i, j, k) * ezamp(i, j, k))

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
                    sar(i, j, k) = eamp(i, j, k) * eamp(i, j, k) * sigma(idbuf) / rho(idbuf)
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
        double precision :: total_weight = 0.d0

        sar_ave_wb = 0.0d0
        total_weight = 0.0d0
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
                    idbuf = idper(i, j, k)
                    sar_ave_wb = sar_ave_wb + sar(i, j, k) * rho(idbuf)
                    total_weight = total_weight + dx * dy * dz * rho(idbuf)
                end do
            end do
        end do
        sar_ave_wb = sar_ave_wb / total_weight
    end subroutine calc_body_sar

    ! this subroutine is called at convergence check times
    subroutine calc_peak_sar_xg(mass_weight)
        implicit none
        double precision, intent(in) :: mass_weight
        double precision :: part_weight
        double precision :: sar_mass, sar_mass_max
        double precision :: part_power
        double precision :: threshold_sar = 1.0d-14
        integer :: i, j, k
        integer :: ci, cj, ck ! cube center i, j, k

        integer :: n, np
        integer :: idbuf = 0
        integer, dimension(nx) :: iair = 0
        integer, dimension(ny) :: jair = 0
        integer, dimension(nz) :: kair = 0

        print *, 'calc. peak sar'

        do ci = npml, nz - npml
            do cj = npml, ny - npml
                do ck = npml, nx - npml
                    if (idper(ci, cj, ck) .ne. 1) then
                        idbuf = idper(ci, cj, ck)
                        n = 0
                        part_weight = rho(idbuf) * dx * dy * dz
                        sar_mass = 0.0d0

                        do while (part_weight < mass_weight)
                            do k = ck - n, ck + n
                                do j = cj - n, cj + n
                                    do i = ci - n, ci + n
                                        if (sar(i, j, k) > threshold_sar) then
                                            idbuf = idper(i, j, k)
                                            part_weight = part_weight + rho(idbuf) * dx * dy * dz
                                            part_power = part_power + sigma(idbuf) * eamp(i, j, k) * eamp(i, j, k)
                                        else
                                            iair(i) = iair(i) + 1
                                            jair(i) = jair(i) + 1
                                            kair(i) = kair(i) + 1
                                        end if
                                    end do
                                end do
                            end do
                            
                        end do

                    end if
                end do
            end do
        end do

    end subroutine calc_peak_sar_xg


    ! this subroutine is called at end of programs
    subroutine calc_field_phase
        implicit none
        integer :: i, j, k

        print *, 'calc. field phase'

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k)
        !$omp do
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
                    
                    exphase(i, j, k) = - omega * t * atan2(etx(i, j, k), etx_sub(i, j, k))
                    eyphase(i, j, k) = - omega * t * atan2(ety(i, j, k), ety_sub(i, j, k))
                    ezphase(i, j, k) = - omega * t * atan2(etz(i, j, k), etz_sub(i, j, k))

                    hxphase(i, j, k) = - omega * t * atan2(hx(i, j, k), hx_sub(i, j, k))
                    hyphase(i, j, k) = - omega * t * atan2(hy(i, j, k), hy_sub(i, j, k))
                    hzphase(i, j, k) = - omega * t * atan2(hz(i, j, k), hz_sub(i, j, k))
                    
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
    subroutine calc_escatter
        implicit none
        integer :: i, j, k
        integer :: n
        integer :: idbuf
        double precision :: rpx, rpy, rpz, rp
        complex(kind(0d0)) :: ppx, ppy, ppz, ds

        print *, 'calc. scatter efield'

        allocate(escatter_x(vertex_num))
        allocate(escatter_y(vertex_num))
        allocate(escatter_z(vertex_num))

        escatter_x = (0.d0, 0.d0)
        escatter_y = (0.d0, 0.d0)
        escatter_z = (0.d0, 0.d0)

        !$omp parallel num_threads(16) &
        !$omp private(i, j, k)
        !$omp do
        do k = npml, nz - npml
            do j = npml, ny - npml
                do i = npml, nx - npml
                    idbuf = idper(i ,j ,k)
                    
                    cur_x(i ,j ,k) = (im * omega * (eps(idbuf) - eps0) + sigma(idbuf)) &
                                   * examp(i, j, k) * exp(im * exphase(i, j, k))* dx * dy * dz

                    cur_y(i ,j ,k) = (im * omega * (eps(idbuf) - eps0) + sigma(idbuf)) &
                                   * eyamp(i, j, k) * exp(im * eyphase(i, j, k))* dx * dy * dz

                    cur_z(i ,j ,k) = (im * omega * (eps(idbuf) - eps0) + sigma(idbuf)) &
                                   * ezamp(i, j, k) * exp(im * ezphase(i, j, k))* dx * dy * dz
                    
                    xi(i, j, k) = i * dx
                    yi(i, j, k) = j * dy
                    zi(i, j, k) = k * dz

                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        !$omp parallel num_threads(16) &
        !$omp private(n, i, j, k, rpx, rpy, rpz, rp, ds, ppx, ppy, ppz)
        !$omp do
        do n = 1, vertex_num
            do k = npml, nz - npml
                do j = npml, ny - npml
                    do i = npml, nx - npml
                                         
                    rpx = (pp(1, n) - xi(i, j, k))
                    rpy = (pp(2, n) - yi(i, j, k))
                    rpz = (pp(3, n) - zi(i, j, k))
             
                    rp = sqrt(rpx ** 2 + rpy ** 2 + rpz ** 2)
                    ds = (1.0d0 + 1.0d0/(im * wn * rp)) / rp
             
                    ppx = (rpx * cur_x(i, j, k) &
                        +  rpy * cur_y(i, j, k) &
                        +  rpz * cur_z(i, j, k)) * rpx / rp**2
             
             
                    ppy = (rpx * cur_x(i, j, k) &
                        +  rpy * cur_y(i, j, k) &
                        +  rpz * cur_z(i, j, k)) * rpy / rp**2
             
             
                    ppz = (rpx * cur_x(i, j, k) &
                        +  rpy * cur_y(i, j, k) &
                        +  rpz * cur_z(i, j, k)) * rpz / rp**2


                    escatter_x(n) = escatter_x(n) - im * wn * z0 * exp(-im * wn * rp)/(4.0d0 * pi * rp) &
                                  * ((cur_x(i, j, k) - ppx) &
                                  + ds / (im * wn) * (cur_x(i, j, k) - 3.0d0 * ppx))
             
                    escatter_y(n) = escatter_y(n) - im * wn * z0 * exp(-im * wn * rp)/(4.0d0 * pi * rp) &
                                  * ((cur_y(i, j, k) - ppy) &
                                  + ds / (im * wn) * (cur_y(i, j, k) - 3.0d0 * ppy))
                    
                    escatter_z(n) = escatter_z(n) - im * wn * z0 * exp(-im * wn * rp)/(4.0d0 * pi * rp) &
                                  * ((cur_z(i, j, k) - ppz) &
                                  + ds / (im * wn) * (cur_z(i, j, k) - 3.0d0 * ppz))
                    
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

    end subroutine calc_escatter


    subroutine calc_return_voltage
    implicit none
        double precision, dimension(3) :: xa, ya, za    
        integer :: i, n
        integer :: mt
        print *, 'calc. return voltage'

        allocate(return_v(vertex_num))
        return_v = (0.d0, 0.d0)

        do n = 1, edge_num
            do i = 1, 2
                mt = tl(i, lrp(n))
                if(mt == freen(1, i)) then
                    xa(1) = pp(1, mt) / lambda
                    ya(1) = pp(2, mt) / lambda
                    za(1) = pp(3, mt) / lambda
                else
                    xa(2) = pp(1, mt) / lambda
                    ya(2) = pp(2, mt) / lambda
                    za(2) = pp(3, mt) / lambda
                endif
            end do

            do i = 1, 2
                mt = tl(i, lrm(n))
                if(mt == freen(2, i)) then
                    xa(3) = pp(1, mt) / lambda
                    ya(3) = pp(2, mt) / lambda
                    za(3) = pp(3, mt) / lambda
                endif
            end do

            return_v(edp(n)) = return_v(edp(n)) &
                             + (escatter_x(edp(n)) * (xa(2) - xa(1)) &
                             +  escatter_y(edp(n)) * (ya(2) - ya(1)) &
                             +  escatter_z(edp(n)) * (za(2) - za(1))) * 0.5d0
            
            return_v(edp(n)) = return_v(edp(n)) &
                             + (escatter_x(edp(n)) * (xa(3) - xa(2)) &
                             +  escatter_y(edp(n)) * (ya(3) - ya(2)) &
                             +  escatter_z(edp(n)) * (za(3) - za(2))) * 0.5d0

        end do

        open(101, file='output/return_voltage.dat', position='append')
            do n = 1, edge_num
                write(101,'(i5, 2e18.10)') edp(n), real(return_v(edp(n))), imag(return_v(edp(n)))
            end do
        close(101)
    end subroutine calc_return_voltage

end module calc_scatter
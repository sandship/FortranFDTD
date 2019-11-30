module load_model
    use setup_parameter
    implicit none
    contains

    ! this subroutine LOAD human model, 
    ! each line format represents (x, y, z, permitivity_ID)
    ! permitivity ID is bound with `tissue_param.csv`
    subroutine load_human_model
        implicit none
        integer :: i, j, k, n, idbuf
        integer :: ios = 1

        open(101, file=human_model, status='old')
            do n = 1, nx * ny * nz
                read(101, *, iostat=ios) i, j, k, idbuf
                idper(i + margin, j + margin, k + margin) = idbuf
                
                idperx(i + margin, j + margin, k + margin) = idbuf
                idpery(i + margin, j + margin, k + margin) = idbuf
                idperz(i + margin, j + margin, k + margin) = idbuf

                ! idperx(i + margin, j + margin, k + margin + 1) = idbuf
                ! idpery(i + margin, j + margin, k + margin + 1) = idbuf
                
                ! idpery(i + margin + 1, j + margin, k + margin) = idbuf
                ! idperz(i + margin + 1, j + margin, k + margin) = idbuf

                ! idperz(i + margin, j + margin + 1, k + margin) = idbuf
                ! idperx(i + margin, j + margin + 1, k + margin) = idbuf

                ! idperx(i + margin, j + margin + 1, k + margin + 1) = idbuf
                ! idpery(i + margin + 1, j + margin, k + margin + 1) = idbuf
                ! idperz(i + margin + 1, j + margin + 1, k + margin) = idbuf
            end do
        close(101)
    
    end subroutine load_human_model

    subroutine load_tissue
        implicit none
        integer :: index, n
        integer :: ios = 1

        open(101, file=tissue_param, status='old')
            read(101, *)
            do n = 1, nmax_per
                read(101, *, iostat=ios) index, sigma(index), eps(index), rho(index)
                eps(index) = eps(index) * eps0
            end do
        close(101)

    end subroutine load_tissue

    
    subroutine make_groud_plane
        implicit none
        integer :: i, j, k
        k = ground_plane_height
        do j = 1, ny
            do i = i, nx
                idpecx(i, j, k) = 0
                idpecy(i, j, k) = 0
            end do
        end do
    end subroutine make_groud_plane


    subroutine set_pml_coefficient
        implicit none
        integer :: n
        double precision :: sigbuf
        double precision :: msigbuf

        double precision :: ce_value, de_value
        double precision :: ch_value, dh_value

        do n = 1, npml
            sigbuf = (((float(npml) - float(n)) / float(npml)) ** dimpml) * ((dimpml + 1) * (-log(refpml))/(2 * npml * dx * z0))
            msigbuf = mu0 / eps0 * sigbuf

            ce_value = (2.0d0 * eps0 - sigbuf * dt)/(2.0d0 * eps0 + sigbuf * dt)
            de_value = (2.0d0 * dt)/(2.0d0 * eps0 + sigbuf * dt)/dx

            ch_value = (2.0d0 * mu0 - msigbuf * dt)/(2.0d0 * mu0 + msigbuf * dt)
            dh_value = (2.0d0 * dt)/(2.0d0 * mu0 + msigbuf * dt)/dx

            ! cex(n, :, :) = ce_value
            cey(n, :, :) = ce_value
            cez(n, :, :) = ce_value

            cex(:, n, :) = ce_value
            ! cey(:, n, :) = ce_value
            cez(:, n, :) = ce_value

            cex(:, :, n) = ce_value
            cey(:, :, n) = ce_value
            ! cez(:, :, n) = ce_value

            ! cex(nx - (n - 1), :, :) = ce_value
            cey(nx - (n - 1), :, :) = ce_value
            cez(nx - (n - 1), :, :) = ce_value

            cex(:, ny - (n - 1), :) = ce_value
            ! cey(:, ny - (n - 1), :) = ce_value
            cez(:, ny - (n - 1), :) = ce_value

            cex(:, :, nz - (n - 1)) = ce_value
            cey(:, :, nz - (n - 1)) = ce_value
            ! cez(:, :, nz - (n - 1)) = ce_value
            
            !###
            ! dex(n, :, :) = de_value
            dey(n, :, :) = de_value
            dez(n, :, :) = de_value

            dex(:, n, :) = de_value
            ! dey(:, n, :) = de_value
            dez(:, n, :) = de_value

            dex(:, :, n) = de_value
            dey(:, :, n) = de_value
            ! dez(:, :, n) = de_value

            ! dex(nx - (n - 1), :, :) = de_value
            dey(nx - (n - 1), :, :) = de_value
            dez(nx - (n - 1), :, :) = de_value

            dex(:, ny - (n - 1), :) = de_value
            ! dey(:, ny - (n - 1), :) = de_value
            dez(:, ny - (n - 1), :) = de_value

            dex(:, :, nz - (n - 1)) = de_value
            dey(:, :, nz - (n - 1)) = de_value
            ! dez(:, :, nz - (n - 1)) = de_value
            
            !###
            ! chx(n, :, :) = ch_value
            chy(n, :, :) = ch_value
            chz(n, :, :) = ch_value

            chx(:, n, :) = ch_value
            ! chy(:, n, :) = ch_value
            chz(:, n, :) = ch_value

            chx(:, :, n) = ch_value
            chy(:, :, n) = ch_value
            ! chz(:, :, n) = ch_value

            ! chx(nx - (n - 1), :, :) = ch_value
            chy(nx - (n - 1), :, :) = ch_value
            chz(nx - (n - 1), :, :) = ch_value

            chx(:, ny - (n - 1), :) = ch_value
            ! chy(:, ny - (n - 1), :) = ch_value
            chz(:, ny - (n - 1), :) = ch_value

            chx(:, :, nz - (n - 1)) = ch_value
            chy(:, :, nz - (n - 1)) = ch_value
            ! chz(:, :, nz - (n - 1)) = ch_value
            
            !###
            ! dhx(n, :, :) = dh_value
            dhy(n, :, :) = dh_value
            dhz(n, :, :) = dh_value

            dhx(:, n, :) = dh_value
            ! dhy(:, n, :) = dh_value
            dhz(:, n, :) = dh_value

            dhx(:, :, n) = dh_value
            dhy(:, :, n) = dh_value
            ! dhz(:, :, n) = dh_value

            ! dhx(nx - (n - 1), :, :) = dh_value
            dhy(nx - (n - 1), :, :) = dh_value
            dhz(nx - (n - 1), :, :) = dh_value

            dhx(:, ny - (n - 1), :) = dh_value
            ! dhy(:, ny - (n - 1), :) = dh_value
            dhz(:, ny - (n - 1), :) = dh_value

            dhx(:, :, nz - (n - 1)) = dh_value
            dhy(:, :, nz - (n - 1)) = dh_value
            ! dhz(:, :, nz - (n - 1)) = dh_value

        end do
    end subroutine set_pml_coefficient

    
    subroutine set_em_coefficient
        implicit none
        integer :: i, j, k
        integer :: idx, idy, idz

        dt = 0.99d0 / (cv * sqrt(1.d0 / dx ** 2 + 1.d0 / dy ** 2 + 1.d0 / dz ** 2))
        hdt = dt / 2.0d0
        cphase = exp(im * omega * dt)

        do k = 1, nz
            do j = 1, ny
                do i = 1, nz
                    idx = idperx(i, j, k)
                    idy = idpery(i, j, k)
                    idz = idperz(i, j, k)

                    cex(i, j, k) = (2.0d0 * eps(idx) - sigma(idx) * dt)/(2.0d0 * eps(idx) + sigma(idx) * dt)
                    cey(i, j, k) = (2.0d0 * eps(idy) - sigma(idy) * dt)/(2.0d0 * eps(idy) + sigma(idy) * dt)
                    cez(i, j, k) = (2.0d0 * eps(idz) - sigma(idz) * dt)/(2.0d0 * eps(idz) + sigma(idz) * dt)

                    dex(i, j, k) = (2.0d0 * dt)/((2.0d0 * eps(idx) + sigma(idx) * dt) * dx)
                    dey(i, j, k) = (2.0d0 * dt)/((2.0d0 * eps(idy) + sigma(idy) * dt) * dy)
                    dez(i, j, k) = (2.0d0 * dt)/((2.0d0 * eps(idz) + sigma(idz) * dt) * dz)
                    
                    chx(i, j, k) = (2.0d0 * mu(idx) - msigma(idx) * dt)/(2.0d0 * mu(idx) + msigma(idx) * dt)
                    chy(i, j, k) = (2.0d0 * mu(idy) - msigma(idy) * dt)/(2.0d0 * mu(idy) + msigma(idy) * dt)
                    chz(i, j, k) = (2.0d0 * mu(idz) - msigma(idz) * dt)/(2.0d0 * mu(idz) + msigma(idz) * dt)

                    dhx(i, j, k) = (2.0d0 * dt)/((2.0d0 * mu(idx) + msigma(idx) * dt) * dx)
                    dhy(i, j, k) = (2.0d0 * dt)/((2.0d0 * mu(idy) + msigma(idy) * dt) * dy)
                    dhz(i, j, k) = (2.0d0 * dt)/((2.0d0 * mu(idz) + msigma(idz) * dt) * dz)

                end do
            end do
        end do
    end subroutine set_em_coefficient
end module load_model



module load_field
    use setup_parameter
    implicit none

    contains

    subroutine load_efield
        implicit none
        integer :: i, j, k, n
        double precision :: px, py, pz
        double precision :: einx_re, einx_imag
        double precision :: einy_re, einy_imag
        double precision :: einz_re, einz_imag
        open(101, file=inc_efield, status='old')
        
        do n = 1, nx * ny * nz
            read(101, *, iostat=ios) px, py, pz, &
                                     einx_re, einx_imag, &
                                     einy_re, einy_imag, &
                                     einz_re, einz_imag
            
            i = int((px + cent_x * dx) / dx) + 1
            j = int((py + cent_y * dy) / dy) + 1
            k = int((pz + cent_z * dz) / dz) + 1

            einx(i, j, k) = einx_re + im * einx_imag
            einy(i, j, k) = einy_re + im * einy_imag
            einz(i, j, k) = einz_re + im * einz_imag

            einx_sub(i, j, k) = im * einx(i, j, k)
            einy_sub(i, j, k) = im * einy(i, j, k)
            einz_sub(i, j, k) = im * einz(i, j, k)

        end do
        close(101)
    end subroutine load_efield

end module load_field
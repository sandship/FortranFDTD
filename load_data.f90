module load_model
    use setup_parameter
    implicit none
    contains

    ! this subroutine LOAD human model, 
    ! each line format represents (x, y, z, permitivity_ID)
    ! permitivity ID is bound with `tissue_param.csv`
    subroutine load_human
        implicit none
        integer :: i, j, k
    
    end subroutine load_human


    ! this subroutine MAKE test Mie sphere model
    ! permitivity ID is bound with `tissue_param.csv`
    subroutine make_mie_model
        implicit none
        integer :: i, j, k

        do k = cent_z - mie_radius, cent_z + mie_radius
            do j = cent_y - mie_radius, cent_y + mie_radius
                do i = cent_x - mie_radius, cent_x + mie_radius

                    if((k - cent_z) ** 2 + (j - cent_y) ** 2 + (i - cent_x) ** 2 <= mie_radius ** 2) then
                        idper(i, j, k) = mie_per

                        idperx(i, j, k) = mie_per
                        ! idperx(i, j, k + 1) = mie_per
                        ! idperx(i, j + 1, k) = mie_per
                        ! idperx(i, j + 1, k + 1) = mie_per

                        idpery(i, j, k) = mie_per
                        ! idpery(i, j, k + 1) = mie_per
                        ! idpery(i + 1, j, k) = mie_per
                        ! idpery(i + 1, j, k + 1) = mie_per

                        idperz(i, j, k) = mie_per
                        ! idperz(i, j + 1, k) = mie_per
                        ! idperz(i + 1, j, k) = mie_per
                        ! idperz(i + 1, j + 1, k) = mie_per
                        
                    end if

                end do
            end do
        end do

        eps(mie_per) = 77.9 * eps0
        sigma(mie_per) = 0.408d0
        rho(mie_per) = 1040.0d0

    end subroutine make_mie_model


    subroutine make_dipole_antenna
        implicit none
        integer :: i, j, k
        i = feedx
        j = feedy
        do k = feedz - int(mie_dipole_len / 2), feedz + int(mie_dipole_len / 2)
            idpecz(i, j, k) = 0
        end do
    end subroutine make_dipole_antenna


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

    end subroutine load_efield

    subroutine load_plane_wave
        implicit none
        integer :: i, j, k

        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nz - 1
                    einx(i, j, k) = exp(im * (k * dz * klambda))
                    einy(i, j, k) = 0.0d0
                    einz(i, j, k) = 0.0d0

                    einx_sub(i, j, k) = einx(i, j, k) * im
                    einy_sub(i, j, k) = einy(i, j, k) * im
                    einz_sub(i, j, k) = einz(i, j, k) * im
                end do
            end do
        end do

        open(101, file='output/load_planewave.dat', status='replace')
        i = cent_x
        do j = 1, ny
            do k = 1, nz
                write(101, *) j, k, real(einx(i, j, k)), imag(einx(i, j, k))
            end do
        end do
        close(101)

    end subroutine load_plane_wave

end module load_field
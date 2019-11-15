program main
    use fileIO
    use setup_parameter
    use load_model
    use load_field
    use field_update
    use calc_sar
    use calc_scatter
    implicit none

    integer :: i, j, k
    character(len=8) :: n

    call make_mie_model
    call make_dipole_antenna

    !call load_efield

    call set_em_coefficient
    call set_pml_coefficient

    do step = 1, nt
        call update_feed
        call update_efield
        t = t + hdt

        call update_feed
        call update_hfield
        t = t + hdt

        if (mod(step, check_interval) == 0) then
            call calc_ave_sar_wb
            call calc_peak_sar_xg

            write(n, '(I8)') step
            open(101, file='./output/'// n //'_Eamp_Instantaneous_YZplane.dat', status='replace')
            i = cent_x
            do k = npml, nz - npml
                do j = npml, ny - npml
                    write(101, *) j, k, ex(i, j, k), ey(i, j, k), ez(i, j, k), &
                                  sqrt(ex(i, j, k) ** 2 + ey(i, j, k) ** 2 + ez(i, j, k) ** 2)
                end do
            end do
            close(101)

        end if
        
    end do

    call scatter_efield
    call calc_return_voltage

    open(101, file='./output/CE_YZplane.dat', status='replace')
    i = cent_x
    do k = npml, nz - npml
        do j = npml, ny - npml
            write(101, *) j, k, cex(i, j, k), cey(i, j, k), cez(i, j, k) &
                              , dex(i, j, k), dey(i, j, k), dez(i, j, k)
        end do
    end do
    close(101)

    open(101, file='./output/CE_XZplane.dat', status='replace')
    j = cent_y
    do k = npml, nz - npml
        do i = npml, nx - npml
            write(101, *) i, k, cex(i, j, k), cey(i, j, k), cez(i, j, k) &
                              , dex(i, j, k), dey(i, j, k), dez(i, j, k)
        end do
    end do
    close(101)

    open(101, file='./output/Eamp_Instantaneous_Yline.dat', status='replace')
    i = cent_x
    k = cent_z
    do j = npml, ny - npml
        write(101, *) j, ex(i, j, k), ey(i, j, k), ez(i, j, k), &
                      sqrt(ex(i, j, k) ** 2 + ey(i, j, k) ** 2 + ez(i, j, k) ** 2)
    end do
    close(101)

end program main
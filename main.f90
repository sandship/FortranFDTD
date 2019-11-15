program main
    use fileIO
    use setup_parameter
    use load_model
    use load_field
    use field_update
    use calc_sar
    use calc_scatter
    implicit none

    integer :: i

    ! call make_mie_model
    call make_dipole_antenna

    !call load_efield

    call set_em_coefficient
    call set_pml_coefficient

    do i = 1, nz
        print *, cex(cent_x, cent_y, i), dex(cent_x, cent_y, i), chx(cent_x, cent_y, i), dhx(cent_x, cent_y, i)
    end do

    do step = 1, nt
        call update_feed
        call update_efield
        t = t + hdt

        call update_feed
        call update_hfield
        t = t + hdt

        if (mod(step, check_interval) == 0) then
            
            do i = 1, nz
                print *, step, i * dx, ex(cent_x + 10, cent_y, i)
            end do

            call calc_ave_sar_wb
            call calc_peak_sar_xg
        end if
        

    end do

    call scatter_efield
    call calc_return_voltage

end program main
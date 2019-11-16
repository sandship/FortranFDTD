program main
    use fileIO
    use setup_parameter
    use load_model
    use load_field
    use field_update
    use calc_amp
    use calc_scatter
    use output
    implicit none

    call make_mie_model
    call make_dipole_antenna

    !call load_efield
    !call load_plane_wave

    call set_em_coefficient
    call set_pml_coefficient

    call output_at_first

    do step = 1, nt
        call update_Vfeed
        call update_Ifeed
        call update_efield
        t = t + hdt

        call update_Vfeed
        call update_Ifeed
        call update_hfield
        t = t + hdt

        if (mod(step, check_interval) == 0) then
            call calc_field_amp
            call calc_ave_sar_wb
            call calc_peak_sar_xg

            call output_at_checkpoint
        end if
        
    end do

    call calc_field_phase
    call scatter_efield
    call calc_return_voltage
    
    call output_at_end
end program main
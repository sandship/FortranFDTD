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

    !call load_efield
    call load_plane_wave

    call set_em_coefficient
    call set_pml_coefficient

    call output_at_first

    do step = 1, nt
        call update_efield
        call update_inc_field
        call update_scatter_efield
        t = t + hdt

        call update_hfield
        t = t + hdt

        if (mod(step, check_interval) == 0) then
            print *, step, ' step'
            call calc_total_field
            call calc_field_amp
            
            call calc_sar
            call calc_body_sar
            call calc_peak_sar_xg

            call output_at_checkpoint
        end if
        
    end do

    call calc_field_phase
    call calc_return_voltage
    
    call output_at_end
end program main